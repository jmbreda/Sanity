#include "gene_fit.hpp"
#include "fit_frac.hpp"
#include "helpers/math/digamma_trigamma.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// "SI" is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: c = cells, b = grid points of the gene variance v_g.
// Where the SI has n_g + 1 (Eqs. 16, 19, 26, 31, 44, 46) this code uses n_g.
// This is intentional: the prior over alpha_g differs from the one in the SI.
//
// The fit of one gene (get_gene_expression_level, at the end of this file) has two stages:
//   1. One pass over the grid of candidate variances v_b. At each grid point fit_fractions solves for q_g(v_b) and
//      the fractions f_gc, from which the log likelihood L(v_b) and the mean LTQ mu_g(v_b) follow. Only q_g, L and
//      mu_g are kept per grid point.
//   2. Either one grid point is selected and the cell estimates are made there (MLE, MAP, EAP), or a second pass over
//      the grid makes the cell estimates at every grid point and averages them over the posterior of v_g (MARG).
// Both branches make the cell estimates with the same two functions: reconstruct_fractions (fit_frac.cpp), which
// gives f_gc back from q_g, and compute_cell_estimates, which turns f_gc into delta_gc and its error bar.

GeneWorkspace::GeneWorkspace(int num_cells, int num_bins, VarianceMethod variance_method)
    : lambert_log_arg_base_c(num_cells), fractions_c(num_cells), log_fold_change_c(num_cells),
      inv_curv_diag_rank1_c(num_cells), inv_curv_diag_num_c(num_cells), inv_curv_diag_den2_c(num_cells),
      log_weighted_total_b(num_bins), log_likelihood_b(num_bins), mean_ltq_b(num_bins)
{
    if (variance_method == VarianceMethod::MARG)
    {
        log_fold_change_bc.resize(static_cast<std::size_t>(num_bins) * num_cells);
        log_fold_change_error_sq_bc.resize(static_cast<std::size_t>(num_bins) * num_cells);
    }
}

// Zero-count error bar (SI Eq. 38): the squared displacement epsilon_gc^2 above delta*_gc at which the log-likelihood has dropped by 1/2.
// log_fold_change = delta*_gc, gene_variance = v_g, gene_total = n_g and fraction = f_gc.
static double zero_count_log_fold_change_error_sq(double log_fold_change, double gene_variance, double gene_total, double fraction)
{
    // Find the positive displacement whose log-likelihood drop is 1/2.
    // Keep log(1 + f*(exp(e)-1)) rather than linearizing the logarithm.
    auto log_likelihood_drop = [&](double displacement) { // sigma_gc in SI Eq. 38
        double logarithm; // log(1 + f*(exp(e)-1)), written to stay finite when exp(e) overflows.
        if (displacement < 700)
        {
            logarithm = std::log1p(fraction * std::expm1(displacement));
        }
        else
        {
            logarithm = displacement + std::log(fraction + (1 - fraction) * std::exp(-displacement));
        }
        return displacement * (2 * log_fold_change + displacement) / (2 * gene_variance) + gene_total * logarithm;
    };
    double lo = 0, hi = 1;
    // The upper bound for the linearized formula need not bracket this root.
    while (log_likelihood_drop(hi) < 0.5)
    {
        hi *= 2;
        if (!std::isfinite(hi)) throw std::runtime_error("Cannot bracket zero-count uncertainty");
    }
    for (int iteration = 0; iteration < 200; ++iteration)
    {
        const double displacement = (lo + hi) / 2;
        const double value = log_likelihood_drop(displacement);
        if (!std::isfinite(value)) throw std::runtime_error("Nonfinite zero-count likelihood");
        if (std::fabs(value - 0.5) <= 1e-7) return displacement * displacement;
        if (value < 0.5) lo = displacement; else hi = displacement;
        if (hi - lo <= 1e-14 * std::max(1.0, hi)) return displacement * displacement;
    }
    throw std::runtime_error("Zero-count uncertainty did not converge");
}

// The log fold change of every cell at one v_g: delta_gc = log(f_gc) - log(N_c) + q_g (SI Eq. 28).
static void compute_log_fold_changes(const std::vector<double>& fractions_c, double log_weighted_total,
                                     const std::vector<double>& log_cell_totals_c, double* log_fold_change_c)
{
    const int num_cells = static_cast<int>(fractions_c.size());
    for (int c = 0; c < num_cells; ++c)
    {
        log_fold_change_c[c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total;
    }
}

// The variance of the log fold change of every cell at one v_g, var(delta_gc), from the diagonal of the inverse of
// M_g in closed form (SI Eqs. 35-37): var(delta_gc) = num_c / (den1 * den2_c). Cells without counts get the
// asymmetric error bar of SI Eq. 38 instead.
static void compute_log_fold_change_errors(const std::vector<double>& umi_counts_c, const std::vector<double>& fractions_c,
                                           const double* log_fold_change_c, double gene_total, double gene_variance,
                                           double* log_fold_change_error_sq_c, GeneWorkspace& workspace)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    double* rank1_c = workspace.inv_curv_diag_rank1_c.data();
    double* num_c = workspace.inv_curv_diag_num_c.data();
    double* den2_c = workspace.inv_curv_diag_den2_c.data();
    const double inv_variance = 1.0 / gene_variance;
    // rank1_c = n_g f_gc^2 / (n_g f_gc + 1/v_g) for each cell c.
    for (int c = 0; c < num_cells; ++c)
    {
        rank1_c[c] = gene_total * fractions_c[c] * fractions_c[c] / (gene_total * fractions_c[c] + inv_variance);
    }
    // den1 is one minus the full sum of the rank-one terms (the matrix determinant lemma, SI Eq. 33); den2_c is the diagonal of M_g.
    double den1 = 1.0;
    for (int c = 0; c < num_cells; ++c)
    {
        den1 -= rank1_c[c];
        den2_c[c] = gene_total * fractions_c[c] + inv_variance;
    }
    // The numerator removes the term of cell c from the sum.
    for (int c = 0; c < num_cells; c++)
    {
        num_c[c] = den1 + rank1_c[c];
    }
    for (int c = 0; c < num_cells; ++c)
    {
        log_fold_change_error_sq_c[c] = num_c[c] / (den1 * den2_c[c]);
    }
    // fix computation of asymmetric variance for zero count (SI Eq. 38)
    for (int c = 0; c < num_cells; ++c)
    {
        if (umi_counts_c[c] <= 0.5)
        {
            log_fold_change_error_sq_c[c] = zero_count_log_fold_change_error_sq(log_fold_change_c[c], gene_variance, gene_total, fractions_c[c]);
        }
    }
}

// The per-cell estimates at one v_g, given the fractions f_gc and q_g at that v_g: delta_gc and its variance.
// Called once at the selected grid point (MLE, MAP, EAP) or once per grid point (MARG).
static void compute_cell_estimates(const std::vector<double>& umi_counts_c, const std::vector<double>& fractions_c,
                                   double log_weighted_total, double gene_variance, double gene_total,
                                   const std::vector<double>& log_cell_totals_c,
                                   double* log_fold_change_c, double* log_fold_change_error_sq_c, GeneWorkspace& workspace)
{
    compute_log_fold_changes(fractions_c, log_weighted_total, log_cell_totals_c, log_fold_change_c);
    compute_log_fold_change_errors(umi_counts_c, fractions_c, log_fold_change_c, gene_total, gene_variance, log_fold_change_error_sq_c, workspace);
}

// L(v_b), the log marginal likelihood at one grid point (SI Eqs. 20, 32, 33), from f_gc, delta_gc and q_g at that point.
// The constant C*log(n_g) of log det(M_g) is dropped because it does not depend on v_g.
static double bin_log_likelihood(const std::vector<double>& umi_counts_c, const std::vector<double>& fractions_c,
                                 const std::vector<double>& log_fold_change_c, double gene_total, double bin_variance,
                                 double log_weighted_total)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const double inv_total_times_variance = 1.0 / (gene_total * bin_variance); // 1/(n_g v_b)
    double sum_log_fold_change_sq = 0; // sum_c delta_gc^2, from the Gaussian prior term of Eq. 20.
    double log_likelihood = -0.5 * ((double)num_cells) * std::log(bin_variance); // 1st term of SI Eq. 20
    for (int c = 0; c < num_cells; ++c)
    {
        log_likelihood += umi_counts_c[c] * log_fold_change_c[c];                 // data term sum_c n_gc*delta_gc, SI Eq. 20
        sum_log_fold_change_sq += log_fold_change_c[c] * log_fold_change_c[c];
    }
    log_likelihood -= sum_log_fold_change_sq / (2 * bin_variance); // 2nd term of SI Eq. 20
    log_likelihood -= gene_total * log_weighted_total;     // last term of SI Eq. 20, using sum_c N_c exp(delta_gc) = exp(q_g)

    // get the determinant of the matrix (SI Eq. 33)
    double log_det_curvature = 0.0; // log det(M_g/n_g)
    for (int c = 0; c < num_cells; ++c)
    {
        log_det_curvature += (fractions_c[c] * fractions_c[c]) / (fractions_c[c] + inv_total_times_variance);
    }
    log_det_curvature = std::log(1 - log_det_curvature);
    for (int c = 0; c < num_cells; ++c)
    {
        log_det_curvature += std::log(fractions_c[c] + inv_total_times_variance);
    }
    log_likelihood -= 0.5 * log_det_curvature; // SI Eq. 32
    return log_likelihood;
}

// P(v_b | n_g), the normalised posterior over the grid (SI Eq. 34), from L(v_b) and its maximum over the grid.
static void normalise_posterior(const std::vector<double>& log_likelihood_b, double max_log_likelihood,
                                std::vector<double>& variance_likelihood_b)
{
    const int num_bins = static_cast<int>(log_likelihood_b.size());
    double likelihood_sum = 0.0; // The normalisation, sum_b exp(L_b).
    for (int b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] = log_likelihood_b[b] - max_log_likelihood;
        variance_likelihood_b[b] = std::exp(variance_likelihood_b[b]);
        likelihood_sum += variance_likelihood_b[b];
    }
    for (int b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] /= likelihood_sum;
    }
}

// The grid point whose v_b is reported by the MLE, MAP and EAP methods.
static int select_variance_bin(VarianceMethod variance_method, const std::vector<double>& variance_likelihood_b,
                               const std::vector<double>& variance_grid_b, int max_likelihood_bin)
{
    const int num_bins = static_cast<int>(variance_grid_b.size());
    int selected_bin = 0;
    if (variance_method == VarianceMethod::MLE)
    {
        selected_bin = max_likelihood_bin;
    }
    else if (variance_method == VarianceMethod::MAP)
    {
        // Dividing the normalised likelihood by v applies an implicit 1/v (Jeffreys) prior on the
        // gene variance. This is the only prior in effect anywhere in this computation.
        double max_posterior = -1e+100;
        for (int b = 0; b < num_bins; b++)
        {
            const double posterior = variance_likelihood_b[b] / variance_grid_b[b];
            if (posterior > max_posterior)
            {
                max_posterior = posterior;
                selected_bin = b;
            }
        }
    }
    else if (variance_method == VarianceMethod::EAP)
    {
        double posterior_mean_variance = 0.0;
        for (int b = 0; b < num_bins; b++)
        {
            posterior_mean_variance += variance_likelihood_b[b] * variance_grid_b[b];
        }
        double min_distance = 1e+100;
        for (int b = 0; b < num_bins; b++)
        {
            const double distance = std::fabs(variance_grid_b[b] - posterior_mean_variance);
            if (distance < min_distance)
            {
                min_distance = distance;
                selected_bin = b;
            }
        }
    }
    else
    {
        throw std::logic_error("select_variance_bin: the marginalising method does not select a grid point");
    }
    return selected_bin;
}

// MLE, MAP and EAP: report v_b, mu_g(v_b) and the cell estimates at the selected grid point.
static void estimate_at_selected_bin(int selected_bin, const std::vector<double>& umi_counts_c, double gene_total, double trigamma_total,
                                     const std::vector<double>& variance_grid_b, const std::vector<double>& log_cell_totals_c,
                                     GeneWorkspace& workspace, RowComputation& result)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const double gene_variance = variance_grid_b[selected_bin];
    const double log_weighted_total = workspace.log_weighted_total_b[selected_bin];
    result.gene_variance = gene_variance;
    result.mean_ltq = workspace.mean_ltq_b[selected_bin];
    result.mean_ltq_error_sq = trigamma_total; // Only psi_1(n_g) (SI Eq. 47); the spread of mu_g(v_b) is not included.
    reconstruct_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, gene_variance, num_cells, log_weighted_total, log_cell_totals_c.data());
    compute_cell_estimates(umi_counts_c, workspace.fractions_c, log_weighted_total, gene_variance, gene_total, log_cell_totals_c,
                           result.log_fold_change_c.data(), result.log_fold_change_error_sq_c.data(), workspace);
}

// MARG: make the cell estimates at every grid point and average everything over the posterior P(v_b | n_g)
// (SI Eqs. 39, 42, 45 and 47). The reported v_g is the posterior mean of the variance.
static void marginalise_over_variance(const std::vector<double>& umi_counts_c, double gene_total, double trigamma_total, const std::vector<double>& variance_grid_b,
                                      const std::vector<double>& log_cell_totals_c, GeneWorkspace& workspace, RowComputation& result)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const int num_bins = static_cast<int>(variance_grid_b.size());
    const std::vector<double>& variance_likelihood_b = result.variance_likelihood_b;
    const std::vector<double>& mean_ltq_b = workspace.mean_ltq_b;
    std::vector<double>& log_fold_change_bc = workspace.log_fold_change_bc;
    std::vector<double>& log_fold_change_error_sq_bc = workspace.log_fold_change_error_sq_bc;

    // Second pass over the grid: delta_gc(v_b) and var(delta_gc)(v_b) at every grid point.
    for (int b = 0; b < num_bins; ++b)
    {
        const double bin_variance = variance_grid_b[b];
        const double log_weighted_total = workspace.log_weighted_total_b[b];
        reconstruct_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, bin_variance, num_cells, log_weighted_total, log_cell_totals_c.data());
        compute_cell_estimates(umi_counts_c, workspace.fractions_c, log_weighted_total, bin_variance, gene_total, log_cell_totals_c,
                               &log_fold_change_bc[static_cast<std::size_t>(b) * num_cells],
                               &log_fold_change_error_sq_bc[static_cast<std::size_t>(b) * num_cells], workspace);
    }

    // Average the mean LTQ over the variance grid (SI Eq. 45)
    double mean_ltq = 0.0;
    for (int b = 0; b < num_bins; b++)
    {
        mean_ltq += variance_likelihood_b[b] * mean_ltq_b[b];
    }
    // Compute (delta mu_g)^2 = psi_1 + < (mu - <mu>)^2 > (SI Eq. 47)
    double mean_ltq_error_sq = trigamma_total;
    for (int b = 0; b < num_bins; b++)
    {
        mean_ltq_error_sq += variance_likelihood_b[b] * (mean_ltq_b[b] - mean_ltq) * (mean_ltq_b[b] - mean_ltq);
    }
    result.mean_ltq = mean_ltq;
    result.mean_ltq_error_sq = mean_ltq_error_sq;

    // Compute <delta> = int p(v)*delta(v) dv (SI Eq. 39)
    std::vector<double>& log_fold_change_c = result.log_fold_change_c;
    for (int c = 0; c < num_cells; c++)
    {
        log_fold_change_c[c] = 0.0;
        for (int b = 0; b < num_bins; b++)
        {
            log_fold_change_c[c] += variance_likelihood_b[b] * log_fold_change_bc[static_cast<std::size_t>(b) * num_cells + c];
        }
    }
    // Compute epsilon^2 = < (delta - <delta>)^2 > + <var(delta)> (SI Eq. 42)
    std::vector<double>& log_fold_change_error_sq_c = result.log_fold_change_error_sq_c;
    for (int c = 0; c < num_cells; c++)
    {
        log_fold_change_error_sq_c[c] = 0.0;
        for (int b = 0; b < num_bins; b++)
        {
            const double log_fold_change = log_fold_change_bc[static_cast<std::size_t>(b) * num_cells + c];
            const double log_fold_change_error_sq = log_fold_change_error_sq_bc[static_cast<std::size_t>(b) * num_cells + c];
            log_fold_change_error_sq_c[c] += variance_likelihood_b[b] * (log_fold_change - log_fold_change_c[c]) * (log_fold_change - log_fold_change_c[c]) + variance_likelihood_b[b] * log_fold_change_error_sq;
        }
    }
    // The reported v_g is the posterior mean of the variance.
    double gene_variance = 0.0;
    for (int b = 0; b < num_bins; ++b)
    {
        gene_variance += variance_grid_b[b] * variance_likelihood_b[b];
    }
    result.gene_variance = gene_variance;
}

RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, double gene_total, const std::vector<double> &variance_grid_b,
                                         VarianceMethod variance_method, const std::vector<double>& log_cell_totals_c,
                                         double sum_cell_totals, GeneWorkspace& workspace)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const int num_bins = static_cast<int>(variance_grid_b.size());
    const double digamma_total = digamma(gene_total); // psi(n_g) (SI Eq. 44).
    const double trigamma_total = trigamma(gene_total); // psi_1(n_g), the variance of log(alpha_g) (SI Eq. 46).

    // Stage 1: fit every grid point, keeping q_g(v_b), L(v_b) and mu_g(v_b).
    double previous_log_weighted_total = 0.0; // q_g of the previous grid point, the starting guess for the next fit.
    double max_log_likelihood = -1e+100; // Maximum of L(v_b) over the grid.
    int max_likelihood_bin = 0; // Grid point at which L(v_b) is maximal (the MLE).
    for (int b = 0; b < num_bins; ++b)
    {
        const double bin_variance = variance_grid_b[b];
        const double log_weighted_total = fit_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, bin_variance, num_cells, previous_log_weighted_total, log_cell_totals_c.data(), sum_cell_totals, workspace.lambert_log_arg_base_c.data());
        previous_log_weighted_total = log_weighted_total;
        workspace.log_weighted_total_b[b] = log_weighted_total;
        workspace.mean_ltq_b[b] = digamma_total - log_weighted_total; /*** SI Eq. 44 ***/
        compute_log_fold_changes(workspace.fractions_c, log_weighted_total, log_cell_totals_c, workspace.log_fold_change_c.data());
        const double log_likelihood = bin_log_likelihood(umi_counts_c, workspace.fractions_c, workspace.log_fold_change_c, gene_total, bin_variance, log_weighted_total);
        workspace.log_likelihood_b[b] = log_likelihood;
        if (log_likelihood > max_log_likelihood)
        {
            max_log_likelihood = log_likelihood;
            max_likelihood_bin = b;
        }
    }

    RowComputation result;
    result.variance_likelihood_b.resize(num_bins);
    result.log_fold_change_c.resize(num_cells);
    result.log_fold_change_error_sq_c.resize(num_cells);
    normalise_posterior(workspace.log_likelihood_b, max_log_likelihood, result.variance_likelihood_b);

    // Stage 2: the reported estimates, either at one grid point or averaged over the grid.
    if (variance_method == VarianceMethod::MARG)
    {
        marginalise_over_variance(umi_counts_c, gene_total, trigamma_total, variance_grid_b, log_cell_totals_c, workspace, result);
    }
    else
    {
        const int selected_bin = select_variance_bin(variance_method, result.variance_likelihood_b, variance_grid_b, max_likelihood_bin);
        estimate_at_selected_bin(selected_bin, umi_counts_c, gene_total, trigamma_total, variance_grid_b, log_cell_totals_c, workspace, result);
    }
    return result;
}
