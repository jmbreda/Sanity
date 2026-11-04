#include "gene_fit.hpp"
#include "fit_frac.hpp"
#include "helpers/math/digamma_trigamma.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

double zero_count_log_fold_change_error_sq(double log_fold_change, double gene_variance, double gene_total, double fraction);

// "SI" is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, b = grid points of the gene variance v_g.
// Where the SI has n_g + 1 (Eqs. 16, 19, 26, 31, 44, 46) this code uses n_g.
// This is intentional: the prior over alpha_g differs from the one in the SI.
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, const std::vector<double> &cell_totals_c,
                                         double gene_total, const std::vector<double> &variance_grid_b, int num_cells, int num_bins,
                                         int variance_method, const std::vector<double>& log_cell_totals_c,
                                         double sum_cell_totals, std::vector<double>& fitting_workspace_c)
{
    // umi_counts_c = n_gc, the counts of this gene in each cell.
    // cell_totals_c = N_c, the total counts of each cell.
    // gene_total = n_g, the total counts of this gene.
    // variance_grid_b = v_b, the candidate values of the gene variance v_g (SI S1.1.4).
    // num_cells = C and num_bins = B.
    int c, b;
    double inv_total_times_variance; // 1/(n_g v_b); the reciprocal of total_times_variance in fit_frac.cpp.
    double log_likelihood; // L(v_b), the log marginal likelihood of the current grid point (SI Eq. 32).
    double log_det_curvature; // log det(M_g/n_g) (SI Eqs. 31, 33); the constant C*log(n_g) is dropped because it does not depend on v_g.
    double log_weighted_total; // q_g = log sum_c N_c exp(delta_gc) (SI Eq. 22).
    double sum_log_fold_change_sq; // sum_c delta_gc^2, from the Gaussian prior term of Eq. 20.
    double previous_log_weighted_total = 0.0; // q_g of the previous grid point, the starting guess for the next fit.
    double *fractions_c = new double[num_cells]; // f_gc, the fitted fraction of the gene's counts expected in each cell (SI Eq. 23).
    double **log_fold_change_bc = new double *[num_bins]; // delta*_gc(v_b) at every grid point (SI Eq. 28), kept for MARG.
    double *log_fold_change_buffer_c = variance_method > 0 ? new double[num_cells] : nullptr; // Reusable delta*_gc for the methods that keep one bin.
    double **log_fold_change_error_sq_bc = new double *[num_bins](); // var(delta_gc)(v_b) at every grid point (SI Eq. 37), kept for MARG.
    // Fixed-bin methods retain only the fitted scalar; reconstruct cell values after selection.
    std::vector<double> log_weighted_total_b(num_bins); // q_g(v_b) at every grid point.
    std::vector<double> log_likelihood_b(num_bins, -1.0); // L_b = log P(n_g | v_b) (SI S1.1.4).
    std::vector<double> log_fold_change_c(num_cells, 0.0); // The reported log fold change <delta_gc> (SI Eq. 39).
    std::vector<double> log_fold_change_error_sq_c(num_cells, 0.0); // The reported error bar squared, epsilon_gc^2 (SI Eq. 42).
    double mean_ltq_error_sq; // (delta mu_g)^2, the variance of the mean LTQ (SI Eq. 47).
    double gene_variance; // The reported v_g.

    for (b = 0; b < num_bins; ++b)
    {
        log_fold_change_bc[b] = variance_method == 0 ? new double[num_cells] : log_fold_change_buffer_c;
        if (variance_method == 0) log_fold_change_error_sq_bc[b] = new double[num_cells];
    }

    /*** To compute the variance of the log fold change ***/
    // Together these give the diagonal of the inverse of M_g in closed form (SI Eqs. 35-37):
    // var(delta_gc) = inv_curv_diag_num_c / (inv_curv_diag_den1 * inv_curv_diag_den2_c).
    double *inv_curv_diag_rank1_c = new double[num_cells]; // n_g f_gc^2 / (n_g f_gc + 1/v_g), one term of each sum over cells in Eq. 37.
    double *inv_curv_diag_num_c = new double[num_cells]; // The numerator of Eq. 37: den1 + rank1_c, i.e. the sum with cell c removed.
    double *inv_curv_diag_den2_c = new double[num_cells]; // n_g f_gc + 1/v_g, the diagonal of M_g before the rank-one correction (Eq. 31).
    double inv_curv_diag_den1; // 1 - sum_c rank1_c, the factor from the matrix determinant lemma (Eq. 33).

    double *mean_ltq_b = new double[num_bins]; // mu_g(v_b) = psi(n_g) - q_g(v_b) at every grid point (SI Eq. 44).
    double max_log_likelihood = -1e+100; // Maximum of L_b over the grid.
    int max_likelihood_bin = 0; // Grid point at which L_b is maximal (the MLE).
    double bin_variance; // v_b, the gene variance at the current grid point.

    auto compute_log_fold_change_error_sq = [&](const double* given_fractions_c, const double* given_log_fold_change_c,
                                                double given_variance, double* error_sq_out_c) {
        // rank1_c = n_g f_gc^2 / (n_g f_gc + 1/v_g) for each cell c.
        const double inv_variance = 1.0 / given_variance;
        for (c = 0; c < num_cells; ++c)
        {
            inv_curv_diag_rank1_c[c] = gene_total * given_fractions_c[c] * given_fractions_c[c] / (gene_total * given_fractions_c[c] + inv_variance);
        }
        // den1 is one minus the full sum of the rank-one terms; den2_c is the diagonal of M_g.
        inv_curv_diag_den1 = 1.0;
        for (c = 0; c < num_cells; ++c)
        {
            inv_curv_diag_den1 -= inv_curv_diag_rank1_c[c];
            inv_curv_diag_den2_c[c] = gene_total * given_fractions_c[c] + inv_variance;
        }
        // The numerator removes the term of cell c from the sum.
        for (c = 0; c < num_cells; c++)
        {
            inv_curv_diag_num_c[c] = inv_curv_diag_den1 + inv_curv_diag_rank1_c[c];
        }
        // The variance of the log fold change of each cell.
        for (c = 0; c < num_cells; ++c)
        {
            error_sq_out_c[c] = inv_curv_diag_num_c[c] / (inv_curv_diag_den1 * inv_curv_diag_den2_c[c]);
        }

        // fix computation of asymmetric variance for zero count (SI Eq. 38)
        for (c = 0; c < num_cells; ++c)
        {
            if (umi_counts_c[c] <= 0.5)
            {
                error_sq_out_c[c] = zero_count_log_fold_change_error_sq(given_log_fold_change_c[c], given_variance, gene_total, given_fractions_c[c]);
            }
        }
    };

    const double digamma_total = digamma(gene_total); // psi(n_g) (SI Eq. 44).
    const double trigamma_total = trigamma(gene_total); // psi_1(n_g), the variance of log(alpha_g) (SI Eq. 46).
    for (b = 0; b < num_bins; ++b)
    {
        bin_variance = variance_grid_b[b];
        inv_total_times_variance = 1.0 / (gene_total * bin_variance);
        log_weighted_total = fit_fractions(fractions_c, umi_counts_c, gene_total, bin_variance, num_cells, cell_totals_c, previous_log_weighted_total, log_cell_totals_c.data(), sum_cell_totals, fitting_workspace_c.data());
        previous_log_weighted_total = log_weighted_total;
        log_weighted_total_b[b] = log_weighted_total;
        mean_ltq_b[b] = digamma_total - log_weighted_total; /*** SI Eq. 44 ***/

        sum_log_fold_change_sq = 0;
        log_likelihood = -0.5 * ((double)num_cells) * std::log(bin_variance); // 1st term of SI Eq. 20
        for (c = 0; c < num_cells; ++c)
        {
            log_fold_change_bc[b][c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total; // SI Eq. 28
            log_likelihood += umi_counts_c[c] * log_fold_change_bc[b][c];                 // data term sum_c n_gc*delta_gc, SI Eq. 20
            sum_log_fold_change_sq += log_fold_change_bc[b][c] * log_fold_change_bc[b][c];
        }
        log_likelihood -= sum_log_fold_change_sq / (2 * bin_variance); // 2nd term of SI Eq. 20
        log_likelihood -= gene_total * log_weighted_total;     // last term of SI Eq. 20, using sum_c N_c exp(delta_gc) = exp(q_g)

        // get the determinant of the matrix (SI Eq. 33)
        log_det_curvature = 0.0;
        for (c = 0; c < num_cells; ++c)
        {
            log_det_curvature += (fractions_c[c] * fractions_c[c]) / (fractions_c[c] + inv_total_times_variance);
        }

        log_det_curvature = std::log(1 - log_det_curvature);
        for (c = 0; c < num_cells; ++c)
        {
            log_det_curvature += std::log(fractions_c[c] + inv_total_times_variance);
        }
        log_likelihood -= 0.5 * log_det_curvature; // SI Eq. 32
        log_likelihood_b[b] = log_likelihood;

        if (log_likelihood > max_log_likelihood)
        {
            max_log_likelihood = log_likelihood;
            max_likelihood_bin = b;
        }

        if (variance_method == 0)
        {
            compute_log_fold_change_error_sq(fractions_c, log_fold_change_bc[b], bin_variance, log_fold_change_error_sq_bc[b]);
        }
    } // end v bins loop

    // get normalized likelihood from loglikelihood (SI Eq. 34)
    std::vector<double> variance_likelihood_b(num_bins); // P(v_b | n_g), the posterior over the grid (SI Eq. 34); written to likelihood.txt.
    double likelihood_sum; // The normalisation, sum_b exp(L_b).
    likelihood_sum = 0.0;
    for (b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] = log_likelihood_b[b] - max_log_likelihood;
        variance_likelihood_b[b] = std::exp(variance_likelihood_b[b]);
        likelihood_sum += variance_likelihood_b[b];
    }
    for (b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] /= likelihood_sum;
    }

    int selected_bin = 0; // Grid point chosen by the MLE, MAP or EAP methods.
    if(variance_method ==1){
        selected_bin = max_likelihood_bin;
    }
    else if(variance_method == 2){
        // Dividing the normalised likelihood by v applies an implicit 1/v (Jeffreys) prior on the
        // gene variance. This is the only prior in effect anywhere in this computation.
        double max_posterior = -1e+100;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];

            if (variance_likelihood_b[b]/grid_variance > max_posterior)
            {
                max_posterior = variance_likelihood_b[b]/grid_variance;
                selected_bin = b;
            }
        }
    }
    else if(variance_method == 3){
        double posterior_mean_variance = 0.0;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];
            posterior_mean_variance += variance_likelihood_b[b] * grid_variance;
        }
        double min_distance = 1e+100;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];
            if (std::fabs(grid_variance - posterior_mean_variance) < min_distance)
            {
                min_distance = std::fabs(grid_variance - posterior_mean_variance);
                selected_bin = b;
            }
        }
    }

    // Average the mean LTQ over the variance grid (SI Eq. 45)
    double mean_ltq = 0.0; // mu_g, the mean LTQ log(alpha_g) of the gene (SI S1.1.7).
    for (b = 0; b < num_bins; b++)
    {
        mean_ltq += variance_likelihood_b[b] * mean_ltq_b[b];
    }

    // Compute (delta mu_g)^2 = psi_1 + < (mu - <mu>)^2 > (SI Eq. 47)
    mean_ltq_error_sq = trigamma_total;
    for (b = 0; b < num_bins; b++)
    {
        mean_ltq_error_sq += variance_likelihood_b[b] * (mean_ltq_b[b] - mean_ltq) * (mean_ltq_b[b] - mean_ltq);
    }

    if (variance_method == 0)  // MARG
    {
        // Compute <delta> = int p(v)*delta(v) dv (SI Eq. 39)
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_c[c] = 0.0;
            for (b = 0; b < num_bins; b++)
            {
                log_fold_change_c[c] += variance_likelihood_b[b] * log_fold_change_bc[b][c];
            }
        }
        // Compute epsilon^2 = < (delta - <delta>)^2 > + <var(delta)> (SI Eq. 42)
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_error_sq_c[c] = 0.0;
            for (b = 0; b < num_bins; b++)
            {
                log_fold_change_error_sq_c[c] += variance_likelihood_b[b] * (log_fold_change_bc[b][c] - log_fold_change_c[c]) * (log_fold_change_bc[b][c] - log_fold_change_c[c]) + variance_likelihood_b[b] * log_fold_change_error_sq_bc[b][c];
            }
        }
        // The reported v_g is the posterior mean of the variance.
        gene_variance = 0.0;
        for (b = 0; b < num_bins; ++b)
        {
            bin_variance = variance_grid_b[b];
            gene_variance += bin_variance * variance_likelihood_b[b];
        }
    }

    if (variance_method > 0)
    {
        // Store the gene-variance that maximizes the likelihood:
        gene_variance = variance_grid_b[selected_bin];
        // And then also the corresponding values for the LTQs etc.
        mean_ltq = mean_ltq_b[selected_bin];
        mean_ltq_error_sq = trigamma_total;
        reconstruct_fractions(fractions_c, umi_counts_c, gene_total, gene_variance, num_cells, cell_totals_c, log_weighted_total_b[selected_bin], log_cell_totals_c.data());
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_c[c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total_b[selected_bin];
        }
        compute_log_fold_change_error_sq(fractions_c, log_fold_change_c.data(), gene_variance, log_fold_change_error_sq_c.data());
    }

    delete[] fractions_c;
    for (b = 0; b < num_bins; ++b)
    {
        if (variance_method == 0) delete[] log_fold_change_bc[b];
        delete[] log_fold_change_error_sq_bc[b];
    }
    delete[] log_fold_change_buffer_c;
    delete[] log_fold_change_bc;
    delete[] log_fold_change_error_sq_bc;
    delete[] inv_curv_diag_rank1_c;
    delete[] inv_curv_diag_num_c;
    delete[] inv_curv_diag_den2_c;
    delete[] mean_ltq_b;

    RowComputation result;
    result.mean_ltq = mean_ltq;
    result.mean_ltq_error_sq = mean_ltq_error_sq;
    result.log_fold_change_c = log_fold_change_c;
    result.log_fold_change_error_sq_c = log_fold_change_error_sq_c;
    result.gene_variance = gene_variance;
    result.variance_likelihood_b = variance_likelihood_b;
    
    return result;
}

// Zero-count error bar (SI Eq. 38): the squared displacement epsilon_gc^2 above delta*_gc at which the log-likelihood has dropped by 1/2.
// log_fold_change = delta*_gc, gene_variance = v_g, gene_total = n_g and fraction = f_gc.
double zero_count_log_fold_change_error_sq(double log_fold_change, double gene_variance, double gene_total, double fraction)
{
    // Find the positive displacement whose log-likelihood drop is 1/2.
    // Keep log(1 + f*(exp(e)-1)) rather than linearizing the logarithm.
    auto log_likelihood_drop = [&](double displacement) { // sigma_gc in SI Eq. 38
        const double logarithm = displacement < 700 ? std::log1p(fraction * std::expm1(displacement))
                                        : displacement + std::log(fraction + (1 - fraction) * std::exp(-displacement));
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
