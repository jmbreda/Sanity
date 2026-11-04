#include "fit_frac.hpp"
#include "helpers/math/FukushimaLambertW.h"
#include <iostream>
#include <iomanip> // Required for setprecision
#include <cmath>
#include <optional>
#include <stdexcept>

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: c = cells.

// W_0(exp(log_arg)), the principal branch of the Lambert W function of exp(log_arg).
// At large arguments the LambertW routine breaks, so the asymptotic series is used instead.
static inline double lambert_w0_of_exp(double log_arg)
{
    if (log_arg > 50.0)
    {
        return lambert_w0_of_exp_asymptotic(log_arg);
    }
    else
    {
        return Fukushima::LambertW(0, std::exp(log_arg));
    }
}

// The residual n_g v_g - sum_c W[e^(Q_c - q)] of the consistency equation for q_g (SI Eq. 27, multiplied by n_g v_g),
// which is zero at the solution and decreases with q_g, and the Newton step residual/derivative at the same q_g.
// Both come from one evaluation of the Lambert W function per cell.
struct ResidualAndStep
{
    double residual;
    double newton_step;
};

static ResidualAndStep residual_and_step(const double *lambert_log_arg_base_c, int num_cells, double total_times_variance, double log_weighted_total)
{
    double residual = total_times_variance;
    double derivative = 0; // Minus the derivative of the residual with respect to q_g.
    for (int c = 0; c < num_cells; ++c)
    {
        const double lambert_log_arg = lambert_log_arg_base_c[c] - log_weighted_total;
        if (std::isnan(lambert_log_arg))
        {
            throw std::runtime_error("Error: the Lambert W argument is NaN in fit_fractions.");
        }
        const double lambert_w = lambert_w0_of_exp(lambert_log_arg);
        residual -= lambert_w;
        derivative += lambert_w / (1.0 + lambert_w);
    }
    return {residual, residual / derivative};
}

// Solves the optimum of the log-likelihood for one gene at one value of its variance v_g (SI Eqs. 21-28).
// Fills fractions_c with f_gc (SI Eq. 23) and returns log_weighted_total, q_g (SI Eq. 22).
double fit_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells,
                     std::optional<double> previous_log_weighted_total, const double* log_cell_totals_c, double sum_cell_totals, double* workspace_c)
{
    // Q_c = y_gc + log(n_g v_g) with y_gc = v_g n_gc + log(N_c) (SI Eq. 24), which is the log of the Lambert-W argument before q_g is subtracted.
    double *lambert_log_arg_base_c = workspace_c;
    // n_g v_g; the SI has v_g (n_g + 1) here (Eq. 26), because the prior over alpha_g differs.
    const double total_times_variance = gene_total * gene_variance;
    const double inv_total_times_variance = 1.0 / total_times_variance;
    const double log_total_times_variance = std::log(total_times_variance);
    for (int c = 0; c < num_cells; ++c)
    {
        lambert_log_arg_base_c[c] = log_cell_totals_c[c] + umi_counts_c[c] * gene_variance + log_total_times_variance;
    }

    // Newton's method on the residual converges if it starts where the residual is negative: the residual is a strictly
    // increasing and concave function of q, so each linear guess lies above the curve and zero is approached from below.
    // The fit of the previous grid point is the preferred start, after at most one Newton step to make the residual
    // negative. Otherwise start from log(sum_c N_c) + v_g/2: at the optimum q_g = log[sum_c N_c exp(delta_gc)] (SI Eq. 22),
    // and the expectation of exp(delta_gc) is exp(v_g/2).
    double log_weighted_total = 0.0;
    ResidualAndStep at_q{0.0, 0.0};
    bool started = false;
    if (previous_log_weighted_total)
    {
        log_weighted_total = *previous_log_weighted_total;
        at_q = residual_and_step(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total);
        if (at_q.residual > 0)
        {
            log_weighted_total -= at_q.newton_step;
            at_q = residual_and_step(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total);
        }
        started = at_q.residual <= 0;
    }
    if (!started)
    {
        log_weighted_total = std::log(sum_cell_totals) + 0.5 * gene_variance;
        at_q = residual_and_step(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total);
        while (at_q.residual > 0)
        {
            log_weighted_total -= 0.5 * gene_variance; // A rough step size: lowering q by the expected average of e^delta_gc.
            at_q = residual_and_step(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total);
        }
    }

    const double newton_tolerance = 1e-7;
    double step = at_q.newton_step;
    log_weighted_total -= step;
    while (std::fabs(step) > newton_tolerance && step < 0)
    {
        step = residual_and_step(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total).newton_step;
        log_weighted_total -= step;
    }

    double residual = total_times_variance;
    for (int c = 0; c < num_cells; ++c)
    {
        const double lambert_w = lambert_w0_of_exp(lambert_log_arg_base_c[c] - log_weighted_total);
        fractions_c[c] = lambert_w * inv_total_times_variance;
        residual -= lambert_w;
    }
    if (std::fabs(step / log_weighted_total) > 1e-6)
    {
        std::cerr << std::fixed << std::setprecision(12);
        std::cerr << "fitfrac: WARNING: dq = " << step << " is not close to zero at final q = " << log_weighted_total
            << " funcq = " << residual << std::endl;
    }
    return log_weighted_total;
}

// Asymptotic series for W_0(exp(log_arg)), used at large log_arg where the Lambert W routine breaks down.
// L2 is ln(ln(z)) for z = exp(log_arg), the standard variable of the series.
double lambert_w0_of_exp_asymptotic(const double& log_arg) {
    double L2 = std::log(log_arg);
    double log_arg_sq = log_arg*log_arg;
    double L2sq = L2*L2;
    double L2cub = L2sq*L2;
    double lambert_w = log_arg - L2 + L2/log_arg +
        L2*(L2 - 2.0)/(2.0*log_arg_sq) +
        L2*(6.0 - 9.0*L2 + 2.0*L2sq)/(6.0*log_arg_sq*log_arg) +
        L2*(-12.0 + 36.0*L2 - 22.0*L2sq + 3.0*L2cub)/(12.0*log_arg_sq*log_arg_sq);

    return lambert_w;
}

// Reconstruct fractions at the saved fitted root without repeating the root search.
void reconstruct_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, double log_weighted_total, const double* log_cell_totals_c) {
    const double total_times_variance = gene_total * gene_variance;
    const double inv_total_times_variance = 1.0 / total_times_variance;
    const double log_total_times_variance = std::log(total_times_variance);
    for (int c = 0; c < num_cells; ++c) {
        const double lambert_log_arg_base = log_cell_totals_c[c] + umi_counts_c[c] * gene_variance + log_total_times_variance;
        const double lambert_log_arg = lambert_log_arg_base - log_weighted_total;
        const double lambert_w = lambert_w0_of_exp(lambert_log_arg);
        fractions_c[c] = lambert_w * inv_total_times_variance;
    }
}
