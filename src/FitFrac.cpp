#include "FitFrac.h"
#include "FukushimaLambertW.h"
#include <iostream>
#include <iomanip> // Required for setprecision
#include <cmath>
#include <stdexcept>

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: c = cells.

// Solves the optimum of the log-likelihood for one gene at one value of its variance v_g (SI Eqs. 21-28).
// Fills fractions_c with f_gc (SI Eq. 23) and returns log_weighted_total, q_g (SI Eq. 22).
double fit_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, const std::vector<double>& cell_totals_c, double previous_log_weighted_total, const double* log_cell_totals_c, double sum_cell_totals, double* workspace_c) {
    double log_weighted_total, lambert_w, lambert_log_arg, step, residual;
    int c;
    // Q_c = y_gc + log(n_g v_g) with y_gc = v_g n_gc + log(N_c) (SI Eq. 24), which is the log of the Lambert-W argument before q_g is subtracted.
    double *lambert_log_arg_base_c = workspace_c ? workspace_c : new double [num_cells];
    // n_g v_g; the SI has v_g (n_g + 1) here (Eq. 26), because the prior over alpha_g differs.
    double total_times_variance = gene_total * gene_variance;
    double inv_total_times_variance = 1.0/total_times_variance;
    double log_total_times_variance = std::log(total_times_variance);

    log_weighted_total = 0.0;
    for(c=0;c<num_cells;++c){
        lambert_log_arg_base_c[c] = (log_cell_totals_c ? log_cell_totals_c[c] : std::log(cell_totals_c[c])) + umi_counts_c[c]*gene_variance + log_total_times_variance;
        if (!log_cell_totals_c) log_weighted_total += cell_totals_c[c]; // Accumulate sum_c N_c; the log is taken after the loop.
    }
    // Initial guess for q_g is log(sum_c N_c) + v_g/2.
    // At the optimum q_g = log[sum_c N_c exp(delta_gc)] (SI Eq. 22) and the expectation of exp(delta_gc) is exp(v_g/2).
    log_weighted_total = std::log(log_cell_totals_c ? sum_cell_totals : log_weighted_total) + 0.5*gene_variance;

    if (previous_log_weighted_total != 0.0) {
        residual = fraction_sum_residual(lambert_log_arg_base_c, num_cells, total_times_variance, previous_log_weighted_total);
        if (residual > 0) {
            // do a single Newton step
            step = newton_step(lambert_log_arg_base_c, num_cells, total_times_variance, previous_log_weighted_total);
            previous_log_weighted_total = previous_log_weighted_total - step;
            residual = fraction_sum_residual(lambert_log_arg_base_c, num_cells, total_times_variance, previous_log_weighted_total);
            if (residual > 0) {
                // no luck
                previous_log_weighted_total = 0.0;
            }
            else {
                log_weighted_total = previous_log_weighted_total;
            }
        }
        else {
            log_weighted_total = previous_log_weighted_total;
        }
    }

    if (previous_log_weighted_total == 0.0) {
        // We are looking for the zero of the residual n_g v_g - sum_c W[e^(Q_c-q)].
        // The residual is a strictly increasing and concave function of q.
        // Consequently, if we start at a value of q where it is negative, Newton's method is guaranteed to converge,
        // because the linear guess is always above the curve, so we approach zero from below.
        // So we check whether the residual is negative at our guess and, if not, keep lowering q until it is.
        residual = fraction_sum_residual(lambert_log_arg_base_c, num_cells, total_times_variance, log_weighted_total);
        while(residual > 0){
            log_weighted_total -= 0.5*gene_variance; // A rough step size: lowering q by the expected average of e^delta_gc.
            residual = fraction_sum_residual(lambert_log_arg_base_c , num_cells, total_times_variance, log_weighted_total);
        }
    }

    double newton_tolerance = 1e-7;
    step = -1.0;
    while(std::fabs(step) > newton_tolerance && step < 0){
        // step = residual/residual'.
        step = newton_step(lambert_log_arg_base_c,num_cells,total_times_variance,log_weighted_total);
        log_weighted_total = log_weighted_total - step; // Newton's method update.
    }

    residual = total_times_variance;

    for(c=0;c<num_cells;++c){
        lambert_log_arg = lambert_log_arg_base_c[c] -log_weighted_total;
        if(lambert_log_arg > 50){
            // At large arguments the LambertW routine breaks, so we use an approximation.
            lambert_w = lambert_w0_of_exp_asymptotic(lambert_log_arg);
        }
        else{
            lambert_w = Fukushima::LambertW(0,std::exp(lambert_log_arg));
        }
        fractions_c[c] = lambert_w*inv_total_times_variance;
        residual -= lambert_w;
    }
    if(std::fabs(step/log_weighted_total) > 1e-6){
        std::cerr << std::fixed << std::setprecision(12);
        std::cerr << "fitfrac: WARNING: dq = " << step << " is not close to zero at final q = " << log_weighted_total
            << " funcq = " << residual << std::endl;
    }
  if (!workspace_c) delete[] lambert_log_arg_base_c;
  return log_weighted_total;
}


// Residual of the consistency equation for q_g (SI Eq. 27, multiplied by n_g v_g): n_g v_g - sum_c W_c.
// It is zero at the solution and a decreasing function of q_g.
double fraction_sum_residual(double *lambert_log_arg_base_c,int num_cells,double total_times_variance, double log_weighted_total){

    double residual = total_times_variance;
    double lambert_w, lambert_log_arg;
    int c;

    for(c=0;c<num_cells;++c){
        lambert_log_arg = lambert_log_arg_base_c[c] - log_weighted_total;
        if(std::isnan(lambert_log_arg)){
            throw std::runtime_error("Error: x is NaN in fq function.");
        }
        if(lambert_log_arg > 50.0){
            lambert_w = lambert_w0_of_exp_asymptotic(lambert_log_arg);
        }
        else{
            lambert_w = Fukushima::LambertW(0,std::exp(lambert_log_arg));
        }
        residual -= lambert_w;
    }

    return residual;
}


// Newton step for q_g: the residual divided by its derivative with respect to q_g.
double newton_step(double *lambert_log_arg_base_c,int num_cells,double total_times_variance, double log_weighted_total){
    double numerator = total_times_variance;
    double denominator = 0;
    double lambert_log_arg, lambert_w;
    int c;

    for(c=0;c<num_cells;++c){
        lambert_log_arg = lambert_log_arg_base_c[c] - log_weighted_total;
        if(std::isnan(lambert_log_arg)){
            throw std::runtime_error("Error: x is NaN in deltaq function.");
        }
        if(lambert_log_arg > 50.0){
            lambert_w = lambert_w0_of_exp_asymptotic(lambert_log_arg);
        }
        else{
            lambert_w = Fukushima::LambertW(0,std::exp(lambert_log_arg));
        }
        numerator -= lambert_w;
        denominator += lambert_w/(1.0+lambert_w);
    }

    double ratio = numerator/denominator;
    return ratio;
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
void reconstruct_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, const std::vector<double>& cell_totals_c, double log_weighted_total, const double* log_cell_totals_c) {
    const double total_times_variance = gene_total * gene_variance;
    const double inv_total_times_variance = 1.0 / total_times_variance;
    const double log_total_times_variance = std::log(total_times_variance);
    for (int c = 0; c < num_cells; ++c) {
        const double lambert_log_arg_base = (log_cell_totals_c ? log_cell_totals_c[c] : std::log(cell_totals_c[c])) + umi_counts_c[c] * gene_variance + log_total_times_variance;
        const double lambert_log_arg = lambert_log_arg_base - log_weighted_total;
        const double lambert_w = lambert_log_arg > 50 ? lambert_w0_of_exp_asymptotic(lambert_log_arg) : Fukushima::LambertW(0, std::exp(lambert_log_arg));
        fractions_c[c] = lambert_w * inv_total_times_variance;
    }
}
