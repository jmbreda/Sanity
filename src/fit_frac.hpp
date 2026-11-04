#ifndef SANITY_FIT_FRAC_HPP
#define SANITY_FIT_FRAC_HPP

#include <vector>

// Index letters at the end of a name give what an array runs over: c = cells.
// umi_counts_c = n_gc, gene_total = n_g, gene_variance = v_g, log_cell_totals_c = log(N_c), sum_cell_totals = sum_c N_c.

// Solves for q_g at one v_g and fills fractions_c with f_gc. previous_log_weighted_total is the starting guess (0 if none);
// workspace_c is scratch space of num_cells doubles.
double fit_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, double previous_log_weighted_total, const double* log_cell_totals_c, double sum_cell_totals, double* workspace_c);
// Fills fractions_c with f_gc at a known q_g (log_weighted_total), without repeating the root search.
void reconstruct_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, double log_weighted_total, const double* log_cell_totals_c);
double newton_step(double *lambert_log_arg_base_c, int num_cells, double total_times_variance, double log_weighted_total);
double fraction_sum_residual(double *lambert_log_arg_base_c, int num_cells, double total_times_variance, double log_weighted_total);
double lambert_w0_of_exp_asymptotic(const double& log_arg);
#endif
