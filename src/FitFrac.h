#ifndef _FitFrac_h_
#define _FitFrac_h_

#include <vector>

double fit_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, const std::vector<double>& cell_totals_c, double previous_log_weighted_total, const double* log_cell_totals_c = nullptr, double sum_cell_totals = 0, double* workspace_c = nullptr);
void reconstruct_fractions(double *fractions_c, const std::vector<double>& umi_counts_c, double gene_total, double gene_variance, int num_cells, const std::vector<double>& cell_totals_c, double log_weighted_total, const double* log_cell_totals_c = nullptr);
double newton_step(double *lambert_log_arg_base_c, int num_cells, double total_times_variance, double log_weighted_total);
double fraction_sum_residual(double *lambert_log_arg_base_c, int num_cells, double total_times_variance, double log_weighted_total);
double lambert_w0_of_exp_asymptotic(const double& log_arg);
#endif
