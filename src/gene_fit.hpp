#ifndef SANITY_GENE_FIT_HPP
#define SANITY_GENE_FIT_HPP

#include <vector>

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, b = grid points of the gene variance v_g.
// The results of the fit of one gene.
struct RowComputation
{
    double mean_ltq;                                    // mu_g, the mean log transcription quotient (LTQ) of the gene (SI Eq. 45).
    double mean_ltq_error_sq;                           // (delta mu_g)^2, the variance of mu_g (SI Eq. 47).
    std::vector<double> log_fold_change_c;              // delta_gc, the log fold change of each cell from the mean LTQ (SI Eqs. 14 and 39).
    std::vector<double> log_fold_change_error_sq_c;     // epsilon_gc^2, the variance of delta_gc (SI Eq. 42).
    double gene_variance;                               // v_g, the variance of the log fold changes across cells (SI Eq. 17).
    std::vector<double> variance_likelihood_b;          // P(v_b | n_g), the posterior of v_g on the grid of variances (SI Eq. 34).
};

// Fits one gene: the log fold change of every cell, its error bar, the mean LTQ and the gene variance.
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, const std::vector<double> &cell_totals_c, double gene_total, const std::vector<double> &variance_grid_b, int num_cells, int num_bins, int variance_method, const std::vector<double>& log_cell_totals_c, double sum_cell_totals, std::vector<double>& fitting_workspace_c);

#endif
