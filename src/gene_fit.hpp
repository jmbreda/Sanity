#ifndef SANITY_GENE_FIT_HPP
#define SANITY_GENE_FIT_HPP

#include <vector>

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, b = grid points of the gene variance v_g.

// How the gene variance v_g is estimated. The values are the indices of the method names in the main program.
enum class VarianceMethod
{
    MARG = 0, // Average every estimate over the posterior of v_g.
    MLE = 1,  // Use the v_g that maximises the likelihood.
    MAP = 2,  // Use the v_g that maximises the posterior under a 1/v_g prior (the default).
    EAP = 3   // Use the grid point closest to the posterior mean of v_g.
};

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

// Scratch space of one thread, allocated once and reused for every gene the thread fits.
struct GeneWorkspace
{
    GeneWorkspace(int num_cells, int num_bins, VarianceMethod variance_method);

    std::vector<double> lambert_log_arg_base_c;         // Q_c, the per-cell constants of the fit at one v_g (SI Eq. 24), used by fit_fractions.
    std::vector<double> fractions_c;                    // f_gc, the fraction of the gene's counts expected in each cell (SI Eq. 23).
    std::vector<double> log_fold_change_c;              // delta_gc at one grid point.
    std::vector<double> inv_curv_diag_rank1_c;          // n_g f_gc^2 / (n_g f_gc + 1/v_g), one term of each sum over cells in SI Eq. 37.
    std::vector<double> inv_curv_diag_num_c;            // The numerator of SI Eq. 37, the sum with cell c removed.
    std::vector<double> inv_curv_diag_den2_c;           // n_g f_gc + 1/v_g, the diagonal of M_g before the rank-one correction (SI Eq. 31).
    std::vector<double> log_weighted_total_b;           // q_g(v_b) at every grid point (SI Eq. 22).
    std::vector<double> log_likelihood_b;               // L(v_b), the log marginal likelihood at every grid point (SI Eq. 32).
    std::vector<double> mean_ltq_b;                     // mu_g(v_b) = psi(n_g) - q_g(v_b) at every grid point (SI Eq. 44).
    std::vector<double> log_fold_change_bc;             // delta_gc at every grid point, row b at [b * num_cells]; MARG only.
    std::vector<double> log_fold_change_error_sq_bc;    // var(delta_gc) at every grid point (SI Eq. 37), same layout; MARG only.
};

// Fits one gene: the log fold change of every cell, its error bar, the mean LTQ and the gene variance.
// umi_counts_c = n_gc, the counts of this gene in each cell (SI Eq. 13); gene_total = n_g, their sum (SI Eq. 15).
// log_cell_totals_c = log(N_c), the log of the total counts of each cell (SI Eq. 9); sum_cell_totals = sum_c N_c.
// variance_grid_b = v_b, the candidate values of the gene variance (SI Section S1.1.4).
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, double gene_total, const std::vector<double> &variance_grid_b,
                                         VarianceMethod variance_method, const std::vector<double>& log_cell_totals_c,
                                         double sum_cell_totals, GeneWorkspace& workspace);

#endif
