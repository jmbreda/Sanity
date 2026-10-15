# Proposal: clearer variable names in Sanity

Status: implemented (see [Implementation notes](#implementation-notes) at the end for where the code differs from this proposal).

Changes since draft 1:

- `mu` becomes `mean_ltq` and `delta` becomes `log_fold_change` everywhere.
- Every array name now ends in index letters that say what it loops over and whether it is a vector or a matrix (`log_fold_change_c`, `log_fold_change_gc`, `mean_ltq_g`).
- The table of misleading names now lists the proposed replacement for each entry.

## Why

Much of the numerical code uses one- and two-letter names (`f`, `q`, `Q`, `L`, `C`, `G`, `n`, `N_c`, `beta`, `ldet`, `lik`) that follow the notation in the Supplementary Information of Breda, Zavolan & van Nimwegen (2021).
A reader without the SI open next to the code cannot tell what they hold.
Several names are also reused for different quantities in different places (for example `beta`, `q`, `f`, `var_delta`), which is confusing even for someone who knows the maths.

This document lists every name we consider problematic, groups them by file, and proposes a replacement.
Each table keeps the paper symbol, so the maths stays traceable after the rename.

## Scope

In scope:

- `src/calc_true_variation_parallel_prior_mu_sigma.cpp` (the `Sanity` binary)
- `src/FitFrac.cpp`, `src/FitFrac.h`
- `src/Digamma_Trigamma.cpp`, `src/Digamma_Trigamma.h`
- `src/ReadInputFiles.cpp`, `src/ReadInputFiles.h`
- `src/FormatOutput.hpp`
- `src/compute_distance.cpp` (the `Sanity_distance` binary)
- `src/compute_gene_correlations.cpp` (the `Sanity_gene_correlations` binary)

Out of scope:

- Command-line flags, output file names and output formats.
These stay exactly as they are; every rename below is internal.
- `reproducibility/`: frozen scripts for the paper figures.
- `src/FukushimaLambertW.*` and `reproducibility/scripts/bhtsne/`: vendored third-party code.
- `src/DenseRowReader.*`, `src/FastMtxReader.hpp`, `src/Writer.hpp`, `src/npy_writer.*`: newer code that already uses descriptive names (the two small consistency fixes are listed under [I/O](#io-readinputfilescpp-readinputfilesh-fastmtxreaderhpp-tests)).

## Where the paper symbols come from

The "Paper symbol" column was built from the code comments and the prose of the bioRxiv preprint ([10.1101/2019.12.28.889956](https://www.biorxiv.org/content/10.1101/2019.12.28.889956v1.full)).
The preprint names the following symbols in its text: $n_{gc}$, $N_c$, $n_g$, $\alpha_{gc}$, $\mu_g$, $\delta_{gc}$, $v_g$, $f_{gc}$ (normalised quantities summing to 1 over cells), $q_g$, $W$ (Lambert W), $L^*(v_g)$, $M_g$ (curvature matrix), $\epsilon_{gc}$, $\delta\mu_g$, $\psi$, $\psi_1$.
Its equations are images, and the published SI numbers its equations differently, so the equation numbers quoted in the code comments (19, 56, 67, 85) have not been checked against the published SI.
Symbols marked *(code)* are not named in the preprint text and were derived from what the code computes.
Please correct any symbol that does not match the published SI.

## Naming conventions this proposal follows

### Words

1. **The two central quantities are spelled out.**
$\mu_g$ is `mean_ltq` (mean log transcription quotient of gene $g$) and $\delta_{gc}$ is `log_fold_change` (deviation of cell $c$ from that mean).
The output file names (`mu.txt`, `delta.txt`, …) do not change.
2. **Counts of things use a `num_` prefix**: `num_cells`, `num_genes`, `num_bins`.
In the paper, `n` and `N` already mean UMI counts ($n_{gc}$, $N_c$, $n_g$).
The current code has `N_c` (UMI total of cell $c$) next to `N_cell` (how many cells there are), which read almost identically but mean different things.
3. **UMI quantities**: `umi_counts` for $n_{gc}$, `cell_totals` for $N_c$ and `gene_totals` (or `gene_total` for a single gene) for $n_g$.
`cell_totals` and `gene_totals` are already used inside `ReadInputFiles.cpp`, `FastMtxReader.hpp` and the test drivers.
4. **The gene variance $v_g$ is always `gene_variance`**, and the grid of candidate values is `variance_grid`.
The bare word "variance" is not used for anything else.
5. **Uncertainties use `_error` (standard deviation) and `_error_sq` (variance)**, for example `log_fold_change_error` for $\epsilon_{gc}$ and `log_fold_change_error_sq` for $\epsilon_{gc}^2$.
This matches the existing output names (`ltq_error`, `d_delta.txt`, `d_mu.txt`).
We avoid `var_delta` because $v_g$ *is* the prior variance of $\delta$, so "variance of delta" is ambiguous.
6. **Functions are named after what they return.**

### Index suffix (new in draft 2)

Every array name ends in letters that say what its entries loop over.
The number of letters says whether it is a vector (one letter) or a matrix (two letters).
Scalars have no suffix.

| Letter | Loops over | Typical size |
|---|---|---|
| `g` | genes | `num_genes` |
| `c` | cells | `num_cells` |
| `b` | grid points of the gene variance $v_g$ ("bins") | `num_bins` |
| `a` | grid points of the scale $\alpha$ (only in `Sanity_distance`) | `num_alpha_bins` |

Rules:

- For a matrix the letters are in storage order: the first letter is the first index in `x[i][j]`.
- A suffix such as `_i_c` means "for the fixed gene `i`, a vector over cells" (used in `Sanity_gene_correlations`, which compares two genes `i` and `j`).
- The same letters are used as loop variables: `g` for genes, `c` for cells, `b` for grid points (currently `i`, `k` and `kk`).
- Inside the per-gene fit the gene is fixed, so a quantity that is a matrix over genes and cells elsewhere is a vector over cells there.

Examples:

| Name | Shape |
|---|---|
| `gene_variance` | scalar (one gene) |
| `mean_ltq_g` | vector over genes |
| `log_fold_change_c` | vector over cells (one gene) |
| `log_fold_change_gc` | matrix, genes × cells |
| `log_fold_change_cg` | matrix, cells × genes |
| `log_fold_change_bc` | matrix, variance-grid points × cells (one gene) |
| `distance_cc` | matrix, cells × cells |

The suffix also makes a storage difference between the two analysis tools visible.
`Sanity_distance` stores $\delta$ as `[cell][gene]` (`log_fold_change_cg`), while `Sanity_gene_correlations` stores it as `[gene][cell]` (`log_fold_change_gc`).

### Before and after

The inner loop of the variance grid (lines 590–599 of `calc_true_variation_parallel_prior_mu_sigma.cpp`), as it stands:

```cpp
delsq = 0;
L = -0.5 * ((double)C) * std::log(v);
for (i = 0; i < C; ++i)
{
    delta_v[k][i] = std::log(f[i]) - log_N_c[i] + q;
    L += n_c[i] * delta_v[k][i];
    delsq += delta_v[k][i] * delta_v[k][i];
}
L -= delsq / (2 * v);
L -= n * q;
```

and with the proposed names:

```cpp
sum_log_fold_change_sq = 0;
log_likelihood = -0.5 * num_cells * std::log(bin_variance);
for (c = 0; c < num_cells; ++c)
{
    log_fold_change_bc[b][c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total;
    log_likelihood += umi_counts_c[c] * log_fold_change_bc[b][c];
    sum_log_fold_change_sq += log_fold_change_bc[b][c] * log_fold_change_bc[b][c];
}
log_likelihood -= sum_log_fold_change_sq / (2 * bin_variance);
log_likelihood -= gene_total * log_weighted_total;
```

## Priority levels

- **High**: the name can only be understood with the SI at hand, or it is misleading.
- **Medium**: the name is guessable but unclear, or inconsistent with the rest of the code.
- **Low**: cosmetic or consistency-only.

---

## Names that are misleading, not only short

These are the most important to fix, because they actively suggest the wrong thing.

| Where | Current | Problem | Proposed |
|---|---|---|---|
| `FitFrac.cpp:12` | `beta` | Means $n_g v_g$ here, but $1/(n_g v_g)$ in `get_gene_expression_level`. | `total_times_variance` |
| `calc_true_variation…cpp:584` | `beta` | The reciprocal of the `FitFrac.cpp` one. | `inv_total_times_variance` |
| `calc_true_variation…cpp:534` | `sig2_delta_c` | Not a variance of $\delta$. It holds $n f_c^2/(n f_c + 1/v)$, one term in the closed-form (Sherman–Morrison) diagonal of the *inverse* curvature matrix $M_g^{-1}$. | `inv_curv_diag_rank1_c` |
| `calc_true_variation…cpp:535–537` | `sig2_delta_num`, `sig2_delta_den1`, `sig2_delta_den2` | Same family: the numerator and the two denominator factors of that diagonal, not variances. The final ratio *is* $\epsilon_{gc}^2$. | `inv_curv_diag_num_c`, `inv_curv_diag_den1`, `inv_curv_diag_den2_c` |
| `calc_true_variation…cpp:401`, `RowComputation` | `var_delta` | Posterior variance of $\delta$, but in the plain NPY branch it is overwritten with the *standard deviation*, so after that line the name is wrong. | `log_fold_change_error_sq_c`; the NPY branch writes into a separate `log_fold_change_error_c` |
| `calc_true_variation…cpp:523` | `lik` | First holds log-likelihoods, then (from line 633) normalised likelihoods, in the same array. | two arrays: `log_likelihood_b` and `variance_likelihood_b` |
| `FitFrac.cpp:140` | `LambertW0_approximation(x)` | Returns $W_0(e^x)$, not $W_0(x)$. The argument is the *log* of the Lambert-W argument. | `lambert_w0_of_exp_asymptotic(log_arg)` |
| `FitFrac.cpp:90` | `fq` | Named after "$f(q)$" in the comments, but `f` is also the fitted fractions in the same file. Returns $nv - \sum_c W_c$. | `fraction_sum_residual` |
| `FitFrac.cpp:115` | `deltaq` | Not a difference of anything; it is the Newton step residual / derivative. | `newton_step` |
| `compute_distance.cpp:252` | `q` | Unrelated to $q_g$ in the inference: here the total variance of a cell-pair difference, $\epsilon^2_{\text{pair}} + \alpha v_g$. | `pair_total_variance` |
| `compute_distance.cpp:252` | `f` | Unrelated to $f_{gc}$: here the shrinkage factor $\alpha v_g/(\alpha v_g + \epsilon^2_{\text{pair}})$. | `shrinkage` |
| `compute_distance.cpp:74` | `var_delta` | The sample variance of $\delta$ across cells, not the posterior variance that `var_delta` means in `Sanity`. | `log_fold_change_sample_var` |
| `compute_gene_correlations.cpp:238` | `vi`, `vj` | Sample variance of $\delta$ in `get_correlation_simple`, but the gene variance $v_g$ in `get_correlation_full`. | `sample_var_i`, `sample_var_j` in the first; `gene_variance_i`, `gene_variance_j` in the second |
| `compute_gene_correlations.cpp:289` | `cov` | A per-cell product $z_i z_j$, not a covariance. | `z_product_c` |
| `compute_gene_correlations.cpp:270`, `:358` | `denom` | Holds the reciprocal $1/(\phi_c - r^2)$, not the denominator. | `inv_denominator` |
| `compute_distance.cpp:202`, `compute_gene_correlations.cpp:182` | `factor` | Rescales `delta` and `epsilon2` in place, after which they no longer hold the posterior values. | `prior_removal_factor`, see the rescaling item below |
| `compute_distance.cpp:200–213`, `compute_gene_correlations.cpp:181–195` | `delta`, `epsilon2` (after rescaling) | The arrays are overwritten with the rescaled values, so the same name holds two different things. | Rescale in place (no extra memory), but call the routine `remove_prior_shrinkage` and name its parameters `unshrunk_log_fold_change_cg` and `unshrunk_log_fold_change_error_sq_cg` (`_gc` in the correlation tool) |

---

## Inference core: `calc_true_variation_parallel_prior_mu_sigma.cpp`

### `RowComputation` (lines 44–52) and `FormatOutput.hpp`

`RowComputation` describes one gene, so there is no gene index in the names.

| Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|
| `mu` | `mean_ltq` | Mean LTQ of the gene | $\mu_g$ | High |
| `var_mu` | `mean_ltq_error_sq` | Posterior variance of $\mu_g$ | $\delta\mu_g^2$ | Medium |
| `delta` | `log_fold_change_c` | Log fold change of each cell from the mean LTQ | $\delta_{gc}$ | High |
| `var_delta` | `log_fold_change_error_sq_c` | Posterior variance of each $\delta_{gc}$ | $\epsilon_{gc}^2$ | Medium |
| `var_gene` | `gene_variance` | Variance of the LTQ across cells | $v_g$ | Medium |
| `lik` | `variance_likelihood_b` | Normalised likelihood of each grid value of $v_g$ (written to `likelihood.txt`) | $P(n_g \mid v_g)$ | High |

`FormattedOutput` (strings holding one formatted output line) follows:

| Current | Proposed |
|---|---|
| `ltq`, `ltq_error` | keep |
| `delta` | `log_fold_change` |
| `delta_error` | `log_fold_change_error` |
| `mu` | `mean_ltq` |
| `mu_error` | `mean_ltq_error` |
| `variance` | `gene_variance` |
| `likelihood` | `variance_likelihood` |

### Problem size, counts and the variance grid (`main`)

| Line | Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|---|
| 114 | `G` | `num_genes` | Number of genes with at least one count | – | High |
| 114 | `C` | `num_cells` | Number of cells | – | High |
| 116 | `N_rows` | `num_input_rows` | Number of genes in the input, including all-zero rows | – | Medium |
| 118 | `gene_idx` | `row_to_gene` | Map from MTX row to kept-gene index (−1 if dropped) | – | Medium |
| 119 | `mtx_rows` | `mtx_rows_g` | Index block of each kept gene | – | Low |
| 120 | `tsv_offsets` | `tsv_offsets_g` | File offset of each kept gene's row | – | Low |
| 121 | `N_c` | `cell_totals_c` | Total UMI count of each cell | $N_c$ | High |
| 122 | `n` | `gene_totals_g` | Total UMI count of each gene | $n_g$ | High |
| 123–124 | `gene_names`, `cell_names` | `gene_names_g`, `cell_names_c` | | – | Low |
| 169 | `mean_N_c` | `mean_cell_total` | Mean of $N_c$ (used with `-no_norm`) | – | Medium |
| 184 | `log_N_c` | `log_cell_totals_c` | $\log N_c$, cached | $\log N_c$ | Medium |
| 185 | `library_total` | `sum_cell_totals` | $\sum_c N_c$ | – | Medium |
| 87–88 | `vmin`, `vmax` | `variance_min`, `variance_max` | Range of the variance grid (matches `--variance_min/max`) | $v_{\min}$, $v_{\max}$ | Low |
| 89 | `numbin` | `num_bins` | Number of variance-grid points | $N_b$ | Medium |
| 194 | `deltav` | `log_variance_step` | Spacing of the grid in $\log v$ | $\Delta v$ | Medium |
| 195 | `v_grid` | `variance_grid_b` | Candidate values of $v_g$ | $v_b$ | Medium |
| 196 | `kk` | `b` | Grid index | $b$ | Low |
| 91 | `v_method` | `variance_method` | Which estimate of $v_g$ to report (0 MARG, 1 MLE, 2 MAP, 3 EAP) | – | Medium |
| 256 | `vm_name` | `variance_method_names` | Labels for `variance_method` | – | Low |
| 85 | `N_threads` | `num_threads` | Thread count | – | Low |
| 362 | `n_c_g` | `umi_counts_c` | Counts of the current gene in every cell | $n_{gc}$ | High |
| 362 | `fitting_workspace` | `fitting_workspace_c` | Reusable scratch vector | – | Low |
| 368 | `binary_delta`, `binary_error` | `binary_log_fold_change_c`, `binary_log_fold_change_error_c` | float32 rows for Bonsai output | – | Low |

### Output streams (`main`, lines 271–277)

| Current | Proposed | File | Priority |
|---|---|---|---|
| `out_exp_lev` | `out_ltq` | `log_transcription_quotients.txt` | Medium |
| `out_d_exp_lev` | `out_ltq_error` | `ltq_error_bars.txt` | Medium |
| `out_mu` | `out_mean_ltq` | `mu.txt` | Low |
| `out_dmu` | `out_mean_ltq_error` | `d_mu.txt` | Low |
| `out_delta` | `out_log_fold_change` | `delta.txt` | Low |
| `out_ddelta` | `out_log_fold_change_error` | `d_delta.txt` | Low |
| `out_var_gene` | `out_gene_variance` | `variance.txt` | Low |
| `out_lik` | `out_variance_likelihood` | `likelihood.txt` | Low |
| `out_delta_npy`, `out_ddelta_npy`, `out_mu_npy`, `out_var_gene_npy` | `out_log_fold_change_npy`, `out_log_fold_change_error_npy`, `out_mean_ltq_npy`, `out_gene_variance_npy` | `delta.npy`, `d_delta.npy`, `mu.npy`, `variance.npy` | Low |

### `get_gene_expression_level` (lines 504–758)

Function name: keep, or `fit_gene` (Low).

Parameters:

| Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|
| `n_c` | `umi_counts_c` | Counts of this gene in each cell | $n_{gc}$ | High |
| `N_c` | `cell_totals_c` | Total UMI count of each cell | $N_c$ | High |
| `n` | `gene_total` | Total UMI count of this gene | $n_g$ | High |
| `v_grid` | `variance_grid_b` | Candidate values of $v_g$ | $v_b$ | Medium |
| `C` | `num_cells` | | | High |
| `numbin` | `num_bins` | | $N_b$ | Medium |
| `v_method` | `variance_method` | | | Medium |
| `log_N_c` | `log_cell_totals_c` | | $\log N_c$ | Medium |
| `library_total` | `sum_cell_totals` | | $\sum_c N_c$ | Medium |
| `fitting_workspace` | `fitting_workspace_c` | | | Low |

Locals:

| Line | Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|---|
| 514 | `i` | `c` | Cell index | $c$ | Low |
| 514 | `k` | `b` | Grid index | $b$ | Low |
| 515 | `beta` | `inv_total_times_variance` | $1/(n_g v_b)$; the opposite of `beta` in `FitFrac.cpp` | – *(code)* | High |
| 515 | `L` | `log_likelihood` | Log marginal likelihood of $v_b$ (Laplace approximation) | $L^*(v_g)$ | High |
| 515 | `ldet` | `log_det_curvature` | $\log\det(M_g/n_g)$; the constant $C\log n_g$ is dropped because it does not depend on $v$ | $\log\det M_g$ | High |
| 515 | `q` | `log_weighted_total` | $\log\sum_c N_c e^{\delta_c}$, fixed so that the fractions sum to 1 | $q_g$ | High |
| 515 | `delsq` | `sum_log_fold_change_sq` | $\sum_c \delta_c^2$ | – | Medium |
| 516 | `prev_q` | `previous_log_weighted_total` | $q_g$ from the previous bin, used as a starting guess | – | Medium |
| 517 | `f` | `fractions_c` | Fitted fraction of the gene's counts expected in each cell, $N_c e^{\delta_c}/\sum_{c'} N_{c'} e^{\delta_{c'}}$ (the lambda at line 544 already calls it `fractions`) | $f_{gc}$ | High |
| 518 | `delta_v` | `log_fold_change_bc` | $\delta_{gc}$ at every grid value (MARG) | – | Medium |
| 519 | `delta_scratch` | `log_fold_change_buffer_c` | Single reusable array for MLE/MAP/EAP | – | Low |
| 520 | `sig2_delta_v` | `log_fold_change_error_sq_bc` | $\epsilon_{gc}^2$ at every grid value (MARG) | – | High |
| 522 | `fitted_q` | `log_weighted_total_b` | $q_g$ at every grid value | – | Medium |
| 523 | `lik` | `log_likelihood_b`, then `variance_likelihood_b` | See [misleading names](#names-that-are-misleading-not-only-short) | – | High |
| 524 | `delta`, `var_delta` | `log_fold_change_c`, `log_fold_change_error_sq_c` | Reported values; copied into `RowComputation` | $\delta_{gc}$, $\epsilon_{gc}^2$ | Medium |
| 525 | `var_mu` | `mean_ltq_error_sq` | Posterior variance of $\mu_g$ | $\delta\mu_g^2$ | Medium |
| 525 | `var_gene` | `gene_variance` | Reported $v_g$ | $v_g$ | Medium |
| 539 | `mu_v` | `mean_ltq_b` | $\mu_g$ at every grid value | – | Medium |
| 678 | `mu` | `mean_ltq` | Reported $\mu_g$ | $\mu_g$ | High |
| 540 | `Lmax` | `max_log_likelihood` | | – | Medium |
| 541 | `Lmax_ind` | `max_likelihood_bin` | MLE bin | – | Medium |
| 542, 583 | `v` | `bin_variance` | $v_b$ in the current bin | $v_b$ | Medium |
| 580 | `psi0` | `digamma_total` | $\psi(n_g)$ | $\psi(n_g)$ | Medium |
| 580 | `psi1` | `trigamma_total` | $\psi_1(n_g)$, the posterior variance of $\mu_g$ | $\psi_1(n_g)$ | Medium |
| 629 | `sum_L` | `likelihood_sum` | Normalising constant over the grid | – | Medium |
| 642 | `vindex` | `selected_bin` | Bin chosen by MLE/MAP/EAP | – | Medium |
| 649 | `mapmax` | `max_posterior` | | – | Low |
| 651 | `curv` | `bin_variance` | Reads as "curvature"; it is $v_b$ | $v_b$ | Medium |
| 661 | `postmean` | `posterior_mean_variance` | EAP value of $v_g$ | – | Low |
| 666 | `mindist` | `min_distance` | | – | Low |

Uncertainty lambda (lines 544–578):

| Line | Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|---|
| 544 | `calculate_uncertainty` | `compute_log_fold_change_error_sq` | Diagonal of $M_g^{-1}$, with the zero-count correction | $\epsilon_{gc}^2$ | Medium |
| 544 | `fractions` | `fractions_c` | | $f_{gc}$ | Medium |
| 544 | `bin_delta` | `log_fold_change_c` | | $\delta_{gc}$ | Medium |
| 545 | `variance` | `bin_variance` | | $v_b$ | Low |
| 545 | `output` | `log_fold_change_error_sq_c` | | $\epsilon_{gc}^2$ | Medium |
| 547 | `inv_v` | `inv_bin_variance` | $1/v_b$ | – | Low |
| 534 | `sig2_delta_c` | `inv_curv_diag_rank1_c` | $n f_c^2/(n f_c + 1/v)$, the rank-one part | – *(code)* | High |
| 536 | `sig2_delta_den2` | `inv_curv_diag_den2_c` | $n f_c + 1/v$, the diagonal of $M_g$ before the rank-one correction | – *(code)* | High |
| 537 | `sig2_delta_den1` | `inv_curv_diag_den1` | $1 - \sum_c$ `inv_curv_diag_rank1_c` | – *(code)* | High |
| 535 | `sig2_delta_num` | `inv_curv_diag_num_c` | `inv_curv_diag_den1` + `inv_curv_diag_rank1_c`; the ratio `num / (den1 * den2)` is $\epsilon_{gc}^2$ | – *(code)* | High |

### `get_epsilon_2` (lines 760–786)

| Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|
| `get_epsilon_2` | `zero_count_log_fold_change_error_sq` | $\epsilon_{gc}^2$ for a cell with $n_{gc}=0$, from a log-likelihood drop of ½ | $\epsilon_{gc}^2$ | High |
| `d` | `log_fold_change` | $\delta_{gc}$ | $\delta_{gc}$ | High |
| `v` | `gene_variance` | | $v_g$ | Medium |
| `n` | `gene_total` | | $n_g$ | High |
| `f` | `fraction` | | $f_{gc}$ | High |
| `e` | `displacement` | Upward shift of $\delta$ being tested | – | Medium |
| `drop` | `log_likelihood_drop` | | – | Low |
| `lo`, `hi` | keep | Bisection bracket | – | – |

### `parse_argv`

`j == 11`, `j == 12` and `t == 11` refer to rows of the `to_find` table by position.
This is a magic-number problem rather than a naming one; see [related suggestions](#related-suggestions-not-renames).

---

## Fraction fitting: `FitFrac.cpp`, `FitFrac.h`

Background for reviewers: for each candidate $v$, the fractions satisfy $W_c = n v f_c$ with $W_c = W_0\!\left(e^{Q_c - q}\right)$ and $Q_c = \log N_c + n_{gc} v + \log(n v)$, and $q$ is found by Newton's method so that $\sum_c W_c = n v$ (equivalently $\sum_c f_c = 1$).

| Line | Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|---|
| 8 | `fitfrac` | `fit_fractions` | Solves for the fractions; returns $q_g$ | – | Medium |
| 155 | `fractions_at_q` | `reconstruct_fractions` | Recomputes fractions from a saved $q_g$ without the root search | – | Medium |
| 90 | `fq` | `fraction_sum_residual` | $n v - \sum_c W_c$; zero at the solution | – | High |
| 115 | `deltaq` | `newton_step` | Ratio residual / derivative | – | High |
| 140 | `LambertW0_approximation` | `lambert_w0_of_exp_asymptotic` | Asymptotic series for $W_0(e^x)$ at large $x$ | – | High |
| 8 | `f` | `fractions_c` | | $f_{gc}$ | High |
| 8 | `n_c` | `umi_counts_c` | | $n_{gc}$ | High |
| 8 | `n` | `gene_total` | | $n_g$ | High |
| 8 | `v` | `gene_variance` | | $v_g$ | Medium |
| 8 | `C` | `num_cells` | | | High |
| 8 | `N_c` | `cell_totals_c` | | $N_c$ | High |
| 8 | `prev_q` | `previous_log_weighted_total` | | – | Medium |
| 8 | `log_N_c` | `log_cell_totals_c` | | | Medium |
| 8 | `library_total` | `sum_cell_totals` | | | Medium |
| 8 | `workspace` | `workspace_c` | | | Low |
| 9 | `q` | `log_weighted_total` | | $q_g$ | High |
| 9 | `W` | `lambert_w` | $W_c = n v f_c$ | $W$ | Medium |
| 9 | `x` | `lambert_log_arg` | $Q_c - q$ | – *(code)* | Medium |
| 9 | `dq` | `step` | Newton step for $q$ | – | Medium |
| 9 | `funcq` | `residual` | | – | Medium |
| 11 | `Q` | `lambert_log_arg_base_c` | $Q_c = \log N_c + n_{gc} v + \log(n v)$ | – *(code)* | High |
| 12 | `beta` | `total_times_variance` | $n_g v$; the reciprocal of `beta` in `get_gene_expression_level` | – *(code)* | High |
| 13 | `inv_beta` | `inv_total_times_variance` | $1/(n_g v)$ | – | High |
| 14 | `logbeta` | `log_total_times_variance` | | – | Medium |
| 57 | `tol` | `newton_tolerance` | | – | Low |
| 136 | `rat` | `step` | | – | Low |
| 142–145 | `L2`, `xsq`, `L2sq`, `L2cub` | keep `L2` (standard notation in the asymptotic series for Lambert W, $L_2 = \ln\ln z$) and add a comment; `xsq` → `x_sq` | – | – | Low |

## Digamma and trigamma: `Digamma_Trigamma.cpp`, `Digamma_Trigamma.h`

| Current | Proposed | Meaning | Paper symbol | Priority |
|---|---|---|---|---|
| `Psi_0` | `digamma` | $\psi(x)$ | $\psi$ | Medium |
| `Psi_1` | `trigamma` | $\psi_1(x)$ | $\psi_1$ | Medium |
| parameter `n` in the header, `x` in the source | `x` in both | | | Low |

---

## I/O: `ReadInputFiles.cpp`, `ReadInputFiles.h`, `FastMtxReader.hpp`, tests

| Where | Current | Proposed | Meaning | Priority |
|---|---|---|---|---|
| `.h:50`, `.cpp:215` | `Get_G_C_UMIcountMatrix` | `index_tsv_counts` | Scans a TSV, returns dimensions, totals and row offsets | Medium |
| `.h:62`, `.cpp:711` | `Get_G_C_MTX` | `index_mtx_counts` | Same for MTX | Medium |
| `.cpp:520` | `Get_G_C_MTX_legacy` | `index_mtx_counts_serial` | Single-threaded fallback | Low |
| `.h`, `.cpp` | `G`, `C`, `N_rows` | `num_genes`, `num_cells`, `num_input_rows` | | High |
| `.h`, `.cpp` | `N_c`, `n` | `cell_totals_c`, `gene_totals_g` | | High |
| `.h`, `.cpp` | `cell_names`, `gene_names` | `cell_names_c`, `gene_names_g` | | Low |
| `.h`, `.cpp` | `tsv_offsets`, `mtx_rows` | `tsv_offsets_g`, `mtx_rows_g` | | Low |
| `.h`, `.cpp` | `gene_idx` | `row_to_gene` | | Medium |
| `.cpp:227`, `:523` | `infp` | `reader` | It is a `FileReader`, not a `FILE*` | Low |
| `.cpp:238` | `ss` | `line` | | Medium |
| `.cpp:364` | `thread_ss` | `thread_line` | | Medium |
| `.cpp:365` | `thread_sc` | `thread_line_copy` | | Medium |
| `.cpp:275` | `N_fields_first_row` | `num_fields_first_row` | | Low |
| `.cpp:608` | `g_idx`, `c_idx` | `row_index`, `cell_index` | | Low |
| `.cpp:693` | `cur` | `next_gene` | | Low |
| `.cpp:348–351`, `:595–598` | `gene_totals_tmp`, `cell_totals`, `gene_totals` | add suffix: `gene_totals_tmp_g`, `cell_totals_c`, `gene_totals_g` (index over input rows where the zero-count genes are not yet dropped) | | Low |
| `FastMtxReader.hpp:96` | `libraries`, `totals` | `cell_totals_c`, `gene_totals_g` | Already the names used inside the same function | Low |
| `tests/mtx_index.cpp:22`, `tests/read_rows.cpp` | `libraries`, `totals` | `cell_totals_c`, `gene_totals_g` | | Low |

---

## Cell-cell distance: `compute_distance.cpp`

Background for reviewers: before computing distances, each $\delta_{gc}$ and $\epsilon_{gc}^2$ is multiplied by $v_g/(v_g - \epsilon_{gc}^2)$.
If the posterior is Gaussian, that factor undoes the shrinkage from the $\mathcal{N}(0, v_g)$ prior, so the rescaled values are the likelihood-only estimate of $\delta$ and its variance.
This is our reading of the code; please confirm.
The distance then averages over a scale $\alpha \in [0, 2]$ applied to $v_g$.

Matrices are stored `[cell][gene]` here, hence `_cg`.

| Line | Current | Proposed | Meaning | Priority |
|---|---|---|---|---|
| 36–37 | `N_gene`, `N_cell` | `num_genes`, `num_cells` | | Medium |
| 35 | `N_threads` | `num_threads` | | Low |
| 54 | `delta_all` | `log_fold_change_all_cg` | Before filtering | Medium |
| 58 | `epsilon2_all` | `log_fold_change_error_sq_all_cg` | Before filtering | Medium |
| 62 | `variance_all` | `gene_variance_all_g` | Before filtering | Medium |
| 69 | `idx_gene` | `keep_gene_g` | Gene passes the signal-to-noise filter | Medium |
| 70 | `N_gene_cutoff` | `num_genes_kept` | | Medium |
| 73 | `mean_delta` | `mean_log_fold_change` | Mean of $\delta$ across cells, one gene | Medium |
| 74 | `var_delta` | `log_fold_change_sample_var` | Sample variance of $\delta$ across cells | High |
| 75 | `mean_epsilon2` | `mean_log_fold_change_error_sq` | | Medium |
| 108 | `delta` | `log_fold_change_cg` | After filtering | Medium |
| 112 | `epsilon2` | `log_fold_change_error_sq_cg` | $\epsilon_{gc}^2$ | Medium |
| 116 | `variance` | `gene_variance_g` | $v_g$ | Medium |
| 155 | `D` | `distance_cc` | Cell-by-cell distance matrix | Medium |
| 202 | `factor` | `prior_removal_factor` | $v_g/(v_g - \epsilon_{gc}^2)$ | High |
| 217 | `N_bin` | `num_alpha_bins` | | Medium |
| 218 | `da` | `alpha_step` | | Medium |
| 219 | `a` | `alpha` | | Medium |
| 221 | `alpha_v` | `alpha_times_variance_ga` | $\alpha_a v_g$ for every gene and grid point | High |
| 233 | `N_parloop` | `num_row_pairs` | Rows are processed in pairs $(i, N-1-i)$ to balance load | Medium |
| 249 | `get_Di_errorbar` | `distance_row_errorbar` | Fills row $i$ of the distance matrix | Medium |
| 249 | `D` (parameter) | `distance_row_c` | Row $i$ of `distance_cc` | Medium |
| 331 | `get_Di_euclidean` | `distance_row_euclidean` | | Medium |
| 252 | `q` | `pair_total_variance` | $\epsilon^2_{\text{pair}} + \alpha v_g$ | High |
| 252 | `f` | `shrinkage` | $\alpha v_g / (\alpha v_g + \epsilon^2_{\text{pair}})$ | High |
| 253 | `x2` | `log_fold_change_diff_sq_g` | $(\delta_{gi} - \delta_{gj})^2$ | High |
| 254 | `eps2` | `pair_error_sq_g` | $\epsilon_{gi}^2 + \epsilon_{gj}^2$ | High |
| 255 | `lik` | `alpha_likelihood_a` | Log-likelihood, then normalised posterior, of each $\alpha$ | High |
| 256 | `lik_max`, `lik_tot` | `max_log_likelihood`, `likelihood_sum` | | Medium |
| 350 | `my_file` | `path` | | Low |
| 351 | `infp` | `file` | | Low |
| 352 | `ss` | `line` | | Low |
| 356 | `epsilon` | `log_fold_change_error` | $\epsilon_{gc}$ read from `d_delta.txt` | Medium |

## Gene-gene correlation: `compute_gene_correlations.cpp`

Matrices are stored `[gene][cell]` here, hence `_gc`.
Functions that compare two genes use `_i_c` and `_j_c` for the vectors over cells of gene `i` and gene `j`.

| Line | Current | Proposed | Meaning | Priority |
|---|---|---|---|---|
| 37–38 | `N_gene`, `N_cell` | `num_genes`, `num_cells` | | Medium |
| 48 | `delta_all` | `log_fold_change_all_gc` | Before filtering | Medium |
| 54 | `epsilon2_all` | `log_fold_change_error_sq_all_gc` | Before filtering | Medium |
| 58 | `variance_all` | `gene_variance_all_g` | Before filtering | Medium |
| 67 | `s2n_gene` | `signal_to_noise_g` | Per-gene $\sqrt{\langle \delta^2/\epsilon^2 \rangle}$ | Medium |
| 68 | `N_gene_cutoff` | `num_genes_kept` | | Medium |
| 99 | `delta` | `log_fold_change_gc` | After filtering | Medium |
| 103 | `epsilon2` | `log_fold_change_error_sq_gc` | | Medium |
| 107 | `variance` | `gene_variance_g` | | Medium |
| 108 | `orignum` | `original_gene_index_g` | Row in `delta.txt` before filtering | Medium |
| 157 | `Rs` | `corr_simple_gg` | Pearson correlation of $\delta$ | High |
| 158 | `Rfull` | `corr_full_gg` | Correlation accounting for error bars | High |
| 159 | `Rerr` | `corr_error_gg` | Error bar on `corr_full_gg` | High |
| 182 | `factor` | `prior_removal_factor` | As in `compute_distance.cpp` | High |
| 236 | `deli`, `delj` | `log_fold_change_i_c`, `log_fold_change_j_c` | $\delta$ of genes $i$ and $j$ over cells | High |
| 238 | `mui`, `muj` | `mean_log_fold_change_i`, `mean_log_fold_change_j` | Sample mean of $\delta$, not $\mu_g$ | High |
| 238 | `vi`, `vj` | `sample_var_i`, `sample_var_j` | Sample variance of $\delta$, not $v_g$ | High |
| 238 | `r` | `corr` | | Low |
| 267 | `deriv` | `log_likelihood_slope` | $\partial \log L / \partial r$ | High |
| 267 | `cov`, `phi`, `V` (parameters) | `z_product_c`, `variance_inflation_product_c`, `cross_term_c` | See below | High |
| 269 | `deriv` (local) | `slope` | Shadows the function name | Medium |
| 270, 358 | `denom` | `inv_denominator` | $1/(\phi_c - r^2)$ | High |
| 284 | `rfull`, `rerr`, `rs` | `corr_full`, `corr_error`, `corr_simple` | Scalars for one gene pair | High |
| 284 | `eps2i`, `eps2j` | `log_fold_change_error_sq_i_c`, `log_fold_change_error_sq_j_c` | | High |
| 284 | `vi`, `vj` | `gene_variance_i`, `gene_variance_j` | $v_g$ here | High |
| 287 | `phi` | `variance_inflation_product_c` | $\phi_{i,c}\,\phi_{j,c}$ | High |
| 288 | `V` | `cross_term_c` | $z_i^2\phi_j + z_j^2\phi_i$ | High |
| 289 | `cov` | `z_product_c` | $z_{i,c}\,z_{j,c}$ | High |
| 292–293 | `sigi`, `sigj` | `gene_sd_i`, `gene_sd_j` | $\sqrt{v_g}$ | Medium |
| 297–298 | `phi_i`, `phi_j` | `variance_inflation_i`, `variance_inflation_j` | $1 + \epsilon^2/v_g$ | High |
| 299–300 | `zi`, `zj` | `z_i`, `z_j` | $\delta/\sqrt{v_g}$ | Low |
| 308 | `ders` | `slope_at_simple` | Slope at `corr_simple` | Medium |
| 310 | `der` | `slope_at_bound` | Slope at $r = \pm 1$ | Medium |
| 310, 340 | `rmin`, `rmax`, `rmid` | `corr_low`, `corr_high`, `corr_mid` | Bisection bracket | Low |
| 339 | `rtol` | `corr_tolerance` | | Low |
| 340 | `curder` | `slope_at_mid` | | Medium |
| 353 | `rfinal` | `corr_estimate` | | Low |
| 357 | `secder` | `log_likelihood_curvature` | Second derivative at the optimum | Medium |

---

## Same name, different meaning

A summary of the collisions above, which is the main argument for renaming rather than only adding comments.

| Name | Meanings in the current code |
|---|---|
| `beta` | $n_g v$ (`FitFrac.cpp`); $1/(n_g v)$ (`calc_true_variation…cpp`) |
| `q` | $\log\sum_c N_c e^{\delta_c}$ (inference); total variance of a cell-pair difference (distance) |
| `f` | Fitted cell fraction (inference); shrinkage factor (distance); "$f(q)$" residual in comments (`FitFrac.cpp`) |
| `var_delta` | Posterior variance of $\delta_{gc}$ (`Sanity`); sample variance across cells (`Sanity_distance`); standard deviation after line 401 (NPY) |
| `vi`, `vj` | Sample variance of $\delta$; gene variance $v_g$ (both in `compute_gene_correlations.cpp`) |
| `lik` | Log-likelihood, then normalised likelihood (both files) |
| `delta`, `epsilon2` | Posterior values, then values with the prior shrinkage removed (distance and correlation tools) |
| `N_…` | UMI total ($N_c$) or number of items (`N_cell`, `N_gene`, `N_rows`, `N_threads`) |
| `variance` | $v_g$ (distance, correlation, `variance.txt`); bin value $v_b$ (lambda parameter) |

## Related suggestions (not renames)

These came up while reading and would make the same code easier to follow.
They are listed separately because they go beyond renaming.

1. Replace the integer `v_method` (0 MARG, 1 MLE, 2 MAP, 3 EAP) with an `enum class VarianceMethod`.
2. Replace the positional option indices in `parse_argv` (`j == 11`, `t == 12`, …) with named constants.
3. Split `lik` into two arrays: `log_likelihood_b` and `variance_likelihood_b`.
4. In the NPY branch, write standard deviations into a separate `log_fold_change_error_c` vector instead of overwriting the variance.
5. Rename `calc_true_variation_parallel_prior_mu_sigma.cpp` to something like `sanity_main.cpp` (needs a Makefile change).
6. `compute_gene_correlations.cpp:321`: the comment says 0.156 is a z-score, but the code uses 0.159, which is the one-sided tail probability beyond 1σ.

## How we would implement it

- Pure renames, no change in behaviour.
- Three separate PRs so each can be reviewed on its own: (1) `FitFrac`, `Digamma_Trigamma` and the inference core; (2) I/O; (3) `Sanity_distance` and `Sanity_gene_correlations`.
- Each PR is checked for bit-identical output against the current `master` binary with `tests/compare.py` and `tests/compare_inference.py` for all four `-v_m` methods, plus `compare_mtx.py`, `compare_mtx_index.py` and `compare_rows.py` for the I/O PR.
- Renames done with the compiler's help (rename one symbol at a time and rebuild), not with a global search and replace, because single letters such as `n`, `f`, `q` and `C` occur inside other identifiers.

## Open questions for reviewers

1. `cell_totals` / `gene_totals`, or `library_sizes` / `gene_totals`?
`cell_totals` is already used in the I/O code and tests; "library size" is the more common term in the field.
2. `_error_sq`, or `_posterior_var` (for example `log_fold_change_posterior_var`)?
3. Is `log_weighted_total` a good name for $q_g$, or does the SI give it a name we should reuse?
4. Are the names with the index suffix too long in the formulas (for example `log_fold_change_error_sq_all_cg`)?
If so, a short alias such as `lfc` could be used for `log_fold_change` throughout.
5. Does anyone call `fitfrac`, `Psi_0`/`Psi_1` or `Get_G_C_*` from outside this repository?
If so, we would keep the old names as thin wrappers for one release.

## Implementation notes

The renames were applied as proposed, with these differences.

- **Sherman–Morrison terms**: named `inv_curv_diag_rank1_c`, `inv_curv_diag_num_c`, `inv_curv_diag_den1` and `inv_curv_diag_den2_c`, because together they give the diagonal of the inverse curvature matrix $M_g^{-1}$ (SI Eqs. 35–37).
- **Local names that would have clashed with a renamed parameter** got an `input_`, `kept_` or `exact_` prefix in `ReadInputFiles.cpp` and `FastMtxReader.hpp` (for example `input_cell_totals_c`, `kept_gene_totals_g`).
- **Lambda parameters** of the uncertainty routine are `given_fractions_c`, `given_log_fold_change_c`, `given_variance` and `error_sq_out_c`, so they do not shadow the outer variables.
- **Rescaling in the distance and correlation tools** is still done in place and inline.
No separate `remove_prior_shrinkage` function was added.
In `Sanity_distance`, only the parameters of `distance_row_errorbar` carry the `unshrunk_` prefix.
- **Two small structural changes**, both without effect on the output: `lik` is split into `log_likelihood_b` and `variance_likelihood_b`, and the plain NPY branch writes standard deviations into its own vector `npy_log_fold_change_error_c` instead of overwriting the variance.
- **Not done** (they go beyond renaming): the `VarianceMethod` enum, named constants for the option indices in `parse_argv`, renaming `calc_true_variation_parallel_prior_mu_sigma.cpp`, and the comment fix at `compute_gene_correlations.cpp:321`.
- **Strings are untouched.**
Log and error messages that mention old function or variable names (for example `fitfrac: WARNING` and `x is NaN in fq function`) keep their text so that the program output does not change.

