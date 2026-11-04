#ifndef SANITY_INPUT_HPP
#define SANITY_INPUT_HPP

#include <ios>
#include <string>
#include <vector>

#include "helpers/io/read_input_files.hpp"
#include "sanity_options.hpp"

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// The result of one scan of the count matrix: its size, the totals and where each gene's row is in the file.
// Only genes with at least one count are kept; index g runs over those.
struct CountIndex
{
    int num_genes = 0;                          // G, the number of genes with at least one count.
    int num_cells = 0;                          // C, the number of cells.
    int num_input_rows = 0;                     // Number of rows (genes) in the input file, including rows without counts.
    std::vector<RowBlock> mtx_rows_g;           // Position and size of the entries of each kept gene in the MTX file.
    std::vector<std::streampos> tsv_offsets_g;  // Position of the row of each kept gene in the TSV file.
    std::vector<double> cell_totals_c;          // N_c, the total counts of each cell (SI Eq. 9).
    std::vector<double> gene_totals_g;          // n_g, the total counts of each gene (SI Eq. 15).
    std::vector<std::string> gene_names_g;
    std::vector<std::string> cell_names_c;
};

// Scans the input named in `options` once. MTX input without (readable) name files gets the
// default names Gene_1, Gene_2, ... and Cell_1, Cell_2, ...
CountIndex load_count_index(const SanityOptions &options);

#endif
