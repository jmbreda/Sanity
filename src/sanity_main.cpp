#include <cmath>
#include <ctime>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

#include "helpers/io/dense_row_reader.hpp"
#include "gene_fit.hpp"
#include "logging.hpp"
#include "sanity_input.hpp"
#include "sanity_options.hpp"
#include "sanity_output.hpp"
#include "Version.h"

int main(int argc, char **argv)
{
    SanityOptions options;
    const ParseResult parse_res = parse_options(argc, argv, options);
    if (parse_res == ParseResult::HelpRequested)
    {
        show_usage();
        return 0;
    }
    else if (parse_res == ParseResult::VersionRequested)
    {
        std::cout << "Sanity version " << SANITY_VERSION << "\n";
        return 0;
    }
    else if (parse_res == ParseResult::Error)
    {
        show_usage();
        return 1;
    }
    const int num_threads = options.num_threads;
    const int num_bins = options.num_bins;

    // counts will contain the index of the input: gene and cell counts and names, cell and gene totals, and where each gene's counts sit in the file.
    const CountIndex counts = load_count_index(options);
    const int num_genes = counts.num_genes; // G
    const int num_cells = counts.num_cells; // C

    // Cache library constants, preserving the sum order.
    std::vector<double> log_cell_totals_c(num_cells); // log(N_c).
    double sum_cell_totals = 0; // sum_c N_c, the total counts of the dataset.
    for (int c = 0; c < num_cells; ++c) {
        log_cell_totals_c[c] = std::log(counts.cell_totals_c[c]);
        sum_cell_totals += counts.cell_totals_c[c];
    }

    // Grid of candidate gene variances, computed once here and shared with
    // get_gene_expression_level. The grid is equally spaced in log(v_g), as in 
    // SI Section S1.1.4.
    double log_variance_step = std::log(options.variance_max / options.variance_min) / ((double)num_bins - 1);
    std::vector<double> variance_grid_b(num_bins); // v_b, the candidate values of v_g.
    for (int b = 0; b < num_bins; ++b)
        variance_grid_b[b] = options.variance_min * std::exp(log_variance_step * b);

    if (!prepare_output_folder(options))
        return 1;
    write_command_metadata(options, argc, argv);
    // Open every output file and write its header and the gene, cell and variance-grid identifiers.
    SanityOutput output(options, counts.gene_names_g, counts.cell_names_c, variance_grid_b);

    std::exception_ptr output_failure;
    logging_debug("Fit gene expression levels");
    const std::clock_t begin = std::clock();
    #pragma omp parallel num_threads(num_threads)
    {
        std::unique_ptr<sanity::DenseRowReader> thread_reader;
        std::vector<double> umi_counts_c; // n_gc, the counts of gene g in each cell (SI Eq. 13).
        GeneWorkspace workspace(num_cells, num_bins, options.variance_method); // Scratch space reused by every fit of this thread.
        #pragma omp for schedule(dynamic) ordered
        for (int g = 0; g < num_genes; ++g)
        {
            SanityOutput::PreparedRow prepared;
            std::exception_ptr local_failure;
            try
            {
                if (!thread_reader) {
                    thread_reader = std::make_unique<sanity::DenseRowReader>(options.in_file, options.in_file_extension,
                                                                             counts.mtx_rows_g, counts.tsv_offsets_g, num_cells);
                    umi_counts_c.resize(num_cells);
                }
                thread_reader->read_into(g, umi_counts_c.data());
                const RowComputation result = get_gene_expression_level(umi_counts_c, counts.gene_totals_g[g], variance_grid_b,
                                                                         options.variance_method, log_cell_totals_c, sum_cell_totals, workspace);
                prepared = output.prepare_row(result, counts.gene_names_g[g]);
            }
            catch (...) { local_failure = std::current_exception(); }
            // The `ordered` clause is what makes the writes below safe: it serialises this block
            // across threads and runs it in ascending g, so the output rows stay in gene order and
            // the shared output files are never written concurrently.
            #pragma omp ordered
            {
                // output esimated running time
                if (g == (3 * num_threads - 1))
                {
                    std::clock_t end = std::clock();
                    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    double estimated_running_time = (elapsed_secs / (3. * num_threads)) * num_genes / num_threads;
                    logging_debug("Estimated running time: " + format_duration(estimated_running_time) + " ");
                }

                if (g % 1000 == 0 && g > 0)
                {
                    logging_debug("Finished " + std::to_string(g) + " genes out of " + std::to_string(num_genes));
                }

                if (local_failure && !output_failure) output_failure = local_failure;
                if (!output_failure)
                {
                    try { output.write_row(prepared); }
                    catch (...) { output_failure = std::current_exception(); }
                }
            }
        }
    }

    if (output_failure) std::rethrow_exception(output_failure);
    output.close();
    logging_debug("Finished fitting all genes");

    return 0;
}
