#include <algorithm>
#include <cerrno>
#include <cmath>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/stat.h> // For mkdir on Unix-like systems
#include <vector>
#include <omp.h>

#ifdef _WIN32
    #include <direct.h> // For _mkdir on Windows
#endif

#include "helpers/io/read_input_files.hpp"
#include "helpers/io/dense_row_reader.hpp"
#include "gene_fit.hpp"
#include "logging.hpp"
#include "sanity_options.hpp"
#include "helpers/io/writer.hpp"
#include "helpers/io/format_output.hpp"
#include "helpers/io/npy_writer.hpp"
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
    const std::string &in_file = options.in_file;
    const std::string &in_file_extension = options.in_file_extension;
    const std::string &out_folder = options.out_folder;
    const int num_threads = options.num_threads;
    const bool print_extended_output = options.extended_output;
    const int num_bins = options.num_bins;
    const VarianceMethod variance_method = options.variance_method;
    const bool npy_output = options.npy_output;
    const bool bonsai_output_mode = options.bonsai_output_mode;
    std::string out_suffix = "";

    // count Number of genes and cells
    int num_genes; // G, the number of genes with at least one count.
    int num_cells; // C, the number of cells.
    int num_input_rows; // Number of rows (genes) in the input file, including rows without counts.
    std::map<int, int> row_to_gene; // Maps each MTX row to the index of its gene among the kept genes (-1 if the row has no counts).
    std::vector<RowBlock> mtx_rows_g; // Position and size of the entries of each kept gene in the MTX file.
    std::vector<std::streampos> tsv_offsets_g; // Position of the row of each kept gene in the TSV file.
    std::vector<double> cell_totals_c; // N_c, the total counts of each cell (SI Eq. 9).
    std::vector<double> gene_totals_g; // n_g, the total counts of each gene (SI Eq. 15).
    std::vector<std::string> gene_names_g;
    std::vector<std::string> cell_names_c;
    logging_debug("Reading input file to get gene and cell counts");
    if (in_file_extension == "mtx")
    {
        index_mtx_counts(in_file, num_input_rows, num_genes, num_cells, row_to_gene, mtx_rows_g, cell_totals_c, gene_totals_g, num_threads);
        gene_names_g = Read_GeneNames(options.gene_name_file, row_to_gene, num_genes);
        if (static_cast<int>(gene_names_g.size()) == 0)
        {
            logging_debug("Warning: gene names file not found or empty, generating default gene names.");
            gene_names_g.clear();
            for (int g = 0; g < num_genes; ++g)
            {
                gene_names_g.push_back("Gene_" + std::to_string(g + 1));
            }
        }
        cell_names_c = Read_CellNames(options.cell_name_file);
        if (static_cast<int>(cell_names_c.size()) == 0)
        {
            logging_debug("Warning: cell names file not found or empty, generating default cell names.");
            for (int c = 0; c < num_cells; ++c)
            {
                cell_names_c.push_back("Cell_" + std::to_string(c + 1));
            }
        }
        logging_debug("There were " + std::to_string(num_input_rows) + " rows");
    }
    else
    {
        index_tsv_counts(in_file, num_input_rows, num_genes, num_cells, tsv_offsets_g, cell_totals_c, gene_totals_g, cell_names_c, gene_names_g, num_threads);
        logging_debug("There were " + std::to_string(num_input_rows) + " rows");
    }
    std::string genes_message = "The input contained " + std::to_string(num_input_rows) + " rows for " + std::to_string(num_cells) +
                                 " cells, and " + std::to_string(num_genes) + " of these genes had at least one count in the dataset.";
    if (num_input_rows != num_genes)
    {
        genes_message += " " + std::to_string(num_input_rows - num_genes) + " genes were thus removed because no counts were detected.";
    }
    logging_debug(genes_message);

    // Cache library constants, preserving the sum order.
    std::vector<double> log_cell_totals_c(num_cells); // log(N_c).
    double sum_cell_totals = 0; // sum_c N_c, the total counts of the dataset.
    for (int c = 0; c < num_cells; ++c) {
        log_cell_totals_c[c] = std::log(cell_totals_c[c]);
        sum_cell_totals += cell_totals_c[c];
    }

    // Grid of candidate gene variances, computed once here and shared with
    // get_gene_expression_level rather than being recomputed there. Each value must stay exactly
    // variance_min * exp(log_variance_step * b), with log_variance_step formed exactly as below, to keep the output bit-identical.
    // The grid is equally spaced in log(v_g), as in SI Section S1.1.4.
    double log_variance_step = std::log(options.variance_max / options.variance_min) / ((double)num_bins - 1);
    std::vector<double> variance_grid_b(num_bins); // v_b, the candidate values of v_g.
    for (int b = 0; b < num_bins; ++b)
        variance_grid_b[b] = options.variance_min * std::exp(log_variance_step * b);

    // create output folder if it does not exist
    int mkdir_result = 0;
    #ifdef _WIN32
        mkdir_result = _mkdir(out_folder.c_str());
    #else
        mkdir_result = mkdir(out_folder.c_str(), 0755);
    #endif

    // Check for errors (ignore if directory already exists)
    if (mkdir_result != 0 && errno != EEXIST)
    {
        if (errno == ENOENT)
        {
            // mkdir only creates the final component, so a missing parent is the likely cause.
            logging_debug("Error: Failed to create output folder: " + out_folder +
                          "\nIts parent directory must already exist; create it first.");
        }
        else
        {
            logging_debug("Error: Failed to create output folder: " + out_folder);
        }
        return 1;
    }

    if (bonsai_output_mode)
    {
        for (const char* filename : {"log_transcription_quotients.txt", "ltq_error_bars.txt",
                                    "d_mu.txt", "likelihood.txt", "delta.txt", "d_delta.txt",
                                    "mu.txt", "variance.txt", "log_transcription_quotients.txt.gz",
                                    "ltq_error_bars.txt.gz", "d_mu.txt.gz", "likelihood.txt.gz",
                                    "delta.txt.gz", "d_delta.txt.gz", "mu.txt.gz", "variance.txt.gz"})
        {
            struct stat previous;
            if (stat((out_folder + filename).c_str(), &previous) == 0)
            {
                logging_debug("Bonsai output mode requires a folder without existing text-output files: " + out_folder);
                return 1;
            }
        }
    }

    // save version and command with parameters to a file
    std::ofstream cmd_file(out_folder + "sanity_command.txt", std::ios::trunc);
    std::time_t now = std::time(NULL);
    std::tm tm_snapshot;
    #ifdef _WIN32
        localtime_s(&tm_snapshot, &now);
    #else
        localtime_r(&now, &tm_snapshot);
    #endif

    char timestamp[20]; // "YYYY-MM-DD HH:MM:SS" + '\0' = 20
    std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &tm_snapshot);
    cmd_file << "# Timestamp: " << timestamp << "\n";

    cmd_file << "# Sanity version: " << SANITY_VERSION << "\n";

    cmd_file << "# Method: " << variance_method_name(variance_method) << "\n";
    if (bonsai_output_mode)
        cmd_file << "# Output mode: bonsai\n"
                    "# Output: Bonsai NPY; delta/d_delta: <f4; mu/variance: <f8; "
                    "C-order genes x cells; gene vectors: genes x 1; d_delta: standard deviation\n";
    cmd_file  << argv[0];
    for (int i = 1; i < argc; ++i)
    {
        cmd_file << " " << argv[i];
    }
    cmd_file << "\n";
    cmd_file.close();

    // -- Open files for writing --
    std::unique_ptr<npy::NpyStreamWriter> out_log_fold_change_npy;
    std::unique_ptr<npy::NpyStreamWriter> out_log_fold_change_error_npy;
    std::unique_ptr<npy::NpyStreamWriter> out_mean_ltq_npy;
    std::unique_ptr<npy::NpyStreamWriter> out_gene_variance_npy;

    sanity::Writer out_ltq, out_ltq_error, out_mean_ltq, out_mean_ltq_error, out_gene_variance,
        out_log_fold_change, out_log_fold_change_error, out_variance_likelihood, out_gene, out_cell;

    if (npy_output) {
        // The per-cell matrices are float32 in Bonsai output mode and float64 otherwise.
        npy::NpyDType cell_matrix_dtype = npy::NpyDType::Float64;
        if (bonsai_output_mode)
        {
            cell_matrix_dtype = npy::NpyDType::Float32;
        }
        out_log_fold_change_npy.reset(new npy::NpyStreamWriter(out_folder + "delta.npy", num_cells, cell_matrix_dtype));
        out_log_fold_change_error_npy.reset(new npy::NpyStreamWriter(out_folder + "d_delta.npy", num_cells, cell_matrix_dtype));
        out_mean_ltq_npy.reset(new npy::NpyStreamWriter(out_folder + "mu.npy", 1));
        out_gene_variance_npy.reset(new npy::NpyStreamWriter(out_folder + "variance.npy", 1));

        out_gene.open(out_folder + "geneID.txt");
        // save gene names
        for (int g = 0; g < num_genes; g++)
        {
            out_gene << gene_names_g[g].c_str() << "\n";
        }
        out_cell.open(out_folder + "cellID.txt");
        // save cell names
        for (int c = 0; c < num_cells; c++)
        {
            out_cell << cell_names_c[c].c_str() << "\n";
        }
    }
    else {
        if (options.gzip_output) {out_suffix = ".gz";}

        out_ltq.open(out_folder + "log_transcription_quotients.txt" + out_suffix);
        out_ltq << std::fixed << std::setprecision(6);

        out_ltq_error.open(out_folder + "ltq_error_bars.txt" + out_suffix);
        out_ltq_error << std::fixed << std::setprecision(6);

        out_ltq << "GeneID";
        out_ltq_error << "GeneID";
        for (int c = 0; c < num_cells; c++)
        {
            out_ltq << "\t" << cell_names_c[c].c_str();
            out_ltq_error << "\t" << cell_names_c[c].c_str();
        }
        out_ltq << "\n";
        out_ltq_error << "\n";
        if (print_extended_output)
        {

            out_gene.open(out_folder + "geneID.txt" + out_suffix);
            out_cell.open(out_folder + "cellID.txt" + out_suffix);
            out_mean_ltq.open(out_folder + "mu.txt" + out_suffix);
            out_mean_ltq << std::fixed << std::setprecision(6);

            out_mean_ltq_error.open(out_folder + "d_mu.txt" + out_suffix);
            out_mean_ltq_error << std::fixed << std::setprecision(6);

            out_gene_variance.open(out_folder + "variance.txt" + out_suffix);
            out_gene_variance << std::fixed << std::setprecision(6);

            out_log_fold_change.open(out_folder + "delta.txt" + out_suffix);
            out_log_fold_change << std::fixed << std::setprecision(6);

            out_log_fold_change_error.open(out_folder + "d_delta.txt" + out_suffix);
            out_log_fold_change_error << std::fixed << std::setprecision(6);

            out_variance_likelihood.open(out_folder + "likelihood.txt" + out_suffix);
            out_variance_likelihood << std::fixed << std::setprecision(6);

            out_variance_likelihood << "Variance";
            for (int k = 0; k < (num_bins); ++k)
            {
                out_variance_likelihood << "\t" << variance_grid_b[k];
            }
            out_variance_likelihood << "\n";

            // save cell names
            for (int c = 0; c < num_cells; c++)
            {
                out_cell << cell_names_c[c].c_str() << "\n";
            }
        }
    }

    std::exception_ptr output_failure;
    logging_debug("Fit gene expression levels");
    const std::clock_t begin = std::clock();
    #pragma omp parallel num_threads(num_threads)
    {
        std::unique_ptr<sanity::DenseRowReader> thread_reader;
        std::vector<double> umi_counts_c; // n_gc, the counts of gene g in each cell (SI Eq. 13).
        GeneWorkspace workspace(num_cells, num_bins, variance_method); // Scratch space reused by every fit of this thread.
        std::vector<double> npy_log_fold_change_error_c; // epsilon_gc (standard deviation) of the current gene, written to d_delta.npy.
        #pragma omp for schedule(dynamic) ordered
        for (int g = 0; g < num_genes; ++g)
        {
            RowComputation result;
            sanity::FormattedOutput formatted;
            std::vector<float> binary_log_fold_change_c, binary_log_fold_change_error_c; // float32 delta_gc and epsilon_gc for the Bonsai output.
            std::exception_ptr local_failure;
            try
            {
                if (!thread_reader) {
                    thread_reader = std::make_unique<sanity::DenseRowReader>(in_file, in_file_extension, mtx_rows_g, tsv_offsets_g, num_cells);
                    umi_counts_c.resize(num_cells);
                }
                thread_reader->read_into(g, umi_counts_c.data());
                result = get_gene_expression_level(umi_counts_c, gene_totals_g[g], variance_grid_b, variance_method, log_cell_totals_c, sum_cell_totals, workspace);
                if (!npy_output)
                    formatted = sanity::format_output_row(result, gene_names_g[g], print_extended_output);
                else
                {
                    if (!std::isfinite(result.mean_ltq) || !std::isfinite(result.gene_variance))
                        throw std::runtime_error("Nonfinite NPY gene output");
                    if (bonsai_output_mode)
                    {
                        binary_log_fold_change_c.resize(num_cells);
                        binary_log_fold_change_error_c.resize(num_cells);
                    }
                    else npy_log_fold_change_error_c.resize(num_cells);
                    for (int c = 0; c < num_cells; ++c)
                    {
                        const double log_fold_change_error = std::sqrt(result.log_fold_change_error_sq_c[c]);
                        if (!std::isfinite(result.log_fold_change_c[c]) || !std::isfinite(log_fold_change_error))
                            throw std::runtime_error("Nonfinite NPY cell output");
                        if (bonsai_output_mode)
                        {
                            binary_log_fold_change_c[c] = static_cast<float>(result.log_fold_change_c[c]);
                            binary_log_fold_change_error_c[c] = static_cast<float>(log_fold_change_error);
                            if (!std::isfinite(binary_log_fold_change_c[c]) || !std::isfinite(binary_log_fold_change_error_c[c]))
                                throw std::runtime_error("NPY output exceeds float32 range");
                        }
                        else npy_log_fold_change_error_c[c] = log_fold_change_error;
                    }
                }
            }
            catch (...) { local_failure = std::current_exception(); }
            // The `ordered` clause is what makes the writes below safe: it serialises this block
            // across threads and runs it in ascending g, so the output rows stay in gene order and
            // the shared output streams are never written concurrently.
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
                    try
                    {
                        // write output
                        if (npy_output) {
                            if (bonsai_output_mode)
                            {
                                out_log_fold_change_npy->write_row(binary_log_fold_change_c);
                                out_log_fold_change_error_npy->write_row(binary_log_fold_change_error_c);
                            }
                            else
                            {
                                out_log_fold_change_npy->write_row(result.log_fold_change_c);
                                out_log_fold_change_error_npy->write_row(npy_log_fold_change_error_c);
                            }
                            out_mean_ltq_npy->write_row(&result.mean_ltq, 1);
                            out_gene_variance_npy->write_row(&result.gene_variance, 1);
                        }
                        else
                        {
                            out_ltq.write_raw(formatted.ltq);
                            out_ltq_error.write_raw(formatted.ltq_error);
                            if (print_extended_output)
                            {
                                out_log_fold_change.write_raw(formatted.log_fold_change);
                                out_log_fold_change_error.write_raw(formatted.log_fold_change_error);
                                out_mean_ltq.write_raw(formatted.mean_ltq);
                                out_gene_variance.write_raw(formatted.gene_variance);
                                out_mean_ltq_error.write_raw(formatted.mean_ltq_error);
                                out_variance_likelihood.write_raw(formatted.variance_likelihood);
                                out_gene.write_raw(gene_names_g[g] + "\n");
                            }
                        }
                    }
                    catch (...) { output_failure = std::current_exception(); }
                }
            }
        }
    }

    if (output_failure) std::rethrow_exception(output_failure);
    for (auto* writer : {out_log_fold_change_npy.get(), out_log_fold_change_error_npy.get(), out_mean_ltq_npy.get(), out_gene_variance_npy.get()})
        if (writer) writer->close();
    for (auto* writer : {&out_ltq, &out_ltq_error, &out_mean_ltq, &out_mean_ltq_error, &out_gene_variance,
                         &out_log_fold_change, &out_log_fold_change_error, &out_variance_likelihood, &out_gene, &out_cell})
        writer->close();
    logging_debug("Finished fitting all genes");

    return 0;
}
