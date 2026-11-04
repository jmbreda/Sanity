#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <memory>
#include <exception>
#include <cerrno>
#include <cstdlib>
#include <sys/stat.h> // For mkdir on Unix-like systems
#include <map>
#include <cmath>
#include <stdexcept>

#ifdef _WIN32
    #include <direct.h> // For _mkdir on Windows
#endif

#include "helpers/io/read_input_files.hpp"
#include "helpers/io/dense_row_reader.hpp"
#include "gene_fit.hpp"
#include "helpers/io/writer.hpp"
#include "helpers/io/format_output.hpp"
#include "helpers/io/npy_writer.hpp"
#include "Version.h"

enum ParseResult
{
    CONTINUE,
    VERSION_REQUESTED,
    HELP_REQUESTED,
    ERROR
};


/***Function declarations ****/
ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &num_threads, bool &print_extended_output, double &variance_min, double &variance_max, int &num_bins, bool &no_norm, VarianceMethod &variance_method, bool &gzip_output, bool &npy_output, bool &bonsai_output_mode);
static void show_usage(void);

// Function to print a timestamped message to stderr
void logging_debug(const std::string &msg)
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::tm buf;
    #ifdef _WIN32
        localtime_s(&buf, &in_time_t);
    #else
        localtime_r(&in_time_t, &buf);
    #endif
    std::cerr << "[" << std::put_time(&buf, "%Y-%m-%d %H:%M:%S");
    std::cerr << ',' << std::setw(3) << std::setfill('0') << ms.count() << "] ";
    std::cerr << msg << std::endl;
}

int main(int argc, char **argv)
{
    std::string in_file("");
    std::string gene_name_file("none");
    std::string cell_name_file("none");
    std::string in_file_extension("");
    std::string out_folder("./");
    int num_threads(4);
    bool print_extended_output(false);
    double variance_min = 0.001; // v_min, the smallest candidate gene variance (SI Section S1.1.4).
    double variance_max = 50.0; // v_max, the largest candidate gene variance.
    int num_bins = 160; // B, the number of grid points of the gene variance between v_min and v_max.
    bool no_norm(false);
    VarianceMethod variance_method = VarianceMethod::MAP; // How to estimate v_g (see gene_fit.hpp).
    bool gzip_output = false;
    bool npy_output = false;
    bool bonsai_output_mode = false;
    std::string out_suffix = "";

    ParseResult parse_res = parse_argv(argc, argv, in_file, gene_name_file, cell_name_file, in_file_extension, out_folder, num_threads, print_extended_output, variance_min, variance_max, num_bins, no_norm, variance_method, gzip_output, npy_output, bonsai_output_mode);
    if (parse_res == HELP_REQUESTED)
    {
        show_usage();
        return 0;
    }
    else if (parse_res == VERSION_REQUESTED)
    {
        std::cout << "Sanity version " << SANITY_VERSION << "\n";
        return 0;
    }
    else if (parse_res == ERROR)
    {
        show_usage();
        return 1;
    }
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
        gene_names_g = Read_GeneNames(gene_name_file, row_to_gene, num_genes);
        if (static_cast<int>(gene_names_g.size()) == 0)
        {
            logging_debug("Warning: gene names file not found or empty, generating default gene names.");
            gene_names_g.clear();
            for (int g = 0; g < num_genes; ++g)
            {
                gene_names_g.push_back("Gene_" + std::to_string(g + 1));
            }
        }
        cell_names_c = Read_CellNames(cell_name_file);
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

    // Remove the total UMI correction if no cell size normalization option is true
    if (no_norm)
    {
        logging_debug("No cell size normalization will be performed");

        // get mean count per cell
        double mean_cell_total = 0;
        for (int c = 0; c < num_cells; ++c)
        {
            mean_cell_total += cell_totals_c[c];
        }
        mean_cell_total /= num_cells;

        // Now replace every N_c by the mean cell total
        for (int c = 0; c < num_cells; ++c)
        {
            cell_totals_c[c] = mean_cell_total;
        }
    }

    // Cache library constants after optional normalisation, preserving the sum order.
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
    double log_variance_step = std::log(variance_max / variance_min) / ((double)num_bins - 1);
    std::vector<double> variance_grid_b(num_bins); // v_b, the candidate values of v_g.
    for (int b = 0; b < num_bins; ++b)
        variance_grid_b[b] = variance_min * std::exp(log_variance_step * b);

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

    static const char* variance_method_names[] = {"MARG", "MLE", "MAP", "EAP"};
    cmd_file << "# Method: " << variance_method_names[static_cast<int>(variance_method)] << "\n";
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
        if (gzip_output) {out_suffix = ".gz";}

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

                    std::string estimated_running_time_str = "Estimated running time: ";
                    if (estimated_running_time > 86400)
                    {
                        std::stringstream stream;
                        stream << std::fixed << std::setprecision(2) << estimated_running_time / (3600. * 24.);
                        estimated_running_time_str += stream.str() + " days ";
                    }
                    else if (estimated_running_time > 3600)
                    {
                        std::stringstream stream;
                        stream << std::fixed << std::setprecision(2) << estimated_running_time / 3600.;
                        estimated_running_time_str += stream.str() + " hours ";
                    }
                    else if (estimated_running_time > 60)
                    {
                        std::stringstream stream;
                        stream << std::fixed << std::setprecision(2) << estimated_running_time / 60.;
                        estimated_running_time_str += stream.str() + " minutes ";
                    }
                    else
                    {
                        std::stringstream stream;
                        stream << std::fixed << std::setprecision(2) << estimated_running_time;
                        estimated_running_time_str += stream.str() + " seconds ";
                    }
                    logging_debug(estimated_running_time_str);
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

ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &num_threads, bool &print_extended_output, double &variance_min, double &variance_max, int &num_bins, bool &no_norm, VarianceMethod &variance_method, bool &gzip_output, bool &npy_output, bool &bonsai_output_mode)
{

    // Running with no arguments at all is a usage error, not a help request: main prints the usage
    // for both, but ERROR exits 1 whereas an explicit -h/--help still exits 0.
    if (argc < 2)
    {
        logging_debug("Error in argument parsing :\n"
                      "Not enough arguments provided.\n");
        return ERROR;
    }
    int i;

    std::string get_help[2];
    get_help[0] = "-h";
    get_help[1] = "--help";

    for (i = 1; i < argc; i++)
    {
        if (argv[i] == get_help[0] || argv[i] == get_help[1])
            return HELP_REQUESTED;
    }
    std::string get_version[2];
    get_version[0] = "-v";
    get_version[1] = "--version";
    for (i = 1; i < argc; i++)
    {
        if (argv[i] == get_version[0] || argv[i] == get_version[1])
        {
            return VERSION_REQUESTED;
        }
    }

    std::string extended_output("false");
    std::string no_norm_str("false");
    std::string variance_method_str("MAP");
    std::string to_find[][2] = {{"-f", "--file"},
                             {"-d", "--destination"},
                             {"-n", "--n_threads"},
                             {"-e", "--extended_output"},
                             {"-vmin", "--variance_min"},
                             {"-vmax", "--variance_max"},
                             {"-nbin", "--number_of_variance_bins"},
                             {"-mtx_genes", "--mtx_gene_name_file"},
                             {"-mtx_cells", "--mtx_cell_name_file"},
                             {"-no_norm", "--no_cell_size_normalization"},
                             {"-v_m", "--v_method"},
                             {"-gz", "--gzip_output"},
                             {"-npy", "--npy_output"},
                             {"--bonsai-output-mode", "--bonsai-output-mode"}};
    // Derived from the table itself so the count cannot drift from the array.
    const int N_param = sizeof(to_find) / sizeof(to_find[0]);

    // Returns the to_find index matched by `arg`, or -1 if `arg` is not a known option. Both the
    // matching loop below and the validation pass after it go through this, so the set of accepted
    // spellings cannot drift between them.
    auto match_option = [&](const std::string &arg) -> int
    {
        for (int t = 0; t < N_param; t++)
        {
            // The two development-only flags also keep their original spellings working.
            const bool legacy_spelling =
                (t == 11 && (arg == "--gz" || arg == "--gzip-output")) ||
                (t == 12 && (arg == "--npy" || arg == "--npy-output"));
            if (arg == to_find[t][0] || arg == to_find[t][1] || legacy_spelling)
                return t;
        }
        return -1;
    };
    // Format flags and Bonsai output mode take no argument; other known options take
    // exactly one. -h/--help and -v/--version are handled earlier and never reach here.
    auto takes_no_argument = [](int option_index) { return option_index == 11 || option_index == 12 || option_index == 13; };

    int j;
    int idx;
    for (j = 0; j < N_param; j++)
    {
        idx = 0;
        for (i = 1; i < argc; i++)
        {
            if (match_option(argv[i]) == j)
            {
                // Both flags are experimental and deliberately undocumented: they appear neither
                // in README.md nor in show_usage(). Neither takes an argument.
                //   -npy silently overrides -e and -gz (see the block at the end of this function)
                //   and writes a different SET of outputs, not merely a different format:
                //   delta.npy, d_delta.npy, mu.npy, variance.npy plus geneID.txt and cellID.txt,
                //   and NOT log_transcription_quotients, ltq_error_bars, likelihood or d_mu.
                //   d_delta.npy holds standard deviations, not variances.
                //   -gz output cannot be read back by Sanity_distance or Sanity_gene_correlation,
                //   which open delta.txt with a plain fopen and have no zlib path.
                if (j == 11)
                {
                    gzip_output = true;
                    continue; // no argument expected for -gz
                }
                if (j == 12)
                {
                    npy_output = true;
                    continue; // no argument expected for -npy
                }
                if (j == 13)
                {
                    bonsai_output_mode = true;
                    continue;
                }
                idx = i;
                if (idx + 1 > argc - 1)
                {
                    logging_debug("Error in argument parsing :\n" + std::string(argv[i]) + " option missing\n");
                    return ERROR;
                }
                if (j == 0)
                    in_file = argv[idx + 1];
                if (j == 1)
                    out_folder = argv[idx + 1];
                if (j == 2)
                    num_threads = std::atoi(argv[idx + 1]);
                if (j == 3)
                    extended_output = argv[idx + 1];
                if (j == 4)
                    variance_min = std::atof(argv[idx + 1]);
                if (j == 5)
                    variance_max = std::atof(argv[idx + 1]);
                if (j == 6)
                    num_bins = std::atoi(argv[idx + 1]);
                if (j == 7)
                    gene_name_file = argv[idx + 1];
                if (j == 8)
                    cell_name_file = argv[idx + 1];
                if (j == 9)
                    no_norm_str = argv[idx + 1];
                if (j == 10)
                    variance_method_str = argv[idx + 1];
                // add '/' to out_folder if not already. `-d ""` gives an empty string, on which
                // back() would be undefined behaviour, so treat it as the current directory.
                if (j == 1)
                {
                    if (out_folder.empty())
                        out_folder = "./";
                    else if (out_folder.back() != '/')
                        out_folder = out_folder + '/';
                }
            }
        }
        if (idx == 0 && j == 0)
        {
            logging_debug("Error in argument parsing :\n"
                          "missing input file name\n");
            return ERROR;
        }
    }

    // Reject anything on the command line that was not consumed above, either as a known option
    // or as the value of one. Without this, a typo such as "-e1" or "--bogus 7" is silently
    // ignored and the run quietly produces something other than what was asked for.
    std::vector<bool> consumed(argc, false);
    for (i = 1; i < argc; i++)
    {
        const int option_index = match_option(argv[i]);
        if (option_index < 0)
            continue;
        consumed[i] = true;
        if (!takes_no_argument(option_index) && i + 1 < argc)
        {
            // Consume the next token as this option's value, whatever it looks like: a value may
            // legitimately begin with '-', and no option here takes a negative number anyway.
            consumed[i + 1] = true;
            i++;
        }
    }
    for (i = 1; i < argc; i++)
    {
        if (!consumed[i])
        {
            logging_debug("Error in argument parsing :\nunrecognised option '" + std::string(argv[i]) + "'\n");
            return ERROR;
        }
    }

    if (extended_output == "true" || extended_output == "1")
        print_extended_output = true;

    if (no_norm_str == "true" || no_norm_str == "1")
        no_norm = true;

    if(variance_method_str == "MLE" || variance_method_str == "mle" || variance_method_str == "MaxLikelihood" || variance_method_str == "maxlikelihood" || variance_method_str == "max_likelihood"){
        variance_method = VarianceMethod::MLE;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the likelihood (MLE).");
    }
    else if(variance_method_str == "MAP" || variance_method_str == "map" || variance_method_str == "MaxAPosterior" || variance_method_str == "maxaposterior" || variance_method_str == "max_a_posterior"){
        variance_method = VarianceMethod::MAP;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the posterior (MAP).");
    }
    else if(variance_method_str == "EAP" || variance_method_str == "eap" || variance_method_str == "ExpectedAPosterior" || variance_method_str == "expectedaposterior" || variance_method_str == "expected_a_posterior"){
        variance_method = VarianceMethod::EAP;
        logging_debug("Outputting results for the expected value of the prior variance (v_g) over the posterior (EAP).");
    }
    else if(variance_method_str == "MARG" || variance_method_str == "marg" || variance_method_str == "marginalizing"){
        variance_method = VarianceMethod::MARG;
        logging_debug("Outputting results for the prior variance (v_g) marginalizing over v.");
    }
    else {
        logging_debug("Invalid method specified for variance estimation.");
        return ERROR;
    }

    // -npy silently overrides both -e and -gz rather than rejecting the combination: the npy path
    // writes its own fixed set of files (delta/d_delta/mu/variance .npy plus geneID.txt and
    // cellID.txt), so there is no extended output to print and no text stream to compress.
    // Like -gz this flag is experimental and intentionally absent from README.md and show_usage().
    if (bonsai_output_mode)
    {
        if (gzip_output || npy_output || variance_method != VarianceMethod::MAP)
        {
            logging_debug("Bonsai output mode writes NPY and requires MAP; omit gzip/NPY flags and use -v_m MAP.");
            return ERROR;
        }
        npy_output = true;
    }
    if(npy_output){
        gzip_output = false;
        print_extended_output = false;
    }

    // Get input file extension
    if (in_file.size() >= 3)
    {
        in_file_extension = in_file.substr(in_file.size() - 3);
    }
    else
    {
        in_file_extension = in_file;
    }

    // Check for .mtx.gz or .tsv.gz extensions
    if (in_file.size() >= 7)
    {
        std::string last_7 = in_file.substr(in_file.size() - 7);
        if (last_7 == ".mtx.gz")
        {
            in_file_extension = "mtx";
        }
        else if (last_7 == ".tsv.gz")
        {
            in_file_extension = "tsv";
        }
    }

    // Handle uncompressed extensions
    if (in_file_extension != "mtx" && in_file_extension != "tsv")
    {
        for (char &ch : in_file_extension)
        {
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
        if (in_file_extension != "mtx" && in_file_extension != "tsv")
        {
            in_file_extension = "tsv"; // default to tsv
        }
    }

    logging_debug("File type : " + in_file_extension + "\n");
    return CONTINUE;
}

static void show_usage(void)
{
    std::cerr << "Usage: Sanity <option(s)> SOURCES\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-v,--version\t\tShow the current version\n"
              << "\t-f,--file\t\tSpecify the input transcript count file: tab-separated values (.tsv) or Matrix Market File Format (.mtx)\n"
              << "\t-mtx_genes,--mtx_gene_name_file\tSpecify the gene name text file (only needed if .mtx input file)\n"
              << "\t-mtx_cells,--mtx_cell_name_file\tSpecity the cell name text file (only needed if .mtx input file)\n"
              << "\t-d,--destination\tSpecify the destination path (default: pwd)\n"
              << "\t-n,--n_threads\t\tSpecify the number of threads to be used (default: 4)\n"
              << "\t--bonsai-output-mode\tWrite only the four numerical files and IDs/metadata required by Bonsai (NPY, MAP)\n"
              << "\t-e,--extended_output\tOption to print extended output (default: false, choice: false,0,true,1)\n"
              << "\t-vmin,--variance_min\tMinimal value of variance in log transcription quotient (default: 0.001)\n"
              << "\t-vmax,--variance_max\tMaximal value of variance in log transcription quotient (default: 50)\n"
              << "\t-nbin,--number_of_variance_bins\tNumber of bins for the variance in log transcription quotient  (default: 160)\n"
              << "\t-no_norm,--no_cell_size_normalization\tOption to skip cell size normalization (default: false, choice: false,0,true,1)\n"
              << "\t-v_m,--v_method\t\tOption to specify the method for variance estimation (default: MAP, choice: MAP, EAP, MLE, MARG)\n";
}
