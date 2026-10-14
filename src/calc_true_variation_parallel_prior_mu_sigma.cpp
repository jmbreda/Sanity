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

#include "ReadInputFiles.h"
#include "DenseRowReader.hpp"
#include "FitFrac.h"
#include "Digamma_Trigamma.h"
#include "Writer.hpp"
#include "FormatOutput.hpp"
#include "npy_writer.hpp"
#include "Version.h"

enum ParseResult
{
    CONTINUE,
    VERSION_REQUESTED,
    HELP_REQUESTED,
    ERROR
};

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

/***Function declarations ****/
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, const std::vector<double> &cell_totals_c, double gene_total, const std::vector<double> &variance_grid_b, int num_cells, int num_bins, int variance_method, const std::vector<double>& log_cell_totals_c, double sum_cell_totals, std::vector<double>& fitting_workspace_c);
double zero_count_log_fold_change_error_sq(double log_fold_change, double gene_variance, double gene_total, double fraction);
ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &num_threads, bool &print_extended_output, double &variance_min, double &variance_max, int &num_bins, bool &no_norm, int &variance_method, bool &gzip_output, bool &npy_output, bool &bonsai_output_mode);
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
    int variance_method = 2; // How to estimate v_g: 0 = MARG, 1 = MLE, 2 = MAP (the default), 3 = EAP.
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
    cmd_file << "# Method: " << variance_method_names[variance_method] << "\n";
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
        out_log_fold_change_npy.reset(new npy::NpyStreamWriter(out_folder + "delta.npy", num_cells,
            bonsai_output_mode ? npy::NpyDType::Float32 : npy::NpyDType::Float64));
        out_log_fold_change_error_npy.reset(new npy::NpyStreamWriter(out_folder + "d_delta.npy", num_cells,
            bonsai_output_mode ? npy::NpyDType::Float32 : npy::NpyDType::Float64));
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
        std::vector<double> fitting_workspace_c(num_cells); // Scratch space reused by every fit of this thread.
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
                result = get_gene_expression_level(umi_counts_c, cell_totals_c, gene_totals_g[g], variance_grid_b, num_cells, num_bins, variance_method, log_cell_totals_c, sum_cell_totals, fitting_workspace_c);
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

// "SI" is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, b = grid points of the gene variance v_g.
// Where the SI has n_g + 1 (Eqs. 16, 19, 26, 31, 44, 46) this code uses n_g.
// This is intentional: the prior over alpha_g differs from the one in the SI.
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, const std::vector<double> &cell_totals_c,
                                         double gene_total, const std::vector<double> &variance_grid_b, int num_cells, int num_bins,
                                         int variance_method, const std::vector<double>& log_cell_totals_c,
                                         double sum_cell_totals, std::vector<double>& fitting_workspace_c)
{
    // umi_counts_c = n_gc, the counts of this gene in each cell.
    // cell_totals_c = N_c, the total counts of each cell.
    // gene_total = n_g, the total counts of this gene.
    // variance_grid_b = v_b, the candidate values of the gene variance v_g (SI S1.1.4).
    // num_cells = C and num_bins = B.
    int c, b;
    double inv_total_times_variance; // 1/(n_g v_b); the reciprocal of total_times_variance in FitFrac.cpp.
    double log_likelihood; // L(v_b), the log marginal likelihood of the current grid point (SI Eq. 32).
    double log_det_curvature; // log det(M_g/n_g) (SI Eqs. 31, 33); the constant C*log(n_g) is dropped because it does not depend on v_g.
    double log_weighted_total; // q_g = log sum_c N_c exp(delta_gc) (SI Eq. 22).
    double sum_log_fold_change_sq; // sum_c delta_gc^2, from the Gaussian prior term of Eq. 20.
    double previous_log_weighted_total = 0.0; // q_g of the previous grid point, the starting guess for the next fit.
    double *fractions_c = new double[num_cells]; // f_gc, the fitted fraction of the gene's counts expected in each cell (SI Eq. 23).
    double **log_fold_change_bc = new double *[num_bins]; // delta*_gc(v_b) at every grid point (SI Eq. 28), kept for MARG.
    double *log_fold_change_buffer_c = variance_method > 0 ? new double[num_cells] : nullptr; // Reusable delta*_gc for the methods that keep one bin.
    double **log_fold_change_error_sq_bc = new double *[num_bins](); // var(delta_gc)(v_b) at every grid point (SI Eq. 37), kept for MARG.
    // Fixed-bin methods retain only the fitted scalar; reconstruct cell values after selection.
    std::vector<double> log_weighted_total_b(num_bins); // q_g(v_b) at every grid point.
    std::vector<double> log_likelihood_b(num_bins, -1.0); // L_b = log P(n_g | v_b) (SI S1.1.4).
    std::vector<double> log_fold_change_c(num_cells, 0.0); // The reported log fold change <delta_gc> (SI Eq. 39).
    std::vector<double> log_fold_change_error_sq_c(num_cells, 0.0); // The reported error bar squared, epsilon_gc^2 (SI Eq. 42).
    double mean_ltq_error_sq; // (delta mu_g)^2, the variance of the mean LTQ (SI Eq. 47).
    double gene_variance; // The reported v_g.

    for (b = 0; b < num_bins; ++b)
    {
        log_fold_change_bc[b] = variance_method == 0 ? new double[num_cells] : log_fold_change_buffer_c;
        if (variance_method == 0) log_fold_change_error_sq_bc[b] = new double[num_cells];
    }

    /*** To compute the variance of the log fold change ***/
    // Together these give the diagonal of the inverse of M_g in closed form (SI Eqs. 35-37):
    // var(delta_gc) = inv_curv_diag_num_c / (inv_curv_diag_den1 * inv_curv_diag_den2_c).
    double *inv_curv_diag_rank1_c = new double[num_cells]; // n_g f_gc^2 / (n_g f_gc + 1/v_g), one term of each sum over cells in Eq. 37.
    double *inv_curv_diag_num_c = new double[num_cells]; // The numerator of Eq. 37: den1 + rank1_c, i.e. the sum with cell c removed.
    double *inv_curv_diag_den2_c = new double[num_cells]; // n_g f_gc + 1/v_g, the diagonal of M_g before the rank-one correction (Eq. 31).
    double inv_curv_diag_den1; // 1 - sum_c rank1_c, the factor from the matrix determinant lemma (Eq. 33).

    double *mean_ltq_b = new double[num_bins]; // mu_g(v_b) = psi(n_g) - q_g(v_b) at every grid point (SI Eq. 44).
    double max_log_likelihood = -1e+100; // Maximum of L_b over the grid.
    int max_likelihood_bin = 0; // Grid point at which L_b is maximal (the MLE).
    double bin_variance; // v_b, the gene variance at the current grid point.

    auto compute_log_fold_change_error_sq = [&](const double* given_fractions_c, const double* given_log_fold_change_c,
                                                double given_variance, double* error_sq_out_c) {
        // rank1_c = n_g f_gc^2 / (n_g f_gc + 1/v_g) for each cell c.
        const double inv_variance = 1.0 / given_variance;
        for (c = 0; c < num_cells; ++c)
        {
            inv_curv_diag_rank1_c[c] = gene_total * given_fractions_c[c] * given_fractions_c[c] / (gene_total * given_fractions_c[c] + inv_variance);
        }
        // den1 is one minus the full sum of the rank-one terms; den2_c is the diagonal of M_g.
        inv_curv_diag_den1 = 1.0;
        for (c = 0; c < num_cells; ++c)
        {
            inv_curv_diag_den1 -= inv_curv_diag_rank1_c[c];
            inv_curv_diag_den2_c[c] = gene_total * given_fractions_c[c] + inv_variance;
        }
        // The numerator removes the term of cell c from the sum.
        for (c = 0; c < num_cells; c++)
        {
            inv_curv_diag_num_c[c] = inv_curv_diag_den1 + inv_curv_diag_rank1_c[c];
        }
        // The variance of the log fold change of each cell.
        for (c = 0; c < num_cells; ++c)
        {
            error_sq_out_c[c] = inv_curv_diag_num_c[c] / (inv_curv_diag_den1 * inv_curv_diag_den2_c[c]);
        }

        // fix computation of asymmetric variance for zero count (SI Eq. 38)
        for (c = 0; c < num_cells; ++c)
        {
            if (umi_counts_c[c] <= 0.5)
            {
                error_sq_out_c[c] = zero_count_log_fold_change_error_sq(given_log_fold_change_c[c], given_variance, gene_total, given_fractions_c[c]);
            }
        }
    };

    const double digamma_total = digamma(gene_total); // psi(n_g) (SI Eq. 44).
    const double trigamma_total = trigamma(gene_total); // psi_1(n_g), the variance of log(alpha_g) (SI Eq. 46).
    for (b = 0; b < num_bins; ++b)
    {
        bin_variance = variance_grid_b[b];
        inv_total_times_variance = 1.0 / (gene_total * bin_variance);
        log_weighted_total = fit_fractions(fractions_c, umi_counts_c, gene_total, bin_variance, num_cells, cell_totals_c, previous_log_weighted_total, log_cell_totals_c.data(), sum_cell_totals, fitting_workspace_c.data());
        previous_log_weighted_total = log_weighted_total;
        log_weighted_total_b[b] = log_weighted_total;
        mean_ltq_b[b] = digamma_total - log_weighted_total; /*** SI Eq. 44 ***/

        sum_log_fold_change_sq = 0;
        log_likelihood = -0.5 * ((double)num_cells) * std::log(bin_variance); // 1st term of SI Eq. 20
        for (c = 0; c < num_cells; ++c)
        {
            log_fold_change_bc[b][c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total; // SI Eq. 28
            log_likelihood += umi_counts_c[c] * log_fold_change_bc[b][c];                 // data term sum_c n_gc*delta_gc, SI Eq. 20
            sum_log_fold_change_sq += log_fold_change_bc[b][c] * log_fold_change_bc[b][c];
        }
        log_likelihood -= sum_log_fold_change_sq / (2 * bin_variance); // 2nd term of SI Eq. 20
        log_likelihood -= gene_total * log_weighted_total;     // last term of SI Eq. 20, using sum_c N_c exp(delta_gc) = exp(q_g)

        // get the determinant of the matrix (SI Eq. 33)
        log_det_curvature = 0.0;
        for (c = 0; c < num_cells; ++c)
        {
            log_det_curvature += (fractions_c[c] * fractions_c[c]) / (fractions_c[c] + inv_total_times_variance);
        }

        log_det_curvature = std::log(1 - log_det_curvature);
        for (c = 0; c < num_cells; ++c)
        {
            log_det_curvature += std::log(fractions_c[c] + inv_total_times_variance);
        }
        log_likelihood -= 0.5 * log_det_curvature; // SI Eq. 32
        log_likelihood_b[b] = log_likelihood;

        if (log_likelihood > max_log_likelihood)
        {
            max_log_likelihood = log_likelihood;
            max_likelihood_bin = b;
        }

        if (variance_method == 0)
        {
            compute_log_fold_change_error_sq(fractions_c, log_fold_change_bc[b], bin_variance, log_fold_change_error_sq_bc[b]);
        }
    } // end v bins loop

    // get normalized likelihood from loglikelihood (SI Eq. 34)
    std::vector<double> variance_likelihood_b(num_bins); // P(v_b | n_g), the posterior over the grid (SI Eq. 34); written to likelihood.txt.
    double likelihood_sum; // The normalisation, sum_b exp(L_b).
    likelihood_sum = 0.0;
    for (b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] = log_likelihood_b[b] - max_log_likelihood;
        variance_likelihood_b[b] = std::exp(variance_likelihood_b[b]);
        likelihood_sum += variance_likelihood_b[b];
    }
    for (b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] /= likelihood_sum;
    }

    int selected_bin = 0; // Grid point chosen by the MLE, MAP or EAP methods.
    if(variance_method ==1){
        selected_bin = max_likelihood_bin;
    }
    else if(variance_method == 2){
        // Dividing the normalised likelihood by v applies an implicit 1/v (Jeffreys) prior on the
        // gene variance. This is the only prior in effect anywhere in this computation.
        double max_posterior = -1e+100;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];

            if (variance_likelihood_b[b]/grid_variance > max_posterior)
            {
                max_posterior = variance_likelihood_b[b]/grid_variance;
                selected_bin = b;
            }
        }
    }
    else if(variance_method == 3){
        double posterior_mean_variance = 0.0;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];
            posterior_mean_variance += variance_likelihood_b[b] * grid_variance;
        }
        double min_distance = 1e+100;
        for (b = 0; b < num_bins; b++)    {
            double grid_variance = variance_grid_b[b];
            if (std::fabs(grid_variance - posterior_mean_variance) < min_distance)
            {
                min_distance = std::fabs(grid_variance - posterior_mean_variance);
                selected_bin = b;
            }
        }
    }

    // Average the mean LTQ over the variance grid (SI Eq. 45)
    double mean_ltq = 0.0; // mu_g, the mean LTQ log(alpha_g) of the gene (SI S1.1.7).
    for (b = 0; b < num_bins; b++)
    {
        mean_ltq += variance_likelihood_b[b] * mean_ltq_b[b];
    }

    // Compute (delta mu_g)^2 = psi_1 + < (mu - <mu>)^2 > (SI Eq. 47)
    mean_ltq_error_sq = trigamma_total;
    for (b = 0; b < num_bins; b++)
    {
        mean_ltq_error_sq += variance_likelihood_b[b] * (mean_ltq_b[b] - mean_ltq) * (mean_ltq_b[b] - mean_ltq);
    }

    if (variance_method == 0)  // MARG
    {
        // Compute <delta> = int p(v)*delta(v) dv (SI Eq. 39)
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_c[c] = 0.0;
            for (b = 0; b < num_bins; b++)
            {
                log_fold_change_c[c] += variance_likelihood_b[b] * log_fold_change_bc[b][c];
            }
        }
        // Compute epsilon^2 = < (delta - <delta>)^2 > + <var(delta)> (SI Eq. 42)
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_error_sq_c[c] = 0.0;
            for (b = 0; b < num_bins; b++)
            {
                log_fold_change_error_sq_c[c] += variance_likelihood_b[b] * (log_fold_change_bc[b][c] - log_fold_change_c[c]) * (log_fold_change_bc[b][c] - log_fold_change_c[c]) + variance_likelihood_b[b] * log_fold_change_error_sq_bc[b][c];
            }
        }
        // The reported v_g is the posterior mean of the variance.
        gene_variance = 0.0;
        for (b = 0; b < num_bins; ++b)
        {
            bin_variance = variance_grid_b[b];
            gene_variance += bin_variance * variance_likelihood_b[b];
        }
    }

    if (variance_method > 0)
    {
        // Store the gene-variance that maximizes the likelihood:
        gene_variance = variance_grid_b[selected_bin];
        // And then also the corresponding values for the LTQs etc.
        mean_ltq = mean_ltq_b[selected_bin];
        mean_ltq_error_sq = trigamma_total;
        reconstruct_fractions(fractions_c, umi_counts_c, gene_total, gene_variance, num_cells, cell_totals_c, log_weighted_total_b[selected_bin], log_cell_totals_c.data());
        for (c = 0; c < num_cells; c++)
        {
            log_fold_change_c[c] = std::log(fractions_c[c]) - log_cell_totals_c[c] + log_weighted_total_b[selected_bin];
        }
        compute_log_fold_change_error_sq(fractions_c, log_fold_change_c.data(), gene_variance, log_fold_change_error_sq_c.data());
    }

    delete[] fractions_c;
    for (b = 0; b < num_bins; ++b)
    {
        if (variance_method == 0) delete[] log_fold_change_bc[b];
        delete[] log_fold_change_error_sq_bc[b];
    }
    delete[] log_fold_change_buffer_c;
    delete[] log_fold_change_bc;
    delete[] log_fold_change_error_sq_bc;
    delete[] inv_curv_diag_rank1_c;
    delete[] inv_curv_diag_num_c;
    delete[] inv_curv_diag_den2_c;
    delete[] mean_ltq_b;

    RowComputation result;
    result.mean_ltq = mean_ltq;
    result.mean_ltq_error_sq = mean_ltq_error_sq;
    result.log_fold_change_c = log_fold_change_c;
    result.log_fold_change_error_sq_c = log_fold_change_error_sq_c;
    result.gene_variance = gene_variance;
    result.variance_likelihood_b = variance_likelihood_b;
    
    return result;
}

// Zero-count error bar (SI Eq. 38): the squared displacement epsilon_gc^2 above delta*_gc at which the log-likelihood has dropped by 1/2.
// log_fold_change = delta*_gc, gene_variance = v_g, gene_total = n_g and fraction = f_gc.
double zero_count_log_fold_change_error_sq(double log_fold_change, double gene_variance, double gene_total, double fraction)
{
    // Find the positive displacement whose log-likelihood drop is 1/2.
    // Keep log(1 + f*(exp(e)-1)) rather than linearizing the logarithm.
    auto log_likelihood_drop = [&](double displacement) { // sigma_gc in SI Eq. 38
        const double logarithm = displacement < 700 ? std::log1p(fraction * std::expm1(displacement))
                                        : displacement + std::log(fraction + (1 - fraction) * std::exp(-displacement));
        return displacement * (2 * log_fold_change + displacement) / (2 * gene_variance) + gene_total * logarithm;
    };
    double lo = 0, hi = 1;
    // The upper bound for the linearized formula need not bracket this root.
    while (log_likelihood_drop(hi) < 0.5)
    {
        hi *= 2;
        if (!std::isfinite(hi)) throw std::runtime_error("Cannot bracket zero-count uncertainty");
    }
    for (int iteration = 0; iteration < 200; ++iteration)
    {
        const double displacement = (lo + hi) / 2;
        const double value = log_likelihood_drop(displacement);
        if (!std::isfinite(value)) throw std::runtime_error("Nonfinite zero-count likelihood");
        if (std::fabs(value - 0.5) <= 1e-7) return displacement * displacement;
        if (value < 0.5) lo = displacement; else hi = displacement;
        if (hi - lo <= 1e-14 * std::max(1.0, hi)) return displacement * displacement;
    }
    throw std::runtime_error("Zero-count uncertainty did not converge");
}

ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &num_threads, bool &print_extended_output, double &variance_min, double &variance_max, int &num_bins, bool &no_norm, int &variance_method, bool &gzip_output, bool &npy_output, bool &bonsai_output_mode)
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
        variance_method = 1;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the likelihood (MLE).");
    }
    else if(variance_method_str == "MAP" || variance_method_str == "map" || variance_method_str == "MaxAPosterior" || variance_method_str == "maxaposterior" || variance_method_str == "max_a_posterior"){
        variance_method = 2;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the posterior (MAP).");
    }
    else if(variance_method_str == "EAP" || variance_method_str == "eap" || variance_method_str == "ExpectedAPosterior" || variance_method_str == "expectedaposterior" || variance_method_str == "expected_a_posterior"){
        variance_method = 3;
        logging_debug("Outputting results for the expected value of the prior variance (v_g) over the posterior (EAP).");
    }
    else if(variance_method_str == "MARG" || variance_method_str == "marg" || variance_method_str == "marginalizing"){
        variance_method = 0;
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
        if (gzip_output || npy_output || variance_method != 2)
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
    in_file_extension = (in_file.size() >= 3) ? in_file.substr(in_file.size() - 3) : in_file;

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
