#include <cctype>
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
#include "FitFrac.h"
#include "Digamma_Trigamma.h"
#include "Writer.hpp"

std::string VERSION("2.0");
enum ParseResult
{
    CONTINUE,
    VERSION_REQUESTED,
    HELP_REQUESTED,
    ERROR
};

struct RowComputation
{
    double mu;
    double var_mu;
    std::vector<double> delta;
    std::vector<double> var_delta;
    double var_gene;
    std::vector<double> lik;
};

/***Function declarations ****/
RowComputation get_gene_expression_level(const std::vector<double> &n_c, const std::vector<double> &N_c, double n, double vmin, double vmax, int C, int numbin, double a, double b, int v_method);
double get_epsilon_2(double &d, double &v, double &n, double &f, double &a);
ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &N_threads, bool &print_extended_output, double &vmin, double &vmax, int &numbin, bool &no_norm, int &v_method, bool &gzip_output);
static void show_usage(void);
std::vector<double> fetch_row(int g, FileReader &infile, const std::string &in_file_extension, const std::vector<RowBlock> &mtx_rows, const std::vector<std::streampos> &tsv_offsets, const int &C);

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
    int N_threads(4);
    bool print_extended_output(false);
    double vmin = 0.001;
    double vmax = 50.0;
    int numbin = 160;
    bool no_norm(false);
    int v_method = 2; // default is to output MAP: the maximum a posteriori estimate of v
    bool gzip_output = false;
    std::string out_suffix = "";

    ParseResult parse_res = parse_argv(argc, argv, in_file, gene_name_file, cell_name_file, in_file_extension, out_folder, N_threads, print_extended_output, vmin, vmax, numbin, no_norm, v_method, gzip_output);
    if (parse_res == HELP_REQUESTED)
    {
        show_usage();
        return 0;
    }
    else if (parse_res == VERSION_REQUESTED)
    {
        std::cout << "Sanity version " << VERSION << "\n";
        return 0;
    }
    else if (parse_res == ERROR)
    {
        show_usage();
        return 1;
    }
    // count Number of genes and cells
    int G, C;
    // Number of rows in file
    int N_rows;
    // Gene idx map for mtx
    std::map<int, int> gene_idx;
    std::vector<RowBlock> mtx_rows;
    std::vector<std::streampos> tsv_offsets;
    std::vector<double> N_c; // total counts per cell
    std::vector<double> n;   // total counts per gene
    std::vector<std::string> gene_names;
    std::vector<std::string> cell_names;
    logging_debug("Reading input file to get gene and cell counts");
    if (in_file_extension == "mtx")
    {
        Get_G_C_MTX(in_file, N_rows, G, C, gene_idx, mtx_rows, N_c, n);
        gene_names = Read_GeneNames(gene_name_file, gene_idx, G);
        if (static_cast<int>(gene_names.size()) == 0)
        {
            logging_debug("Warning: gene names file not found or empty, generating default gene names.");
            gene_names.clear();
            for (int g = 0; g < G; ++g)
            {
                gene_names.push_back("Gene_" + std::to_string(g + 1));
            }
        }
        cell_names = Read_CellNames(cell_name_file);
        if (static_cast<int>(cell_names.size()) == 0)
        {
            logging_debug("Warning: cell names file not found or empty, generating default cell names.");
            for (int c = 0; c < C; ++c)
            {
                cell_names.push_back("Cell_" + std::to_string(c + 1));
            }
        }
        logging_debug("There were " + std::to_string(N_rows) + " rows");
    }
    else
    {
        Get_G_C_UMIcountMatrix(in_file, N_rows, G, C, tsv_offsets, N_c, n, cell_names, gene_names, N_threads);
        logging_debug("There were " + std::to_string(N_rows) + " rows");
    }
    logging_debug("There were " + std::to_string(G) + " genes and " + std::to_string(C) + " cells");

    int g, c, k;

    // Remove the total UMI correction if no cell size normalization option is true
    if (no_norm)
    {
        logging_debug("No cell size normalization will be performed");

        // get mean count per cell
        double mean_N_c = 0;
        for (c = 0; c < C; ++c)
        {
            mean_N_c += N_c[c];
        }
        mean_N_c /= C;

        // Now replce N_c by N
        for (c = 0; c < C; ++c)
        {
            N_c[c] = mean_N_c;
        }
    }

    // alpha and beta of gamma prior on mu
    double a = 1.0;
    double b = 0.0;
    double deltav = std::log(vmax / vmin) / ((double)numbin - 1);

    // create output folder if it does not exist
    if (out_folder == "/") {
        out_folder = "./";
    }
    int mkdir_result = 0;
    #ifdef _WIN32
        mkdir_result = _mkdir(out_folder.c_str());
    #else
        mkdir_result = mkdir(out_folder.c_str(), 0755);
    #endif
    
    // Check for errors (ignore if directory already exists)
    if (mkdir_result != 0 && errno != EEXIST)
    {
        logging_debug("Error: Failed to create output folder: " + out_folder);
        return 1;
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

    cmd_file << "# Sanity version: " << VERSION << "\n";

    static const char* vm_name[] = {"MARG", "MLE", "MAP", "EAP"};
    cmd_file << "# Method: " << vm_name[v_method] << "\n";
    cmd_file  << argv[0];
    for (int i = 1; i < argc; ++i)
    {
        cmd_file << " " << argv[i];
    }
    cmd_file << "\n";
    cmd_file.close();

    // -- Open files for writing --
    if (gzip_output) {out_suffix = ".gz";}

    sanity::Writer out_exp_lev, out_d_exp_lev, out_mu, out_dmu, out_var_gene,
        out_delta, out_ddelta, out_lik, out_gene, out_cell;

    out_exp_lev.open(out_folder + "log_transcription_quotients.txt" + out_suffix);
    out_exp_lev << std::fixed << std::setprecision(6);

    out_d_exp_lev.open(out_folder + "ltq_error_bars.txt" + out_suffix);
    out_d_exp_lev << std::fixed << std::setprecision(6);

    out_exp_lev << "GeneID";
    out_d_exp_lev << "GeneID";
    for (c = 0; c < C; c++)
    {
        out_exp_lev << "\t" << cell_names[c].c_str();
        out_d_exp_lev << "\t" << cell_names[c].c_str();
    }
    out_exp_lev << "\n";
    out_d_exp_lev << "\n";
    if (print_extended_output)
    {

        out_gene.open(out_folder + "geneID.txt" + out_suffix);
        out_cell.open(out_folder + "cellID.txt" + out_suffix);
        out_mu.open(out_folder + "mu.txt" + out_suffix);
        out_mu << std::fixed << std::setprecision(6);

        out_dmu.open(out_folder + "d_mu.txt" + out_suffix);
        out_dmu << std::fixed << std::setprecision(6);

        out_var_gene.open(out_folder + "variance.txt" + out_suffix);
        out_var_gene << std::fixed << std::setprecision(6);

        out_delta.open(out_folder + "delta.txt" + out_suffix);
        out_delta << std::fixed << std::setprecision(6);

        out_ddelta.open(out_folder + "d_delta.txt" + out_suffix);
        out_ddelta << std::fixed << std::setprecision(6);

        out_lik.open(out_folder + "likelihood.txt" + out_suffix);
        out_lik << std::fixed << std::setprecision(6);

        out_lik << "Variance";
        for (k = 0; k < (numbin); ++k)
        {
            out_lik << "\t" << vmin * std::exp(deltav * k);
        }
        out_lik << "\n";

        // save cell names
        for (c = 0; c < C; c++)
        {
            out_cell << cell_names[c].c_str() << "\n";
        }
    }

    logging_debug("Fit gene expression levels");
    const std::clock_t begin = std::clock();
    #pragma omp parallel num_threads(N_threads)
    {
        FileReader thread_reader(in_file);
        #pragma omp for schedule(dynamic) ordered
        for (int g = 0; g < G; ++g)
        {
            std::vector<double> n_c_g = fetch_row(g, thread_reader, in_file_extension, mtx_rows, tsv_offsets, C);
            RowComputation result = get_gene_expression_level(n_c_g, N_c, n[g], vmin, vmax, C, numbin, a, b, v_method);
            #pragma omp ordered
            {
                // output esimated running time
                if (g == (3 * N_threads - 1))
                {
                    std::clock_t end = std::clock();
                    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    double estimated_running_time = (elapsed_secs / (3. * N_threads)) * G / N_threads;

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
                    logging_debug("Finished " + std::to_string(g) + " genes out of " + std::to_string(G));
                }

                out_exp_lev << gene_names[g];
                out_d_exp_lev << gene_names[g];
                for (c = 0; c < C; c++)
                {
                    out_exp_lev << "\t" << result.mu + result.delta[c];
                    out_d_exp_lev << "\t" << std::sqrt(result.var_mu + result.var_delta[c]);
                    if (print_extended_output)
                    {
                        out_delta << result.delta[c];
                        out_ddelta << std::sqrt(result.var_delta[c]);
                        if (c < C - 1)
                        {
                            out_delta << "\t";
                            out_ddelta << "\t";
                        }
                    }
                }
                out_exp_lev << "\n";
                out_d_exp_lev << "\n";

                if (print_extended_output)
                {
                    out_delta << "\n";
                    out_ddelta << "\n";
                    // Write gene names
                    out_gene << gene_names[g].c_str() << "\n";
                    // print best fit to file : mu, delta
                    //  Print diagonal of invM : variance of mu, delta
                    out_mu << result.mu << "\n";
                    out_dmu << std::sqrt(result.var_mu) << "\n";
                    out_var_gene << result.var_gene << "\n";
                    // Write likelihood
                    out_lik << gene_names[g];
                    for (k = 0; k < numbin; ++k)
                    {
                        out_lik << "\t" << result.lik[k];
                    }
                    out_lik << "\n";
                }
            }
        }
    }

    logging_debug("Finished fitting all genes");

    return 0;
}

RowComputation get_gene_expression_level(const std::vector<double> &n_c, const std::vector<double> &N_c,
                                         double n, double vmin, double vmax, int C, int numbin, double a, double b,
                                         int v_method)
{
    // n = total counts for the gene
    // n_c = counts for the gene in each cell
    // N_c = total counts for each cell
    // C = number of cells
    // numbin = number of bins for v
    int i, k;
    double beta, L, ldet, q, delsq, inv_v;
    double prev_q = 0.0;
    double *f = new double[C];
    double **delta_v = new double *[numbin];
    double **sig2_delta_v = new double *[numbin];
    std::vector<double> lik(numbin, -1.0);
    std::vector<double> delta(C, 0.0), var_delta(C, 0.0);
    double var_mu, var_gene;

    for (k = 0; k < numbin; ++k)
    {
        delta_v[k] = new double[C];
        sig2_delta_v[k] = new double[C];
    }

    /*** To compute var of delta ***/
    double *sig2_delta_c = new double[C];
    double *sig2_delta_num = new double[C];
    double *sig2_delta_den2 = new double[C];
    double sig2_delta_den1;

    double *mu_v = new double[numbin];
    double Lmax = -1e+100;
    int Lmax_ind = 0;
    double v;
    double deltav;
    deltav = std::log(vmax / vmin) / ((double)numbin - 1);

    for (k = 0; k < numbin; ++k)
    {
        v = vmin * std::exp(deltav * k);
        beta = 1.0 / (n * v);
        q = fitfrac(f, n_c, n, v, C, N_c, a, b, prev_q);
        prev_q = q;
        mu_v[k] = Psi_0(n) - q; /*** equation (85) ***/

        delsq = 0;
        L = -0.5 * ((double)C) * std::log(v); // (56) 1st term
        for (i = 0; i < C; ++i)
        {
            delta_v[k][i] = std::log(f[i]) - std::log(N_c[i]) + q; // equation (67)
            L += n_c[i] * delta_v[k][i];                 // Bug fix: remove a term as in Equation 19 of Sanity paper SI
            delsq += delta_v[k][i] * delta_v[k][i];
        }
        L -= delsq / (2 * v); // (56) 2nd term
        L -= n * q;     // 4th term in equation (56)

        // get the determinant of the matrix
        ldet = 0.0;
        for (i = 0; i < C; ++i)
        {
            ldet += (f[i] * f[i]) / (f[i] + beta);
        }

        ldet = std::log(1 - ldet);
        for (i = 0; i < C; ++i)
        {
            ldet += std::log(f[i] + beta);
        }
        L -= 0.5 * ldet;
        // substract prior with a = 1, b = 1 ( log(v^a*exp(-b*v)) = alog(v) - bv
        lik[k] = L;

        if (L > Lmax)
        {
            Lmax = L;
            Lmax_ind = k;
        }

        /* compute nf^2/(nf+1/sigma^2) for each c */
        inv_v = 1.0 / v;
        for (i = 0; i < C; ++i)
        {
            sig2_delta_c[i] = n * f[i] * f[i] / (n * f[i] + inv_v);
        }
        /* Compute the full sum of the denominator in Delta_delta and the second tern in the denominator*/
        sig2_delta_den1 = 1.0;
        for (i = 0; i < C; ++i)
        {
            sig2_delta_den1 -= sig2_delta_c[i];
            sig2_delta_den2[i] = n * f[i] + inv_v;
        }
        /* compute the different terms in the numerator : remove the \tilde{c} terms */
        for (i = 0; i < C; i++)
        {
            sig2_delta_num[i] = sig2_delta_den1 + sig2_delta_c[i];
        }
        /* compute sig2_delta */
        for (i = 0; i < C; ++i)
        {
            sig2_delta_v[k][i] = sig2_delta_num[i] / (sig2_delta_den1 * sig2_delta_den2[i]);
        }

        // fix computation of asymmetric sig2_delta for zero count
        for (i = 0; i < C; ++i)
        {
            if (n_c[i] <= 0.5)
            {
                sig2_delta_v[k][i] = get_epsilon_2(delta_v[k][i], v, n, f[i], a);
            }
        }
    } // end v bins loop

    // get normalized likelihood from loglikelihood
    double sum_L;
    sum_L = 0.0;
    for (k = 0; k < numbin; k++)
    {
        lik[k] -= Lmax;
        lik[k] = std::exp(lik[k]);
        sum_L += lik[k];
    }
    for (k = 0; k < numbin; k++)
    {
        lik[k] /= sum_L;
    }

    int vindex = 0;
    if(v_method ==1){
        vindex = Lmax_ind;
    }
    else if(v_method == 2){
        double mapmax = -1e+100;
        for (k = 0; k < numbin; k++)    {
            double curv = vmin * std::exp(deltav * k);

            if (lik[k]/curv > mapmax)
            {
                mapmax = lik[k]/curv;
                vindex = k;
            }
        }
    }
    else if(v_method == 3){
        double postmean = 0.0;
        for (k = 0; k < numbin; k++)    {
            double curv = vmin * std::exp(deltav * k);
            postmean += lik[k] * curv;
        }
        double mindist = 1e+100;
        for (k = 0; k < numbin; k++)    {
            double curv = vmin * std::exp(deltav * k);
            if (std::fabs(curv - postmean) < mindist)
            {
                mindist = std::fabs(curv - postmean);
                vindex = k;
            }
        }
    }

    // Average delta, and mu
    double mu = 0.0;
    for (k = 0; k < numbin; k++)
    {
        mu += lik[k] * mu_v[k];
    }

    // Compute var_delta = < (mu - <mu>)^2 > + <d_mu>
    var_mu = Psi_1((double)n);
    for (k = 0; k < numbin; k++)
    {
        var_mu += lik[k] * (mu_v[k] - mu) * (mu_v[k] - mu);
    }

    if (v_method == 0)  // MARG
    {
        // Compute <delta> = int p(v)*delta(v) dv
        for (i = 0; i < C; i++)
        {
            delta[i] = 0.0;
            for (k = 0; k < numbin; k++)
            {
                delta[i] += lik[k] * delta_v[k][i];
            }
        }
        // Compute var_delta = < (delta - <delta>)^2 > + <sig2_delta^2>
        for (i = 0; i < C; i++)
        {
            var_delta[i] = 0.0;
            for (k = 0; k < numbin; k++)
            {
                var_delta[i] += lik[k] * (delta_v[k][i] - delta[i]) * (delta_v[k][i] - delta[i]) + lik[k] * sig2_delta_v[k][i];
            }
        }
        // var_gene
        var_gene = 0.0;
        for (k = 0; k < numbin; ++k)
        {
            v = vmin * std::exp(deltav * k);
            var_gene += v * lik[k];
        }
    }

    if (v_method > 0)
    {
        // Store the gene-variance that maximizes the likelihood:
        var_gene = vmin * std::exp(deltav * vindex);
        // And then also the corresponding values for the LTQs etc.
        mu = mu_v[vindex];
        var_mu = Psi_1((double)n);
        for (i = 0; i < C; i++)
        {
            delta[i] = delta_v[vindex][i];
        }
        for (i = 0; i < C; i++)
        {
            var_delta[i] = sig2_delta_v[vindex][i];
        }
    }

    delete[] f;
    for (k = 0; k < numbin; ++k)
    {
        delete[] delta_v[k];
        delete[] sig2_delta_v[k];
    }
    delete[] delta_v;
    delete[] sig2_delta_v;
    delete[] sig2_delta_c;
    delete[] sig2_delta_num;
    delete[] sig2_delta_den2;
    delete[] mu_v;

    RowComputation result;
    result.mu = mu;
    result.var_mu = var_mu;
    result.delta = delta;
    result.var_delta = var_delta;
    result.var_gene = var_gene;
    result.lik = lik;
    
    return result;
}

double get_epsilon_2(double &d, double &v, double &n, double &f, double &a)
{

    double e;
    double dL;
    double e_low = 0.0;
    double e_high = 0.0;
    double vnf = v * n * f;
    e_high = (-(d + vnf) + std::sqrt((d + vnf) * (d + vnf) + v * (1.0 + vnf))) / (1.0 + vnf);

    // bisection method :
    double tol = 0.0000001;
    double diff = 1.0;
    while (diff > tol)
    {
        e = (e_high + e_low) / 2.0;
        dL = e * (2.0 * d + e) / (2.0 * v) + n * f * (std::exp(e) - 1.0);
        if (dL < 0.5)
        {
            e_low = e;
        }
        else
        {
            e_high = e;
        }
        diff = std::fabs(dL - 0.5);
    }
    return e * e;
}

ParseResult parse_argv(int argc, char **argv, std::string &in_file, std::string &gene_name_file, std::string &cell_name_file, std::string &in_file_extension, std::string &out_folder, int &N_threads, bool &print_extended_output, double &vmin, double &vmax, int &numbin, bool &no_norm, int &v_method, bool &gzip_output)
{

    if (argc < 2)
    {
        logging_debug("Error in argument parsing :\n"
                      "Not enough arguments provided.\n");
        return HELP_REQUESTED;
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

    int N_param(12);
    std::string extended_output("false");
    std::string no_norm_str("false");
    std::string v_method_str("MAP");
    std::string to_find[12][2] = {{"-f", "--file"},
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
                             {"--gz", "--gzip-output"}};

    int j;
    int idx;
    for (j = 0; j < N_param; j++)
    {
        idx = 0;
        for (i = 1; i < argc; i++)
        {
            if (argv[i] == to_find[j][0] || argv[i] == to_find[j][1])
            {
                if (j == 11)
                {
                    gzip_output = true;
                    continue; // no argument expected for --gz
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
                    N_threads = std::atoi(argv[idx + 1]);
                if (j == 3)
                    extended_output = argv[idx + 1];
                if (j == 4)
                    vmin = std::atof(argv[idx + 1]);
                if (j == 5)
                    vmax = std::atof(argv[idx + 1]);
                if (j == 6)
                    numbin = std::atoi(argv[idx + 1]);
                if (j == 7)
                    gene_name_file = argv[idx + 1];
                if (j == 8)
                    cell_name_file = argv[idx + 1];
                if (j == 9)
                    no_norm_str = argv[idx + 1];
                if (j == 10)
                    v_method_str = argv[idx + 1];
                // add '/' to out_folder if not already
                if (j == 1 && out_folder.back() != '/')
                    out_folder = out_folder + '/';
            }
        }
        if (idx == 0 && j == 0)
        {
            logging_debug("Error in argument parsing :\n"
                          "missing input file name\n");
            return ERROR;
        }
    }

    if (extended_output == "true" || extended_output == "1")
        print_extended_output = true;

    if (no_norm_str == "true" || no_norm_str == "1")
        no_norm = true;

    if(v_method_str == "MLE" || v_method_str == "mle" || v_method_str == "MaxLikelihood" || v_method_str == "maxlikelihood" || v_method_str == "max_likelihood"){
        v_method = 1;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the likelihood (MLE).");
    }
    else if(v_method_str == "MAP" || v_method_str == "map" || v_method_str == "MaxAPosterior" || v_method_str == "maxaposterior" || v_method_str == "max_a_posterior"){
        v_method = 2;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the posterior (MAP).");
    }
    else if(v_method_str == "EAP" || v_method_str == "eap" || v_method_str == "ExpectedAPosterior" || v_method_str == "expectedaposterior" || v_method_str == "expected_a_posterior"){
        v_method = 3;
        logging_debug("Outputting results for the expected value of the prior variance (v_g) over the posterior (EAP).");
    }
    else if(v_method_str == "MARG" || v_method_str == "marg" || v_method_str == "marginalizing"){
        v_method = 0;
        logging_debug("Outputting results for the prior variance (v_g) marginalizing over v.");
    }
    else {
        logging_debug("Invalid method specified for variance estimation.");
        return ERROR;
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
              << "\t-f,--file\t\tSpecify the input transcript count text file (.mtx for Matrix Market File Format)\n"
              << "\t-mtx_genes,--mtx_gene_name_file\tSpecify the gene name text file (only needed if .mtx input file)\n"
              << "\t-mtx_cells,--mtx_cell_name_file\tSpecity the cell name text file (only needed if .mtx input file)\n"
              << "\t-d,--destination\tSpecify the destination path (default: pwd)\n"
              << "\t-n,--n_threads\t\tSpecify the number of threads to be used (default: 4)\n"
              << "\t-e,--extended_output\tOption to print extended output (default: false, choice: false,0,true,1)\n"
              << "\t-vmin,--variance_min\tMinimal value of variance in log transcription quotient (default: 0.001)\n"
              << "\t-vmax,--variance_max\tMaximal value of variance in log transcription quotient (default: 50)\n"
              << "\t-nbin,--number_of_bins\tNumber of bins for the variance in log transcription quotient  (default: 160)\n"
              << "\t-no_norm,--no_cell_size_normalization\tOption to skip cell size normalization (default: false, choice: false,0,true,1)\n"
              << "\t-v_m,--v_method\t\tOption to specify the method for variance estimation (default: MAP, choice: MAP, EAP, MLE, MARG)\n";
}

std::vector<double> fetch_row(int g, FileReader &thread_reader, const std::string &in_file_extension, const std::vector<RowBlock> &mtx_rows, const std::vector<std::streampos> &tsv_offsets, const int &C)
{
    std::vector<double> row_data(C, 0.0);
    if (!thread_reader.is_open())
    {
        throw std::runtime_error("RowReader: cannot open input file.");
    }

    if (in_file_extension == "mtx")
    {
        // Fetch from Matrix Market format
        const RowBlock &row_block = mtx_rows[g];
        thread_reader.seekg(row_block.offset);
        char *line = nullptr;
        char *saveptr = nullptr;
        char *token = nullptr;
        for (size_t idx = 0; idx < row_block.nnz; ++idx)
        {
            line = thread_reader.getline();
            token = strtok_r(line, " \t\r\n", &saveptr);
            int g_idx = std::stoi(token) - 1; // Convert to 0-based index
            token = strtok_r(NULL, " \t\r\n", &saveptr);
            int c_idx = std::stoi(token) - 1; // Convert to 0-based index
            token = strtok_r(NULL, " \t\r\n", &saveptr);
            double value = std::stod(token);
            row_data[c_idx] = value;
        }
    }
    else
    {
        // Fetch from TSV format
        thread_reader.seekg(tsv_offsets[g]);
        char *line = nullptr;
        char *saveptr = nullptr;
        char *token = nullptr;
        line = thread_reader.getline();
        token = strtok_r(line, "\t\r\n", &saveptr); // skip gene name
        for (int c = 0; c < C; ++c)
        {
            token = strtok_r(NULL, "\t\r\n", &saveptr);
            row_data[c] = std::stod(token);
        }
    }

    return row_data;
}
