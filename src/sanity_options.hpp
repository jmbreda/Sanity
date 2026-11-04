#ifndef SANITY_OPTIONS_HPP
#define SANITY_OPTIONS_HPP

#include <string>

#include "gene_fit.hpp"

// Everything the command line of Sanity can set, with its default.
struct SanityOptions
{
    std::string in_file;                    // The count matrix (-f).
    std::string in_file_extension;          // "mtx" or "tsv", derived from in_file.
    std::string gene_name_file = "none";    // Gene names of an MTX input (-mtx_genes).
    std::string cell_name_file = "none";    // Cell names of an MTX input (-mtx_cells).
    std::string out_folder = "./";          // Output folder, always ending in '/' (-d).
    int num_threads = 4;                    // -n
    bool extended_output = false;           // -e
    double variance_min = 0.001;            // v_min, the smallest candidate gene variance (SI Section S1.1.4).
    double variance_max = 50.0;             // v_max, the largest candidate gene variance.
    int num_bins = 160;                     // B, the number of grid points of the gene variance between v_min and v_max.
    VarianceMethod variance_method = VarianceMethod::MAP; // How to estimate v_g (see gene_fit.hpp).
    bool gzip_output = false;               // -gz
    bool npy_output = false;                // -npy
    bool bonsai_output_mode = false;        // --bonsai-output-mode
};

enum class ParseResult
{
    Continue,
    VersionRequested,
    HelpRequested,
    Error
};

// Fills `options` from the command line. Errors are logged before Error is returned.
ParseResult parse_options(int argc, char **argv, SanityOptions &options);

void show_usage();

// The name of a variance method as written to sanity_command.txt ("MARG", "MLE", "MAP" or "EAP").
const char *variance_method_name(VarianceMethod variance_method);

#endif
