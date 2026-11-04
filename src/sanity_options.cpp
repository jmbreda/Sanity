#include "sanity_options.hpp"

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "logging.hpp"

namespace
{

enum class Option
{
    File,
    Destination,
    Threads,
    ExtendedOutput,
    VarianceMin,
    VarianceMax,
    NumBins,
    GeneNameFile,
    CellNameFile,
    VarianceMethod,
    GzipOutput,
    NpyOutput,
    BonsaiOutputMode
};

struct OptionSpec
{
    Option option;
    std::vector<std::string> spellings;
    bool takes_value; // Flags take no value; every other option takes exactly one.
};

// Every option the command line accepts. -h/--help and -v/--version are handled before this table is used.
const std::vector<OptionSpec> option_table = {
    {Option::File, {"-f", "--file"}, true},
    {Option::Destination, {"-d", "--destination"}, true},
    {Option::Threads, {"-n", "--n_threads"}, true},
    {Option::ExtendedOutput, {"-e", "--extended_output"}, true},
    {Option::VarianceMin, {"-vmin", "--variance_min"}, true},
    {Option::VarianceMax, {"-vmax", "--variance_max"}, true},
    {Option::NumBins, {"-nbin", "--number_of_variance_bins"}, true},
    {Option::GeneNameFile, {"-mtx_genes", "--mtx_gene_name_file"}, true},
    {Option::CellNameFile, {"-mtx_cells", "--mtx_cell_name_file"}, true},
    {Option::VarianceMethod, {"-v_m", "--v_method"}, true},
    // -gz and -npy also accept their older spellings.
    {Option::GzipOutput, {"-gz", "--gzip_output", "--gz", "--gzip-output"}, false},
    {Option::NpyOutput, {"-npy", "--npy_output", "--npy", "--npy-output"}, false},
    {Option::BonsaiOutputMode, {"--bonsai-output-mode"}, false},
};

// Returns the table entry matched by `arg`, or nullptr if `arg` is not a known option.
const OptionSpec *match_option(const std::string &arg)
{
    for (const OptionSpec &spec : option_table)
        for (const std::string &spelling : spec.spellings)
            if (arg == spelling)
                return &spec;
    return nullptr;
}

bool parse_variance_method(const std::string &name, VarianceMethod &variance_method)
{
    if (name == "MLE" || name == "mle" || name == "MaxLikelihood" || name == "maxlikelihood" || name == "max_likelihood")
    {
        variance_method = VarianceMethod::MLE;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the likelihood (MLE).");
    }
    else if (name == "MAP" || name == "map" || name == "MaxAPosterior" || name == "maxaposterior" || name == "max_a_posterior")
    {
        variance_method = VarianceMethod::MAP;
        logging_debug("Outputting results for the prior variance (v_g) that maximizes the posterior (MAP).");
    }
    else if (name == "EAP" || name == "eap" || name == "ExpectedAPosterior" || name == "expectedaposterior" || name == "expected_a_posterior")
    {
        variance_method = VarianceMethod::EAP;
        logging_debug("Outputting results for the expected value of the prior variance (v_g) over the posterior (EAP).");
    }
    else if (name == "MARG" || name == "marg" || name == "marginalizing")
    {
        variance_method = VarianceMethod::MARG;
        logging_debug("Outputting results for the prior variance (v_g) marginalizing over v.");
    }
    else
    {
        logging_debug("Invalid method specified for variance estimation.");
        return false;
    }
    return true;
}

// "mtx" for .mtx and .mtx.gz files (in any letter case), and "tsv" for everything else.
std::string detect_input_type(const std::string &in_file)
{
    std::string extension = in_file.size() >= 3 ? in_file.substr(in_file.size() - 3) : in_file;
    if (in_file.size() >= 7)
    {
        const std::string last_7 = in_file.substr(in_file.size() - 7);
        if (last_7 == ".mtx.gz")
            return "mtx";
        if (last_7 == ".tsv.gz")
            return "tsv";
    }
    for (char &ch : extension)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return extension == "mtx" ? "mtx" : "tsv";
}

} // namespace

ParseResult parse_options(int argc, char **argv, SanityOptions &options)
{
    // Running with no arguments at all is a usage error, not a help request: main prints the usage
    // for both, but Error exits 1 whereas an explicit -h/--help still exits 0.
    if (argc < 2)
    {
        logging_debug("Error in argument parsing :\n"
                      "Not enough arguments provided.\n");
        return ParseResult::Error;
    }
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
            return ParseResult::HelpRequested;
    }
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "-v" || arg == "--version")
            return ParseResult::VersionRequested;
    }

    bool have_input_file = false;
    std::string extended_output("false");
    std::string variance_method_str("MAP");
    for (int i = 1; i < argc; i++)
    {
        // Anything that is not a known option, or the value of one, is rejected. Without this, a
        // typo such as "-e1" or "--bogus 7" is silently ignored and the run quietly produces
        // something other than what was asked for.
        const OptionSpec *spec = match_option(argv[i]);
        if (spec == nullptr)
        {
            logging_debug("Error in argument parsing :\nunrecognised option '" + std::string(argv[i]) + "'\n");
            return ParseResult::Error;
        }
        if (!spec->takes_value)
        {
            if (spec->option == Option::GzipOutput)
                options.gzip_output = true;
            else if (spec->option == Option::NpyOutput)
                options.npy_output = true;
            else if (spec->option == Option::BonsaiOutputMode)
                options.bonsai_output_mode = true;
            continue;
        }
        // The next token is this option's value, whatever it looks like: a value may legitimately
        // begin with '-', and no option here takes a negative number anyway.
        if (i + 1 >= argc)
        {
            logging_debug("Error in argument parsing :\n" + std::string(argv[i]) + " option missing\n");
            return ParseResult::Error;
        }
        const std::string value = argv[++i];
        switch (spec->option)
        {
        case Option::File:
            options.in_file = value;
            have_input_file = true;
            break;
        case Option::Destination:
            // `-d ""` means the current directory; otherwise make sure the folder ends in '/'.
            options.out_folder = value.empty() ? "./" : (value.back() == '/' ? value : value + '/');
            break;
        case Option::Threads:
            options.num_threads = std::atoi(value.c_str());
            break;
        case Option::ExtendedOutput:
            extended_output = value;
            break;
        case Option::VarianceMin:
            options.variance_min = std::atof(value.c_str());
            break;
        case Option::VarianceMax:
            options.variance_max = std::atof(value.c_str());
            break;
        case Option::NumBins:
            options.num_bins = std::atoi(value.c_str());
            break;
        case Option::GeneNameFile:
            options.gene_name_file = value;
            break;
        case Option::CellNameFile:
            options.cell_name_file = value;
            break;
        case Option::VarianceMethod:
            variance_method_str = value;
            break;
        default:
            break;
        }
    }
    if (!have_input_file)
    {
        logging_debug("Error in argument parsing :\n"
                      "missing input file name\n");
        return ParseResult::Error;
    }

    options.extended_output = (extended_output == "true" || extended_output == "1");
    if (!parse_variance_method(variance_method_str, options.variance_method))
        return ParseResult::Error;

    // Bonsai output mode writes one fixed set of files, whatever -e, -gz and -npy say. It does need
    // MAP: that is a choice of inference, not of output, so it is not overridden silently.
    if (options.bonsai_output_mode)
    {
        if (options.variance_method != VarianceMethod::MAP)
        {
            logging_debug("Bonsai output mode requires MAP; use -v_m MAP.");
            return ParseResult::Error;
        }
        if (options.extended_output || options.gzip_output || options.npy_output)
            logging_debug("Bonsai output mode ignores -e, -gz and -npy.");
    }

    options.in_file_extension = detect_input_type(options.in_file);
    logging_debug("File type : " + options.in_file_extension + "\n");
    return ParseResult::Continue;
}

void show_usage()
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
              << "\t-e,--extended_output\tOption to print extended output (default: false, choice: false,0,true,1)\n"
              << "\t-gz,--gzip_output\tCompress every text output file with gzip (.gz)\n"
              << "\t-npy,--npy_output\tWrite the numerical output files as float64 NumPy arrays (.npy) instead of text\n"
              << "\t--bonsai-output-mode\tWrite only the four numerical files and IDs/metadata required by Bonsai (NPY, MAP); ignores -e, -gz and -npy\n"
              << "\t-vmin,--variance_min\tMinimal value of variance in log transcription quotient (default: 0.001)\n"
              << "\t-vmax,--variance_max\tMaximal value of variance in log transcription quotient (default: 50)\n"
              << "\t-nbin,--number_of_variance_bins\tNumber of bins for the variance in log transcription quotient  (default: 160)\n"
              << "\t-v_m,--v_method\t\tOption to specify the method for variance estimation (default: MAP, choice: MAP, EAP, MLE, MARG)\n";
}

const char *variance_method_name(VarianceMethod variance_method)
{
    static const char *names[] = {"MARG", "MLE", "MAP", "EAP"};
    return names[static_cast<int>(variance_method)];
}
