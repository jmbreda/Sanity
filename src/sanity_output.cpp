#include "sanity_output.hpp"

#include <cerrno>
#include <charconv>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h> // For mkdir and stat on Unix-like systems
#include <system_error>

#ifdef _WIN32
    #include <direct.h> // For _mkdir on Windows
#endif

#include "logging.hpp"
#include "Version.h"

namespace
{

// The base name of each field's file; the extension depends on the output format.
const char *field_file_name(OutputField field)
{
    switch (field)
    {
    case OutputField::Ltq: return "log_transcription_quotients";
    case OutputField::LtqError: return "ltq_error_bars";
    case OutputField::LogFoldChange: return "delta";
    case OutputField::LogFoldChangeError: return "d_delta";
    case OutputField::MeanLtq: return "mu";
    case OutputField::MeanLtqError: return "d_mu";
    case OutputField::GeneVariance: return "variance";
    case OutputField::VarianceLikelihood: return "likelihood";
    }
    throw std::logic_error("Unknown output field");
}

// Text lines of these fields start with the gene name; the others hold only numbers.
bool text_line_starts_with_gene(OutputField field)
{
    return field == OutputField::Ltq || field == OutputField::LtqError || field == OutputField::VarianceLikelihood;
}

// Fills `values` with the numbers of one gene's row of `field`, in file order.
void field_values(OutputField field, const RowComputation &result, std::vector<double> &values)
{
    const std::size_t num_cells = result.log_fold_change_c.size();
    values.clear();
    switch (field)
    {
    case OutputField::Ltq:
        for (std::size_t c = 0; c < num_cells; ++c)
            values.push_back(result.mean_ltq + result.log_fold_change_c[c]);
        break;
    case OutputField::LtqError:
        for (std::size_t c = 0; c < num_cells; ++c)
            values.push_back(std::sqrt(result.mean_ltq_error_sq + result.log_fold_change_error_sq_c[c]));
        break;
    case OutputField::LogFoldChange:
        values = result.log_fold_change_c;
        break;
    case OutputField::LogFoldChangeError:
        for (double error_sq : result.log_fold_change_error_sq_c)
            values.push_back(std::sqrt(error_sq));
        break;
    case OutputField::MeanLtq:
        values.push_back(result.mean_ltq);
        break;
    case OutputField::MeanLtqError:
        values.push_back(std::sqrt(result.mean_ltq_error_sq));
        break;
    case OutputField::GeneVariance:
        values.push_back(result.gene_variance);
        break;
    case OutputField::VarianceLikelihood:
        values = result.variance_likelihood_b;
        break;
    }
}

// Appends `value` with six decimals, the precision of all text output.
void append_fixed_number(std::string &output, double value)
{
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite output value");
    char buffer[512];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::fixed, 6);
    if (result.ec != std::errc{}) throw std::runtime_error("Could not format output value");
    output.append(buffer, result.ptr);
}

std::string format_text_line(OutputField field, const std::vector<double> &values, const std::string &gene_name)
{
    std::string line;
    line.reserve(values.size() * 12 + gene_name.size() + 1);
    const bool with_gene = text_line_starts_with_gene(field);
    if (with_gene) line = gene_name;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (with_gene || i > 0) line.push_back('\t');
        append_fixed_number(line, values[i]);
    }
    line.push_back('\n');
    return line;
}

void write_names(const std::string &path, const std::vector<std::string> &names)
{
    sanity::Writer writer(path);
    for (const std::string &name : names)
        writer.write_raw(name + "\n");
    writer.close();
}

} // namespace

bool prepare_output_folder(const SanityOptions &options)
{
    const std::string &out_folder = options.out_folder;
    int mkdir_result = 0;
    #ifdef _WIN32
        mkdir_result = _mkdir(out_folder.c_str());
    #else
        mkdir_result = mkdir(out_folder.c_str(), 0755);
    #endif

    // An existing folder is fine.
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
        return false;
    }

    if (options.bonsai_output_mode)
    {
        for (const char *filename : {"log_transcription_quotients.txt", "ltq_error_bars.txt",
                                     "d_mu.txt", "likelihood.txt", "delta.txt", "d_delta.txt",
                                     "mu.txt", "variance.txt", "log_transcription_quotients.txt.gz",
                                     "ltq_error_bars.txt.gz", "d_mu.txt.gz", "likelihood.txt.gz",
                                     "delta.txt.gz", "d_delta.txt.gz", "mu.txt.gz", "variance.txt.gz"})
        {
            struct stat previous;
            if (stat((out_folder + filename).c_str(), &previous) == 0)
            {
                logging_debug("Bonsai output mode requires a folder without existing text-output files: " + out_folder);
                return false;
            }
        }
    }
    return true;
}

void write_command_metadata(const SanityOptions &options, int argc, char **argv)
{
    std::ofstream cmd_file(options.out_folder + "sanity_command.txt", std::ios::trunc);
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
    cmd_file << "# Method: " << variance_method_name(options.variance_method) << "\n";
    if (options.bonsai_output_mode)
        cmd_file << "# Output mode: bonsai\n"
                    "# Output: Bonsai NPY; delta/d_delta: <f4; mu/variance: <f8; "
                    "C-order genes x cells; gene vectors: genes x 1; d_delta: standard deviation\n";
    cmd_file << argv[0];
    for (int i = 1; i < argc; ++i)
        cmd_file << " " << argv[i];
    cmd_file << "\n";
}

SanityOutput::SanityOutput(const SanityOptions &options, const std::vector<std::string> &gene_names_g,
                           const std::vector<std::string> &cell_names_c, const std::vector<double> &variance_grid_b)
    : out_folder_(options.out_folder)
{
    const std::size_t num_cells = cell_names_c.size();
    if (options.npy_output)
    {
        // The per-cell matrices are float32 in Bonsai output mode and float64 otherwise.
        const FileFormat matrix_format = options.bonsai_output_mode ? FileFormat::Float32Npy : FileFormat::Float64Npy;
        add_field(OutputField::LogFoldChange, matrix_format, "", num_cells);
        add_field(OutputField::LogFoldChangeError, matrix_format, "", num_cells);
        add_field(OutputField::MeanLtq, FileFormat::Float64Npy, "", 1);
        add_field(OutputField::GeneVariance, FileFormat::Float64Npy, "", 1);
        write_names(out_folder_ + "geneID.txt", gene_names_g);
        write_names(out_folder_ + "cellID.txt", cell_names_c);
        return;
    }

    const std::string suffix = options.gzip_output ? ".gz" : "";
    add_field(OutputField::Ltq, FileFormat::Text, suffix, num_cells);
    add_field(OutputField::LtqError, FileFormat::Text, suffix, num_cells);
    if (options.extended_output)
    {
        add_field(OutputField::MeanLtq, FileFormat::Text, suffix, 1);
        add_field(OutputField::MeanLtqError, FileFormat::Text, suffix, 1);
        add_field(OutputField::GeneVariance, FileFormat::Text, suffix, 1);
        add_field(OutputField::LogFoldChange, FileFormat::Text, suffix, num_cells);
        add_field(OutputField::LogFoldChangeError, FileFormat::Text, suffix, num_cells);
        add_field(OutputField::VarianceLikelihood, FileFormat::Text, suffix, variance_grid_b.size());
        write_names(out_folder_ + "geneID.txt" + suffix, gene_names_g);
        write_names(out_folder_ + "cellID.txt" + suffix, cell_names_c);
    }

    // Header lines of the text tables: the cell names above the per-cell values, the grid above the likelihoods.
    for (FieldFile &file : files_)
    {
        std::ostringstream header;
        if (file.field == OutputField::Ltq || file.field == OutputField::LtqError)
        {
            header << "GeneID";
            for (const std::string &cell_name : cell_names_c)
                header << "\t" << cell_name;
        }
        else if (file.field == OutputField::VarianceLikelihood)
        {
            header << "Variance" << std::fixed << std::setprecision(6);
            for (double variance : variance_grid_b)
                header << "\t" << variance;
        }
        else
            continue;
        header << "\n";
        file.text->write_raw(header.str());
    }
}

void SanityOutput::add_field(OutputField field, FileFormat format, const std::string &text_suffix, std::size_t num_columns)
{
    FieldFile file;
    file.field = field;
    file.format = format;
    const std::string path = out_folder_ + field_file_name(field);
    if (format == FileFormat::Text)
        file.text = std::make_unique<sanity::Writer>(path + ".txt" + text_suffix);
    else
        file.npy = std::make_unique<npy::NpyStreamWriter>(path + ".npy", num_columns,
            format == FileFormat::Float32Npy ? npy::NpyDType::Float32 : npy::NpyDType::Float64);
    files_.push_back(std::move(file));
}

SanityOutput::PreparedRow SanityOutput::prepare_row(const RowComputation &result, const std::string &gene_name) const
{
    PreparedRow row;
    row.text.resize(files_.size());
    row.float64.resize(files_.size());
    row.float32.resize(files_.size());
    std::vector<double> values; // Reused for every file of this gene.
    for (std::size_t f = 0; f < files_.size(); ++f)
    {
        const FieldFile &file = files_[f];
        field_values(file.field, result, values);
        if (file.format == FileFormat::Text)
        {
            row.text[f] = format_text_line(file.field, values, gene_name);
            continue;
        }
        for (double value : values)
            if (!std::isfinite(value))
                throw std::runtime_error("Nonfinite NPY output");
        if (file.format == FileFormat::Float64Npy)
            row.float64[f] = values;
        else
        {
            row.float32[f].assign(values.begin(), values.end());
            for (float value : row.float32[f])
                if (!std::isfinite(value))
                    throw std::runtime_error("NPY output exceeds float32 range");
        }
    }
    return row;
}

void SanityOutput::write_row(const PreparedRow &row)
{
    for (std::size_t f = 0; f < files_.size(); ++f)
    {
        FieldFile &file = files_[f];
        if (file.format == FileFormat::Text)
            file.text->write_raw(row.text[f]);
        else if (file.format == FileFormat::Float64Npy)
            file.npy->write_row(row.float64[f]);
        else
            file.npy->write_row(row.float32[f]);
    }
}

void SanityOutput::close()
{
    for (FieldFile &file : files_)
    {
        if (file.text) file.text->close();
        if (file.npy) file.npy->close();
    }
}
