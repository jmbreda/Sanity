#ifndef SANITY_OUTPUT_HPP
#define SANITY_OUTPUT_HPP

#include <memory>
#include <string>
#include <vector>

#include "gene_fit.hpp"
#include "helpers/io/npy_writer.hpp"
#include "helpers/io/writer.hpp"
#include "sanity_options.hpp"

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.

// Creates the output folder (its parent must exist) and, in Bonsai output mode, checks that it holds no
// earlier text output. Logs the reason and returns false if the run cannot write there.
bool prepare_output_folder(const SanityOptions &options);

// Writes sanity_command.txt: timestamp, version, method, output mode and the command line.
void write_command_metadata(const SanityOptions &options, int argc, char **argv);

// One numerical result per gene, each written to its own file.
enum class OutputField
{
    Ltq,                // x_gc = mu_g + delta_gc, the log transcription quotients (SI Eq. 45 and 39).
    LtqError,           // sqrt((delta mu_g)^2 + epsilon_gc^2), the error bars on x_gc (SI Eqs. 42, 47).
    LogFoldChange,      // delta_gc (SI Eq. 39).
    LogFoldChangeError, // epsilon_gc, the standard deviation of delta_gc (SI Eq. 42).
    MeanLtq,            // mu_g (SI Eq. 45).
    MeanLtqError,       // delta mu_g, the standard deviation of mu_g (SI Eq. 47).
    GeneVariance,       // v_g (SI Eq. 17).
    VarianceLikelihood  // P(v_b | n_g) on the grid of variances (SI Eq. 34).
};

// The output files of one run. The constructor opens every file and writes the headers and identifier
// files; prepare_row may then run in any thread, while write_row must be called once per gene, in gene order.
class SanityOutput
{
public:
    // The content of one gene's row for every output file, prepared outside the ordered write section.
    struct PreparedRow
    {
        // Indexed like the output files; only the entry matching each file's format is filled.
        std::vector<std::string> text;              // A complete text line, for text files.
        std::vector<std::vector<double>> float64;   // The values, for float64 NPY files.
        std::vector<std::vector<float>> float32;    // The values, for float32 NPY files.
    };

    SanityOutput(const SanityOptions &options, const std::vector<std::string> &gene_names_g,
                 const std::vector<std::string> &cell_names_c, const std::vector<double> &variance_grid_b);

    PreparedRow prepare_row(const RowComputation &result, const std::string &gene_name) const;
    void write_row(const PreparedRow &row);

    // Finishes every file; errors surface here as exceptions rather than in the destructors.
    void close();

private:
    enum class FileFormat { Text, Float64Npy, Float32Npy };

    struct FieldFile
    {
        OutputField field;
        FileFormat format;
        std::unique_ptr<sanity::Writer> text;
        std::unique_ptr<npy::NpyStreamWriter> npy;
    };

    void add_field(OutputField field, FileFormat format, const std::string &text_suffix, std::size_t num_columns);

    std::string out_folder_;
    std::vector<FieldFile> files_;
};

#endif
