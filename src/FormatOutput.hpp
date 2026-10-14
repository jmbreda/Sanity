#ifndef SANITY_FORMAT_OUTPUT_HPP
#define SANITY_FORMAT_OUTPUT_HPP

#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>
#include <system_error>

namespace sanity
{

// One formatted line (one gene) for each output file.
// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
struct FormattedOutput
{
    std::string ltq, ltq_error;                          // mu_g + delta_gc and its error bar sqrt((delta mu_g)^2 + epsilon_gc^2) (SI Eqs. 42, 45, 47).
    std::string log_fold_change, log_fold_change_error;  // delta_gc and its error bar epsilon_gc (SI Eqs. 39, 42).
    std::string mean_ltq, mean_ltq_error;                // mu_g and its error bar delta mu_g (SI Eqs. 45, 47).
    std::string gene_variance;                           // v_g (SI Eq. 17).
    std::string variance_likelihood;                     // P(v_b | n_g) on the grid of variances (SI Eq. 34).
};

inline void append_fixed_number(std::string& output, double value)
{
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite output value");
    char buffer[512];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value,
                                      std::chars_format::fixed, 6);
    if (result.ec != std::errc{}) throw std::runtime_error("Could not format output value");
    output.append(buffer, result.ptr);
}

// Called by inference workers; shared file streams are used only by the caller's
// ordered write section. Each worker holds at most one formatted gene row.
template <typename Row>
FormattedOutput format_output_row(const Row& result, const std::string& gene,
                                  bool extended)
{
    FormattedOutput row;
    const std::size_t cells = result.log_fold_change_c.size();
    row.ltq.reserve(cells * 12 + gene.size() + 1);
    row.ltq_error.reserve(cells * 12 + gene.size() + 1);
    row.ltq = gene;
    row.ltq_error = gene;
    if (extended)
    {
        row.log_fold_change.reserve(cells * 12 + 1);
        row.log_fold_change_error.reserve(cells * 12 + 1);
    }
    for (std::size_t c = 0; c < cells; ++c)
    {
        row.ltq.push_back('\t');
        append_fixed_number(row.ltq, result.mean_ltq + result.log_fold_change_c[c]);
        row.ltq_error.push_back('\t');
        append_fixed_number(row.ltq_error, std::sqrt(result.mean_ltq_error_sq + result.log_fold_change_error_sq_c[c]));
        if (extended)
        {
            if (c) { row.log_fold_change.push_back('\t'); row.log_fold_change_error.push_back('\t'); }
            append_fixed_number(row.log_fold_change, result.log_fold_change_c[c]);
            append_fixed_number(row.log_fold_change_error, std::sqrt(result.log_fold_change_error_sq_c[c]));
        }
    }
    row.ltq.push_back('\n');
    row.ltq_error.push_back('\n');
    if (extended)
    {
        row.log_fold_change.push_back('\n');
        row.log_fold_change_error.push_back('\n');
        append_fixed_number(row.mean_ltq, result.mean_ltq); row.mean_ltq.push_back('\n');
        append_fixed_number(row.mean_ltq_error, std::sqrt(result.mean_ltq_error_sq)); row.mean_ltq_error.push_back('\n');
        append_fixed_number(row.gene_variance, result.gene_variance); row.gene_variance.push_back('\n');
        row.variance_likelihood = gene;
        for (double value : result.variance_likelihood_b)
        {
            row.variance_likelihood.push_back('\t');
            append_fixed_number(row.variance_likelihood, value);
        }
        row.variance_likelihood.push_back('\n');
    }
    return row;
}

} // namespace sanity

#endif
