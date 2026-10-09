#ifndef SANITY_FORMAT_OUTPUT_HPP
#define SANITY_FORMAT_OUTPUT_HPP

#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>
#include <system_error>

namespace sanity
{

enum class TextOutputMode { Standard, Extended, Bonsai };

struct FormattedOutput
{
    std::string ltq, ltq_error, delta, delta_error;
    std::string mu, mu_error, variance, likelihood;
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
                                  TextOutputMode mode)
{
    const bool extended = mode != TextOutputMode::Standard;
    const bool bonsai = mode == TextOutputMode::Bonsai;
    FormattedOutput row;
    const std::size_t cells = result.delta.size();
    if (!bonsai)
    {
        row.ltq.reserve(cells * 12 + gene.size() + 1);
        row.ltq_error.reserve(cells * 12 + gene.size() + 1);
        row.ltq = gene;
        row.ltq_error = gene;
    }
    if (extended)
    {
        row.delta.reserve(cells * 12 + 1);
        row.delta_error.reserve(cells * 12 + 1);
    }
    for (std::size_t c = 0; c < cells; ++c)
    {
        if (!bonsai)
        {
            row.ltq.push_back('\t');
            append_fixed_number(row.ltq, result.mu + result.delta[c]);
            row.ltq_error.push_back('\t');
            append_fixed_number(row.ltq_error, std::sqrt(result.var_mu + result.var_delta[c]));
        }
        if (extended)
        {
            if (c) { row.delta.push_back('\t'); row.delta_error.push_back('\t'); }
            append_fixed_number(row.delta, result.delta[c]);
            append_fixed_number(row.delta_error, std::sqrt(result.var_delta[c]));
        }
    }
    if (!bonsai) { row.ltq.push_back('\n'); row.ltq_error.push_back('\n'); }
    if (extended)
    {
        row.delta.push_back('\n');
        row.delta_error.push_back('\n');
        append_fixed_number(row.mu, result.mu); row.mu.push_back('\n');
        if (!bonsai) { append_fixed_number(row.mu_error, std::sqrt(result.var_mu)); row.mu_error.push_back('\n'); }
        append_fixed_number(row.variance, result.var_gene); row.variance.push_back('\n');
        if (!bonsai)
        {
            row.likelihood = gene;
            for (double value : result.lik)
            {
                row.likelihood.push_back('\t');
                append_fixed_number(row.likelihood, value);
            }
            row.likelihood.push_back('\n');
        }
    }
    return row;
}

} // namespace sanity

#endif
