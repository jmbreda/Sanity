#ifndef SANITY_FAST_MTX_READER_HPP
#define SANITY_FAST_MTX_READER_HPP

#include "ReadInputFiles.h"
#include <algorithm>
#include <atomic>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <omp.h>
#include <sstream>
#include <stdexcept>

// Internal implementation of the uncompressed, integral-count MTX reader.
namespace fast_mtx {
constexpr std::uint64_t exact_limit = std::uint64_t{1} << 53;
constexpr std::size_t buffer_size = 8 * 1024 * 1024;

struct Unsupported : std::runtime_error { using std::runtime_error::runtime_error; };

inline bool whitespace(char value) { return value == ' ' || value == '\t' || value == '\r'; }
inline void skip_space(const char*& cursor, const char* ending) {
    while (cursor < ending && whitespace(*cursor)) ++cursor;
}
inline std::uint64_t parse_integer(const char*& cursor, const char* ending) {
    skip_space(cursor, ending);
    if (cursor < ending && *cursor == '+') ++cursor;
    std::uint64_t value = 0;
    const auto parsed = std::from_chars(cursor, ending, value);
    if (parsed.ec != std::errc{} || parsed.ptr == cursor ||
        (parsed.ptr != ending && !whitespace(*parsed.ptr)))
        throw std::runtime_error("Invalid integer field in MTX input");
    cursor = parsed.ptr;
    return value;
}
inline std::uint64_t parse_count(const char*& cursor, const char* ending) {
    skip_space(cursor, ending);
    if (cursor < ending && *cursor == '+') ++cursor;
    const char* beginning = cursor;
    while (cursor < ending && !whitespace(*cursor)) ++cursor;
    std::uint64_t integral = 0;
    const auto integer_result = std::from_chars(beginning, cursor, integral);
    if (integer_result.ec == std::errc{} && integer_result.ptr == cursor && cursor != beginning) {
        if (integral > exact_limit) throw Unsupported("count exceeds exact integer range");
        return integral;
    }
    double value = 0;
    const auto real_result = std::from_chars(beginning, cursor, value);
    if (real_result.ec != std::errc{} || real_result.ptr != cursor || cursor == beginning)
        throw std::runtime_error("Invalid count in MTX input");
    if (!std::isfinite(value) || value < 0)
        throw std::runtime_error("MTX counts must be finite and nonnegative");
    if (value > double(exact_limit) || std::floor(value) != value)
        throw Unsupported("counts are outside the nonnegative integral fast path");
    return static_cast<std::uint64_t>(value);
}
inline void add_exact(std::uint64_t& total, std::uint64_t value) {
    if (value > exact_limit - total) throw Unsupported("total exceeds exact integer range");
    total += value;
}

struct Partial {
    std::vector<std::uint64_t> cell_totals, gene_totals, entries;
    std::vector<long long> offsets;
    int first_row = -1, last_row = -1;
    std::uint64_t nonzeros = 0;
};

// Move each byte partition to a line boundary, including CRLF and missing final newline.
inline long long next_line(std::ifstream& input, long long offset, long long ending) {
    if (offset >= ending) return ending;
    input.clear();
    input.seekg(offset - 1);
    char previous = 0;
    input.get(previous);
    if (!input) throw std::runtime_error("Cannot seek in MTX input");
    if (previous == '\n') return offset;
    char buffer[4096];
    while (offset < ending) {
        const auto length = std::min<long long>(sizeof(buffer), ending - offset);
        input.read(buffer, length);
        if (input.gcount() != length) throw std::runtime_error("Cannot read MTX partition boundary");
        const auto found = static_cast<const char*>(std::memchr(buffer, '\n', length));
        if (found) return offset + (found - buffer) + 1;
        offset += length;
    }
    return ending;
}

// Same outputs as index_mtx_counts (see ReadInputFiles.h).
// The exact_*_totals vectors are the integer totals before they are converted to double.
inline void index(const std::string& filename, int num_threads, int& num_input_rows, int& num_genes, int& num_cells,
                  std::map<int, int>& row_to_gene, std::vector<RowBlock>& mtx_rows_g,
                  std::vector<double>& cell_totals_c, std::vector<double>& gene_totals_g, int& actual_threads) {
    if (num_threads < 1) throw std::runtime_error("Thread count must be positive");
    std::ifstream input(filename, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open MTX input");
    input.seekg(0, std::ios::end);
    const long long file_size = static_cast<long long>(input.tellg());
    if (file_size < 0) throw std::runtime_error("Cannot determine MTX input size");
    input.seekg(0);
    std::string line, magic, object, layout, field, symmetry;
    if (!std::getline(input, line)) throw std::runtime_error("Missing MTX header");
    std::istringstream header(line);
    header >> magic >> object >> layout >> field >> symmetry;
    if (magic != "%%MatrixMarket" || object != "matrix" || layout != "coordinate" ||
        (field != "real" && field != "integer") || symmetry != "general")
        throw Unsupported("MTX header is outside the general real/integer fast path");
    std::uint64_t rows = 0, columns = 0, nonzeros = 0;
    bool dimensions_found = false;
    while (std::getline(input, line)) {
        const char* fields = line.data();
        const char* ending = fields + line.size();
        skip_space(fields, ending);
        if (fields == ending || *fields == '%') continue;
        rows = parse_integer(fields, ending);
        columns = parse_integer(fields, ending);
        nonzeros = parse_integer(fields, ending);
        skip_space(fields, ending);
        if (fields != ending) throw std::runtime_error("Unexpected MTX dimension fields");
        dimensions_found = true;
        break;
    }
    if (!dimensions_found || !rows || !columns || rows > std::uint64_t(std::numeric_limits<int>::max()) ||
        columns > std::uint64_t(std::numeric_limits<int>::max()))
        throw std::runtime_error("Invalid MTX dimensions");
    const long long body = input.eof() ? file_size : static_cast<long long>(input.tellg());
    if (body < 0) throw std::runtime_error("Cannot locate MTX entries");

    std::vector<Partial> partials;
    std::vector<std::exception_ptr> failures;
    std::vector<long long> boundaries;
    std::atomic<bool> cancelled{false};
#pragma omp parallel num_threads(num_threads)
    {
#pragma omp single
        {
            try {
                actual_threads = omp_get_num_threads();
                partials.resize(actual_threads);
                failures.resize(actual_threads);
                boundaries.resize(actual_threads + 1);
                boundaries.front() = body;
                boundaries.back() = file_size;
                const long long length = file_size - body;
                for (int part = 1; part < actual_threads; ++part) {
                    const long long offset = body + (length / actual_threads) * part +
                                             ((length % actual_threads) * part) / actual_threads;
                    boundaries[part] = offset == body ? body : next_line(input, offset, file_size);
                }
            } catch (...) {
                // No exceptions may escape an OpenMP region.
                if (failures.empty()) failures.resize(1);
                failures[0] = std::current_exception();
                cancelled.store(true);
            }
        }
        const int part = omp_get_thread_num();
        if (!cancelled.load()) {
            try {
                Partial& partial = partials[part];
                partial.cell_totals.assign(columns, 0);
                partial.gene_totals.assign(rows, 0);
                partial.entries.assign(rows, 0);
                partial.offsets.assign(rows, -1);
                std::ifstream stream(filename, std::ios::binary);
                stream.seekg(boundaries[part]);
                if (!stream) throw std::runtime_error("Cannot open MTX partition");
                std::vector<char> buffer(buffer_size);
                std::string pending;
                long long position = boundaries[part];
                auto consume = [&](const char* beginning, const char* ending, long long offset) {
                    const char* fields = beginning;
                    skip_space(fields, ending);
                    if (fields == ending || *fields == '%') return;
                    const auto row_number = parse_integer(fields, ending);
                    const auto cell_number = parse_integer(fields, ending);
                    const auto count = parse_count(fields, ending);
                    skip_space(fields, ending);
                    if (fields != ending) throw std::runtime_error("Unexpected MTX entry fields");
                    if (!row_number || row_number > rows || !cell_number || cell_number > columns)
                        throw std::runtime_error("MTX index out of range");
                    const int row = static_cast<int>(row_number - 1);
                    if (partial.last_row > row) throw std::runtime_error("MTX input is not sorted by row; use scripts/sort_mtx_by_row.py");
                    if (partial.first_row < 0) partial.first_row = row;
                    partial.last_row = row;
                    if (!partial.entries[row]) partial.offsets[row] = offset;
                    ++partial.entries[row];
                    ++partial.nonzeros;
                    add_exact(partial.gene_totals[row], count);
                    add_exact(partial.cell_totals[cell_number - 1], count);
                };
                long long pending_offset = position;
                while (position < boundaries[part + 1] && !cancelled.load()) {
                    const auto length = std::min<long long>(buffer.size(), boundaries[part + 1] - position);
                    stream.read(buffer.data(), length);
                    if (stream.gcount() != length) throw std::runtime_error("Cannot read MTX partition");
                    const char* cursor = buffer.data();
                    const char* ending = cursor + length;
                    while (cursor < ending) {
                        const auto newline = static_cast<const char*>(std::memchr(cursor, '\n', ending - cursor));
                        if (!newline) {
                            if (pending.empty()) pending_offset = position + (cursor - buffer.data());
                            pending.append(cursor, ending);
                            if (pending.size() > buffer_size) throw Unsupported("MTX line exceeds fast reader buffer");
                            break;
                        }
                        if (pending.empty()) consume(cursor, newline, position + (cursor - buffer.data()));
                        else {
                            pending.append(cursor, newline);
                            consume(pending.data(), pending.data() + pending.size(), pending_offset);
                            pending.clear();
                        }
                        cursor = newline + 1;
                    }
                    position += length;
                }
                if (!pending.empty() && !cancelled.load())
                    consume(pending.data(), pending.data() + pending.size(), pending_offset);
            } catch (...) {
                failures[part] = std::current_exception();
                cancelled.store(true);
            }
        }
    }
    for (const auto& failure : failures) if (failure) std::rethrow_exception(failure);
    std::vector<std::uint64_t> exact_gene_totals_g(rows, 0), exact_cell_totals_c(columns, 0), entries(rows, 0);
    std::vector<long long> offsets(rows, -1);
    int last_row = -1;
    std::uint64_t observed = 0;
    for (const auto& partial : partials) {
        if (partial.first_row >= 0 && last_row > partial.first_row)
            throw std::runtime_error("MTX input is not sorted across partitions; use scripts/sort_mtx_by_row.py");
        if (partial.last_row >= 0) last_row = partial.last_row;
        observed += partial.nonzeros;
        for (std::size_t row = 0; row < rows; ++row) {
            if (offsets[row] < 0 && partial.entries[row]) offsets[row] = partial.offsets[row];
            entries[row] += partial.entries[row];
            add_exact(exact_gene_totals_g[row], partial.gene_totals[row]);
        }
        for (std::size_t cell = 0; cell < columns; ++cell) add_exact(exact_cell_totals_c[cell], partial.cell_totals[cell]);
    }
    if (observed != nonzeros) throw std::runtime_error("MTX entry count differs from header");

    // Publish only a complete index; an unsupported input leaves the caller untouched.
    num_input_rows = static_cast<int>(rows);
    num_cells = static_cast<int>(columns);
    num_genes = 0;
    mtx_rows_g.clear();
    gene_totals_g.clear();
    row_to_gene.clear();
    cell_totals_c.assign(exact_cell_totals_c.begin(), exact_cell_totals_c.end());
    for (int row = 0; row < num_input_rows; ++row) {
        row_to_gene[row] = entries[row] ? num_genes : -1;
        if (entries[row]) {
            RowBlock block;
            block.offset = offsets[row];
            block.nnz = entries[row];
            block.row_index = row;
            mtx_rows_g.push_back(block);
            gene_totals_g.push_back(double(exact_gene_totals_g[row]));
            ++num_genes;
        }
    }
}
} // namespace fast_mtx
#endif
