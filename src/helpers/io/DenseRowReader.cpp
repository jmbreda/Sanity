#include "helpers/io/DenseRowReader.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string_view>

namespace sanity {
namespace {
constexpr std::size_t input_buffer_size = 1024 * 1024;
constexpr std::size_t maximum_mtx_line_size = 8 * 1024 * 1024;

bool whitespace(char value) { return value == ' ' || value == '\t' || value == '\r'; }

std::string_view field(std::string_view& line) {
    while (!line.empty() && whitespace(line.front())) line.remove_prefix(1);
    const auto length = line.find_first_of(" \t\r");
    const auto result = line.substr(0, length);
    line.remove_prefix(result.size());
    if (result.empty()) throw std::runtime_error("Missing MTX field");
    return result;
}

int index_value(std::string_view token) {
    if (!token.empty() && token.front() == '+') token.remove_prefix(1);
    int value = 0;
    const auto parsed = std::from_chars(token.data(), token.data()+token.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != token.data()+token.size() || value <= 0)
        throw std::runtime_error("Invalid MTX index");
    return value-1;
}

double count_value(std::string_view token) {
    if (!token.empty() && token.front() == '+') token.remove_prefix(1);
    std::uint64_t integer = 0;
    const auto integral = std::from_chars(token.data(), token.data()+token.size(), integer);
    double value = 0;
    if (integral.ec == std::errc{} && integral.ptr == token.data()+token.size() &&
        integer <= (std::uint64_t{1} << 53)) {
        value = static_cast<double>(integer);
    } else {
        const auto parsed = std::from_chars(token.data(), token.data()+token.size(), value);
        if (parsed.ec != std::errc{} || parsed.ptr != token.data()+token.size())
            throw std::runtime_error("Invalid MTX count");
    }
    if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid count");
    return value;
}

class BufferedLines {
    std::ifstream stream;
    std::vector<char> buffer;
    std::string continuation;
    long long beginning = 0;
    std::size_t position = 0, available = 0;
public:
    void open(const std::string& filename) {
        stream.open(filename, std::ios::binary);
        if (!stream) throw std::runtime_error("Cannot open row input");
        buffer.resize(input_buffer_size);
    }
    void seek(long long offset) {
        if (offset < 0) throw std::runtime_error("Invalid row offset");
        if (offset >= beginning && offset < beginning + static_cast<long long>(available)) {
            position = static_cast<std::size_t>(offset-beginning);
            return;
        }
        stream.clear();
        stream.seekg(offset);
        if (!stream) throw std::runtime_error("Cannot seek to input row");
        beginning = offset;
        position = available = 0;
    }
    std::string_view line() {
        continuation.clear();
        for (;;) {
            if (position == available) {
                beginning += static_cast<long long>(available);
                stream.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                available = static_cast<std::size_t>(stream.gcount());
                position = 0;
                if (stream.bad()) throw std::runtime_error("Error reading row input");
                if (!available) {
                    if (!continuation.empty()) return continuation;
                    throw std::runtime_error("Unexpected end of row input");
                }
            }
            const char* cursor = buffer.data()+position;
            const char* newline = static_cast<const char*>(std::memchr(cursor, '\n', available-position));
            const std::size_t length = newline ? static_cast<std::size_t>(newline-cursor) : available-position;
            position += length + (newline ? 1 : 0);
            if (continuation.size()+length > maximum_mtx_line_size)
                throw std::runtime_error("MTX line exceeds bounded row reader limit");
            if (newline && continuation.empty()) return std::string_view(cursor, length);
            continuation.append(cursor, length);
            if (newline) return continuation;
        }
    }
};
} // namespace

struct DenseRowReader::Impl {
    std::string filename, format;
    const std::vector<RowBlock>& blocks;
    const std::vector<std::streampos>& offsets;
    int cells;
    bool gzip, opened = false;
    BufferedLines buffered;
    FileReader existing;

    Impl(const std::string& name, const std::string& kind,
         const std::vector<RowBlock>& rows, const std::vector<std::streampos>& positions, int columns)
        : filename(name), format(kind), blocks(rows), offsets(positions), cells(columns),
          gzip(name.size() >= 3 && name.substr(name.size()-3) == ".gz") {
        if (cells <= 0) throw std::runtime_error("Invalid row width");
    }
    void open() {
        if (opened) return;
        if (format == "mtx" && !gzip) buffered.open(filename);
        else if (!existing.open(filename)) throw std::runtime_error("Cannot open row input");
        opened = true;
    }
    void read(int gene, double* destination) {
        if (!destination || gene < 0) throw std::runtime_error("Invalid row destination");
        open();
        if (format == "mtx") {
            const auto& block = blocks.at(static_cast<std::size_t>(gene));
            std::fill(destination, destination+cells, 0.0);
            if (gzip) {
                if (!existing.seekg(block.offset)) throw std::runtime_error("Cannot seek to MTX row");
            } else buffered.seek(block.offset);
            std::size_t entries = 0;
            while (entries < block.nnz) {
                std::string_view line;
                if (gzip) {
                    const char* text = existing.getline();
                    if (!text) throw std::runtime_error("Unexpected end of MTX row");
                    line = text;
                } else line = buffered.line();
                while (!line.empty() && whitespace(line.front())) line.remove_prefix(1);
                if (line.empty() || line.front() == '%') continue;
                const int row = index_value(field(line));
                const int column = index_value(field(line));
                if (row != block.row_index) throw std::runtime_error("MTX row does not match indexed gene");
                if (column >= cells) throw std::runtime_error("MTX column exceeds matrix dimensions");
                const auto token = field(line);
                // Preserve the existing conversion for compressed input.
                const double value = gzip ? std::stod(std::string(token)) : count_value(token);
                if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid count");
                while (!line.empty() && whitespace(line.front())) line.remove_prefix(1);
                if (!line.empty()) throw std::runtime_error("Unexpected MTX entry fields");
                destination[column] = value;
                ++entries;
            }
        } else {
            if (!existing.seekg(static_cast<long long>(offsets.at(static_cast<std::size_t>(gene)))))
                throw std::runtime_error("Cannot seek to TSV row");
            char* line = existing.getline();
            if (!line) throw std::runtime_error("Unexpected end of TSV row");
            char* state = nullptr;
            if (!strtok_r(line, "\t\r\n", &state)) throw std::runtime_error("Missing TSV gene ID");
            for (int column = 0; column < cells; ++column) {
                const char* token = strtok_r(nullptr, "\t\r\n", &state);
                if (!token) throw std::runtime_error("Missing TSV count");
                const double value = std::stod(token);
                if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid count");
                destination[column] = value;
            }
        }
    }
};

DenseRowReader::DenseRowReader(const std::string& filename, const std::string& format,
                             const std::vector<RowBlock>& blocks,
                             const std::vector<std::streampos>& offsets, int cells)
    : implementation(std::make_unique<Impl>(filename, format, blocks, offsets, cells)) {}
DenseRowReader::~DenseRowReader() = default;
void DenseRowReader::read_into(int gene, double* destination) {
    implementation->read(gene, destination);
}
} // namespace sanity
