#include "npy_writer.hpp"

#include <algorithm>
#include <cstring>
#include <istream>
#include <limits>
#include <sstream>
#include <string>

namespace npy {

namespace {

// --- Header parsing helpers (reader) ---------------------------------------

// NumPy stores header lengths as little-endian integers.  These helpers read
// those integer fields explicitly rather than relying on the machine
// endianness.
std::uint16_t read_uint16_little_endian(std::istream& in) {
    unsigned char bytes[2];
    in.read(reinterpret_cast<char*>(bytes), 2);
    if (!in) throw std::runtime_error("Could not read .npy header length");
    return static_cast<std::uint16_t>(bytes[0]) | (static_cast<std::uint16_t>(bytes[1]) << 8);
}

std::uint32_t read_uint32_little_endian(std::istream& in) {
    unsigned char bytes[4];
    in.read(reinterpret_cast<char*>(bytes), 4);
    if (!in) throw std::runtime_error("Could not read .npy header length");
    return static_cast<std::uint32_t>(bytes[0]) |
           (static_cast<std::uint32_t>(bytes[1]) << 8) |
           (static_cast<std::uint32_t>(bytes[2]) << 16) |
           (static_cast<std::uint32_t>(bytes[3]) << 24);
}

std::string trim(std::string s) {
    const auto first = s.find_first_not_of(" \t\n\r");
    const auto last = s.find_last_not_of(" \t\n\r");
    if (first == std::string::npos) return {};
    return s.substr(first, last - first + 1);
}

// The .npy header is a small Python-literal dictionary.  We parse only the
// fields we need rather than implementing a general Python literal parser.
bool header_contains_true_fortran(const std::string& header) {
    const auto field_position = header.find("'fortran_order'");
    if (field_position == std::string::npos) return false;
    const auto true_position = header.find("True", field_position);
    const auto false_position = header.find("False", field_position);
    return true_position != std::string::npos &&
           (false_position == std::string::npos || true_position < false_position);
}

// Extract the two shape entries from a header fragment such as
//     'shape': (65, 1131)
// One- and higher-dimensional arrays are rejected: this reader is for matrices.
std::pair<std::size_t, std::size_t> parse_shape_2d(const std::string& header) {
    const auto shape_key_position = header.find("'shape'");
    if (shape_key_position == std::string::npos) throw std::runtime_error(".npy header lacks shape");
    const auto open = header.find('(', shape_key_position);
    const auto close = header.find(')', open);
    if (open == std::string::npos || close == std::string::npos) {
        throw std::runtime_error("Could not parse .npy shape");
    }
    const auto shape = header.substr(open + 1, close - open - 1);
    const auto comma = shape.find(',');
    if (comma == std::string::npos) throw std::runtime_error("Expected a 2D .npy shape");

    const auto rows_s = trim(shape.substr(0, comma));
    auto cols_s = trim(shape.substr(comma + 1));
    const auto trailing_comma = cols_s.find(',');
    if (trailing_comma != std::string::npos) {
        cols_s = trim(cols_s.substr(0, trailing_comma));
    }
    if (rows_s.empty() || cols_s.empty()) throw std::runtime_error("Could not parse .npy shape values");
    return {static_cast<std::size_t>(std::stoull(rows_s)), static_cast<std::size_t>(std::stoull(cols_s))};
}

// Supporting only float64 keeps this short and prevents silent conversion.
void validate_float64_descr(const std::string& header) {
    const bool little_f8 = header.find("'<f8'") != std::string::npos || header.find("\"<f8\"") != std::string::npos;
    const bool native_f8 = header.find("'f8'") != std::string::npos || header.find("\"f8\"") != std::string::npos;
    if (!little_f8 && !native_f8) {
        throw std::runtime_error("Only little-endian float64 .npy files are supported");
    }
}

// --- Header building helpers (writers) -------------------------------------

// The Python-literal dictionary describing a C-order little-endian float64
// matrix of the given shape, without any padding or trailing newline.
std::string shape_dict(std::size_t rows, std::size_t columns) {
    std::ostringstream dict;
    dict << "{'descr': '<f8', 'fortran_order': False, 'shape': ("
         << rows << ", " << columns << "), }";
    return dict.str();
}

// The .npy preamble preceding the header dictionary for format version 1.0:
// 6-byte magic, 2-byte version, 2-byte uint16 header-length field.
constexpr std::size_t kPreamble = 6 + 2 + 2;

// Pad `dict` with spaces and a trailing newline so that the whole header block
// (preamble + header) is a multiple of 64 bytes, as NumPy requires.  If
// `fixed_length` is non-zero the header is padded to exactly that many bytes
// instead (used by the streaming writer, whose slot size is decided up front).
std::string pad_header(std::string dict, std::size_t fixed_length = 0) {
    std::size_t header_length = fixed_length;
    if (header_length == 0) {
        const std::size_t unpadded_total = kPreamble + dict.size() + 1; // +1 for newline
        const std::size_t padding = (64 - (unpadded_total % 64)) % 64;
        header_length = dict.size() + padding + 1;
    }
    if (header_length > 0xFFFF) {
        throw std::runtime_error(".npy header too large for format version 1.0");
    }
    if (dict.size() + 1 > header_length) {
        throw std::runtime_error(".npy header does not fit its reserved slot");
    }
    dict.append(header_length - dict.size() - 1, ' ');
    dict.push_back('\n');
    return dict;
}

// Write the magic, version, little-endian uint16 length field, then `header`.
void write_header_block(std::ostream& out, const std::string& header) {
    out.write("\x93NUMPY", 6);
    const unsigned char version[2] = {1, 0};
    out.write(reinterpret_cast<const char*>(version), 2);
    const auto header_len = static_cast<std::uint16_t>(header.size());
    const unsigned char length_bytes[2] = {
        static_cast<unsigned char>(header_len & 0xFF),
        static_cast<unsigned char>((header_len >> 8) & 0xFF),
    };
    out.write(reinterpret_cast<const char*>(length_bytes), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
}

} // namespace

// --- Reader -----------------------------------------------------------------

NpyArray read_npy_float64_2d(const std::filesystem::path& path) {
    // File layout, simplified:
    //   magic string, version, header length, header dictionary, raw data.
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Could not open .npy file: " + path.string());

    char magic[6];
    in.read(magic, 6);
    if (!in || std::memcmp(magic, "\x93NUMPY", 6) != 0) {
        throw std::runtime_error("Not a NumPy .npy file: " + path.string());
    }

    unsigned char version[2];
    in.read(reinterpret_cast<char*>(version), 2);
    if (!in) throw std::runtime_error("Could not read .npy version");

    std::size_t header_len = 0;
    if (version[0] == 1) {
        header_len = read_uint16_little_endian(in);
    } else if (version[0] == 2 || version[0] == 3) {
        header_len = read_uint32_little_endian(in);
    } else {
        throw std::runtime_error("Unsupported .npy version");
    }

    std::string header(header_len, '\0');
    in.read(header.data(), static_cast<std::streamsize>(header.size()));
    if (!in) throw std::runtime_error("Could not read .npy header");

    validate_float64_descr(header);
    const bool fortran_order = header_contains_true_fortran(header);
    const auto [rows, columns] = parse_shape_2d(header);

    Matrix matrix(rows, columns);
    std::vector<double> raw(rows * columns);
    in.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size() * sizeof(double)));
    if (!in) throw std::runtime_error("Could not read .npy data block");

    if (!fortran_order) {
        // C-order .npy storage already matches Matrix's row-major layout.
        matrix.data = std::move(raw);
    } else {
        // Fortran-order files are column-major; convert into row-major layout.
        for (std::size_t c = 0; c < columns; ++c) {
            for (std::size_t r = 0; r < rows; ++r) {
                matrix(r, c) = raw[c * rows + r];
            }
        }
    }

    return NpyArray{std::move(matrix), fortran_order};
}

// --- Whole-matrix writer ----------------------------------------------------

void write_npy_float64_2d(const Matrix& matrix, const std::filesystem::path& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open .npy file for writing: " + path.string());

    const std::string header = pad_header(shape_dict(matrix.rows, matrix.columns));
    write_header_block(out, header);

    // Matrix is row-major, which already matches C-order .npy storage.  The
    // doubles are written in the host's native byte order; on the little-endian
    // platforms this targets that matches the declared '<f8'.
    if (!matrix.data.empty()) {
        out.write(reinterpret_cast<const char*>(matrix.data.data()),
                  static_cast<std::streamsize>(matrix.data.size() * sizeof(double)));
    }
    if (!out) throw std::runtime_error("Could not write .npy data block: " + path.string());
}

// --- Streaming writer -------------------------------------------------------

NpyStreamWriter::NpyStreamWriter(const std::filesystem::path& path, std::size_t columns)
    : path_(path), out_(path, std::ios::binary), columns_(columns) {
    if (!out_) throw std::runtime_error("Could not open .npy file for writing: " + path_.string());
}

void NpyStreamWriter::write_header() {
    // Reserve a header slot wide enough for the largest row count the shape
    // field could ever need, so the final rewrite always fits without moving
    // the data block.  std::size_t's maximum has at most 20 decimal digits.
    constexpr std::size_t kMaxRowDigits = std::numeric_limits<std::size_t>::digits10 + 1;
    std::string widest = shape_dict(0, columns_);
    widest.append(kMaxRowDigits, '0'); // stand in for the widest possible row count
    const std::size_t unpadded_total = kPreamble + widest.size() + 1;
    const std::size_t padding = (64 - (unpadded_total % 64)) % 64;
    header_length_ = widest.size() + padding + 1;

    const std::string header = pad_header(shape_dict(rows_, columns_), header_length_);
    write_header_block(out_, header);
    if (!out_) throw std::runtime_error("Could not write .npy header: " + path_.string());
    header_written_ = true;
}

void NpyStreamWriter::ensure_header_written(std::size_t row_length) {
    if (header_written_) return;
    if (columns_ == 0) columns_ = row_length; // infer from the first row
    write_header();
}

void NpyStreamWriter::write_row(const double* row, std::size_t length) {
    if (closed_) throw std::runtime_error("Cannot write to a closed .npy stream: " + path_.string());
    ensure_header_written(length);
    if (length != columns_) {
        throw std::runtime_error("Row length does not match the .npy column count");
    }
    out_.write(reinterpret_cast<const char*>(row), static_cast<std::streamsize>(length * sizeof(double)));
    if (!out_) throw std::runtime_error("Could not write .npy row: " + path_.string());
    ++rows_;
}

void NpyStreamWriter::close() {
    if (closed_) return;
    closed_ = true;

    // An empty stream still needs a valid (0, columns) header on disk.
    if (!header_written_) write_header();

    // Seek back over the reserved slot and rewrite the header with the true row
    // count.  The slot size is fixed, so the length field is unchanged and only
    // the dictionary text (the row count) differs.
    out_.seekp(static_cast<std::streamoff>(kPreamble));
    if (!out_) throw std::runtime_error("Could not seek to rewrite .npy header: " + path_.string());
    const std::string header = pad_header(shape_dict(rows_, columns_), header_length_);
    out_.write(header.data(), static_cast<std::streamsize>(header.size()));
    out_.flush();
    if (!out_) throw std::runtime_error("Could not finalise .npy header: " + path_.string());
}

NpyStreamWriter::~NpyStreamWriter() {
    // Best-effort finalisation.  Callers that need to observe write errors
    // should call close() explicitly; a throwing destructor is not safe.
    try {
        close();
    } catch (...) {
        // Swallowed deliberately; the explicit close() path reports errors.
    }
}

} // namespace npy
