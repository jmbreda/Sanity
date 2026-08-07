#pragma once

// Self-contained NumPy ".npy" reader/writer for two-dimensional float64 arrays.
//
// This is a standalone copy of the Bonsai port's general_utils/npy.{hpp,cpp},
// carved out so it can be dropped into an unrelated project without pulling in
// the rest of that codebase.  Everything it needs lives here: a minimal
// row-major `Matrix`, the reader, the whole-matrix writer, and -- the reason
// this copy exists -- a streaming writer that appends rows one at a time.
//
// Only the little-endian, C-order, float64 (`<f8`) subset of the .npy format is
// supported.  That is enough for numeric matrices produced on the platforms
// this targets, and NumPy loads the result directly with `numpy.load`.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace npy {

/// Small row-major dense matrix, self-contained so this writer carries no
/// external linear-algebra dependency.  Element (row_index, column_index) is
/// stored at data[row_index * columns + column_index].
struct Matrix {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<double> data;

    Matrix() = default;
    Matrix(std::size_t row_count, std::size_t column_count)
        : rows(row_count), columns(column_count), data(row_count * column_count) {}

    /// Bounds-checked element access, kept simple for reviewability.
    double& operator()(std::size_t row_index, std::size_t column_index) {
        return data.at(row_index * columns + column_index);
    }
    const double& operator()(std::size_t row_index, std::size_t column_index) const {
        return data.at(row_index * columns + column_index);
    }

    /// Pointer to the first value in row row_index.
    double* row_pointer(std::size_t row_index) { return data.data() + row_index * columns; }
    const double* row_pointer(std::size_t row_index) const { return data.data() + row_index * columns; }

    /// Append one complete row.  The first appended row fixes the column count.
    void append_row(const double* row, std::size_t row_length) {
        if (columns == 0) columns = row_length;
        if (row_length != columns) throw std::runtime_error("Cannot append row with wrong column count");
        data.insert(data.end(), row, row + row_length);
        ++rows;
    }
};

/// Result of reading a NumPy .npy array.
struct NpyArray {
    Matrix matrix;
    bool fortran_order = false; ///< true if the file was stored in column-major order
};

/// Read a NumPy .npy file containing a two-dimensional float64 array.
///
/// Supported subset: format versions 1.0/2.0/3.0; little-endian or
/// native-endian float64 (`<f8` or `f8`); two-dimensional arrays; C-order and
/// Fortran-order storage.  The returned Matrix is always row-major.
NpyArray read_npy_float64_2d(const std::filesystem::path& path);

/// Write a two-dimensional float64 array as a NumPy .npy file.
///
/// Produces format version 1.0, C-order (row-major), little-endian `<f8`
/// layout that `read_npy_float64_2d` reads back and that NumPy loads directly.
void write_npy_float64_2d(const Matrix& matrix, const std::filesystem::path& path);

/// Streaming NumPy .npy writer: append rows one at a time.
///
/// The .npy format records the array shape in a header at the very start of the
/// file, so an ordinary writer needs the whole matrix (and thus the final row
/// count) up front.  This writer instead reserves a fixed-size header slot,
/// streams each row straight to disk as it arrives, and -- on `close()` (or in
/// the destructor) -- seeks back and rewrites the header with the true row
/// count.  The reserved slot is sized for the widest possible count, so the
/// data block never has to move.  This is the same technique NumPy's own
/// incremental writers use.
///
/// Typical use:
/// @code
///   npy::NpyStreamWriter writer("out.npy", columns);
///   for (...) writer.write_row(row.data(), row.size());
///   writer.close();          // or just let `writer` go out of scope
/// @endcode
///
/// The number of columns may be given up front, or left as 0 and inferred from
/// the first row.  Every subsequent row must have the same length.
class NpyStreamWriter {
public:
    /// Open `path` for writing.  If `columns` is 0 the column count is taken
    /// from the first row written.  The header is written lazily (on the first
    /// row, or on close for an empty array), so no bytes hit disk until then.
    explicit NpyStreamWriter(const std::filesystem::path& path, std::size_t columns = 0);

    /// Append one row of exactly `columns()` values (contiguous, row-major).
    void write_row(const double* row, std::size_t length);

    /// Convenience overload for a contiguous vector.
    void write_row(const std::vector<double>& row) { write_row(row.data(), row.size()); }

    /// Finalise the file: patch the header with the true row count and flush.
    /// Idempotent; safe to call explicitly before the object is destroyed so
    /// that write errors surface as exceptions rather than in the destructor.
    void close();

    std::size_t rows_written() const { return rows_; }
    std::size_t columns() const { return columns_; }

    ~NpyStreamWriter();

    // Non-copyable (owns a file handle); movable is not needed here.
    NpyStreamWriter(const NpyStreamWriter&) = delete;
    NpyStreamWriter& operator=(const NpyStreamWriter&) = delete;

private:
    void write_header();       ///< reserve and write the placeholder header
    void ensure_header_written(std::size_t row_length);

    std::filesystem::path path_;
    std::ofstream out_;
    std::size_t columns_ = 0;
    std::size_t rows_ = 0;
    std::size_t header_length_ = 0;   ///< reserved header slot size, in bytes
    bool header_written_ = false;
    bool closed_ = false;
};

} // namespace npy
