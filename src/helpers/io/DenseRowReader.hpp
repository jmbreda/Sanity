#pragma once

#include "helpers/io/ReadInputFiles.h"
#include <memory>

namespace sanity {

// Each worker owns a reader. The index vectors must outlive it.
// Reads into caller-owned storage, preserving input order for duplicate entries.
class DenseRowReader {
public:
    DenseRowReader(const std::string& filename, const std::string& format,
                   const std::vector<RowBlock>& blocks,
                   const std::vector<std::streampos>& offsets, int cells);
    ~DenseRowReader();
    DenseRowReader(const DenseRowReader&) = delete;
    DenseRowReader& operator=(const DenseRowReader&) = delete;

    // MTX rows are cleared before filling; TSV overwrites every element.
    void read_into(int gene, double* destination);

private:
    struct Impl;
    std::unique_ptr<Impl> implementation;
};

} // namespace sanity
