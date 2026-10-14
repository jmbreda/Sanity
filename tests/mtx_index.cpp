#include "ReadInputFiles.h"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

template<class T>
void save(const std::filesystem::path& path, const std::vector<T>& values) {
    std::ofstream output(path, std::ios::binary);
    output.write(reinterpret_cast<const char*>(values.data()), values.size()*sizeof(T));
    output.close();
    if (!output) throw std::runtime_error("Cannot write test index");
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("Usage: mtx_index matrix threads output");
        int rows=0, genes=0, cells=0;
        std::map<int,int> index;
        std::vector<RowBlock> blocks;
        std::vector<double> cell_totals_c, gene_totals_g;
        index_mtx_counts(argv[1], rows, genes, cells, index, blocks, cell_totals_c, gene_totals_g, std::stoi(argv[2]));
        const std::filesystem::path directory(argv[3]);
        std::filesystem::create_directory(directory);
        std::vector<std::int64_t> offsets;
        std::vector<std::uint64_t> entries;
        std::vector<std::int32_t> row_indices;
        for (const auto& block : blocks) {
            offsets.push_back(block.offset);
            entries.push_back(block.nnz);
            row_indices.push_back(block.row_index);
        }
        save(directory/"offsets.bin", offsets);
        save(directory/"entries.bin", entries);
        save(directory/"row_indices.bin", row_indices);
        save(directory/"gene_totals.bin", gene_totals_g);
        save(directory/"cell_totals.bin", cell_totals_c);
        std::ofstream metadata(directory/"dimensions.txt");
        metadata<<rows<<' '<<genes<<' '<<cells<<'\n';
        metadata.close();
        if (!metadata) throw std::runtime_error("Cannot write test dimensions");
    } catch (const std::exception& failure) {
        std::cerr<<failure.what()<<'\n';
        return 1;
    }
}
