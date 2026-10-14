#include "helpers/io/DenseRowReader.hpp"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::runtime_error("Usage: read_rows mode matrix format output threads");
        const std::string mode=argv[1], filename=argv[2], format=argv[3];
        const int threads=std::stoi(argv[5]);
        int rows=0, genes=0, cells=0;
        std::map<int,int> index;
        std::vector<RowBlock> blocks;
        std::vector<std::streampos> offsets;
        std::vector<double> cell_totals_c, gene_totals_g;
        std::vector<std::string> gene_ids, cell_ids;
        if (mode.substr(0,3)=="raw") {
            genes=1;cells=2;blocks.resize(1);
            blocks[0].offset=0;blocks[0].row_index=0;
            blocks[0].nnz=mode=="raw_truncated" ? 2 : 1;
        } else if (format=="mtx") index_mtx_counts(filename,rows,genes,cells,index,blocks,cell_totals_c,gene_totals_g,threads);
        else index_tsv_counts(filename,rows,genes,cells,offsets,cell_totals_c,gene_totals_g,cell_ids,gene_ids,threads);
        if (mode=="stale") blocks.at(0).row_index+=1;
        std::vector<double> values(static_cast<size_t>(genes)*cells);
        if (mode=="parallel" || mode=="raw_parallel") {
            std::vector<std::unique_ptr<sanity::DenseRowReader>> readers;
            for(int i=0;i<threads;++i) readers.push_back(std::make_unique<sanity::DenseRowReader>(filename,format,blocks,offsets,cells));
            std::vector<std::exception_ptr> failures(threads);
#pragma omp parallel num_threads(threads)
            {
                const int worker=omp_get_thread_num();
#pragma omp for schedule(static)
                for(int gene=0;gene<genes;++gene) {
                    try { readers[worker]->read_into(gene,values.data()+static_cast<size_t>(gene)*cells); }
                    catch (...) { failures[worker]=std::current_exception(); }
                }
            }
            for(const auto& failure:failures) if(failure) std::rethrow_exception(failure);
        } else {
            sanity::DenseRowReader reader(filename,format,blocks,offsets,cells);
            for(int position=0;position<genes;++position) {
                const int gene=mode=="reverse" ? genes-position-1 : position;
                reader.read_into(gene,values.data()+static_cast<size_t>(position)*cells);
            }
        }
        std::ofstream output(argv[4],std::ios::binary);
        output.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(double));
        output.close();
        if(!output) throw std::runtime_error("Cannot write row test output");
    } catch(const std::exception& failure) {
        std::cerr<<failure.what()<<'\n';
        return 1;
    }
}
