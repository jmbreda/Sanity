#ifndef _ReadInputFiles_h_
#define _ReadInputFiles_h_

#include <cstddef>
#include <cstdio>
#include <ios>
#include <map>
#include <string>
#include <vector>
#include <zlib.h>


struct RowBlock
{
    long long offset;
    std::size_t nnz;
    int row_index; // 0-based row index in the mtx file, which differs from the position in
                   // mtx_rows_g because zero-count genes are dropped when that vector is built
    RowBlock() : offset(-1), nnz(0), row_index(-1) {}
};

// Wrapper class for transparent gzip/plain file reading
class FileReader
{
private:
    gzFile gz_fp;
    FILE *plain_fp;
    std::vector<char> plain_buffer;
    bool is_gzipped;
    bool is_open_flag;
    std::string filename;
    std::size_t max_line_length;
    char *line_buffer;

public:
    FileReader();
    explicit FileReader(const std::string &fname);
    ~FileReader();

    bool open(const std::string &fname);
    bool is_open() const;
    char *getline();
    long long tellg();
    bool seekg(long long offset);
    void close();
    bool gzipped() const;
};

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells.
// Both index functions scan the input once and return:
//   num_input_rows: the number of genes (rows) in the input, including genes without counts.
//   num_genes: G, the number of genes with at least one count; only these are kept.
//   num_cells: C, the number of cells.
//   cell_totals_c: N_c, the total counts of each cell (SI Eq. 9).
//   gene_totals_g: n_g, the total counts of each kept gene (SI Eq. 15).

// For UMI count matrix
// tsv_offsets_g holds the file position of the row of each kept gene.
void index_tsv_counts(std::string in_file,
                      int &num_input_rows,
                      int &num_genes,
                      int &num_cells,
                      std::vector<std::streampos> &tsv_offsets_g,
                      std::vector<double> &cell_totals_c,
                      std::vector<double> &gene_totals_g,
                      std::vector<std::string> &cell_names_c,
                      std::vector<std::string> &gene_names_g,
                      int num_threads);

// For mtx file
// row_to_gene maps each (0-based) MTX row to the index of its gene among the kept genes, or -1 if the row has no counts.
// mtx_rows_g holds the position and size of the entries of each kept gene.
void index_mtx_counts(std::string in_file,
                      int &num_input_rows,
                      int &num_genes,
                      int &num_cells,
                      std::map<int, int> &row_to_gene,
                      std::vector<RowBlock> &mtx_rows_g,
                      std::vector<double> &cell_totals_c,
                      std::vector<double> &gene_totals_g,
                      int num_threads = 4);

std::vector<std::string> Read_CellNames(const std::string &filename);

std::vector<std::string> Read_GeneNames(const std::string &filename,
                                        const std::map<int, int> &row_to_gene,
                                        const int num_genes);

#endif
