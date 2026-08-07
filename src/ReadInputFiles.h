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
                   // mtx_rows because zero-count genes are dropped when that vector is built
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

// For UMI count matrix
void Get_G_C_UMIcountMatrix(std::string in_file,
                            int &N_rows,
                            int &G,
                            int &C,
                            std::vector<std::streampos> &tsv_offsets,
                            std::vector<double> &N_c,
                            std::vector<double> &n,
                            std::vector<std::string> &cell_names,
                            std::vector<std::string> &gene_names,
                            int N_threads);

// For mtx file
void Get_G_C_MTX(std::string in_file,
                 int &N_rows,
                 int &G,
                 int &C,
                 std::map<int, int> &gene_idx,
                 std::vector<RowBlock> &mtx_rows,
                 std::vector<double> &N_c,
                 std::vector<double> &n);

std::vector<std::string> Read_CellNames(const std::string &filename);

std::vector<std::string> Read_GeneNames(const std::string &filename,
                                        const std::map<int, int> &gene_idx,
                                        const int G);

#endif
