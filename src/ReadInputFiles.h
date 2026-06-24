#ifndef _ReadInputFiles_h_
#define _ReadInputFiles_h_

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <cstring>
#include <zlib.h>

using namespace std;

struct RowBlock
{
    long long offset;
    std::size_t nnz;
    RowBlock() : offset(-1), nnz(0) {}
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
    size_t max_line_length;
    char *line_buffer;

public:
    FileReader() : gz_fp(nullptr), plain_fp(nullptr), is_gzipped(false), is_open_flag(false), max_line_length(0), line_buffer(nullptr) {}

    explicit FileReader(const std::string &fname) : gz_fp(nullptr), plain_fp(nullptr), is_gzipped(false), is_open_flag(false), max_line_length(0), line_buffer(nullptr)
    {
        open(fname);
    }

    ~FileReader()
    {
        close();
        delete[] line_buffer;
    }

    bool open(const std::string &fname)
    {
        close();
        filename = fname;
        // Check if file ends with .gz
        is_gzipped = (fname.size() >= 3 && fname.substr(fname.size() - 3) == ".gz");

        if (is_gzipped)
        {
            gz_fp = gzopen(fname.c_str(), "rb");
            is_open_flag = (gz_fp != nullptr);
        }
        else
        {
            plain_fp = fopen(fname.c_str(), "rb");
            is_open_flag = (plain_fp != nullptr);
            if (is_open_flag)
            {
                // Allocate a reusable buffer to reduce system calls on plain files
                const std::size_t buf_size = 1u << 20; // 1 MB
                plain_buffer.resize(buf_size);
                if (!plain_buffer.empty())
                {
                    setvbuf(plain_fp, plain_buffer.data(), _IOFBF, plain_buffer.size());
                }
            }
        }
        return is_open_flag;
    }

    bool is_open() const
    {
        return is_open_flag;
    }

    char *getline()
    {
        if (!is_open_flag)
            return nullptr;

        if (max_line_length == 0)
        {
            max_line_length = 1024;
            line_buffer = new char[max_line_length];
        }

        size_t current_pos = 0;

        while (true)
        {
            char *res = nullptr;
            size_t remaining = max_line_length - current_pos;

            if (is_gzipped)
            {
                res = gzgets(gz_fp, line_buffer + current_pos, static_cast<int>(remaining));
            }
            else
            {
                res = fgets(line_buffer + current_pos, static_cast<int>(remaining), plain_fp);
            }

            if (res == nullptr)
            {
                if (current_pos > 0)
                {
                    // We have partial data but hit EOF/error
                    line_buffer[current_pos] = '\0';
                    break;
                }
                return nullptr;
            }

            std::size_t chunk_len = std::strlen(line_buffer + current_pos);
            if (chunk_len == 0)
            {
                if (current_pos > 0)
                {
                    break;
                }
                return nullptr;
            }

            current_pos += chunk_len;
            bool has_newline = (line_buffer[current_pos - 1] == '\n');

            if (has_newline)
            {
                line_buffer[current_pos - 1] = '\0';
                current_pos -= 1;
                break;
            }
            else if (chunk_len == remaining - 1)
            {
                // Buffer is full but no newline yet - need to expand
                size_t new_length = max_line_length * 2;
                char *new_buffer = new char[new_length];
                std::memcpy(new_buffer, line_buffer, current_pos);
                delete[] line_buffer;
                line_buffer = new_buffer;
                max_line_length = new_length;
                continue; // Read more of the line
            }
            else
            {
                // Read less than buffer size but no newline - must be EOF
                break;
            }
        }

        // Trim trailing carriage returns if present
        while (current_pos > 0 && line_buffer[current_pos - 1] == '\r')
        {
            line_buffer[current_pos - 1] = '\0';
            --current_pos;
        }

        return line_buffer;
    }

    long long tellg()
    {
        if (!is_open_flag)
            return -1;

        if (is_gzipped)
        {
            return static_cast<long long>(gztell(gz_fp));
        }
        else
        {
            if (plain_fp == nullptr)
            {
                return -1;
            }
#ifdef _WIN32
            return static_cast<long long>(_ftelli64(plain_fp));
#else
            return static_cast<long long>(ftello(plain_fp));
#endif
        }
    }

    bool seekg(long long offset)
    {
        if (!is_open_flag)
            return false;

        if (is_gzipped)
        {
            return (gzseek(gz_fp, offset, SEEK_SET) >= 0);
        }
        else
        {
            if (plain_fp == nullptr)
            {
                return false;
            }
#ifdef _WIN32
            return (_fseeki64(plain_fp, offset, SEEK_SET) == 0);
#else
            return (fseeko(plain_fp, offset, SEEK_SET) == 0);
#endif
        }
    }

    void close()
    {
        if (is_open_flag)
        {
            if (is_gzipped && gz_fp != nullptr)
            {
                gzclose(gz_fp);
                gz_fp = nullptr;
            }
            else if (!is_gzipped && plain_fp != nullptr)
            {
                fclose(plain_fp);
                plain_fp = nullptr;
            }
            plain_buffer.clear();
            is_open_flag = false;
        }
    }

    bool gzipped() const
    {
        return is_gzipped;
    }
};

// For UMI count matrix
void Get_G_C_UMIcountMatrix(string in_file,
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
void Get_G_C_MTX(string in_file,
                 int &N_rows,
                 int &G,
                 int &C,
                 map<int, int> &gene_idx,
                 vector<RowBlock> &mtx_rows,
                 vector<double> &N_c,
                 vector<double> &n);

std::vector<std::string> Read_CellNames(const std::string &filename);

std::vector<std::string> Read_GeneNames(const std::string &filename,
                                        const std::map<int, int> &gene_idx,
                                        const int G);

#endif
