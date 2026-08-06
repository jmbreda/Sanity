#include "ReadInputFiles.h"
#include <omp.h>
#include <atomic>
#include <cstring>


/*** FileReader class implementation ***/
FileReader::FileReader() : gz_fp(nullptr), plain_fp(nullptr), is_gzipped(false), is_open_flag(false), max_line_length(0), line_buffer(nullptr) {}

FileReader::FileReader(const std::string &fname) : gz_fp(nullptr), plain_fp(nullptr), is_gzipped(false), is_open_flag(false), max_line_length(0), line_buffer(nullptr)
{
    open(fname);
}

FileReader::~FileReader()
{
    close();
    delete[] line_buffer;
}

bool FileReader::open(const std::string &fname)
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
            // Allocate a reusable buffer to reduce system calls on plain files.
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

bool FileReader::is_open() const
{
    return is_open_flag;
}

char *FileReader::getline()
{
    if (!is_open_flag)
        return nullptr;

    if (max_line_length == 0)
    {
        max_line_length = 1024;
        line_buffer = new char[max_line_length];
    }

    std::size_t current_pos = 0;

    while (true)
    {
        char *res = nullptr;
        std::size_t remaining = max_line_length - current_pos;

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
                // We have partial data but hit EOF/error.
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
            // Buffer is full but no newline yet, need to expand.
            std::size_t new_length = max_line_length * 2;
            char *new_buffer = new char[new_length];
            std::memcpy(new_buffer, line_buffer, current_pos);
            delete[] line_buffer;
            line_buffer = new_buffer;
            max_line_length = new_length;
            continue; // Read more of the line.
        }
        else
        {
            // Read less than buffer size but no newline, must be EOF.
            break;
        }
    }

    // Trim trailing carriage returns if present.
    while (current_pos > 0 && line_buffer[current_pos - 1] == '\r')
    {
        line_buffer[current_pos - 1] = '\0';
        --current_pos;
    }

    return line_buffer;
}

long long FileReader::tellg()
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

bool FileReader::seekg(long long offset)
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

void FileReader::close()
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

bool FileReader::gzipped() const
{
    return is_gzipped;
}
/*** end of FileReader class implementation ***/


void Get_G_C_UMIcountMatrix(std::string in_file,
                            int &N_rows,
                            int &G,
                            int &C,
                            std::vector<std::streampos> &tsv_offsets,
                            std::vector<double> &N_c,
                            std::vector<double> &n,
                            std::vector<std::string> &cell_names,
                            std::vector<std::string> &gene_names,
                            int N_threads)
{

    FileReader infp(in_file);
    if (!infp.is_open())
    {
        fprintf(stderr, "Cannot open input file %s\n", in_file.c_str());
        exit(EXIT_FAILURE);
    }

    // Count cell. First line should have the names of the columns (cell names)
    char *ss = infp.getline();
    if (ss == nullptr)
    {
        fprintf(stderr, "Error: Unable to read header line from %s\n", in_file.c_str());
        exit(EXIT_FAILURE);
    }
    char *saveptr = nullptr;

    char *token = strtok_r(ss, " \t,", &saveptr);
    cell_names.clear();
    if (token == NULL)
    {
        fprintf(stderr, "Error: header line missing first field in %s\n", in_file.c_str());
        exit(EXIT_FAILURE);
    }
    // remaining tokens are cell names
    token = strtok_r(NULL, " \t,", &saveptr);
    while (token)
    {
        std::string name(token);
        std::size_t pos = name.find('\r');
        if (pos != std::string::npos)
        {
            name.erase(pos, 1);
        }
        pos = name.find('\n');
        if (pos != std::string::npos)
        {
            name.erase(pos, 1);
        }
        cell_names.push_back(name);
        token = strtok_r(NULL, " \t,", &saveptr);
    }
    C = static_cast<int>(cell_names.size());

    // PASS 1: Sequential scan to record all line offsets
    std::vector<long long> line_offsets;
    N_rows = 0;

    while (true)
    {
        long long row_offset = infp.tellg();
        ss = infp.getline();
        if (ss == nullptr)
        {
            break;
        }
        if (!ss || ss[0] == '\0' || ss[0] == '\n')
        {
            continue;
        }
        else
        {
            line_offsets.push_back(row_offset);
            ++N_rows;
        }
    }
    infp.close();

    if (N_rows == 0)
    {
        G = 0;
        N_c.assign(static_cast<std::size_t>(C), 0.0);
        n.clear();
        gene_names.clear();
        tsv_offsets.clear();
        return;
    }

    // PASS 2: Parallel processing of rows
    std::vector<std::string> gene_names_tmp(N_rows);
    std::vector<double> gene_totals_tmp(N_rows, 0.0);
    std::vector<char> is_expressed(N_rows, 0);
    std::vector<double> cell_totals(static_cast<size_t>(C), 0.0);
    std::atomic<bool> parse_error(false);
    std::string parse_error_msg;
    std::string parse_error_line;
    int parse_error_row = -1;

    #pragma omp parallel num_threads(N_threads)
    {
        // Thread-local accumulators for cell totals
        std::vector<double> local_cell_totals(static_cast<size_t>(C), 0.0);

        // Each thread opens its own FileReader instance once and reuses it
        FileReader thread_reader(in_file);
        char *thread_ss = nullptr;
        char *thread_sc = nullptr;
        char *thread_token = nullptr;
        if (!thread_reader.is_open())
        {
            #pragma omp critical
            {
                fprintf(stderr, "Error: Thread cannot open file %s\n", in_file.c_str());
            }
        }

        #pragma omp for schedule(static)
        for (int row_idx = 0; row_idx < N_rows; ++row_idx)
        {
            if (parse_error.load())
            {
                continue;
            }
            if (!thread_reader.is_open())
            {
                continue;
            }
            char *thread_saveptr = nullptr;
            // Seek to the line offset
            thread_reader.seekg(line_offsets[row_idx]);
            thread_ss = thread_reader.getline();
            if (thread_ss == nullptr)
            {
                continue;
            }
            thread_sc = new char[static_cast<std::size_t>(strlen(thread_ss) + 1)];
            std::strncpy(thread_sc, thread_ss, static_cast<std::size_t>(strlen(thread_ss) + 1));
            thread_sc[strlen(thread_ss)] = '\0';

            thread_token = strtok_r(thread_sc, " \t,", &thread_saveptr);
            if (thread_token == NULL)
            {
                delete[] thread_sc;
                continue;
            }
            std::string gene_id(thread_token);
            std::size_t pos = gene_id.find('\r');
            if (pos != std::string::npos)
            {
                gene_id.erase(pos, 1);
            }
            pos = gene_id.find('\n');
            if (pos != std::string::npos)
            {
                gene_id.erase(pos, 1);
            }

            double row_sum = 0.0;
            bool row_parse_failed = false;
            for (int c = 0; c < C; ++c)
            {
                thread_token = strtok_r(NULL, " \t,", &thread_saveptr);
                if (thread_token != NULL)
                {
                    double value = std::stod(thread_token); /**total count this gene**/
                    row_sum += value;
                    local_cell_totals[static_cast<std::size_t>(c)] += value;
                }
                else
                {
                    #pragma omp critical
                    {
                        if (!parse_error.load())
                        {
                            parse_error = true;
                            parse_error_row = row_idx + 1;
                            parse_error_msg = "Error: not enough fields on line number " + std::to_string(parse_error_row) + " in " + in_file;
                            parse_error_line = thread_sc ? thread_sc : "";
                        }
                    }
                    row_parse_failed = true;
                    break;
                }
            }
            delete[] thread_sc;
            thread_sc = nullptr;
            if (row_parse_failed)
            {
                continue;
            }

            // Record gene if total count is greater than 0
            if (row_sum > 0.0)
            {
                is_expressed[row_idx] = 1;
                gene_totals_tmp[row_idx] = row_sum;
                gene_names_tmp[row_idx] = gene_id;
            }
        }

        // Close thread's FileReader
        thread_reader.close();

        // Reduce cell totals from all threads
        #pragma omp critical
        {
            for (std::size_t c = 0; c < static_cast<std::size_t>(C); ++c)
            {
                cell_totals[c] += local_cell_totals[c];
            }
        }
    }

    // PASS 3: Build final results (only expressed genes)
    if (parse_error.load())
    {
        fprintf(stderr, "%s\n", parse_error_msg.c_str());
        if (!parse_error_line.empty())
        {
            fprintf(stderr, "--%s--\n", parse_error_line.c_str());
        }
        exit(EXIT_FAILURE);
    }

    gene_names.clear();
    tsv_offsets.clear();
    std::vector<double> gene_totals;
    G = 0;

    for (int row_idx = 0; row_idx < N_rows; ++row_idx)
    {
        if (is_expressed[row_idx])
        {
            ++G;
            tsv_offsets.push_back(static_cast<std::streampos>(line_offsets[row_idx]));
            gene_totals.push_back(gene_totals_tmp[row_idx]);
            gene_names.push_back(gene_names_tmp[row_idx]);
        }
    }

    // Assign results to output vectors
    N_c.assign(cell_totals.begin(), cell_totals.end());
    n.assign(gene_totals.begin(), gene_totals.end());
}

void Get_G_C_MTX(std::string in_file, int &N_rows, int &G, int &C, std::map<int, int> &gene_idx, std::vector<RowBlock> &mtx_rows, std::vector<double> &N_c, std::vector<double> &n)
{

    FileReader infp(in_file);
    if (!infp.is_open())
    {
        fprintf(stderr, "Cannot open input file %s\n", in_file.c_str());
        exit(EXIT_FAILURE);
    }

    // Check first line for MTX format - must be "coordinate"
    char *line_buffer = nullptr;
    char *token = NULL;

    for (int i = 0; i < 100; ++i)
    {
        line_buffer = infp.getline();
        if (line_buffer == nullptr)
        {
            fprintf(stderr, "Error: unable to read from file %s\n", in_file.c_str());
            exit(EXIT_FAILURE);
        }
        else if (line_buffer[0] == '\0')
        {
            continue;
        }
        else if (i > 99)
        {
            fprintf(stderr, "Error: no valid MTX header found in first 100 lines of %s\n", in_file.c_str());
            exit(EXIT_FAILURE);
        }
        else
        {
            if (line_buffer[0] == '%')
            {
                token = strtok(line_buffer, " ");
                token = strtok(NULL, " ");
                token = strtok(NULL, " ");
                if (std::string(token) != "coordinate")
                {
                    fprintf(stderr, "Error: only MTX coordinate format is supported in %s\n", in_file.c_str());
                    exit(EXIT_FAILURE);
                }
                break;
            }
        }
    }

    // Read header: first non-comment line contains N_rows C and (ignored third value)
    while (true)
    {
        line_buffer = infp.getline();
        if (line_buffer == nullptr)
        {
            fprintf(stderr, "Error: no valid header found in %s\n", in_file.c_str());
            exit(EXIT_FAILURE);
        }

        if (line_buffer[0] == '\0' || line_buffer[0] == '\n' || line_buffer[0] == '%')
        {
            continue;
        }
        else
        {
            token = strtok(line_buffer, " ");
            // Highest number of genes expressed and not
            N_rows = std::stoi(token);
            token = strtok(NULL, " \t");
            // Number of cells
            C = std::stoi(token);
            break;
        }
    }

    // Prepare accumulators
    std::vector<char> expressed_genes(N_rows, 0);
    std::vector<RowBlock> row_blocks_tmp(N_rows);
    std::vector<double> gene_totals(N_rows, 0.0);
    std::vector<double> cell_totals(static_cast<std::size_t>(C), 0.0);
    for (int g = 0; g < N_rows; ++g)
    {
        row_blocks_tmp[g].offset = -1;
        row_blocks_tmp[g].nnz = 0;
    }

    // Ensure file sorted by row index
    int last_row_index = -1; // zero-based
    G = 0;
    int g_idx, c_idx;
    double count;
    while (true)
    {
        long long offset = infp.tellg();
        line_buffer = infp.getline();
        if (line_buffer == nullptr)
        {
            break;
        }
        if (line_buffer[0] == '\0' || line_buffer[0] == '\n' || line_buffer[0] == '%')
        {
            continue;
        }

        // Read values gene idx and add as expressed
        token = strtok(line_buffer, " \t\r\n");
        if (token == NULL)
        {
            continue;
        }
        g_idx = std::stoi(token) - 1;

        token = strtok(NULL, " \t\r\n");
        if (token == NULL)
        {
            fprintf(stderr, "Error: missing cell index for gene %d in MTX file %s\n", g_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }
        c_idx = std::stoi(token) - 1;

        token = strtok(NULL, " \t\r\n");
        if (token == NULL)
        {
            fprintf(stderr, "Error: missing count for gene %d cell %d in MTX file %s\n", g_idx + 1, c_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }
        count = std::stod(token);

        if (g_idx < 0 || g_idx >= N_rows)
        {
            fprintf(stderr, "Error: gene index %d out of range in MTX file %s\n", g_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }

        // Sorted-by-row check: smaller row index cannot follow larger
        if (last_row_index > g_idx)
        {
            fprintf(stderr, "Error: MTX file %s is not sorted by row. Row %d appears after row %d.\nRun scripts/sort_mtx_by_row.py on the file to sort it, e.g.: python3 sort_mtx_by_row.py -i %s\n", in_file.c_str(), g_idx + 1, last_row_index + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }
        last_row_index = g_idx;

        if (c_idx < 0 || c_idx >= C)
        {
            fprintf(stderr, "Error: cell index %d out of range in MTX file %s\n", c_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }

        if (!expressed_genes[g_idx])
        {
            expressed_genes[g_idx] = 1;
            G++;
            row_blocks_tmp[g_idx].offset = offset;
            row_blocks_tmp[g_idx].nnz = 1;
            row_blocks_tmp[g_idx].row_index = g_idx;
        }
        else
        {
            row_blocks_tmp[g_idx].nnz += 1;
        }
        gene_totals[g_idx] += count;
        cell_totals[static_cast<std::size_t>(c_idx)] += count;
    }
    infp.close();

    // Build gene_idx map, mtx_rows blocks and totals n (only for non-zero genes)
    mtx_rows.clear();
    mtx_rows.resize(G);
    n.assign(static_cast<std::size_t>(G), 0.0);
    int cur = 0;
    for (int g = 0; g < N_rows; ++g)
    {
        if (expressed_genes[g])
        {
            gene_idx[g] = cur;
            mtx_rows[cur] = row_blocks_tmp[g];
            n[static_cast<std::size_t>(cur)] = gene_totals[g];
            cur++;
        }
        else
        {
            gene_idx[g] = -1;
        }
    }
    N_c.assign(cell_totals.begin(), cell_totals.end());
}

std::vector<std::string> Read_CellNames(const std::string &filename)
{
    std::vector<std::string> cell_names;
    if (filename == "none" || filename.empty())
    {
        return cell_names;
    }
    FileReader infp(filename);
    if (!infp.is_open())
    {
        fprintf(stderr, "Cannot open input file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    char *line = nullptr;
    while ((line = infp.getline()) != nullptr)
    {
        if (line[strlen(line) - 1] == '\r')
        {
            line[strlen(line) - 1] = '\0';
        }
        cell_names.push_back(std::string(line));
    }
    return cell_names;
}

std::vector<std::string> Read_GeneNames(const std::string &filename,
                                        const std::map<int, int> &gene_idx,
                                        const int G)
{
    std::vector<std::string> gene_names;
    if (filename == "none" || filename.empty())
    {
        return gene_names;
    }
    FileReader infp(filename);
    if (!infp.is_open())
    {
        fprintf(stderr, "Cannot open input file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    char *line = nullptr;
    int gene_index = 0;
    while ((line = infp.getline()) != nullptr)
    {
        if (line[strlen(line) - 1] == '\r')
        {
            line[strlen(line) - 1] = '\0';
        }
        // Only add gene if its index in gene_idx is not -1
        auto it = gene_idx.find(gene_index);
        if (it != gene_idx.end() && it->second != -1)
        {
            gene_names.push_back(std::string(line));
        }
        gene_index++;
    }
    // Check that the length of gene_names equals G
    if (static_cast<int>(gene_names.size()) != G)
    {
        fprintf(stderr, "Error: Number of gene names (%d) does not match expected G (%d)\n",
                static_cast<int>(gene_names.size()), G);
        exit(EXIT_FAILURE);
    }
    return gene_names;
}