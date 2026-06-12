#include <ReadInputFiles.h>
#include <omp.h>
#include <atomic>

void Get_G_C_UMIcountMatrix(string in_file,
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
        size_t pos = name.find('\r');
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
        N_c.assign(static_cast<size_t>(C), 0.0);
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
            size_t pos = gene_id.find('\r');
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
            double value;
            for (int c = 0; c < C; ++c)
            {
                thread_token = strtok_r(NULL, " \t,", &thread_saveptr);
                if (thread_token != NULL)
                {
                    double value = stod(thread_token); /**total count this gene**/
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
                    delete[] thread_sc;
                    continue;
                }
            }
            delete[] thread_sc;
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
            for (size_t c = 0; c < static_cast<size_t>(C); ++c)
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

    return;
}

void Get_G_C_MTX(string in_file, int &N_rows, int &G, int &C, map<int, int> &gene_idx, vector<RowBlock> &mtx_rows, vector<double> &N_c, vector<double> &n)
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
            N_rows = stoi(token);
            token = strtok(NULL, " \t");
            // Number of cells
            C = stoi(token);
            break;
        }
    }

    // Prepare accumulators
    vector<char> expressed_genes(N_rows, 0);
    vector<RowBlock> row_blocks_tmp(N_rows);
    vector<double> gene_totals(N_rows, 0.0);
    vector<double> cell_totals(static_cast<size_t>(C), 0.0);
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
        g_idx = stoi(token) - 1;

        token = strtok(NULL, " \t\r\n");
        if (token == NULL)
        {
            fprintf(stderr, "Error: missing cell index for gene %d in MTX file %s\n", g_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }
        c_idx = stoi(token) - 1;

        token = strtok(NULL, " \t\r\n");
        if (token == NULL)
        {
            fprintf(stderr, "Error: missing count for gene %d cell %d in MTX file %s\n", g_idx + 1, c_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }
        count = stod(token);

        if (g_idx < 0 || g_idx >= N_rows)
        {
            fprintf(stderr, "Error: gene index %d out of range in MTX file %s\n", g_idx + 1, in_file.c_str());
            infp.close();
            exit(EXIT_FAILURE);
        }

        // Sorted-by-row check: smaller row index cannot follow larger
        if (last_row_index > g_idx)
        {
            fprintf(stderr, "Error: MTX file %s is not sorted by row. Row %d appears after row %d.\n", in_file.c_str(), g_idx + 1, last_row_index + 1);
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
        }
        else
        {
            row_blocks_tmp[g_idx].nnz += 1;
        }
        gene_totals[g_idx] += count;
        cell_totals[static_cast<size_t>(c_idx)] += count;
    }
    infp.close();

    // Build gene_idx map, mtx_rows blocks and totals n (only for non-zero genes)
    mtx_rows.clear();
    mtx_rows.resize(G);
    n.assign(static_cast<size_t>(G), 0.0);
    int cur = 0;
    for (int g = 0; g < N_rows; ++g)
    {
        if (expressed_genes[g])
        {
            gene_idx[g] = cur;
            mtx_rows[cur] = row_blocks_tmp[g];
            n[static_cast<size_t>(cur)] = gene_totals[g];
            cur++;
        }
        else
        {
            gene_idx[g] = -1;
        }
    }
    N_c.assign(cell_totals.begin(), cell_totals.end());

    return;
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
        cell_names.push_back(string(line));
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
            gene_names.push_back(string(line));
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