#ifndef _ReadInputFiles_h_
#define _ReadInputFiles_h_

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cstring>

struct RowBlock {
	long long offset;
	std::size_t nnz;
	RowBlock() : offset(-1), nnz(0) {}
};

using namespace std;
// For UMI count matrix
void Get_G_C_UMIcountMatrix(string in_file,
							int &N_rows,
							int &G,
							int &C,
							int N_char,
							std::vector<std::streampos> &tsv_offsets,
							std::vector<double> &N_c,
							std::vector<double> &n,
							std::vector<std::string> &cell_names,
							std::vector<std::string> &gene_names);

void ReadUMIcountMatrix(string in_file, double **n_c, double *N_c, double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, int N_char);

// For mtx file
void Get_G_C_MTX(string in_file, int &N_rows, int &G, int &C, map<int,int> &gene_idx, vector<RowBlock> &mtx_rows, vector<double> &N_c, vector<double> &n);
void Get_G_C_MTX_2(string in_file, int &N_rows, int &G, int &C, map<int,int> &gene_idx, vector<RowBlock> &mtx_rows, vector<double> &N_c, vector<double> &n);

void ReadMTX(string mtx_file, string gene_name_file, string cell_name_file, double **n_c, double *N_c,    double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, map<int,int> gene_idx);

std::vector<std::string> Read_CellNames(const std::string &filename);
std::vector<std::string> Read_GeneNames(const std::string &filename,
									   const std::map<int,int> &gene_idx,
									   const int G);

#endif

