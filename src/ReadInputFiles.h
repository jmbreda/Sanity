#ifndef _ReadInputFiles_h_
#define _ReadInputFiles_h_

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;
// For UMI count matrix
void Get_G_C_UMIcountMatrix(string in_file, int &N_rows, int &G, int &C, int N_char);

void ReadUMIcountMatrix(string in_file, double **n_c, double *N_c, double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, int N_char);

// For mtx file
void Get_G_C_MTX(string in_file, int &N_rows, int &G, int &C, map<int,int> &gene_idx);

void ReadMTX(string mtx_file, string gene_name_file, string cell_name_file, double **n_c, double *N_c,    double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, map<int,int> gene_idx);

#endif

