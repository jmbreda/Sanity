#ifndef _FitFrac_h_
#define _FitFrac_h_

#include <math.h>
#include <iostream>
#include <FukushimaLambertW.h>
#include <vector>

using namespace std;

double fitfrac(double *f, const std::vector<double>& n_c, double n, double &v, int C, const std::vector<double>& N_c, double a, double b);
double normalization(double *Q,int C,double q);

#endif
