#ifndef _FitFrac_h_
#define _FitFrac_h_

#include <vector>

double fitfrac(double *f, const std::vector<double>& n_c, double n, double v, int C, const std::vector<double>& N_c, double prev_q, const double* log_N_c = nullptr, double library_total = 0, double* workspace = nullptr);
void fractions_at_q(double *f, const std::vector<double>& n_c, double n, double v, int C, const std::vector<double>& N_c, double q, const double* log_N_c = nullptr);
double deltaq(double *Q,int C,double beta, double q);
double fq(double *Q,int C,double beta, double q);
double LambertW0_approximation(const double& x);
#endif
