#include "FitFrac.h"
#include "FukushimaLambertW.h"
#include <iostream>
#include <iomanip> // Required for setprecision
#include <cmath>

double fitfrac(double *f, const std::vector<double>& n_c, double n, double &v, int C, const std::vector<double>& N_c, double a, double b, double prev_q) {
    double q, W, x, dq, funcq;
    int i;
    double *Q = new double [C];
    double beta = n * v;
    double inv_beta = 1.0/beta;
    double logbeta = std::log(beta);

    q = 0.0;
    for(i=0;i<C;++i){
        Q[i] = std::log( N_c[i]) + n_c[i]*v + logbeta;//Q[i] = y_gc +log(nv) = log(Nc)+nc v + log(beta)
        q += N_c[i]; //q = log[sum_c N_c]
    }
    //Initial guess for q is q = log(sum_c Nc))+v/2 because at optimum q = log[ sum_c Nc exp(delta_c)] and expectation of exp(delta_c) is exp(v/2)
    q = std::log(q) + 0.5*v;

    if (prev_q != 0.0) {
        funcq = fq(Q, C, beta, prev_q);
        if (funcq > 0) {
            // do a single Newton step
            dq = deltaq(Q, C, beta, prev_q);
            prev_q = prev_q - dq;
            funcq = fq(Q, C, beta, prev_q);
            if (funcq > 0) {
                // no luck
                prev_q = 0.0;
            }
            else {
                q = prev_q;
            }
        }
        else {
            q = prev_q;
        }
    }

    if (prev_q == 0.0) {
        //We are looking for the zero of the function f(q) = beta -sum_c W[e^(Q_c-q)]
        //This function f(q) is a strictly increasing and  concave function of q
        //Consequently, if we start at a value of q with f(q) < 0. Then Newton's method is guaranteed to converge,
        //because the linear guess is always above the curve of f(q), so we will approach zero from below.
        //So we check if fq is negative at our guess.. if not.. we keep lowering q until it is below zero
        funcq = fq(Q, C, beta, q);
        while(funcq > 0){
            q -= 0.5*v; //This is just a rough guess for the stepsize, i.e lowering q by expected average of e^delta_c
            funcq = fq(Q , C, beta, q);
        }
    }

    double tol = 1e-7;
    dq = -1.0;
    while(std::fabs(dq) > tol && dq < 0){
        //dq = f(q)/f'(q)
        dq = deltaq(Q,C,beta,q);
        q = q - dq;//Newton's method update
    }

    funcq = beta;

    for(i=0;i<C;++i){
        x = Q[i] -q;
        if(x > 50){
            // at large x LambertW breaks so we use an approximation
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,std::exp(x));
        }
        f[i] = W*inv_beta;
        funcq -= W;
    }
    if(std::fabs(dq/q) > 1e-6){
        std::cerr << std::fixed << std::setprecision(12);
        std::cerr << "fitfrac: WARNING: dq = " << dq << " is not close to zero at final q = " << q
            << " funcq = " << funcq << std::endl;
    }
  delete[] Q;
  return q;
}


////calculates the function f(q) = beta -sum_c W[e^{Q_c-q)]
double fq(double *Q,int C,double beta, double q){

    double funcq = beta;
    double W, x;
    int i;

    for(i=0;i<C;++i){
        x = Q[i] - q;
        if(std::isnan(x)){
            //cout << "x = nan !\n";
            throw std::runtime_error("Error: x is NaN in fq function.");
        }
        if(x > 50.0){
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,std::exp(x));
        }
        funcq -= W;
    }

    return funcq;
}


//Estimate the ratio f(q)/f'(q)
double deltaq(double *Q,int C,double beta, double q){
    double num = beta;
    double denom = 0;
    double x, W;
    int i;

    for(i=0;i<C;++i){
        x = Q[i] - q;
        if(std::isnan(x)){
            throw std::runtime_error("Error: x is NaN in deltaq function.");
        }
        if(x > 50.0){
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,std::exp(x));
        }
        num -= W;
        denom += W/(1.0+W);
    }

    double rat = num/denom;
    return rat;
}

double LambertW0_approximation(const double& x) {
    // our Lambert approximation
    double L2 = std::log(x);
    double xsq = x*x;
    double L2sq = L2*L2;
    double L2cub = L2sq*L2;
    double W = x - L2 + L2/x +
        L2*(L2 - 2.0)/(2.0*xsq) +
        L2*(6.0 - 9.0*L2 + 2.0*L2sq)/(6.0*xsq*x) +
        L2*(-12.0 + 36.0*L2 - 22.0*L2sq + 3.0*L2cub)/(12.0*xsq*xsq);

    return W;
}
