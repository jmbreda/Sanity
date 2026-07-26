#include <FitFrac.h>
#include <iostream>
#include <iomanip> // Required for setprecision

double fitfrac(double *f, const std::vector<double>& n_c, double n, double &v, int C, const std::vector<double>& N_c, double a, double b){
    double q,W,x;
    int i;
    double *Q = new double [C];
    double beta = n * v;
    double inv_beta = 1.0/beta;
    double logbeta = log(beta);
    q=0;

    //std::cerr << "fitfrac: beta = " << beta << " logbeta = " << logbeta << std::endl;

    for(i=0;i<C;++i){
        Q[i] = log( N_c[i]) + n_c[i]*v + logbeta;//Q[i] = y_gc +log(nv) = log(Nc)+nc v + log(beta)
        q += N_c[i]; //q = log[sum_c N_c]
        //std::cerr << "fitfrac: N_c[" << i << "] = " << N_c[i] << " n_c[" << i << "] = " << n_c[i]
        //<< " Q[" << i << "] = " << Q[i] << std::endl;
    }
    //Initial guess for q is q = log(sum_c Nc))+v/2 because at optimum q = log[ sum_c Nc exp(delta_c)] and expectation of exp(delta_c) is exp(v/2)
    q = log(q)+0.5*v;
    //std::cerr << "fitfrac: Initial q = " << q << std::endl;

    //We are looking for the zero of the function f(q) = beta -sum_c W[e^(Q_c-q)]
    //This function f(q) is a strictly increasing and  concave function of q
    //Consequently, if we start at a value of q with f(q) < 0. Then Newton's method is guaranteed to converge,
    //because the linear guess is always above the curve of f(q), so we will approach zero from below.
    //So we check if fq is negative at our guess.. if not.. we keep lowering q until it is below zero
    double funcq = fq(Q,C,beta,q);
    //std::cerr << "fitfrac: funcq = " << funcq << std::endl;
    while(funcq > 0){
        q -= 0.5*v;//This is just a rough guess for the stepsize, i.e lowering q by expected average of e^delta_c
        funcq = fq(Q,C,beta,q);
        //std::cerr << "fitfrac: q = " << q << " funcq = " << funcq  << std::endl;
    }

    double tol = 1e-7;
    double dq = -1.0;
    while(fabs(dq) > tol && dq < 0){
        //dq = f(q)/f'(q)
        dq = deltaq(Q,C,beta,q);
        q = q - dq;//Newton's method update
        //std::cerr << "fitfrac: dq =  " << dq << " new q = " << q  << std::endl;
        //funcq = fq(Q,C,beta,q);
        //std::cerr << "fitfrac update: funcq = " << funcq << std::endl;
    }

    // commented this block
    // funcq is calculated along with f[i] in the next
    /* funcq = fq(Q,C,beta,q);
    if(fabs(dq/q) > 1e-6){
        std::cerr << std::fixed << std::setprecision(12);
        std::cerr << "fitfrac: WARNING: dq = " << dq << " is not close to zero at final q = " << q
            << " funcq = " << funcq << std::endl;
    } */

    funcq = beta;

    for(i=0;i<C;++i){
        x = Q[i] -q;
        if(x > 50){
            // at large x LambertW breaks so we use an approximation
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,exp(x));//But now I worry about taking exponent when x is large and feeding this function such a large value.
        }
        f[i] = W*inv_beta;
        funcq -= W;
    }
    if(fabs(dq/q) > 1e-6){
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
        if(isnan(x)){
            //cout << "x = nan !\n";
            throw std::runtime_error("Error: x is NaN in fq function.");
        }
        if(x > 50.0){
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,exp(x));
        }
        //double res = W + log(W);
        //std::cerr << "fq: x =  " << x << " res = " << res <<  " x - res =  " << x - res << std::endl;
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
        if(isnan(x)){
            throw std::runtime_error("Error: x is NaN in deltaq function.");
        }
        if(x > 50.0){
            W = LambertW0_approximation(x);
        }
        else{
            W = Fukushima::LambertW(0,exp(x));
        }
        //double res = W + log(W);
        //std::cerr << "deltaq: x =  " << x << " res = " << res << " x - res =  " << x - res << std::endl;
        num -= W;
        denom += W/(1.0+W);
    }

    double rat = num/denom;
    return rat;
}

double LambertW0_approximation(const double& x) {
    // our Lambert approximation
    double L2 = log(x);
    double xsq = x*x;
    double L2sq = L2*L2;
    double L2cub = L2sq*L2;
    double W = x - L2 + L2/x +
        L2*(L2 - 2.0)/(2.0*xsq) +
        L2*(6.0 - 9.0*L2 + 2.0*L2sq)/(6.0*xsq*x) +
        L2*(-12.0 + 36.0*L2 - 22.0*L2sq + 3.0*L2cub)/(12.0*xsq*xsq);

    return W;
}
