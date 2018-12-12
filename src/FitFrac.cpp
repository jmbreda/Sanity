#include <FitFrac.h>

double fitfrac(double *f, double *n_c, double n, double &v, int C, double *N_c, double a, double b){
    double q,qmin,qmax,nor;
    int i;
	double *Q = new double [C];
	double beta = 1.0/((n+a)*v);
	double logbeta = log(beta);
	q=0;
	for(i=0;i<C;++i){
        Q[i] = log( N_c[i]) + n_c[i]*v -logbeta;
		q += N_c[i];
    }
	q =log(q+b);//initial guess for q
	double dq = log(2.0);

    /**** get initial nor ****/
    nor = beta*normalization(Q,C,q)+b*exp(-q);
    /**** initial c too small. Set as lower bound and increase by factor 2 until a value too large is found****/
    if(nor < 1){
        while(nor < 1){
            qmax = q;
            q = q-dq;
            nor = beta*normalization(Q,C,q)+b*exp(-q);
        }
        qmin = q;
    }
    /**** initial c too large. Set as upper bound and decrease by factor 2 until a value too small is found****/
    else{
        while(nor > 1){
            qmin = q;
            q = q+dq;
            nor = beta*normalization(Q,C,q)+b*exp(-q);
        }
        qmax = q;
    }
    /***now do the bisection*****/
    double tol = 0.0000001;
    double diff = 1.0;
    while(diff > tol){
        q = (qmax+qmin)/2;
    	nor = beta*normalization(Q,C,q)+b*exp(-q);
        if(nor > 1){
            qmin = q;
        }else{
            qmax = q;
        }
        diff = fabs(nor-1.0);
    }

    q = (qmax+qmin)/2;
	for(i=0;i<C;++i){
        double x = Q[i] -q;
        if(x > 50){
            f[i] = beta*(x-log(x)+log(x)/x);
        }else{
            //f[i] = beta*gsl_sf_lambert_W0(exp(x));
			f[i] = beta*Fukushima::LambertW(0,exp(x));
        }
    }
	
	delete[] Q;
    return q;
}

double normalization(double *Q,int C,double q){
    double nor = 0;
    double x;
    int i;
    for(i=0;i<C;++i){
        x = Q[i]-q;
		if(isnan(x)){
			cout << "x = nan !\n";
		}
        if(x > 50){
            nor += (x-log(x)+log(x)/x);
        }else{
            //nor += gsl_sf_lambert_W0(exp(x));
			nor += Fukushima::LambertW(0,exp(x));
        }
    }
    return nor;
}

