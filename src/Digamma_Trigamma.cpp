#include "Digamma_Trigamma.h"
#include <cmath>

double Psi_0(double x){

    if ( x == 1.0 )
        return -0.577215664901532;
    else
        return std::log(x) - 1.0/(2.0*x) - 1.0/(12.0*std::pow(x,2)) + 1.0/(120.0*std::pow(x,4)) - 1.0/(252.0*std::pow(x,6)) + 1.0/(240.0*std::pow(x,8)) - 5.0/(660.0*std::pow(x,10)) + 691.0/(32760.0*std::pow(x,12)) - 1.0/(12.0*std::pow(x,14));

}

double Psi_1(double x){
	
    if ( x == 1.0 )
    	return 1.644934066848226;
    else
		return 1.0/x + 1.0/(2.0*std::pow(x,2)) + 1.0/(6.0*std::pow(x,3)) - 1.0/(30.0*std::pow(x,5)) + 1.0/(42.0*std::pow(x,7)) - 1.0/(30.0*std::pow(x,9)) - 5.0/(66.0*std::pow(x,11)) - 691.0/(2730.0*std::pow(x,13)) + 7.0/(6.0*std::pow(x,15));

}
