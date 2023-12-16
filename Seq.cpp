#include "Seq.h"

#ifndef PI 
#define PI 3.141592653589793238462643
#endif

double n(const double& z) {  // normal distribution function    
    return (1.0/sqrt(2.0*PI))*exp(-0.5*z*z);
};


#ifdef __cplusplus
extern "C"{
#endif
	int seqa( const double &a, const double &b, const size_t &leny, double *y){
		
		std::vector<Real> dVy(leny); 
		double step( (b - a)/(leny - 1) );
		generate(dVy.begin(), dVy.end(), aritGen<Real>( a, step ) );
		copy(dVy.begin(), dVy.end(), y);
		
		return 0;
	}
#ifdef __cplusplus
}
#endif
