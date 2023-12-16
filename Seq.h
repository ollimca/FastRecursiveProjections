

#ifndef _UTILS_H_
#define _UTILS_H_

#ifndef Real
#define Real double
#endif  

#ifndef PI 
#define PI 3.141592653589793238462643
#endif

#include <algorithm>
#include <vector>
#include <math.h>

double n(const double& z) ;


template< class  T  > 
class aritGen 
{
public:
  aritGen (T start=0, T step=1 ) : current(start), step(step) { }
  T operator() () { T tmp= current; current += step ; return tmp; }
private:
  T current;
  T step ;
};


#ifdef __cplusplus
extern "C"{
#endif
  int seqa( const double &a, const double &b, const size_t &leny, double *y);
#ifdef __cplusplus
}
#endif

#endif
