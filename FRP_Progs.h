#include <iostream>
#include <time.h>
#include <newmatio.h>
#include <math.h>
#include "Seq.h"
#include "FRP_Merton_class.h"
#include "FRP_Heston_class.h"
#include "Model.h"

ReturnMatrix AmericanCall_1D(const double &K, const double &rate, 
							   const double &maturity,	const double *dividends,
						      const double *dividends_times,const int &ndiv, Real *x, Real *y, 
							   const int &lenx, const int &leny, const int &iJ,const double dt,Model_1D *MDL );
ReturnMatrix AmericanPut_1D(const double &K, const double &rate, const double &ov_rate,
							   const double &maturity,const double *dividends,
						     const double *dividends_times,const int &ndiv, Real *x, Real *y, 
							   const int &lenx, const int &leny, const int &iJ,const double dt,Model_1D *MDL );
							   
ReturnMatrix AmericanCall_2D(const double &K, 
                          const double &rate, 
                          const double &ov_rate, 
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  const double &x0, const double &V0, Real *yy,
						  double *vv, double *grid, double *grid2, 
						  const int &ylen, const int &Wlen,
						  const int &iJ1 ,const int &iJ2, 
						  const double dt,
						  Heston_FRP *HST,const double &shift, const int &Nrowstruncation);

ReturnMatrix AmericanPut_2D(const double &K, 
                          const double &rate, 
                          const double &ov_rate,
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  const double &x0, const double &V0, Real *yy, 
						  double *vv, double *grid, double *grid2, 
						  const int &ylen, const int &Wlen,
						  const int &iJ ,const int &iJ2, 
						  const double dt,
						  Heston_FRP *HST,const double &shift, const int &Nrowstruncation);

double payoff_call( const double &s, const double &k );
double payoff_put ( const double &s, const double &k );
double payoff_binary_call(const double& S, const double& K);
double payoff_binary_put(const double& S, const double& K);
