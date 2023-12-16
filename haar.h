#include <cmath>
#include <vector>
#include <newmatio.h>
#include <algorithm>
#include "Seq.h"


//using namespace NEWMAT;


#ifndef _HAAR_H_
#define _HAAR_H_

////////////////////////////////////////
// HAAR Lifting Scheme
////////////////////////////////////////

void haarm2dr( Matrix& mM, int iJ );

void haarm2dc( Matrix& mM, int iJ );

ReturnMatrix haarl1l2sq( int iJ );

ReturnMatrix haarl1l2( int iJ );

ReturnMatrix make_limits( const Matrix &mM1, int iJ);

	
////////////////////////////////////////
// Fourier Transform of the Haar basis
////////////////////////////////////////


// Analytical definitions, functions defined x \in [d,u]

double ReEhat_du( double x, double d, double u );
double ImEhat_du( double x, double d, double u );
double RePsihat_du(double x, double d, double u);
double ImPsihat_du(double x, double d, double u);

//Matrix Forms

ReturnMatrix ReEhat_du_matrix( const size_t k, const Real *yy, double d, double u );
ReturnMatrix ImEhat_du_matrix( const size_t k, const Real *yy, double d, double u );
ReturnMatrix RePsihat_du_matrix( const size_t k, Real *yy, double d, double u );
ReturnMatrix ImPsihat_du_matrix( const size_t k, Real *yy, double d, double u );


#endif /*_HAAR_H_*/
