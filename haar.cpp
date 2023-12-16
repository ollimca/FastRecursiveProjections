#include "haar.h"

// READ!!
// in the following codes j=0 is the finest resolution level, j=iJ is the coarsest (digital - processing convention)

void haarm2dr( Matrix& mM, int iJ )
{ 	  
	  for( int j=0; j <iJ; j++){
	  	int len = int( ldexp ( 1, iJ-j) ); int step = int( ldexp ( 1, j) );
      
	    vector<int> iVec(len);
	    generate(iVec.begin(), iVec.end(), aritGen<int>( 1, step ) );
		  
	    for( vector<int>::size_type i=0; i< (iVec.size() / 2) ; ++i)
		{
			mM.Column( iVec[2*i + 1] )   = mM.Column( iVec[2*i + 1] ) - mM.Column( iVec[2*i]);	
		    mM.Column( iVec[2*i] )       = mM.Column( iVec[2*i] )     + mM.Column( iVec[2*i + 1] ) / 2;
		}
	  }
}

void haarm2dc( Matrix& mM, int iJ )
{ 	  
	  for( int j=0; j <iJ; j++){
	  	int len = int( ldexp ( 1, iJ-j) ); int step = int( ldexp ( 1, j) );
      
	    vector<int> iVec(len);
	    generate(iVec.begin(), iVec.end(), aritGen<int>( 1, step ) );
		  
	    for( vector<int>::size_type i=0; i< (iVec.size() / 2) ; ++i)
		{
			mM.Row( iVec[2*i + 1] )   = mM.Row( iVec[2*i + 1] ) - mM.Row( iVec[2*i]);	
		    mM.Row( iVec[2*i] )       = mM.Row( iVec[2*i] )     + mM.Row( iVec[2*i + 1] ) / 2;
		}
	  }
}

ReturnMatrix haarl1l2sq( int iJ )
{
		int veclen = int( ldexp ( 1, iJ ) );
		Matrix mConver(1, veclen);
			
	  	for( int j=0; j < iJ; j++){ // here with respect to the R code, we start from j=0,
	  								// so the normalising factor is (j+1) / 2 instead of j/2
	  	int len = int( ldexp ( 1, iJ-j) ); int step = int( ldexp ( 1, j) );
      
	    vector<int> iVec(len);
	    generate(iVec.begin(), iVec.end(), aritGen<int>( 1, step ) );
		  
	    for( vector<int>::size_type i=0; i< (iVec.size()) ; ++i)
				mConver(1, iVec[i] ) = pow( 2.0, (j - 1) );
	  }
	  mConver(1,1) *= 4;
	  mConver.Release();
	  return mConver;
} 

ReturnMatrix haarl1l2( int iJ )
{
		int veclen = int( ldexp ( 1, iJ ) );
		Matrix mConver(1, veclen);
			
	  	for( double j=0.0; j < (double)iJ; j++){ // here with respect to the R code, we start from j=0,
	  								// so the normalising factor is (j+1) / 2 instead of j/2
	  	int len = int( pow ( 2.0, double(iJ-j) ) ); int step = int( pow ( 2.0, j) );
      
	    vector<int> iVec(len);
	    generate(iVec.begin(), iVec.end(), aritGen<int>( 1, step ) );
		  
	    for( vector<int>::size_type i=0; i< (iVec.size()) ; ++i)
				mConver(1, iVec[i] ) = pow( 2.0, (j/2 - 0.5) );
	  	}

	mConver(1,1) *= 2;
	mConver.Release();
	return mConver;
} 

ReturnMatrix make_limits( const Matrix &mM1, int iJ){
	
	size_t len = (size_t)pow(2.0, iJ);
	int step, l;
	Matrix mL(2, len), mO(1,1), mM(mM1); 
	
	mO << -mM(1,2) + 2*mM(1,1);
	mM = mO | mM;	
	
	for(int j=1; j!=iJ+1; ++j){
		step = (int)pow(2.0, j);
		for( size_t i=0; i!= (len/step ); i++ ){
			l = ( (int)pow(2.0, j-1) + 1 +i*step ) ;
			mL(1,l)= mM(1, step*i + 1);
			mL(2,l)= mM(1, step*(i+1) + 1 );
		}
	}
	mL(1,1) = mM(1,1);
	mL(2,1) = mM(1, mM.Ncols() );
	
	mL.Release(); return mL;	
}

//*******************************************************************************
// Fourier Transform of the Haar basis
//*******************************************************************************

//	Analytical form


double ReEhat_du( double x, double d, double u ){
	
	double y = x/2.0;
	
	return x != 0.0 
	? ( 1 / sqrt(u-d) * cos( ( u + d )* y ) * sin( (u-d)*y ) / y ) 
	: ( sqrt(u-d) ) ;
}

double ImEhat_du( double x, double d, double u ){
	
	double y = x/2.0;
	
	return x != 0.0 
	? (-1 /sqrt(u-d) * sin( ( u + d) * y )  * sin( (u-d) * y ) / y ) 
	: 0.0;
}

double RePsihat_du(double x, double d, double u){
	
	double m = ( d + u )/2;
	
	return 1/sqrt(2.0) * ( -ReEhat_du( x, d, m ) + ReEhat_du( x, m, u ) );
	
}

double ImPsihat_du(double x, double d, double u){
	
	double m = ( d + u )/2;
	
	return 1/sqrt(2.0) * ( -ImEhat_du( x, d, m ) + ImEhat_du( x, m, u ) );
	
}

//Matrix Forms


ReturnMatrix ReEhat_du_matrix( const size_t k, const Real *yy, double d, double u ){
	
	Matrix mE; mE.resize(1,k) ; //initialize the ehat matrix (vector)
	
	for(size_t j=0; j!= k; ++j )	//Build hathaar matrix
		mE( 1, j + 1 ) = ReEhat_du( yy[j], d, u ) ;
	
	mE.release(); return mE;
}

ReturnMatrix ImEhat_du_matrix( const size_t k, const Real *yy, double d, double u ){
	
	Matrix mE; mE.resize(1,k) ; //initialize the ehat matrix (vector)
	
	for(size_t j=0; j!= k; ++j )	//Build hathaar matrix
		mE( 1, j + 1 ) = ImEhat_du( yy[j], d, u ) ;
	
	mE.release(); return mE;
}

ReturnMatrix RePsihat_du_matrix( const size_t k, Real *yy, double d, double u ){
	
	Matrix mE; mE.resize(1,k) ; //initialize the ehat matrix (vector)
	
	for(size_t j=0; j!= k; ++j )	//Build hathaar matrix
		mE( 1, j + 1 ) = RePsihat_du( yy[j], d, u ) ;
	
	mE.release(); return mE;
}

ReturnMatrix ImPsihat_du_matrix( const size_t k, Real *yy, double d, double u ){
	
	Matrix mE; mE.resize(1,k) ; //initialize the ehat matrix (vector)
	
	for(size_t j=0; j!= k; ++j )	//Build hathaar matrix
		mE( 1, j + 1 ) = ImPsihat_du( yy[j], d, u ) ;

	  	mE.release(); return mE;
}



