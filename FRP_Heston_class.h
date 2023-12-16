#ifndef FRP_Heston_class
#define FRP_Heston_class

#include <newmatio.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <complex>
#include "Model.h"
#include "Seq.h"
#include "fftw3.h"
#include "haar.h"
//#include "nr.h"

using namespace std;
typedef enum { Y, W } Dimension;

class Heston_FRP
	{
	public:
		Heston_FRP(const double rate,//rate
		          const double dyield,
				  const double longrun_var_Heston, 
				  const double mean_rev_Heston,
				  const double volvol_Heston,
				  const double rho,
				  double lambda,
				  double muj,
				  double sigmaj,
				  const int Nrowstruncation,const int ylen, const int Wlen,
				  const int glen, const int glen2,
				  const double *yy, const double *vv,
				  const double* grid,const double* grid2,
				  const double Delta, const double Delta2,const double shift,
				  const char Method);
		~Heston_FRP();
		void Make_Green_0( const double &x0, const double &v0, const double& t);
		void Make_Green(const double &t);
		void Make_Green_maturity(const double &t);
		void Make_Green_maturity(const double &t,const double &x0, const double &v0);
		
		
		double project_back_start( Matrix &mV );
		 double project_back_start( Matrix &mV ,const int flag);

		ReturnMatrix project_back( Matrix &mV );
		ReturnMatrix project_back_maturity( Matrix &mV,int *iIndDrop,vector<double> &vSmallDelta,vector<int>&vSmallDeltaSign,
		                                                    const double &K,double (*pt2Payoff)(const double &s, const double &k) );
		ReturnMatrix project_back_dividend( Matrix &mV,
                                                                      int *iIndDrop,vector<double> &vSmallDelta,
                                                                      vector<int>&vSmallDeltaSign,const double &K,
                                                                      double (*pt2Payoff)(const double &s, const double &k) );

		complex<double> Psic( const double &u1,const double& u2,const double& r,const double& dyield,const double &tau,const double& ka,const double &theta, 
					 const double &alpha,const double& rho,const double& x0,const double &v0);
		void PsiFn (double T,	double U,	double V,double lambda,	double muJ,	double sigmaJ,	double a,double b,		double alpha,	double rho,	double x1,		double x2,		double *RePsi,		double *ImPsi);
		
		ReturnMatrix TransitionM( const double &t, const double &x0, const double &v0, int &Truncation);
		ReturnMatrix TransitionM_1D( const double &t, const double &x0, const double &v0);
        
        double dotproduct_sparse(Matrix &A, Matrix &Sparse, Matrix &Coordinates, int & vmTrowsStart,int &vmTrowsEnd);		
		void Set_Rate(double Rate);
		
		double Get_dyield();
		//Sampling functions
		void SampleHaar( const Dimension WhichDimension =Y);
		
	private:
		double  rate,dyield, longrun_var_Heston, mean_rev_Heston,volvol_Heston,rho,lambda, muj, sigmaj;
		int Nrowstruncation,ylen,Wlen,glen,glen2;		
		const double *yy,*vv,*grid, *grid2;
		double Delta, Delta2,shift;
		char method;
		int ik;
		
		bool is_0,is_mG,is_maturity_mG;		
		Matrix mG_0;
		vector<Matrix> vmT,maturity_mG;
		vector<Matrix> vmT_nonZero;
		
		//Parameters for the sampling method
		
		Matrix mRePsic, mImPsic;
		Matrix mReTheta, mImTheta; 
		Matrix mReEta, mImEta;
		Matrix mLim,mLim2;
		
	
	};
#endif 


