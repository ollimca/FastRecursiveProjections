#ifndef FRP_Merton_class
#define FRP_Merton_class

#include <newmatio.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include "Model.h"
#include "Seq.h"

using namespace std;

class Merton_FRP: public Model_1D
	{
	public:	
		Merton_FRP(const Real rate,//rate
		          const Real dyield,
				  const Real sigma, //vol
				  const Real lambda,//mean jumps per year
				  const Real muj,//jump mean
				  const Real sigmaj,//jump vol
				  const int n  ,//Maximum number of jumps that can "resonably" happen between two steps (for approximation)
				  const int lower=0,
				  const int upper=0);
		~Merton_FRP();
		void Make_Green( const int m, const int k, 
						  Real *xx, Real *yy, const int J, Real t);
		void Make_Green_Call_Maturity( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
	   void Make_Green_Call_dyield( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);	
	   void Make_Green_Call_DD( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
	   void Make_Green_Call_0d( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
	   void Make_Green_Put_last( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
	void Make_Green_Put_yydtSI( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
		void Make_Green_Put_yyintdt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
		void Make_Green_Put_yintydt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
		void Make_Green_Put_yintytime( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t);
		void Get_Green_Put_last(Matrix &mM) const;				  
		void Get_Green_Put_yydtSI(Matrix &mM) const;
		void Get_Green_Put_yyintdt(Matrix &mM) const;
		void Get_Green_Put_yintydt(Matrix &mM) const;
		void Get_Green_Put_yintytime(Matrix &mM) const;
						
		void Make_Green_SpaceInvariance( const int lenx,//number of starting points
                              const int leny, //number of arrivals  points
							  Real *xx, //grid of starting points
							  Real *yy, //grid of arrivals points
							  const int J,//precision
							  Real t);//time step
		void Get_Green( Matrix & mM ) const;
		void Get_Green_Call_mat(Matrix & mM) const;  
		void Get_Green_Call_DD(Matrix & mM) const;    
		void Get_Green_Call_0d(Matrix & mM) const;  
      
		bool Get_IsBand_mat() const;
		bool Get_IsGreen_mat() const ;
		ColumnVector project_back( ColumnVector &Y );
		ColumnVector project_back( ColumnVector &Y ,Matrix &mM);
		ColumnVector project_back_Call_maturity( ColumnVector &Y );
		ColumnVector project_back_Call_dyield( ColumnVector &Y );
		ColumnVector project_back_Call_maturity_start( ColumnVector &Y );
	    ColumnVector project_back_Call_DD( ColumnVector &Y );
	    ColumnVector project_back_Call_0d( ColumnVector &Y );
	     ColumnVector project_back_Put_last( ColumnVector &Y );
	     ColumnVector project_back_Put_yydtSI( ColumnVector &Y );
	     ColumnVector project_back_Put_yyintdt( ColumnVector &Y );
	     ColumnVector project_back_Put_yintydt( ColumnVector &Y );
	     ColumnVector project_back_Put_yintytime( ColumnVector &Y );
	     	     
		void Set_Rate(double Rate);
		void Set_dL(double DL);
		double Get_dyield();
		
		private:
		Real  rate,dyield, sigma, time_to_mat;
		int iJ; 		   
		Real Muj, Sigmaj, Lambda;
		bool isGreen_mat;
		int N;
		int lower,upper;
		Real kappa;
		bool isBand_mat;
		Real *x, *y;
		double dL;
		Matrix mG,Call_dyield_mG,Call_maturity_mG,Call_DD_mG,Call_0d_mG;
        BandMatrix bG,Call_dyield_bG,Call_maturity_bG,Call_DD_bG,Call_0d_bG;
        bool is_Call_mat_mG,is_Call_dyield_mG,is_Call_DD_mG,is_Call_0d_mG;
        Matrix Put_last_mG,Put_yydtSI_mG,Put_yyintdt_mG,Put_yintydt_mG,Put_yintytime_mG;
        BandMatrix Put_last_bG,Put_yydtSI_bG,Put_yyintdt_bG,Put_yintydt_bG,Put_yintytime_bG;
		bool is_Put_last_mG,is_Put_yydtSI_mG,is_Put_yyintdt_mG,is_Put_yintydt_mG,is_Put_yintytime_mG;
		int lenx, leny;
	
	};
#endif 

