#ifndef Model_class_H
#define Model_class_H
#include <newmatio.h>
class Model_1D
	{
	public:	
		Model_1D( );
		virtual ~Model_1D(); 
		virtual void Make_Green( const int m, const int k, 
						  Real *xx, Real *yy, const int J,Real t)=0;
		virtual void Make_Green_Call_Maturity( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
	   virtual void Make_Green_Call_DD( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
	   virtual void Make_Green_Call_0d( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
	   virtual void Make_Green_Call_dyield( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;						  
	    virtual void Make_Green_Put_last( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
	 virtual void Make_Green_Put_yydtSI( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
		 virtual void Make_Green_Put_yyintdt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
		 virtual void Make_Green_Put_yintydt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
		 virtual void Make_Green_Put_yintytime( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)=0;
		 virtual void Get_Green_Put_last(Matrix &mM) const=0;				  
		 virtual void Get_Green_Put_yydtSI(Matrix &mM) const=0;
		 virtual void Get_Green_Put_yyintdt(Matrix &mM) const=0;
		 virtual void Get_Green_Put_yintydt(Matrix &mM) const=0;
		 virtual void Get_Green_Put_yintytime(Matrix &mM) const=0;	
		 virtual bool Get_IsBand_mat() const=0;					  
		virtual	void Make_Green_SpaceInvariance( const int lenx,//number of starting points
                              const int leny, //number of arrivals  points
							  Real *xx, //grid of starting points
							  Real *yy, //grid of arrivals points
							  const int J,//precision
							  Real t)=0;//time step
		virtual void Get_Green( Matrix & mM ) const=0;
		virtual void Get_Green_Call_mat(Matrix & mM) const=0;  
		virtual void Get_Green_Call_DD(Matrix & mM) const=0;    
		virtual void Get_Green_Call_0d(Matrix & mM) const=0;  
		virtual void Set_Rate(double Rate)=0;
		virtual void Set_dL(double DL)=0;
		virtual double Get_dyield()=0;
		virtual ColumnVector project_back( ColumnVector &Y )=0;
        virtual ColumnVector project_back( ColumnVector &Y ,Matrix &mM)=0;
        virtual ColumnVector project_back_Call_maturity( ColumnVector &Y )=0;
        virtual ColumnVector project_back_Call_dyield( ColumnVector &Y )=0;
        virtual ColumnVector project_back_Call_maturity_start( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Call_DD( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Call_0d( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Put_last( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Put_yydtSI( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Put_yyintdt( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Put_yintydt( ColumnVector &Y )=0;
	    virtual ColumnVector project_back_Put_yintytime( ColumnVector &Y )=0;
		
	};
#endif 
