
#include "FRP_BS_class.h"


FRP_BS_class::FRP_BS_class(const Real rate,//rate
                  const Real dyield,
				  const Real sigma, //vol				  
				  const int lower,
				  const int upper):
 Model_1D(),
 rate(rate),dyield(dyield), sigma(sigma), lower(lower),upper(upper)
 {
	 	is_Call_DD_mG=0;
	 	is_Call_0d_mG=0;
	 	is_Call_mat_mG=0;	 
	 	is_Call_dyield_mG=0;	
	    is_Put_last_mG=0;
	    is_Put_yydtSI_mG=0;
	    is_Put_yyintdt_mG=0;
	    is_Put_yintydt_mG=0;
	    is_Put_yintytime_mG=0;
        if (lower!=0 || upper!=0)
	        isBand_mat=1;
	    else
	         isBand_mat=0;
}
 
 FRP_BS_class::~FRP_BS_class()
 {
	 mG.release();
	 bG.release();
	 Call_maturity_mG.release();
	 Call_DD_mG.release();
	 Call_0d_mG.release();
     Call_maturity_bG.release();
     Call_DD_bG.release();
     Call_0d_bG.release();
     Put_last_mG.release();
     Put_yydtSI_mG.release();
     Put_yyintdt_mG.release();
     Put_yintydt_mG.release();
     Put_yintytime_mG.release();
     Put_last_bG.release();
     Put_yydtSI_bG.release();
     Put_yyintdt_bG.release();
     Put_yintydt_bG.release();
     Put_yintytime_bG.release();     
	  }
 
  void FRP_BS_class::Make_Green_Call_Maturity( const int lenx,//number of starting points
                              const int leny, //number of arrivals  points
							  Real *xx, //grid of starting points
							  Real *yy, //grid of arrivals points
							  const int J,//precision
							  Real t)//time step
{
	if (!is_Call_mat_mG){
	Make_Green(lenx,leny,xx,yy,J,t);
	if (isBand_mat && lenx==leny){
		Call_maturity_bG.resize(lenx,lower,upper);Call_maturity_bG=0;
	    Call_maturity_bG=bG;}
	else{
		Call_maturity_mG.resize(lenx,leny);Call_maturity_mG=0;
	    Call_maturity_mG=mG;}
	 is_Call_mat_mG=1;
	}
}
 
 void FRP_BS_class::Make_Green_Call_dyield( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Call_dyield_mG){
	Make_Green_SpaceInvariance(lenx,leny,xx,yy,J,t);
	if (isBand_mat && lenx==leny){
		Call_dyield_bG.resize(lenx,lower,upper);Call_dyield_bG=0;
	    Call_dyield_bG=bG;}
	else{
		Call_dyield_mG.resize(lenx,leny);Call_dyield_mG=0;
	    Call_dyield_mG=mG;}
	 is_Call_dyield_mG=1;
	}	
}
 
 
void FRP_BS_class::Make_Green_Call_DD( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
						  
{
	if (!is_Call_DD_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
		  Call_DD_bG.resize(lenx,lower,upper);Call_DD_bG=0;
	       Call_DD_bG=bG;}
	      else{
			Call_DD_mG.resize(lenx,leny);Call_DD_mG=0;
	       Call_DD_mG=mG;}
	       
	       is_Call_DD_mG=1;
		   
	   }
	   
}
						  
void FRP_BS_class::Make_Green_Call_0d( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Call_0d_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Call_0d_bG.resize(lenx,lower, upper);Call_0d_bG=0;   
	       Call_0d_bG=bG;}
	      else{
			  Call_0d_mG.resize(lenx,leny);Call_0d_mG=0;   
	       Call_0d_mG=mG;}
		  is_Call_0d_mG=1;
	   }	
}

 void FRP_BS_class::Get_Green_Call_mat(Matrix & mM) const  
{
	if(isBand_mat)
	mM = Call_maturity_bG;
	else
	mM = Call_maturity_mG;
	}

void FRP_BS_class::Get_Green_Call_DD(Matrix & mM) const  
{
	if(isBand_mat)
	mM = Call_DD_bG;
	else
	mM = Call_DD_mG;
	}

void FRP_BS_class::Get_Green_Call_0d(Matrix & mM) const  
{	
	mM = Call_0d_mG;
	}
 
 void FRP_BS_class::Make_Green_Put_last( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Put_last_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Put_last_bG.resize(lenx,lower,upper);Put_last_bG=0;   
	       Put_last_bG=bG;}
	      else{
			  Put_last_mG.resize(lenx,leny);Put_last_mG=0;   
	       Put_last_mG=mG;}
		  is_Put_last_mG=1;
	   }	
}

void FRP_BS_class::Make_Green_Put_yydtSI( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Put_yydtSI_mG){
		  Make_Green_SpaceInvariance(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Put_yydtSI_bG.resize(lenx,lower,upper);Put_yydtSI_bG=0;   
	       Put_yydtSI_bG=bG;}
	      else{
			  Put_yydtSI_mG.resize(lenx,leny);Put_yydtSI_mG=0;   
	       Put_yydtSI_mG=mG;}
		  is_Put_yydtSI_mG=1;
	   }	
}
						  
						  
						  
void FRP_BS_class::Make_Green_Put_yyintdt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Put_yyintdt_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Put_yyintdt_bG.resize(lenx,lower,upper);Put_yyintdt_bG=0;   
	       Put_yyintdt_bG=bG;}
	      else{
			  Put_yyintdt_mG.resize(lenx,leny);Put_yyintdt_mG=0;   
	       Put_yyintdt_mG=mG;}
		  is_Put_yyintdt_mG=1;
	   }	
}
void FRP_BS_class::Make_Green_Put_yintydt( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
						 
{
	if (!is_Put_yintydt_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Put_yintydt_bG.resize(lenx,lower,upper);Put_yintydt_bG=0;   
	       Put_yintydt_bG=bG;}
	      else{
			  Put_yintydt_mG.resize(lenx,leny);Put_yintydt_mG=0;   
	       Put_yintydt_mG=mG;}
		  is_Put_yintydt_mG=1;
	   }	
}

void FRP_BS_class::Make_Green_Put_yintytime( const int lenx, const int leny, 
						  Real *xx, Real *yy, const int J,Real t)
{
	if (!is_Put_yintytime_mG){
		  Make_Green(lenx,leny,xx,yy,J,t);
		   if (isBand_mat && lenx==leny){
			 Put_yintytime_bG.resize(lenx,lower,upper);Put_yintytime_bG=0;   
	       Put_yintytime_bG=bG;}
	      else{
			  Put_yintytime_mG.resize(lenx,leny);Put_yintytime_mG=0;   
	       Put_yintytime_mG=mG;}
		  is_Put_yintytime_mG=1;
	   }	
}						  
void FRP_BS_class::Get_Green_Put_last(Matrix &mM) const
{
mM = Put_last_mG;	
}		  
void FRP_BS_class::Get_Green_Put_yydtSI(Matrix &mM) const
{
	if(isBand_mat)
	mM = Put_yydtSI_bG;
	else
	mM = Put_yydtSI_mG;
}
void FRP_BS_class::Get_Green_Put_yyintdt(Matrix &mM) const
{
	if(isBand_mat)
	mM = Put_yyintdt_bG;
	else
	mM = Put_yyintdt_mG;
}
void FRP_BS_class::Get_Green_Put_yintydt(Matrix &mM) const
{
	if(isBand_mat)
	mM = Put_yintydt_bG;
	else
	mM = Put_yintydt_mG;
}
void FRP_BS_class::Get_Green_Put_yintytime(Matrix &mM) const
{
	if(isBand_mat)
	mM = Put_yintytime_bG;
	else
	mM = Put_yintytime_mG;
}

 
 
 void FRP_BS_class::Make_Green_SpaceInvariance( const int lenx,//number of starting points
                              const int leny, //number of arrivals  points
							  Real *xx, //grid of starting points
							  Real *yy, //grid of arrivals points
							  const int J,//precision
							  Real t)//time step
{
	//This function computes the discounted green function as a matrix m X k,
	//by exploiting the space translation invariance that occurs when the vectors xx 
	//and yy are the same equally spaced vectors
	double time_to_mat=t;
	double s_t = sqrt(time_to_mat);
	double xydiff;
	mG.resize(lenx, leny);
	mG =0; //initialize the Green matrix
	if (isBand_mat && lenx==leny)
    {
		bG.resize(leny,lower, upper); bG =0; 
		
		////////////////////////
		//Compute the first row of the green matrix
		////////////////////////////

		int v;//row index
		v=1;				
			for(int i=1; i<=min(upper+1,leny); ++i)//column index
			{
				xydiff=xx[v-1]-yy[i-1];
				bG(v, i) = exp(-rate * time_to_mat) * (n( (xydiff + (rate -dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) )* dL;
				
				for (int j=i+1;j<=leny;++j)//uses space translation invariance
					 bG(j-i+1,j)=bG(v,i);
			}
		////////////////////////
		//Compute the first column of the green matrix
		////////////////////////////
		int ii;//column index
		ii=1;				
			for(int vv=2; vv<=min(lower+1,lenx); ++vv)//row index
			{
				
				xydiff=xx[vv-1]-yy[ii-1];

				bG(vv,ii)= exp(-rate * time_to_mat) * n( (xydiff + (rate -dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) * dL;
				for (int j=vv+1;j<=lenx;++j)//uses space translation invariance
					 bG(j,j-vv+1)=bG(vv,ii);
			  }
	}
	else
	{
	////////////////////////
    //Compute the first row of the green matrix
    ////////////////////////////

    int v;//row index
    v=1;				
		for(int i=1; i<=leny; ++i)//column index
		{
			
            xydiff=xx[v-1]-yy[i-1];
		    mG(v,i)= exp(-rate * time_to_mat) * n( (xydiff + (rate-dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) * dL;
		    for (int j=i+1;j<=leny;++j)//uses space translation invariance
		         mG(j-i+1,j)=mG(v,i);

		}
	////////////////////////
    //Compute the first column of the green matrix
    ////////////////////////////
	int ii;//column index
    ii=1;				
		for(int vv=2; vv<=lenx; ++vv)//row index
		{
			
            xydiff=xx[vv-1]-yy[ii-1];
		    mG(vv,ii)= exp(-rate * time_to_mat) * n( (xydiff + (rate-dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) * dL;
		    for (int j=vv+1;j<=lenx;++j)//uses space translation invariance
		         mG(j,j-vv+1)=mG(vv,ii);
	      }
	
	if ( ! mG.IsZero() )
		isGreen_mat = 1;
  
}
		////Uncomment if you want to print the green matrix
		//ofstream outfile;
		//outfile.open ("Green_Matrix.txt");
        //outfile<<mG;
        //outfile.close();
		
}

 void FRP_BS_class::Make_Green( const int lenx,//number of starting points
                              const int leny, //number of arrivals  points
							  Real *xx, //grid of starting points
							  Real *yy, //grid of arrivals points
							  const int J,//precision
							  Real t)//time step
{
	//This function computes the discounted green function as a matrix m X k,

	
	double time_to_mat=t;
	double s_t = sqrt(time_to_mat);
	double xydiff;
	mG.resize(lenx, leny); mG =0; //initialize the Green matrix

    if (isBand_mat && lenx==leny)
    {
	bG.resize(leny,lower, upper); bG =0; 
	int indexstart, indexend;
	for(int v=1; v<=lenx; ++v )	//Build 'direct' Green matrix
		{	
			indexstart=max(1,v-lower);
			indexend=min(leny,v+upper);		
			for(int i=indexstart; i<=indexend; ++i)
			{
				
				xydiff=xx[v-1]-yy[i-1];
				bG(v,i)= exp(-rate * time_to_mat) * n( (xydiff + (rate-dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) * dL;
				//if (v==1 && i==1){
				    //cout<<bG(v,i)<<endl;
				    //cout<<disc<<endl;
				    //cout<<temp<<endl;
				    //cout<<dL<<endl;
				    //cout<<coeff<<endl;
				    //}
			}
		}
		////Uncomment if you want to print the green matrix
		//ofstream outfile;
		//outfile.open ("Green_Matrix.txt");
        //outfile<<bG;
        //outfile.close();	
    }	
	else
	{
		for(int v=1; v<=lenx; ++v )	//Build 'direct' Green matrix
		{				
			for(int i=1; i<=leny; ++i)
			{
				
				xydiff=xx[v-1]-yy[i-1];
				mG(v,i)= exp(-rate * time_to_mat) * n( (xydiff + (rate-dyield - 0.5*sigma*sigma)*time_to_mat ) / (s_t * sigma ) ) / (s_t*sigma) * dL;
			}
		}
		
		if ( ! mG.IsZero() )
			isGreen_mat = 1;	

    }	
}

void FRP_BS_class::Get_Green(Matrix & mM) const  
{
	if(isBand_mat)
	mM=bG;
	else
	mM = mG;
	 }

bool FRP_BS_class::Get_IsBand_mat() const 
{
	return isBand_mat;
}

bool FRP_BS_class::Get_IsGreen_mat() const {return isGreen_mat;}


ColumnVector FRP_BS_class::project_back( ColumnVector &vY )
// projects back the value function
{
	if(isBand_mat)
	   vY=bG*vY;
	else
	{	
		
		if(isGreen_mat)
			vY = mG * vY;
		else
		{
			cout << "Make Green function first" << endl;
			terminate();
		}
     }
	return vY;
	
}

ColumnVector FRP_BS_class::project_back_Call_maturity( ColumnVector &Y )
{
	if(isBand_mat)
	   Y=Call_maturity_bG*Y;
	else
	   Y = Call_maturity_mG* Y;

  	return Y;
	
}

ColumnVector FRP_BS_class::project_back_Call_dyield( ColumnVector &Y )
{
	if(isBand_mat)
	   Y=Call_dyield_bG*Y;
	else
	   Y = Call_dyield_mG* Y;

  	return Y;
	
}

ColumnVector FRP_BS_class::project_back_Call_maturity_start( ColumnVector &Y )
{

	   Y = Call_maturity_mG* Y;

  	return Y;
	
}

	    ColumnVector FRP_BS_class::project_back_Call_DD( ColumnVector &Y )
{
	if(isBand_mat)
	   Y=Call_DD_bG*Y;
	else
	   Y = Call_DD_mG* Y;

  	return Y;
	
}
 ColumnVector FRP_BS_class::project_back_Call_0d( ColumnVector &Y )
 {
	   Y = Call_0d_mG* Y;

  	return Y;
	
}
ColumnVector FRP_BS_class::project_back_Put_last( ColumnVector &Y )
 {

	   Y = Put_last_mG* Y;

  	return Y;
	
}
	 ColumnVector FRP_BS_class::project_back_Put_yydtSI( ColumnVector &Y )
 {
	if(isBand_mat)
	   Y=Put_yydtSI_bG*Y;
	else
	   Y = Put_yydtSI_mG* Y;

  	return Y;
	
}	 
	 ColumnVector FRP_BS_class::project_back_Put_yyintdt( ColumnVector &Y )
 {

	if(isBand_mat)
	   Y=Put_yyintdt_bG*Y;
	else
	   Y = Put_yyintdt_mG* Y;

  	return Y;
	
}	     
	     ColumnVector FRP_BS_class::project_back_Put_yintydt( ColumnVector &Y )
 {
	if(isBand_mat)
	   Y=Put_yintydt_bG*Y;
	else
	   Y = Put_yintydt_mG* Y;

  	return Y;
	
}	      
	     ColumnVector FRP_BS_class::project_back_Put_yintytime( ColumnVector &Y )
	     
 {
	if(isBand_mat)
	   Y=Put_yintytime_bG*Y;
	else
	   Y = Put_yintytime_mG* Y;

  	return Y;
	
}	 

 


ColumnVector FRP_BS_class::project_back( ColumnVector &vY ,Matrix &mM)

{
	vY = mM * vY;
	return vY;
}

void FRP_BS_class::Set_Rate(double Rate){rate=Rate;}
void FRP_BS_class::Set_dL(double DL){dL=DL;}
double FRP_BS_class::Get_dyield(){return dyield;}

