
 #include "FRP_Pricing.h"
#define OUT_PREC 6
#define COL_WIDTH OUT_PREC +5
#define NUM_STOCK_PRICES 3

#ifdef __cplusplus
extern "C"{
#endif
void PriceMerton( double rate, double ov_rate,double dyield, double sigma, double lambda, double muj, double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ,int lenx, double S0,
                      double *result, double *error,char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend,int n,int Noptions)
	{		    				
		int leny( pow(2.0, iJ));	
		double dividend;
		ndiv==0 ? dividend=0 : dividend=dividends[0];
		for (int i=0;i<Noptions;i++)		   
		   dividend == 0 ? ystart[i] = ystart[i] : ystart[i] = max(ystart[i],log(dividend)+0.01);
		
		
		Real y[leny];	
		double gridstart=*std::min_element(ystart, ystart+Noptions);
        double gridend=*std::max_element(yend, yend+Noptions); 

		seqa(gridstart, gridend, leny, y);//Creating the grid
		
		////////////////////////////////////////////////////
		//PRICE CALCULATION
		Model_1D *MRT;
		MRT=new Merton_FRP(rate,dyield,sigma,lambda,muj,sigmaj,n ,lower,upper);//model
		ColumnVector Z;//Prices	
		if (flag =='P'){
			for (int i=0;i<Noptions;i++){		
				   Z = AmericanPut_1D( K[i], rate,ov_rate,  maturity, dividends,dividends_times,ndiv,x, y, lenx, leny, iJ ,dt,MRT); 
		    for (int i=1;i<=lenx;i++)
		          *result++ = Z(i);
		       
		          
		          }
	          }
			  
		else{
		    for (int i=0;i<Noptions;i++)	{
		       Z = AmericanCall_1D( K[i], rate,  maturity, dividends, dividends_times,ndiv, x, y, lenx, leny, iJ ,dt,MRT);
		    for (int i=1;i<=lenx;i++)
		          *result++ = Z(i);
		         }
			  }
		          
		delete MRT;
		Z.release();

	}
#ifdef __cplusplus
}
#endif



#ifdef __cplusplus
extern "C"{
#endif
void PriceBS(double rate,double ov_rate, double dyield, double sigma, double *K, double maturity,const double *dividends,const double *dividends_times,
             int ndiv,int iJ,int lenx, double S0, double *result, double *error,
              char flag,int lower, int upper,double dt,double *x, double *ystart, double *yend,int Noptions)
	{		
	    int leny=pow(2.0, iJ);	
		double dividend;
		ndiv==0 ? dividend=0 : dividend=dividends[0];
		Real y[leny];	
		for (int i=0;i<Noptions;i++)
		     dividend == 0 ? ystart[i] = ystart[i] : ystart[i] = max(ystart[i],log(dividend)+0.01);	
		double gridstart=*std::min_element(ystart, ystart+Noptions);
        double gridend=*std::max_element(yend, yend+Noptions); 
		seqa(gridstart, gridend, leny, y);//Creating the grid
		////////////////////////////////////////////////////
		//PRICE CALCULATION
		Model_1D *FRP_BS;
		FRP_BS=new FRP_BS_class(rate,dyield,sigma,lower,upper);//model
		
		ColumnVector Z;//Prices	
		if (flag =='P'){
			  for (int i=0;i<Noptions;i++){
			       Z = AmericanPut_1D( K[i], rate,ov_rate,  maturity, dividends,dividends_times,ndiv,x, y, lenx, leny, iJ ,dt,FRP_BS);
		    for (int i=1;i<=lenx;i++)
		          *result++ = Z(i);          
		                
		    }
	          }
		else{
			 for (int i=0;i<Noptions;i++){
		         Z = AmericanCall_1D( K[i], rate,  maturity,dividends, dividends_times,ndiv, x, y, lenx, leny, iJ ,dt,FRP_BS); 
		     for (int i=1;i<=lenx;i++)
		          *result++ = Z(i);
		          }
			  }
		   


		delete FRP_BS;   
		Z.release();

	}
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"{
#endif
void PriceHeston(double rate, double ov_rate,double dyield, double longrun_var_Heston, double mean_rev_Heston, double volvol_Heston, double rho,double lambda,double muj,double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ1,int iJ2, double S0, double V0,
                      double *resultHESTON, double *errorHESTON,char flag,int Nrowstruncation,double dt, double ymin, double ymax,
                      double vmin, double vmax,double gmin,double gmax,double gmin2,double gmax2,int ij,int ij2,int Noptions, double dYShift, char method)
 	{	
		
		//cout<<"Nrowstruncation "<<Nrowstruncation<<endl;
		//cout<<"ymin "<<ymin<<endl;
		//cout<<"ymax "<<ymax<<endl;
		//cout<<"rate "<<rate<<endl;
		//cout<<"ov_rate "<<ov_rate<<endl;
		//cout<<"dyield "<<dyield<<endl;
		//cout<<"longrun_var_Heston "<<longrun_var_Heston<<endl;
		//cout<<"mean_rev_Heston "<<mean_rev_Heston<<endl;
		//cout<<"volvol_Heston "<<volvol_Heston<<endl;
		//cout<<"rho "<<rho<<endl;
		//cout<<"lambda "<<lambda<<endl;
		//cout<<"muj "<<muj<<endl;
		//cout<<"sigmaj "<<sigmaj<<endl;
		//cout<<"maturity "<<maturity<<endl;
		//cout<<"ndiv "<<ndiv<<endl;
		//cout<<"iJ1 "<<iJ1<<endl;
		//cout<<"iJ2 "<<iJ2<<endl;
		//cout<<"S0 "<<S0<<endl;
		//cout<<"V0 "<<V0<<endl;		
						
        int ylen,Wlen,glen,glen2;
        ylen=( pow(2.0, iJ1));
        Wlen =pow(2.0, (double)iJ2);	
        glen=( pow(2.0, ij));
        glen2=( pow(2.0, ij2));
        
        double x0;
        double  yy[ylen], vv[Wlen],grid[glen],grid2[glen2];
        double Delta, Delta2,Delta3,shift;
        
        
        switch(method){
           
           case 'F':
        
        // *******************************************************************
		// *** Initializations and grids
		// *******************************************************************
		//1^ option: set the grid for direct space and deduce the grid for Fourier space
		/////////////////////////////////////////////////////////////////////////////////////////////////////  
        shift=ymin+dYShift;//change the interval to have [0,ymax]
		ymax=ymax-ymin;
		ymin=0;
		
        x0=log(S0)-shift;	
		    
		    ylen=( pow(2.0, iJ1));
		    Wlen = (size_t)pow(2.0, (double)iJ2);	
		//	GRID on y							
			seqa(ymin, ymax, ylen, yy);	
			Delta = yy[2] - yy[1];
		//	GRID on vv									
			seqa(vmin, vmax, Wlen, vv);	
			Delta2 = vv[2] - vv[1];			
		//	Grid on transform of y, u1
			glen = ylen;							 
			//double gmin(-(M_PI*(ylen-1))/(ylen*Delta)+M_PI/(ylen*Delta));	
			gmin=(-(M_PI*(ylen-1))/(ylen*Delta));	
			 gmax=((M_PI*(ylen-1))/(ylen*Delta));	
			//cout<<"gmin= "<<gmin<<endl;
			seqa( gmin, gmax, glen , grid);	
			 Delta3=grid[2]-grid[1];
		//	Grid on transform of vv, u2
			 glen2 = Wlen;				
			 gmin2=(-M_PI*(Wlen-1)/(Wlen*Delta2));	 gmax2=(-gmin2);
			//cout<<"gmin2= "<<gmin2<<endl;		
			seqa( gmin2, gmax2, glen2 , grid2);		
		//2^ option: set the grid for Fourier space and deduce the grid for direct space 
		///////////////////////////////////////////////////////////////////////////////////////////////////// 
		////	Grid on transform of y, u1
			//int glen = ylen;				
			//double grid[glen];
			//double gmin(-800);	double gmax=(-gmin);	
			//seqa( gmin, gmax, glen , grid);		
			//double Delta3=grid[2]-grid[1];
			
		////	Grid on transform of vv, u2
			//int glen2 = Wlen;				
			//double grid2[glen2];
			//double gmin2(-1750);	double gmax2=(-gmin2);		
			//seqa( gmin2, gmax2, glen2 , grid2);	
			//double Delta4=grid2[2]-grid2[1];
		
		////	GRID on y
			//double yy[ylen];								
			//seqa(0, +2*M_PI*(ylen-1)/ylen/Delta3, ylen, yy);	
			//double Delta = yy[2] - yy[1];
			//cout<<Delta<<endl;
		////	GRID on vv
			
			//double vv[Wlen];									
			//seqa(0, +2*M_PI*(Wlen-1)/Wlen/Delta4, Wlen, vv);	
			//double Delta2 = vv[2] - vv[1];		
			//cout<<Delta2<<endl;	
		  break;
		  
		case 'S':
		
			shift=dYShift;
		
            x0=log(S0)-shift;	
		  		
		    ylen=( pow(2.0, iJ1));
		     Wlen = (size_t)pow(2.0, (double)iJ2);	
		//	GRID on y								
			seqa(ymin, ymax, ylen, yy);	
			Delta = yy[2] - yy[1];
		//	GRID on vv										
			seqa(vmin, vmax, Wlen, vv);	
			 Delta2 = vv[2] - vv[1];	
		  //GRID on the transform of y
		   glen=( pow(2.0, ij));
		  seqa( gmin, gmax, glen , grid);	
		  //GRID on the transform of volatility
			glen2 = ( pow(2.0, ij2));						
			seqa( gmin2, gmax2, glen2 , grid2);			  
		  break;	
		}		
		
		
		////////////////////////////////////////////////////
		//PRICE CALCULATION
		//******************************************************/
		//*************
		Heston_FRP *HST;
		HST=new Heston_FRP(rate,dyield,longrun_var_Heston,mean_rev_Heston,
		                                 volvol_Heston,rho, lambda, muj, sigmaj,Nrowstruncation,ylen, Wlen,glen,glen2,
		                                 yy,vv,grid,grid2,Delta,Delta2,shift,method);//model                           
		ColumnVector Z;//Prices	
		if (flag =='P'){
			for (int i=0;i<Noptions;i++){		
				   Z = AmericanPut_2D( K[i], rate,ov_rate, maturity, dividends,dividends_times,
				                                      ndiv,x0,V0, yy, vv, grid, grid2,  ylen, Wlen, iJ1,iJ2 ,dt,HST,shift,Nrowstruncation); 
                      *resultHESTON++=Z.element(0);
		          }
	          }
			  
		else{
		    for (int i=0;i<Noptions;i++)	{
		       Z = AmericanCall_2D( K[i], rate,ov_rate, maturity, dividends, dividends_times,
		                              ndiv, x0,V0,  yy,vv, grid, grid2,  ylen, Wlen, iJ1,
		                              iJ2 ,dt,HST,shift,Nrowstruncation);
		        *resultHESTON++=Z.element(0);
		         }
			  }    
		delete HST;
		Z.release();

	}
#ifdef __cplusplus
}
#endif
