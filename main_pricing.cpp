#include "FRP_Pricing.h"
#define OUT_PREC 6
#define COL_WIDTH OUT_PREC +5
#define NUM_STOCK_PRICES 3

void PriceMerton_Internal( double rate, double ov_rate, double dyield, double sigma, double lambda, double muj, double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ, int lenx,double S0,
                      double *result, double *error, char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend,int n, int Noptions);
                      
 void PriceBS_Internal( double rate, double ov_rate, double dyield, double sigma, double *K, double maturity,const double *dividends,
              const double *dividends_times,int ndiv, int iJ,int lenx, double S0, double *result, double *error,  
              char flag,int lower, int upper,double dt,double *x, double *ystart, double *yend, int Noptions);

void PriceHeston_Internal( double rate, double ov_rate,double dyield, double longrun_var_Heston, double mean_rev_Heston, double volvol_Heston, double rho,double lambda,double muj,double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ1,int iJ2, double S0, double V0,
                      double *resultHESTON, double *errorHESTON,char flag,int Nrowstruncation,double dt,  double ymin, double ymax,
                      double vmin, double vmax,double gmin,double gmax,double gmin2,double gmax2,int ij,int ij2,int Noptions, double dYShift, char method);



int main( int argc, char *argv[] )
{
	///////////////////////////////////////
	/*INPUT PARAMETERS*/
	////////////////////////////////////////

	//Option-stock parameters
	double	maturity( 1), //maturity of the option
			S0(100);//stock value
			
	double	rate(  0.05), ov_rate(   0.05);//rate to maturity and overnight rate
	char option='C';// call ('C') or put ('P') option
	double K[]={100};//strike of the option
	int ndiv=3;//number of discrete dividends
    double dividends_times[]={0.25,0.5,0.75};//time of the payments of the dividends
	double dividends[]={2,2,2};//dividend values. NB: in this version the dividends need to have the same value.
	double dyield(0);//dividend yield

	//Black Scholes Merton model parameters
	double sigma_BS(0.2);//Black-Scholes diffusive volatility
	double lambda( 0);//mean jumps arrival per unit of time
	double muj(0);//mean jump size
	double sigmaj( sqrt(0.05));//std dev of jump
	double sigma(0.2);//Merton diffusive vol
    //Heston model parameters
    //NB: Heston model uses as well the jump parameters of Merton
    //In order to have the usual Heston model (no jump) one should set lambda=0.
    double longrun_var_Heston(0.2*0.2);//Heston long run variance
    double mean_rev_Heston(2);//Heston mean reversion speed
    double volvol_Heston(0.2);//Heston volatility of volatility
    double rho(0);//Heston correlation between stock brownian motion and variance brownian motion
    double V0(0.2*0.2);//Heston variance starting value
 	///////////////////////////////////////
	/*END INPUT PARAMETERS*/
	////////////////////////////////////////	
 	double dt= double(7)/double(365);//time discretization	
    int Noptions=1;
	 	/*1Dimensional discretization*/		
	int iJ(9);//Resolution level	
	int lenx(3);
	Real x[]= {log(0.8*S0), log(S0), log(1.2*S0) };
	int lower;
	int upper;   
	if (option=='P'){
	    lower=max(double(0),pow(2,(iJ-10))*200);
	    upper=max(double(0),pow(2,(iJ-10))*200);}
	else{
		lower=max(double(0),pow(2,(iJ-10))*300);
	    upper=max(double(0),pow(2,(iJ-10))*300);
	}
	 

    
    char method='F'; 
    //double longrun_vol_Heston(0.2), mean_rev_Heston(2), volvol_Heston(0.2), rho(0), V0(0.2*0.2);
    
    /*2Dimensional discretization*/
    int iJ1(8), iJ2(7);	
	double ymin(log(K[0]) -15*sqrt(V0)*sqrt(maturity)), ymax(log(K[0])+15*sqrt(V0)*sqrt(maturity));		
	double vmin( 0 ), vmax( 0.5);	
	
	double gmin(-100),gmin2(-400);
	double gmax(-gmin),gmax2(-gmin2);
	int ij(8),ij2(7);
	
	int ylen = (size_t)pow(2.0, (double)iJ1);
	int NRowstruncation=ylen/2;
    double dYshift;
    dYshift=0;//dYshift=yminreal-ymin;
	//Declaring the vectors to store the results and utils
	double total_var_log,max_mean_log,min_mean_log;
    double ystart[Noptions];
    double yend[Noptions];		
    double resultHESTON[NUM_STOCK_PRICES*Noptions],errorHESTON[NUM_STOCK_PRICES*Noptions];
    double resultMRT[NUM_STOCK_PRICES*Noptions],resultBS[NUM_STOCK_PRICES*Noptions];
    double errorMRT[NUM_STOCK_PRICES*Noptions],errorBS[NUM_STOCK_PRICES*Noptions];	
    double t,t2,t_init;
	


try{
	///*///////////////////////*/
	///* HESTON */
	///*///////////////////////*/
	///*Varianza equivalente di Black-Scholes*/
	//double beta=mean_rev_Heston;
   //double alpha=mean_rev_Heston*longrun_var_Heston;
   //double fun=alpha*(volvol_Heston*volvol_Heston+beta*volvol_Heston*(volvol_Heston*maturity-double(4)*rho)+double(2)*beta*beta*(double(1)-rho*volvol_Heston*maturity))+V0*beta*beta*(volvol_Heston*(double(2)*rho-volvol_Heston*maturity)+double(2)*beta*(rho*volvol_Heston*maturity-double(1)));
   //double g=double(2)*V0*beta*(double(4)*beta*(beta-rho*volvol_Heston)+volvol_Heston*volvol_Heston)+alpha*(double(8)*beta*beta*beta*maturity-double(5)*volvol_Heston*volvol_Heston+double(2)*beta*volvol_Heston*(double(8)*rho+volvol_Heston*maturity)-double(8)*beta*beta*(double(1)+rho*volvol_Heston*maturity));
   //double Var=exp(-double(2)*beta*maturity)/(maturity*double(8)*beta*beta*beta*beta)*((alpha-double(2)*V0*beta)*volvol_Heston*volvol_Heston+double(4)*exp(beta*maturity)*fun+exp(double(2)*beta*maturity)*g);
	
    t = clock();
    PriceHeston_Internal(rate, ov_rate,dyield,longrun_var_Heston, mean_rev_Heston, volvol_Heston,rho,lambda,muj,sigmaj,
                     K, maturity,dividends,dividends_times,ndiv,  iJ1,iJ2, S0,  V0,
                      resultHESTON, errorHESTON,option,NRowstruncation,dt,  ymin, ymax,vmin,
                       vmax, gmin, gmax,gmin2,gmax2, ij,ij2,Noptions,dYshift,method);
	t2 = clock();
	cout<< "Seconds elapsed end Heston: ";
	cout<<fixed<<setprecision(6)<<(t2-t)/CLOCKS_PER_SEC<<endl;   
	cout<< "Heston price: "<<resultHESTON[0]<<endl<<endl;
	
	//cout<<"NII: "<<(S0-K[0])-resultHESTON[0]<<endl;
	
	///*///////////////////////*/
	///* MERTON */
	///*///////////////////////*/
	double kappa=exp(muj+sigmaj*sigmaj/2)-1;//expected value of the jump size	
	total_var_log=sigma*sigma+lambda*(muj*muj+sigmaj*sigmaj); 	

	max_mean_log=x[lenx-1]+(rate-dyield -sigma*sigma/2-lambda*kappa)*maturity+lambda*muj*maturity;
	min_mean_log=x[0]+(rate-dyield -sigma*sigma/2-lambda*kappa)*maturity+lambda*muj*maturity;

	for (int i=0;i<Noptions;i++){
		 ystart[i]=min(min_mean_log-5*sqrt(total_var_log),log(K[i])-10*sigma);//Discretization interval for log(S)
		 yend[i]= max(max_mean_log+5*sqrt(total_var_log),log(K[i])+10*sigma);
		 }
	int n=int(10*lambda);	
	t_init = clock();
    PriceMerton_Internal(rate,ov_rate,dyield,sigma,lambda,muj,sigmaj,K,maturity,dividends, 
             dividends_times,ndiv, iJ,lenx, S0, resultMRT, errorMRT,option,lower,upper,dt,x,ystart,yend,n,Noptions);       
		t = clock();
		cout<< "Seconds elapsed end Merton: ";
		cout<<fixed<<setprecision(3)<<(t-t_init)/CLOCKS_PER_SEC<<endl;   
		cout<< "Merton price: "<<resultMRT[1]<<endl<<endl;
		//cout<<resultMRT[0]<<" "<<resultMRT[1]<<" "<<resultMRT[2]<<endl;//" "<<resultMRT[4]<<" "<<resultMRT[7]<<" "<<resultMRT[10]<<" "<<resultMRT[13]<<endl;  
	/////*///////////////////////*/
	/////* BLACK-SCHOLES */
	/////*///////////////////////*/    
    
    total_var_log=sigma_BS*sigma_BS;
	max_mean_log=x[lenx-1]+(rate-dyield -sigma_BS*sigma_BS/2)*maturity;
	min_mean_log=x[0]+(rate -dyield -sigma_BS*sigma_BS/2)*maturity;
	for (int i=0;i<Noptions;i++){
		ystart[i]=min(min_mean_log-5*sqrt(total_var_log),log(K[i])-10*sigma_BS);//Discretization interval for log(S)
		yend[i]=max(max_mean_log+5*sqrt(total_var_log),log(K[i])+10*sigma_BS);}
    t=clock();
	    PriceBS_Internal( rate,ov_rate,dyield, sigma_BS, K, maturity, dividends,
             dividends_times,ndiv, iJ,lenx, S0, resultBS, errorBS,option,lower,upper,dt,x,ystart,yend,Noptions);
     t2 = clock();
		cout<< "Seconds elapsed end BS: ";
		cout<<fixed<<setprecision(3)<<(t2-t)/CLOCKS_PER_SEC<<endl;  
		cout<< "Black-Scholes price: "<<resultBS[1]<<endl; 
		//cout<<setprecision(6)<<resultBS[0]<<" "<<resultBS[1]<<" "<<resultBS[2]<<endl;//<<"	"<<resultBS[4]<<" "<<resultBS[7]<<" "<<resultBS[10]<<" "<<resultBS[13]<<" "<<resultBS[16]<<" "<<resultBS[19]<<" "<<resultBS[22]<<endl;  
        }
      catch(BaseException) { cout << BaseException::what() << endl; }
      catch(exception& e){cout << "Error: standard exception " << e.what() << endl;}
      catch(...){cout << "Error: non standard exception " <<  endl;}
    
	return 0; 	
}

void PriceMerton_Internal( double rate, double ov_rate,double dyield, double sigma, double lambda, double muj, double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ,int lenx, double S0,
                      double *result, double *error,char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend,int n,
                      int Noptions)
	{
		cout<<"Starting the MRT Pricing"<<endl;	

        PriceMerton(rate,ov_rate,dyield,sigma,lambda,muj,sigmaj,K,maturity,dividends, 
             dividends_times,ndiv, iJ,lenx, S0, result, error,flag,lower,upper,dt,x,ystart,yend,n,Noptions); 

		

	   
	
	}


void PriceBS_Internal(double rate,double ov_rate, double dyield, double sigma, double *K, double maturity,const double *dividends,const double *dividends_times,
             int ndiv,int iJ,int lenx, double S0, double *result, double *error,
              char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend, int Noptions)

	{		

	    cout<<"Starting the BS Pricing"<<endl;
	    PriceBS( rate,ov_rate,dyield, sigma, K, maturity, dividends, 
             dividends_times,ndiv, iJ,lenx, S0, result, error,flag,lower,upper,dt,x,ystart,yend,Noptions);				


		  
	   

	}
void PriceHeston_Internal( double rate, double ov_rate,double dyield, double longrun_var_Heston, double mean_rev_Heston, double volvol_Heston, double rho,double lambda,double muj,double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ1,int iJ2, double S0, double V0,
                      double *resultHESTON, double *errorHESTON,char flag,int Nrowstruncation,double dt,double ymin, double ymax,
                      double vmin, double vmax,double gmin,double gmax,double gmin2,double gmax2,int ij,int ij2,int Noptions, double dYShift,char method)
	{
	cout<<"Starting the Heston Pricing"<<endl;	
    PriceHeston(rate, ov_rate,dyield,longrun_var_Heston, mean_rev_Heston, volvol_Heston,rho, lambda, muj, sigmaj,
                     K, maturity,dividends,dividends_times,ndiv,  iJ1,iJ2,S0,  V0,
                      resultHESTON, errorHESTON,flag,Nrowstruncation,dt, ymin, ymax,vmin,
                       vmax, gmin, gmax, gmin2, gmax2, ij, ij2,Noptions,dYShift,method);
}


