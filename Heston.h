#ifdef __cplusplus
extern "C"{
#endif
void PriceMerton(int test, double rate, double ov_rate, double dyield, double sigma, double lambda, double muj, double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ, int lenx,double S0,
                      double *result, double *error, char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend,int n,int Noptions);
#ifdef __cplusplus
}
#endif  

#ifdef __cplusplus
extern "C"{
#endif
void PriceHeston(double rate, double ov_rate,double dyield, double longrun_var_Heston, double mean_rev_Heston, double volvol_Heston, double rho,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ1,int iJ2,int lenx, double S0, double V0,
                      double *resultHESTON, double *errorHESTON,char flag,int Nrowstruncation,double dt, double *x, double ymin, double ymax,
                      double vmin, double vmax,int Noptions, double *dYShift);
#ifdef __cplusplus
}
#endif  

                   
#ifdef __cplusplus
extern "C"{
#endif                     
void PriceBS( double rate, double ov_rate, double dyield, double sigma, double *K, double maturity,const double *dividends,
              const double *dividends_times,int ndiv, int iJ, int lenx, double S0, double *result, double *error,  
              char flag,int lower, int upper,double dt,double *x, double *ystart, double *yend, int Noptions);            
              
#ifdef __cplusplus
}
#endif
