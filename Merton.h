#ifdef __cplusplus
extern "C"{
#endif
void PriceMerton(int test, double rate,double ov_rate, double dyield, double sigma, double lambda, double muj, double sigmaj,
                     double *K, double maturity,const double *dividends,const double *dividends_times,int ndiv, int iJ, int lenx, double S0,
                      double *result, double *error, char flag,int lower, int upper,double dt, double *x, double *ystart, double *yend,int n,int Noptions);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"{
#endif                    
void PriceBS( double rate,double ov_rate, double dyield, double sigma, double *K, double maturity,const double *dividends,
              const double *dividends_times,int ndiv, int iJ, int lenx, double S0, double *result, double *error,  
              char flag,int lower, int upper,double dt,double *x, double *ystart, double *yend,int Noptions);

#ifdef __cplusplus
}
#endif
