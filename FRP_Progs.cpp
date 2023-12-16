
#include "FRP_Progs.h"
#include "FRP_Merton_class.h"
#include "FRP_Heston_class.h"

double payoff_call( const double &s, const double &k ){ return max(s-k, 0.0); };
double payoff_put ( const double &s, const double &k ){ return max(k-s, 0.0); };
double payoff_binary_call(const double& S, const double& K){
    if (S>=K) return 1;
    return 0;
};
double payoff_asset_call(const double& S, const double& K){
    if (S>=K) return 1;
    return 0;
};

double payoff_binary_put(const double& S, const double& K){
    if (S<=K) return 1;
    return 0;
};

ReturnMatrix AmericanCall_1D(const double &K, 
                          const double &rate, 
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  Real *x, Real *y, 
						  const int &lenx, 
						  const int &leny, 
						  const int &iJ ,
						  const double dt,
						  Model_1D *MDL)
{
	double dL = y[1] - y[0];
	MDL->Set_dL(dL);
	ColumnVector vP;
	int numdiv=ndiv;
	int i1=numdiv-1;
	while (dividends_times[i1]>maturity && i1>-1)
	i1-=1;
	numdiv=i1+1;
	ColumnVector vZ(leny),//Continuation value
	             vV(leny); //Intrinsic value
    /*PAYOFF*/
	
	for(int it=0; it != leny; ++it ) 
		vV(it+1) = payoff_call( exp( y[it] ), K ) ;//payoff valutato a maturity
		
	vZ=vV;
	double dyield=MDL->Get_dyield();
	int nsteps;
	if (dyield!=0)
	{
		maturity>dt? nsteps=int(maturity/dt) : nsteps=1;
		MDL->Make_Green_Call_dyield(leny, leny, y, y, iJ, maturity/nsteps);
		for (int i=nsteps; i>1;i--)
		{
	     vZ = MDL->project_back_Call_dyield(vZ) ;	// Continuation value
         for(int it=0; it != leny; ++it )
	       vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;  			
		}
				//Last step 1->0
		MDL->Make_Green_Call_0d( lenx, leny, x, y, iJ, maturity/nsteps);//Build green matrix for the last step
		vP.resize(lenx);
		vP=MDL->project_back_Call_0d(vZ);	
		vP.release();
		vZ.release();
		vV.release();
		return vP;	
   }	
	
	if (numdiv==0)
	{
		MDL->Make_Green_Call_Maturity( lenx, leny, x, y, iJ, maturity);
		vP.resize(lenx);
		vP=MDL->project_back_Call_maturity_start(vZ);	
		vP.release(); return vP;
    }
 	double yint[leny];//grid of the log(S_t) after dividend payment	
	for(int i=0; i<leny; ++i)	
	    yint[i] = log( exp(y[i]) - dividends[0] );   


    //T->T-1

    MDL->Make_Green_Call_Maturity( leny, leny, yint, y, iJ, maturity-dividends_times[numdiv-1]);//Build green matrix
    vZ = MDL->project_back_Call_maturity(vZ) ;	// Continuation value
    for(int it=0; it != leny; ++it )
	    vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;
	    
	//T-1->T-2,...2->1   

	int flag=0;
	for (int it=numdiv-1;it>0;it--)
	{
	    if (flag==0)
	    {
	        MDL->Make_Green_Call_DD( leny, leny, yint, y, iJ, dividends_times[it]-dividends_times[it-1]);//Build green matrix
	        flag=1;
	     }

	   vZ = MDL->project_back_Call_DD(vZ) ;	// Continuation value

       for(int it=0; it != leny; ++it )
	      vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;   
	}
	//Last step 1->0
	MDL->Make_Green_Call_0d( lenx, leny, x, y, iJ, dividends_times[0]);//Build green matrix for the last step
	vP.resize(lenx);
	vP=MDL->project_back_Call_0d(vZ);	
	vP.release();
	vZ.release();
	vV.release();
	return vP;	
	
}	

ReturnMatrix AmericanPut_1D(const double &K, 
                          const double &rate, 
                          const double &ov_rate,
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  Real *x, Real *y, 
						  const int &lenx, 
						  const int &leny, 
						  const int &iJ ,
						  const double dt,
						  Model_1D *MDL)
{	
    double Tt;
    int flag=0;
	ColumnVector vP;
	int numdiv=ndiv;
	int i1=numdiv-1;
	while (dividends_times[i1]>maturity && i1>-1)
	i1-=1;
	numdiv=i1+1;
	ColumnVector vZ(leny),//Continuation value
	             vV(leny), //Intrinsic value
	             dvV(leny);//Intrinsic value at ex-dividend date
    /*PAYOFF*/
	double dL = y[1] - y[0];
	MDL->Set_dL(dL);
	for(int it=0; it != leny; ++it ) 
		vV(it+1) = payoff_put( exp( y[it] ), K ) ;//payoff valutato a maturity

	vZ=vV;	
	Matrix mG;
	int nsteps;
	if (numdiv==0)
	{
		maturity>dt? nsteps=int(maturity/dt) : nsteps=1;
		if (nsteps>1){
		MDL->Make_Green_Put_yydtSI( leny, leny, y, y, iJ, maturity/nsteps);}
		  
		for (int i=nsteps; i>1;i--)
		{
	     vZ = MDL->project_back_Put_yydtSI(vZ) ;	// Continuation value
         for(int it=0; it != leny; ++it )
	       vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;  			
		}
		MDL->Set_Rate(ov_rate);
		MDL->Make_Green_Put_last( lenx, leny, x, y, iJ, maturity/nsteps);//Build green matrix for the last step		
		vP.resize(lenx);
		vP=MDL->project_back_Put_last(vZ);	
		//~ for(int it=0; it != lenx; ++it )
	       //~ vP(it+1) = max( payoff_put( exp( x[it] ), K ) , vP(it+1) ) ;  
		vP.release(); 
		vZ.release();
	    vV.release();
	    dvV.release();
		return vP;	
    }	
	double yint[leny];//grid of the log(S_t) after dividend payment	
	for(int i=0; i<leny; ++i)	
	    yint[i] = log( exp(y[i]) - dividends[0] );  
	for(int it=0; it != leny; ++it ) 
	    dvV(it+1) = payoff_put( exp( yint[it] ), K ) ;	//payoff valutato at ex-dividends dates
	if ((maturity-dividends_times[numdiv-1])<dt)
	    {
         MDL->Make_Green_Put_yintytime( leny, leny, yint, y, iJ, (maturity-dividends_times[numdiv-1]));//Build green matrix
	    
	    vZ = MDL->project_back_Put_yintytime(vZ) ;
	    for(int it=0; it != leny; ++it )
	       vZ(it+1) = max( dvV(it + 1), vZ(it+1) ) ;//ex-dividend date     
	     
	     }
	 else
	 {
		 nsteps=int((maturity-dividends_times[numdiv-1])/dt);
		 Tt=(maturity-dividends_times[numdiv-1])/nsteps;
		 
		 if (nsteps==1)
		    {
	        MDL->Make_Green_Put_yintytime( leny, leny, yint, y, iJ, (maturity-dividends_times[numdiv-1]));//Build green matrix
	        vZ = MDL->project_back_Put_yintytime(vZ) ;
	        for(int it=0; it != leny; ++it )
	            vZ(it+1) = max( dvV(it + 1), vZ(it+1) ) ;//ex-dividend date  				
			}
		else
		    { 	
			 MDL->Make_Green_Put_yydtSI( leny, leny, y, y, iJ,Tt);//Build green matrix for the days in which there are no dividends
		    for (int i=nsteps;i>1;i--){
			   vZ = MDL->project_back_Put_yydtSI(vZ) ;	// Continuation value
               for(int it=0; it != leny; ++it )
	               vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;   
			    }
		    MDL->Make_Green_Put_yintydt( leny, leny, yint, y, iJ, Tt);//Build green matrix

	        vZ = MDL->project_back_Put_yintydt(vZ) ;
	        
	        for(int it=0; it != leny; ++it )
	            vZ(it+1) = max( dvV(it + 1), vZ(it+1) ) ;//ex-dividend date  		   
		    }
	 }
	///////////////////////////////////////////////////
	//Codice per il caso multidividend
	
	if (numdiv>1)
	{
      nsteps=int((dividends_times[numdiv-1]-dividends_times[numdiv-2])/dt); 
	  Tt=(dividends_times[numdiv-1]-dividends_times[numdiv-2])/nsteps;		
	  MDL->Make_Green_Put_yintydt( leny, leny, yint, y, iJ, Tt);
	  MDL->Make_Green_Put_yydtSI( leny, leny, y, y, iJ,Tt);
	for (int it=numdiv;it>1;it--)
		  {   
		   vZ = MDL->project_back_Put_yydtSI(vZ) ;
		   for(int it=0; it != leny; ++it )
			   vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ; 

		   for (int j=1; j<nsteps-1;j++)
			   {
				vZ = MDL->project_back_Put_yydtSI(vZ) ;	// Continuation value
				for(int it=0; it != leny; ++it )
					vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;  			   
			   }

		   vZ = MDL->project_back_Put_yintydt(vZ) ;
		   for(int it=0; it != leny; ++it )
			   vZ(it+1) = max( dvV(it + 1), vZ(it+1) ) ; 
		  }
      }

	
	/////////////////////////////////////////////////
	if (dividends_times[0]<dt)
	    {
	    MDL->Set_Rate(ov_rate);
	    MDL->Make_Green_Put_last( lenx, leny, x, y, iJ, dividends_times[0]);//Build green matrix
	    vP.resize(lenx);
		vP=MDL->project_back_Put_last(vZ);	
		//~ for(int it=0; it != lenx; ++it )
	       //~ vP(it+1) = max( payoff_put( exp( x[it] ), K ) , vP(it+1) ) ;  
		vP.release(); 
		vZ.release();
	    vV.release();
	    dvV.release();
		return vP;	  
	 }
	 else
	 {
		 nsteps=int((dividends_times[0])/dt);
		 double Tt;
		 Tt=(dividends_times[0])/nsteps;
		 
		 if (nsteps==1)
		    {
			MDL->Set_Rate(ov_rate);
			MDL->Make_Green_Put_last( lenx, leny, x, y, iJ, dividends_times[0]);//Build green matrix
			vP.resize(lenx);
			vP=MDL->project_back_Put_last(vZ);	
			//~ for(int it=0; it != lenx; ++it )
			   //~ vP(it+1) = max( payoff_put( exp( x[it] ), K ) , vP(it+1) ) ;  
			vP.release(); 
			vZ.release();
	        vV.release();
	        dvV.release();
			return vP;			
			}
		else
		    {
			 MDL->Make_Green_Put_yydtSI( leny, leny, y, y, iJ,Tt);
			 vZ = MDL->project_back_Put_yydtSI(vZ) ;	// Continuation value
             for(int it=0; it != leny; ++it )
	             vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;  
		    
		    for (int i=nsteps-1;i>1;i--){
			   if (flag==0)
			      {
					MDL->Make_Green_Put_yydtSI( leny, leny, y, y, iJ,Tt);//Build green matrix for the days in which there are no dividends				
				    flag=1;
				  }
			   vZ = MDL->project_back_Put_yydtSI(vZ) ;	// Continuation value
               for(int it=0; it != leny; ++it )
	               vZ(it+1) = max( vV(it + 1), vZ(it+1) ) ;   
			    }
				MDL->Set_Rate(ov_rate);
				MDL->Make_Green_Put_last( lenx, leny, x, y, iJ, Tt);//Build green matrix
				vP.resize(lenx);
				vP=MDL->project_back_Put_last(vZ);	
				//~ for(int it=0; it != lenx; ++it )
				    //~ vP(it+1) = max( payoff_put( exp( x[it] ), K ) , vP(it+1) ) ;  
				vP.release(); 
			    vZ.release();
	            vV.release();
	            dvV.release();
				return vP;			 
		        
		    }
	 }
	     
}	

ReturnMatrix AmericanCall_2D(const double &K, 
                          const double &rate, 
                          const double &ov_rate, 
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  const double &x0, const double &V0, Real *yy,
						  double *vv, double *grid, double *grid2, 
						    const int &ylen, const int &Wlen, 						  
						  const int &iJ1 ,const int &iJ2, 
						  const double dt,
						  Heston_FRP *HST,const double &shift, const int &Nrowstruncation)
{
	ColumnVector vP;//Output
	double Delta=yy[1]-yy[0];
	
    double (*pt2Payoff)(const double &, const double &) = &payoff_call;
	int numdiv=ndiv;
	int i1=numdiv-1;
	while (dividends_times[i1]>maturity && i1>-1)
	     i1-=1;
	numdiv=i1+1;
    Matrix vV, //intrinsic value
              vZ,//continuation value
             mZ;//value approximated
    vV.resize(ylen, Wlen);	vZ.resize(ylen, Wlen);
    
    /*PAYOFF*/		
for(int j=0; j!=Wlen; ++j)
		for(int i=0; i!=ylen; ++i) 
		vV.element(i,j)=pt2Payoff(exp(yy[i] + shift), K);
    vZ=vV;
	double dyield=HST->Get_dyield();
	int nsteps;
	/* CASO DIVIDENDYIELD!=0 */	
	if (dyield!=0)
	{
		maturity>dt? nsteps=int(maturity/dt) : nsteps=1;
		double t_dyield=maturity/nsteps;
		if (nsteps>1){
			HST->Make_Green( t_dyield);		
			for (int i=nsteps; i>1;i--)
			{
			 vZ = HST->project_back(vZ) ;	// Continuation value	 
		     //Uncomment if you want to print the green matrix
			//ofstream outfile;
			//outfile.open ("vZ.txt");
			//outfile<<vZ;
			//outfile.close();
  	 
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );			
			}
	     }
	     ////Uncomment if you want to print the green matrix
		//ofstream outfile;
		//outfile.open ("vZ.txt");
        //outfile<<vZ;
        //outfile.close();
		//Last step 1->0
		HST->Set_Rate(ov_rate);
		HST->Make_Green_0( x0,V0,t_dyield);//Build green matrix for the last step
		vP.resize(1);
		vP=HST->project_back_start(vZ);	
		vP.release();		vV.release();		vZ.release();		mZ.release();
		return vP;	
	
    }	
	//*********************//
    /* CASO NO DIVIDENDS */		
	if (numdiv==0)
	{
		//HST->Make_Green_0( x[1], V0, maturity);  //uses 2 dimensional transform
		HST->Make_Green_maturity(maturity,x0,V0);//uses 1 dimensional transform
		vP.resize(1);
		vP= HST->project_back_start(vZ,1);
		vV.release();		vZ.release();		mZ.release();
		vP.release(); return vP;
    }
    
 	//**************************//
    /* CASO DIVIDENDI DISCRETI */	
    //Questa parte è stata effettuata utilizzando l'approssimazione 
    //del payoff con espansione di taylor di secondo grado
     double dividend=dividends[0];
    
 		vector<double>  vY_m_d(ylen), vSmallDelta(ylen);		int iIndDrop[ylen];		vector<int>vSmallDeltaSign(ylen);			
		for(int i=0; i<ylen; ++i){									//	y grid, not translated, I need it to compute the dividend drop
			vY_m_d[i] = log( exp(yy[i]+ shift) - dividend);						//	ex dividend y grid
			iIndDrop[i] = int((yy[i] + shift- vY_m_d[i])/Delta);	
			if(i-iIndDrop[i]>0){											//	make sure that the ex dividend stock price is within the range
				int k = i-iIndDrop[i];										//	smallest integer such that yy[k] >= vY_m_d[i], includes dividend=0 case
				vSmallDelta[i] = vY_m_d[i] - (yy[k]+ shift);							//	Mind that vSmallDelta[i]<=0 !! Price to pay to include the dividend=0.0 case
				if(abs(vSmallDelta[i])/Delta< 0.5)
					vSmallDeltaSign[i] = 1;	
				else
					vSmallDeltaSign[i] = -1;}						//	delta to be used in the 1st order approx
			else { iIndDrop[i]=i;  vSmallDelta[i]=0;	vSmallDeltaSign[i]=0.0; };
		}


    ////T->T-1

    HST->Make_Green_maturity( maturity-dividends_times[numdiv-1]);//Build green matrix
    vZ = HST->project_back_maturity(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;	// Continuation value
	 for(int i=1; i<= ylen; ++i)	
		for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
			  vZ(i,w ) = max( vZ(i,w), vV(i,w) );	
	    
	////T-1->T-2,...2->1   

	int flag=0;
	for (int it=numdiv-1;it>0;it--)
	{
	    if (flag==0)
	    {
	        HST->Make_Green( dividends_times[it]-dividends_times[it-1]);//Build green matrix
	        flag=1;
	     }

	   vZ = HST->project_back_dividend(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;	// Continuation value

	 for(int i=1; i<= ylen; ++i)	
		for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
			  vZ(i,w ) = max( vZ(i,w), vV(i,w) );	 
	}
	//Last step 1->0
	HST->Make_Green_0( x0, V0, dividends_times[0]);  //It is done only for S0, not for 0.8*S0 or 1.2*S0
	vP.resize(1);
	vP= HST->project_back_start(vZ);	
	vV.release();		vZ.release();		mZ.release();
	vP.release(); return vP;
	
}	

ReturnMatrix AmericanPut_2D(const double &K, 
                          const double &rate, 
                          const double &ov_rate,
						  const double &maturity, 
						  const double *dividends,
						  const double *dividends_times,
						  const int &ndiv,
						  const double  &x0, const double &V0, Real *yy, 
						  double *vv, double *grid, double *grid2, 
						   const int &ylen, const int &Wlen,
						  const int &iJ ,const int &iJ2, 
						  const double dt,
						  Heston_FRP *HST,const double &shift, const int &Nrowstruncation)
{	
	ColumnVector vP;//Output
	double (*pt2Payoff)(const double &, const double &) = &payoff_put;
    double Delta=yy[1]-yy[0];	
    double Tt;
    int flag=0;
	int numdiv=ndiv;
	int i1=numdiv-1;
	while (dividends_times[i1]>maturity && i1>-1)
	     i1-=1;
	numdiv=i1+1;
    Matrix vV, //intrinsic value
              vZ,//continuation value
             mZ;//value approximated
    vV.resize(ylen, Wlen);	vZ.resize(ylen, Wlen);	
    
    /*PAYOFF*/
    for(int j=0; j!=Wlen; ++j)
		for(int i=0; i!=ylen; ++i) 
		vV.element(i,j)=payoff_put(exp(yy[i] + shift), K);
    
    vZ=vV;
	double dyield=HST->Get_dyield();
	int nsteps;  
	//*********************//
    /* CASO NO DIVIDENDS */	
	if (numdiv==0)
	{		
		maturity>dt? nsteps=int(maturity/dt) : nsteps=1;
		double t_dyield=maturity/nsteps;
		if (nsteps>1){
		HST->Make_Green( t_dyield);}
		for (int i=nsteps; i>1;i--)
		{
	     vZ = HST->project_back(vZ) ;	// Continuation value
         for(int i=1; i<= ylen; ++i)	
			for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
				  vZ(i,w ) = max( vZ(i,w), vV(i,w) );			
		}	
		//////////////////////
		//Last step 1->0
		HST->Set_Rate(ov_rate);
		HST->Make_Green_0( x0,V0,t_dyield);//Build green matrix for the last step
		vP.resize(1);
		vP=HST->project_back_start(vZ);	
		vP.release();		vV.release();		vZ.release();		mZ.release();
		return vP;					
    }	
  	//**************************//
    /* CASO DIVIDENDI DISCRETI */	
    //Questa parte è stata effettuata utilizzando l'approssimazione 
    //del payoff con espansione di taylor di secondo grado
    
        double dividend=dividends[0];
        Matrix dvV(ylen,Wlen);
		double yint[ylen];//grid of the log(S_t) after dividend payment	
		for(int i=0; i<ylen; ++i)	
			yint[i] = log( exp(yy[i]+shift) - dividends[0] );  

		for(int it=0; it != ylen; ++it ) 
		for(int j=0; j!=Wlen; ++j)
			dvV(it+1,j+1) = payoff_put( exp( yint[it] ), K ) ;	//payoff valutato at ex-dividends dates


		vector<double>  vY_m_d(ylen),vSmallDelta(ylen);
			int iIndDrop[ylen];		vector<int>vSmallDeltaSign(ylen);			
		
		for(int i=0; i<ylen; ++i){									//	y grid, not translated, I need it to compute the dividend drop
			vY_m_d[i] = log( exp(yy[i]+ shift) - dividend);						//	ex dividend y grid
			
			
			iIndDrop[i] = int((yy[i] + shift- vY_m_d[i])/Delta);				
			if(i-iIndDrop[i]>0){											//	make sure that the ex dividend stock price is within the range
				int k = i-iIndDrop[i];											//	smallest integer such that yy[k] >= vY_m_d[i], includes dividend=0 case
				vSmallDelta[i] = vY_m_d[i] - (yy[k]+ shift);							//	Mind that vSmallDelta[i]<=0 !! Price to pay to include the dividend=0.0 case
				if(abs(vSmallDelta[i])/Delta< 0.5)
					vSmallDeltaSign[i] = 1;	
				else
					vSmallDeltaSign[i] = -1;}						//	delta to be used in the 1st order approx
			else { iIndDrop[i]=i;  vSmallDelta[i]=0;	vSmallDeltaSign[i]=0.0; };
			
		}
	//////////////
	if ((maturity-dividends_times[numdiv-1])<dt)
	    {
         HST->Make_Green_maturity( maturity-dividends_times[numdiv-1]);//Build green matrix
	    vZ = HST->project_back_maturity(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;

         for(int i=1; i<= ylen; ++i)	
			for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
				  vZ(i,w ) = max( vZ(i,w), dvV(i,w) );	 
	     
	     }
	 else
	 {
		 nsteps=int((maturity-dividends_times[numdiv-1])/dt);
		 Tt=(maturity-dividends_times[numdiv-1])/nsteps;
		 
		 if (nsteps==1)
		    {
	        HST->Make_Green_maturity(maturity-dividends_times[numdiv-1]);//Build green matrix
	        vZ = HST->project_back_maturity(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), dvV(i,w) );	 
					  
		//Uncomment if you want to print the matrix
		//ofstream outfile;
		//outfile.open ("vZ.txt");
        //outfile<<vZ;
        //outfile.close();			    	
			}
		else
		    { 	
			 HST->Make_Green(Tt);//Build green matrix for the days in which there are no dividends
		    for (int i=nsteps;i>1;i--){
			   vZ = HST->project_back(vZ) ;	// Continuation value
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );	 	
			    }		    
	        vZ = HST->project_back_dividend(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;	// Continuation value
	        
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), dvV(i,w) );	 		   
		    }
	 }
	/////////////////////////////////////////////////////
	////Codice per il caso multidividend
	if (numdiv>1)
	{
      nsteps=int((dividends_times[numdiv-1]-dividends_times[numdiv-2])/dt); 
	  Tt=(dividends_times[numdiv-1]-dividends_times[numdiv-2])/nsteps;		
	  HST->Make_Green(Tt);
	for (int it=numdiv;it>1;it--)
		  { 
		   vZ=HST->project_back(vZ);
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );	 

		   for (int j=1; j<nsteps-1;j++)
			   {
				vZ = HST->project_back(vZ) ;	// Continuation value	
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );			   
			   }

		   vZ = HST->project_back_dividend(vZ,iIndDrop,vSmallDelta,vSmallDeltaSign,K,pt2Payoff) ;
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), dvV(i,w) );	
		  }
      }

	/////////////////////////////////////////////////
	if (dividends_times[0]<dt)
	    {
	    HST->Set_Rate(ov_rate);
	    HST->Make_Green_0( x0, V0, dividends_times[0]);
	    vP.resize(1);
		vP=HST->project_back_start(vZ) ;
		vP.release();		vV.release();		vZ.release();		mZ.release();
		return vP;	 
	 }
	 else
	 {
		 nsteps=int((dividends_times[0])/dt);
		 double Tt;
		 Tt=(dividends_times[0])/nsteps;
		 
		 if (nsteps==1)
		    {	
			HST->Set_Rate(ov_rate);
			HST->Make_Green_0( x0, V0, dividends_times[0]);
			vP.resize(1);
			vP=HST->project_back_start(vZ);
			vP.release();		vV.release();		vZ.release();		mZ.release();
			return vP;	 		
			}
		else
		    {
			 HST->Make_Green( Tt);
			vZ = HST->project_back(vZ) ;	
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );	
		    
		    for (int i=nsteps-1;i>1;i--){
			   if (flag==0)
			      {
					HST->Make_Green(Tt);//Build green matrix for the days in which there are no dividends				
				    flag=1;
				  }
			   vZ = HST->project_back(vZ) ;	// Continuation value
			 for(int i=1; i<= ylen; ++i)	
				for(int w=1; w<=Wlen; ++w)			//	Compare time and intrinsic value
					  vZ(i,w ) = max( vZ(i,w), vV(i,w) );
			    }
				HST->Set_Rate(ov_rate);
				HST->Make_Green_0( x0, V0, Tt);//Build green matrix
				vP.resize(1);
				vP=HST->project_back_start(vZ);	
				//~ for(int it=0; it != lenx; ++it )
				    //~ vP(it+1) = max( payoff_put( exp( x[it] ), K ) , vP(it+1) ) ;  
			vP.release();		vV.release();		vZ.release();		mZ.release();
			return vP;	 			 
		        
		    }
	 }
	     
}	

