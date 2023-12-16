#include "FRP_Heston_class.h"

Heston_FRP::Heston_FRP(  const double rate,
		          const double dyield,
				  const double longrun_var_Heston,
				  const double mean_rev_Heston,
				  const double volvol_Heston,
				  const double rho,
				  double lambda,
				  double muj,
				  double sigmaj,
				  const int Nrowstruncation,
				  const int ylen, const int Wlen,
				  const int glen, const int glen2,
				  const double *yy, const double *vv,				  
				  const double* grid,const double* grid2,
				  const double Delta, const double Delta2,const double shift,
				  const char Method) :
 rate(rate), dyield(dyield), longrun_var_Heston(longrun_var_Heston), mean_rev_Heston(mean_rev_Heston), 
 volvol_Heston(volvol_Heston), rho(rho),lambda(lambda), muj(muj), sigmaj(sigmaj), Nrowstruncation( Nrowstruncation ),ylen(ylen),Wlen(Wlen),glen(glen),glen2(glen2),
 yy(yy),vv(vv),grid(grid),grid2(grid2),
 Delta(Delta), Delta2(Delta2), shift(shift),
 method(Method)
 {
	 ik=ylen/2;
	 is_0=0;  is_mG=0;  is_maturity_mG=0;
	if (method=='S')
	    {
	     SampleHaar();
	     SampleHaar( W );//SampleHaar function computes mReTheta, mImTheta, mReEta, mImEta, 
	                     //which are the transforms of the basis functions.
          }
	 }
 
 
 Heston_FRP::~Heston_FRP()
 { 
	 if (is_0)
	    mG_0.release();
	 if (is_mG){
		 for(int w=0; w<Wlen; ++w){
					vmT[w].release();		
					vmT_nonZero[w].release();		   
				}
		}
	if (is_maturity_mG){
		for(int w=0; w<Wlen; ++w)
	         maturity_mG[w].release();
    }
    	mReEta.release(); mImEta.release();
		mRePsic.release(); mImPsic.release();
		mReTheta.release(); mImTheta.release(); 
 }

complex<double> Heston_FRP::Psic(const double &u1,const double& u2,const double& r,const double& dyield,const double &tau,
                      const double& ka,const double &theta, 
					 const double &alpha,const double& rho,const double& x0,const double &v0)
{
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%PARAMETERS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%Heston models according to the dynamics of S_t(stock) e v_t(variance):
//%dS_t=(rate-delta)S_tdt+sqrt(v_t)S_tDW_1t
//%dv_t=ka(theta-v_t)dt+alpha*sqrt(v_t)dW_1t
//%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%alpha: volvol
//%ka: mean rev speed
//% theta= long run variance
//%%%%%%%%%%%%%%%%%%%%%%%%%%%
complex<double> Psi0,EspPsiJ;
double beta(ka),gamma(alpha), alpha2(ka*theta);
double rate=r;
complex<double> i(0,1);
complex<double> A,B,p,q;
if (u1==0 && u2==0){
    Psi0=double(1);
    return Psi0;
}
A=gamma*gamma*(double(1)-rho*rho)*u1*u1+(gamma*gamma-double(2)*rho*gamma*beta)*i*u1+beta*beta;
B=-(rho*gamma*i*u1-beta-sqrt(A)+gamma*gamma*i*u2)/(rho*gamma*i*u1-beta+sqrt(A)+gamma*gamma*i*u2);
q=(beta-rho*gamma*i*u1-sqrt(A)*(B*exp(tau*sqrt(A))-double(1))/(B*exp(tau*sqrt(A))+double(1)))/(gamma*gamma);
p=tau*(rate-dyield -alpha2*rho/gamma)*i*u1+alpha2/(gamma*gamma)*(tau*beta+tau*sqrt(A)+double(2)*log((B+double(1))/(B*exp(tau*sqrt(A))+double(1))));

EspPsiJ=lambda*tau*(exp(i*muj*u1-double(0.5)*sigmaj*sigmaj*u1*u1)-double(1))-i*lambda*tau*u1*(-double(1)+exp(muj+double(0.5)*sigmaj*sigmaj));

Psi0=exp(p+q*v0+i*u1*x0+EspPsiJ);

return Psi0;
	
}

ReturnMatrix  Heston_FRP::TransitionM( const double &t, const double &x0, const double &v0, int &Truncation){
	//OUTPUT
	//Matrix mT(ylen,Wlen). 
	//mT(i,j)=probability of reaching a state (S(i),v(j)) in a time t starting from the state (x0,v0).
	//INPUT
	//Truncation is an integer number (from 0 to glen/2) stating how many diagonals of the green matrix we consider different than zero.
	//If Truncation=0, the green matrix will be a zero matrix (we don't consider any diagonal)
	//If Truncation=glen/2 we consider the full green matrix (no Truncation)
	//the suitable value for Truncation should be set accordingly to the variable t (delta time of the green matrix),
	//as the higher is t the more diagonals we should consider.
	double b (mean_rev_Heston), theta(longrun_var_Heston), alpha(volvol_Heston);
	double Disc( exp(-rate*t) );
	Matrix mT(ylen, Wlen);
	double Delta3= (grid[1]-grid[0]);	
	double Delta4= (grid2[1]-grid2[0]);	
	complex<double> Psi0;		
		if (method=='F')//FFT method for inverting the Fourier transform of the green function
		
		  {
	        complex<double> I(0.0, 1.0), ComplexDens; 
			double DeltaF = 2*M_PI/glen/(grid[1]-grid[0]);	double DeltaF2 = 2*M_PI/glen2/(grid2[1]-grid2[0]);
			complex<double>*  in=new complex<double> [glen*glen2];
			complex<double>* out=new complex<double> [glen*glen2];
			fftw_plan plan;
									
			plan = fftw_plan_dft_2d(glen, glen2, reinterpret_cast<fftw_complex*>(in), 
									reinterpret_cast<fftw_complex*>(out), 
									FFTW_FORWARD, FFTW_ESTIMATE);				//Here we set which kind of transformation we want to perform
			
			
			for(int j=0; j!=glen2; ++j){
				for(int i=0; i!=glen; ++i){
					Psi0=Psic( grid[i], grid2[j], rate,dyield, t, b, theta, alpha, rho, x0, v0);
					in[j + glen2*i] =	Psi0*exp( -I* (double)i * Delta3 * yy[0] )*exp( -I*(double)j * Delta4 * vv[0] );	
				}
			}		

				
			fftw_execute(plan); //Execution of FFT
			mT=0;
			for(int j=0; j!=glen2; ++j)
				for(int i=glen/2-Truncation; i!=glen/2+Truncation; ++i){
					ComplexDens = exp(-I*grid[0]*DeltaF*(double)i)*exp(-I*grid2[0]*DeltaF2*(double)j)*out[j + glen2*i];
					mT.element(i,j)=1/(4*PI*PI)*Delta3*Delta4*real(ComplexDens)*Disc*Delta*Delta2;
				}
			
					////Uncomment if you want to print the green matrix
				//ofstream outfile;
				//outfile.open ("Green_Matrix.txt");
				//outfile<<mT;
				//outfile.close();
				
				delete [] in;	delete [] out;
				
			}
		else// sampling method
		    {
                 mRePsic.resize( glen, glen2 );
	             mImPsic.resize( glen, glen2 );
                //Transform of the green matrix
				for(int k=0; k!= glen; ++k){	
					for (int l=0; l!= glen2; ++l) {			
						Psi0 = Psic( grid[k], grid2[l], rate,dyield, t, b, theta, alpha, rho, x0, v0) ;
						mRePsic(k+1, l+1) = real(Psi0);	mImPsic(k+1,l+1) = imag(Psi0);
					}
				}
	           //Parseval multiplications in the Fourier space with the basis functions 
	           //which give us the green function
	           mT= Disc*( mReEta * (mRePsic*mReTheta - mImPsic*mImTheta) - mImEta * (mRePsic*mImTheta + mImPsic*mReTheta) ) *
	                  Delta3 * Delta4 * sqrt(Delta) * sqrt(Delta2) / 4 / PI / PI;			       	            		        
			}
		        	
	mT.Release();	return mT;
	
}

ReturnMatrix  Heston_FRP::TransitionM_1D( const double &t, const double &x0, const double &v0){
	//OUTPUT
	//Matrix mT(ylen,1). 
	//mT(i)=probability of reaching a state (S(i)) in a time t starting from the state (x0,v0).
    /////////////////////////////////////////////////////////////////////////////////////////////
	double b (mean_rev_Heston), theta(longrun_var_Heston), alpha(volvol_Heston);
	double RePsi, ImPsi;
	double Disc( exp(-rate*t) );
	Matrix mT(ylen,1);
	double Delta3= (grid[1]-grid[0]);	
	complex<double> Psi0;		
		if (method=='F')//FFT method for inverting the Fourier transform of the green function
		
		  {
	        complex<double> I(0.0, 1.0), ComplexDens; 
			double DeltaF = 2*M_PI/glen/(grid[1]-grid[0]);	
			complex<double>*  in=new complex<double> [glen];
			complex<double>* out=new complex<double> [glen];
			fftw_plan plan;
									
			plan = fftw_plan_dft_1d(glen,  reinterpret_cast<fftw_complex*>(in), 
									reinterpret_cast<fftw_complex*>(out), 
									FFTW_FORWARD, FFTW_ESTIMATE);				//Here we set which kind of transformation we want to perform
			
			
			
				for(int i=0; i!=glen; ++i){
					PsiFn (	 t, 0.0, grid[i],lambda,muj,sigmaj, longrun_var_Heston*mean_rev_Heston, mean_rev_Heston, volvol_Heston, rho, x0+rate*t, v0, &RePsi, &ImPsi);
					Psi0=exp(RePsi+I*ImPsi);
					in[i] =	Psi0*exp( -I* (double)i * Delta3 * yy[0] );	
				}
					

				
			fftw_execute(plan); //Execution of FFT
			mT=0;
			
				for(int i=0; i!=glen; ++i){
					ComplexDens = exp(-I*grid[0]*DeltaF*(double)i)*out[i];
					mT.element(i,0)=1/(2*PI)*Delta3*real(ComplexDens)*Disc*Delta;
				}
			
					////Uncomment if you want to print the green matrix
				//ofstream outfile;
				//outfile.open ("Green_Matrix.txt");
				//outfile<<mT;
				//outfile.close();
				
				delete [] in;	delete [] out;
				
			}
		else// sampling method
		    {
				//////////////////////////////////////////
	
	                     
		Matrix ReMatrix( 1, glen );
		Matrix ImMatrix( 1, glen );        

		for (int l=0; l!= glen; ++l) {
			PsiFn (	 t, 0.0, grid[l],lambda,muj,sigmaj, longrun_var_Heston*mean_rev_Heston, mean_rev_Heston, volvol_Heston, rho, x0+rate*t, v0, &RePsi, &ImPsi);
			ReMatrix(+1, l+1) = exp(RePsi)*cos(ImPsi) ;	ImMatrix(1,l+1) = exp(RePsi)*sin(ImPsi) ;
		
	        }
	     mT = Disc*( ReMatrix * mReEta.t() - ImMatrix * mImEta.t() ).t() * Delta3 *  sqrt(Delta) / 2 / PI;          
	               
		       	            		        
			}
		        	
	mT.Release();	return mT;
	
}




void Heston_FRP::Make_Green(const double &t)
{
	if (!is_mG){
		   vmT.resize(Wlen);
		   vmT_nonZero.resize(Wlen);
			for(int w=0; w<Wlen; ++w){
				vmT[w] = TransitionM( t, yy[ik], vv[w],Nrowstruncation);				   
			}	
	 is_mG=1;
	}	
}

void Heston_FRP::Make_Green_maturity(const double &t)
{
	if (!is_maturity_mG){
		maturity_mG.resize(Wlen);
			for(int w=0; w<Wlen; ++w){
				maturity_mG[w] = TransitionM( t, yy[ik], vv[w],Nrowstruncation);
			}		
	 is_maturity_mG=1;
	}	
}	

void Heston_FRP::Make_Green_maturity(const double &t,const double &x0, const double &v0)
{
	if (!is_0){
		mG_0.resize(ylen,1);	
		mG_0=TransitionM_1D( t, x0, v0); 
		is_0=1;
	}	
}	

void Heston_FRP::Make_Green_0( const double &x0, const double &v0, const double& t)
{
	int Truncation=glen/2;
	if (!is_0){
		mG_0.resize(ylen,Wlen);
		  mG_0=TransitionM( t, x0, v0, Truncation); 
		  is_0=1;
	   }	
}

double Heston_FRP::project_back_start( Matrix &mV )
{
	return dotproduct(mG_0,mV);
		   
	//int rowstart=1; 
	//int rowend=ylen;
	//return dotproduct_sparse(mV, mG_0, vmT_nonZero[0], rowstart,rowend);
}

double Heston_FRP::project_back_start( Matrix &mV ,const int flag)
{
	return dotproduct(mG_0,mV.columns(1,1));
		   

}

double Heston_FRP::dotproduct_sparse(Matrix &A, Matrix &Sparse, Matrix &Coordinates, int & vmTrowsStart,int &vmTrowsEnd)
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%
	//Dotproduct A dot Sparse where Sparse is a sparse matrix
	//Dimension of A and Sparse is (ylen x Wlen)
	//%%%%%%%%%%%%%%%%%%%%%%%%%%
	double result=0;
	int index=0;

	if (Coordinates.element(index,2)==0)
		return result;
	while (index<ylen*Wlen){
		if (Coordinates.element(index,0)<vmTrowsStart-1){
		   if (Coordinates.element(index,2)==0)
		       return result;
	       index+=1;}
	    else
	     break;
   }
	while (index<ylen*Wlen)
	{
		if (Coordinates.element(index,0)<vmTrowsEnd)
		{
		   if (Coordinates.element(index,2)==0)
		       return result;
		   result+=Coordinates.element(index,2)*A.element(Coordinates.element(index,0)-vmTrowsStart+1,Coordinates.element(index,1));
		   index+=1;
		  }
		  else
		    return result;
	}
	return result;
	
}

ReturnMatrix Heston_FRP::project_back( Matrix &mV ) //senza approssimazione
{
	Matrix mU,mZ;  	mU.resize(ylen,Wlen);
	int id,iu,vmTrowsStart,vmTrowsEnd;
	for(int l=0; l<ik; ++l){// l:	loop on the conditioning log(S)
		id = 1; iu = ylen - ik + l;											// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )									                           // submatrix of mV actually used (instead of padding the shifted transition matrix with zero					
	   vmTrowsStart=max(ylen - iu + 1,ylen/2-Nrowstruncation);			
	   vmTrowsEnd=min(ylen,ylen/2+Nrowstruncation);
	   mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );	
		for(int m = 0; m< Wlen; ++m)										// m:	loop on the conditioning variance, interpolation to be done only in the log(S) dimension
			mU.element(l,m) = dotproduct( mZ, (vmT[m]).rows(vmTrowsStart, vmTrowsEnd) );	
		    //mU.element(l,m) =dotproduct_sparse(mZ, vmT[m], vmT_nonZero[m],vmTrowsStart,vmTrowsEnd);
		
	}//end for l
	
		
	for(int l=ik; l<ylen; ++l){										// loop on the conditioning log(S)
		id = l + 1 - ik; iu =ylen;	
		vmTrowsStart=max(1,ylen/2-Nrowstruncation);					
	   vmTrowsEnd=min(ylen - id +1,ylen/2+Nrowstruncation);										// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )
		mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	
		mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );					// submatrix of mV actually used (instead of padding the shifted transition matrix with zeros
		for(int m = 0; m< Wlen; ++m)										// loop on the conditioning variance
			mU.element(l,m) =dotproduct( mZ, (vmT[m]).rows(vmTrowsStart, vmTrowsEnd) );
			//mU.element(l,m) =dotproduct_sparse(mZ, vmT[m], vmT_nonZero[m],vmTrowsStart,vmTrowsEnd);
	}//end for l
	mZ.release();
	mU.release();
	return mU;
}

ReturnMatrix Heston_FRP::project_back_dividend( Matrix &mV,
                                                                      int *iIndDrop,vector<double> &vSmallDelta,
                                                                      vector<int>&vSmallDeltaSign,const double &K,
                                                                      double (*pt2Payoff)(const double &s, const double &k) )//con approssimazione
{
	Matrix mU,mZ;  	mU.resize(ylen,Wlen);
	int id,iu,vmTrowsStart,vmTrowsEnd;
	double gamma;
	for(int l=0; l<ik; ++l){// l:	loop on the conditioning log(S)
		id = 1; iu = ylen - ik + l;											// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )									                           // submatrix of mV actually used (instead of padding the shifted transition matrix with zero					
	   vmTrowsStart=max(ylen - iu + 1,ylen/2-Nrowstruncation);					
	   vmTrowsEnd=min(ylen,ylen/2+Nrowstruncation);
	   mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );	
	   //INTERPOLATION
		for(int j=0; j<Wlen; ++j){	
			for(int i=vmTrowsStart-ik+l-1; i<vmTrowsEnd-ik+l; ++i){								// i si an index of vY, vY_m_d, so i = (row of mZ) -1
				int  q; q=i-iIndDrop[l];
				if(q>2 && q<(ylen-1))	{ 
					if(vSmallDeltaSign[l]==1){
						gamma = ( mV.element(q+1,j) -2*mV.element(q,j) + mV.element(q-1,j) ) / (Delta*Delta);					
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q,j) + ( mV.element(q+1,j) - mV.element(q-1,j) )/2/Delta * vSmallDelta[l] 
						+ 0.5*gamma*vSmallDelta[l]*vSmallDelta[l];																
					}
					else{
						gamma = ( mV.element(q,j) -2*mV.element(q-1,j) + mV.element(q-2,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q-1,j) + ( mV.element(q,j) - mV.element(q-2,j) )/2/Delta * (Delta + vSmallDelta[l]) 
						+ 0.5*gamma*(Delta + vSmallDelta[l])*(Delta + vSmallDelta[l]);
					}

				}//end if(q>2 && q<(ylen-1))
				else if (q <= 2)	mZ.element(i-id+1,j) = 0.0;
				else	mZ.element(i-(vmTrowsStart-ik+l)+1,j) =pt2Payoff( exp(yy[ylen-1] + shift), K )  ;

			}//end for i
		
		}// end for j                  
	   /////////////////////	   
		for(int m = 0; m< Wlen; ++m)										// m:	loop on the conditioning variance, interpolation to be done only in the log(S) dimension
			mU.element(l,m) = dotproduct( mZ, (vmT[m]).rows(vmTrowsStart, vmTrowsEnd) );
	}//end for l
	for(int l=ik; l<ylen; ++l){										// loop on the conditioning log(S)
		id = l + 1 - ik; iu =ylen;	
		vmTrowsStart=max(1,ylen/2-Nrowstruncation);					
	   vmTrowsEnd=min(ylen - id +1,ylen/2+Nrowstruncation);										// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )
		mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	
		mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );					// submatrix of mV actually used (instead of padding the shifted transition matrix with zeros
	  //INTERPOLATION
		for(int j=0; j<Wlen; ++j){										// Computing elements of mZ by interpolation
			for(int i=vmTrowsStart-ik+l-1; i<vmTrowsEnd-ik+l ; ++i){									// i si an index of vY, vY_m_d, so it is row of mZ -1
				int  q; q=i-iIndDrop[l];
				if(q>2 && q<(ylen-1))	{ 
					if(vSmallDeltaSign[l]==1){
						gamma = ( mV.element(q+1,j) -2*mV.element(q,j) + mV.element(q-1,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q,j) + ( mV.element(q+1,j) - mV.element(q-1,j) )/2/Delta * vSmallDelta[l] + 0.5*gamma*vSmallDelta[l]*vSmallDelta[l];
					}
					else{
						gamma = ( mV.element(q,j) -2*mV.element(q-1,j) + mV.element(q-2,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q-1,j) + ( mV.element(q,j) - mV.element(q-2,j) )/2/Delta * (Delta + vSmallDelta[l]) + 0.5*gamma*(Delta + vSmallDelta[l])*(Delta + vSmallDelta[l]);
					}
				}//end if(q>2 && q<(ylen-1))
				else if (q <= 2)	mZ.element(i-id+1,j) = 0.0;
				else	mZ.element(i-(vmTrowsStart-ik+l)+1,j) = pt2Payoff( exp(yy[ylen-1] + shift), K ) ;
			}//end for i
		}//end for j                 
	  //////////////////////////		
		for(int m = 0; m< Wlen; ++m)										// loop on the conditioning variance
			mU.element(l,m) =dotproduct( mZ, (vmT[m]).rows(vmTrowsStart, vmTrowsEnd) );
	}//end for l	
	mZ.release();
	mU.release();
	return mU;
}


ReturnMatrix Heston_FRP::project_back_maturity( Matrix &mV,
                                                                      int *iIndDrop,vector<double> &vSmallDelta,
                                                                      vector<int>&vSmallDeltaSign,const double &K,
                                                                      double (*pt2Payoff)(const double &s, const double &k) )//con approssimazione
{
	Matrix mU,mZ;  	mU.resize(ylen,Wlen);
	int id,iu,vmTrowsStart,vmTrowsEnd;
	double gamma;
	for(int l=0; l<ik; ++l){// l:	loop on the conditioning log(S)
		id = 1; iu = ylen - ik + l;											// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )									                           // submatrix of mV actually used (instead of padding the shifted transition matrix with zero					
	   vmTrowsStart=max(ylen - iu + 1,ylen/2-Nrowstruncation);					
	   vmTrowsEnd=min(ylen,ylen/2+Nrowstruncation);
	   mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );	
	   //INTERPOLATION
		for(int j=0; j<Wlen; ++j){	
			for(int i=vmTrowsStart-ik+l-1; i<vmTrowsEnd-ik+l; ++i){								// i si an index of vY, vY_m_d, so i = (row of mZ) -1
				int  q; q=i-iIndDrop[l];
				if(q>2 && q<(ylen-1))	{ 
					if(vSmallDeltaSign[l]==1){
						gamma = ( mV.element(q+1,j) -2*mV.element(q,j) + mV.element(q-1,j) ) / (Delta*Delta);					
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q,j) + ( mV.element(q+1,j) - mV.element(q-1,j) )/2/Delta * vSmallDelta[l] 
						+ 0.5*gamma*vSmallDelta[l]*vSmallDelta[l];																
					}
					else{
						gamma = ( mV.element(q,j) -2*mV.element(q-1,j) + mV.element(q-2,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q-1,j) + ( mV.element(q,j) - mV.element(q-2,j) )/2/Delta * (Delta + vSmallDelta[l]) 
						+ 0.5*gamma*(Delta + vSmallDelta[l])*(Delta + vSmallDelta[l]);
					}

				}//end if(q>2 && q<(ylen-1))
				else if (q <= 2)	mZ.element(i-id+1,j) = 0.0;
				else	mZ.element(i-(vmTrowsStart-ik+l)+1,j) =pt2Payoff( exp(yy[ylen-1] + shift), K)  ;

			}//end for i
		
		}// end for j                  
	   /////////////////////	   
		for(int m = 0; m< Wlen; ++m)										// m:	loop on the conditioning variance, interpolation to be done only in the log(S) dimension
			mU.element(l,m) = dotproduct( mZ, (maturity_mG[m]).rows(vmTrowsStart, vmTrowsEnd) );
	}//end for l
	for(int l=ik; l<ylen; ++l){										// loop on the conditioning log(S)
		id = l + 1 - ik; iu =ylen;	
		vmTrowsStart=max(1,ylen/2-Nrowstruncation);					
	   vmTrowsEnd=min(ylen - id +1,ylen/2+Nrowstruncation);										// rows of the "mother" transition matrix to use (when l=ik, vmT[m] )
		mZ.resize(vmTrowsEnd-vmTrowsStart+1, Wlen);	
		mZ = mV.rows( vmTrowsStart-ik+l, vmTrowsEnd-ik+l );					// submatrix of mV actually used (instead of padding the shifted transition matrix with zeros
	  //INTERPOLATION
		for(int j=0; j<Wlen; ++j){										// Computing elements of mZ by interpolation
			for(int i=vmTrowsStart-ik+l-1; i<vmTrowsEnd-ik+l ; ++i){									// i si an index of vY, vY_m_d, so it is row of mZ -1
				int  q; q=i-iIndDrop[l];
				if(q>2 && q<(ylen-1))	{ 
					if(vSmallDeltaSign[l]==1){
						gamma = ( mV.element(q+1,j) -2*mV.element(q,j) + mV.element(q-1,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q,j) + ( mV.element(q+1,j) - mV.element(q-1,j) )/2/Delta * vSmallDelta[l] + 0.5*gamma*vSmallDelta[l]*vSmallDelta[l];
					}
					else{
						gamma = ( mV.element(q,j) -2*mV.element(q-1,j) + mV.element(q-2,j) ) / (Delta*Delta);
						mZ.element(i-(vmTrowsStart-ik+l)+1,j) = mV.element(q-1,j) + ( mV.element(q,j) - mV.element(q-2,j) )/2/Delta * (Delta + vSmallDelta[l]) + 0.5*gamma*(Delta + vSmallDelta[l])*(Delta + vSmallDelta[l]);
					}
				}//end if(q>2 && q<(ylen-1))
				else if (q <= 2)	mZ.element(i-id+1,j) = 0.0;
				else	mZ.element(i-(vmTrowsStart-ik+l)+1,j) =pt2Payoff( exp(yy[ylen-1] + shift), K) ;
			}//end for i
		}//end for j                 
	  //////////////////////////		
		for(int m = 0; m< Wlen; ++m)										// loop on the conditioning variance
			mU.element(l,m) =dotproduct( mZ, (maturity_mG[m]).rows(vmTrowsStart, vmTrowsEnd) );
	}//end for l	
	mZ.release();
	mU.release();
	return mU;
}

double Heston_FRP::Get_dyield(){return dyield;}

void Heston_FRP::Set_Rate(double Rate){rate=Rate;};

//////////////////////////////////////////////////////
//Functions for the sampling method
//////////////////////////////////////////////////////
void Heston_FRP::SampleHaar( const Dimension WhichDimension ){
	
	ColumnVector yy_t;
	yy_t.resize(ylen);	yy_t << yy;
	Matrix	mLim = (( yy_t - 0.5*Delta ) | ( yy_t + 0.5*Delta )).t();
		
	ColumnVector vv_t;
	vv_t.resize(Wlen); vv_t<<vv;
	Matrix mLim2 = (( vv_t - 0.5*Delta2 ) | ( vv_t + 0.5*Delta2 )).t();
	
	switch (WhichDimension){
		
		case W:
			
			mReTheta.resize(glen2, Wlen);	mImTheta.resize(glen2, Wlen);
			
			for( int l = 1; l<= Wlen; ++l ){
				mReTheta.Column(l) = ( ReEhat_du_matrix( glen2, grid2, mLim2(1,l), mLim2(2,l) ) ).t();	//Direct wav matrix
				mImTheta.Column(l) = ( ImEhat_du_matrix( glen2, grid2, mLim2(1,l), mLim2(2,l) ) ).t();	//Direct wav matrix
			}
			
			break;
			
		default:
			
			mReEta.resize(ylen, glen);	mImEta.resize(ylen, glen);
			
			for( int l = 1; l<= ylen; ++l ){
				mReEta.Row(l) = ReEhat_du_matrix( glen, grid, mLim(1,l), mLim(2,l) );	//Direct wav matrix
				mImEta.Row(l) = ImEhat_du_matrix( glen, grid, mLim(1,l), mLim(2,l) );	//Direct wav matrix
			}
			
			break;
			
	}
	mLim.release();	yy_t.release();	
	mLim2.release(); vv_t.release();    
	
}

void Heston_FRP::PsiFn (
			double T,
			double U,
			double V,
            double lambda,
	        double muJ,	
	        double sigmaJ,
			double a,
			double b,
			double alpha,
			double rho,
			double x1,
			double x2,
			double *RePsi,
			double *ImPsi
			)

{
	double	Phi1,Phi2,Phi3,Phi4,Psi1,Psi2,Psi3,Psi4,NumRe,NumIm,DenRe,DenIm, Chi, Den, 
	PhiG, PsiG, PhiG2, PsiG2, Phi3G, Psi3G, SinT, Sin2T, CosT, Cos2T, U2=U*U, U3=U*U*U, 
	V2=V*V, V3=V*V*V, r=0.0, Phi2G,Psi2G,Esp,Esp2,ReJumps,ImJumps;
	
	//  Defines the Re and Im part of the square of the function gamma 
	
	PhiG2=		b*b + U*alpha*alpha - U2*alpha*alpha + 
	V2*alpha*alpha - 2*b*U*alpha*rho + 
	U2*alpha*alpha*rho*rho - 
	V2*alpha*alpha*rho*rho;
	
	PsiG2=	V*alpha*alpha - 2*U*V*alpha*alpha - 2*b*V*alpha*rho + 
	2*U*V*alpha*alpha*rho*rho;
	
	//  Defines the Re and Im part of the function gamma 
	
	PhiG=pow(PhiG2*PhiG2+PsiG2*PsiG2,0.25)*cos(atan(PsiG2/PhiG2)/2.);
	PsiG=pow(PhiG2*PhiG2+PsiG2*PsiG2,0.25)*sin(atan(PsiG2/PhiG2)/2.);
	
	Phi2G=PhiG*PhiG;
	Psi2G=PsiG*PsiG;
	Phi3G=PhiG*PhiG*PhiG;
	Psi3G=PsiG*PsiG*PsiG;
	
	SinT=sin(T*PsiG);
	Sin2T=SinT*SinT;
	CosT=cos(T*PsiG);
	Cos2T=CosT*CosT;
	Esp=exp(-T*PhiG);
	Esp2=exp(-2.*T*PhiG);
	
	//  Defines the L1 function 
	Phi1=x1*U;
	Psi1=x1*V;
	
	//  Defines the L2 function 
	NumRe=
	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
	U3*x2*alpha*rho + V2*x2*alpha*rho - 
	U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
	V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
	
	(2*b*U*x2*CosT)*Esp - 
	(2*b*U2*x2*CosT)*Esp +
	(2*b*V2*x2*CosT)*Esp - 
	(2*U2*x2*alpha*rho*CosT)*Esp + 
	(2*U3*x2*alpha*rho*CosT)*Esp - 
	(2*V2*x2*alpha*rho*CosT)*Esp + 
	(2*U*V2*x2*alpha*rho*CosT)*Esp - 
	
	(b*U*x2*Cos2T)*Esp2 + 
	(b*U2*x2*Cos2T)*Esp2 - 
	(b*V2*x2*Cos2T)*Esp2 + 
	(U2*x2*alpha*rho*Cos2T)*Esp2 - 
	(U3*x2*alpha*rho*Cos2T)*Esp2 + 
	(V2*x2*alpha*rho*Cos2T)*Esp2 - 
	(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
	(U*x2*PhiG*Cos2T)*Esp2 - 
	(U2*x2*PhiG*Cos2T)*Esp2 + 
	(V2*x2*PhiG*Cos2T)*Esp2 + 
	(V*x2*PsiG*Cos2T)*Esp2 - 
	(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
	
	(2*V*x2*PhiG*SinT)*Esp - 
	(4*U*V*x2*PhiG*SinT)*Esp - 
	(2*U*x2*PsiG*SinT)*Esp + 
	(2*U2*x2*PsiG*SinT)*Esp - 
	(2*V2*x2*PsiG*SinT)*Esp - 
	
	(b*U*x2*Sin2T)*Esp2 + 
	(b*U2*x2*Sin2T)*Esp2 - 
	(b*V2*x2*Sin2T)*Esp2 + 
	(U2*x2*alpha*rho*Sin2T)*Esp2 - 
	(U3*x2*alpha*rho*Sin2T)*Esp2 + 
	(V2*x2*alpha*rho*Sin2T)*Esp2 - 
	(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
	(U*x2*PhiG*Sin2T)*Esp2 - 
	(U2*x2*PhiG*Sin2T)*Esp2 + 
	(V2*x2*PhiG*Sin2T)*Esp2 + 
	(V*x2*PsiG*Sin2T)*Esp2 - 
	(2*U*V*x2*PsiG*Sin2T)*Esp2;
	
	
	NumIm=
	-(b*V*x2) + 2*b*U*V*x2 - 
	U2*V*x2*alpha*rho -
	V3*x2*alpha*rho - 
	V*x2*PhiG + 2*U*V*x2*PhiG + 
	U*x2*PsiG - U2*x2*PsiG + 
	V2*x2*PsiG + 
	
	(2*b*V*x2*CosT)*Esp - 
	(4*b*U*V*x2*CosT)*Esp + 
	(2*U2*V*x2*alpha*rho*CosT)*Esp +
	(2*V3*x2*alpha*rho*CosT)*Esp - 
	
	(b*V*x2*Cos2T)*Esp2 +  
	(2*b*U*V*x2*Cos2T)*Esp2 - 
	(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
	(V3*x2*alpha*rho*Cos2T)*Esp2 + 
	(V*x2*PhiG*Cos2T)*Esp2 - 
	(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
	(U*x2*PsiG*Cos2T)*Esp2 + 
	(U2*x2*PsiG*Cos2T)*Esp2 - 
	(V2*x2*PsiG*Cos2T)*Esp2 - 
	
	(2*U*x2*PhiG*SinT)*Esp + 
	(2*U2*x2*PhiG*SinT)*Esp - 
	(2*V2*x2*PhiG*SinT)*Esp - 
	(2*V*x2*PsiG*SinT)*Esp + 
	(4*U*V*x2*PsiG*SinT)*Esp - 
	(b*V*x2*Sin2T)*Esp2 + 
	(2*b*U*V*x2*Sin2T)*Esp2 - 
	(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
	(V3*x2*alpha*rho*Sin2T)*Esp2 + 
	(V*x2*PhiG*Sin2T)*Esp2 - 
	(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
	(U*x2*PsiG*Sin2T)*Esp2 + 
	(U2*x2*PsiG*Sin2T)*Esp2 - 
	(V2*x2*PsiG*Sin2T)*Esp2;
	
	
	DenRe=
	b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
	(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
	(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;
	
	DenIm=
	-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
	(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
	(PhiG*SinT)*Esp;
	
	Phi2=NumRe/(DenRe*DenRe+DenIm*DenIm);
	Psi2=NumIm/(DenRe*DenRe+DenIm*DenIm);
	
	//  Defines the L3 function 
	Phi3=T/alpha/alpha*(r*alpha*alpha*(U-1)+a*b-a*U*alpha*rho-a*PhiG);
	Psi3=T/alpha/alpha*(r*V*alpha*alpha-a*V*alpha*rho-a*PsiG);
	
	
	//  Defines the L4 function 
	
	Chi=PhiG*PhiG+PsiG*PsiG;
	Den=2*Chi;
	
	
	NumRe=
	b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
	rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
	(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
													  PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
																							PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp;
	
	NumIm=
	-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
													   alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
	(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
	(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
													  PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;
	
	Phi4=-a/alpha/alpha*log(NumRe*NumRe/Den/Den+NumIm*NumIm/Den/Den);
	Psi4=-2.*a/alpha/alpha*atan(NumIm/NumRe);
	
	//Re and Im part of the jump component
	ReJumps=-lambda*T+exp(-0.5*sigmaJ*sigmaJ*V*V)*lambda*T*cos(muJ*V);
    ImJumps= lambda*T*V-exp(muJ+0.5*sigmaJ*sigmaJ)*lambda*T*V+exp(-0.5*sigmaJ*sigmaJ*V*V)*lambda*T*sin(muJ*V);
	// computes Re and Im part of the Psi function 
	
	*RePsi=Phi1+Phi2+Phi3+Phi4+ReJumps;
	*ImPsi=Psi1+Psi2+Psi3+Psi4+ImJumps;	
	
}	
