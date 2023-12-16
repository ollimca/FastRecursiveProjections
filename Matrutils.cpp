/*
 *  MatrUtils.h
 *  WavHeston
 *
 *  Created by Antonio Cosma on 25/February/2011.
 *  Copyright 2011 Université du Luxembourg. All rights reserved.
 *
 */

#include "Matrutils.h"


int find_index(const Matrix &mX, const Real X){
	double val = X;
	int col = mX.Ncols();
	int row = mX.Nrows();
	int ires(-1);
	int nmax( max(col, row) );
	if(col>row && row>1) cout << "find_index: argument is not a vector" << endl;
	if(col<row && col>1) cout << "find_index: argument is not a vector" << endl;
	if(nmax==col){
		int i=1;
		while(ires==-1 && (i<=nmax) ){
			if(mX(1,i)>= (val-1E-7)) ires = i; 
			++i;
		}
	}
	if(nmax==row){
		int i=1;
		while(ires==-1 && (i<=nmax) ){
			if(mX(i,1)>=(val-1E-7)) ires = i;
			++i;
		}
	}
	return ires;
}

int find_closest_index(const Matrix &mX, const Real X){
	double val = X;
	int col = mX.Ncols();
	int row = mX.Nrows();
	int ires(-1);
	int nmax( max(col, row) );
	if(col>row && row>1) cout << "find_index: argument is not a vector" << endl;
	if(col<row && col>1) cout << "find_index: argument is not a vector" << endl;
	if(nmax==col){
		int i=1;
		while(ires==-1 && (i<=nmax) ){
			if(mX(1,i)>= (val-1E-7)) ires = i;
			++i;
		}
		if(ires > 1) if ( fabs(mX(1, ires) - val ) >  fabs(mX(1, ires - 1) - val ) ) --ires;
	}
	if(nmax==row){
		int i=1;
		while(ires==-1 && (i<=nmax) ){
			if(mX(i,1)>=(val-1E-7)) ires = i;
			++i;
		}
		if(ires > 1) if ( fabs(mX(ires, 1) - val ) >  fabs(mX(ires - 1, 1) - val ) ) --ires;
	}
	return ires;
}


void print_table( vector<string> vLables, Matrix mM, int iWidth, int iPrecision, string outDevice, const char sep, string fileMode ){
	
	switch (sep) {
			
		case 't' :
			if( outDevice == "screen"){
				for(vector<string>::size_type sv = 0; sv != vLables.size(); ++ sv)
					cout << left <<  setw(iWidth+1) << vLables[sv]  << "\t" ;
				cout << endl;
				for( int i=1; i <= mM.Nrows(); ++i ){
					for( int j=1; j <= mM.Ncols(); ++j )
						cout << left << setw(iWidth) << setprecision(iPrecision) << mM(i,j) << "\t";
					cout << endl;
				}
			}
			
			else
				
			{		
				ofstream outputFile;
				if(fileMode == "app")
					outputFile.open( outDevice.c_str(), ios:: app );
				else
					outputFile.open( outDevice.c_str() );
				for(vector<string>::size_type sv = 0; sv != vLables.size(); ++ sv)
					outputFile << left <<  setw(iWidth+1) << vLables[sv] << "\t" ;
				outputFile << endl;
				for( int i=1; i <= mM.Nrows(); ++i ){
					for( int j=1; j <= mM.Ncols(); ++j )
						outputFile << left << setw(iWidth) << setprecision(iPrecision) << mM(i,j) << "\t";
					outputFile << endl;
				}
				outputFile.close();		
			}
			
			break;
			
		default:
			if( outDevice == "screen"){
				for(vector<string>::size_type sv = 0; sv != vLables.size(); ++ sv)
					cout << left <<  setw(iWidth+1) << vLables[sv] ;
				cout << endl;
				cout << setw(iWidth) << setprecision(iPrecision) << mM << endl;
			}
			
			else
				
			{		
				ofstream outputFile;
				if(fileMode == "app")
					outputFile.open( outDevice.c_str(), ios:: app );
				else
					outputFile.open( outDevice.c_str() );
				for(vector<string>::size_type sv = 0; sv != vLables.size(); ++ sv)
					outputFile << left <<  setw(iWidth+1) << vLables[sv] ;
				outputFile << endl;
				outputFile << setw(iWidth) << setprecision(iPrecision) << mM << endl;
				outputFile.close();		
			}
			
			break;
			
	}
	
}

void print_table( Matrix mM, int iWidth, int iPrecision, string outDevice, const char sep, string fileMode ){
	
	switch (sep) {
			
		case 't' : 
			if( outDevice == "screen"){
				for( int i=1; i <= mM.Nrows(); ++i ){
					for( int j=1; j <= mM.Ncols(); ++j )
						cout << left << setw(iWidth) << setprecision(iPrecision) << mM(i,j) << "\t";
					cout << endl;
				}
			}
			
			else
				
			{		
				ofstream outputFile;
				if(fileMode == "app")
					outputFile.open( outDevice.c_str(), ios:: app );
				else
					outputFile.open( outDevice.c_str() );
				
				for( int i=1; i <= mM.Nrows(); ++i ){
					for( int j=1; j <= mM.Ncols(); ++j )
						outputFile << left << setw(iWidth) << setprecision(iPrecision) << mM(i,j) << "\t";
					outputFile << endl;
				}
				outputFile.close();		
			}
			
			break;
			
			
		default : 
			
			
			if( outDevice == "screen")
				cout << setw(iWidth) << setprecision(iPrecision) << mM << endl;
			
			else
				
			{		
				ofstream outputFile;
				if(fileMode == "app")
					outputFile.open( outDevice.c_str(), ios:: app );
				else
					outputFile.open( outDevice.c_str() );
				outputFile << setw(iWidth) << setprecision(iPrecision) << mM << endl;
				outputFile.close();		
			}
			
			break;
			
	}
	
}

void print_table( vector< vector<double> > dTable, string outDevice, string fileMode ){
	
	if( outDevice == "screen"){
		for( vector< vector<double> >::iterator il= dTable.begin(); il!= dTable.end(); ++il ){
			for( vector<double>::iterator im= (*il).begin(); im!= (*il).end(); ++im )
				cout << left << (*im) << "\t" ;
			cout << endl;
		}
	}
	
	else
		
	{		
		ofstream outputFile;
		if(fileMode == "app")
			outputFile.open( outDevice.c_str(), ios:: app );
		else
			outputFile.open( outDevice.c_str() );
		for( vector< vector<double> >::iterator il= dTable.begin(); il!= dTable.end(); ++il ){
			for( vector<double>::iterator im= (*il).begin(); im!= (*il).end(); ++im )
				outputFile << left << (*im) << "\t" ;
			outputFile << endl;
		}			
		outputFile.close();		
	}
	
}

ReturnMatrix ReadTable(size_t nRows, size_t nCols, string inputFileName){
	
	double dData[nRows*nCols];
	Matrix mM(nRows, nCols);
	
	ifstream inputFile( inputFileName.c_str() ); 
	for(size_t i=0; i<(nRows*nCols); ++i)	inputFile >> dData[i];
	inputFile.close();
	
	mM << dData;
	
	mM.Release();	return mM;
	
}


ReturnMatrix compare_prices( const ColumnVector &mPrices, string inputFileName ){
	
	int np = mPrices.Nrows();
	Matrix mErr(np, 3);
	
	// read from file
	
	ifstream inputFile ( inputFileName.c_str() ); 
	for(int i=0; i<np; ++i){
		inputFile >> mErr.element(i,0) ;
	}
	inputFile.close();
	
	
	// Compute errors
	for( int l=0; l < np; ++l ){
		mErr.element(l,1) = 10000 * ( mPrices.element(l) - mErr.element(l,0) );
		mErr.element(l,2) = 10000 * abs( mPrices.element(l) / mErr.element(l,0) - 1 ); 
	}
	
	return mErr;
}

ReturnMatrix compare_prices( const ColumnVector &mPrices, const ColumnVector &mTruePrices ){
	
	int np = mPrices.Nrows();
	Matrix mErr(np, 3);

	mErr.Column(1) = mTruePrices;

	// Compute errors
	for( int l=0; l < np; ++l ){
		mErr.element(l,1) = 10000 * ( mPrices.element(l) - mErr.element(l,0) );
		mErr.element(l,2) = 10000 * abs( mPrices.element(l) / mErr.element(l,0) - 1 ); 
	}
	
	return mErr;
}



