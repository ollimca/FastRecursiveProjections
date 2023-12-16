/*
 *  MatrUtils.h
 *  WavHeston
 *
 *  Created by Antonio Cosma on 25/February/2011.
 *  Copyright 2011 Université du Luxembourg. All rights reserved.
 *
 */

#include "newmatio.h"
#include "Seq.h"
#include "cmath"
#include<fstream>

#ifndef _MATRUTILS_
#define _MATRUTILS_

//using namespace NEWMAT;
using namespace std;


int find_index(const Matrix &mX, const Real X);
int find_closest_index(const Matrix &mX, const Real X);

void print_table( vector< string > vLables, Matrix mM, int iWidth, int iPrecision, string outDevice = "screen", const char sep = 'm', string fileMode = " " );
void print_table( Matrix mM, int iWidth, int iPrecision, string outDevice = "screen", const char sep =  'm', string fileMode = " " );
void print_table( vector< vector<double> > dTable, string outDevice = "screen", string fileMode = " " );
ReturnMatrix ReadTable(size_t nRows, size_t nCols, string inputFileName);
ReturnMatrix compare_prices( const ColumnVector &mPrices, string inputFileName );
ReturnMatrix compare_prices( const ColumnVector &mPrices, const ColumnVector &mTruePrices );


#endif

//print_table( vector<string> vLables, Matrix mM, int iWidth, int iPrecision, string outDevice, string sep, string fileMode )
