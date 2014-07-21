#include<iostream>
#include<vector>
#include<string>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_rmultinom (const SEXP proba_matrix, const SEXP z_matrix, const SEXP seed);
}



SEXP C_rmultinom (const SEXP proba_matrix,  const SEXP z_matrix, const SEXP seed) {

  const unsigned int  Seed = INTEGER(seed)[0];
  srand((unsigned)Seed); 

  const double * pmatrix = REAL(proba_matrix);
  double * zmatrix =  REAL(z_matrix);

  const int * dims = INTEGER(getAttrib(proba_matrix, R_DimSymbol));

  int nrow = dims[0];
  int ncol = dims[1];

  //cout<<"Number of rows "<<nrow<<" and nb of columns "<<ncol<<endl;

  for (int i = 0; i != nrow; i++) {    

    double sum = 0.;
    for (int j = 0; j != ncol; j++) sum += pmatrix[ j*nrow + i ];
    
    double randN = rand() / double(RAND_MAX);
    
    int k = 0;
    double current = pmatrix[ k*nrow + i ] / sum;
    while (current < randN) {k++;current += pmatrix[ k*nrow + i ]/sum;}
    
    for (int j = 0; j != ncol; j++) {
      if (j == k) {zmatrix[  j*nrow + i ] = 1;} else zmatrix[  j*nrow + i ] = 0;
    }

  }

  return z_matrix;
}

