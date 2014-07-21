#include<iostream>
#include<vector>
#include<string>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_rmultinom_weight_v2 (const SEXP proba_matrix, const SEXP z_matrix, const SEXP seed, const SEXP weights_col, const SEXP weights_row);
}



SEXP C_rmultinom_weight_v2 (const SEXP proba_matrix, const SEXP seed, const SEXP weights_col, const SEXP weights_row, const SEXP sumCols) {

  const unsigned int  Seed = INTEGER(seed)[0];
  srand((unsigned)Seed); 

  const double * pmatrix = REAL(proba_matrix);
  double * cweights = REAL(weights_col);
  double * rweights = REAL(weights_row);
  const int * dims = INTEGER(getAttrib(proba_matrix, R_DimSymbol));

  int nrow = dims[0];
  int ncol = dims[1];

  float * norm = new float [ ncol ];
  
  /////////// This will store the sum of the counts
  double * sumCol_c = REAL( sumCols );
  for (int i = 0; i != ncol; i++) { sumCol_c[ i ] = 0;}


  for (int i = 0; i != nrow; i++) {    

    float sum = 0.;
    for (int j = 0; j != ncol; j++) {
      norm[ j ] =  pmatrix[ j*nrow + i ]*cweights[ j ];
      sum += norm[ j ];
    }
    
    float randN = sum*(rand() / float(RAND_MAX));
    
    int k = 0;
    float current = norm[ 0 ];   //probability for species 1
    while (current < randN) {k++;current += norm[ k ];}  //move to the next species in an iterative manner
    sumCol_c[ k ] += rweights[ i ];
  }
  
  delete [] norm;

  return (sumCols);
}

