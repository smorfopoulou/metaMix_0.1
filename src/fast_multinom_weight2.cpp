#include<iostream>
#include<vector>
#include<string>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_rmultinom_weight2 (const SEXP proba_matrix, const SEXP z_matrix, const SEXP seed, const SEXP weights);
}



SEXP C_rmultinom_weight2 (const SEXP proba_matrix,  const SEXP z_matrix, const SEXP seed, const SEXP weights) {

  const unsigned int  Seed = INTEGER(seed)[0];
  srand((unsigned)Seed); 

  const double * pmatrix = REAL(proba_matrix);
  double * zmatrix =  REAL(z_matrix);
  double * cweights = REAL(weights);
  const int * dims = INTEGER(getAttrib(proba_matrix, R_DimSymbol));

  int nrow = dims[0];
  int ncol = dims[1];

  float * norm = new float [ncol];
  
  //cout<<"Number of rows "<<nrow<<" and nb of columns "<<ncol<<endl;

  for (int i = 0; i != nrow; i++) {    

    float sum = 0.;
    for (int j = 0; j != ncol; j++) {
      norm[ j ] =  pmatrix[ j*nrow + i ]*cweights[ j ];
      sum += norm[ j ];
    }
    
    //    float randN = sum*(rand() / float(RAND_MAX));
    
    float randN = rand() / float(RAND_MAX);
    int k = 0;

    //    float current = norm[ 0 ] ;   //probability for species 1
    float current = norm[ 0 ] / sum ;   //probability for species 1
    //    while (current < randN) {k++;current += norm[ k ];}  //move to the next species in an iterative manner

        while (current < randN) {k++;current += norm[ k ]/sum;}  //move to the next species in an iterative manner
    
    // now we assign the result to the output matrix
    for (int j = 0; j != ncol; j++) {
      if (j == k) {zmatrix[  j*nrow + i ] = 1;} else zmatrix[  j*nrow + i ] = 0;
    }

  }

  return z_matrix;
}

