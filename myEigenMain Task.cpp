
#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <fstream>
#include <string>

#include "myEigenFunctions.h"

using namespace std;

#define PI 3.14159265358979323846
#define TOL 1E-6


int main(void)
{
  // Defining local variables to be used:

  double **spectra;

  const int N = 56;   // rows
  const int M = 286;  // columns
  Eigenpair eigenPP(M);
  double *v_long = NULL;
  v_long = new double[M];

  double allvaules[56];
  double allvectors[56][M];

  spectra = ReadData("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv", N, M);


  double **cov  = CovarianceMatrix(spectra, N, M);

  //Output eigenpair as two csv files

 //setup both vectors with zeros to avoid any unforseen issues later
  for(int i = 0;i <56;i++){
        allvaules[i] = 0.0;
    for(int j = 0;j < M;j++){
        allvectors[i][j] = 0.0;

    }
  }

  //start looping over the cov matrix to find each eigenpair.//
  for (int i = 0; i < 56;i++){

  //initialize / clean the eigenPP and v_long variable
  for(int j = 0; j < M + 1; j++){
    eigenPP.vector[j] = 1.0;
    v_long[j] = 1.0;
  }

  eigenPP.normalize();

  //get the eigenpair using the power method
  eigenPP = power_method(cov,v_long,M,TOL);




  //record infromation
  allvaules[i] += eigenPP.value;
  for(int k =0;k < M + 1; k++){
    allvectors[i][k] += eigenPP.vector[k];

  }

  //deflate cov matrix
  deflate(cov,eigenPP);




  }
  //Note: both csv files contain a column for indexing, this is done to ensure that the data is recorded in the
  //proper order and is removed when imported into matlab.

  //save eigenValues
  ofstream myfile;
  myfile.open ("eigenValues.csv");
  myfile << "index,";
  myfile << "eigenValue\n";
  for (int i = 0;i < 56; i++){
    myfile << i+1;
    myfile << ",";
    myfile << allvaules[i];
    myfile << "\n";
  }
  myfile.close();

  //save eigenVectors
  //NOTE: each row represents an eigenvector.

  myfile.open ("eigenVectors.csv");
  myfile << "index,";
  myfile << "eigenVectors\n";
  for (int i = 0;i < 56; i++){
    myfile << i+1;
    myfile << ",";
    for (int j = 0; j < M; j++) {
        myfile << allvectors[i][j];
        myfile << ",";
    }
    myfile << "\n";
  }
  myfile.close();

  return 0;

}
// end main

