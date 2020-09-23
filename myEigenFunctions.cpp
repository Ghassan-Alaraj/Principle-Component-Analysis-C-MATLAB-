/*****************************************************************************************************

 This file is where you'll put the source information for the functions you write.  Make sure it's
 included in your project (shift-alt-A) and that the functions you add are declared in the header
 file, myEigenFunctions.h.  You can add any extra comments or checks that could be relevant to these
 functions as you need to.

====================================================================================================*/

#include <string>

#include "myEigenFunctions.h"


using namespace std;

double DotProduct(double **A, double **B, int n, int m)
{
	// This is a function that takes two matrices A and B of identical dimensions (n*m) and
	//  calculates and returns their dot product.
	//
	double dot = 0.0;
/*  */
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			dot += A[i][j]*B[i][j];
		}
	}
	return dot;
}

double DotProduct(double *A, double *B, int n)
{
	//
	//	This is a function that takes two vectors A and B of identical length (n) and
	//  calculates and returns their dot product.
	//

	double dot = 0.0;

	for(int i = 0; i < n; i++) {
		dot += A[i]*B[i];
	}
	return dot;
}

double* DotProduct(double **A, double *v, int n)
    {
          //
          //  This is a function that takes a nxn-matrix A and an n-dimensional vector v stores
          //  the product A.v at the original location of v
          //

    double *result = new double[n]; // pointer to result vector

    for(int i = 0; i < n; i++) {
      result[i] = 0.0; // initialize ith element of result v
      for(int j = 0; j <n; j++){
        result[i] += A[i][j] * v[j];
      }
    }

    return result;
  }



double** ReadData(string inputFileName, int n, int m)
    {
          //
          //  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
          //  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
          // The code for parsing the CSV file has been adapted from Bernardo Trindade
          // https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c

      double **spectra = new double* [n];


      for(int sample = 0; sample < n; sample++)
	spectra[sample] = new double [m];

    vector<vector<double> > data;
    cout << inputFileName << endl;
    ifstream inputFile(inputFileName);
    int l = 0;

    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
	// cout << s << endl;
        if (l!=0) // ignore first line, which contains wavelengths
	  {
            istringstream ss(s);
            vector<double> record;

	    int column = 0;
            while (ss) {
	      // cout << "Row " << l << " " << ss << endl;
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
		  // cout << "Row " << l << " Column " << column << " line " << line << endl;
                    spectra[l-1][column] = stod(line);
		    column++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
		}

        }
	l++;
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return spectra;
  }

void print_matrix(double **A, int n, int m)
{
  //
  // This procedure prints an nxm matrix.
  //
  for (int i=0; i < n; i++)
    {
    for (int j=0; j < m; j++)
      std::cout << A[i][j] << '\t';
    std::cout << '\n';
    }
}

Eigenpair power_method(double **A, double *v, int n, double tol)
{
  /*
      ''' This function computes the largest eigenvalue and the corresponding eigenvector

        Parameters:
        -----------
        A : n by n matrix
            matrix
        v : vector
            intial eigenvector guess
        n : integer
            length of A/v matrix/vector
        tol : float
            tolerance for error

        Returns:
        --------
        eigenvector : struct
            contains the eigenpair


    ''' */


    int k, index;
    double lam[1000];
    double val = 0.0;

    double *t = NULL;
    t = new double [n];

    double *Av;
    Av = new double [n];


    Eigenpair eigenpair(n);

    // intilize the eigenvector
    for(int i = 0; i < n; i++)
        eigenpair.vector[i] = v[i];

    //norm the vector
    eigenpair.normalize();

    //intialize
    lam[0] = 1;
    lam[1] = eigenpair.value;
    k = 1;

    while((abs(lam[k-1]-lam[k])) > tol){
        k++;

        t = DotProduct(A, v, n);

        //argmax
        for(int i = 0; i < n; i++){
            if(t[i] > val){
                index = i;
                val = t[i];
            }
        }
        //calculate norm of t vector
        lam[k] = abs(sqrt(DotProduct(t, t, n)));

        //norm guess vector
        for (int i =0; i < n; i++){
            v[i] = t[i]/lam[k];
        }




    }

    //negativity test
    Av = DotProduct(A,v,n);

    if (Av[index] / v[index] < 0)
        lam[k] = lam[k] * -1;

    //record results within the eigenpair struct
    for(int i = 0; i < n; i++)
        eigenpair.vector[i] = v[i];

    eigenpair.normalize();
    eigenpair.value = lam[k];



  return eigenpair;
}

void deflate(double **A, Eigenpair eigenpair)
{
    /*
    '''
        This procedure removes eigenpair.vector from transformation matrix A in place

        Parameters:
        -----------
        A : n by n matrix
            matrix
        eigenpair : struct
            struct containing the eigen pair

        Returns:
        --------
        None

    ''' */
    //record the eigenvector length
    int n = eigenpair.length;
    double R[n][1];
    double R_trans[1][n];
    double p[n][n];

    //record the eigenvalue
    double lam = eigenpair.value;


    //NOTE: the eigenvector is normalized here !!
    eigenpair.normalize();

    for(int i =0;i<n;i++){
        R[i][0] = eigenpair.vector[i];
        R_trans[0][i] = eigenpair.vector[i];
    }

    //preform dot product
    for(int i = 0; i < n;i++){
        for(int j = 0; j < n; j++){
            p[i][j] = 0;
            for(int k =0;k <1;k++){
                p[i][j] += R[i][k] * R_trans[k][j];
            }

        }
    }
    //deflate the matrix
    for(int i = 0; i < n;i++){
        for(int j = 0; j < n; j++){
            A[i][j] -= (lam * p[i][j]);

        }
    }



}

double** CenterMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and subtracts the emperical mean of each column.
          //  The function returns the point to the centered matrix
          //


      double **result = new double* [n];

      for(int row = 0; row < n; row++)
            result[row] = new double [m];

      double *means = NULL;
      means = new double [m];
      double temp;

      //first lets find the sample mean or the emperical mean for each column and
      //store it in a vector.
      for (int j =0 ; j < m; j++){
        means[j] = 0;
        temp = 0;
        for(int i = 0;i < n; i++){
            temp += A[i][j];
        }
        means[j] = (temp/n);
      }

      //column mean minus observation and record the result.
      for (int j =0 ; j < m; j++){
        for(int i = 0;i < n; i++){
            result[i][j] = A[i][j] - means[j];
        }
      }

    return result;
  }


double** CovarianceMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and subtracts the emperical mean of each column.
          //  The function returns the point to the centered matrix
          //




      double **cov = new double* [m];

      double cAT[m][n];
      double **cA = CenterMatrix(A, n, m);

      for(int i = 0; i < m; i++){
            cov[i] = new double [m];
      }


    // Computing transpose of the matrix
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) {
            cAT[j][i] = cA[i][j];
      }



    for(int i =0;i<m;i++){
        for(int j=0;j<m;j++){
            cov[i][j] = 0;
        }
    }


    //preform dot product
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < m; ++j)
            for(int k = 0; k < n; ++k)
            {
                cov[i][j] += cAT[i][k] * cA[k][j];
            }

      // 1/(num of rows) - 1
      double R = pow((n - 1),-1);


      for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) {
            cov[i][j] *= R;
      }


    return cov;
  }
