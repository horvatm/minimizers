/*
  Testing Brent's MINFIT algorithm which computes the singular value 
  decomposition and compare C++ and fortran implementation.
  
  Compile: 
    make test_minfit
  
  Author: Martin Horvat, June 2013
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

#include "praxis.h"
#include "matrix.h"

typedef double myreal;

extern "C" void minfit_fort_ (int *m, int *n, double *tol, double *ab, double *q);

int main(){
  
  int n = 4,
      i, j;
  
  myreal tol = 1e-16, 
         *ab = new myreal [n*n], 
         *abF = new myreal [n*n],
         *q = new myreal [n],
         *qF = new myreal [n],
         *work = new myreal [n];
  
  
  std::cout.precision(std::numeric_limits<myreal>::digits10+1);
  std::cout << std::scientific;
    
  srand(0);
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      abF[i+n*j] = ab[i + n*j] = rand();
  
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) std::cerr << int(ab[i + n*j]) << ' ';
    std::cerr << '\n';
  }
  
  minfit_(n, tol, ab, q, work);
  
  minfit_fort_(&n, &n, &tol, abF, qF);
  
  sort_(n, q, ab);
  
  sort_(n, qF, abF);
  
  for (i = 0; i < n; ++i) std::cout << q[i] << ' '<<  qF[i] << '\n';
  
  std::cout << '\n';
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) std::cout << ab[i + n*j] << ' ';
    std::cout << '\n';
  }
  
  return 0;
}
