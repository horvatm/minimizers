/*
  Testing powell's + brent's algorithm for minimization of functions.
  
  Compile:
    make test_powell
    
  Author: Martin Horvat, June 2013
*/ 
 
#include <cmath>
#include <iostream>
#include <cstdlib>

typedef double myreal;

#include "powell.h"
#include "matrix.h"


void print(std::string label, int nfunk, myreal *v, const myreal & y){
  
  std::cerr 
    << label 
    << ": Nr. eval:" << nfunk 
    << " Vec:";
  
  for (int i = 0; i < 3; ++i) std::cerr << ' ' << v[i];
  
  std::cerr << " Value:" <<  y << '\n';
}

myreal f(myreal *x) {
  
  myreal xx = x[0] - 1,
         yy = x[1] - 2,
         zz = x[2] - 3;
  
  return 4 + (xx*xx + xx*yy + yy*yy + zz*zz);//*std::abs(std::sin(xx));
}

int main(){
  
  int i, j, 
      nfunk,
      iter, 
      n = 3,
      seed = 0;
  
  myreal fret, 
         tol[2] = {1e-16,  1e-16},
         *p = new myreal [n], 
         **xi = matrix <myreal> (n, n);          
  
  srand(seed);
  
  // initial point
  for (i = 0; i < n; ++i) p[i] = 2*myreal(rand())/RAND_MAX-1;
  
  // directions
  for (i = 0; i < n; ++i) {

    for (j = 0; j < n; ++j) xi[j][i] = 0;
    
    xi[i][i] = 1;
  }
  
  std::cout.precision(std::numeric_limits<myreal>::digits10+1);
  std::cout << std::scientific;
      
  std::cerr.precision(std::numeric_limits<myreal>::digits10+1);
  std::cerr << std::scientific;
  
  std::cout << "Init:\n";
  
  
  std::cout  << " Value:" << f(p) << " Vec:";

  for (i = 0; i < n; ++i) std::cout << ' ' << p[i];

  std::cout << std::endl;
    
  
  if (!powell(p, xi, n, tol, fret, f, iter, nfunk, print)) {
    std::cerr << "Powell failed to converge\n";
    return EXIT_FAILURE;
  } 
  
  std::cout << "------------------\nResult:\n";
  
  
  std::cout << " Value:" << fret << " Vec:";
   
  for (i = 0; i < n; ++i) std::cout << ' ' << p[i];
     
  std::cout << '\n';
  
  std::cout << " Iter:" << iter << " Nr F eval:" << nfunk << '\n';
  
  
  free_matrix (xi);
  
  delete [] p;
  
  return EXIT_SUCCESS;
}
