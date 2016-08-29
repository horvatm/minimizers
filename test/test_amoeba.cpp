/*
  Testing amoeba/downhill simplex algorithm for minimization of functions.
  
  Compile: 
    make test_amoeba

  Author: Martin Horvat, May 2013
*/ 
 
#include <cmath>
#include <iostream>
#include <cstdlib>

typedef double myreal;

#include "amoeba.h"
#include "matrix.h"

void print(int nfunk, myreal *v, const myreal & y){
  
  std::cerr 
    << " Nr. eval:" << nfunk 
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
  
  int i, j, nfunk, 
      ndim = 3,
      seed = 0;
  
  myreal ftol = 1e-12,  //fractional tolerance 
         dx = 0.1,
         *y = new myreal [ndim+1],
         **p = matrix <myreal> (ndim +1, ndim);
  
  srand(seed);
  
  // initial point
  for (i = 0; i < ndim; ++i) p[0][i] = myreal(rand())/RAND_MAX;
  y[0] = f(p[0]);
  
  // neighboring points around initial point -- creating a simplex
  for (i = 1; i <= ndim; ++i) {

    for (j = 0; j < ndim; ++j) p[i][j] = p[0][j];
    
    p[i][i-1] += dx;
    
    y[i] = f(p[i]);
  }
  
  std::cout.precision(16);
  std::cout << std::scientific;
  
    
  std::cerr.precision(16);
  std::cerr << std::scientific;
  
  std::cout << "Init:\n";
  for (i = 0; i <= ndim; ++i){
     std::cout 
      << "Nr: " << i 
      << " Value:" << y[i]
      << " Vec:";
     
     for (j = 0; j < ndim; ++j) std::cout << ' ' << p[i][j];
     
     std::cout << '\n';
  }
  
  if (!amoeba<myreal>(p, y, ndim, ftol, f, nfunk, print)) {
    std::cerr << "Amoeba failed to converge\n";
    return EXIT_FAILURE;
  } 
  
  std::cout << "------------------\nResult:\n";
  for (i = 0; i <= ndim; ++i){
     std::cout 
      << "Nr: " << i 
      << " Value:" << y[i]
      << " Vec:";
     
     for (j = 0; j < ndim; ++j) std::cout << ' ' << p[i][j];
     
     std::cout << '\n';
  }
  
  std::cout << "Nr. F eval =" << nfunk << '\n';
  
  free_matrix (p);
  
  delete [] y;
  
  return EXIT_SUCCESS;
}
