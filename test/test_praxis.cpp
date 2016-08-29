/*
  Testing Brent's Praxis algorithm for minimization of functions.
  
  Compile: 
    make test_praxis
    
  Author: Martin Horvat, June 2013
*/ 
 
#include <cmath>
#include <iostream>
#include <cstdlib>

typedef double myreal;

#include "praxis.h"
#include "matrix.h"

void print(std::string label, int nfunk, myreal *v, const myreal & y){
  
  std::cerr 
    << label 
    << ": Nr. eval:" << nfunk 
    << " Vec:";
  
  for (int i = 0; i < 3; ++i) std::cerr << ' ' << v[i];
  
  std::cerr << " Value:" <<  y << '\n';
}

myreal sqr(const myreal &x){return x*x;}

myreal x0[3] = {1,2,3};

myreal f(myreal *x) {
  
  myreal xx = x[0] - x0[0],
          yy = x[1] - x0[1],
          zz = x[2] - x0[2];
  
  return 4 + (xx*xx + 2*xx*yy + 3*yy*yy + 4*zz*zz);//*std::abs(std::sin(xx));
}

int main(){
  
  int i,
      ret, 
      nfunk,
      n = 3,      // dimension
      
      h0 = 2,     // maximum step size
      seed = 3;
  
  myreal fret,        // returned value
         t0 = 1e-8,   // tolerance (abs. precision, accuracy)
         x[3] = {0.1, 0.2, 1};
  
  std::cout.precision(std::numeric_limits<myreal>::digits10+1);
  std::cout << std::scientific;
      
  std::cerr.precision(std::numeric_limits<myreal>::digits10+1);
  std::cerr << std::scientific;
  
  std::cout << "Init:\n";
  
  std::cout  << " Value:" << f(x) << " Vec:";
  
  for (i = 0; i < n; ++i) std::cout << ' ' << x[i];

  std::cout << std::endl;
   
  if ((ret = praxis<myreal> (seed, t0, h0, n, x, f, fret, nfunk, print)) != 0) {
    std::cerr << "Praxis failed to converge:" << ret << '\n';
    return EXIT_FAILURE;
  } 
  
  std::cout << "------------------\nResult:\n";
  
  
  std::cout << " Value:" << fret << " Vec:";
   
  for (i = 0; i < n; ++i) std::cout << ' ' << x[i];
     
  std::cout << '\n';
  
  std::cout << " Nr F eval:" << nfunk << '\n';
  
  myreal sum = 0, sum1 = 0;
  for (i = 0; i < n; ++i) {
    sum += sqr(x[i] - x0[i]);
    sum1 += sqr(x0[i]);
  }
  
  std::cout << " EPS:" << sqrt(sum) << " " << sqrt(sum1) << '\n';
    
  return EXIT_SUCCESS;
}

