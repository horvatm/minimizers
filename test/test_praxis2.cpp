/*
  Testing Brent's Praxis algorithm for minimization of functions.
  
  Compile: 
    g++ test_praxis.cpp -o test_praxis -O2 -Wall
    icc test_praxis.cpp -o  test_praxis -Qoption,cpp,--extended_float_type -DQUAD_TYPE -O2 -lifcore ../intel_quad/libquad/libquad.a
    
  Author: Martin Horvat, June 2013
*/ 
 
#include <cmath>
#include <iostream>
#include <cstdlib>

typedef double myreal;

#include "praxis.h"
#include "matrix.h"

int n = 2,      // dimension
    seed = 3,
    ret, 
    nfunk;

myreal fret,            // returned value
       t0 = 1e-8,       // tolerance (abs. precision, accuracy)
       h0 = 1,          // maximum step size
       x[2] = {-1.2, 1},
       x0[2]= {1, 1};
       
void print(std::string label, int nfunk, myreal *v, const myreal & y){
  
  std::cerr 
    << label 
    << ": Nr. eval:" << nfunk 
    << " Vec:";
  
  for (int i = 0; i < n; ++i) std::cerr << ' ' << v[i];
  
  std::cerr << " Value:" <<  y << '\n';
}

template <class T> T sqr(const T & x) {return x*x;}

myreal f(myreal *x) {
  
  return sqr(1-x[0]) + 100*sqr(x[1]-sqr(x[0]));
}

int main(){
 
  int i;
  
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

