/*
  Testing BFGS algorithm for minimization of functions.
   
  Using calculation of derivative via num_diff.h
  
  export OMP_NUM_THREADS=<N>    N - is number of threads
  export OMP_NESTED=TRUE
  
  Compile: 
    make test_bfgs_diff

  Author: Martin Horvat, November 2013
*/ 

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <omp.h>

typedef double myreal;

#include "bfgs.h"
#include "num_diff.h"
#include "matrix.h"

Tnum_diff type = central_7; 

myreal step[3] = {myreal(1)/10, myreal(1)/10, myreal(1)/10};
  
template <class T>
void print(int iter, int nfunk[2], T *v, const T & y){
  
  std::cerr 
    << " Iter:" << iter
    << " Nr. fun:" << nfunk[0] << ":" << nfunk[1]
    << " Vec:";
  
  for (int i = 0; i < 3; ++i) std::cerr << ' ' << v[i];
  
  std::cerr << " Value:" <<  y << '\n';
}

template <class T>
T func(T *x) {
  
  T xx = x[0] - 1,
    yy = x[1] - 2,
    zz = x[2] - 3;
  
  return 4 + (xx*xx + xx*yy + yy*yy + zz*zz);
}
template <class T>
void dfunc(T *x, T *df){
  
  T xx = x[0] - 1,
    yy = x[1] - 2,
    zz = x[2] - 3;
    
  df[0] = 2*xx + yy;
  df[1] = 2*yy + xx;
  df[2] = 2*zz;
}


template <class T>
void dfunc_num(T *x, T *df){
  
  T h[3];
  
  for (int i = 0; i < 3; ++i) h[i] = step[i];
  
  //num_diff_kahan(type, 3, h, x, func, df);  // serial + Kahan summation
  //num_diff(type, 3, h, x, func, df);        // serial
  //num_diff_omp(type, 3, h, x, func, df);    // openmp nested  
  num_diff_omp_f(type, 3, h, x, func, df);    // openmp function calls
}

int main(){
  
  int i,
      nfunk[2],
      iter, 
      n = 3,
      seed = 0;
  
  myreal fret, 
         gtol = std::numeric_limits<myreal>::epsilon(),
         stepmax = 1,
         p[3];
  
  #if defined(_OPENMP)
  std::cout 
    << "OpenMP: Max. threads:" 
    << omp_get_max_threads() << '\n';  // OMP_NUM_THREADS
  #endif
  
  srand(seed);

  std::cout.precision(std::numeric_limits<myreal>::digits10+1);
  std::cout << std::scientific;
      
  std::cerr.precision(std::numeric_limits<myreal>::digits10+1);
  std::cerr << std::scientific;
  
  // initial point
  for (i = 0; i < n; ++i) p[i] = 2*myreal(rand())/RAND_MAX-1;
  
  std::cout << "Original h[]:\n";
  for (i = 0; i < n; ++i) std::cout << step[i] << '\n';
  
  //num_diff_estimate_step(type, n, step, p, func);
  //num_diff_estimate_step_omp(type, n, step, p, func);
  num_diff_estimate_step_omp_f(type, n, step, p, func);
  
  std::cout << "Optimal h[]:\n";
  for (i = 0; i < n; ++i) std::cout << step[i] << '\n';
    
  // value and derivative
  myreal df1[n], df0[n];
  
  dfunc(p, df0);
  dfunc_num(p, df1);
  
  std::cout << "\nNum. Derivative:" << std::endl;
  for (i = 0; i < n; ++i) 
    std::cout << " err_df[" << i << "]=" << df1[i] - df0[i] << '\n';
  
  std::cout << "\nInit:\n Value:" << func(p) << " Vec:";

  for (i = 0; i < n; ++i) std::cout << ' ' << p[i];

  std::cout << "\n\nSimulation:" << std::endl;
    
  int nerr = bfgs(p, n, gtol, stepmax, fret, func, dfunc_num, iter, nfunk, print);
  
  std::cout << "\nStatus:" << std::endl;
  switch (nerr) {
    case -1: 
      std::cout << "BFGS: no error: zeroing gradient\n"; 
      break;
    case 0: 
      std::cout << "BFGS: no error: converged in x\n"; 
      break;
    case 1: 
      std::cout << "BFGS: maximal number of iterations exceeded \n";
      return EXIT_FAILURE;
    case -2:
      std::cout << "BFGS: roundoff problem in lnsrch\n"; 
      return EXIT_FAILURE;
  }
  
  std::cout << "\nResults:\n Value:" << fret << " Vec:";
   
  for (i = 0; i < n; ++i) std::cout << ' ' << p[i];
     
  std::cout << "\n Iter:" << iter 
            << " Nr F eval:" << nfunk[0] << ":" << nfunk[1] << '\n';
  
  return EXIT_SUCCESS;
}
