/*
  Testing the calculation of the derivative and checking the error.
  
  Compile:

    export OMP_NUM_THREADS=<N>    N - is number of threads
    export OMP_NESTED=TRUE
  
    make test_diff
    
  Gnuplot:
    
    set log
    set format '10^{%L}'
    plot '< ./test_diff' u 1:(abs($2)) w l, 1e-17/x + x**6/140
    
  Author: Martin Horvat, November 2013
*/ 
#include <iostream>
#include <cmath>
#include <omp.h>

typedef double myreal;

#include "num_diff.h"


template <class T>
T func(T *x){
  //return std::sin(x[0]);
  return x[0]*x[0];
}

template <class T>
T dfunc(T *x){
  //return std::cos(x[0]);  
  return 2*x[0];
}

int main(){

  const int d = 1;

  myreal x[d] = { 0.2 },
         h[d] = { 0.1 },
         df[d], df0[d];

  Tnum_diff type = central_7;


  #if defined(_OPENMP)
  std::cout 
    << "OpenMP: Max. threads:" 
    << omp_get_max_threads() << '\n';  // OMP_NUM_THREADS
  #endif
  
  std::cout.precision(std::numeric_limits<myreal>::digits10+1);
  std::cout << std::scientific;

  for (int i = 0; i < 200; ++i){
    
    num_diff(type, d, h, x, func, df);
    num_diff_kahan(type, d, h, x, func, df0);
    
    std::cout 
      << h[0] << ' ' 
      << df[0] - dfunc(x) << ' '
      << df[0] - df0[0] << '\n';
    
    h[0] *= 0.9;
  }
  
  //
  // estimate optimal step
  //
  h[0] = 0.1;
  
  num_diff_estimate_step(type, d, h, x, func);
   
  for (int i = 0; i < d; ++i)
    std::cout << "#h[" << i << "]=" << h[i] << '\n'; 
  
  return EXIT_SUCCESS;
}
