#if !defined(__amoeba_h)
#define __amoeba_h

/*
  Downhill simplex (or Nelder Mead method) optimization
  
  Using source code given in NRC.
  
  Ref: 
   * http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
   * NRC 2nd ed. - 10.4 Downhill Simplex Method in Multidimensions called amoeba
     File: c10-4.pdf
  
  Author: Martin Horvat, May 2013
*/


#include <cmath>
#include <limits>

/*
Extrapolates by a factor fac through the face of the simplex across from the high point, tries it, and replaces the high point if the new point is better.*/
 
template <class T> 
  T amotry(T **p, T *y, T *psum, const int & n,
           T func(T *), const int &ihi, const T &fac, T *ptry){
  int j;
  
  T fac1 = (1 - fac)/n, fac2 = fac1 - fac, ytry;
         
  for (j = 0; j < n; ++j) ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
  ytry = func(ptry);
  
  if (ytry < y[ihi]) {
    y[ihi] = ytry;
    for (j = 0; j < n; ++j) {
      psum[j] += ptry[j] - p[ihi][j];
      p[ihi][j] = ptry[j];
    }
  }
  
  return ytry;
}
  
/*
  Multidimensional minimization of the function funk(x) where x[0..n-1] is a vector in n dimensions, by the downhill simplex method of Nelder and Mead. 

  Input:
    p[0..n][0..n-1] - starting simplex  (n + 1 vectors of dimension n)
    y[n+1] = {f(p_i): p_i in T^n, i = 0, ..., n} 
              - values of func evaluated at the n+1 vertices (rows) of p;
    n - dimension of the domain
    ftol - fractional convergence tolerance
            = relative difference between values of function
    T func(T *) - function which is minimized
  
  Output:
    p, y - set to new points all within ftol of a minimum function value
    nfunk - number of function evaluations
  
  Return:
    true - if OK
    false - NMAX exceeded
*/ 

template <class T>
bool amoeba(
  T **p, T *y, 
  const int & n, 
  const T& ftol, 
  T func(T *), 
  int &nfunk,
  void print(int, T * v, const T &y) = 0){

  // constants
  const int NMAX = 10000;
  
  const T TINY = std::numeric_limits<T>::epsilon()/100;
  
  bool ok = true;
  
  int i, ihi, ilo, inhi, j;
  
  T rtol, sum, swap, ysave, ytry,
    *ptry = new T [n],
    *psum = new T [n];

  nfunk = 0;
  
  //GET_PSUM
  for (j = 0; j < n; ++j) {
    sum = 0;
    for (i = 0; i <= n; ++i) sum += p[i][j];
    psum[j] = sum;
  }
  
  do {
    ilo = 0;
    
    /* First we must determine which point is the highest (worst), next-highest, and lowest (best), by looping over the points in the simplex.*/
    if (y[0] > y[1]){
      ihi = 0;
      inhi = 1;
    } else {
      ihi = 1;
      inhi = 0;
    }
    
    for (i = 0; i <= n; ++i) {
      
      if (y[i] <= y[ilo]) ilo = i;
      
      if (y[i] > y[ihi]) {
        inhi = ihi;
        ihi = i;
      } else if (y[i] > y[inhi] && i != ihi) inhi=i;
    }
    
    /* Compute the fractional range from highest to lowest and return if satisfactory.*/
    rtol = 2*std::abs(y[ihi] - y[ilo])/
          (std::abs(y[ihi]) + std::abs(y[ilo]) + TINY);
    
    if (rtol < ftol) {    // If returning, put best point and value in slot 0.
      //SWAP(y[0],y[ilo]);
      swap = y[0];
      y[0] = y[ilo];
      y[ilo] = swap;
          
      for (i = 0; i< n; ++i)  {
        //SWAP(p[1][i],p[ilo][i])
        swap = p[0][i];
        p[0][i] = p[ilo][i];
        p[ilo][i] = swap;
      }
      break;
    }
    
    if (nfunk >= NMAX) {
      ok = false;
      break;
    }

    /*
      Printing vector with the lowest found value
    */ 
    
    if (print) print(nfunk, p[ilo], y[ilo]);    
    
    nfunk += 2;
    
    /* Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across from the high point, i.e., reflect the simplex from the high point. */

    ytry = amotry(p, y, psum, n, func, ihi, T(-1), ptry);
    
    if (ytry <= y[ilo])
      /* Gives a result better than the best point, so try an additional extrapolation by a factor 2.*/
      ytry = amotry(p, y, psum, n, func, ihi, T(2), ptry);
      
    else if (ytry >= y[inhi]) {
      /* The reflected point is worse than the second-highest, so look for an 
      intermediate lower point, i.e., do a one-dimensional contraction. */

      ysave = y[ihi];
      ytry = amotry(p, y, psum, n, func, ihi, T(0.5), ptry);
      
      if (ytry >= ysave) {
        
        for (i = 0; i <= n; ++i) {
          if (i != ilo) {
            for (j = 0; j < n; ++j) 
              p[i][j] = psum[j] = (p[i][j] + p[ilo][j])/2;
             
            y[i] = func(psum);
          }
        }
        
        nfunk += n;  //Keep track of function evaluations.
        
        //GET_PSUM
        for (j = 0; j < n; ++j) {
          sum = 0;
          for (i = 0; i <= n; ++i) sum += p[i][j];
          psum[j] = sum;
        }
      }
    } else --nfunk;
        
  } while (1);


  delete [] ptry;
  delete [] psum;
  
  return ok;
}


#endif
