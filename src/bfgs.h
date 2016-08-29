#if !defined(__bfgs_h)
#define __bfgs_h
/*
  Library describing BFGS - Broyden–Fletcher–Goldfarb–Shanno algorithm

  Using source code given in NRC++ 3rd ed.
  
  Ref: http://en.wikipedia.org/wiki/BFGS
  
  Author: Martin Horvat, November 2013
*/

#include <iostream>
#include <cmath>
#include <limits>

#include "matrix.h"

/*
  Approximative line search of a function given a starting point and a direction, a new point is calculated where the function has decreased "sufficiently".

  The line search algorithms are based on the Newton method of approximating root values of monotone decreasing functions. When the derivative f' of the function f at a starting point  x on the X-axis can be calculated, the intersection x' of the tangent at point f(x) (with gradient f'(x)) with the X-axis can be used to get a better approximation of the minima at x_{min}.

  This function is based on this idea: Given a nonlinear function f, a n-dimensional starting point x_{old} and a direction p (known as Newton direction), a new point x_{new} is calculated as

    x_{new} = x_{old} + lambda p,    0 < lambda <=1 

  in a way that f(x_{new}) has decreased sufficiently. Sufficiently means that

    f(x_{new}) <= f(x_{old}) + alpha nabla f (x_{new} - x_{old})

  where  0 < alpha < 1 is a fraction of the initial rate of decrease \nabla f p.

  Input: 
    xold[n] - starting point
    n - dimension of the domain
    fold -  value of function func at point xold
    g[n] - gradient of function func at point xold
    p[n] - point to specify the search direction (the Newton direction)
    stepmax  - maximal step size
    T func(T *) - function which is minimized
    nfunk - number of function evaluation 
    max_nfunk - maximal nfunk  
    
  Output:
    p[n] - direction
    x[n] - new position
    f  - returned value
    nfunk - number of function evaluation after leaving line search

  Return:
   -1 - failed to converge
   0 - no error
   1 - roundoff problem in lnsrch
   2 - number of function evaluations exceeded
  
  Ref:  NRC++ 3rd ed., p. 478
*/
template <class T>
int lnsrch (
  T *xold, 
  const int &n, 
  const T fold, 
  T *g, 
  T *p,
	T *x, 
  T &f, 
  const T stpmax, 
  T func(T *),
  int & nfunk,
  const int &max_nfunk ) {
    
	const T ALF = 1.0e-4, 
          TOLX = std::numeric_limits<T>::epsilon();
          
  // ALF -- ensures sufficient decrease in function value
  // TOLX  -- convergence criterion on delta x
      
	int i;
  
	T a, alam, alam2 = 0, alamin, b, disc, f2 = 0,
	  rhs1, rhs2, slope, sum, temp, test, tmplam;

	sum = 0;
	for (i = 0; i < n; ++i) sum += p[i]*p[i];
	sum = std::sqrt(sum);
  
	if (sum > stpmax) {
		temp = stpmax/sum;
    for (i = 0; i < n; ++i) p[i] *= temp;
  }
  
	slope = 0;
	for (i = 0; i < n; ++i) slope += g[i]*p[i];
  
	if (slope >= T(0)) return 1;
  
	test = 0;
	for (i = 0; i < n; ++i) {
		temp = std::abs(p[i])/std::max(std::abs(xold[i]), T(1));
		if (temp > test) test = temp;
	}
  
	alamin = TOLX/test;
	alam = 1;
  
	do {
    
		for (i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
		
    f = func(x);
    
		if (++nfunk > max_nfunk) return 2;
    
    if (alam < alamin) {
		
    	for (i = 0; i < n; ++i) x[i] = xold[i];
  
			return 0;
		
    } else if (f <= fold + ALF*alam*slope) {
    
      return -1;
    
    } else {
      
			if (alam == T(1)) {
			
      	tmplam = -slope/(2*(f - fold - slope));
			
      } else {
			
      	rhs1 = f - fold - alam*slope;
				rhs2 = f2 - fold - alam2*slope;
        
				a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
				b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))/(alam - alam2);
				
        if (a == T(0)) {
      
          tmplam = -slope/(2*b);
			
      	} else {
			
      		disc = b*b - 3*a*slope;
				
        	if (disc < T(0)) {
            tmplam = alam/2;
				  } else if (b <= T(0)) {
            tmplam = (-b + std::sqrt(disc))/(3*a);
          } else {
            tmplam = -slope/(b + std::sqrt(disc));
          }
				}
				
        if (2*tmplam > alam) tmplam = alam/2;
			}
		}
		
    alam2 = alam;
		f2 = f;
		alam = std::max(tmplam, T(0.1)*alam);
  
  } while (true);
}


/*
  Broyden–Fletcher–Goldfarb–Shanno algorithm optimization of function 
    
    f(x)    x in T^n
    
  Input:
    p[n] - initial starting point
    n - dimension of the domain
    gtol - convergence requirement on zeroing the gradient
    stpmx -  scaled maximum step length allowed in line searches
    T func(T *) - function which is minimized: f(x) = func(x)
    void dfunc(T *, T*) - derivative of the function: dfunc(x,df) df = nabla f(x)
    void print(int, int, T *, const T &) - printing intermediate results
    
  Output: 
    p[n] - final point
    fret - value of the optimum
    iter - number of iterations 
    nfunk[2] - total number of evaluations of functions and of derivative
  
  Return:
    -1 - no error: zeroing gradient
    0 - no error: converged in x
    1 - maximal number of iterations exceeded 
    2 - roundoff problem in lnsrch  
  
  Requirements: routine lnsrch.
  
  Ref:  NRC++ 3rd ed. p. 521.

*/ 

template <class T>
int bfgs(
  T* p, const int& n, const T & gtol, const T & stpmx,
  T &fret, T func(T *), void dfunc(T *, T*),
  int &iter, int nfunk[2],
  void print(int, int[2], T *, const T &) = 0) {

  // constants
  const int ITMAX = 1000,   // maximal number of iterations
            NMAX = 100000;  // maximal number of evaluations

	const T EPS = std::numeric_limits<T>::epsilon(),
          TOLX = 4*EPS;     // convergence criterion on  x

  int i, j, err, ok = 1;
  
	T den, fac, fad, fae, fp, stpmax, sum = 0, 
    sumdg, sumxi, temp, test,
    
    *dg   = new T [5*n],
    *g    = dg + n,
    *hdg  = dg + 2*n,
    *pnew = dg + 3*n,
    *xi   = dg + 4*n,
     
    **hess = matrix <T> (n, n);

  nfunk[0] = nfunk[1] = 1;
  
  fp = func(p);
  
	dfunc(p, g);
	  
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) hess[i][j] = 0;
		hess[i][i] = 1;
    
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
  
	//stpmax = stpmx*std::max(std::sqrt(sum), T(n));  // NR way of doing it
	stpmax = stpmx;
  
  for (iter = 0; iter < ITMAX; ++iter) {
    
    err = lnsrch(p, n, fp, g, xi, pnew, fret, stpmax, func, nfunk[0], NMAX);
    
    if (err > 0) {
      std::cerr << "Lnsrch: err=" << err << "\n";
      ok = 2;
      break;
    }
    
    if (print) print(iter, nfunk, pnew, fret);
    
		fp = fret;
		for (i = 0; i < n; ++i) {
			xi[i] = pnew[i] - p[i];
			p[i] = pnew[i];
		}
		
    test = 0;
		for (i = 0; i < n; ++i) {
			temp = std::abs(xi[i])/std::max(std::abs(p[i]), T(1));
			if (temp > test) test = temp;
		}
		
    if (test < TOLX) {
      ok = 0; 
      break;
    }
  		
		for (i = 0; i < n; ++i) dg[i] = g[i];
    
    ++nfunk[1];
		dfunc(p, g);
    
		test = 0;
		den = std::max(fret, T(1));
		for (i = 0; i < n; ++i) {
			temp = std::abs(g[i])*std::max(std::abs(p[i]), T(1))/den;
			if (temp > test) test = temp;
		}
    
		if (test < gtol) {
      ok = -1; 
      break;
    }
			
		for (i = 0; i < n; ++i) dg[i] = g[i] - dg[i];
    
		for (i = 0; i < n; ++i) {
			temp = 0;
			for (j = 0; j < n; ++j) temp += hess[i][j]*dg[j];
      hdg[i] = temp; 
		}
    
		fac = fae = sumdg = sumxi = 0;
    
		for (i = 0; i < n; ++i) {
      
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			
      sumdg += dg[i]*dg[i];
			sumxi += xi[i]*xi[i];
		}
    
		if (fac > std::sqrt(EPS*sumdg*sumxi)) {
      
			fac = T(1)/fac;
			fad = T(1)/fae;
      
			for (i = 0; i < n; ++i) dg[i] = fac*xi[i] - fad*hdg[i];
      
			for (i = 0; i < n; ++i) {
				for (j = i; j < n; ++j) {
					hess[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					hess[j][i] = hess[i][j];
				}
			}
		}
    
		for (i = 0; i < n; ++i) {
			temp = 0;
			for (j = 0; j < n; ++j) temp -= hess[i][j]*g[j];
      xi[i] = temp;
		}
	}
  
  delete [] dg;
     
  free_matrix(hess);
  
  return ok;
}

#endif
