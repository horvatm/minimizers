#if !defined(__powell_h)
#define __powell_h

/*
  Powell's method for optimization
  
  Using source code given in NRC.
  
  Ref: 
   * http://en.wikipedia.org/wiki/Powell%27s_method
   * Press et al - Numerical recipes 2nd and 3rd ed.
     Section 10.7. Direction Set (Powell's) Methods in Multidimensions 
     File: c10-7.pdf
   * Richard Brent, Algorithms for Minimization Without Derivatives, Dover, 2002
   
  Author: Martin Horvat, June 2013
*/


#include <iostream>
#include <cmath>
#include <limits>

template <class T>
class Tlinmin {
  
  int n, nfunk;
  
  T (*func)(T *);
  
  void (*print)(std::string, int, T *, const T &);
   
  T TINY,
    GOLD,
    GLIMIT,
    ITMAX,
    CGOLD,
    RTOL,
    ATOL,
    *p_,
    *xi_, 
    *xt;
  
  std::string label;
  
  // Univariate function: func from p along xi direcion
  T F(const T & x){
    
    for (int i = 0; i < n; ++i) xt[i] = p_[i] + x*xi_[i];
    T val = func(xt);
    
    ++nfunk;
    
    if (print) print(label, nfunk, xt, val);
            
    return val;
  }
   
  void SHFT(T &a, T &b, T &c, const T& d){
   a = b;
   b = c;
   c = d;
  }

  T SIGN(const T &a, const T &b){
   return ((b)  > T(0) ? std::abs(a) : -std::abs(a));
  }
   
  /* Bracket a minimum along a line. 
    
  Global in  class:
  
    GOLD - default ratio by which successive intervals are magnified
    GLIMIT - maximum magnification allowed for a parabolic-fit step.
    F(x) - function minimized
  
  Input:
    ax,bx -- distinct initial points
    
  Output:
   ax, bx, cx - bracket a minimum of the function. 
   fa, fb, and fc - function values at the three point

  
  */
  void mnbrak(T &ax, T &bx, T&cx, T &fa, T &fb, T &fc){
    
    T ulim, u, r, q, fu, dum;

    label = "mnbrak";
    
    fa = F(ax);
    fb = F(bx);
    
    /*
     Switch roles of a and b so that we can go downhill in 
     the direction from a to b.
    */

    if (fb > fa){
      SHFT(dum, ax, bx, dum);
      SHFT(dum, fb, fa, dum);
    }
    
    cx = bx + GOLD*(bx - ax);   // First guess for c.

    fc = F(cx);
    
    while (fb > fc){            // Keep returning here until we bracket.
      
      r = (bx - ax)*(fb - fc);
      q = (bx - cx)*(fb - fa);
      
      dum = q -r;
      /*
       Compute u by parabolic extrapolation from a, b, c. 
       TINY is used to prevent any possible division by zero. */

      u = bx - ((bx - cx)*q - (bx - ax)*r)/
          (2*SIGN(std::max(std::abs(dum), TINY), dum));
      
      ulim = bx + GLIMIT*(cx - bx);
      
      // We won’t go farther than this. Test various possibilities

      if ((bx - u)*(u - cx) > T(0)){  //Parabolic u is between b and c: try it.

        fu = F(u);
        
        if (fu < fc){         // Got a minimum between b and c.
          ax = bx;
          bx = u;
          fa = fb;
          fb = fu;
          return;
        } else if (fu > fb){  // Got a minimum between between a and u.
          cx = u;
          fc = fu;
          return;
        }

        // Parabolic fit was no use. Use default magnification.
        u = cx + GOLD*(cx - bx);
        fu = F(u);
     
      } else if ((cx - u)*(u - ulim) > T(0)){
        // Parabolic fit is between c and its allowed limit.
        
        fu = F(u);
        
        if (fu < fc){

          SHFT(bx, cx, u, cx + GOLD*(cx - bx));
          SHFT(fb, fc, fu, F(u));
        }
        
      } else if ((u - ulim)*(ulim - cx) >= T(0)){
        
        // Limit parabolic u to maximum allowed value.
        u = ulim;
        fu = F(u);
        
      } else{
        
        // Reject parabolic u, use default magnification.
        u = cx + GOLD*(cx - bx);
        fu = F(u);
        
      }
      // Eliminate oldest point and continue.

      SHFT(ax, bx, cx, u);
      SHFT(fa, fb, fc, fu);
    }
  }
   
  /*  Brent’s method for minimization of univariate function F(x).
   
  The values rtol and atol define a tolerance 
  
     TOL(x) = RTOL * abs (x) + ATOL.
  
  F is never evaluated at two points closer than TOL. 
  
  RTOL -- relative error tolerance
  ATOL -  positive absolute error tolerance 
  
  Def: F(x) is unimodal on the interval if and only if it is monotonic on 
        either side of the single optimal point x* in the interval.
  
  If F is a delta-unimodal function near minimum and delta < RTOL, then xmin 
  approximates the abscissa of the global minimum of F on the interval [A,B]:
    
    |argmin_x {f(x)} - xmin| < 2TOL(xmin) + delta <  3 TOL(xmin)
  
  If F is not unimodal, then xmin may approximate a local, but perhaps non-global, 
  minimum to the same accuracy. (Brent, p. 79.)
  
  It was shown in (Brent, p.64) that
  
    delta > [2 |F_0| eps/F_0'']^(1/2) (1 - eps + M delta/(6 F''_0))
    
  with eps = rel. machine precision and
     
    F_0 = F(xmin)
    F''_0 = F''(xmin)  
  
  From here follows that RTOL should not be smaller than 
  
    [2 |F_0| eps/(xmin^2 F'')]^(1/2)
  
    
  Global in class:
  
    ITMAX - maximum allowed number of iterations;
    CGOLD - square of the inverse of the golden ratio
    ATOL > 0 
    RTOL  in  [2 eps_machine, sqrt(eps_machine)]
    F(x) - function minimized
  
  Input:
    ax, bx, cx - bracketing triplet of abscissas: 
                 bx in [ax, cx], f(bx) < f(ax), f(cx)
  
  Output:
    xmin - abscissa of the minimum  
    fx - minimum function value
    
  Return:
    true  - everything OK
 */

  bool brent(const T & ax, const T & bx, const T & cx, T &xmin, T &fx){
    
    label = "brent";
    
    int iter;
    
    T a, b, d = 0, etemp,
      fu, fv, fw,
      p, q, r, tol1, tol2,
      u, v, w, x, xm,
      e = 0;  // This will be the distance moved on the step before last.

    // a and b must be in ascending order,but input abscissas need not be.
    a = ((ax < cx) ? ax : cx);    
    b = ((ax > cx) ? ax : cx);
    
    x = w = v = bx;
    
    fw = fv = fx = F(x);
    
    for (iter = 0; iter < ITMAX; ++iter){             // Main program loop.

      xm = (a + b)/2;
      tol1 = RTOL*std::abs(x) + ATOL;
      tol2 = 2*tol1;

      if (std::abs(x - xm) <= (tol2 - (b - a)/2)) {   // Test for done here.
        xmin = x;
        return true;
      }
      
      if (std::abs(e) > tol1){              // Construct a trial parabolic fit.
        
        r = (x - w)*(fx - fv);
        q = (x - v)*(fx - fw);
        p = (x - v)*q - (x - w)*r;
        
        q = 2*(q - r);
        
        if (q > T(0)) p = -p;
        
        q = std::abs(q);
        etemp = e;
        e = d;
        
        if (std::abs(p) >= std::abs(q*etemp/2) || 
            p <= q*(a - x) || 
            p >= q*(b - x)) {
              
          e = (x >= xm ? a - x : b - x);
          d = CGOLD*e;
          
      /* The above conditions determine the acceptability of the parabolic fit. 
      Here we take the golden section step into the larger of the two segments.*/

        } else {
          
          // Take the parabolic step.
          d = p / q;
          u = x + d;
          
          if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
        }
      } else {
        
        //  A golden-section step.
        e = (x >= xm ? a - x : b - x);
        d = CGOLD*e;
      }
      
      u = (std::abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));		
      
      fu = F(u);
      
      if (fu <= fx){ 
        
        if (u >= x) a = x; else b = x;
        
        SHFT(v, w, x, u);
        SHFT(fv, fw, fx, fu);
        
      } else {
    
        if (u < x) a = u; else b = u;
    
        if (fu <= fw || w == x){
          v = w;
          fv = fw;
          w = u;
          fw = fu;

        } else if (fu <= fv || v == x || v == w){
          v = u;
          fv = fu;
        }
      }
    }
    
    std::cerr << "Too many iterations in BRENT\n";
    
    xmin = x;
    
    return false;
  }

    
  public:
  
  
  /* 
  
  Initialization of Powell's method for minimization using the Brent's line 
  minimization.

  Input:
   n - dimension
   tol[2] - absolute and relative tolerance for Brent's methods
   func - function which is minimized
  
  Note: 
    tol[0] - absolute tolerance just need to be positive, NR suggests
             = std::numeric_limits<T>::epsilon()*1e-3;
    
            
    tol[1] - should not be smaller than sqrt(machine precision)
            = sqrt(std::numeric_limits<T>::epsilon())
  */
  
  Tlinmin(const int &n, const T tol[2], 
          T func(T *),
          void print(std::string, int, T *, const T &) = 0
  ): n(n), nfunk(0), func(func), print(print) {
          
    TINY = std::numeric_limits<T>::epsilon()*1e-10; // some small value
    
    // mnbrak
    GOLD = (std::sqrt(T(5)) + 1)/2;
    GLIMIT = 1000;
    
    // brent
    ATOL = tol[0];
    RTOL = tol[1];
    
    ITMAX = 1000;
    CGOLD = (3 - std::sqrt(T(5)))/2;
  
    xt = new T [n];
  }
  
  ~Tlinmin(){
    delete [] xt;
  }
  
  /*
    Finding minimum from p along xi direction.
     
    Input:
      p[n] -  n-dimensional point
      xi[n] - n-dimensional direction 
      
    Output:
      p[n] -  point where F(p) takes on a minimum along the direction xi
      xi[n] - replaced by actual vector displacement that p was moved.
      fret - value of func at the returned location p
      
    Returns:
      true - everything OK
  */ 
  
  bool get_min(T *p,  T *xi, T &fret){
   
    p_ = p, xi_ = xi;
    
    T ax = 0, xx = 1, bx = 2,   //  Initial guess for brackets.
      fa, fx, fb, xmin;
    
    mnbrak(ax, xx, bx, fa, fx, fb);
    
    bool ok = brent(ax, xx, bx, xmin, fret);
    
    for(int i = 0; i < n; ++i){     // Construct the vector results to return.
      xi[i] *= xmin;
      p[i] += xi[i];
    }
         
    return ok;
  }
  
  int get_nfunk() const {return nfunk;}
  
  int & get_nfunk() {return nfunk;}
};

/*
  Powell Minimization of a function func of n variables. 

  Input:
    p[n] - initial starting point
    xi[n][n] - initial matrix
               columns contain the initial set of directions (usually n unit vectors)
    tol[2] - absolute tolerance of position (Brent method) and 
             relative tolerance of values (Powell method)
    T func(T *) - function which is minimized
           
  Output:
    p - best point found, 
    xi - current direction set, 
    fret - returned function value at p
    iter - number of iterations taken
    nfunk - nr. of functions evaluations
  
  Return:
    true - everything OK
     
  Requirements: class Tlinmin. 
*/

template <class T>
bool powell(T *p, T **xi, int n, T tol[2], 
            T &fret, T func(T *), 
            int &iter, int &nfunk, 
            void print(std::string, int, T *, const T &) = 0) {
  
  // constants
  const int ITMAX = 1000; // maximal number of iterations
  const int NMAX = 100000; // maximal number of evaluations
  
  const T TINY = std::numeric_limits<T>::epsilon()*1e-10;   // some tiny number
    
  int i, ibig, j;
  
	T t, fptt, fp, del,
    t1, t2,
    tol_brent[2],
    ftol_powell,
    *pt = new T [n],
    *ptt = new T [n], 
    *xit = new T [n];
  
  bool ok = true;
    
  //
  // tolerances for methods
  //
  tol_brent[0] = tol[0]/3;
  tol_brent[1] = sqrt(std::numeric_limits<T>::epsilon());
    
  ftol_powell = tol[1];
  
  
  //
  // setup line minimization 
  //
  Tlinmin<T> linmin(n, tol_brent, func, print);
    
  fret = func(p);
	
  nfunk = 1;
  
  for (j = 0; j < n; ++j) pt[j] = p[j]; // Save the initial point.
    
	for (iter = 1; ; ++iter){
    
    /* 
      Will be the biggest function decrease. 
      In each iteration, loop over all directions in the set. 
    */

		fp = fret;
    
		ibig = 0;
		
    del = 0;
    
		for (i = 0; i < n; ++i){
      
      // Copy the direction.
			for(j = 0; j < n; ++j) xit[j] = xi[j][i]; 
      
			fptt = fret;
			
      // Minimize along it and record it, if it is the largest decrease so far.

      ok = linmin.get_min(p, xit, fret);
      
      if (!ok) {
        std::cerr << "powell::Line search nr 1. failed\n";
        break;
      }
            
			if (std::abs(fptt - fret) > del){
        
				del = std::abs(fptt - fret);
				
        ibig = i;
			}
		}
    
    // Print best point
    
    if (print) print("powell", linmin.get_nfunk(), p, fret);
    
    // Termination criterion.
    
		if (2*std::abs(fp - fret) <= ftol_powell*(std::abs(fp) + std::abs(fret)) + TINY) break;
     
		if (iter >= ITMAX || linmin.get_nfunk() > NMAX) {
      std::cerr << "powell::Too many iterations in routine POWELL\n";
      ok = false;
      break;
    }
    
    /* 
    Construct the extrapolated point and the average direction moved. 
    Save the old starting point. 
    */

		for (j = 0; j < n; ++j){
			ptt[j] = 2*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
    
		fptt = func(ptt); //  Function value at extrapolated point.
    ++linmin.get_nfunk();
    
		if (fptt < fp){
      t1 = fp - fret - del;
      t2 = fp - fptt;
			t = 2*(fp - 2*fret + fptt)*t1*t1 - del*t2*t2;
      
			if (t < T(0)){
        
        // Move to the minimum of the new direction and save the new direction.
				ok = linmin.get_min(p, xit, fret);
        
        if (!ok) {
           std::cerr << "powell::Line search nr. 2 failed\n";
          break;
        }
        
        for (j = 0; j < n; ++j) xi[j][ibig] = xit[j];
			}
		}
	}
  
  delete [] xit;
	delete [] ptt;
  delete [] pt;
  
  nfunk = linmin.get_nfunk();
  
  return ok;
}


#endif // #if !defined(__powell_h)

