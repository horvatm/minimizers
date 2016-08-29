#if !defined(__praxis_h)
#define __praxis_h

/*
  Praxis's method for optimization
  
  Using NLOPT version by S. G. Johnson and FORTRAN90 version by John Burkardt.
  
  The minfit was checked using EISPACK.
  
  Ref: 
   * Richard Brent, Algorithms for Minimization Without Derivatives, Dover, 2002
   
  Author: Martin Horvat, June 2013
*/

#include <cmath>
#include <cstdlib>
#include <string>
#include <limits>

#include "matrix.h"
#include "MersenneTwister.h"

// for using BRENTs MINFIT to do SVD
#define BRENT_MINFIT

template <typename T> 
struct Tq {
  T *v,                 // size n x n 
    *q0, *q1, *t_flin,  // size n 
    qa, qb, qc, 
    qd0, qd1, qf1,

    fbest, *xbest,      // size n

    fx, ldt, dmin;
    
  int nf, nl, max_nf;
  
  T machep, small, vsmall,
    large, vlarge,  m2, m4;
};

/*
  Calculates x^n for integer x and n for n>=0
 
  Input:
    x
    n
  Return
    value
*/ 
template <class T> T pow_ii(T x, int n){
  
  T p = 1;
  
  while (n > 0) {
    if (n & 1) {
      --n; 
      p *= x;
    } else {
      n >>= 1; 
      x *= x; 
    }
  }
  
  return p;
}
/*
  Calculates sqrt(x^2+y^2)
  
  Input:
    x
    y
  
  Return: 
    value
*/
template <class T> T hypot_(const T & x, const T & y){
  
  T a = std::abs(x), 
    b = std::abs(y);
  
  if (a < b){
    T t = a;
    a = b;
    b = t;
  } 
  
  if (a == T(0)) return 0;
  
  b /= a;
  
  return a*std::sqrt(T(1) + b*b);
}

/*
  Copy array
*/ 

template <class T> void cpy(T * dest, T *src, const int & n){
  for (int i = 0; i < n; ++i) dest[i] = src[i];
}

/* 
  The function of one real variable l that is minimized by the subroutine minny.
    
  Input:
    n - number of variables.
    j -  indicates the kind of search. 
         If J is nonzero, then the search is linear in the direction of V(*,J).
         If J is zero, then the search is parabolic, based on X, Q0 and Q1.
    l -  parameter determining the particular point at which F is to be evaluated.
         For a linear search ( j >=0 ), l is the size of the step along the 
         direction V(*,J). For a quadratic search ( j < 0), l is a parameter 
         which specifies  a point in the plane of X, Q0 and Q1.   
   f  - name of the function to be minimized.
   x[n] - base point of the search.
   
   global - struct with global variables
   q  - variable containing 
        v[n][n] - a matrix whose columns constitute search directions.
        
        v - stored in column first format
        
        If j !=0, then a linear search is being carried out in the direction of 
        v(*,j).
        
        q0[n],q1[n], two auxiliary points used to  determine the plane when a 
        quadratic search is performed.
  
  
  Output:
    nf -  the function evaluation counter
    ret - return status
          0 - OK
          1 - too many evaluations
   Return:  
    value of the function at the given point.
*/
template <class T>
T flin_(
  const int &n, 
  const int &j, 
  const T &l, 
  T f(T*), 
  T *x, 
  Tq<T> &q,
  int &ret) {
  
  // System generated locals
  T ret_val;

  // Local variables 
  int i;
  T *t = q.t_flin;
  
  // Function Body 
  if (j >= 0) {         //  THE SEARCH IS LINEAR
    
    for (i = 0; i < n; ++i) t[i] = x[i] + l * q.v[i + j * n];
    
  } else {              // THE SEARCH IS ALONG A PARABOLIC SPACE CURVE

    q.qa = l * (l - q.qd1) / (q.qd0 * (q.qd0 + q.qd1));
    q.qb = (l + q.qd0) * (q.qd1 - l) / (q.qd0 * q.qd1);
    q.qc = l * (l + q.qd0) / (q.qd1 * (q.qd0 + q.qd1));

    for (i = 0; i < n; ++i)
      t[i] = q.qa * q.q0[i] + q.qb * x[i] + q.qc * q.q1[i];
  }
  
  ret_val = f(t);

  // Save the best approximation of the minimum
  if (ret_val < q.fbest) {
    q.fbest = ret_val;
    cpy(q.xbest, t, n);
  }
  
  // The function evaluation counter nf is incremented
  // If to many evalutions
  if (++q.nf > q.max_nf) 
    ret = 1;
  else 
    ret = 0;
  
  return ret_val;
}
/* 
  The subroutine min minimizes f from x in the direction v(*,j) unless 
  j is less than 1, when a quadratic search is made in the plane 
  defined by q0,q1,x.
    
  Input:
   n - number of variables
   j - indicates the kind of search.
       If j is positive, then the search is linear in the direction of v(*,j).
       If j is zero, then the search is parabolic, based on x, q0 and q1.
   nits - maximum number of times the interval may be halved to retry 
          calculation
   d2 - is either zero or an approximation to half f
   x1 - estimate of the distance from x to the minimum along v(*,j) (or, if j=0, a curve).
   f1 -
   fk - if fk == true, then on input f1 contains the value FLIN(X1).
        otherwise x1 and f1 are ignored on entry unless final fx 
        is greater than f1
   f  - name of the function to be minimized.
   x  - base point of the search
   t  -
   h  - 
   global - struct with global variable
   q - struct containing
      v[n*n] - matrix whose columns are direction vectors along which 
               function may be minimized.
      q0[n], q1[N]
      
  Output:
    d2 -
    x1 - distance between X and the minimizer that was found.
    f1 -
    x - update base point for linear but not parabolic search
    global -
    q -
*/

template <class T>
int min_(
  const int &n, 
  const int &j, 
  const int &nits, 
  T & d2, 
  T & x1, 
  T & f1, 
  const bool &fk, 
  T f(T*), 
  T *x, 
  const T &t, 
  const T &h, 
  Tq<T> &q) {
    
  // Local variables
  bool ok, dz;

  int i, k, ret;

  T s, d1, f0, f2, t2, x2, fm,
    xm, sf1, sx1, temp;
  
  // Function Body

  sf1 = f1;
  sx1 = x1;
  k = 0;
  xm = 0;

  fm = f0 = q.fx;

  dz = (d2 < q.machep);

  // find the step size
  
  s = 0;
  for (i = 0; i < n; ++i) s += x[i]*x[i];
  s = std::sqrt(s);

  temp = (dz ? q.dmin : d2);
    
  t2 = q.m4 * std::sqrt(std::abs(q.fx)/temp + s * q.ldt) + q.m2 * q.ldt;       
  s = q.m4 * s + t;
  
  if (dz && t2 > s) t2 = s;

  if (t2 < q.small) t2 = q.small;
  temp = h/100;
  if (t2 > temp) t2 = temp;

  if (fk && f1 <= fm){
   xm = x1;
   fm = f1;  
  }

  if (!fk || std::abs(x1) < t2) {
    
    x1 = (x1 >= 0 ? t2 : -t2);
    
    f1 = flin_(n, j, x1, f, x, q, ret);
    if (ret != 0) return ret;
  }

  if (f1 <= fm) {
    xm = x1;
    fm = f1;
  }
  
  ok = true;
  do {
  
    if (dz){
      
      // Evaluate flin at another point and estimate the second derivative
      if (f0 < f1)    
        x2 = -x1;
      else
        x2 = 2*x1;
      
        
      f2 = flin_(n, j, x2, f, x, q, ret);
      if (ret != 0) return ret;
      
      if (f2 <= fm){
        xm = x2;
        fm = f2;
      }

      d2 = (x2 * (f1 - f0) - x1 * (f2 - f0)) / (x1 * x2 * (x1 - x2));
    }
    
    // Estimate the first derivative at 0
    d1 = (f1 - f0) / x1 - x1 * d2;
    dz = true;
    
    // Predict the minimum
    if (d2 <= q.small) 
      x2 = (d1 < 0 ? h : -h);
    else
      x2 = -d1/(2*d2);
    
    if (h < std::abs(x2))  x2 = (x2 > 0 ? h : -h);
      
  
    // evaluate f at the predicted minimum

    do {
    
      f2 = flin_(n, j, x2, f, x, q, ret);
      
      if (ret != 0) return ret;
      
      if (nits <= k || f2 <= f0) {
        ok = false;
        break;
      }
      
      // no success, so try again
      ++k;
      
      if (f0 < f1 && x1 * x2 > 0) break;

      x2 /= 2;
     
    } while (1);
    

  } while (ok);


  // increment the one-dimensional search counter

  ++q.nl;

  if (f2 > fm)
    x2 = xm;
  else
    fm = f2;

  // get a new estimate of the second derivative

  if (std::abs(x2 *(x2 - x1)) > q.small) {
    
     d2 = (x2 * (f1 - f0) - x1 * (fm - f0)) / (x1 * x2 * (x1 - x2));
   
  } else {
    if (k > 0) d2 = 0;
  }
  
  if (d2 <= q.small) d2 = q.small;

  x1 = x2;
  q.fx = fm;

  if (sf1 < q.fx) {
    q.fx = sf1;
    x1 = sx1;
  }

  // update x for linear but not parabolic search

  if (j >= 0) 
    for (i = 0; i < n; ++i)
      x[i] += x1 * q.v[i + j * n];
  
  return 0;
}

/* 
  Looks for the minimum of f along a curve defined by q0,q1 and x.

  Input:
    n - number of variables
    f - name of the function to be minimized
    x[n] - base point of the curve along which we search for minimum
    t - absolute tolerance 
    h - step size
    q - containing global variables
          v(N,N) -  matrix of search directions.
          q0(N), q1(N) - two auxiliary points used to define a curve through x
        
  Output:
    q      -
    x[n]   - updated base point
*/
template <class T>
void quad_(
  const int &n, 
  T f(T*),
  T *x, 
  const T &t, 
  const T &h,
  Tq<T> &q){
             
  // Local variables 
  int i, nits;
  T l, s, value;
  
  // Function Body
  s = q.fx;
  q.fx = q.qf1;
  q.qf1 = s;
  
  q.qd1 = 0;
  for (i = 0; i < n; ++i) {
    
    s = x[i];
    x[i] = l = q.q1[i];
    q.q1[i] = s;
        
    s -= l;
    q.qd1 += s * s;
  }

  l = q.qd1 = std::sqrt(q.qd1);
  
  s = 0;

  if (q.qd0 <= 0 || q.qd1 <= 0 || q.nl < 3*n*n) {

    q.fx = q.qf1;
    q.qb = q.qa = 0;
    q.qc = 1;

  } else {
    
    nits = 2;
    value = q.qf1;
    min_(n, -1, nits, s, l, value, true, f, x, t, h, q);
    
    q.qa = l * (l - q.qd1) / (q.qd0 * (q.qd0 + q.qd1));
    q.qb = (l + q.qd0) * (q.qd1 - l) / (q.qd0 * q.qd1);
    q.qc = l * (l + q.qd0) / (q.qd1 * (q.qd0 + q.qd1));
  }
  
  q.qd0 = q.qd1;
  for (i = 0; i < n; ++i) {
    s = q.q0[i];
    q.q0[i] = x[i];
    x[i] = q.qa * s + q.qb * x[i] + q.qc * q.q1[i];
  }
}
  
/* 
  Sorts the elements of d(n) into descending order and moves the 
  corresponding columns of v(n,n).  
  m is the row dimension of v as declared in the calling program.
  
  Input:
    m - the row dimension of v, which must be at least n
    n - the length of d, and the order of v
    d[n] - the vector to be sorted.  
    v[n*n] - matrix nxn stored in column first notation, to be adjusted 
             as d is sorted.  In particular, if the value that was in d(I) 
             on input is moved to d(J) on output, then the input column v(*,I) 
             is moved to the output column v(*,J).
  Output:
    d[n] - the entries of d are in descending order
    v[n*n] - sorted columns of the matrix according to d
*/

template <class T> void sort_(const int &n, T *d, T *v) {
    
  // Local variables 
  
  int i, j, k;
  
  T s;

  // Function Body
  if (n == 1) return;
    
  for (i = 0; i < n-1; ++i) {
    k = i;
    s = d[i];
    
    for (j = i+1; j < n; ++j) {
      
      if (d[j] > s) {     
        k = j;
        s = d[j];
      }
    }
    
    if (k > i) {
      
      d[k] = d[i];
      d[i] = s;
  
      for (j = 0; j < n; ++j) {
        s = v[j + i * n];
        v[j + i * n] = v[j + k * n];
        v[j + k * n] = s;
      }
    }
  }
}
/* 
  An improved version of minfit (see Golub and Reinsch, 1969)
  restricted to m=n,p=0 !!! 
  
  The singular values of the array ab are returned in q and ab is 
  overwritten with the orthogonal matrix v such that u.diag(q) = ab.v,
  where u is another orthogonal matrix.


  Input:
    m - leading dimension of AB, which must be at least N.
    n - order of the matrix AB
    tol - a tolerance which determines when a vector (a column or part of 
          a column of the matrix) may be considered "essentially" equal to zero.
    
    ab[m*n] - an N by N matrix in column first format whose singular value 
              decomposition is desired. 
    e[n] - working area
  Output:
    q[n] - singular values
    ab[m*n] - orthogonal matrix v in column first format
  
*/

#if defined(BRENT_MINFIT)
#warning "Using BRENT'S MINFIT"

template <class T>
void minfit_(
  const int &n,  
  const T& tol, 
  T *ab, 
  T *q,
  T *e) {
  
  // Local variables
  int i, j, k, l, kt;
 
  // Constants
  const int kt_max = 30;
  
  T c, f, g, h,
    s, x, y, z,
    eps, temp;
  
  bool ok;
  
  // Function Body
  
  // Householder's reduction to bidiagonal form.
  
  if (n == 1) {
    q[0] = ab[0];
    ab[0] = 1;
    return;
  }
  
  eps = std::numeric_limits<T>::epsilon();
  g = x = 0;
    
  for (i = 0; i < n; ++i) {
    
    e[i] = g;
    l = i + 1;
  
    s = 0;
    for (j = i; j < n; ++j) {
      temp = ab[j + i * n];
      s += temp * temp;
    }
  
    if (s < tol) {
      
      g = 0;

    } else {
  
      f = ab[i + i * n];

      g = (f < 0 ? std::sqrt(s) : -std::sqrt(s));
        
      h = f * g - s;
      ab[i + i * n] = f - g;
      
      for (j = l; j < n; ++j) {
        
        f = 0;
        
        for (k = i; k < n; ++k) f += ab[k + i * n] * ab[k + j * n];
      
        f /= h;
        
        for (k = i; k < n; ++k) ab[k + j * n] += f * ab[k + i * n];
        
      }
    }

    q[i] = g;
       
    if (i != n - 1) {
      
      s = 0;
      for (j = l; j < n; ++j) {
        temp = ab[i + j * n];
        s += temp * temp;
      }
      
      if ( s < tol) {
        
        g = 0;
      
      } else  {
    
        f = ab[i + l*n];

        g = (f < 0 ? std::sqrt(s) : -std::sqrt(s));
        
        h = f * g - s;
        
        ab[i + l * n] = f - g;
            
        for (j = l; j < n; ++j) e[j] = ab[i + j * n] / h;
             
        for (j = l; j < n; ++j) {
          s = 0;
          for (k = l; k < n; ++k) s += ab[j + k * n] * ab[i + k * n];
          for (k = l; k < n; ++k) ab[j + k * n] += s * e[k];
        }
      }
    }

    y = std::abs(q[i]) + std::abs(e[i]);
  
    if (y > x) x = y;
  }
  
  // Accumulation of right-hand transformations
  l = n - 1;
  ab[l + l * n] = 1;
  g = e[l];
   
  for (i = n-2; i >= 0; --i) {
    
    
    if (g != 0) {
      
      h = ab[i + (i + 1) * n] * g;
    
      for (j = l; j < n; ++j) ab[j + i * n] = ab[i + j * n] / h;
    
      for (j = l; j < n; ++j) {
        
        s = 0;
        for (k = l; k < n; ++k) s += ab[i + k * n] * ab[k + j * n];
        
        for (k = l; k < n; ++k) ab[k + j * n] += s * ab[k + i * n];
      }
    }
    
    for (j = l; j < n; ++j) ab[i + j * n] = ab[j + i * n] = 0;
    
    l = i;
    ab[l + l * n] = 1;
    g = e[l];
  }
  
      
  // Diagonalization of the bidiagonal form

  eps *= x;

  for (k = n-1; k >= 0; --k) {
    
    kt = 0;
        
    do {
      
      ++kt;
      
      if (kt > kt_max){
        e[k] = 0;
        std::cerr 
          << "minfit: QR  algorithm failed to converge in k=" << k << "\n";
      }
      
      ok = true;
      
      for (l = k; l >= 0; --l) {
          
          if (std::abs(e[l]) <= eps) {
            ok = false;
            break;
          }
          
          if (l> 0 && std::abs(q[l - 1]) <= eps) break;
      }
      
      // Cancellation of e(l) if l>1
      if (ok) {
        
        c = 0;
        s = 1;
        
        for (i = l; i <= k; ++i) {
          
          f = s * e[i];
          e[i] *= c;
          
          if (std::abs(f) <= eps) break;
          
          g = q[i];
          
          // q(i) = h = sqrt(g*g + f*f)
          q[i] = h = hypot_(f, g);
          
          if (h == 0) g = h = 1;
          
          c = g / h;
          s = -f / h;
        }
      }
      
      // Test for convergence
      z = q[k];
      
      if (l == k) {
        
        if (z < 0) {
          
          // Convergence:  q(k) is made non-negative
          q[k] = -z;
          for (j = 0; j < n; ++j) ab[j + k * n] = -ab[j + k * n];
        }
        
        break;
      }
    
      // Shift from bottom 2*2 minor
      x = q[l];
      y = q[k - 1];
      g = e[k - 1];
      h = e[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (h * 2 * y);
      
      g = hypot_(f, T(1));
      
      
      temp = (f < 0 ? f - g : f + g);
      
      f = ((x - z) * (x + z) + h * (y / temp - h)) / x;
      
      // Next QR transformation
      c = s = 1;
      
      for (i = l+1; i <= k; ++i) {
        
        g = e[i];
        y = q[i];
        h = s * g;
        g *= c;
        
        e[i - 1] = z = hypot_(f, h);
        
        if (z == 0) f = z = 1;
        
        c = f / z;
        s = h / z;
        f =  x * c + g * s;
        g = -x * s + g * c;
        h = y * s;
        y *= c;
        
        for (j = 0; j < n; ++j) {
          x = ab[j + (i - 1) * n];
          z = ab[j + i * n];
          
          ab[j + (i - 1) * n] = x * c + z * s;
          ab[j + i * n] = -x * s + z * c;
        }
        
        z = hypot_(f, h);
          
        q[i - 1] = z;
        if (z == 0) f = z = 1;
        
        c = f / z;
        s = h / z;
               
        f = c * g + s * y;
        x = -s * g + c * y;
      }

      e[l] = 0;
      e[k] = f;
      q[k] = x;
            
    } while (1);
  }
} 

#elif defined(LAPACK_MINFIT)
#warning "Using LAPACK MINFIT"

extern "C" {
  int dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
    double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, 
    double *work, int *lwork, int *info);
}
  
void minfit_(
  int n,  
  double tol, 
  double *ab, 
  double *q,
  double *e) {

  char jobu = 'N',
       jobvt= 'A';
  
  int i, j,
      lwork,
      info;
  
  double *vt = new double [n*n],
         dlwork[1];
        
  lwork = -1;
  dgesvd_(&jobu, &jobvt, &n, &n, ab, &n, q, 0, &n, vt, &n, dlwork, &lwork, &info);
  
  lwork = int(*dlwork);
  double *work = new double [lwork];
  
  dgesvd_(&jobu, &jobvt, &n, &n, ab, &n, q, 0, &n, vt, &n, work, &lwork, &info);
 
  if (info != 0) {
    if (info < 0) 
      std::cerr << " The " << -info << "-th argument had an illegal value.\n";
    else
      std::cerr << " DBDSQR did not converge info=" <<  info << "\n";
  }
  
  for (i = 0; i < n; ++i)  
    for (j = 0; j < n; ++j)
      ab[i + j*n] = vt[j + i*n];
  
  delete [] work;
  delete [] vt;
}

#elif defined(NR_MINFIT)
#warning "Using NR MINFIT"

template<class T>
inline T SIGN(const T &a, const T &b)
        {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

void minfit_(
  int n,  
  double tol, 
  double *ab, 
  double *w,
  double *e) {
    
	bool flag;
  
	int i, its, j, jj, k, l, nm, m = n;
  
	double eps,anorm, c, f, g, h, s, scale, x, y, z;
	
  double *rv1 = e,
         **u = matrix<double>(n,n),
         **v = matrix<double>(n,n); 
    
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j) u[i][j] = ab[i + n*j];
  
  eps = std::numeric_limits<double>::epsilon();
	g = scale = anorm = 0.0;
  
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += std::abs(u[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}
				f=u[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
					f=s/h;
					for (k=i;k<m;k++) u[k][j] += f*u[k][i];
				}
				for (k=i;k<m;k++) u[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += abs(u[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}
				f=u[i][l-1];
				g = -SIGN(std::sqrt(s),f);
				h=f*g-s;
				u[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=u[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += u[j][k]*u[i][k];
					for (k=l-1;k<n;k++) u[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) u[i][k] *= scale;
			}
		}
		anorm=std::max(anorm,(std::abs(w[i])+ std::abs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(u[i][j]/u[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=std::min(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) u[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
				f=(s/u[i][i])*g;
				for (k=i;k<m;k++) u[k][j] += f*u[k][i];
			}
			for (j=i;j<m;j++) u[j][i] *= g;
		} else for (j=i;j<m;j++) u[j][i]=0.0;
		++u[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (l == 0 || std::abs(rv1[l]) <= eps*anorm) {
					flag=false;
					break;
				}
				if (std::abs(w[nm]) <= eps*anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (std::abs(f) <= eps*anorm) break;
					g=w[i];
					h=hypot_(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u[j][nm];
						z=u[j][i];
						u[j][nm]=y*c+z*s;
						u[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) {
        std::cerr << "SVD failed\n";
        exit(1);
      }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=hypot_(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=hypot_(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=hypot_(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
  
  for (i = 0; i< n; ++i)
    for (j = 0; j < n; ++j) ab[i + n*j] = v[i][j];
  
  free_matrix(v);
  free_matrix(u);
}


#endif

/*
  Praxis returns the minimum of the function f(n,x) of n variables
  using the principal axis method.  the gradient of the function is
  not required.

  For a description of the algorithm, see chapter seven of
  "algorithms for finding zeros and extrema of functions without
  calculating derivatives" by R. P. Brent.

  The approximating quadratic form is
          q(x') = f(n,x) + (1/2) * (x'-x)-transpose * a * (x'-x)
  where x is the best estimate of the minimum and a is
          inverse(v-transpose) * d * inverse(v)
  (v(*,*) is the matrix of search directions; d(*) is the array
  of second differences).  If f has continuous second derivatives
  near x0, a will tend to the hessian of f at x0 as x approaches x0.

  It is assumed that on floating-point underflow the result is set
  to zero. The user should observe the comment on heuristic numbers after
  the initialization of machine dependent numbers.

  Input:
    seed - seed for the random generator 
    t0  - tolerance.  praxis attempts to return praxis=f(x) such that if x0 is 
          the true local minimum near x, then
              norm(x-x0) < t0 + sqrt(machep)*norm(x).
    h0 - maximum step size.  h0 should be set to about the maximum distance 
         from the initial guess to the minimum. if h0 is set too large or too 
         small, the initial rate of convergence may be slow.
    n - number of variables upon which the function depends; at least 2
    
    x[n] - containing initial guess of the point of minimum
    
    T func(T *) - function which is minimized
    
  Output:
   x[n] - estimated point of minimum.
   minf - estimated value of function in x
   nfunk - nr. of functions evaluations
  
  Return:
    0 - OK
    1 - too many evaluation
    2 - out of many memory    
*/

template <class T>
int praxis (
  const unsigned &seed, 
  const T &t0,  
  const T &h0, 
  const int &n, 
  T *x, 
  T f(T*), 
  T &minf,
  int &nfunk,
  void print(std::string, int, T *, const T &) = 0 
  ) {

  // Global variables
  Tq<T> q;

  int ret, i, j, k, k2, kl, 
      kt, ktm, nits;
  
  bool illc;
  
  // Local variables 
  T *d, *y, *z, *e_minfit, *work, 
    h, s, t, f1,t2, df, dn,
    sf, sl, dni, lds, scbd,
    ldfac, value, temp;
    
  q.machep = std::numeric_limits<T>::epsilon();
  q.small = q.machep * q.machep;
  
  q.vsmall = q.small * q.small;
  q.large = 1 / q.small;
  q.vlarge = 1 / q.vsmall;
  q.m2 = std::sqrt(q.machep);
  q.m4 = std::sqrt(q.m2);

  work = new (std::nothrow) T [n*n + n*8];
  if (!work) return 2;  // OUT_OF_MEMORY

  // Function Body
  q.v = work;
  q.q0 = q.v + n*n;
  q.q1 = q.q0 + n;
  q.t_flin = q.q1 + n;
  q.xbest = q.t_flin + n;
  d = q.xbest + n;
  y = d + n;
  z = y + n;
  e_minfit = z + n;
  
  /*
    Heuristic numbers:
    If the axes may be badly scaled (which is to be avoided if possible), 
    then set scbd=10, otherwise set scbd=1.
   
     If the problem is known to be ill-conditioned, set illc=false
     otherwise set illc=true.
     
     ktm is the number of iterations without improvement before the algorithm 
     terminates. ktm=4 is very cautious; usually ktm=1 is satisfactory.
  */
  scbd = 1;
  illc = false;
  ktm = 1;

  if (illc) 
    ldfac = .1;
  else 
    ldfac = .01;
  
  kt = 0;
  q.nl = 0;
  q.nf = 1;
  q.fbest = q.fx = f(x);

  // maximal number of evaluation
  q.max_nf = 15000;
  
  // initial return status
  ret = 0; // OK
  
  cpy(q.xbest, x, n);

  // prepare random generator
  MTRand rnd(seed);
    
  q.qf1 = q.fx;

  t2 = t = q.small + std::abs(t0);
  q.dmin = q.small;

  h = h0;
  temp = t *100; 
  if (h < temp) h = temp;

  q.ldt = h;

  // The first set of search directions v is the identity matrix
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) q.v[i + j * n] = 0;
    q.v[i + i * n] = 1;
  }
  
  d[0] = q.qd0 = 0;
  
  for (i = 0; i < n; ++i) 
    q.q0[i] = q.q1[i] = x[i];

  
  if (print) print("praxis", q.nf, x, q.fx);
    
  // The main loop starts here
  do {

    sf = d[0];
    d[0] = s = 0;

    // Minimize along the first direction v(*,1). 
    // fx must be passed to min by value. 
    nits = 2;
    value = q.fx;
    ret = min_(n, 0, nits, d[0], s, value, false, f, x, t, h, q);

    if (ret != 0) goto done;

    if (s <= 0) for (i = 0; i < n; ++i) q.v[i] = -q.v[i];
    
    if (sf <= d[0]*.9 || d[0] <= sf*.9) for (i = 1; i < n; ++i) d[i] = 0;

    // The inner loop starts here
    
    for (k = 1; k < n; ++k) {
    
      cpy(y, x, n);

      sf = q.fx;
      
      if (kt > 0) illc = true; 
      
      do {
        
        kl = k;
        df = 0;

        // A random step follows (to avoid resolution valleys).
        // praxis assumes that random returns a random number uniformly
        // distributed in (0,1).

        if (illc) {

          for (i = 0; i < n; ++i) {
            
            temp = rnd() - 0.5;
           
            z[i] = s = (q.ldt*T(0.1) + t2 * pow_ii(T(10), kt))* temp;
           
            for (j = 0; j < n; ++j) x[j] += s * q.v[j + i * n];
          }
          
          q.fx = f(x);
          ++q.nf;
        }
        
        // Minimize along the "non-conjugate" directions v(*,k),...,v(*,n)

        for (k2 = k; k2 < n; ++k2) {
          
          s = 0;
          nits = 2;
          value = sl = q.fx;
          
          ret = min_(n, k2, nits, d[k2], s, value, false, f, x, t, h, q);
                    
          if (ret != 0) goto done;
          
          if (illc) {
            temp = s + z[k2];
            s = d[k2] * (temp * temp);
          } else 
            s = sl - q.fx;
          
          if (df <= s){
            df = s;
            kl = k2;
          }
        }
        
        // If there was not much improvement on the first try, set
        // illc=true and start the inner loop again
        
        if (illc || df > std::abs(100*q.machep*q.fx)) break;
        
        illc = true;
      } while (1);

      // Minimize along the "conjugate" directions v(*,1),...,v(*,k-1) 
      
      for (k2 = 0; k2 < k; ++k2) {
        
        s = 0;
        nits = 2;
        value = q.fx;
        
        ret = min_(n, k2, nits, d[k2], s, value, false, f, x, t, h, q);
             
        if (ret != 0) goto done;
      }
      
      f1 = q.fx;
      q.fx = sf;
      lds = 0;
      
      for (i = 0; i < n; ++i) {
        sl = x[i];
        x[i] = y[i];
        sl -= y[i];

        y[i] = sl;
        lds += sl * sl;
      }
      
      lds = std::sqrt(lds);

      // Discard direction v(*,kl).
      // If no random step was taken, v(*,kl) is the "non-conjugate"
      // direction along which the greatest improvement was made.

      if (q.small < lds) {
      
    
        for (i = kl - 1; i>= k; --i) {
          
          for (j = 0; j < n; ++j) q.v[j + (i + 1) * n] = q.v[j + i * n];
            
          d[i+1] = d[i];
        }

        d[k] = 0;
        
        for (i = 0; i < n; ++i) q.v[i + k * n] = y[i] / lds;


        // Minimize along the new "conjugate" direction v(*,k), which is 
        // the normalized vector:  (new x) - (old x)
        
        nits = 4;
        value = f1;
        
        ret = min_(n, k, nits, d[k], lds, value, true, f, x, t, h, q);
                 
        if (ret != 0) goto done;
      
        if (lds <= 0) {

          lds = -lds;
          for (i = 0; i < n; ++i) q.v[i + k * n] = -q.v[i + k * n];
        }
      }
      
      q.ldt *= ldfac;
      if (q.ldt < lds) q.ldt = lds;
    
    
      if (print) print("praxis", q.nf, q.xbest, q.fbest);
    
      t2 = 0;
      for (i = 0; i < n; ++i) t2 += x[i]*x[i];
      t2 = q.m2 * std::sqrt(t2) + t;

      // See whether the length of the step taken since starting the
      // inner loop exceeds half the tolerance 
      
      if (t2 < 2*q.ldt)
        kt = 0;
       else 
        ++kt;
      
      if (kt > ktm) {
        ret = 0;     
        goto done;
      }
      
      // The inner loop ends here.
    }
    

    // Try quadratic extrapolation in case we are in a curved valley.

    quad_(n, f, x, t, h, q);
      
    dn = 0;
    for (i = 0; i < n; ++i) {
      d[i] = 1/std::sqrt(d[i]);
      
      if (dn < d[i]) dn = d[i];  
    }
      
    for (j = 0; j < n; ++j) {
      s = d[j] / dn;
        
      for (i = 0; i < n; ++i) q.v[i + j * n] *= s;
    }

    
    // Scale the axes to try to reduce the condition number
    
    if (scbd > 1) {
      
      s = q.vlarge;
      
      for (i = 0; i < n; ++i) {
      
        sl = 0;
        for (j = 0; j < n; ++j) {
          temp = q.v[i + j * n];
          sl += temp*temp;
        }
    
        z[i] = std::sqrt(sl);
      
        if (z[i] < q.m4) z[i] = q.m4;
    
        if (s > z[i]) s = z[i];
      }
      
      for (i = 0; i < n; ++i) {
        
        sl = s / z[i];
        z[i] = 1 / sl;
      
        if (z[i] > scbd){
          sl = 1 / scbd;
          z[i] = scbd;
        }
        
        for (j = 0; j < n; ++j) q.v[i + j * n] *= sl;
      }
    }

    // Calculate a new set of orthogonal directions before repeating the 
    // main loop.
    
    // Transpose v for minfit

    for (i = 1; i < n; ++i) {
      for (j = 0; j < i; ++j) {
        s = q.v[i + j * n];
        q.v[i + j * n] = q.v[j + i * n];
        q.v[j + i * n] = s;
      }
    }

    // Call minfit to find the singular value decomposition of v.
    // this gives the principal values and principal directions of the
    // approximating quadratic form without squaring the condition 
    // number

    minfit_(n, q.vsmall, q.v, d, e_minfit);

    // Unscale the axes

    if (scbd > 1) {
            
      for (i = 0; i < n; ++i) {
        s = z[i];
          
        for (j = 0; j < n; ++j) q.v[i + j * n] *= s;
      }

      for (i = 0; i < n; ++i) {
        
        s = 0;
        
        for (j = 0; j < n; ++j) {
          temp = q.v[j + i * n];
          s += temp * temp;
        }
        
        s = std::sqrt(s);
        
        d[i] *= s;
                
        for (j = 0; j < n; ++j) q.v[j + i * n] /= s;
      }
    }
    
    for (i = 0; i < n; ++i) {
        
      dni = dn * d[i];
      
      if (dni > q.large)
        d[i] = q.vsmall;
      else if (dni < q.small) 
        d[i] = q.vlarge;
      else 
        d[i] = 1 / (dni * dni);
    }

    // Sort the eigenvalues and eigenvectors

    sort_(n, d, q.v);
    
    // Determine the smallest eigenvalue
    
    q.dmin = d[n-1];
    if (q.dmin < q.small) q.dmin = q.small;
    
    //  The ratio of the smallest to largest eigenvalue determines whether
    // the system is ill conditioned.
    
    illc = (q.dmin < q.m2 * d[0]);
      
    // the main loop ends here
  } while (1);

done:

  if (print) print("praxis", q.nf, q.xbest, q.fbest);

  if (ret != 2) {
    minf = q.fbest;
    cpy(x, q.xbest, n);
  }
  nfunk = q.nf;
    
  delete [] work;
  return ret;
}

#endif
