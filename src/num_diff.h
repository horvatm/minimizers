#if !defined(__num_diff_h)
#define __num_diff_h

/*
  Library for numerical differentiation.


  Derivative in j direction:
  
    partial_j f'(x) = 
        1/h_j sum_{i=1}^n a_i f(x + b_i*h_j e_j) + A h_j^p f^{(p+1)}(c) 
    
  direction:
    
    x in T^d
    
    e_j = ( delta_{kj} )_{k=1}^d
    
  Ref:
    * http://en.wikipedia.org/wiki/Numerical_differentiation
    * http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
    * http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/lanczos-low-noise-differentiators

  Author: Martin Horvat, November 2013
*/

#include <cmath>
#include <vector>
#include <limits>

//#define DEBUG 

enum Tnum_diff {
  central_3 = 0,
  central_5,
  central_7,
  central_9
};

template <class T>
struct Tnum_diff_scheme {
  int n, p;
  T A;
  std::vector<T> a, b;
};

/*
  Finite difference coefficient:
  
  Ref:
    * http://en.wikipedia.org/wiki/Finite_difference_coefficients
  
    * Fornberg, Bengt (1988), "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", Mathematics of Computation 51 (184): 699–706, doi:10.1090/S0025-5718-1988-0935077-0, ISSN 0025-5718
    
    * http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
*/
Tnum_diff_scheme<double> const num_diff_schemes_d [] = {
  
  // approximated via 3 point central differences, error: O(h^2)
  // Tnum_diff type = central_3
  { 2,
    2,
    -1./6,
    {-1./2, 1./2}, 
    {-1, 1} 
  },
  // approximated via 5 point central differences, error: O(h^4)
  // Tnum_diff type = central_5
  { 4, 
    4,
    +1./30,
    {1./12, -8./12,  8./12, -1./12},  
    {-2, -1, 1,  2}
  },
  // approximated via 7 point central differences, error: O(h^6)
  // Tnum_diff type = central_7
  { 6, 
    6,
    -1./140,
    {-1./60, +9./60, -45./60, +45./60, -9./60, +1./60}, 
    {-3, -2,  -1,   1,  2,  3}
  },
  // approximated via 9 point central differences, error: O(h^8)
  // Tnum_diff type = central_9
  { 8,
    8,
    1./630,
    { 3./840, -32./840, 168./840, -672./840, 672./840, -168./840, 32./840, -3./840},
    {-4,  -3,  -2,   -1,   1,   2,  3,  4}
  }
};

void num_diff_set_scheme(const Tnum_diff &type, int &n, double *&a, double *&b){
  n = num_diff_schemes_d[type].n;
  a = const_cast<double *> (num_diff_schemes_d[type].a.data());
  b = const_cast<double *> (num_diff_schemes_d[type].b.data()); 
}

void num_diff_set_scheme(const Tnum_diff &type, int &n, int &p, double &A, double *&a, double *&b){
  n = num_diff_schemes_d[type].n;
  p = num_diff_schemes_d[type].p;
  A = num_diff_schemes_d[type].A;
  a = const_cast<double *> (num_diff_schemes_d[type].a.data());
  b = const_cast<double *> (num_diff_schemes_d[type].b.data());
}

Tnum_diff_scheme<long double> const num_diff_schemes_l [] = {
  
  // approximated via 3 point central differences, error: O(h^2)
  // Tnum_diff type = central_3
  { 2, 
    2,
    -1.L/6,
    {-1.L/2, 1.L/2}, 
    {-1L, 1L} 
  },
  // approximated via 5 point central differences, error: O(h^4)
  // Tnum_diff type = central_5
  { 4, 
    4,
    1.L/30,
    {1.L/12, -8.L/12,  8.L/12, -1.L/12},  
    {-2L, -1L, 1L,  2L}
  },
  // approximated via 7 point central differences, error: O(h^6)
  // Tnum_diff type = central_7
  { 6,
    6,
    -1.L/140,
    {-1.L/60, +9.L/60, -45.L/60, +45.L/60, -9.L/60, +1.L/60}, 
    {-3L, -2L,  -1L,   1L,  2L,  3L}
  },
  // approximated via 9 point central differences, error: O(h^8)
  // Tnum_diff type = central_9
  { 8,
    8,
    1.L/630,
    { 3.L/840, -32.L/840, 168.L/840, -672.L/840, 672.L/840, -168.L/840, 32.L/840, -3.L/840},
    {-4L,  -3L,  -2L,   -1L,   1L,   2L,  3L,  4L}
  }
};

void num_diff_set_scheme(const Tnum_diff &type, int &n, long double *&a, long double *&b){
  n = num_diff_schemes_l[type].n;
  a = const_cast<long double *> (num_diff_schemes_l[type].a.data());
  b = const_cast<long double *> (num_diff_schemes_l[type].b.data()); 
}

void num_diff_set_scheme(const Tnum_diff &type, int &n, int &p, long double &A, long double *&a, long double *&b){
  n = num_diff_schemes_l[type].n;
  p = num_diff_schemes_l[type].p;
  A = num_diff_schemes_l[type].A;
  a = const_cast<long double *> (num_diff_schemes_l[type].a.data());
  b = const_cast<long double *> (num_diff_schemes_l[type].b.data());
}

#if defined(_QUAD)
Tnum_diff_scheme<quad> const num_diff_schemes_q [] = {
  
  // approximated via 3 point central differences, error: O(h^2)
  // Tnum_diff type = central_3
  { 2,
    2,
    -1.Q/6,
    {-1.Q/2, 1.Q/2}, 
    {-1.Q, 1.Q} 
  },
  // approximated via 5 point central differences, error: O(h^4)
  // Tnum_diff type = central_5
  { 4,
    4,
    1.Q/30,
    {1.Q/12, -8.Q/12,  8.Q/12, -1.Q/12},  
    {-2.Q, -1.Q, 1.Q,  2.Q}
  },
  // approximated via 7 point central differences, error: O(h^6)
  // Tnum_diff type = central_7
  { 6,
    6,
    -1.Q/140,
    {-1.Q/60, +9.Q/60, -45.Q/60, +45.Q/60, -9.Q/60, +1.Q/60}, 
    {-3.Q, -2.Q,  -1.Q,   1.Q,  2.Q,  3.Q}
  },
  // approximated via 9 point central differences, error: O(h^8)
  // Tnum_diff type = central_9
  { 8,
    8,
    1.Q/630,
    { 3.Q/840, -32.Q/840, 168.Q/840, -672.Q/840, 672.Q/840, -168.Q/840, 32.Q/840, -3.Q/840},
    {-4.Q,  -3.Q,  -2.Q,   -1.Q,   1.Q,   2.Q,  3.Q,  4.Q}
  }
};

void num_diff_set_scheme(const Tnum_diff &type, int &n, quad *&a, quad *&b){
  n = num_diff_schemes_q[type].n;
  a = const_cast<quad *> (num_diff_schemes_q[type].a.data());
  b = const_cast<quad *> (num_diff_schemes_q[type].b.data()); 
}

void num_diff_set_scheme(const Tnum_diff &type, int &n, int &p, quad &A, quad *&a, quad *&b){
  n = num_diff_schemes_q[type].n;
  p = num_diff_schemes_q[type].p;
  A = num_diff_schemes_q[type].A;
  a = const_cast<quad *> (num_diff_schemes_q[type].a.data());
  b = const_cast<quad *> (num_diff_schemes_q[type].b.data());
}

#endif


/*
  Routine for numerical differentiation. Serial code.
  
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
   df[d] - derivative
*/

template <class T>
void num_diff(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x), 
  T *df){

  int n;
  
  T *a, *b;
  
  num_diff_set_scheme(type, n, a, b);
  
  int i, j, k;
  
  T *y = new T [d], sum;

  for (i = 0; i < d; ++i){
    
    sum = 0;
    
    for (k = 0; k < n; ++k) {

      for (j = 0; j < d; ++j) y[j] = x[j];
  
      y[i] += b[k]*h[i];

      sum += a[k]*func(y);
      
    }

    df[i] = sum/h[i];
  }
  
  delete [] y;
}
/*
  Routine for numerical differentiation. Parallelized first by dimension and 
  then by function calls. Nested OpenMP.
  
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
   df[d] - derivative
*/
template <class T>
void num_diff_omp(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x), 
  T *df){

  int n;
  
  T *a, *b;
  
  num_diff_set_scheme(type, n, a, b);
  
  #pragma omp parallel for shared(x, a, b, df)  
  for (int i = 0; i < d; ++i){
    
    T sum = 0;
    
    #if defined(DEBUG)
      std::cerr << "num_diff::i=" << i << std::endl;
      #if defined(_OPENMP)
      std::cerr << "num_diff::OpenMP L1 id=" << omp_get_thread_num() << std::endl;
      #endif
    #endif
    
    #pragma omp parallel for shared(x, a, b) reduction(+:sum) 
    for (int k = 0; k < n; ++k) {

      #if defined(DEBUG)
        std::cerr << "num_diff::i,k=" << i << "," << k << std::endl;
        #if defined(_OPENMP)
        std::cerr << "num_diff::OpenMP L2 id=" << omp_get_thread_num() << std::endl;
        #endif
      #endif
            
      T *y = new T [d];
        
      for (int j = 0; j < d; ++j) y[j] = x[j];
  
      y[i] += b[k]*h[i];

      sum += a[k]*func(y);
      
      delete [] y;
    }

    df[i] = sum/h[i];
  }
}

/*
  Routine for numerical differentiation. Openmp parallelization of only 
  function evaluation without nesting.
  
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
   df[d] - derivative
*/

template <class T>
void num_diff_omp_f(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x), 
  T *df){

  int n;
  
  T *a, *b;
  
  num_diff_set_scheme(type, n, a, b);
  
  int i, j, k, len = n*d;
  
  T *f = new T [len],
    *y = new T [len*d],
    *p = y;
  
  for (i = 0; i < d; ++i)
    for (j = 0; j < n; ++j) {
      for (k = 0; k < d; ++k) p[k] = x[k];
      p[i] += b[j]*h[i];
      p += d;
    }
    
  #pragma omp parallel for shared(y,f) private(i) 
  for (i = 0; i < len; ++i) f[i] = func(y + i*d);
  
  T sum;
  p = f;
  for (i = 0; i < d; ++i) {
    sum = 0;
    for (j = 0; j < n; ++j) sum += a[j]*p[j];
    df[i] = sum/h[i];
    p += n;
  }
  
  delete [] f;
  delete [] y;
}


/*
  Routine for numerical differentiation. Using Kahan summation to reduce 
  rounding errors.
  
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
   df[d] - derivative
*/

template <class T>
void num_diff_kahan(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x), 
  T *df){

  int n;
  
  T *a, *b;
  
  num_diff_set_scheme(type, n, a, b);
  
  #pragma omp parallel for shared(x, a, b, df)  
  for (int i = 0; i < d; ++i){
    
    T sum = 0,   // sum
      c = 0,     // correction
      sum_,      // partial sum
      t,         // term of summation
      
      *y = new T [d];

    #if defined(DEBUG)
      std::cerr << "num_diff_kahan::i=" << i << std::endl;
      #if defined(_OPENMP)
      std::cerr << "num_diff_kahan::OpenMP L1 id=" << omp_get_thread_num() << std::endl; 
      #endif
    #endif
    
    for (int k = 0; k < n; ++k) {  // parallelize !!!

      for (int j = 0; j < d; ++j) y[j] = x[j];

      y[i] += b[k]*h[i];
      
      t = a[k]*func(y) - c;
    
      sum_ = sum + t;
      
      c = static_cast<T>(sum_ - sum) - t;
       
      sum = sum_;
    }

    df[i] = sum/h[i];
    
    delete [] y;
  }
}
/*
  Estimate optimal step size in all dimension. Serial code rewritten. 
  Using OpenMP nested if enabled.
  
  The upper bound of the differential error is given by

  | Df(x) - f'(x) | <= |S|/h + |A| h^p |f^{p+1}|
     ^        ^
    numerical  exact
  
  with summation error
  
    S  propto eps_0

  where  eps_0 is machine precision
   
  Meaning 
      | Df(x) - f'(x) | <= eps_0 |B| /h + |A| h^p |f^{p+1}|
        
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- initial step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
    h[d] -- optimal step sizes in each dimension
*/


template <class T>
void num_diff_estimate_step(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x)){

  int n, p;
  
  T A, *a, *b, eps = std::numeric_limits<T>::epsilon();
   
  num_diff_set_scheme(type, n, p, A, a, b);
  
  int i, j, k, r = p + 1;
  
  T dS, 
    *y = new T [d],
    *f = new T [r + 1],
    t, fmax;
   
  for (i = 0; i < d; ++i) {
    
    // Estimate rounding error of the scheme
    // Ref: Higham, 1993
    dS = 0;
    for (j = 0; j < n; ++j) {
          
      for (k = 0; k < d; ++k) y[k] = x[k];
      y[i] += b[j]*h[i];
      
      dS += std::abs(a[j]*func(y)); 
    }
    
    dS *= eps*(n - 1);  // sequential summation
    // dS *= 2*eps;     // Kahan summation
    
    //  Estimate r = p+1-th derivative using central differences:
    //
    //    f^{(r)}(x) = h^{-r} delta_h^r[f](x) + O(h^2)
    //
    //  with
    //    delta_h^r[f](x) = sum_{i=0}^r (-1)^i Binom(r,i) f(x + (r/2-i)h)
    //    
    //    delta_h[f](x) = f(x+h/2) - f(x-h/2)
    //  Ref: 
    //    http://en.wikipedia.org/wiki/Finite_difference
    for (j = 0; j <= r; ++j) {
      for (k = 0; k < d; ++k) y[k] = x[k];
      y[i] += (T(r)/2 - j)*h[i];
      
      f[j] = func(y);
    }
    
    fmax = std::abs(f[0]);
    for (j = 1; j <= r; ++j) {
      t = std::abs(f[j]);
      if (fmax < t) fmax  = t;
    }
    
    // divided differences scheme
    for (j = 0; j < r; ++j)
      for (k = 0; k < r-j; ++k) f[k] -= f[k + 1];
      
    if (std::abs(f[0]) < r*(r+1)*fmax*eps/2) {
      h[i] = std::abs(dS/(fmax*eps));              // dS/h ~ fmax*eps = abs. error of f
    } else {
      h[i] *= exp(log(std::abs(dS/(p*A*f[0])))/r); // min dS/h + A f^(r)/h^p
    }
    
  }
  
  delete [] y;
  delete [] f;  
}

/*
  Estimate optimal step size in all dimension. Serial code rewritten. 
  Using OpenMP nested if enabled.
  
  The upper bound of the differential error is given by

  | Df(x) - f'(x) | <= |S|/h + |A| h^p |f^{p+1}|
     ^        ^
    numerical  exact
  
  with summation error
  
    S  propto eps_0

  where  eps_0 is machine precision
   
  Meaning 
      | Df(x) - f'(x) | <= eps_0 |B| /h + |A| h^p |f^{p+1}|
        
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- initial step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
    h[d] -- optimal step sizes in each dimension
*/


template <class T>
void num_diff_estimate_step_omp(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x)){

  int n, p;
  
  T A, *a, *b, eps = std::numeric_limits<T>::epsilon();
   
  num_diff_set_scheme(type, n, p, A, a, b);
  
  int r = p + 1;
      
  #pragma omp parallel for shared(x, a, b, h)  
  for (int i = 0; i < d; ++i) {
    
    
    #if defined(DEBUG)
      std::cerr << "num_diff_estimate_step::i=" << i << std::endl;
      #if defined(_OPENMP)
      std::cerr << "num_diff_estimate_step::OpenMP L2 id=" << omp_get_thread_num() << std::endl;
      #endif
    #endif
    
        
    // Estimate rounding error of the scheme
    // Ref: Higham, 1993
    T dS = 0;
    #if defined(DEBUG)
      std::cerr << "num_diff_estimate_step::estimate rounding error" << std::endl;
    #endif
    
    #pragma omp parallel for shared(x, a, b, h) reduction(+:dS) 
    for (int j = 0; j < n; ++j) {
      
      #if defined(DEBUG)
        std::cerr << "num_diff_estimate_step::i,j=" << i << "," << j << std::endl;
        #if defined(_OPENMP)
        std::cerr << "num_diff_estimate_step::OpenMP L2 id=" << omp_get_thread_num() << std::endl;
        #endif
      #endif
      
      T *y = new T [d];
          
      for (int k = 0; k < d; ++k) y[k] = x[k];
      y[i] += b[j]*h[i];
      
      dS += std::abs(a[j]*func(y)); 
      
      delete [] y;
    }
    
    dS *= eps*(n - 1);  // sequential summation
    // dS *= 2*eps;     // Kahan summation
    
    //  Estimate r = p+1-th derivative using central differences:
    //
    //    f^{(r)}(x) = h^{-r} delta_h^r[f](x) + O(h^2)
    //
    //  with
    //    delta_h^r[f](x) = sum_{i=0}^r (-1)^i Binom(r,i) f(x + (r/2-i)h)
    //    
    //    delta_h[f](x) = f(x+h/2) - f(x-h/2)
    //  Ref: 
    //    http://en.wikipedia.org/wiki/Finite_difference
    T *f = new T [r + 1];

    #if defined(DEBUG)
      std::cerr << "num_diff_estimate_step::estimate f^{(p+1)}" << std::endl;
    #endif
    
    #pragma omp parallel for shared(x, a, b, f, h)
    for (int j = 0; j <= r; ++j) {
      
      
      #if defined(DEBUG)
        std::cerr << "num_diff_estimate_step::i,j=" << i << "," << j << std::endl;
        #if defined(_OPENMP)
        std::cerr << "num_diff_estimate_step::OpenMP L2 id=" << omp_get_thread_num() << std::endl;
        #endif
      #endif
      
      
      T *y = new T [d];
      
      for (int k = 0; k < d; ++k) y[k] = x[k];
      y[i] += (T(r)/2 - j)*h[i];
      
      f[j] = func(y);
      
      delete [] y;
    }
    
    T t, fmax = std::abs(f[0]);
    for (int j = 1; j <= r; ++j) {
      t = std::abs(f[j]);
      if (fmax < t) fmax  = t;
    }
    
    // divided differences scheme
    for (int j = 0; j < r; ++j)
      for (int k = 0; k < r-j; ++k) f[k] -= f[k + 1];
      
    /*
    std::cerr << "#dS="   << dS 
              << " A="    << A 
              << " f[0]=" << f[0]
              << " df="   << f[0]/exp(log(h[i])*r) << '\n';
    */
    if (std::abs(f[0]) < r*(r+1)*fmax*eps/2) {
      h[i] = std::abs(dS/(fmax*eps));              // dS/h ~ fmax*eps = abs. error of f
    } else {
      h[i] *= exp(log(std::abs(dS/(p*A*f[0])))/r); // min dS/h + A f^(r)/h^p
    }
    
    delete [] f;
  }
}

/*
  Estimate optimal step size in all dimension. Intended for use with OpenMP 
  without nesting, if enabled. 
  
  The upper bound of the differential error is given by

  | Df(x) - f'(x) | <= |S|/h + |A| h^p |f^{p+1}|
     ^        ^
    numerical  exact
  
  with summation error
  
    S  propto eps_0

  where  eps_0 is machine precision
   
  Meaning 
      | Df(x) - f'(x) | <= eps_0 |B| /h + |A| h^p |f^{p+1}|
        
  Input:
    type in Tnum_diff 
    d -- dimension of domain
    h[d] -- initial step sizes in each dimension
    x[d] -- point at which derivative is calculated
    T func(T *x) -- function which is differentiated
    
  Output:
    h[d] -- optimal step sizes in each dimension
*/


template <class T>
void num_diff_estimate_step_omp_f(
  const Tnum_diff &type, 
  const int &d,
  T *h,
  T *x, 
  T func(T *x)){

  int n, p;
  
  T A, *a, *b, eps = std::numeric_limits<T>::epsilon();
   
  num_diff_set_scheme(type, n, p, A, a, b);
  
  
  // Estimate rounding error of the scheme
  // Ref: Higham, 1993
  
  int i, j, k, len = d*n;
  
  T *y = new T [len*d],
    *q = y, tmp;
    
  for (i = 0; i < d; ++i) 
    for (j = 0; j < n; ++j) {
      for (k = 0; k < d; ++k) q[k] = x[k];
      q[i] += b[j]*h[i];
      q += d;
    }
  
  T *f = new T [len];
  #pragma omp parallel for shared(y, f) private(i)
  for (i = 0; i < len; ++i) f[i] = func(y + i*d);
  
  delete [] y;
  
  q = f;
  
  T *dS = new T [d];
    
  for (i = 0; i < d; ++i) {
    tmp = 0;
    for (j = 0; j < n; ++j) tmp += std::abs(a[j]*q[j]);
    
    dS[i] = eps*(n - 1)*tmp;  // sequential summation
    // dS[i] *= 2*eps;        // Kahan summation

    q += n;
  }
  delete [] f;
  
  //  Estimate r = p+1-th derivative using central differences:
  //
  //    f^{(r)}(x) = h^{-r} delta_h^r[f](x) + O(h^2)
  //
  //  with
  //    delta_h^r[f](x) = sum_{i=0}^r (-1)^i Binom(r,i) f(x + (r/2-i)h)
  //    
  //    delta_h[f](x) = f(x+h/2) - f(x-h/2)
  //
  //  and calculate optimal step size
  //
  //
  //  Ref: 
  //    http://en.wikipedia.org/wiki/Finite_difference

  int r = p + 1;
  len = d*(r + 1);
  
  y = new T [len*d];
  f = new T [len];

  q = y;
  for (i = 0; i < d; ++i)
    for (j = 0; j <= r; ++j) {
      for (k = 0; k < d; ++k) q[k] = x[k];
      q[i] += (T(r)/2 - j)*h[i];
      q += d;  
    }
  
  #pragma omp parallel for shared(y, f) private(i)
  for (i = 0; i < len; ++i) f[i] = func(y + i*d);
  
  delete [] y;
    
  T fmax, thresh;
  
  q = f;
  
  for (i = 0; i < d; ++i) {
    
    fmax = std::abs(q[0]);
    for (j = 1; j <= r; ++j) {
      tmp = std::abs(q[j]);
      if (fmax < tmp) fmax  = tmp;
    }
    
    // divided differences scheme
    for (j = 0; j < r; ++j)
      for (k = 0; k < r - j; ++k) 
        q[k] -= q[k + 1];
    
    thresh = fmax*r*(r+1)*eps/2;                      // threshold for computing r+1 
                                                      // central difference  
    
    if (std::abs(q[0]) < thresh)
      h[i] = std::abs(dS[i]/(fmax*eps));               // dS/h ~ fmax*eps = abs. error of f
    else
      h[i] *= exp(log(std::abs(dS[i]/(p*A*q[0])))/r);  // min dS/h + A f^(r)/h^p
    
    q += r + 1;
  }
  
  delete [] dS;
  delete [] f;
}



#endif
