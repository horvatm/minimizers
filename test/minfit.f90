function r8_hypot ( x, y )

!*****************************************************************************80
!
!! R8_HYPOT returns the value of sqrt ( X^2 + Y^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_hypot
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( abs ( x ) < abs ( y ) ) then
    a = abs ( y )
    b = abs ( x )
  else
    a = abs ( x )
    b = abs ( y )
  end if
!
!  A contains the larger value.
!
  if ( a == 0.0D+00 ) then
    c = 0.0D+00
  else
    c = a * sqrt ( 1.0D+00 + ( b / a )**2 )
  end if

  r8_hypot = c

  return
end

subroutine minfit_fort ( m, n, tol, ab, q )

!*****************************************************************************80
!
!! MINFIT computes the singular value decomposition of an N by N array.
!
!  Discussion:
!
!    This is an improved version of the EISPACK routine MINFIT
!    restricted to the case M = N and P = 0.
!
!    The singular values of the array AB are returned in Q.  AB is
!    overwritten with the orthogonal matrix V such that u.diag(q) = ab.v,
!    where U is another orthogonal matrix.
!
!    Thanks to Andreas Zuend for pointing out a potential for overflow in
!    the computation z = sqrt ( f*f + 1 ), 22 March 2012.
!
!  Modified:
!
!    22 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973,
!    Reprinted by Dover, 2002.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, Jack Dongarra, Burton Garbow, Y Ikebe, 
!    Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the leading dimension of AB, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix AB.
!
!    Input, real ( kind = 8 ) TOL, a tolerance which determines when a vector
!    (a column or part of a column of the matrix) may be considered
!    "essentially" equal to zero.
!
!    Input/output, real ( kind = 8 ) AB(M,N).  On input, an N by N array whose
!    singular value decomposition is desired.  On output, ?
!
!    Input, real ( kind = 8 ) Q(N), ?
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) ab(m,n)
  real ( kind = 8 ) c
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kt
  integer ( kind = 4 ), parameter :: kt_max = 30
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) r8_hypot
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Householder's reduction to bidiagonal form.
!
  if ( n == 1 ) then
    q(1) = ab(1,1)
    ab(1,1) = 1.0D+00
    return
  end if

  eps = epsilon ( eps )
  g = 0.0D+00
  x = 0.0D+00

  do i = 1, n

    e(i) = g
    l = i + 1

    s = sum ( ab(i:n,i)**2 )

    g = 0.0D+00

    if ( tol <= s ) then

      f = ab(i,i)

      g = sqrt ( s )
      if ( 0.0D+00 <= f ) then
        g = -g
      end if

      h = f * g - s
      ab(i,i) = f - g

      do j = l, n

        f = dot_product ( ab(i:n,i), ab(i:n,j) ) / h

        ab(i:n,j) = ab(i:n,j) + f * ab(i:n,i)

      end do 

    end if

    q(i) = g

    s = sum ( ab(i,l:n)**2 )

    g = 0.0D+00

    if ( tol <= s ) then

      if ( i /= n ) then
        f = ab(i,i+1)
      end if

      g = sqrt ( s )
      if ( 0.0D+00 <= f ) then
        g = - g
      end if

      h = f * g - s

      if ( i /= n ) then

        ab(i,i+1) = f - g
        e(l:n) = ab(i,l:n) / h

        do j = l, n

          s = dot_product ( ab(j,l:n), ab(i,l:n) )

          ab(j,l:n) = ab(j,l:n) + s * e(l:n)

        end do

      end if

    end if
    
    y = abs ( q(i) ) + abs ( e(i) )
    
    x = max ( x, y )
    
  end do
  
!
!  Accumulation of right-hand transformations.
!
  ab(n,n) = 1.0D+00
  g = e(n)
  l = n

  do ii = 2, n

    i = n - ii + 1

    if ( g /= 0.0D+00 ) then

      h = ab(i,i+1) * g

      ab(l:n,i) = ab(i,l:n) / h

      do j = l, n

        s = dot_product ( ab(i,l:n), ab(l:n,j) )

        ab(l:n,j) = ab(l:n,j) + s * ab(l:n,i)

      end do

    end if

    ab(i,l:n) = 0.0D+00
    ab(l:n,i) = 0.0D+00
    ab(i,i) = 1.0D+00

    g = e(i)

    l = i
  end do

!
!  Diagonalization of the bidiagonal form.
!
  eps = eps * x

  do kk = 1, n

    k = n - kk + 1
    kt = 0

10  continue

    kt = kt + 1

    if ( kt_max < kt ) then
      e(k) = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINFIT - Warning!'
      write ( *, '(a)' ) '  The QR algorithm failed to converge.'
    end if

    do ll = 1, k

      l = k - ll + 1
      
      if ( abs ( e(l) ) <= eps ) then
        go to 20
      end if

      if ( l /= 1  .and. abs ( q(l-1) ) <= eps ) then
        exit
      end if

    end do
!
!  Cancellation of E(L) if 1 < L.
!
    c = 0.0D+00
    s = 1.0D+00

    do i = l, k

      f = s * e(i)
      e(i) = c * e(i)
      if ( abs ( f ) <= eps ) then
        exit
      end if
      g = q(i)
!
!  q(i) = h = sqrt(g*g + f*f).
!
      h = r8_hypot ( f, g )

      q(i) = h

      if ( h == 0.0D+00 ) then
        g = 1.0D+00
        h = 1.0D+00
      end if

      c = g / h
      s = - f / h

    end do
!
!  Test for convergence.
!
20  continue

    z = q(k)

    if ( l == k ) then
      if ( z < 0.0D+00 ) then
        q(k) = - z
        ab(1:n,k) = - ab(1:n,k)
      end if
      cycle
    end if
!
!  Shift from bottom 2*2 minor.
!
    x = q(l)
    y = q(k-1)
    g = e(k-1)
    h = e(k)
    f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0D+00 * h * y )

    g = r8_hypot ( f, 1.0D+00 )

    if ( f < 0.0D+00 ) then
      temp = f - g
    else
      temp = f + g
    end if

    f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
!
!  Next QR transformation.
!
    c = 1.0D+00
    s = 1.0D+00

    do i = l + 1, k

      g = e(i)
      y = q(i)
      h = s * g
      g = g * c

      z = r8_hypot ( f, h )

      e(i-1) = z

      if ( z == 0.0D+00 ) then
        f = 1.0D+00
        z = 1.0D+00
      end if

      c = f / z
      s = h / z
      f =   x * c + g * s
      g = - x * s + g * c
      h = y * s
      y = y * c

      do j = 1, n
        x = ab(j,i-1)
        z = ab(j,i)
        ab(j,i-1) = x * c + z * s
        ab(j,i) = - x * s + z * c
      end do

      z = r8_hypot ( f, h )

      q(i-1) = z

      if ( z == 0.0D+00 ) then
        f = 1.0D+00
        z = 1.0D+00
      end if

      c = f / z
      s = h / z
      f = c * g + s * y
      x = - s * g + c * y

    end do

    e(l) = 0.0D+00
    e(k) = f
    q(k) = x
    go to 10

  end do

  return
end
