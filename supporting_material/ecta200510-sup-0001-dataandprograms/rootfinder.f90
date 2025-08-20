module rootfinder
	
	implicit none
	
	contains
	
		subroutine zero_rc ( a, b, t, arg, status, value )

	!*****************************************************************************80
	!
	!! ZERO_RC seeks the root of a function F(X) using reverse communication.
	!
	!  Discussion:
	!
	!    The interval [A,B] must be a change of sign interval for F.
	!    That is, F(A) and F(B) must be of opposite signs.  Then
	!    assuming that F is continuous implies the existence of at least
	!    one value C between A and B for which F(C) = 0.
	!
	!    The location of the zero is determined to within an accuracy
	!    of 6 * MACHEPS * abs ( C ) + 2 * T.
	!
	!    The routine is a revised version of the Brent zero finder 
	!    algorithm, using reverse communication.
	!
	!    Thanks to Thomas Secretin for pointing out a transcription error in the
	!    setting of the value of P, 11 February 2013.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license. 
	!
	!  Modified:
	!
	!    11 February 2013
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Reference:
	!
	!    Richard Brent,
	!    Algorithms for Minimization Without Derivatives,
	!    Dover, 2002,
	!    ISBN: 0-486-41998-3,
	!    LC: QA402.5.B74.
	!
	!  Parameters:
	!
	!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
	!    sign interval.
	!
	!    Input, real ( kind = 8 ) T, a positive error tolerance.
	!
	!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
	!    does not need to initialize this value.  On return with STATUS positive,
	!    the user is requested to evaluate the function at ARG, and return
	!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
	!    estimate for the function's zero.
	!
	!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between 
	!    the user and the routine.  The user only sets STATUS to zero on the first 
	!    call, to indicate that this is a startup call.  The routine returns STATUS
	!    positive to request that the function be evaluated at ARG, or returns
	!    STATUS as 0, to indicate that the iteration is complete and that
	!    ARG is the estimated zero
	!
	!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as requested
	!    by the routine on the previous call.
	!
	  implicit none

	  real ( kind = 8 ) a
	  real ( kind = 8 ) arg
	  real ( kind = 8 ) b
	  real ( kind = 8 ), save :: c
	  real ( kind = 8 ), save :: d
	  real ( kind = 8 ), save :: e
	  real ( kind = 8 ), save :: fa
	  real ( kind = 8 ), save :: fb
	  real ( kind = 8 ), save :: fc
	  real ( kind = 8 ) m
	  real ( kind = 8 ), save :: machep
	  real ( kind = 8 ) p
	  real ( kind = 8 ) q
	  real ( kind = 8 ) r
	  real ( kind = 8 ) s
	  real ( kind = 8 ), save :: sa
	  real ( kind = 8 ), save :: sb
	  integer ( kind = 4 ) status
	  real ( kind = 8 ) t
	  real ( kind = 8 ) tol
	  real ( kind = 8 ) value
	!
	!  Input STATUS = 0.
	!  Initialize, request F(A).
	!
	  if ( status == 0 ) then

		machep = epsilon ( a )

		sa = a
		sb = b
		e = sb - sa
		d = e

		status = 1
		arg = a
		return
	!
	!  Input STATUS = 1.
	!  Receive F(A), request F(B).
	!
	  else if ( status == 1 ) then

		fa = value

		status = 2
		arg = sb
		return
	!
	!  Input STATUS = 2
	!  Receive F(B).
	!
	  else if ( status == 2 ) then

		fb = value

		if ( 0.0D+00 < fa * fb ) then
		  status = -1
		  return
		end if

		c = sa
		fc = fa

	  else

		fb = value

		if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
			 ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
		  c = sa
		  fc = fa
		  e = sb - sa
		  d = e
		end if

	  end if
	!
	!  Compute the next point at which a function value is requested.
	!
	  if ( abs ( fc ) < abs ( fb ) ) then

		sa = sb
		sb = c
		c = sa
		fa = fb
		fb = fc
		fc = fa

	  end if

	  tol = 2.0D+00 * machep * abs ( sb ) + t
	  m = 0.5D+00 * ( c - sb )

	  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
		status = 0
		arg = sb
		return
	  end if

	  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

		e = m
		d = e

	  else

		s = fb / fa

		if ( sa == c ) then

		  p = 2.0D+00 * m * s
		  q = 1.0D+00 - s

		else

		  q = fa / fc
		  r = fb / fc
		  p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
		  q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

		end if

		if ( 0.0D+00 < p ) then
		  q = - q
		else
		  p = - p
		end if

		s = e
		e = d

		if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
		  p < abs ( 0.5D+00 * s * q ) ) then
		  d = p / q
		else
		  e = m
		  d = e
		end if

	  end if

	  sa = sb
	  fa = fb

	  if ( tol < abs ( d ) ) then
		sb = sb + d
	  else if ( 0.0D+00 < m ) then
		sb = sb + tol
	  else
		sb = sb - tol
	  end if

	  arg = sb
	  status = status + 1

	  return
	end
	
	subroutine local_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization 
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real ( kind = 8 ) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between 
!    the user and the routine.  The user only sets STATUS to zero on the first 
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: eps
  real ( kind = 8 ), save :: fu
  real ( kind = 8 ), save :: fv
  real ( kind = 8 ), save :: fw
  real ( kind = 8 ), save :: fx
  real ( kind = 8 ), save :: midpoint
  real ( kind = 8 ), save :: p
  real ( kind = 8 ), save :: q
  real ( kind = 8 ), save :: r
  integer ( kind = 4 ) status
  real ( kind = 8 ), save :: tol
  real ( kind = 8 ), save :: tol1
  real ( kind = 8 ), save :: tol2
  real ( kind = 8 ), save :: u
  real ( kind = 8 ), save :: v
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: w
  real ( kind = 8 ), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
  if ( status == 0 ) then

    if ( b <= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
      write ( *, '(a)' ) '  A < B is required, but'
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      status = -1
      stop 1
    end if

    c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0D+00

    status = 1
    arg = x

    return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
  else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
  else if ( 2 <= status ) then

    fu = value

    if ( fu <= fx ) then

      if ( x <= u ) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        a = u
      else
        b = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end if
!
!  Take the next step.
!
  midpoint = 0.5D+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0D+00
  tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
  if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
    status = 0
    return
  end if
!
!  Is golden-section necessary?
!
  if ( abs ( e ) <= tol1 ) then

    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Consider fitting a parabola.
!
  else

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )
    if ( 0.0D+00 < q ) then
      p = - p
    end if
    q = abs ( q )
    r = e
    e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
    if ( &
      ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
      ( p <= q * ( a - x ) ) .or. &
      ( q * ( b - x ) <= p ) ) then

      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Choose a parabolic interpolation step.
!
    else

      d = p / q
      u = x + d

      if ( ( u - a ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

      if ( ( b - u ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

    end if

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol1 <= abs ( d ) ) then
    u = x + d
  end if

  if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
  end if
!
!  Request value of F(U).
!
  arg = u
  status = status + 1

  return
end	
		
end module