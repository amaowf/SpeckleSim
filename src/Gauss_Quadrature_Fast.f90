module Gauss_Quadrature_Fast
    
    use libmath
	 use lib_sie_constants

!*****************************************************************************80
!
!! MAIN is the main program for LEGENDRE_RULE_FAST.
!
!  Discussion:!
!    This program computes a standard Gauss-Legendre quadrature rule
!    and writes it to a file.
!
!  Usage:!
!    legendre_rule_fast ( n, a, b )! 2D quadrature points
!    where!
!    * n is the number of points in the rule;
!    * a is the left endpoint;
!    * b is the right endpoint.
!
!  Licensing:!
!    This code is distributed under the GNU LGPL license. !
!  Modified:!
!    28 June 2009!
!  Author:!
!    John Burkardt
!
    implicit none
    contains

subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, 
!    which guarantee the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end

subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal. 
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
    digit = 0
 
  else

    digit = -1

  end if
 
  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH 
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if
 
  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine i4_to_s_left ( i4, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an I4 to a left-justified string.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ).
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!        I4  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = i4
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '
 
  return
end
subroutine legendre_compute_glr ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( dp ) X(N), the abscissas.
!
!    Output, real ( dp ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = dp ) p
  real ( dp ) pp
  real ( dp ) w(n)
  real ( dp ) x(n)
!
!  Get the value and derivative of the N-th Legendre polynomial at 0.0.
!
  call legendre_compute_glr0 ( n, p, pp )
!
!  If N is odd, then zero is a root.
!
  if ( mod ( n, 2 ) == 1 ) then

    x((n+1)/2) = 0.0D+00
    w((n+1)/2) = pp
!
!  If N is even, we have to compute a root.
!
  else

    call legendre_compute_glr2 ( p, n, x((n/2)+1), w((n/2)+1) )

  end if
!
!  Get the complete set of roots and derivatives.
!
  call legendre_compute_glr1 ( n, x, w )
!
!  Compute W.
!
  w(1:n) = 2.0D+00 / &
    ( 1.0D+00 - x(1:n) ) / ( 1.0D+00 + x(1:n) ) / w(1:n) / w(1:n)

  w(1:n) = 2.0D+00 * w(1:n) / sum ( w(1:n) )

  return
end

subroutine legendre_compute_glr0 ( n, p, pp )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Output, real ( dp ) P, PP, the value of the N-th Legendre polynomial
!    and its derivative at 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) k
  real ( dp ) p
  real ( dp ) pm1
  real ( dp ) pm2
  real ( dp ) pp
  real ( dp ) ppm1
  real ( dp ) ppm2
  real ( dp ) rk
!
!  Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
!
  pm2 = 0.0D+00
  pm1 = 1.0D+00
  ppm2 = 0.0D+00
  ppm1 = 0.0D+00

  do k = 0, n - 1
    rk = real ( k, dp )
    p = - rk * pm2 / ( rk + 1.0D+00 )
    pp = ( ( 2.0D+00 * rk + 1.0D+00 ) * pm1 &
                     - rk             * ppm2 ) &
         / (           rk + 1.0D+00 )
    pm2 = pm1
    pm1 = p
    ppm2 = ppm1
    ppm1 = pp
  end do

  return
end

subroutine legendre_compute_glr1 ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
!
!  Discussion:
!
!    This routine requires that a starting estimate be provided for one
!    root and its derivative.  This information will be stored in entry
!    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Input/output, real ( dp ) X(N).  On input, a starting value
!    has been set in one entry.  On output, the roots of the Legendre 
!    polynomial.
!
!    Input/output, real ( dp ) W(N).  On input, a starting value
!    has been set in one entry.  On output, the derivatives of the Legendre 
!    polynomial at the zeros.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) M, the number of terms in the Taylor expansion.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 30
  integer ( kind = 4 ) n

  real ( dp ) dk
  real ( dp ) dn
  real ( dp ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n2
  real ( dp ), parameter :: pi = 3.141592653589793D+00
  !real ( dp ) rk2_leg
  integer ( kind = 4 ) s
  !real ( dp ) ts_mult
  real ( dp ) u(m+2)
  real ( dp ) up(m+1)
  real ( dp ) w(n)
  real ( dp ) x(n)
  real ( dp ) xp

  if ( mod ( n, 2 ) == 1 ) then
    n2 = ( n - 1 ) / 2 - 1
    s = 1
  else
    n2 = n / 2 - 1
    s = 0
  end if

  dn = real ( n, dp )

  do j = n2 + 1, n - 2

    xp = x(j+1)

    h = rk2_leg ( pi/2.00, -pi/2.0, xp, n ) - xp

    u(1) = 0.0D+00
    u(2) = 0.0D+00
    u(3) = w(j+1)

    up(1) = 0.0D+00
    up(2) = u(3)

    do k = 0, m - 2

      dk = real ( k, dp )

      u(k+4) = &
      ( &
        2.0D+00 * xp * ( dk + 1.0D+00 ) * u(k+3) &
        + ( dk * ( dk + 1.0D+00 ) - dn * ( dn + 1.0D+00 ) ) * u(k+2) / ( dk + 1.0D+00 ) &
      ) / ( 1.0D+00 - xp ) / ( 1.0D+00 + xp ) / ( dk + 2.0D+00 )

      up(k+3) = ( dk + 2.0D+00 ) * u(k+4)

    end do

    do l = 0, 4
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 )
    end do

    x(j+2) = xp + h
    w(j+2) = ts_mult ( up, h, m - 1 )   

  end do

  do k = 0, n2 + s
    x(k+1) = - x(n-1-k+1)
    w(k+1) = w(n-1-k+1)
  end do

  return
end

subroutine legendre_compute_glr2 ( pn0, n, x1, d1 )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR2 finds the first real root.
!
!  Discussion:
!
!    This routine is only called if N is even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, real ( dp ) PN0, the value of the N-th Legendre polynomial
!    at 0.
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Output, real ( dp ) X1, the first real root.
!
!    Output, real ( dp ) D1, the derivative at X1.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) M, the number of terms in the Taylor expansion.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 30

  real ( dp ) d1
  real ( dp ) eps
  !integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  real ( dp ), parameter :: pi = 3.141592653589793D+00
  real ( dp ) pn0
  real ( dp ) rk
  real ( dp ) rn
  real ( dp ) scale
  real ( dp ) step
  real ( dp ) theta
  real ( dp ) u(m+1)
  real ( dp ) up(m+1)
  real ( dp ) x1
  real ( dp ) x1k(m+1)

  k = ( n + 1 ) / 2

  theta = pi * real ( 4 * k - 1, dp ) / real ( 4 * n + 2, dp )

  x1 = ( 1.0D+00 - real ( n - 1, dp ) &
    / real ( 8 * n * n * n, dp ) &
    - 1.0D+00 / real ( 384 * n * n * n * n, kind = 4 ) &
    * ( 39.0D+00 - 28.0D+00 / ( sin ( theta ) * sin ( theta ) ) ) ) &
    * cos ( theta )
!
!  Scaling.
!
  scale = 1.0D+00 / x1
!
!  Recurrence relation for Legendre polynomials.
!
  u(1:m+1) = 0.0D+00
  up(1:m+1) = 0.0D+00

  rn = real ( n, dp )

  u(1) = pn0

  do k = 0, m - 2, 2

    rk = real ( k, dp )

    u(k+3) = ( rk * ( rk + 1.0D+00 ) - rn * ( rn + 1.0D+00 ) ) * u(k+1) &
      / ( rk + 1.0D+00 ) / ( rk + 2.0D+00 ) / scale / scale

    up(k+2) = ( rk + 2.0D+00 ) * u(k+3) * scale

  end do
!
!  Flip for more accuracy in inner product calculation
!
  u = u(m+1:1:-1)
  up = up(m+1:1:-1)

  x1k(1:m+1) = 1.0D+00

  step = huge ( step )
  l = 0
!
!  Newton iteration.
!
  eps = epsilon ( eps )

  do while ( eps < abs ( step ) .and. l < 10 )
    l = l + 1
    step = dot_product ( u(1:m+1),  x1k(1:m+1) ) &
         / dot_product ( up(1:m+1), x1k(1:m+1) )
    x1 = x1 - step
    x1k(1) = 1.0D+00
    x1k(2) = scale * x1
    do kk = 3, m + 1
      x1k(kk) = x1k(kk-1) * scale * x1
    end do
    x1k(1:m+1) = x1k(m+1:1:-1)
  end do

  d1 = dot_product ( up(1:m+1), x1k(1:m+1) )

  return
end

subroutine legendre_handle ( n, a, b, x, w )

!*****************************************************************************80
!
!! LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( dp ) A, B, the left and right endpoints.
! 
  implicit none

  real ( dp ) a
  real ( dp ) b
  !integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  !character ( len = 80 ) output_r
  !character ( len = 80 ) output_w
  !character ( len = 80 ) output_x
  real ( dp ) r(2)
  real ( dp ) t1
  real ( dp ) t2
  character ( len = 10 ) tag
  real ( dp ), allocatable, dimension ( : ), intent(out) :: w
  real ( dp ), allocatable, dimension ( : ), intent(out) :: x

  r(1) = a
  r(2) = b
!
!  Compute the rule.
!
  allocate ( w(n) )
  allocate ( x(n) )

  call cpu_time ( t1 )
  call legendre_compute_glr ( n, x, w )
  call cpu_time ( t2 )
  
!  write ( *, '(a)' ) ' '
!  write ( *, '(a,g14.6,a)' ) '  Computation required ', t2 - t1, ' seconds.'
!!
!  Rescale the data.
!
  call rescale ( n, a, b, x, w )
!
!  Write the data to files.
!
  call i4_to_s_left ( n, tag )

  !output_w = 'leg_o' // trim ( tag ) // '_w.txt'
  !output_x = 'leg_o' // trim ( tag ) // '_x.txt'
  !output_r = 'leg_o' // trim ( tag ) // '_r.txt'
  !
  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) '  Weight file will be   "' // trim ( output_w ) // '".'
  !write ( *, '(a)' ) '  Abscissa file will be "' // trim ( output_x ) // '".'
  !write ( *, '(a)' ) '  Region file will be   "' // trim ( output_r ) // '".'
  !          
  !call r8mat_write ( output_w, 1, n, w )
  !call r8mat_write ( output_x, 1, n, x )
  !call r8mat_write ( output_r, 1, 2, r )
  !    
  !deallocate ( w )
  !deallocate ( x )

  return
end

subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( dp ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( dp ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end

subroutine rescale ( n, a, b, x, w )

!*****************************************************************************80
!
!! RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( dp ) A, B, the endpoints of the new interval.
!
!    Input/output, real ( dp ) X(N), on input, the abscissas for [-1,+1].
!    On output, the abscissas for [A,B].
!
!    Input/output, real ( dp ) W(N), on input, the weights for [-1,+1].
!    On output, the weights for [A,B].
!
  implicit none

  integer ( kind = 4 ) n

  real ( dp ) a
  real ( dp ) b
  real ( dp ) w(n)
  real ( dp ) x(n)

  x(1:n) = ( ( a + b ) + ( b - a ) * x(1:n) ) / 2.0D+00
  w(1:n) = ( b - a ) * w(1:n) / 2.0D+00

  return
end

function rk2_leg ( t1, t2, x, n )

!*****************************************************************************80
!
!! RK2_LEG advances the value of X(T) using a Runge-Kutta method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( dp ) T1, T2, the range of the integration interval.
!
!    Input, real ( dp ) X, the value of X at T1.
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Output, real ( dp ) RK2_LEG, the value of X at T2.
!
  implicit none

  real ( dp ) f
  real ( dp ) h
  integer ( kind = 4 ) j
  real ( dp ) k1
  real ( dp ) k2
  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) n
  real ( dp ) rk2_leg
  real ( dp ) snn1
  real ( dp ) t
  real ( dp ) t1
  real ( dp ) t2
  real ( dp ) x
  real ( dp ) x2

  x2 = x

  h = ( t2 - t1 ) / real ( m, dp )
  snn1 = sqrt ( real ( n * ( n + 1 ), dp ) )
  t = t1

  do j = 0, m - 1

    f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
    x2 = x2 + k1

    t = t + h

    f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
    x2 = x2 + 0.5D+00 * ( k2 - k1 )

  end do

  rk2_leg = x2

  return
end
!

!
function ts_mult ( u, h, n )

!*****************************************************************************80
!
!! TS_MULT...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( dp ) U(N+1), ...
!
!    Input, real ( dp ) H, ...
!
!    Input, integer ( kind = 4 ) N, ...
!
!    Output, real ( dp ) TS_MULT, ...
!
  implicit none

  integer ( kind = 4 ) n

  real ( dp ) h
  real ( dp ) hk
  integer ( kind = 4 ) k
  real ( dp ) ts
  real ( dp ) ts_mult
  real ( dp ) u(n+1)
  
  ts = 0.0D+00
  hk = 1.0D+00
  do k = 1, n
    ts = ts + u(k+1) * hk
    hk = hk * h
  end do

  ts_mult = ts
  return  
end 

end 