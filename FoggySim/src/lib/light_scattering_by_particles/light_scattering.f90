  module light_scattering
      use file_io
      implicit none

      contains

!     ..................................................................
!     .  Light Scattering by Particles: Computational Methods          .
!     .  by P.W. Barber and S.C. Hill                                  .
!     .  copyright (c) 1990 by World Scientific Publishing Co Pte Ltd  .
!     .                                                                .
!     .  equation numbers in columns 73-80 are references to the text  .
!     ..................................................................
!     ................................................................
!     .  calculate the differential scattering cross section in all  .
!     .    directions by projecting the surface of a spherical       .
!     .    surface onto a rectangular coordinate system.             .
!     .    theta: the theta scattering angle, is projected onto the  .
!     .           radius r in the rectangular plane                  .
!     .           r = theta/180 degrees   (r varies from 0 to 1)     .
!     .           when r > 1, set scattering to backscatter value    .
!     .    phi: the phi scattering angle, is the azimuthal angle in  .
!     .         the spherical and rectangular coordinate systems     .
!     .                                                              .
!     .                               |                              .
!     .                      x     x  |  x      x                    .
!     .                               |                              .
!     .                      x     x  |  x      x                    .
!     .                    ------------------------ y                .
!     .                      x     x  |  x      x                    .
!     .                               |                              .
!     .                      x     x  |  x      x                    .
!     .                               |                              .
!     .                               x                              .
!     .                                                              .
!     .  start at (x,y) = (-1,-1) and continue in the + x direction  .
!     .    then increment y by dlt and continue at (-1,-1+dlt)       .
!     .                                                              .
!     .  inputs: ip = polarization parallel or perpendicular         .
!     .          x = size parameter (ka)                             .
!     .          cm = complex index of refraction, (real,imag)       .
!     .               (imag is positive for absorption)              .
!     .          npnts = number of grid points = npnts by npnts      .
!     .                                                              .
!     .  dimension of arrays f(*), g(*), amat(*) and cnrm(*):        .
!     .    nc = int(x+4.05*x**.3333+2.0), e.g., for x = 200,         .
!     .    nc = 225                                                  .
!     .                                                              .
!     .  dimension of arrays bj(*), by(*), hkl(*), and               .
!     .    pnmllg(*): nc+1, e.g., for x = 200, nc+1 = 226            .
!     .                                                              .
!     .  arrays are set for a maximum size parameter of 200          .
!     ................................................................
      subroutine s2(ip, x, cmr, cmi, npnts)
      implicit none
      ! dummy
      integer, optional :: ip, npnts
      real, optional :: x, cmr, cmi

      ! auxiliary
      complex cm,ci,cim,f(225),g(225),fth,fph
      real cnrm, pi, ac, pnmllg, snorm
      real rn, rf, bkscat, dlt, yg, xg, r
      real theta, phi, sinph, cosph, costh, p1, p2
      integer nc, nci, n, ig, jg, n1, j
      dimension pnmllg(226),cnrm(225),ac(1000)
      common /cfcom/ f,g,cnrm
      pi = 3.14159265358979
      ci = cmplx(0.0,1.0)
      write(6,100)
      write(6,*) 'enter polarization, size parameter, mr, mi, npnts'
!      read(5,*) ip,x,cmr,cmi,npnts
!      ip = 1
!      x = 2*pi * 10 ! sphere radius: a=10um;  wave length: 1um; a / wave length = 10
!      cmr = 1.33    ! water: 1.33
!      cmi = 0       ! water: 0
!      npnts = 10
!     .........................................
!     .  set the complex index of refraction  .
!     .    for an exp(-iwt) time variation    .
!     .........................................
      cm = cmplx(cmr,cmi)
      snorm = 1.0/(pi*x**2)
      open(unit=9,file='s2.dat')
      rewind 9
      call sphere(x,cm,ip,nc)
      nci = nc+1
!     ............................................................
!     .  calculate the logarithm of the differential scattering  .
!     .    cross section in the backscatter direction            .
!     .  note: calculated using the expression for parallel      .
!     .    incident polarization, but the result is the same if  .
!     .    the coefficients for perpendicular polarization are   .
!     .    used because the absolute value is taken              .
!     ............................................................
      fth = 0.0
      do 10 n = 1,nc
        cim = ci**n
        rn = real(n)
        rf = (2.0*rn+1.0)/(rn*(rn+1.0))
        fth = fth+rf*cim*(ci*f(n)+g(n))
10    continue
      bkscat = alog10(snorm*abs(fth)**2/16.)                            !eq 4.24
!     .....................................................
!     .  calculate the differential scattering cross      .
!     .    section over a rectangular grid corresponding  .
!     .    to the surface of a sphere                     .
!     .....................................................
      dlt = 2.0/real(npnts-1)
!     ....................................
!     .  set starting y grid value (yg)  .
!     .  enter loop to vary y values     .
!     ....................................
      yg = -1.0
      do 40 ig = 1,npnts
!     ....................................
!     .  set starting x grid value (xg)  .
!     .  enter loop to vary x values     .
!     ....................................
        xg = -1.0
        do 30 jg = 1,npnts
          r = sqrt(xg**2+yg**2)
          if(r.lt.1.0) then
!     .............................................................
!     .  calculate the theta and phi scattering angles for radii  .
!     .    within the unit circle on the rectangular grid         .
!     .............................................................
            theta = r*pi
            phi = 0.0
            if(abs(xg).gt.0.0.or.abs(yg).gt.0.0) phi = atan2(yg,xg)
!     .............................................
!     .  calculate the logarithm of the           .
!     .    differential scattering cross section  .
!     .    at (xg,yg) and store in array ac(*)    .
!     .............................................
            sinph = sin(phi)
            cosph = cos(phi)
            costh = cos(theta)
            call genlgp(theta,pnmllg,nci)
            fth = 0.0
            fph = 0.0
            do 20 n = 1,nc
              n1 = n+1
              cim = ci**(-n1)
              rn = real(n)
              p1 = rn*costh*pnmllg(n1)-(rn+1.0)*pnmllg(n)
              p2 = pnmllg(n1)
              if(ip.eq.1) then
!     ..........................................
!     .  calculate parallel polarization case  .
!     ..........................................
                fth = fth+cim*cosph*(p2*f(n)+ci*p1*g(n))*cnrm(n)        !eq 4.10a
                fph = fph-cim*sinph*(p1*f(n)+ci*p2*g(n))*cnrm(n)        !eq 4.10b
              else
!     ...............................................
!     .  calculate perpendicular polarization case  .
!     ...............................................
                fth = fth+cim*sinph*(-p2*f(n)+ci*p1*g(n))*cnrm(n)       !eq 4.11a
                fph = fph+cim*cosph*(-p1*f(n)+ci*p2*g(n))*cnrm(n)       !eq 4.11b
              end if
20          continue
!     ............................................................
!     .  calculate the logarithm of the differential scattering  .
!     .    cross section in the (theta,phi) direction            .
!     ............................................................
            ac(jg) = alog10(snorm*(abs(fth)**2+abs(fph)**2))
          else
!     ....................................................
!     .  set the differential scattering cross section   .
!     .    to the backscatter value for all grid points  .
!     .    on or outside the unit circle                 .
!     ....................................................
            ac(jg) = bkscat
          end if
!     .......................
!     .  increment x value  .
!     .......................
          xg = xg+dlt
30      continue
!     ......................................
!     .  write out data for all xg values  .
!     .    (for given yg value)            .
!     ......................................
        write(9,110) (ac(j),j=1,npnts)
!     .......................
!     .  increment y value  .
!     .......................
        yg = yg+dlt
40    continue
      close(unit=9)
!      stop
100   format('.....................................................',/,  &
            '.  calculate scattered intensity in all directions  .',/,  &
            '.  output is written to s2.dat                      .',/,  &
            '.....................................................',//, &
            'polarization: parallel (1) perpendicular (2)',/,           &
            'size parameter: x',/,                                      &
            'index of refraction: real,imaginary (+ for absorption)',/, &
            'npnts: number of grid points (npnts by npnts)',/)
110   format(e14.6)
!      end

      contains

      subroutine sphere(x,cm,ip,nc)
!     ..............................................................
!     .  calculate the scattered field f(n) and g(n) coefficients  .
!     .    the f(n) and g(n) for theta incident polarization are,  .
!     .    within an n-dependent factor, the same as the b(n) and  .
!     .    a(n) coefficients, respectively, defined in C.F.        .
!     .    Bohren and D.R. Huffman, Absorption and Scattering of   .
!     .    Light by Small Particles (Wiley- Interscience,New       .
!     .    York,1983), p.100                                       .
!     ..............................................................
      real x, cnrm, xc, bj, bjm
      integer ip, nc, nmx
      complex b,z,cm,ci,hkl(226),an,amat(225),f(225),g(225)
      logical rv
      common /cfcom/ f,g,cnrm
      complex test_a(255)
      complex test_b(255)
      double complex a_bh(255)
      double complex b_bh(255)
      dimension cnrm(225)
      ci = (0.0,1.0)
!     ......................................................
!     .  set the number of terms required for convergence  .
!     ......................................................
      xc = x+4.05*x**.3333+2.0                                          !eq 4.16
      nc = int(xc)
      nci = nc+1
      z = cm*x
!     ..................................................
!     .  logarithmic derivative calculation - set the  .
!     .    starting order for downward recursion       .
!     ..................................................
      nmx = int(max(xc,abs(z)))+15                                      !eq 4.20
      an = 0.0
      do 10 n = nmx,nc+1,-1
        rn = real(n)
        an = rn/z-1.0/(an+rn/z)
10    continue
      amat(nc) = an
      do 20 n = nc,2,-1
        rn = real(n)
        amat(n-1) = rn/z-1.0/(amat(n)+rn/z)                             !eq 4.19
20    continue
!     ...................................................
!     .  calculate the Bessel functions - the order is  .
!     .    incremented by one in the hkl(*) array       .
!     ...................................................
      call besh(x,hkl,nci)
      bj = real(hkl(1))
!     ................................
!     .  calculate the coefficients  .
!     ................................
      do 30 n = 1,nc
        rn = real(n)
        rf = 2.0*rn*(rn+1.0)
        bjm = bj
        bj = real(hkl(n+1))
!     .......................................................
!     .  scattering coefficients for theta                  .
!     .    (parallel) incident polarization                 .
!     .    f(n) = -ci**n*rf*(Bohren and Huffman's b(n))     .
!     .    g(n) = ci**(n+1)*rf*(Bohren and Huffman's a(n))  .
!     .......................................................
        b = cm*amat(n)+rn/x
        test_b(n) = (b*bj-bjm)/(b*hkl(n+1)-hkl(n))
        f(n) = -ci**n*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n))                 !eq 4.18a
        b_bh(n) = (b*bj-bjm)/(b*hkl(n+1)-hkl(n))
        b = amat(n)/cm+rn/x
        test_a(n) = (b*bj-bjm)/(b*hkl(n+1)-hkl(n))
        g(n) = ci**(n+1)*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n))              !eq 4.18b
        a_bh(n) = (b*bj-bjm)/(b*hkl(n+1)-hkl(n))
        if(ip.eq.2) then
!     ...............................................
!     .  scattering coefficients for phi            .
!     .    (perpendicular) incident polarization    .
!     ...............................................
          f(n) = -f(n)                                                  !eq 4.8a
          g(n) = g(n)                                                   !eq 4.8b
        end if
!        print *, "f n", n,
!     ........................................
!     .  calculate the normalization factor  .
!     .    (used in main program)            .
!     ........................................
        cnrm(n) = (2.0*rn+1.0)/(rf*rn*(rn+1.0))
30    continue

      open(unit=99, file="temp/test_get_coefficients_a_b_barberh_gt.csv", status='unknown')
      rv = write_csv(99, (/ a_bh, b_bh /), (/ "n_", "ab" /))
      close(99)

      return
      end subroutine
      subroutine besh(x,hankel,nc)
!     ...................................................
!     .  calculate Hankel functions                     .
!     .  bj = Bessel function of the first kind         .
!     .  by = Bessel function of the second kind        .
!     .  x = real argument                              .
!     .  nc = number of orders (0 to nc-1)              .
!     .  the order of the functions is incremented by   .
!     .    one in the bj(*),by(*) and hankel(*) arrays  .
!     .                                                 .
!     .  arrays are set for nc = 226 maximum            .
!     ...................................................
      real  x, bj, by, t, a, ri, alpha, rn
      integer nc, nst, n, i, k
      complex hankel(nc)
      dimension bj(226),by(226),t(3)
!     ................................................
!     .  by(*) calculation - obtain the zeroeth and  .
!     .                      first order functions   .
!     ................................................
      a = sin(x)/x                                                      !eq 4.68
      by(1) = -cos(x)/x                                                 !eq 4.69a
      by(2) = by(1)/x-a                                                 !eq 4.69b
!     ...........................................................
!     .  obtain the higher order functions by upward recursion  .
!     ...........................................................
        do 10 n = 3,nc
        rn = real(n-2)
        by(n) = (2.0*rn+1.0)*by(n-1)/x-by(n-2)
10      continue
!     ................................................
!     .  bj(*) calculation - set the starting order  .
!     .                      for downward recursion  .
!     ................................................
      nst = nc+int((101.0+x)**.5)                                       !eq 4.21
!     ....................................................
!     .  the t(*) array is used to recur down to the     .
!     .    two highest order functions that are needed   .
!     .  set starting values for the two highest orders  .
!     .    nst and nst-1                                 .
!     ....................................................
      t(3) = 0.0
      t(2) = 1.0e-35
!     ...................................................
!     .  recur downward to obtain orders nc-1 and nc-2  .
!     ...................................................
        do 20 i = nst-1,nc-1,-1
        ri = real(i)
        t(1) = (2.0*ri+1.0)*t(2)/x-t(3)
        t(3) = t(2)
        t(2) = t(1)
20      continue
!     ...............................................
!     .  continue downward recursion to order zero  .
!     ...............................................
      bj(nc) = t(3)
      bj(nc-1) = t(2)
        do 30 i = nc-2,1,-1
        ri = real(i)
        bj(i) = (2.0*ri+1.0)*bj(i+1)/x-bj(i+2)
30      continue
!     ..................................................
!     .  calculate the scale factor and the functions  .
!     ..................................................
      alpha = a/bj(1)
        do 40 k = 1,nc
        hankel(k) = cmplx(bj(k)*alpha,by(k))
40      continue
      return
      end subroutine
      subroutine genlgp(theta,pnmllg,nc)
      implicit none
!     ........................................................
!     .  calculate associated Legendre functions (argument   .
!     .    cos(theta)) divided by sin(theta) for m = 1       .
!     .  generate first two orders by formula and remaining  .
!     .    orders by recursion                               .
!     .                                                      .
!     .  pnmllg = associated Legendre function/sin(theta)    .
!     .  nc = number of orders (0 to nc-1)                   .
!     .  the order of the associated Legendre functions is   .
!     .    incremented by one in the pnmllg(*) array         .
!     ........................................................
      real theta, pnmllg, costh, rn
      integer nc, n
      dimension pnmllg(nc)
      costh = cos(theta)
!     ..............................
!     .  calculate orders 0 and 1  .
!     ..............................
      pnmllg(1) = 0.0                                                   !eq 4.70a
      pnmllg(2) = 1.0                                                   !eq 4.70b
!     .................................................
!     .  recur upward to obtain all remaining orders  .
!     .................................................
      do 10 n = 3,nc
      rn = real(n-1)
      pnmllg(n) = ((2.0*rn-1.0)*costh*pnmllg(n-1) &
                  -rn*pnmllg(n-2))/(rn-1.0)                            !eq 4.71
10    continue
      return
      end subroutine
      end subroutine
  end module
