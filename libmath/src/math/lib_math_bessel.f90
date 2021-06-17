!    Copyright (C) 2020  Max Daiber-Huppert <max_daiber-huppert@gmx.de>
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
! Created on Thu Jan 30 13:07:51 2020
! 
! @author: Max Daiber-Huppert
!

!#define _DEBUG_

module lib_math_bessel
    use lib_math_factorial
    implicit none

    private

    ! --- interface ---
    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    interface lib_math_bessel_spherical_first_kind
        module procedure lib_math_bessel_spherical_first_kind_real
        module procedure lib_math_bessel_spherical_first_kind_cmplx
    end interface

    ! calculates the the spherical Bessel function of the second kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} Y_{n+1 / 2}(\rho) $$
    interface lib_math_bessel_spherical_second_kind
        module procedure lib_math_bessel_spherical_second_kind_real
        module procedure lib_math_bessel_spherical_second_kind_cmplx
    end interface

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ h_{n}^{(1)}(\rho)=j_{n}(\rho)+i y_{n}(\rho) $$
    interface lib_math_hankel_spherical_1
        module procedure lib_math_bessel_spherical_third_kind_1_real
        module procedure lib_math_bessel_spherical_third_kind_1_cmplx
    end interface

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(2)}(\rho)=j_{n}(\rho) - i y_{n}(\rho) $$
    interface lib_math_hankel_spherical_2
        module procedure lib_math_bessel_spherical_third_kind_2_real
        module procedure lib_math_bessel_spherical_third_kind_2_cmplx
    end interface

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Formula: S = x * j_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    interface lib_math_riccati_s
        module procedure lib_math_bessel_riccati_s_real
        module procedure lib_math_bessel_riccati_s_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Formula: C = -x * y_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    interface lib_math_riccati_c
        module procedure lib_math_bessel_riccati_c_real
        module procedure lib_math_bessel_riccati_c_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Formula: Xi = x * h^(1)_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    interface lib_math_riccati_xi
        module procedure lib_math_bessel_riccati_xi_real
        module procedure lib_math_bessel_riccati_xi_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Formula: Zeta = x * h^(2)_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR.complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    interface lib_math_riccati_zeta
        module procedure lib_math_bessel_riccati_zeta_real
        module procedure lib_math_bessel_riccati_zeta_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n(a x)
    !
    ! Argument
    ! ----
    !   a: integer, optional
    !       coefficient
    !   x: double precision .OR. complex
    !       control variable
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    interface lib_math_bessel_spherical_first_kind_derivative
        module procedure lib_math_bessel_spherical_first_kind_derivative_real
        module procedure lib_math_bessel_spherical_first_kind_derivative_cmplx
        module procedure lib_math_bessel_spherical_first_kind_derivative_coeff_real
        module procedure lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    interface lib_math_bessel_spherical_second_kind_derivative
        module procedure lib_math_bessel_spherical_second_kind_derivative_real
        module procedure lib_math_bessel_spherical_second_kind_derivative_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    interface lib_math_hankel_spherical_1_derivative
        module procedure lib_math_bessel_spherical_third_kind_1_derivative_real
        module procedure lib_math_bessel_spherical_third_kind_1_derivative_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    interface lib_math_hankel_spherical_2_derivative
        module procedure lib_math_bessel_spherical_third_kind_2_derivative_real
        module procedure lib_math_bessel_spherical_third_kind_2_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Formula: S' = [ x * j_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = -x j_{n+1}(x)+(n+1) j_{n}(x) $$
    interface lib_math_riccati_s_derivative
        module procedure lib_math_bessel_riccati_s_derivative_real
        module procedure lib_math_bessel_riccati_s_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Formula: C' = [ -x * y_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = x y_{n+1}(x)-(n+1) y_{n}(x) $$
    interface lib_math_riccati_c_derivative
        module procedure lib_math_bessel_riccati_c_derivative_real
        module procedure lib_math_bessel_riccati_c_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Formula: Xi' = [ x * h^(1)_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision .OR. compelx
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    interface lib_math_riccati_xi_derivative
        module procedure lib_math_bessel_riccati_xi_derivative_real
        module procedure lib_math_bessel_riccati_xi_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Formula: Zeta' = [ x * h^(2)_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    interface lib_math_riccati_zeta_derivative
        module procedure lib_math_bessel_riccati_zeta_derivative_real
        module procedure lib_math_bessel_riccati_zeta_derivative_cmplx
    end interface

    ! --- public functions ---
    public :: lib_math_bessel_spherical_first_kind
    public :: lib_math_bessel_spherical_second_kind
    public :: lib_math_hankel_spherical_1
    public :: lib_math_hankel_spherical_2
    public :: lib_math_bessel_spherical_first_kind_derivative
    public :: lib_math_bessel_spherical_second_kind_derivative
    public :: lib_math_hankel_spherical_1_derivative
    public :: lib_math_hankel_spherical_2_derivative
    public :: lib_math_riccati_s_derivative
    public :: lib_math_riccati_c_derivative
    public :: lib_math_riccati_xi_derivative
    public :: lib_math_riccati_zeta_derivative

    public :: lib_math_bessel_test_functions

    ! --- paraemters ---
    double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

    contains

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    !
    function lib_math_bessel_spherical_first_kind_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n-1)
        real(kind=8) sj(0:fnu+n-1)

        order = fnu+n-1

        call SPHJ( order, x, nm, sj, dj)

        if (order .ne. nm) then
#ifdef _DEBUG_
            print *, "lib_math_bessel_spherical_first_kind_real: WARNING"
            print *, "  calculated highest order / requested: ", nm, " / ", order
            print *, "  x = ", x
            print *, "  sj(highest + 1:requested) = 0"
#endif
            sj(nm+1:order) = 0d0
        end if

        rv = sj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_real

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !   kode: integer, optional (std_value = 1)
    !       A PARAMETER TO INDICATE THE SCALING OPTION
    !           KODE= 1  RETURNS
    !                    CY(I)=J(FNU+I-1,Z), I=1,...,N
    !               = 2  RETURNS
    !                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
    !
    ! Returns
    ! ----
    !   rv: array<cmplx>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    function lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: z

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) dj(0:fnu+n-1)
        complex(kind=8) sj(0:fnu+n-1)
        complex(kind=8) dy(0:fnu+n-1)
        complex(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call CSPHJY( order, z, nm, sj, dj, sy, dy)

        if (order .ne. nm) then
#ifdef _DEBUG_
            print *, "lib_math_bessel_spherical_first_kind_cmplx: WARNING"
            print *, "  calculated highest order / requested: ", nm, " / ", order
            print *, "  z = ", z
            print *, "  sj(highest + 1:requested) = 0"
#endif
            sj(nm+1:order) = 0d0
        end if

        rv = sj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_cmplx


    ! calculates the the spherical Bessel function of the second kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.10)
    !
    function lib_math_bessel_spherical_second_kind_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dy(0:fnu+n-1)
        real(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call SPHY( order, x, nm, sy, dy)

        if (order .ne. nm) then
#ifdef _DEBUG_
            print *, "lib_math_bessel_spherical_second_kind_real: WARNING"
            print *, "  calculated highest order / requested: ", nm, " / ", order
            print *, "  x = ", x
            print *, "  sy(highest + 1:requested) = 0"
#endif
            sy(nm+1:order) = 0d0
        end if

        rv = sy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_real

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !   kode: integer, optional (std_value = 1)
    !       A PARAMETER TO INDICATE THE SCALING OPTION
    !           KODE= 1  RETURNS
    !                    CY(I)=Y(FNU+I-1,Z), I=1,...,N
    !               = 2  RETURNS
    !                    CY(I)=Y(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
    !
    ! Returns
    ! ----
    !   rv: array<cmplx>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    function lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n) result (rv)
    !
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: z

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) dj(0:fnu+n-1)
        complex(kind=8) sj(0:fnu+n-1)
        complex(kind=8) dy(0:fnu+n-1)
        complex(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call CSPHJY( order, z, nm, sj, dj, sy, dy)

        if (order .ne. nm) then
#ifdef _DEBUG_
            print *, "lib_math_bessel_spherical_second_kind_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
            print *, "  z = ", z
            print *, "  sy(highest + 1:requested) = 0"
#endif
            sy(nm+1:order) = 0d0
        end if

        rv = sy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_cmplx

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: double precision
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(z)=i^{-n-1} z^{-1} e^{i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(-2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.16
    !
    function lib_math_bessel_spherical_third_kind_1_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = cmplx(cos(x), sin(x), kind=8) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * cmplx(0.0_8, -2.0_8 * x, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(-i-1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

    end function lib_math_bessel_spherical_third_kind_1_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: cmplx
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(z)=i^{-n-1} z^{-1} e^{i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(-2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.16
    !
    function lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = exp(cmplx(0.0, 1.0, kind=8) * x) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * (-2.0_8) * x * cmplx(0.0, 1.0, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(-i-1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n) &
             + cmplx(0,1, kind=8) * lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

    end function lib_math_bessel_spherical_third_kind_1_cmplx

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{\mathrm{n}}^{(2)}(z)=i^{n+1} z^{-1} e^{-i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.17
    !
    function lib_math_bessel_spherical_third_kind_2_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = cmplx(cos(x), -sin(x), kind=8) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * cmplx(0.0_8, 2.0_8 * x, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(i+1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   -lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

    end function lib_math_bessel_spherical_third_kind_2_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   z: complex
    !       z .ne. (0.0, 0.0)
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{\mathrm{n}}^{(2)}(z)=i^{n+1} z^{-1} e^{-i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.17
    !
    function lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = exp(-cmplx(0.0, 1.0, kind=8) * z) / z
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * (2.0_8) * z * cmplx(0.0, 1.0, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(i+1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n) &
             - cmplx(0,1, kind=8) * lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_spherical_third_kind_2_cmplx

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_s_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rctj( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_s_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = rf(fnu:order)

    end function lib_math_bessel_riccati_s_real

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_s_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_s_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_c_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rcty( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_c_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> riccatiBesselC[2, 4]
        ! sageMath
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> numerical_approx(C(2,4))
        rv = -rf(fnu:order)

    end function lib_math_bessel_riccati_c_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_c_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = -z * lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_c_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (page 101)
    !
    function lib_math_bessel_riccati_xi_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        double precision, dimension(n) :: j
        double precision, dimension(n) :: y
        double precision :: rv_1

        ! nz: integer
        !    number of components of Y set to zero due to
        !    underflow,
        !    NZ=0   , normal return, computation completed
        !    NZ .NE. 0, last NZ components of Y set to zero,
        !             Y(K)=0.0D0, K=N-NZ+1,...,N.
        integer :: nz

        rv_1 = sqrt(PI/2.D0 * x)

        call DBESJ (X, FNU+0.5D0, N, j, nz)

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_riccati_s_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        call DBESY (X, FNU+0.5D0, N, y)

        rv(:) = cmplx(rv_1 * j(:), rv_1 * y(:), kind=8)

    end function lib_math_bessel_riccati_xi_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (page 101)
    !
    function lib_math_bessel_riccati_xi_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_xi_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_zeta_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        double precision, dimension(n) :: j
        double precision, dimension(n) :: y
        double precision :: rv_1

        ! nz: integer
        !    number of components of Y set to zero due to
        !    underflow,
        !    NZ=0   , normal return, computation completed
        !    NZ .NE. 0, last NZ components of Y set to zero,
        !             Y(K)=0.0D0, K=N-NZ+1,...,N.
        integer :: nz

        rv_1 = sqrt(PI * x/2.D0)

        call DBESJ (X, FNU+0.5D0, N, j, nz)

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_riccati_zeta_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        call DBESY (X, FNU+0.5D0, N, y)

        rv(:) = cmplx(rv_1 * j(:), -rv_1 * y(:), kind=8)

    end function lib_math_bessel_riccati_zeta_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_zeta_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_zeta_cmplx

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_real(x, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: j_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n-1)
        real(kind=8) sj(0:fnu+n-1)

        order = fnu+n-1

        call SPHJ( order, x, nm, sj, dj)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = sj(fnu:order)
        rv = dj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n(a x)
    !
    ! Argument
    ! ----
    !   a: double precision
    !       coefficient
    !   x: double precision
    !       control variable
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(a x)\right) = \frac{n j_n(a x)}{x}-a j_{n+1}(a x) $$
    !
    ! Mathematica Source Code:
    !   >>> jdn2[n_, x_, t_] := FullSimplify[D[SphericalBesselJ[n, t y], y] /. y -> x];
    function lib_math_bessel_spherical_first_kind_derivative_coeff_real(a, x, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: a
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: j_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n)
        real(kind=8) sj(0:fnu+n)

        order = fnu+n-1

        call SPHJ( order, a * x, nm, sj, dj)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = sj(fnu:order)
        rv = a * dj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_coeff_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_cmplx(z, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: j_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = csj(fnu:order)
        rv = cdj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   a: double precision
    !       coefficient
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(a x)\right)= j'_{n}(a x) (a x)'$$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx(a, z, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        real(kind=8), intent(in) :: a
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: j_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, a * z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = csj(fnu:order)
        rv = a * cdj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_second_kind_derivative_real(x, fnu, n, y_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: y_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dy(0:fnu+n-1)
        real(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call SPHY( order, x, nm, sy, dy)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        y_n = sy(fnu:order)
        rv = dy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_second_kind_derivative_cmplx(z, fnu, n, y_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: y_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) :: order
        integer(kind=4) :: nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_derivative_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        y_n = csy(fnu:order)
        rv = cdy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
    !
    function lib_math_bessel_spherical_third_kind_1_derivative_real(x, fnu, n, h_1_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i
        complex(kind=8), dimension(n+2) :: buffer_h_1_n

        buffer_h_1_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu-1, n+2)

        h_1_n = buffer_h_1_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_1_n(i) - (buffer_h_1_n(i+1) + x * buffer_h_1_n(i+2))/x)
        end do

    end function lib_math_bessel_spherical_third_kind_1_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
    !
    function lib_math_bessel_spherical_third_kind_1_derivative_cmplx(z, fnu, n, h_1_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        complex(kind=8), dimension(n+2) :: buffer_h_1_n
        integer :: i

        buffer_h_1_n = lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu-1, n+2)
        h_1_n = buffer_h_1_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_1_n(i) - (buffer_h_1_n(i+1) + z * buffer_h_1_n(i+2))/z)
        end do

    end function lib_math_bessel_spherical_third_kind_1_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_third_kind_2_derivative_real(x, fnu, n, h_2_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i
        complex(kind=8), dimension(n+2) :: buffer_h_2_n

        buffer_h_2_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu-1, n+2)

        h_2_n = buffer_h_2_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_2_n(i) - (buffer_h_2_n(i+1) + x * buffer_h_2_n(i+2))/x)
        end do

    end function lib_math_bessel_spherical_third_kind_2_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_third_kind_2_derivative_cmplx(z, fnu, n, h_2_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        complex(kind=8), dimension(n+2) :: buffer_h_2_n
        integer :: i

        buffer_h_2_n = lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu-1, n+2)
        h_2_n = buffer_h_2_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_2_n(i) - (buffer_h_2_n(i+1) + z * buffer_h_2_n(i+2))/z)
        end do

    end function lib_math_bessel_spherical_third_kind_2_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = -x j_{n+1}(x)+(n+1) j_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_s_derivative_real(x, fnu, n, r_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, intent(inout), dimension(n) :: r_n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rctj( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_s_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = df(fnu:order)
        r_n = rf(fnu:order)

    end function lib_math_bessel_riccati_s_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = $$ \frac{1}{2}\left(x j_{n-1}(x)+j_{n}(x)-x j_{n+1}(x)\right) $$
    !
    ! Reference: Wolfram Alpha: https://www.wolframalpha.com/input/?i=derivative%5Bx+SphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[x SphericalBesselJ[n,x], x]
    !
    function lib_math_bessel_riccati_s_derivative_cmplx(z, fnu, n, s_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: s_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: jf
        complex(kind=8), dimension(n+2) :: buffer_r_n

        order = fnu+n-1

        jf = lib_math_bessel_spherical_first_kind_cmplx(z, fnu-1, n+2)

        buffer_r_n = z * jf

        s_n = buffer_r_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_r_n(i-1) + jf(i) - buffer_r_n(i+1)) / 2
        end do

    end function lib_math_bessel_riccati_s_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = x y_{n+1}(x)-(n+1) y_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_c_derivative_real(x, fnu, n, r_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, intent(inout), dimension(n) :: r_n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rcty( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_c_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> derivative[ riccatiBesselC[n,x], x]
        !   >>> N[1/2*(-x sphericalBesselY[n-1,x] - sphericalBesselY[n,x] + x sphericalBesselY[n+1,x] ), 20] where n = 3, x=4
        ! sageMath
        !   >>> var('x')
        !   >>> var('n')
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> C_(n,x) = derivative(C(n,x), x)
        !   >>> n=5
        !   >>> x=4.0
        !   >>> for i in range(1,6):
        !   >>>     value = numerical_approx(C_(i,x))
        !   >>>     print("n = {}: {}".format(i, value))
        rv = -df(fnu:order)
        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> riccatiBesselC[2, 4]
        ! sageMath
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> numerical_approx(C(2,4))
        r_n = -rf(fnu:order)

    end function lib_math_bessel_riccati_c_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   c_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = \frac{1}{2}\left(-x y_{n-1}(x)-y_{n}(x)+x y_{n+1}(x)\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5B-x+SphericalBesselY%5Bn,x%5D,+x%5D
    !               >>> derivative[-x SphericalBesselY[n,x], x]
    !
    function lib_math_bessel_riccati_c_derivative_cmplx(z, fnu, n, c_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: c_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: yf
        complex(kind=8), dimension(n+2) :: buffer_c_n

        order = fnu+n-1

        yf = lib_math_bessel_spherical_second_kind_cmplx(z, fnu-1, n+2)

        buffer_c_n = z * yf

        c_n = -buffer_c_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (-buffer_c_n(i-1) - yf(i) + buffer_c_n(i+1)) / 2
        end do

    end function lib_math_bessel_riccati_c_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_xi_derivative_real(x, fnu, n, xi_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: xi_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_xi_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu-1, n+2)

        buffer_xi_n = x * h_n

        xi_n = buffer_xi_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_xi_n(i-1) + h_n(i) - buffer_xi_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_xi_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_xi_derivative_cmplx(x, fnu, n, xi_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: xi_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_xi_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu-1, n+2)

        buffer_xi_n = x * h_n

        xi_n = buffer_xi_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_xi_n(i-1) + h_n(i) - buffer_xi_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_xi_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_zeta_derivative_real(x, fnu, n, zeta_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: zeta_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_zeta_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu-1, n+2)

        buffer_zeta_n = x * h_n

        zeta_n = buffer_zeta_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_zeta_n(i-1) + h_n(i) - buffer_zeta_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_zeta_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_zeta_derivative_cmplx(z, fnu, n, zeta_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: zeta_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_zeta_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu-1, n+2)

        buffer_zeta_n = z * h_n

        zeta_n = buffer_zeta_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_zeta_n(i-1) + h_n(i) - buffer_zeta_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_zeta_derivative_cmplx

    function lib_math_bessel_test_functions() result (rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        ! test: Bessel functions
        if (.not. test_lib_math_bessel_spherical_first_kind_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_real()) then
            rv = rv + 1
        end if

        ! test: derivatives of the Bessel functions
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_coeff_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_derivative_real()) then
            rv = rv + 1
        end if


        ! test: Bessel functions with a complex argument
        if (.not. test_lib_math_bessel_spherical_first_kind_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_cmplx()) then
            rv = rv + 1
        end if

        ! test: derivatives of Bessel functions with complex argument
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_coeff_c()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_derivative_cmplx()) then
            rv = rv + 1
        end if


        print *, "-------------lib_math_bessel_test_functions----------------"
        if (rv == 0) then
            print *, "lib_math_bessel_test_functions tests: OK"
        else
            print *, rv,"lib_math_bessel_test_functions test(s) FAILED"
        end if
        print *, "-----------------------------------------------------------"

        contains

        function test_lib_math_bessel_spherical_first_kind_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 5
            double precision, dimension(n) :: j
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_j

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_j / -0.1917848549326277_8, &
                                  -0.09508940807917079_8, &
                                  0.1347312100851252_8, &
                                  0.229820618164296_8, &
                                  0.1870176553448892_8, &
                                  0.1068111614565045_8, &
                                  0.0479668998594208_8, &
                                  0.01790277817798953_8, &
                                  0.005741434674547791_8, &
                                  0.001618099715472961_8, &
                                  0.0004073442442494604_8, &
                                  0.00009274611037477283_8, &
                                  0.0000192878634744946_8, &
                                  3.693206997700175d-6, &
                                  6.554543130863411d-7, &
                                  1.084280182006037d-7, &
                                  1.679939975740199d-8, &
                                  2.448020198249413d-9, &
                                  3.367416303439037d-10, &
                                  4.386786629547396d-11, &
                                  5.427726760793208d-12, &
                                  6.394931430303462d-13, &
                                  7.191426926776938d-14, &
                                  7.735280379578157d-15, &
                                  7.973663002653014d-16, &
                                  7.890936302179649d-17, &
                                  7.509202557022752d-18, &
                                  6.881840826446885d-19, &
                                  6.082235206882156d-20, &
                                  5.190730939877244d-21, &
                                  4.282730217299213d-22, &
                                  3.419992522779504d-23, &
                                  2.646036140296298d-24, &
                                  1.98544596056825d-25, &
                                  1.44614468651574d-26, &
                                  1.023370682347077d-27, &
                                  7.04168241710901d-29, &
                                  4.714950550838944d-30, &
                                  3.07434091494052d-31, &
                                  1.953445816945712d-32, &
                                  1.210347583370466d-33   /

            j = lib_math_bessel_spherical_first_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_real:"
            do i=1, n
                buffer = (j(i) - ground_truth_j(i)) / ground_truth_j(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function

        function test_lib_math_bessel_spherical_first_kind_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(2, 2, kind=8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / (0.4779120939483495_8,-1.2325653366101572_8), &
                                  (1.0272168572412931_8,0.0054478926092466_8), &
                                  (0.2965864684395553_8,0.4662386131361224_8), &
                                  (-0.0736855052716959_8,0.2066172882614622_8), &
                                  (-0.06395584820746431_8,0.02429127554690435_8), &
                                  (-0.01555978321456398_8,-0.00806125981413275_8), &
                                  (-0.0010020201214517_8,-0.003670336195718445_8), &
                                  (0.0003746251837610103_8,-0.000610767427234174_8), &
                                  (0.00011648670842733624_8,-0.00002488609551349602_8), &
                                  (0.00001467742112281068_8,9.93301048563694d-6), &
                                  (4.12841712789929d-7,2.350144986920765d-6), &
                                  (-1.717409493295354d-7,2.378317035499526d-7), &
                                  (-3.28198760225293d-8,4.89776713629127d-9), &
                                  (-2.77223120945226d-9,-2.09643380732408d-9), &
                                  (-4.36128407109971d-11,-3.36134671926056d-10), &
                                  (1.906174283362537d-11,-2.434946898509715d-11), &
                                  (2.632963037090804d-12,-3.02219669043536d-13), &
                                  (1.668899527645897d-13,1.342116594888402d-13), &
                                  (1.67607012670835d-15,1.628460288072816d-14), &
                                  (-7.53727445802026d-16,9.172684858430605d-16), &
                                  (-8.154498630825985d-17,7.60745281142975d-18), &
                                  (-4.132272540482617d-18,-3.455984866242127d-18), &
                                  (-2.87808140311425d-20,-3.373603133444852d-19), &
                                  (1.318485750680503d-20,-1.553450103297935d-20), &
                                  (1.1725025985941834d-21,-9.21494979811842d-23), &
                                  (4.94679757042115d-23,4.251284993109568d-23), &
                                  (2.52928255983108d-25,3.47164437395756d-24), &
                                  (-1.173883574976466d-25,1.351386320658078d-25), &
                                  (-8.861980670892263d-27,6.01732539937958d-28), &
                                  (-3.201783684522278d-28,-2.807188114771226d-28), &
                                  (-1.25273306565584d-30,-1.970407455515633d-29), &
                                  (5.870522348423053d-31,-6.64146237759846d-31), &
                                  (3.850251970457037d-32,-2.301388327556d-33), &
                                  (1.216150034178248d-33,1.082732237792418d-33), &
                                  (3.75835093828684d-36,6.66402380933398d-35), &
                                  (-1.774373382688831d-36,1.9803156322458d-36), &
                                  (-1.0287600865064304d-37,5.49192174990615d-39), &
                                  (-2.886203249616963d-39,-2.600902435777421d-39), &
                                  (-7.2229505016509d-42,-1.425314904147324d-40), &
                                  (3.430261976585573d-42,-3.786957549397847d-42), &
                                  (1.782129386084693d-43,-8.5952234401768d-45)      /

            y = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_cmplx

        function test_lib_math_bessel_spherical_second_kind_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 5
            double precision, dimension(n) :: y
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_Y(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.05673243709264525_8, &
                                  0.1804383675140986_8, &
                                  0.1649954576011044_8, &
                                  -0.0154429099129942_8, &
                                  -0.1866155314792963_8, &
                                  -0.3204650467497392_8, &
                                  -0.5184075713701299_8, &
                                  -1.027394638812598_8, &
                                  -2.563776345067666_8, &
                                  -7.689444934417465_8, &
                                  -26.6561144057187_8, &
                                  -104.2662355696011_8, &
                                  -452.9685692144462_8, &
                                  -2160.57661050263_8, &
                                  -11214.14512749976_8, &
                                  -62881.46512899596_8, &
                                  -378650.9386722752_8, &
                                  -2.43621473010802d6, &
                                  -1.667485217208387d7, &
                                  -1.209576913433126d8, &
                                  -9.267951403057543d8, &
                                  -7.478762459163873d9, &
                                  -6.339056200850355d10, &
                                  -5.630362956173681d11, &
                                  -5.229150616794757d12, &
                                  -5.068263974897125d13, &
                                  -5.11733774822712d14, &
                                  -5.373695373371776d15, &
                                  -5.859891533226682d16, &
                                  -6.6265393941447d17, &
                                  -7.760717569758479d18, &
                                  -9.401810041163897d19, &
                                  -1.176867347616893d21, &
                                  -1.520525741860796d22, &
                                  -2.025735820617298d23, &
                                  -2.780310175033264d24, &
                                  -3.927783090341061d25, &
                                  -5.706760210147617d26, &
                                  -8.520862484318015d27, &
                                  -1.306506062374827d29, &
                                  -2.055758716067908d30    /

            y = lib_math_bessel_spherical_second_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_real:"
            do i=1, n
                buffer = 1 - abs(y(i) / ground_truth_y(i))
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_real

        function test_lib_math_bessel_spherical_second_kind_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(5, 2, kind = 8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_Y(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y   / (-0.4238528355573756_8,-0.5260357629568812_8), &
                                    (0.4417021333808715_8,-0.4876486358757115_8), &
                                    (0.551426635055955_8,0.1824170616320225_8), &
                                    (0.096568090505994_8,0.4547576079529878_8), &
                                    (-0.215340784054175_8,0.3198092491014484_8), &
                                    (-0.2322152904581253_8,0.1751579204104719_8), &
                                    (-0.0921891030550149_8,0.1885501306452657_8), &
                                    (0.1946294869478129_8,0.3301067406713438_8), &
                                    (0.9370309562007486_8,0.4638333501104723_8), &
                                    (3.095645174459763_8,-0.069183214858768_8), &
                                    (9.113221782042157_8,-4.746830661870954_8), &
                                    (23.02578238815563_8,-30.3160765900454_8), &
                                    (34.1083455109168_8,-151.9957485695214_8), &
                                    (-138.0686872023444_8,-683.6455043318863_8), &
                                    (-1949.836966415688_8,-2773.398526460756_8), &
                                    (-15157.91319779761_8,-9283.67319514052_8), &
                                    (-98914.44868038918_8,-14439.66102468852_8), &
                                    (-580493.8681778458_8,152242.6229824536_8), &
                                    (-3.03658359691793d6,2.334337240141192d6), &
                                    (-1.283423060283487d7,2.248774274246745d7), &
                                    (-2.27782796320593d7,1.833960014426958d8), &
                                    (3.703833614210729d8,1.338339816760343d9), &
                                    (6.737593691594824d9,8.640400257911148d9), &
                                    (7.871908573584923d10,4.478947486553187d10), &
                                    (7.763381230423557d11,9.91503767845082d10), &
                                    (6.815059433928252d12,-1.830627190242992d12), &
                                    (5.271042678202382d13,-4.01662536792895d13), &
                                    (3.280346028844037d14,-5.578715267375585d14), &
                                    (9.4189811984141d14,-6.494263959083227d15), &
                                    (-1.660062519532185d16,-6.696770136983625d16), &
                                    (-4.42299594473656d17,-6.071815336634587d17), &
                                    (-7.189520872094576d18,-4.4581984104084d18), &
                                    (-9.702087538556945d19,-1.658084844546514d19), &
                                    (-1.1544413342389585d21,2.535595106637828d20), &
                                    (-1.20671471432732d22,8.27994549880366d21), &
                                    (-1.030018788586696d23,1.556722164499901d23), &
                                    (-4.86560861785316d23,2.401716731662492d24), &
                                    (6.07044768130982d24,3.252241374484194d25), &
                                    (2.472031105934223d26,3.867478695483118d26), &
                                    (5.329528500212056d27,3.789155198141813d27), &
                                    (9.298901340286066d28,2.21874176557107d28)     /

            y = lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_cmplx:"
            do i=1, n
!                buffer = abs(y(i) - ground_truth_y(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
                buffer = 1 - abs(real(y(i)) / real(ground_truth_y(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(y(i)) / aimag(ground_truth_y(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_cmplx

        function test_lib_math_bessel_spherical_third_kind_1_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 2
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-14.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.454648713412841_8, + 0.208073418273571_8, kind=8)
            ground_truth_y(2)  = cmplx(0.435397774979992_8, - 0.350612004276055_8, kind=8)
            ground_truth_y(3)  = cmplx(0.198447949057147_8, - 0.733991424687654_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0607220976628748_8, - 1.48436655744308_8, kind=8)
            ground_truth_y(5)  = cmplx(0.0140793927629153_8, - 4.46129152636313_8, kind=8)
            ground_truth_y(6)  = cmplx(0.00263516977024412_8, - 18.5914453111910_8, kind=8)
            ground_truth_y(7)  = cmplx(0.000414040973427324_8, - 97.7916576851873_8, kind=8)
            ground_truth_y(8)  = cmplx(0.0000560965570334895_8, - 617.054329642526_8, kind=8)
            ground_truth_y(9)  = cmplx(6.68320432384702d-6, - 4530.11581463376_8, kind=8)
            ground_truth_y(10) = cmplx(7.10679719210186d-7, - 37888.9300947444_8, kind=8)
            ground_truth_y(11) = cmplx(6.82530086497472d-8, - 355414.720085438_8, kind=8)
            ground_truth_y(12) = cmplx(5.97687161216011d-9, - 3.69396563080236d6, kind=8)
            ground_truth_y(13) = cmplx(4.81014890094075d-10, - 4.21251900341417d7, kind=8)
            ground_truth_y(14) = cmplx(3.58145140158186d-11, - 5.22870909795969d8, kind=8)
            ground_truth_y(15) = cmplx(2.48104911947672d-12, - 7.01663209221144d9, kind=8)
            ground_truth_y(16) = cmplx(1.60698216593841d-13, - 1.01218294427270d11, kind=8)
            ground_truth_y(17) = cmplx(9.77323772781462d-15, - 1.56186693153047d12, kind=8)
            ground_truth_y(18) = cmplx(5.60205915100117d-16, - 2.56695860758255d13, kind=8)
            ground_truth_y(19) = cmplx(3.03657864374240d-17, - 4.47655889395416d14, kind=8)
            ground_truth_y(20) = cmplx(1.56113399222733d-18, - 8.25596436773937d15, kind=8)
            ground_truth_y(21) = cmplx(7.63264110088761d-20, - 1.60543649281522d17, kind=8)
            ground_truth_y(22) = cmplx(3.55743345463290d-21, - 3.28288884590347d18, kind=8)
            ground_truth_y(23) = cmplx(1.58408265731184d-22, - 7.04215665376430d19, kind=8)
            ground_truth_y(24) = cmplx(6.75252431874067d-24, - 1.58120235825106d21, kind=8)
            ground_truth_y(25) = cmplx(2.76055759221868d-25, - 3.70878338523624d22, kind=8)
            ground_truth_y(26) = cmplx(1.08417821950868d-26, - 9.07070727024627d23, kind=8)
            ground_truth_y(27) = cmplx(4.09686752846144d-28, - 2.30932157052756d25, kind=8)
            ground_truth_y(28) = cmplx(1.49167553359982d-29, - 6.11063145462780d26, kind=8)
            ground_truth_y(29) = cmplx(5.24018893806175d-31, - 1.67811432845212d28, kind=8)
            ground_truth_y(30) = cmplx(1.77831374777941d-32, - 4.77651520463390d29, kind=8)
            ground_truth_y(31) = cmplx(5.83661788752249d-34, - 1.40739387103855d31, kind=8)
            ground_truth_y(32) = cmplx(1.85470791494552d-35, - 4.28777479146294d32, kind=8)
            ground_truth_y(33) = cmplx(5.71204455589985d-37, - 1.34924166543979d34, kind=8)
            ground_truth_y(34) = cmplx(1.70656572193294d-38, - 4.38074763788785d35, kind=8)
            ground_truth_y(35) = cmplx(4.95061257549874d-40, - 1.46620121702699d37, kind=8)
            ground_truth_y(36) = cmplx(1.39561661412576d-41, - 5.05401345110523d38, kind=8)
            ground_truth_y(37) = cmplx(3.82640464769116d-43, - 1.79270857392533d40, kind=8)
            ground_truth_y(38) = cmplx(1.02108228151588d-44, - 6.53833228137634d41, kind=8)
            ground_truth_y(39) = cmplx(2.65390799339793d-46, - 2.45008189694220d43, kind=8)
            ground_truth_y(40) = cmplx(6.72295942323183d-48, - 9.42627697094610d44, kind=8)
            ground_truth_y(41) = cmplx(1.66097877863811d-49, - 3.72092932162677d46, kind=8)

            y = lib_math_bessel_spherical_third_kind_1_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_real:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_real

        function test_lib_math_bessel_spherical_third_kind_1_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(5, 2, kind=8)

            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(-0.0250227739998745_8, + 0.00233120915731335_8, kind=8)
            ground_truth_y(2)  = cmplx(-0.00182228917664341_8, + 0.0271504151649199_8, kind=8)
            ground_truth_y(3)  = cmplx(0.0296975379081458_8, + 0.0120891343783301_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0315922819865381_8, - 0.0269692779105477_8, kind=8)
            ground_truth_y(5)  = cmplx(-0.00458857312258838_8, - 0.0598897093673198_8, kind=8)
            ground_truth_y(6)  = cmplx(-0.0758854047150979_8, - 0.0631149498592040_8, kind=8)
            ground_truth_y(7)  = cmplx(-0.187212328816476_8, - 0.00224281954730325_8, kind=8)
            ground_truth_y(8)  = cmplx(-0.345739239467690_8, + 0.225933476709331_8, kind=8)
            ground_truth_y(9)  = cmplx(-0.473216590452379_8, + 0.944214817382837_8, kind=8)
            ground_truth_y(10) = cmplx(0.0657355705905944_8, + 3.09639836994315_8, kind=8)
            ground_truth_y(11) = cmplx(4.74590684093329_8, + 9.11302288820878_8, kind=8)
            ground_truth_y(12) = cmplx(30.3158913122633_8, + 23.0256470073922_8, kind=8)
            ground_truth_y(13) = cmplx(151.995722926319_8, + 34.1083014388942_8, kind=8)
            ground_truth_y(14) = cmplx(683.645503092378_8, - 138.068697574778_8, kind=8)
            ground_truth_y(15) = cmplx(2773.39852701965_8, - 1949.83696832108_8, kind=8)
            ground_truth_y(16) = cmplx(9283.67319536370_8, - 15157.9131980699_8, kind=8)
            ground_truth_y(17) = cmplx(14439.6610247403_8, - 98914.4486804165_8, kind=8)
            ground_truth_y(18) = cmplx(-152242.622982445_8, - 580493.868177847_8, kind=8)
            ground_truth_y(19) = cmplx(-2.33433724014119d6, - 3.03658359691793d6, kind=8)
            ground_truth_y(20) = cmplx(-2.24877427424674d7, - 1.28342306028349d7, kind=8)
            ground_truth_y(21) = cmplx(-1.83396001442696d8, - 2.27782796320593d7, kind=8)
            ground_truth_y(22) = cmplx(-1.33833981676034d9, + 3.70383361421073d8, kind=8)
            ground_truth_y(23) = cmplx(-8.64040025791115d9, + 6.73759369159482d9, kind=8)
            ground_truth_y(24) = cmplx(-4.47894748655319d10, + 7.87190857358492d10, kind=8)
            ground_truth_y(25) = cmplx(-9.91503767845082d10, + 7.76338123042356d11, kind=8)
            ground_truth_y(26) = cmplx(1.83062719024299d12, + 6.81505943392825d12, kind=8)
            ground_truth_y(27) = cmplx(4.01662536792895d13, + 5.27104267820238d13, kind=8)
            ground_truth_y(28) = cmplx(5.57871526737558d14, + 3.28034602884404d14, kind=8)
            ground_truth_y(29) = cmplx(6.49426395908323d15, + 9.41898119841410d14, kind=8)
            ground_truth_y(30) = cmplx(6.69677013698362d16, - 1.66006251953219d16, kind=8)
            ground_truth_y(31) = cmplx(6.07181533663459d17, - 4.42299594473656d17, kind=8)
            ground_truth_y(32) = cmplx(4.45819841040840d18, - 7.18952087209458d18, kind=8)
            ground_truth_y(33) = cmplx(1.65808484454651d19, - 9.70208753855694d19, kind=8)
            ground_truth_y(34) = cmplx(-2.53559510663783d20, - 1.15444133423896d21, kind=8)
            ground_truth_y(35) = cmplx(-8.27994549880366d21, - 1.20671471432732d22, kind=8)
            ground_truth_y(36) = cmplx(-1.55672216449990d23, - 1.03001878858670d23, kind=8)
            ground_truth_y(37) = cmplx(-2.40171673166249d24, - 4.86560861785316d23, kind=8)
            ground_truth_y(38) = cmplx(-3.25224137448419d25, + 6.07044768130982d24, kind=8)
            ground_truth_y(39) = cmplx(-3.86747869548312d26, + 2.47203110593422d26, kind=8)
            ground_truth_y(40) = cmplx(-3.78915519814181d27, + 5.32952850021206d27, kind=8)
            ground_truth_y(41) = cmplx(-2.21874176557107d28, + 9.29890134028607d28, kind=8)

            y = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_cmplx:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_cmplx

        function test_lib_math_bessel_spherical_third_kind_2_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 2
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=2.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.454648713412841_8, - 0.208073418273571_8, kind=8)
            ground_truth_y(2)  = cmplx(0.435397774979992_8, + 0.350612004276055_8, kind=8)
            ground_truth_y(3)  = cmplx(0.198447949057147_8, + 0.733991424687654_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0607220976628748_8, + 1.48436655744308_8, kind=8)
            ground_truth_y(5)  = cmplx(0.0140793927629153_8, + 4.46129152636313_8, kind=8)
            ground_truth_y(6)  = cmplx(0.00263516977024412_8, + 18.5914453111910_8, kind=8)
            ground_truth_y(7)  = cmplx(0.000414040973427324_8, + 97.7916576851873_8, kind=8)
            ground_truth_y(8)  = cmplx(0.0000560965570334895_8, + 617.054329642526_8, kind=8)
            ground_truth_y(9)  = cmplx(6.68320432384702d-6, + 4530.11581463376_8, kind=8)
            ground_truth_y(10) = cmplx(7.10679719210186d-7, + 37888.9300947444_8, kind=8)
            ground_truth_y(11) = cmplx(6.82530086497472d-8, + 355414.720085438_8, kind=8)
            ground_truth_y(12) = cmplx(5.97687161216011d-9, + 3.69396563080236d6, kind=8)
            ground_truth_y(13) = cmplx(4.81014890094075d-10, + 4.21251900341417d7, kind=8)
            ground_truth_y(14) = cmplx(3.58145140158186d-11, + 5.22870909795969d8, kind=8)
            ground_truth_y(15) = cmplx(2.48104911947672d-12, + 7.01663209221144d9, kind=8)
            ground_truth_y(16) = cmplx(1.60698216593841d-13, + 1.01218294427270d11, kind=8)
            ground_truth_y(17) = cmplx(9.77323772781462d-15, + 1.56186693153047d12, kind=8)
            ground_truth_y(18) = cmplx(5.60205915100117d-16, + 2.56695860758255d13, kind=8)
            ground_truth_y(19) = cmplx(3.03657864374240d-17, + 4.47655889395416d14, kind=8)
            ground_truth_y(20) = cmplx(1.56113399222733d-18, + 8.25596436773937d15, kind=8)
            ground_truth_y(21) = cmplx(7.63264110088761d-20, + 1.60543649281522d17, kind=8)
            ground_truth_y(22) = cmplx(3.55743345463290d-21, + 3.28288884590347d18, kind=8)
            ground_truth_y(23) = cmplx(1.58408265731184d-22, + 7.04215665376430d19, kind=8)
            ground_truth_y(24) = cmplx(6.75252431874067d-24, + 1.58120235825106d21, kind=8)
            ground_truth_y(25) = cmplx(2.76055759221868d-25, + 3.70878338523624d22, kind=8)
            ground_truth_y(26) = cmplx(1.08417821950868d-26, + 9.07070727024627d23, kind=8)
            ground_truth_y(27) = cmplx(4.09686752846144d-28, + 2.30932157052756d25, kind=8)
            ground_truth_y(28) = cmplx(1.49167553359982d-29, + 6.11063145462780d26, kind=8)
            ground_truth_y(29) = cmplx(5.24018893806175d-31, + 1.67811432845212d28, kind=8)
            ground_truth_y(30) = cmplx(1.77831374777941d-32, + 4.77651520463390d29, kind=8)
            ground_truth_y(31) = cmplx(5.83661788752249d-34, + 1.40739387103855d31, kind=8)
            ground_truth_y(32) = cmplx(1.85470791494552d-35, + 4.28777479146294d32, kind=8)
            ground_truth_y(33) = cmplx(5.71204455589985d-37, + 1.34924166543979d34, kind=8)
            ground_truth_y(34) = cmplx(1.70656572193294d-38, + 4.38074763788785d35, kind=8)
            ground_truth_y(35) = cmplx(4.95061257549874d-40, + 1.46620121702699d37, kind=8)
            ground_truth_y(36) = cmplx(1.39561661412576d-41, + 5.05401345110523d38, kind=8)
            ground_truth_y(37) = cmplx(3.82640464769116d-43, + 1.79270857392533d40, kind=8)
            ground_truth_y(38) = cmplx(1.02108228151588d-44, + 6.53833228137634d41, kind=8)
            ground_truth_y(39) = cmplx(2.65390799339793d-46, + 2.45008189694220d43, kind=8)
            ground_truth_y(40) = cmplx(6.72295942323183d-48, + 9.42627697094610d44, kind=8)
            ground_truth_y(41) = cmplx(1.66097877863811d-49, + 3.72092932162677d46, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_real:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_real

        function test_lib_math_bessel_spherical_third_kind_2_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(2_8, 2_8, kind=8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=2.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel2(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.910979344197223_8, - 2.44844550451690_8, kind=8)
            ground_truth_y(2)  = cmplx(2.06407896443698_8, + 0.0711231320186915_8, kind=8)
            ground_truth_y(3)  = cmplx(0.690422228144533_8, + 0.953728630203184_8, kind=8)
            ground_truth_y(4)  = cmplx(-0.00889039150233644_8, + 0.258009870554623_8, kind=8)
            ground_truth_y(5)  = cmplx(-0.254463139803031_8, - 0.486653171603505_8, kind=8)
            ground_truth_y(6)  = cmplx(-1.65862130916237_8, - 0.780437442105689_8, kind=8)
            ground_truth_y(7)  = cmplx(-6.45294842618414_8, + 2.90165880600938_8, kind=8)
            ground_truth_y(8)  = cmplx(-9.88306995640557_8, + 31.1829109467346_8, kind=8)
            ground_truth_y(9)  = cmplx(86.3273521399181_8, + 151.095769580766_8, kind=8)
            ground_truth_y(10) = cmplx(1018.93133726931_8, + 244.082863176871_8, kind=8)
            ground_truth_y(11) = cmplx(5912.99009997946_8, - 3831.62602151987_8, kind=8)
            ground_truth_y(12) = cmplx(9908.23007464352_8, - 51403.3175010484_8, kind=8)
            ground_truth_y(13) = cmplx(-244509.742801807_8, - 348709.772538709_8, kind=8)
            ground_truth_y(14) = cmplx(-3.71753020095287d6, - 599846.868354584_8, kind=8)
            ground_truth_y(15) = cmplx(-2.88977854750235d7, + 2.13930722675771d7, kind=8)
            ground_truth_y(16) = cmplx(-5.06916405530332d7, + 3.65208565502209d8, kind=8)
            ground_truth_y(17) = cmplx(2.46640395383114d9, + 3.20183352466055d9, kind=8)
            ground_truth_y(18) = cmplx(4.68136508381095d10, + 5.70208539384046d9, kind=8)
            ground_truth_y(19) = cmplx(4.57046288075731d11, - 3.62928031162014d11, kind=8)
            ground_truth_y(20) = cmplx(8.23780225613767d11, - 7.59046453834298d12, kind=8)
            ground_truth_y(21) = cmplx(-6.64322183371856d13, - 8.16759584174163d13, kind=8)
            ground_truth_y(22) = cmplx(-1.51893259196028d15, - 1.48657871284022d14, kind=8)
            ground_truth_y(23) = cmplx(-1.78601652615391d16, + 1.48121292056872d16, kind=8)
            ground_truth_y(24) = cmplx(-3.27714730363732d16, + 3.67711970627580d17, kind=8)
            ground_truth_y(25) = cmplx(3.95341101195822d18, + 4.69086833384576d18, kind=8)
            ground_truth_y(26) = cmplx(1.05925193459135d20, + 8.66614022249484d18, kind=8)
            ground_truth_y(27) = cmplx(1.45708609342882d21, - 1.24474379710101d21, kind=8)
            ground_truth_y(28) = cmplx(2.70761023288441d21, - 3.58079121897428d22, kind=8)
            ground_truth_y(29) = cmplx(-4.56586238000232d23, - 5.28343689514023d23, kind=8)
            ground_truth_y(30) = cmplx(-1.40379590773110d25, - 9.86735771881784d23, kind=8)
            ground_truth_y(31) = cmplx(-2.21157662787594d26, + 1.93033887444595d26, kind=8)
            ground_truth_y(32) = cmplx(-4.14849614903414d26, + 6.31740787681276d27, kind=8)
            ground_truth_y(33) = cmplx(9.31864502878598d28, + 1.05840021607085d29, kind=8)
            ground_truth_y(34) = cmplx(3.23459501790776d30, + 1.99303126060599d29, kind=8)
            ground_truth_y(35) = cmplx(5.74246074611821d31, - 5.09469792100470d31, kind=8)
            ground_truth_y(36) = cmplx(1.08504492314173d32, - 1.86960917320476d33, kind=8)
            ground_truth_y(37) = cmplx(-3.13170326932691d34, - 3.50605705837511d34, kind=8)
            ground_truth_y(38) = cmplx(-1.21149976429793d36, - 6.64499573280903d34, kind=8)
            ground_truth_y(39) = cmplx(-2.39302402477947d37, + 2.15047444512683d37, kind=8)
            ground_truth_y(40) = cmplx(-4.54792943188346d37, + 8.74689905414290d38, kind=8)
            ground_truth_y(41) = cmplx(1.64008398093830d40, + 1.81518369502779d40, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_cmplx:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_cmplx

        function test_lib_math_bessel_riccati_s_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: s
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_s
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_s /-0.7568024953079283_8, &
                                  0.4644429970366299_8, &
                                  1.105134743085401_8, &
                                  0.9169754318201209_8, &
                                  0.499572262599811_8, &
                                  0.2070621590294538_8, &
                                  0.06984867473118706_8, &
                                  0.01994603384690409_8, &
                                  0.004948952194703267_8, &
                                  0.0010870129805848_8, &
                                  0.0002143594630745305_8, &
                                  0.00003837420055648577_8, &
                                  6.29219012526265d-6, &
                                  9.519877264057912d-7, &
                                  1.337270279764404d-7, &
                                  1.753322642340162d-8, &
                                  2.155476804922155d-9, &
                                  2.494572172061618d-10, &
                                  2.72738456317607d-11, &
                                  2.825854887624648d-12, &
                                  2.782395225796244d-13, &
                                  2.610021881650173d-14, &
                                  2.337829697769169d-15, &
                                  2.003652834014281d-16, &
                                  1.646238219761019d-17, &
                                  1.298898519296812d-18, &
                                  9.857392342415918d-20, &
                                  7.205966073297192d-21, &
                                  5.081100836772061d-22, &
                                  3.460261910299511d-23, &
                                  2.278548091971816d-24, &
                                  1.45239299575079d-25, &
                                  8.970876335678207d-27, &
                                  5.374408796918764d-28, &
                                  3.125839916072345d-29, &
                                  1.766505830603009d-30, &
                                  9.707933247997202d-32, &
                                  5.19198715647987d-33, &
                                  2.704267040255438d-34, &
                                  1.372689601184775d-35, &
                                  6.794922084492532d-37   /

            s = lib_math_bessel_riccati_s_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_real:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_s_real

        function test_lib_math_bessel_riccati_s_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: s
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_s
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0 + I*2
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_s /(-2.847239086848828_8,-2.370674169352002_8), &
                                 (1.652619979612418_8,-2.934227931977672_8), &
                                 (2.958542695022977_8,0.114351416281673_8), &
                                 (1.363098423551396_8,1.569308000747856_8), &
                                 (0.0483106984724765_8,1.1285108882793486_8), &
                                 (-0.2604793668495243_8,0.4185319695297422_8), &
                                 (-0.1609801390587136_8,0.078786748220561_8), &
                                 (-0.05564622201640157_8,-0.0044122433799559_8), &
                                 (-0.01257689206042502_8,-0.00855414533582638_8), &
                                 (-0.001657258059948338_8,-0.003291134259131266_8), &
                                 (0.0000261563402719274_8,-0.0008033745349705883_8), &
                                 (0.000080028165652198_8,-0.0001379671023162519_8), &
                                 (0.0000246488864008038_8,-0.00001533891668422602_8), &
                                 (4.868974641255984d-6,-3.49697106887638d-7), &
                                 (6.993944733818892d-7,3.043207756416122d-7), &
                                 (7.004355371964893d-8,8.651363280151069d-8), &
                                 (3.06782136461726d-9,1.492873119684237d-8), &
                                 (-5.31119763595177d-10,1.892182594411971d-9), &
                                 (-1.630206293416059d-10,1.754661366245485d-10), &
                                 (-2.600818802187669d-11,9.4431451736298d-12), &
                                 (-3.01497105187608d-12,-3.77670984916983d-13), &
                                 (-2.630256416668016d-13,-1.786659372571313d-13), &
                                 (-1.531299666407742d-14,-2.7845816327099d-14), &
                                 (-9.7501781840681d-17,-3.037924698411361d-15), &
                                 (1.182338322416188d-16,-2.524174633165996d-16), &
                                 (1.934776755720683d-17,-1.511222007524695d-17), &
                                 (2.041074458131367d-18,-4.00795992674089d-19), &
                                 (1.634029378129971d-19,4.60879248053659d-20), &
                                 (9.841444241113018d-21,9.047007561629002d-21), &
                                 (3.574696369766228d-22,9.517292228605612d-22), &
                                 (-8.10010991155721d-24,7.432640996354553d-23), &
                                 (-2.899877119993025d-24,4.463649155193212d-24), &
                                 (-3.173521226376781d-25,1.847952478450063d-25), &
                                 (-2.453136330424939d-26,1.47786393677701d-27), &
                                 (-1.466457262857802d-27,-6.31736953723452d-28), &
                                 (-6.473190388008502d-29,-7.727878444180707d-29), &
                                 (-1.415141776235862d-30,-6.025267801605198d-30), &
                                 (8.63789953234888d-32,-3.595904951070299d-31), &
                                 (1.389799278546981d-32,-1.643208992641743d-32), &
                                 (1.1230011393319822d-33,-4.782342079160211d-34), &
                                 (6.737497343894507d-35,4.28044062163919d-36)      /

            s = lib_math_bessel_riccati_s_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_cmplx:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_s_cmplx

        function test_lib_math_bessel_riccati_c_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: c
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_c
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_c /-0.6536436208636119_8, &
                                 -0.9202134005238312_8, &
                                 -0.03651642952926151_8, &
                                 0.8745678636122543_8, &
                                 1.567010190850707_8, &
                                 2.651205065801836_8, &
                                 5.723803740104341_8, &
                                 15.95115708953727_8, &
                                 54.09303534566043_8, &
                                 213.9442431295196_8, &
                                 962.1421195195575_8, &
                                 4837.301884348157_8, &
                                 26852.34371548235_8, &
                                 162989.8463374165_8, &
                                 1.073329119062079d6, &
                                 7.618646266862657d6, &
                                 5.797117944912351d7, &
                                 4.706435841884063d8, &
                                 4.060160182199432d9, &
                                 3.708583810115634d10, &
                                 3.575267613040749d11, &
                                 3.627563465265611d12, &
                                 3.863878049030124d13, &
                                 4.310587170506234d14, &
                                 5.026301144854523d15, &
                                 6.114113030741729d16, &
                                 7.745231102747159d17, &
                                 1.020129008083257d19, &
                                 1.394932155011731d20, &
                                 1.977577030810884d21, &
                                 2.902976798895937d22, &
                                 4.407263848008195d23, &
                                 6.912410792623947d24, &
                                 1.118859489953383d26, &
                                 1.867177234879293d27, &
                                 3.209692135267247d28, &
                                 5.67853176775057d29, &
                                 1.033122355479212d31, &
                                 1.931425884755771d32, &
                                 3.707663604600068d33, &
                                 7.303321360237576d34   /

            c = lib_math_bessel_riccati_c_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_real:"
            do i=1, n
                buffer = 1 - abs(c(i) / ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_c_real

        function test_lib_math_bessel_riccati_c_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: c
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_c
            double precision :: ground_truth_e = 2 * 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0 + I*2.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_c / (-2.459135213917384_8,2.744817006792154_8), &
                                  (-3.064584428953089_8,-1.575797246601832_8), &
                                  (0.14764538256498_8,-2.770920026067327_8), &
                                  (1.826769798484406_8,-1.268945470747985_8), &
                                  (1.521570505789599_8,-0.284342491918936_8), &
                                  (0.6561488692098297_8,-0.6122844699167391_8), &
                                  (-0.751555910436387_8,-1.784447098028703_8), &
                                  (-4.929975463781749_8,-3.050255301390585_8), &
                                  (-18.61375343299474_8,0.02864438952957_8), &
                                  (-58.30809074620009_8,34.79102706188217_8), &
                                  (-136.8540399849895_8,242.9626308634029_8), &
                                  (-6.2573523776097_8,1273.0455065328877_8), &
                                  (3036.074884073627_8,5627.438609656383_8), &
                                  (29255.2282968867_8,19273.96033156496_8), &
                                  (206981.85081434_8,19462.8307792003_8), &
                                  (1.2276817156859659d6,-506636.909173789_8), &
                                  (5.834070367999903d6,-6.966424986283187d6), &
                                  (1.428798025837888d7,-6.472420021469492d7), &
                                  (-1.3235290931078d8,-4.961109074209073d8), &
                                  (-2.829309866615508d9,-3.116790750250133d9), &
                                  (-3.40917479762657d10,-1.278054846472965d10), &
                                  (-3.291232722441548d11,3.80924600421564d10), &
                                  (-2.632570815142193d12,1.75560577547714d12), &
                                  (-1.546378807438846d13,2.760892818739198d13), &
                                  (-1.2965074603367d13,3.304481231356332d14), &
                                  (1.507601860325994d15,3.274311544098311d15), &
                                  (3.20894929248299d16,2.537876013900457d16), &
                                  (4.731484518795951d17,9.56662334277517d16), &
                                  (5.698707761603351d18,-1.575366677771509d18), &
                                  (5.5512529967101d19,-5.053748060116206d19), &
                                  (3.511780103033324d20,-9.222908312218367d20), &
                                  (-1.39711487471965d21,-1.3343596523155573d22), &
                                  (-1.02019483527651d23,-1.584052016498046d23), &
                                  (-2.354489981708474d24,-1.382797381994572d24), &
                                  (-4.071288873072953d25,-2.59599683963069d24), &
                                  (-5.773957526958108d26,2.464769732371248d26), &
                                  (-6.408320289566198d27,7.602078860947059d27), &
                                  (-3.74889047900571d28,1.575246125104232d29), &
                                  (6.25509342266883d29,2.636433894720829d30), &
                                  (2.997087376505044d31,3.562713543073535d31), &
                                  (7.543686660483393d32,3.235024031669992d32)   /

            c = lib_math_bessel_riccati_c_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_cmplx:"
            do i=1, n
!                buffer = abs(c(i) - ground_truth_c(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
                buffer = 1 - abs(real(c(i)) / real(ground_truth_c(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(c(i)) / aimag(ground_truth_c(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_c_cmplx

        function test_lib_math_bessel_riccati_xi_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 0.5_8
            complex(kind=8), dimension(n) :: xi
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_xi
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with mathematica
            !
            ! source code:
            !  >>> xin[n_, x_] := FullSimplify[x hn[n, x]];
            !  >>>
            !  >>> nMax = 40;
            !  >>> x = 0.5;
            !  >>>
            !  >>> For[n = 0, n <= nMax, n++,
            !  >>>  erg = N[xin[n, x], 16];
            !  >>>  Print[FortranForm[erg]];
            !  >>>  ]

            data ground_truth_xi /  (0.47942553860420295_8,-0.8775825618903728_8), &
                                    (0.0812685153180333_8,-2.2345906623849485_8), &
                                    (0.008185553303996706_8,-12.529961412419318_8), &
                                    (0.0005870177219337786_8,-123.06502346180824_8), &
                                    (0.000032694803076194855_8,-1710.380367052896_8), &
                                    (1.4887334377287222d-6,-30663.78158349032_8), &
                                    (5.733255383704708d-8,-672892.8144697341_8), &
                                    (1.9129620345019236d-9,-1.7464549394629598d7), &
                                    (5.6307198010606435d-11,-5.2326358902441823d8), &
                                    (1.4826978586953883d-12,-1.777349747743559d10), &
                                    (3.532061981830939d-14,-6.748696405535281d11), &
                                    (7.681736736061813d-16,-2.832675140577074d13), &
                                    (1.5369167574956977d-17,-1.3023556950249008d15), &
                                    (2.847051416677187d-19,-6.508945799983926d16), &
                                    (4.910075099825432d-21,-3.513528376296296d18), &
                                    (7.921412215629845d-23,-2.037195563671853d20), &
                                    (1.2004738650741575d-24,-1.2627098966389193d22), &
                                    (1.7152938595909842d-26,-8.331848122253196d23), &
                                    (2.3183663953195413d-28,-5.831030975680598d25), &
                                    (2.9727294547738398d-30,-4.314129737191417d27), &
                                    (3.625794040507704d-32,-3.364438091911737d29), &
                                    (4.216584424744342d-34,-2.758407822393905d31), &
                                    (4.685647724476791d-36,-2.371894283449567d33), &
                                    (4.985272847540434d-38,-2.1344290143223713d35), &
                                    (5.087522112610439d-40,-2.0061260840346838d37), &
                                    (4.988228177000061d-42,-1.9657901194525576d39), &
                                    (4.706279296883083d-44,-2.0049053092332057d41), &
                                    (4.278776960091815d-46,-2.1250030487752523d43), &
                                    (3.7535921782660104d-48,-2.3373028631218544d45), &
                                    (3.1812313151664884d-50,-2.6643127636540365d47), &
                                    (2.607736304099903d-52,-3.1436553308254513d49), &
                                    (2.069758353509157d-54,-3.834993072330685d51), &
                                    (1.5922132167825223d-56,-4.8317769056035805d53), &
                                    (1.1882830810796369d-58,-6.280926477977422d55), &
                                    (8.611186421470097d-61,-8.415958302799185d57), &
                                    (6.064508322338716d-63,-1.1613394365215078d60), &
                                    (4.15396251038509d-65,-1.649017840277513d62), &
                                    (2.7694282338524123d-67,-2.407449912861517d64), &
                                    (1.7984039359308915d-69,-3.611009967508248d66), &
                                    (1.138274811658968d-71,-5.560714604971416d68), &
                                    (7.0266490269757956d-74,-8.785567974858086d70) /

            xi = lib_math_bessel_riccati_xi_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_real:"
            do i=1, n
                buffer = 1 - abs(real(xi(i)) / real(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(xi(i)) / aimag(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_real

        function test_lib_math_bessel_riccati_xi_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(0.5, 2)
            complex(kind=8), dimension(n) :: xi
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_xi
            double precision :: ground_truth_e = 3* 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_xi / (0.06488319105786539_8,-0.11876788457694576_8), &
                                   (-0.16702533719458318_8,-0.10938914974120742_8), &
                                   (-0.2782650508788346_8,0.31596042541357816_8), &
                                   (0.7467763082389822_8,0.9499895196993936_8), &
                                   (4.02263454608541_8,-1.9935851304483712_8), &
                                   (-4.930935576635771_8,-20.097884794182935_8), &
                                   (-114.4399548151492_8,1.5094007346791365_8), &
                                   (-160.8608967049671_8,722.5096330222519_8), &
                                   (4930.635781963457_8,2408.996281339652_8), &
                                   (29294.102711349085_8,-35349.6033260506_8), &
                                   (-255516.74181304732_8,-343349.49854692706_8), &
                                   (-4.053671685772039d6,1.7121809965978751d6), &
                                   (7.81859943760688d6,4.8851344558873825d7), &
                                   (6.01771253077249d8,4.998472138118529d7), &
                                   (2.5387900761216574d9,-7.536111092211919d9), &
                                   (-9.47857682517896d10,-6.040844007482326d10), &
                                   (-1.2294806000859526d12,1.1699800664924927d12), &
                                   (1.3490610353448566d13,2.3695794487792234d13), &
                                   (4.470627265602757d14,-1.2579734976179488d14), &
                                   (-2.5780671411997475d14,-8.355435261794524d15), &
                                   (-1.5497675245486906d17,-3.3479394156975948d16), &
                                   (-1.3932324865086664d18,2.837006287280914d18), &
                                   (5.0514516104507154d19,4.257786269210637d19), &
                                   (1.1704718806373283d21,-8.47143956600987d20), &
                                   (-1.231532456916699d22,-3.061486957093495d22), &
                                   (-7.781075706797471d23,1.0833832120141485d23), &
                                   (-2.0562103906753614d24,1.9355226493093356d25), &
                                   (4.70699150609503d26,1.718614977914014d26), &
                                   (7.495936539523732d27,-1.1090111785971161d28), &
                                   (-2.4767977026220627d29,-2.7560832065211702d29), &
                                   (-9.378868891759496d30,4.97480009571606d30), &
                                   (7.574652929466099d31,3.0520582190041642d32), &
                                   (9.619249276476301d33,1.1477129018599958d31), &
                                   (7.38342836490342d34,-2.9445329976279795d35), &
                                   (-8.711570700093484d36,-4.648947842075915d36), &
                                   (-2.2174536166240226d38,2.454257606082996d38), &
                                   (6.356593374785566d39,9.463579737878725d39), &
                                   (3.7991546298753706d41,-1.3733765453381294d41), &
                                   (-1.5013726682075593d42,-1.4630047460478712d43), &
                                   (-5.441046582017436d44,-7.799100030416413d43), &
                                   (-7.954901461926335d45,1.9517663337780713d46) /

            xi = lib_math_bessel_riccati_xi_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_cmplx:"
            do i=1, n
!                buffer = abs(xi(i) - ground_truth_xi(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
                buffer = 1 - abs(real(xi(i)) / real(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(xi(i)) / aimag(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_cmplx

        function test_lib_math_bessel_riccati_zeta_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 0.5_8
            complex(kind=8), dimension(n) :: zeta
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_zeta
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_zeta / (0.47942553860420295_8,0.8775825618903728_8), &
                                     (0.0812685153180333_8,2.2345906623849485_8), &
                                     (0.008185553303996706_8,12.529961412419318_8), &
                                     (0.0005870177219337786_8,123.06502346180824_8), &
                                     (0.000032694803076194855_8,1710.380367052896_8), &
                                     (1.4887334377287222d-6,30663.78158349032_8), &
                                     (5.733255383704708d-8,672892.8144697341_8), &
                                     (1.9129620345019236d-9,1.7464549394629598d7), &
                                     (5.6307198010606435d-11,5.2326358902441823d8), &
                                     (1.4826978586953883d-12,1.777349747743559d10), &
                                     (3.532061981830939d-14,6.748696405535281d11), &
                                     (7.681736736061813d-16,2.832675140577074d13), &
                                     (1.5369167574956977d-17,1.3023556950249008d15), &
                                     (2.847051416677187d-19,6.508945799983926d16), &
                                     (4.910075099825432d-21,3.513528376296296d18), &
                                     (7.921412215629845d-23,2.037195563671853d20), &
                                     (1.2004738650741575d-24,1.2627098966389193d22), &
                                     (1.7152938595909842d-26,8.331848122253196d23), &
                                     (2.3183663953195413d-28,5.831030975680598d25), &
                                     (2.9727294547738398d-30,4.314129737191417d27), &
                                     (3.625794040507704d-32,3.364438091911737d29), &
                                     (4.216584424744342d-34,2.758407822393905d31), &
                                     (4.685647724476791d-36,2.371894283449567d33), &
                                     (4.985272847540434d-38,2.1344290143223713d35), &
                                     (5.087522112610439d-40,2.0061260840346838d37), &
                                     (4.988228177000061d-42,1.9657901194525576d39), &
                                     (4.706279296883083d-44,2.0049053092332057d41), &
                                     (4.278776960091815d-46,2.1250030487752523d43), &
                                     (3.7535921782660104d-48,2.3373028631218544d45), &
                                     (3.1812313151664884d-50,2.6643127636540365d47), &
                                     (2.607736304099903d-52,3.1436553308254513d49), &
                                     (2.069758353509157d-54,3.834993072330685d51), &
                                     (1.5922132167825223d-56,4.8317769056035805d53), &
                                     (1.1882830810796369d-58,6.280926477977422d55), &
                                     (8.611186421470097d-61,8.415958302799185d57), &
                                     (6.064508322338716d-63,1.1613394365215078d60), &
                                     (4.15396251038509d-65,1.649017840277513d62), &
                                     (2.7694282338524123d-67,2.407449912861517d64), &
                                     (1.7984039359308915d-69,3.611009967508248d66), &
                                     (1.138274811658968d-71,5.560714604971416d68), &
                                     (7.0266490269757956d-74,8.785567974858086d70) /

            zeta = lib_math_bessel_riccati_zeta_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_real:"
            do i=1, n
!                buffer = abs(zeta(i) - ground_truth_zeta(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
                buffer = 1 - abs(real(zeta(i)) / real(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(zeta(i)) / aimag(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_real

        function test_lib_math_bessel_riccati_zeta_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: z = cmplx(0.5, 2)
            complex(kind=8), dimension(n) :: zeta
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_zeta
            double precision :: ground_truth_e = 6* 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_zeta / (3.542502200006498_8,6.484506781251242_8), &
                                     (-3.016209213602835_8,2.638325491915351_8), &
                                     (-0.882351816221121_8,-1.295155364900645_8), &
                                     (-0.5502456545287542_8,-1.3240596683954458_8), &
                                     (-3.932399983399321_8,2.017327205846175_8), &
                                     (4.930502073336801_8,20.114923698394403_8), &
                                     (114.43736063470033_8,-1.50884844637336_8), &
                                     (160.86074134309743_8,-722.509956969838_8), &
                                     (-4930.635748639856_8,-2408.996308628145_8), &
                                     (-29294.10270764797_8,35349.60332883239_8), &
                                     (255516.74181286927_8,343349.4985473411_8), &
                                     (4.05367168577199d6,-1.712180996597874d6), &
                                     (-7.818599437606826d6,-4.885134455887371d7), &
                                     (-6.017712530772474d8,-4.998472138118568d7), &
                                     (-2.5387900761216555d9,7.536111092211899d9), &
                                     (9.478576825178926d10,6.04084400748232d10), &
                                     (1.2294806000859502d12,-1.1699800664924888d12), &
                                     (-1.3490610353448512d13,-2.3695794487792188d13), &
                                     (-4.4706272656027456d14,1.2579734976179428d14), &
                                     (2.578067141199645d14,8.355435261794503d15), &
                                     (1.5497675245486864d17,3.3479394156975972d16), &
                                     (1.393232486508666d18,-2.8370062872809047d18), &
                                     (-5.0514516104507d19,-4.2577862692106306d19), &
                                     (-1.1704718806373267d21,8.471439566009837d20), &
                                     (1.2315324569166937d22,3.061486957093489d22), &
                                     (7.781075706797452d23,-1.0833832120141384d23), &
                                     (2.05621039067537d24,-1.9355226493093305d25), &
                                     (-4.7069915060950156d26,-1.7186149779140146d26), &
                                     (-7.495936539523721d27,1.1090111785971126d28), &
                                     (2.4767977026220532d29,2.7560832065211643d29), &
                                     (9.378868891759476d30,-4.974800095716042d30), &
                                     (-7.574652929466055d31,-3.0520582190041577d32), &
                                     (-9.619249276476275d33,-1.1477129018603129d31), &
                                     (-7.38342836490342d34,2.944532997627971d35), &
                                     (8.711570700093456d36,4.648947842075908d36), &
                                     (2.217453616624019d38,-2.4542576060829878d38), &
                                     (-6.356593374785546d39,-9.463579737878706d39), &
                                     (-3.799154629875362d41,1.373376545338122d41), &
                                     (1.5013726682075457d42,1.4630047460478677d43), &
                                     (5.441046582017421d44,7.799100030416439d43), &
                                     (7.954901461926326d45,-1.9517663337780652d46) /

            zeta = lib_math_bessel_riccati_zeta_cmplx(z, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_cmplx:"
            do i=1, n
!                buffer = abs(zeta(i) - ground_truth_zeta(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
                buffer = 1 - abs(real(zeta(i)) / real(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(zeta(i)) / aimag(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_cmplx

        function test_lib_math_bessel_spherical_first_kind_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: y
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(0,5):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.1161107492591575_8, &
                                  -0.2472559984565608_8, &
                                  -0.09110201506935516_8, &
                                  0.04703982781631992_8, &
                                  0.07312752589258929_8, &
                                  0.04724475601390756_8, &
                                  0.02120674456246912_8, &
                                  0.007489151759344721_8, &
                                  0.002202722852205434_8, &
                                  0.0005578549358103171_8, &
                                  0.0001243811142824601_8, &
                                  0.00002480921535126831_8, &
                                  4.48114566234554d-6, &
                                  7.400582707105952d-7, &
                                  1.126278428735349d-7, &
                                  1.589853057070848d-8, &
                                  2.093112500620615d-9, &
                                  2.582298318736067d-10, &
                                  2.997661261382462d-11, &
                                  3.286142798409364d-12, &
                                  3.412743485204051d-13, &
                                  3.367207977221623d-14, &
                                  3.164424513582251d-15, &
                                  2.839094993401503d-16, &
                                  2.436884866659109d-17, &
                                  2.004885455545229d-18, &
                                  1.583811340459344d-19, &
                                  1.203304022776971d-20, &
                                  8.805419916593619d-22, &
                                  6.214761010118569d-23, &
                                  4.235967847553385d-24, &
                                  2.79158423842796d-25, &
                                  1.780739245143344d-26, &
                                  1.100657214574314d-27, &
                                  6.598247175888657d-29, &
                                  3.83996167132409d-30, &
                                  2.171305012908171d-31, &
                                  1.193886362335331d-32, &
                                  6.388316980577045d-34, &
                                  3.328943597676657d-35, &
                                  1.690525218810726d-36  /

            ground_truth_y_n = lib_math_bessel_spherical_first_kind_real(x, fnu, n)

            y = lib_math_bessel_spherical_first_kind_derivative_real(x, fnu, n, y_n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_real:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_real

        function test_lib_math_bessel_spherical_first_kind_derivative_coeff_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            ! auxiliary
            double precision :: a = 1.5
            double precision :: x = 4
            double precision, dimension(n) :: y_d
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y_d
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with mathematica
            !
            ! source code:
            !  >>> dn2[n1_, x1_, t_] := FullSimplify[D[SphericalBesselJ[n1, t y], y] /. y -> x1];
            !  >>> n = {1, 2, 3, 4, 5};
            !  >>> a = 1.5;
            !  >>> x = 4;
            !  >>> N[jdn2[n, x, a], 16]
            data ground_truth_y_d / 0.25168488408754675_8, &
                                    0.01404108681278443_8, &
                                    -0.22369060034052574_8, &
                                    -0.19267373050016112_8, &
                                    -0.04096192445311332_8, &
                                    0.05743394313180334_8, &
                                    0.07361053530996367_8, &
                                    0.05124934489528964_8, &
                                    0.026561338909918065_8, &
                                    0.01125078294800103_8, &
                                    0.004074194866196927_8, &
                                    0.0012957355341196735_8, &
                                    0.00036854587434187995_8, &
                                    0.00009499655184698529_8, &
                                    0.000022416451337026936_8, &
                                    4.881707322527168d-6, &
                                    9.876189688047362d-7, &
                                    1.8664551404040657d-7, &
                                    3.3105196704535166d-8, &
                                    5.5333177921378685d-9, &
                                    8.746266919019054d-10, &
                                    1.3114826465814094d-10, &
                                    1.870724982108646d-11, &
                                    2.5447442062829254d-12, &
                                    3.3085522154100107d-13, &
                                    4.1197752328392323d-14, &
                                    4.9221446387635015d-15, &
                                    5.652181689834273d-16, &
                                    6.247909976337875d-17, &
                                    6.657848323212348d-18, &
                                    6.848424367908085d-19, &
                                    6.808334622863263d-20, &
                                    6.549127770211743d-21, &
                                    6.102166767628009d-22, &
                                    5.512879447508527d-23, &
                                    4.8336525574023176d-24, &
                                    4.1167925372345107d-25, &
                                    3.4087162624577767d-26, &
                                    2.7460742427712334d-27, &
                                    2.1540001137316756d-28, &
                                    1.6462524570794128d-29 /

            ground_truth_y_n = lib_math_bessel_spherical_first_kind_real(a*x, fnu, n)

            y_d = lib_math_bessel_spherical_first_kind_derivative_coeff_real(a, x, fnu, n, y_n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_coeff_real:"
            do i=1, n
                buffer = y_d(i) - ground_truth_y_d(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_coeff_real

        function test_lib_math_bessel_spherical_first_kind_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: y
            complex(kind=8), dimension(n) :: y_n
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / (-0.0371012027247164_8,0.7521075843567761_8), &
                                  (-0.6709341985234971_8,0.1188523491021362_8), &
                                  (-0.2428898097811523_8,-0.4073740884193695_8), &
                                  (0.1884825896869382_8,-0.2432051985674825_8), &
                                  (0.1966017023583236_8,0.0179372642470247_8), &
                                  (0.06895184186820641_8,0.08302030524938263_8), &
                                  (0.00150286328364212_8,0.04813467452814326_8), &
                                  (-0.00955034182416434_8,0.01510751091592448_8), &
                                  (-0.005095215564063652_8,0.002464111624569094_8), &
                                  (-0.001557161824030634_8,-0.00012870280736273_8), &
                                  (-0.0003157118348054715_8,-0.0002158786595714773_8), &
                                  (-0.00003769211596566368_8,-0.00007520886467440467_8), &
                                  (5.7209387829845d-7,-0.00001679663859320406_8), &
                                  (1.546746642267723d-6,-2.659173392929424d-6), &
                                  (4.415052391556155d-7,-2.741445505127181d-7), &
                                  (8.132134146214074d-8,-5.77396157401629d-9), &
                                  (1.094394791427324d-8,4.770836805887377d-9), &
                                  (1.030870596972466d-9,1.274779272135755d-9), &
                                  (4.25618016122201d-11,2.076371190662988d-10), &
                                  (-7.00711553164415d-12,2.492285273739106d-11), &
                                  (-2.040647697000161d-12,2.194804873835421d-12), &
                                  (-3.099383605805499d-13,1.124198621271327d-13), &
                                  (-3.428760333992894d-14,-4.30496694999901d-15), &
                                  (-2.860571972125186d-15,-1.944159526722931d-15), &
                                  (-1.595507370738868d-16,-2.902878317691357d-16), &
                                  (-9.695296997007d-19,-3.039766596929408d-17), &
                                  (1.137920864918278d-18,-2.428215601922184d-18), &
                                  (1.792583488139373d-19,-1.39969210833905d-19), &
                                  (1.823279478187804d-20,-3.57753007915267d-21), &
                                  (1.40919186367387d-21,3.97664352011941d-22), &
                                  (8.203520354433243d-23,7.543118331802168d-23), &
                                  (2.883030130589105d-24,7.67834708128827d-24), &
                                  (-6.33616342233908d-26,5.808654457896537d-25), &
                                  (-2.197880412671735d-26,3.382418652347734d-26), &
                                  (-2.334274405958654d-27,1.358992751189042d-27), &
                                  (-1.752732421572115d-28,1.05454811509191d-29), &
                                  (-1.0186055502200211d-29,-4.388929069400218d-30), &
                                  (-4.374496375137839d-31,-5.223101453837363d-31), &
                                  (-9.30984171979371d-33,-3.964958459451358d-32), &
                                  (5.53962929919473d-34,-2.30553439140007d-33), &
                                  (8.688566712540365d-35,-1.0271694104114872d-34) /

            y = lib_math_bessel_spherical_first_kind_derivative_cmplx(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_cmplx

        function test_lib_math_bessel_spherical_first_kind_derivative_coeff_c() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-5.0_8)

            ! auxiliary
            real(kind=8) :: a = 1.5
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: y_d
            complex(kind=8), dimension(n) :: y_n
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y_d
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y_d / (2.0130070590108713_8,-0.6742293589646774_8), &
                                    (0.8461790601590833_8,1.756883810435237_8), &
                                    (-1.2785972896445386_8,1.0734175202677885_8), &
                                    (-1.1332150020510934_8,-0.607965638904099_8), &
                                    (0.03027219668789972_8,-0.888808920365951_8), &
                                    (0.5082519115384914_8,-0.27754598415735443_8), &
                                    (0.3237504750376738_8,0.14432830477409223_8), &
                                    (0.06783467143800523_8,0.18338987733011186_8), &
                                    (-0.03256390381766831_8,0.08915632912535072_8), &
                                    (-0.033906364977412506_8,0.02216971246550141_8), &
                                    (-0.015331614759037363_8,-0.00030908304433859457_8), &
                                    (-0.004322986159924545_8,-0.0028962537204345424_8), &
                                    (-0.0006680876306840445_8,-0.0014525932232997692_8), &
                                    (0.00004622467764871443_8,-0.000446420509218682_8), &
                                    (0.00006794441436436117_8,-0.00009397770355430132_8), &
                                    (0.000025224705421438888_8,-0.00001159448650579087_8), &
                                    (6.173987753692962d-6,4.7097025681319643d-7), &
                                    (1.0832156218721137d-6,7.058472653791909d-7), &
                                    (1.2089127437225773d-7,2.231567905648795d-7), &
                                    (1.619282852965994d-10,4.701200926001258d-8), &
                                    (-3.853951621082915d-9,7.28806016629315d-9), &
                                    (-1.1351572875813267d-9,7.719051688614596d-10), &
                                    (-2.1622746652111304d-10,2.4375743622698155d-11), &
                                    (-3.0677920168656045d-11,-1.1923512678517219d-11), &
                                    (-3.1350981047479056d-12,-3.5829017710299307d-12), &
                                    (-1.562457850847462d-13,-6.410749600733837d-13), &
                                    (2.0815784721487364d-14,-8.527849809335664d-14), &
                                    (7.37599489439778d-15,-8.465418094745307d-15), &
                                    (1.2853889412886149d-15,-5.10057126670183d-16), &
                                    (1.6345305600578746d-16,1.601088959718294d-17), &
                                    (1.5832336361567132d-17,1.019782462107601d-17), &
                                    (1.0457762919946237d-18,1.806663959492635d-18), &
                                    (1.1533483396542314d-20,2.236387577478277d-19), &
                                    (-9.512463165338626d-21,2.124009501965607d-20), &
                                    (-1.8272255375841344d-21,1.4705999436073604d-21), &
                                    (-2.2459190931098317d-22,4.695997327306761d-23), &
                                    (-2.1044856391269004d-23,-5.71083696348631d-24), &
                                    (-1.4938333369062834d-24,-1.3521990464578699d-24), &
                                    (-6.463646323193738d-26,-1.692461729617684d-25), &
                                    (1.6647540750813945d-27,-1.5774409468746083d-26), &
                                    (7.357137931665447d-28,-1.1358244789533994d-27) /

            y_d = lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx(a, x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_first_kind_cmplx(a * x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx:"
            do i=1, n
                buffer = abs(y_d(i) - ground_truth_y_d(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_coeff_c

        function test_lib_math_bessel_spherical_second_kind_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: y
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_Y(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.2300533501309578_8, &
                                    0.04838423015042407_8, &
                                    0.2232065195942213_8, &
                                    0.227771073285379_8, &
                                    0.2710487187377822_8, &
                                    0.6024493519630117_8, &
                                    1.84136286984519_8, &
                                    6.544627609742551_8, &
                                    26.43954310954968_8, &
                                    120.1918931195346_8, &
                                    607.9866463873159_8, &
                                    3387.440883381229_8, &
                                    20608.20379774237_8, &
                                    135903.0296163689_8, &
                                    965498.587536345_8, &
                                    7.350313987097137d6, &
                                    5.968971659797807d7, &
                                    5.149812373496762d8, &
                                    4.703779320314724d9, &
                                    4.534225758089556d10, &
                                    4.599824146863092d11, &
                                    4.898518074414196d12, &
                                    5.463635608849164d13, &
                                    6.369283804533598d14, &
                                    7.745830859572537d15, &
                                    9.809776146333946d16, &
                                    1.291722466011729d18, &
                                    1.765862686388832d19, &
                                    2.502811305756681d20, &
                                    3.673083628895114d21, &
                                    5.575078122090605d22, &
                                    8.741953276043991d23, &
                                    1.414666566358669d25, &
                                    2.36029538916938d26, &
                                    4.056478714049619d27, &
                                    7.175127873479323d28, &
                                    1.305136240954151d30, &
                                    2.439469264843751d31, &
                                    4.682022535205213d32, &
                                    9.220873364381276d33, &
                                    1.862206939549379d35 /

            y = lib_math_bessel_spherical_second_kind_derivative_real(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_second_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_derivative_real:"
            do i=1, n
                buffer = 1 - abs(y(i) / ground_truth_y(i))
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = 1 - abs(y_n(i) / ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_derivative_real

        function test_lib_math_bessel_spherical_second_kind_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: y
            complex(kind=8), dimension(n+1) :: y_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_y
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_Y(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / (-0.0925935033610706_8,-0.644258003230032_8), &
                                    (0.4512742917535701_8,-0.2583992418287995_8), &
                                    (0.2637440265921878_8,0.1243919192227003_8), &
                                    (-0.0670923351374301_8,0.08950059905227736_8), &
                                    (-0.3047213282882742_8,-0.0586613842326251_8), &
                                    (-0.72747325548997_8,0.0237735452404915_8), &
                                    (-1.830519915862459_8,1.127264704638521_8), &
                                    (-3.724380823782526_8,6.825738734951031_8), &
                                    (0.14386984853883_8,31.8934401958393_8), &
                                    (69.9242398279025_8,127.6044312913219_8), &
                                    (611.8837411905042_8,399.0219813219127_8), &
                                    (3984.284201632267_8,360.707281503365_8), &
                                    (21910.65483258843_8,-9109.74474052318_8), &
                                    (97041.0896414393_8,-116360.1058745882_8), &
                                    (221696.9484172507_8,-1.012096395516843d6), &
                                    (-1.956665755212431d6,-7.295949039841005d6), &
                                    (-3.935608639002741d7,-4.326172186791142d7), &
                                    (-4.488706239772405d8,-1.678213680849475d8), &
                                    (-4.114936847565165d9,4.7936033369543d8), &
                                    (-3.133592088707219d10,2.092574993067198d10), &
                                    (-1.756035904023185d11,3.139164380979988d11), &
                                    (-1.39301140564529d11,3.592802315777198d12), &
                                    (1.572159723181733d13,3.410928660972549d13), &
                                    (3.210561819792464d14,2.537330025101746d14), &
                                    (4.550761659523629d15,9.18670937570185d14), &
                                    (5.277259869276212d16,-1.460395876900965d16), &
                                    (4.956351064190754d17,-4.514334327611531d17), &
                                    (3.026254092661271d18,-7.952888414618691d18), &
                                    (-1.166590463554708d19,-1.1121382438985242d20), &
                                    (-8.23093839066282d20,-1.277532456695498d21), &
                                    (-1.839904034830848d22,-1.08019696553135d22), &
                                    (-3.084803722428571d23,-1.9626423559648d22), &
                                    (-4.245939803984751d24,1.813132137466256d24), &
                                    (-4.577375100820782d25,5.431335801963718d25), &
                                    (-2.60258394019907d26,1.0940929254729155d27), &
                                    (4.22867264599862d27,1.78155865942922d28), &
                                    (1.972215157269966d29,2.343999457416233d29), &
                                    (4.836428125923716d30,2.073579269905943d30), &
                                    (9.257971994931163d31,-5.66557135191576d30), &
                                    (1.448860512137253d33,-8.43384996133509d32) /

            y = lib_math_bessel_spherical_second_kind_derivative_cmplx(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(y(i)) / real(ground_truth_y(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(y(i)) / aimag(ground_truth_y(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(y_n(i) / ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_derivative_cmplx

        function test_lib_math_bessel_spherical_third_kind_1_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 0.5
            complex(kind=8), dimension(n) :: h_1
            complex(kind=8), dimension(n) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1
            complex(kind=8), dimension(n) :: ground_truth_h_1_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel1(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_h_1 / (0.30870295466413966_8,16.121560175298843_8), &
                                    (0.0643103909881061_8,145.89035562426193_8), &
                                    (0.006978823057052953_8,1943.980452564093_8), &
                                    (0.0005201393823436602_8,33961.4772941343_8), &
                                    (0.000029660003646900365_8,732509.9972696619_8), &
                                    (1.3721553680201262d-6,1.8779671241985574d7), &
                                    (5.345032257003262d-8,5.575197949992077d8), &
                                    (1.7988649406220154d-9,1.8802560106089794d10), &
                                    (5.3306481673397344d-11,7.098933719193749d11), &
                                    (1.4112884453851635d-12,2.9658717189400355d13), &
                                    (3.376890330352207d-14,1.3583343281958888d15), &
                                    (7.371506333145997d-16,6.76658426384833d16), &
                                    (1.4794847216521714d-17,3.642404936600949d18), &
                                    (2.7480577734591153d-19,2.1068152366177806d20), &
                                    (4.750446381647761d-21,1.3031024550747268d22), &
                                    (7.679602148755423d-23,8.582352906017307d23), &
                                    (1.1659361512428058d-24,5.996405228229023d25), &
                                    (1.668629258739116d-26,4.4299171718928035d27), &
                                    (2.2585492268200136d-28,3.450137583557998d29), &
                                    (2.8997919155212113d-30,2.825265171258421d31), &
                                    (3.540993787240369d-32,2.426725996088254d33), &
                                    (4.122372942970052d-34,2.181591059209123d35), &
                                    (4.585433515314721d-36,2.0485774748927863d37), &
                                    (4.883023582470525d-38,2.005699198231819d39), &
                                    (4.987286921140749d-40,2.0440204990138535d41), &
                                    (4.89367471336639d-42,2.1649045759479713d43), &
                                    (4.620328398463419d-44,2.379602433566436d45), &
                                    (4.203386993394978d-46,2.710846320611596d47), &
                                    (3.6897067783322705d-48,3.1967078558122196d49), &
                                    (3.1288696132491395d-50,3.8975997476708284d51), &
                                    (2.5661819157080416d-52,4.908162401517112d53), &
                                    (2.0377952608653985d-54,6.377178516782261d55), &
                                    (1.5683614432967147d-56,8.541093654668174d57), &
                                    (1.1710000631534732d-58,1.1781085438623265d60), &
                                    (8.489480858772285d-61,1.6721604694249153d62), &
                                    (5.98115212930763d-63,2.440314135723415d64), &
                                    (4.098394105314449d-65,3.658994063981451d66), &
                                    (2.7333463276526288d-67,5.632694059330295d68), &
                                    (1.7755681732074424d-69,8.896421165960764d70), &
                                    (1.1241791828939346d-71,1.4407219335846267d73) /

            h_1 = lib_math_bessel_spherical_third_kind_1_derivative_real(x, fnu, n, h_1_n)
            ground_truth_h_1_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_derivative_real:"
            do i=1, n
                buffer = 1 - abs(h_1(i) / ground_truth_h_1(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_1_n(i) - ground_truth_h_1_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_derivative_real

        function test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40
            double precision, parameter :: ground_truth_e = 2* 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(0.5, 2, kind=8)
            complex(kind=8), dimension(n) :: h_1
            complex(kind=8), dimension(n) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1
            complex(kind=8), dimension(n) :: ground_truth_h_1_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> f(n,x) = derivative(spherical_hankel1(n, x), x)
            !  >>>
            !  >>> x=4.0 + I*2
            !  >>> diff = 0
            !  >>> for n in range(1,6):
            !  >>>     value = numerical_approx(f(n,x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_h_1 / (-0.09338594701179898_8,-0.12691536888887217_8), &
                                    (-0.34939687137478814_8,0.17008874206691427_8), &
                                    (0.31535384630518887_8,1.2877908427890217_8), &
                                    (5.814370082548807_8,-0.08206209692584164_8), &
                                    (6.744979480000595_8,-30.438887021878138_8), &
                                    (-177.52197553129855_8,-86.55128863186137_8), &
                                    (-919.9355586016592_8,1111.5549749967086_8), &
                                    (7133.3717399221605_8,9575.513863679997_8), &
                                    (101673.35571624429_8,-42990.69107536378_8), &
                                    (-178524.47290051883_8,-1.1132818222871753d6), &
                                    (-1.256575853797051d7,-1.0407031690311907d6), &
                                    (-4.88879584669103d7,1.4521141738794243d8), &
                                    (1.6955166527577977d9,1.0802017423193443d9), &
                                    (2.0517309818282627d10,-1.9529488430255642d10), &
                                    (-2.1108584238533273d11,-3.7067100015306006d11), &
                                    (-6.581237481765829d12,1.8525231795684558d12), &
                                    (3.593067492720328d12,1.1615584184870042d14), &
                                    (2.0408913693295595d15,4.407477267291251d14), &
                                    (1.7426305203252506d16,-3.548995242870003d16), &
                                    (-6.017912105312503d17,-5.0718836461329536d17), &
                                    (-1.3308448842996597d19,9.633049022980035d18), &
                                    (1.339487482123132d20,3.3294834671265323d20), &
                                    (8.109322317166058d21,-1.1293675433724366d21), &
                                    (2.056572065343048d22,-1.9364204526652986d23), &
                                    (-4.527915698259452d24,-1.6530925230288849d24), &
                                    (-6.943148155240989d25,1.0272805285052258d26), &
                                    (2.2122805212894362d27,2.461631410206111d27), &
                                    (8.087865138744946d28,-4.290220481823455d28), &
                                    (-6.3146364016085496d29,-2.544166324552342d30), &
                                    (-7.759723763813355d31,-9.133562885833038d28), &
                                    (-5.769540972649826d32,2.3010529184050733d33), &
                                    (6.60141075880758d34,3.5227430523793565d34), &
                                    (1.630850325781663d36,-1.8050548739971377d36), &
                                    (-4.541526702489605d37,-6.761177345063087d37), &
                                    (-2.6388543519788082d39,9.539643448699718d38), &
                                    (1.0147395998345191d40,9.88714720376597d40), &
                                    (3.580321849632509d42,5.1316476009643184d41), &
                                    (5.100112938337967d43,-1.2513620114416501d44), &
                                    (-4.180150285190999d45,-3.0585660835951755d45), &
                                    (-1.5641985546687648d47,1.3026978419514984d47) /

            h_1 = lib_math_bessel_spherical_third_kind_1_derivative_cmplx(x, fnu, n, h_1_n)
            ground_truth_h_1_n = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(h_1(i)) / real(ground_truth_h_1(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(h_1(i)) / aimag(ground_truth_h_1(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(h_1_n(i) / ground_truth_h_1_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx

        function test_lib_math_bessel_spherical_third_kind_2_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 0.5
            complex(kind=8), dimension(n) :: h_2
            complex(kind=8), dimension(n) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2
            complex(kind=8), dimension(n) :: ground_truth_h_2_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel2(n, x))
            !  >>>     x=0.5
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_h_2 / (0.30870295466413966_8,-16.121560175298843_8), &
                                    (0.0643103909881061_8,-145.89035562426193_8), &
                                    (0.006978823057052953_8,-1943.980452564093_8), &
                                    (0.0005201393823436602_8,-33961.4772941343_8), &
                                    (0.000029660003646900365_8,-732509.9972696619_8), &
                                    (1.3721553680201262d-6,-1.8779671241985574d7), &
                                    (5.345032257003262d-8,-5.575197949992077d8), &
                                    (1.7988649406220154d-9,-1.8802560106089794d10), &
                                    (5.3306481673397344d-11,-7.098933719193749d11), &
                                    (1.4112884453851635d-12,-2.9658717189400355d13), &
                                    (3.376890330352207d-14,-1.3583343281958888d15), &
                                    (7.371506333145997d-16,-6.76658426384833d16), &
                                    (1.4794847216521714d-17,-3.642404936600949d18), &
                                    (2.7480577734591153d-19,-2.1068152366177806d20), &
                                    (4.750446381647761d-21,-1.3031024550747268d22), &
                                    (7.679602148755423d-23,-8.582352906017307d23), &
                                    (1.1659361512428058d-24,-5.996405228229023d25), &
                                    (1.668629258739116d-26,-4.4299171718928035d27), &
                                    (2.2585492268200136d-28,-3.450137583557998d29), &
                                    (2.8997919155212113d-30,-2.825265171258421d31), &
                                    (3.540993787240369d-32,-2.426725996088254d33), &
                                    (4.122372942970052d-34,-2.181591059209123d35), &
                                    (4.585433515314721d-36,-2.0485774748927863d37), &
                                    (4.883023582470525d-38,-2.005699198231819d39), &
                                    (4.987286921140749d-40,-2.0440204990138535d41), &
                                    (4.89367471336639d-42,-2.1649045759479713d43), &
                                    (4.620328398463419d-44,-2.379602433566436d45), &
                                    (4.203386993394978d-46,-2.710846320611596d47), &
                                    (3.6897067783322705d-48,-3.1967078558122196d49), &
                                    (3.1288696132491395d-50,-3.8975997476708284d51), &
                                    (2.5661819157080416d-52,-4.908162401517112d53), &
                                    (2.0377952608653985d-54,-6.377178516782261d55), &
                                    (1.5683614432967147d-56,-8.541093654668174d57), &
                                    (1.1710000631534732d-58,-1.1781085438623265d60), &
                                    (8.489480858772285d-61,-1.6721604694249153d62), &
                                    (5.98115212930763d-63,-2.440314135723415d64), &
                                    (4.098394105314449d-65,-3.658994063981451d66), &
                                    (2.7333463276526288d-67,-5.632694059330295d68), &
                                    (1.7755681732074424d-69,-8.896421165960764d70), &
                                    (1.1241791828939346d-71,-1.4407219335846267d73) /

            h_2 = lib_math_bessel_spherical_third_kind_2_derivative_real(x, fnu, n, h_2_n)
            ground_truth_h_2_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_derivative_real:"
            do i=1, n
                buffer = 1 - abs(h_2(i) / ground_truth_h_2(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_2_n(i) - ground_truth_h_2_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_derivative_real

        function test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40
            double precision, parameter :: ground_truth_e = 10.0_8**(-12.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(0.5, 2, kind=8)
            complex(kind=8), dimension(n) :: h_2
            complex(kind=8), dimension(n) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2
            complex(kind=8), dimension(n) :: ground_truth_h_2_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel2(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_h_2 / (1.6316265049670142_8,-0.4766276791663696_8), &
                                    (0.7673797292340883_8,0.6300130786704627_8), &
                                    (-0.5838076502086155_8,-1.0804193757617195_8), &
                                    (-5.886757145712604_8,0.020171434949572634_8), &
                                    (-6.735050875885715_8,30.42015381399214_8), &
                                    (177.5257446783282_8,86.55228912964594_8), &
                                    (919.9355422208935_8,-1111.5543648968348_8), &
                                    (-7133.371821186327_8,-9575.513846460597_8), &
                                    (-101673.35572056167_8,42990.69106634515_8), &
                                    (178524.4729013523_8,1.11328182228649d6), &
                                    (1.2565758537970562d7,1.0407031690312602d6), &
                                    (4.8887958466910295d7,-1.4521141738794205d8), &
                                    (-1.695516652757793d9,-1.0802017423193424d9), &
                                    (-2.051730981828259d10,1.9529488430255566d10), &
                                    (2.1108584238533194d11,3.7067100015305927d11), &
                                    (6.581237481765814d12,-1.8525231795684456d12), &
                                    (-3.593067492720249d12,-1.1615584184870012d14), &
                                    (-2.0408913693295538d15,-4.407477267291263d14), &
                                    (-1.7426305203252488d16,3.5489952428699932d16), &
                                    (6.017912105312481d17,5.071883646132948d17), &
                                    (1.3308448842996572d19,-9.633049022979998d18), &
                                    (-1.3394874821231246d20,-3.329483467126526d20), &
                                    (-8.109322317166042d21,1.1293675433724277d21), &
                                    (-2.0565720653430608d22,1.9364204526652936d23), &
                                    (4.5279156982594386d24,1.6530925230288838d24), &
                                    (6.943148155240982d25,-1.0272805285052222d26), &
                                    (-2.2122805212894285d27,-2.4616314102061064d27), &
                                    (-8.087865138744927d28,4.290220481823436d28), &
                                    (6.314636401608519d29,2.544166324552336d30), &
                                    (7.759723763813338d31,9.133562885839941d28), &
                                    (5.769540972649818d32,-2.301052918405067d33), &
                                    (-6.60141075880756d34,-3.522743052379351d34), &
                                    (-1.63085032578166d36,1.8050548739971315d36), &
                                    (4.541526702489588d37,6.761177345063074d37), &
                                    (2.6388543519788022d39,-9.53964344869968d38), &
                                    (-1.014739599834507d40,-9.887147203765945d40), &
                                    (-3.5803218496325005d42,-5.131647600964331d41), &
                                    (-5.100112938337964d43,1.2513620114416464d44), &
                                    (4.180150285190985d45,3.05856608359517d45), &
                                    (1.5641985546687618d47,-1.3026978419514933d47) /

            h_2 = lib_math_bessel_spherical_third_kind_2_derivative_cmplx(x, fnu, n, h_2_n)
            ground_truth_h_2_n = lib_math_bessel_spherical_third_kind_2_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(h_2(i)) / real(ground_truth_h_2(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(h_2(i)) / aimag(ground_truth_h_2(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(h_2_n(i) / ground_truth_h_2_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx

        function test_lib_math_bessel_riccati_s_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: s
            double precision, dimension(n) :: s_n

            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_s
            double precision, dimension(n) :: ground_truth_s_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            data ground_truth_s / -0.8729132445670857_8, &
                                    -0.08812437450607047_8, &
                                    0.4174031692203099_8, &
                                    0.4174031692203099_8, &
                                    0.2407445638129937_8, &
                                    0.1022891469326733_8, &
                                    0.03494311549910491_8, &
                                    0.01004812945749755_8, &
                                    0.002503172988387468_8, &
                                    0.0005511143228984732_8, &
                                    0.0001088304115441947_8, &
                                    0.00001949763018069782_8, &
                                    3.198230014443829d-6, &
                                    4.839431284882498d-7, &
                                    6.797742888868432d-8, &
                                    8.911319203712997d-9, &
                                    1.095283631795967d-9, &
                                    1.267249118632387d-10, &
                                    1.385103491554362d-11, &
                                    1.434657274726526d-12, &
                                    1.412133737929903d-13, &
                                    1.32421554787713d-14, &
                                    1.185729318210958d-15, &
                                    1.015909902157669d-16, &
                                    8.344266452005118d-18, &
                                    6.581680170397773d-19, &
                                    4.993365242940314d-20, &
                                    3.649195487556749d-21, &
                                    2.572410951804915d-22, &
                                    1.751350841320649d-23, &
                                    1.152943520264954d-24, &
                                    7.347228888965333d-26, &
                                    4.536989078220226d-27, &
                                    2.717444868257271d-28, &
                                    1.580147314294711d-29, &
                                    8.927918382832613d-31, &
                                    4.905345128253322d-32, &
                                    2.622933468237204d-33, &
                                    1.365894679100282d-34, &
                                    6.931973927355216d-36 /

            s = lib_math_bessel_riccati_s_derivative_real(x, fnu, n, s_n)
            ground_truth_s_n = lib_math_bessel_riccati_s_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_derivative_real:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(s_n(i) - ground_truth_s_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_s_derivative_real

        function test_lib_math_bessel_riccati_s_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: s
            complex(kind=8), dimension(n) :: s_n

            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_s
            complex(kind=8), dimension(n) :: ground_truth_s_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>> S_(n,x) = derivative(S(n,x), x)
            !  >>> # WolframAlpha: S_w = x SphericalBesselJ[-1 + n, x] + SphericalBesselJ[n, x] - x SphericalBesselJ[1 + n, x])/2
            !  >>> S_w(n,x) = 1/2 * (x*spherical_bessel_J(n-1, x) + spherical_bessel_J(n,x) - x*spherical_bessel_J(n+1,x))
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S_w)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2.0
            !  >>> print("WolframAlpha:")
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S_w(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            data ground_truth_s / (-2.884340289573544_8,-1.618566584995225_8), &
                                    (0.446332618346892_8,-2.388259959485746_8), &
                                    (1.669891240667783_8,-0.418303857101622_8), &
                                    (0.8730455094616752_8,0.6858235695133681_8), &
                                    (0.0995240805571297_8,0.5797392353248442_8), &
                                    (-0.1145752489114047_8,0.2273997882298409_8), &
                                    (-0.07998685786978222_8,0.0460115335410182_8), &
                                    (-0.02867987845106044_8,-0.0007871244909737_8), &
                                    (-0.00663180671929987_8,-0.004121635923343606_8), &
                                    (-0.000906196205521605_8,-0.001658228848918162_8), &
                                    (1.858188384969d-6,-0.0004118159276574162_8), &
                                    (0.00003927753831134008_8,-0.00007157503859314495_8), &
                                    (0.00001244415857249217_8,-8.10003717268538d-6), &
                                    (2.484621029888437d-6,-2.22643015949508d-7), &
                                    (3.594933630206764d-7,1.498452078165536d-7), &
                                    (3.63405554379259d-8,4.365020715500273d-8), &
                                    (1.656918150340515d-9,7.592406777729866d-9), &
                                    (-2.600845438895835d-10,9.670673697487057d-10), &
                                    (-8.213149068837107d-11,9.016662772318956d-11), &
                                    (-1.3192961844538401d-11,4.923887009545571d-12), &
                                    (-1.535064888635538d-12,-1.79627895937315d-13), &
                                    (-1.343876604252431d-13,-8.9832938078866d-14), &
                                    (-7.87726166126415d-15,-1.40956168126403d-14), &
                                    (-5.9222264640612d-17,-1.542559677111798d-15), &
                                    (5.92755446437021d-17,-1.284869440473478d-16), &
                                    (9.776249955876346d-18,-7.721287322200135d-18), &
                                    (1.0342611969666951d-18,-2.084828545279721d-19), &
                                    (8.295922889020298d-20,2.29807263353599d-20), &
                                    (5.008105600352978d-21,4.563640016269953d-21), &
                                    (1.830910665553295d-22,4.814704333446164d-22), &
                                    (-3.95818414869941d-24,3.766216612936923d-23), &
                                    (-1.460168328215906d-24,2.265432776544602d-24), &
                                    (-1.603220758209962d-25,9.4087846958255d-26), &
                                    (-1.2411548274156607d-26,7.87720528379952d-28), &
                                    (-7.428581901508817d-28,-3.173471262110996d-28), &
                                    (-3.28519190054081d-29,-3.899136666469875d-29), &
                                    (-7.23861509733668d-31,-3.044695855116268d-30), &
                                    (4.31961918743045d-32,-1.818942390814722d-31), &
                                    (7.003697309552829d-33,-8.322158661277737d-33), &
                                    (5.66879589333865d-34,-2.429778391333543d-34) /

            s = lib_math_bessel_riccati_s_derivative_cmplx(x, fnu, n, s_n)
            ground_truth_s_n = lib_math_bessel_riccati_s_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_derivative_cmplx:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(s_n(i) - ground_truth_s_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_s_derivative_cmplx

        function test_lib_math_bessel_riccati_c_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: c
            double precision, dimension(n) :: c_n
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_c
            double precision, dimension(n) :: ground_truth_c_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> n=5
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            data ground_truth_c / -0.4235902707326541_8, &
                                    -0.9019551857592005_8, &
                                    -0.6924423272384523_8, &
                                    -0.6924423272384523_8, &
                                    -1.746996141401588_8, &
                                    -5.934500544354676_8, &
                                    -22.19072116658589_8, &
                                    -92.23491360178359_8, &
                                    -427.2815116957586_8, &
                                    -2191.411055669374_8, &
                                    -12340.43806243788_8, &
                                    -75719.72926209888_8, &
                                    -502864.6568811213_8, &
                                    -3.59366207037986d6, &
                                    -2.749659438167288d7, &
                                    -2.242660715296314d8, &
                                    -1.942264053351603d9, &
                                    -1.780007723570904d10, &
                                    -1.720975707982932d11, &
                                    -1.750547968419218d12, &
                                    -1.868718143134038d13, &
                                    -2.088857292313912d14, &
                                    -2.439948842550783d15, &
                                    -2.972674815207652d16, &
                                    -3.771057632765035d17, &
                                    -4.973259086478236d18, &
                                    -6.808418493534512d19, &
                                    -9.662512184273791d20, &
                                    -1.419794025787774d22, &
                                    -2.157456828863844d23, &
                                    -3.386599714217392d24, &
                                    -5.485855995619076d25, &
                                    -9.161466684189172d26, &
                                    -1.575912054747865d28, &
                                    -2.789808846010048d29, &
                                    -5.07858166962284d30, &
                                    -9.499596470505203d31, &
                                    -1.824523366963191d33, &
                                    -3.595657755637509d34, &
                                    -7.266244724191576d35 /

            c = lib_math_bessel_riccati_c_derivative_real(x, fnu, n, c_n)
            ground_truth_c_n = lib_math_bessel_riccati_c_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_derivative_real:"
            do i=1, n
                buffer = 1 - abs(c(i) / ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(c_n(i) - ground_truth_c_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_c_derivative_real

        function test_lib_math_bessel_riccati_c_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: c
            complex(kind=8), dimension(n) :: c_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_c
            complex(kind=8), dimension(n) :: ground_truth_c_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>> C_(n,x) = derivative(C(n,x), x)
            !  >>>
            !  >>> n=5
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            data ground_truth_c / (-1.688638603466583_8,2.753518013217212_8), &
                                    (-2.569458576765616_8,-0.437900159661906_8), &
                                    (-0.567732855301268_8,-1.461521804073214_8), &
                                    (0.7232503906203013_8,-0.4328432748969963_8), &
                                    (1.1715638715381388_8,0.656016412602718_8), &
                                    (2.628684220550715_8,1.078118501455872_8), &
                                    (8.285588449831471_8,-0.965072500729107_8), &
                                    (24.82911451738618_8,-17.98708907103369_8), &
                                    (55.02888555447147_8,-115.07248599343842_8), &
                                    (-27.562641639624_8,-587.988274649913_8), &
                                    (-1523.437921940425_8,-2564.620571124321_8), &
                                    (-14045.76340574197_8,-8589.51729575408_8), &
                                    (-98083.66711886624_8,-6453.0614664598_8), &
                                    (-577541.9170741456_8,254552.62528988_8), &
                                    (-2.716107932482874d6,3.380896131829516d6), &
                                    (-6.29506348386062d6,3.112043563573225d7), &
                                    (6.72860778544931d7,2.373854221829236d8), &
                                    (1.38375808713482d9,1.483039829741168d9), &
                                    (1.65409270093034d10,5.97200519696013d9), &
                                    (1.590987789679066d11,-2.01780928438629d10), &
                                    (1.268231829360656d12,-8.63927752354512d11), &
                                    (7.39185560833179d12,-1.347822874537009d13), &
                                    (5.0003194960432d12,-1.608121764576194d14), &
                                    (-7.46306925503746d14,-1.589658241911728d15), &
                                    (-1.573675323647912d16,-1.227210494654094d16), &
                                    (-2.313425377102014d17,-4.52625595741677d16), &
                                    (-2.781210977479913d18,7.86281919704052d17), &
                                    (-2.702858831533894d19,2.487410136143759d19), &
                                    (-1.697152723042125d20,4.525283577135613d20), &
                                    (7.15316961812617d20,6.536741537639855d21), &
                                    (5.037843945534744d22,7.74769515007118d22), &
                                    (1.1584242249816216d24,6.739873467471105d23), &
                                    (2.000084575633036d25,1.19824057987641d24), &
                                    (2.833195426419967d26,-1.221538405569863d26), &
                                    (3.138387973810009d27,-3.748819943934842d27), &
                                    (1.81950264327714d28,-7.755844386802002d28), &
                                    (-3.118314911317091d29,-1.296789001439396d30), &
                                    (-1.480980870595752d31,-1.750243748675373d31), &
                                    (-3.720931342049944d32,-1.583688147813182d32), &
                                    (-7.298988067289661d33,4.65082574288099d32) /

            c = lib_math_bessel_riccati_c_derivative_cmplx(x, fnu, n, c_n)
            ground_truth_c_n = lib_math_bessel_riccati_c_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(c(i)) / real(ground_truth_c(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(c(i)) / aimag(ground_truth_c(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(c_n(i) / ground_truth_c_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_c_derivative_cmplx

        function test_lib_math_bessel_riccati_xi_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            double precision :: x = 0.5
            complex(kind=8), dimension(n) :: xi
            complex(kind=8), dimension(n) :: xi_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            complex(kind=8), dimension(n) :: ground_truth_xi_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Xi(n,x) = x*spherical_hankel1(n,x)
            !  >>> Xi_(n,x) = derivative(Xi(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Xi_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Xi_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            data ground_truth_xi / (0.3168885079681364_8,3.5915987628795243_8), &
                                    (0.04852630210204646_8,47.88525498729233_8), &
                                    (0.004663446972394033_8,725.8601793584301_8), &
                                    (0.00032545929732421977_8,13559.97791296136_8), &
                                    (0.000017807468698907628_8,304927.4354678503_8), &
                                    (8.007427916841572d-7,8.044049992053319d6), &
                                    (3.055108535402016d-8,2.4383079871034467d8), &
                                    (1.0120468663322205d-9,8.35475287499606d9), &
                                    (2.961863655408945d-11,3.193996910048163d11), &
                                    (7.762854623292005d-13,1.3479619313593121d13), &
                                    (1.84207989989734d-14,6.225136612864029d14), &
                                    (3.993136518072138d-16,3.122820992919185d16), &
                                    (7.966833891596294d-18,1.6910235523007962d18), &
                                    (1.4722303887260664d-19,9.831370507829644d19), &
                                    (2.5336514351364773d-21,6.108073162639263d21), &
                                    (4.0798958473925434d-23,4.03863447368087d23), &
                                    (6.172739528132226d-25,2.8315656516694476d25), &
                                    (8.806819572759488d-27,2.0983379664327897d27), &
                                    (1.1887292025054836d-28,1.6387861970351707d29), &
                                    (1.5224118385707597d-30,1.3453438237909758d31), &
                                    (1.8548285821150714d-32,1.1581948415962488d33), &
                                    (2.154899425974562d-34,1.0433576439355702d35), &
                                    (2.3924222146081692d-36,9.816001571599457d36), &
                                    (2.5432622334874715d-38,9.627270774352158d38), &
                                    (2.593408024110376d-40,9.826944471178756d40), &
                                    (2.5409629426208566d-42,1.0423541817893215d43), &
                                    (2.395739738433546d-44,1.1473011558077129d45), &
                                    (2.1767653402628094d-46,1.308677103043361d47), &
                                    (1.908478015469465d-48,1.545067672633029d49), &
                                    (1.6165895327065677d-50,1.885926767218905d51), &
                                    (1.324486124924204d-52,2.3773813393119423d53), &
                                    (1.0507418947683496d-54,3.091953720279059d55), &
                                    (8.079463832699501d-57,4.144928297774538d57), &
                                    (6.027224044196768d-59,5.722223553255649d59), &
                                    (4.366030595832917d-61,8.128534459820275d61), &
                                    (3.0736553148615166d-63,1.1871767110561572d64), &
                                    (2.1045856173342726d-65,1.781348033733495d66), &
                                    (1.4026412425449323d-67,2.7441268303149825d68), &
                                    (9.105495828369005d-70,4.336996290880954d70), &
                                    (5.761428895009189d-72,7.027898308425972d72) /

            xi = lib_math_bessel_riccati_xi_derivative_real(x, fnu, n, xi_n)
            ground_truth_xi_n = lib_math_bessel_riccati_xi_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_derivative_real:"
            do i=1, n
                buffer = 1 - abs(xi(i) / ground_truth_xi(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = 1 - abs(xi_n(i) / ground_truth_xi_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_derivative_real

        function test_lib_math_bessel_riccati_xi_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            complex(kind=8) :: x = cmplx(0.5, 2)
            complex(kind=8), dimension(n) :: xi
            complex(kind=8), dimension(n) :: xi_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            complex(kind=8), dimension(n) :: ground_truth_xi_n
            double precision :: ground_truth_e = 10.0_8**(-12.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Xi(n,x) = x*spherical_hankel1(n,x)
            !  >>> Xi_(n,x) = derivative(Xi(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Xi_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Xi_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            data ground_truth_xi / (0.1360104776648551_8,-0.1844987315224899_8), &
                                    (-0.39892572561234285_8,-0.4456292977244819_8), &
                                    (-1.8829948345387955_8,1.034942441857061_8), &
                                    (2.606402649748666_8,9.460165197247633_8), &
                                    (54.21232557747803_8,-1.7735013141896265_8), &
                                    (71.58837163025785_8,-344.28792243790707_8), &
                                    (-2361.9980074843597_8,-1123.3932509000088_8), &
                                    (-13870.62175065514_8,17017.612547388533_8), &
                                    (123629.3175861611_8,163907.12893346_8), &
                                    (1.9456643803592962d6,-833840.7429793003_8), &
                                    (-3.8726432484970544d6,-2.354282539696414d7), &
                                    (-2.9095811075368094d8,-2.310233214473806d7), &
                                    (-1.2183263184242616d9,3.6538282600955253d9), &
                                    (4.606990773521424d10,2.9188549374767014d10), &
                                    (5.962203110486205d11,-5.70008992737509d11), &
                                    (-6.58973102168053d12,-1.1519989554119629d13), &
                                    (-2.1777705720932097d14,6.170327391849651d13), &
                                    (1.3234709326689375d14,4.0769744542586345d15), &
                                    (7.573075772240302d16,1.6245962615059158d16), &
                                    (6.794935558335805d17,-1.38818511917356d18), &
                                    (-2.4749170389503664d19,-2.081096891176017d19), &
                                    (-5.7294279380468715d20,4.1560929310731515d20), &
                                    (6.052442839943595d21,1.5003486570981156d22), &
                                    (3.817110917006095d23,-5.349588383161185d22), &
                                    (1.001667869178151d24,-9.503463704630947d24), &
                                    (-2.313054117029327d26,-8.425422279064145d25), &
                                    (-3.680870190146903d27,5.454090264782327d27), &
                                    (1.2190673408262197d29,1.3547398178412682d29), &
                                    (4.613763999274765d30,-2.4508797649041507d30), &
                                    (-3.7378261503573083d31,-1.5024128713037226d32), &
                                    (-4.738045259925423d33,-3.1205935954553684d30), &
                                    (-3.6310729866190325d34,1.4511657514664997d35), &
                                    (4.295655038484853d36,2.289786104640177d36), &
                                    (1.0930328255723757d38,-1.2108379313295482d38), &
                                    (-3.137949079756313d39,-4.6675021541981706d39), &
                                    (-1.8746796227305262d41,6.785255227890099d40), &
                                    (7.438978569589538d41,7.222284960903357d42), &
                                    (2.6871160668471134d44,3.841950445495348d43), &
                                    (3.926343064663275d45,-9.642709773061609d45), &
                                    (-3.3050058413677362d47,-2.416651401084738d47) /

            xi = lib_math_bessel_riccati_xi_derivative_cmplx(x, fnu, n, xi_n)
            ground_truth_xi_n = lib_math_bessel_riccati_xi_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(xi(i)) / real(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(xi(i)) / aimag(ground_truth_xi(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(xi_n(i) / ground_truth_xi_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_derivative_cmplx

        function test_lib_math_bessel_riccati_zeta_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            double precision :: x = 0.5
            complex(kind=8), dimension(n) :: zeta
            complex(kind=8), dimension(n) :: zeta_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            complex(kind=8), dimension(n) :: ground_truth_zeta_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Zeta(n,x) = x*spherical_hankel2(n,x)
            !  >>> Zeta_(n,x) = derivative(Zeta(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Zeta_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=0.5
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Zeta_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            data ground_truth_zeta / (0.3168885079681364_8,-3.5915987628795243_8), &
                                        (0.04852630210204646_8,-47.88525498729233_8), &
                                        (0.004663446972394033_8,-725.8601793584301_8), &
                                        (0.00032545929732421977_8,-13559.97791296136_8), &
                                        (0.000017807468698907628_8,-304927.4354678503_8), &
                                        (8.007427916841572d-7,-8.044049992053319d6), &
                                        (3.055108535402016d-8,-2.4383079871034467d8), &
                                        (1.0120468663322205d-9,-8.35475287499606d9), &
                                        (2.961863655408945d-11,-3.193996910048163d11), &
                                        (7.762854623292005d-13,-1.3479619313593121d13), &
                                        (1.84207989989734d-14,-6.225136612864029d14), &
                                        (3.993136518072138d-16,-3.122820992919185d16), &
                                        (7.966833891596294d-18,-1.6910235523007962d18), &
                                        (1.4722303887260664d-19,-9.831370507829644d19), &
                                        (2.5336514351364773d-21,-6.108073162639263d21), &
                                        (4.0798958473925434d-23,-4.03863447368087d23), &
                                        (6.172739528132226d-25,-2.8315656516694476d25), &
                                        (8.806819572759488d-27,-2.0983379664327897d27), &
                                        (1.1887292025054836d-28,-1.6387861970351707d29), &
                                        (1.5224118385707597d-30,-1.3453438237909758d31), &
                                        (1.8548285821150714d-32,-1.1581948415962488d33), &
                                        (2.154899425974562d-34,-1.0433576439355702d35), &
                                        (2.3924222146081692d-36,-9.816001571599457d36), &
                                        (2.5432622334874715d-38,-9.627270774352158d38), &
                                        (2.593408024110376d-40,-9.826944471178756d40), &
                                        (2.5409629426208566d-42,-1.0423541817893215d43), &
                                        (2.395739738433546d-44,-1.1473011558077129d45), &
                                        (2.1767653402628094d-46,-1.308677103043361d47), &
                                        (1.908478015469465d-48,-1.545067672633029d49), &
                                        (1.6165895327065677d-50,-1.885926767218905d51), &
                                        (1.324486124924204d-52,-2.3773813393119423d53), &
                                        (1.0507418947683496d-54,-3.091953720279059d55), &
                                        (8.079463832699501d-57,-4.144928297774538d57), &
                                        (6.027224044196768d-59,-5.722223553255649d59), &
                                        (4.366030595832917d-61,-8.128534459820275d61), &
                                        (3.0736553148615166d-63,-1.1871767110561572d64), &
                                        (2.1045856173342726d-65,-1.781348033733495d66), &
                                        (1.4026412425449323d-67,-2.7441268303149825d68), &
                                        (9.105495828369005d-70,-4.336996290880954d70), &
                                        (5.761428895009189d-72,-7.027898308425972d72) /

            zeta = lib_math_bessel_riccati_zeta_derivative_real(x, fnu, n, zeta_n)
            ground_truth_zeta_n = lib_math_bessel_riccati_zeta_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_derivative_real:"
            do i=1, n
                buffer = 1 - abs(zeta(i) / ground_truth_zeta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = 1 - abs(zeta_n(i) / ground_truth_zeta_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_derivative_real

        function test_lib_math_bessel_riccati_zeta_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 40

            ! auxiliary
            complex(kind=8) :: x = cmplx(0.5,2)
            complex(kind=8), dimension(n) :: zeta
            complex(kind=8), dimension(n) :: zeta_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            complex(kind=8), dimension(n) :: ground_truth_zeta_n
            double precision :: ground_truth_e = 2* 10.0_8**(-12.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Zeta(n,x) = x*spherical_hankel2(n,x)
            !  >>> Zeta_(n,x) = derivative(Zeta(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Zeta_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Zeta_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            data ground_truth_zeta / (2.655785405411372_8,4.754722975801043_8), &
                                        (-1.5896272663501991_8,2.112619162507389_8), &
                                        (1.181113240759068_8,-1.6046575824487295_8), &
                                        (-2.497026285110109_8,-9.675554792839428_8), &
                                        (-54.1619275521725_8,1.7862004969948981_8), &
                                        (-71.58853334963285_8,344.297246746753_8), &
                                        (2361.9966083704394_8,1123.3935581884014_8), &
                                        (13870.621666663039_8,-17017.612720199275_8), &
                                        (-123629.31756853801_8,-163907.1289480185_8), &
                                        (-1.9456643803573349d6,833840.7429807571_8), &
                                        (3.872643248496936d6,2.3542825396964304d7), &
                                        (2.9095811075368017d8,2.3102332144738227d7), &
                                        (1.2183263184242601d9,-3.6538282600955157d9), &
                                        (-4.606990773521408d10,-2.918854937476698d10), &
                                        (-5.962203110486194d11,5.700089927375072d11), &
                                        (6.589731021680504d12,1.1519989554119607d13), &
                                        (2.1777705720932044d14,-6.170327391849623d13), &
                                        (-1.3234709326688866d14,-4.0769744542586245d15), &
                                        (-7.573075772240282d16,-1.624596261505917d16), &
                                        (-6.794935558335805d17,1.3881851191735555d18), &
                                        (2.4749170389503582d19,2.081096891176014d19), &
                                        (5.729427938046863d20,-4.156092931073134d20), &
                                        (-6.052442839943569d21,-1.500348657098113d22), &
                                        (-3.817110917006086d23,5.349588383161136d22), &
                                        (-1.0016678691781554d24,9.50346370463092d24), &
                                        (2.3130541170293197d26,8.42542227906415d25), &
                                        (3.680870190146898d27,-5.45409026478231d27), &
                                        (-1.219067340826215d29,-1.3547398178412652d29), &
                                        (-4.613763999274755d30,2.4508797649041417d30), &
                                        (3.737826150357287d31,1.5024128713037196d32), &
                                        (4.73804525992541d33,3.1205935954568174d30), &
                                        (3.6310729866190316d34,-1.4511657514664953d35), &
                                        (-4.2956550384848394d36,-2.2897861046401736d36), &
                                        (-1.093032825572374d38,1.2108379313295444d38), &
                                        (3.137949079756303d39,4.667502154198161d39), &
                                        (1.874679622730522d41,-6.785255227890062d40), &
                                        (-7.438978569589473d41,-7.22228496090334d42), &
                                        (-2.687116066847106d44,-3.8419504454953607d43), &
                                        (-3.926343064663271d45,9.642709773061577d45), &
                                        (3.305005841367725d47,2.4166514010847342d47) /

            zeta = lib_math_bessel_riccati_zeta_derivative_cmplx(x, fnu, n, zeta_n)
            ground_truth_zeta_n = lib_math_bessel_riccati_zeta_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_derivative_cmplx:"
            do i=1, n
                buffer = 1 - abs(real(zeta(i)) / real(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = 1 - abs(aimag(zeta(i)) / aimag(ground_truth_zeta(i)))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if

                buffer_y_n = 1 - abs(zeta_n(i) / ground_truth_zeta_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_derivative_cmplx

    end function lib_math_bessel_test_functions

end module lib_math_bessel
