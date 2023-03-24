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

module lib_math_legendre
    use lib_math_factorial
    use lib_math_type
    use lib_math_type_operator
    use libmath_generated_functions
    implicit none

    private

    ! --- public ---
    public :: lib_math_associated_legendre_polynomial
    public :: lib_math_associated_legendre_polynomial_with_negative_m
    public :: lib_math_associated_legendre_polynomial_range
    public :: lib_math_associated_legendre_polynomial_theta
    public :: lib_math_legendre_polynomial

    public :: lib_math_legendre_test_functions

    interface lib_math_associated_legendre_polynomial_theta
        module procedure lib_math_associated_legendre_polynomial_theta_wa
    end interface

    contains

        ! calculates the associated Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   m: integer
        !       order of the polynomial, >= 0
        !   fnu: integer
        !       degree of initial function, fnu .GE. 0
        !   n: integer
        !       number of members of the sequence, n .GE. 1
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pm: double precision, dimension(n)
        !       result of the associated Legendre polynomial
        !   pd: double precision, dimension(n)
        !       result of the deriviative of the associated Legendre polynomial
        !
        subroutine lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: fnu
            integer(kind=4), intent(in) :: n

            double precision, intent(inout) :: pm(n)
            double precision, intent(inout) :: pd(n)
            logical, optional :: condon_shortley_phase

            ! auxiliary
            double precision, dimension(0:n+1) :: buffer_pm
            double precision, dimension(0:n+1) :: buffer_pd

            call LPMNS(m, n+1, x, buffer_pm, buffer_pd)

            pm = buffer_pm(fnu:fnu+n-1)
            pd = buffer_pd(fnu:fnu+n-1)

            if (present(condon_shortley_phase)) then
                if (condon_shortley_phase) then
                    if (IAND(m, 1) .eq. 1) then
                        pm = -pm
                        pd = -pd
                    end if
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            end if

        end subroutine lib_math_associated_legendre_polynomial

        ! calculates the associated Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   m: integer
        !       order of the polynomial, >= 0
        !   fnu: integer
        !       degree of initial function, fnu .GE. 0
        !   n: integer
        !       number of members of the sequence, n .GE. 1
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pm: double precision, dimension(2, n)
        !       pm(1,fnu:fnu+n-1): result of the associated Legendre polynomial with a negativ m
        !       pm(2,fnu:fnu+n-1): result of the associated Legendre polynomial with a positiv m
        !   pd: double precision, dimension(2, n)
        !       pd(1,fnu:fnu+n-1): result of the deriviative of the associated Legendre polynomial with a negativ m
        !       pd(2,fnu:fnu+n-1): result of the deriviative of the associated Legendre polynomial with a positiv m
        !
        ! Refrence: http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
        subroutine lib_math_associated_legendre_polynomial_with_negative_m(x, m, fnu, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: fnu
            integer(kind=4), intent(in) :: n

            double precision, allocatable, dimension(:,:), intent(out) :: pm
            double precision, allocatable, dimension(:,:), intent(out) :: pd
            logical, optional :: condon_shortley_phase

            ! auxiliary
            integer :: i
            double precision, dimension(0:fnu+n) :: buffer_pm
            double precision, dimension(0:fnu+n) :: buffer_pd
            real(kind=8) :: buffer

            call LPMNS(m, fnu+n, x, buffer_pm, buffer_pd)

            if (allocated(pm)) deallocate(pm)
            allocate(pm(2, fnu:fnu+n-1))

            if (allocated(pd)) deallocate(pd)
            allocate(pd(2, fnu:fnu+n-1))

            pm(2, fnu:fnu+n-1) = buffer_pm(fnu:fnu+n-1)
            pd(2, fnu:fnu+n-1) = buffer_pd(fnu:fnu+n-1)

            ! calculation with negative m

            do i=fnu, fnu+n-1
                if (i .eq. 0) then
                    pm(1, 0) = pm(2, 0)
                    pd(1, 0) = pd(2, 0)
                else
                    buffer = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(i, m)

                    if ( IAND(m, 1) .eq. 1) then
                        ! m is odd
                        buffer = -buffer
                    end if

                    pm(1, i) = buffer * pm(2, i)
                    pd(1, i) = buffer * pd(2, i)
                end if
            end do

            if (present(condon_shortley_phase)) then
                if (condon_shortley_phase) then
                    if (IAND(m, 1) .eq. 1) then
                        ! m is odd
                        pm(2,:) = -pm(2,:)
                        pd(2,:) = -pd(2,:)
                    end if
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            end if

        end subroutine lib_math_associated_legendre_polynomial_with_negative_m

        ! calculates the associated Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   fnu: integer
        !       degree of initial function, fnu .GE. 0
        !   n: integer
        !       number of members of the sequence, n .GE. 1
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pm: type(list_list_real)
        !       result of the associated Legendre polynomial
        !   pd: double precision, dimension(2, n)
        !       result of the deriviative of the associated Legendre polynomial
        subroutine lib_math_associated_legendre_polynomial_range(x, fnu, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: fnu
            integer(kind=4), intent(in) :: n

            type(list_list_real) :: pm
            type(list_list_real) :: pd

            logical, optional :: condon_shortley_phase


            ! auxiliaray
            integer :: m_n
            integer :: m_m

            double precision, dimension(:, :), allocatable :: buffer_pm
            double precision, dimension(:, :), allocatable :: buffer_pd

            logical :: m_condon_shortley_phase

            m_condon_shortley_phase = .false.
            if (present(condon_shortley_phase)) m_condon_shortley_phase = condon_shortley_phase


            call init_list(pm, fnu, n)
            call init_list(pd, fnu, n)

!            allocate (buffer_pm(2, 0:fnu+n-1))
!            allocate (buffer_pd(2, 0:fnu+n-1))
            call lib_math_associated_legendre_polynomial_with_negative_m(x, 0, &
                                                                         0, fnu+n, &
                                                                         buffer_pm, buffer_pd, &
                                                                         m_condon_shortley_phase)
            do m_n=fnu, fnu+n-1
                pm%item(m_n)%item(0) = buffer_pm(1,m_n)

                pd%item(m_n)%item(0) = buffer_pd(1,m_n)
            end do

            do m_m=1, fnu+n-1
!                allocate (buffer_pm(2, m_m:fnu+n-1))
!                allocate (buffer_pd(2, m_m:fnu+n-1))
                call lib_math_associated_legendre_polynomial_with_negative_m(x, m_m, &
                                                                             m_m, fnu+n-m_m, &
                                                                             buffer_pm, buffer_pd, &
                                                                             m_condon_shortley_phase)
                do m_n=m_m, fnu+n-1
                    if (m_n .ge. fnu) then
                        pm%item(m_n)%item(-m_m) = buffer_pm(1,m_n)
                        pm%item(m_n)%item(m_m) = buffer_pm(2,m_n)

                        pd%item(m_n)%item(-m_m) = buffer_pd(1,m_n)
                        pd%item(m_n)%item(m_m) = buffer_pd(2,m_n)
                    end if
                end do
            end do

            deallocate (buffer_pm)
            deallocate (buffer_pd)
        end subroutine lib_math_associated_legendre_polynomial_range

        ! calculates the associated Legendre polynomial
        !
        ! Formula
        ! ----
        !   pm = P_1l(cos(theta)) / sin(theta)
        !
        !   pd = derivation[P_1l(cos(theta)), theta ]
        !
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       input value
        !   n_max: integer
        !       maximum degree of the polynomial, > 0
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pi_mn: double precision, dimension(n, -n:n)
        !       result of the associated Legendre polynomial (cos(theta)) / sin(theta)
        !       1st-dimension: degree n of the associated Legendre function
        !       2nd-dimension: order m of the associated Legendre function
        !   tau_mn: double precision, dimension(n, -n:n)
        !       result of the deriviative of the associated Legendre polynomial
        !       1st-dimension: degree n of the deriviative of the associated Legendre function
        !       2nd-dimension: order m of the deriviative of the associated Legendre function
        !
        !
        ! Formula
        ! ----
        !   pi_mn
        !       LaTeX: $$ \pi_{m n}(\cos \theta)=\frac{m}{\sin \theta} P_{n}^{m}(\cos \theta) $$
        !
        !   tau_mn
        !       LaTeX: $$ \tau_{m n}(\cos \theta)=\frac{d}{d \theta} P_{n}^{m}(\cos \theta) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, Appendix A
        !
        !
        subroutine lib_math_associated_legendre_polynomial_theta_Xu(theta, n_max, pi_nm, tau_nm, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            integer(kind=4), intent(in) :: n_max

            real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: pi_nm
            real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm
            logical, optional :: condon_shortley_phase

            ! auxiliary
            logical :: func_rv
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8), dimension(0:n_max+1, -n_max-1:n_max+1) :: m_pi_nm

            real(kind=8) :: x

            x = cos(theta)

            ! Recurrence Formulas for pi_mn, tau_mn
            ! where -1 .le. x .le. 1
            ! Special Values, eq. A6
            m_pi_nm(0,0) = 0.0_8
            m_pi_nm(1,0) = 0
            m_pi_nm(0,1) = 0
            m_pi_nm(1,1) = 1

            tau_nm(0,0) = 0
            tau_nm(1,0) = -sqrt(1.0_8 - x*x)
            tau_nm(0,1) = 0
            tau_nm(1,1) = x

            if (x .eq. 1.0_8) then
                ! Special Values, eq. 7
                m_pi_nm = 0
                tau_nm = 0
                do n = 1, n_max
                    m_pi_nm(n, -1) = 0.5_8
                    m_pi_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8

                    tau_nm(n, -1) = -0.5_8
                    tau_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                end do
            else if (x .eq. -1.0_8) then
                ! Special Values, eq. 7
                m_pi_nm = 0
                tau_nm = 0
                do n = 1, n_max
                    m_pi_nm(n, -1) = 0.5_8
                    tau_nm(n, -1) = -0.5_8
                    ! check if n+1 is even or odd
                    if (IAND(n+1_4, 1_4) .eq. 1) then
                        ! n+1 is odd
                        m_pi_nm(n, -1) = -m_pi_nm(n, -1)
                    else
                        ! n+1 is even
                        ! --> n is odd
                        tau_nm(n, -1) = -tau_nm(n, -1)
                    end if

                    m_pi_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                    tau_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                    ! check if n+1 is even or odd
                    if (IAND(n+1_4, 1_4) .eq. 1) then
                        ! is odd
                        m_pi_nm(n, -1) = -m_pi_nm(n, -1)
                    else
                        ! n+1 is even
                        ! --> n is odd
                        tau_nm(n, 1) = -tau_nm(n, 1)
                    end if
                end do
            else
                ! Recurrence Relations eq. A2, A3, A4, A5
                ! eq. A3: Electromagnetic scattering by an aggregate of spheres: errata

                ! Calculation of pi_mn
                ! ----
                !

                ! Step 1
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |  ^  |     |
                !  ^       |-----------------/-----------|
                !  |       |     |     |  /  |     |     |
                !          |-----------/-----------------|
                !          |  0  |  1  |     |     |     |
                !          |-----------------------------|
                !      0 --|  0  |  0  |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |   pi_m,n+1
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=2, n_max+1
                    func_rv = set_pi_n_n(x, n, m_pi_nm, n_max+1)
                end do

                ! Step 2
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |  x  |-->  |
                !  ^       |-----------------------------|
                !  |       |     |     |  x  |-------->  |
                !          |-----------------------------|
                !          |  0  |  1  |-------------->  |
                !          |-----------------------------|
                !      0 --|  0  |  0  |-------------->  |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |   pi_m,n+1
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=1, n_max
                    do m=0, n
                        func_rv = set_pi_m_n_plus_1(x, m, n, m_pi_nm, n_max+1)
                    end do
                end do

                ! Step 3
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |  x  |  x  |
                !  ^       |-----------------------------|
                !  |       |     |     |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |  0  |  1  |  x  |  x  |  x  |
                !          |-----------------------------|
                !      0 --|  0  |  0  |  x  |  x  |  x  |
                !          |--------|-----|-----|-----|--|
                !          |     |  V  |  |  |  |  |  |  |
                !          |--------------|-----|-----|--|
                !          |     |     |  V  |  |  |  |  |   pi_-m,n
                !          |--------------------|-----|--|
                !          |     |     |     |  V  |  |  |
                !  -n_max__---------------------------V---             x: calculated element
                !
                do n=1, n_max+1
                    do m=1, n
                        func_rv = set_pi_minus_m_n(m, n, m_pi_nm, n_max+1)
                    end do
                end do



                ! Calculation of tau_mn
                ! ----
                !
                ! Step 1
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |     |     |
                !  ^       |-----------------------------|
                !  |       |     |     |     |     |     |
                !          |-----------------------------|
                !          |  0  |  x  |     |     |     |
                !          |-----------------------------|
                !      0 --|  0  |  x  -------------->   |   tau_0,n+1
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=2, n_max-1
                    func_rv = set_tau_0_n_plus_1(x, n, tau_nm, n_max)
                end do

                ! Step 2
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |  ^  |  ^  |
                !  ^       |--------------------|-----|--|
                !  |       |     |     |  ^  |  |  |  |  |   tau_mn
                !          |--------------|-----|-----|--|
                !          |  0  |  x  |  |  |  |  |  |  |
                !          |-----------------------------|
                !      0 --|  0  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=2, n_max-1
                    do m=1, n
                        func_rv = set_tau_m_n(x, m, n, m_pi_nm, tau_nm, n_max+1, n_max)
                    end do
                end do

                ! Step 3
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |  x  |  x  |
                !  ^       |-----------------------------|
                !  |       |     |     |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |  0  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !      0 --|  0  |  x  |  x  |  x  |  x  |
                !          |--------|--------------------|
                !          |     |  V  |  |  |  |  |  |  |
                !          |--------------|-----|-----|--|
                !          |     |     |  V  |  |  |  |  |
                !          |--------------------|-----|--|
                !          |     |     |     |  V  |  V  |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=1, n_max
                    do m=1, n
                        func_rv = set_tau_minus_m_n(m, n, tau_nm, n_max)
                    end do
                end do
            end if

            pi_nm = m_pi_nm(0:n_max, -n_max:n_max)

            if (present(condon_shortley_phase)) then
                if (condon_shortley_phase) then
                    do n=0, n_max
                        do m=-n, n
                            if (IAND(abs(m), 1) .eq. 1) then
                                pi_nm(n, m) = -pi_nm(n, m)
                                tau_nm(n, m) = -tau_nm(n, m)
                            end if
                        end do
                    end do
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            end if

            ! reference: Absorption and Scattering of Light by smal particles eq. 4.47
            !  --> m=1
!            pm(0) = 0
!            pm(1) = 1
!            do i=2, n
!                pm(i) = (2*n - 1)/(n - 1) * cos_theta * pm(i-1) - n/(n - 1) * pm(i-2)
!                pd(i) = i * cos_theta * pm(i) - (n + 1) * pm(i-1)
!            end do

            contains

                ! first line of eq. A1
                ! Restriction
                ! ----
                !   - pi_nm(n, m) has to be known
                !   - pi_nm(n-1, m) has to be known
                function set_pi_m_n_plus_1(x, m, n, pi_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: pi_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    ! auxiliaray
                    real(kind=8) :: denominator

                    denominator = n - m + 1

                    pi_nm(n+1, m) = 0.0_8

                    if (pi_nm(n,m) .ne. 0.0_8) then
                        pi_nm(n+1, m) = (2*n + 1) / denominator * x * pi_nm(n,m)
                    end if
                    if (pi_nm(n-1,m) .ne. 0.0_8) then
                        pi_nm(n+1, m) = pi_nm(n+1, m) - (n + m) / denominator * pi_nm(n-1, m)
                    end if

!                    pi_nm(n+1, m) = (2*n + 1) / denominator * x * pi_nm(n,m) &
!                                    -(n + m) / denominator * pi_nm(n-1, m)

                    rv = .true.
                end function

                ! second line of eq. A1
                ! Restriction
                ! ----
                !   - m .ne. 1
                !   - pi_nm(n, m) has to be known
                !   - pi_nm(n, m-1) has to be known
                function set_pi_m_plus_1_n(x, m, n, pi_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: pi_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    pi_nm(n, m+1) = 0.0_8

                    if (pi_nm(n,m) .ne. 0.0_8) then
                        pi_nm(n, m+1) = 2*(m + 1)*x / sqrt(1 - x*x) * pi_nm(n,m)
                    end if
                    if (pi_nm(n,m-1) .ne. 0.0_8) then
                        pi_nm(n, m+1) = pi_nm(n, m+1) - (m + 1)*(n + m)*(n - m + 1)/(m - 1) * pi_nm(n,m-1)
                    end if

!                    pi_nm(n, m+1) = 2*(m + 1)*x / sqrt(1 - x*x) * pi_nm(n,m) &
!                                    -(m - 1)*(n + m)*(n - m + 1)/(m - 1) * pi_nm(n,m-1)

                    rv = .true.
                end function

                ! eq. A2
                ! Restriction
                ! ----
                !   - n .nq. 1
                !   - pi_nm(n-1, n-1) has to be known
                function set_pi_n_n(x, n, pi_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: pi_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    pi_nm(n,n) = sqrt(1 - x*x) * n*(2*n - 1)/(n - 1) * pi_nm(n-1,n-1)

                    rv = .true.
                end function

                ! first line of eq. A4
                !
                ! Restriction
                ! ----
                !   - pi_nm(n, m) has to be known
                function set_pi_minus_m_n(m, n, pi_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    if (pi_nm(n, m) .eq. 0.0_8) then
                        pi_nm(n, -m) = 0.0_8
                    else
                        pi_nm(n, -m) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m) * pi_nm(n,m)
                        ! check if m+1 is even or odd
                        if (IAND(m+1_4, 1_4) .eq. 1) then
                            ! m+1 is odd
                            pi_nm(n, -m) = -pi_nm(n, -m)
                        end if
                    end if

                    rv = .true.
                end function

                ! first line eq. A3
                !
                ! Restriction
                ! ----
                !   - m .ne. 0
                !   - pi_nm(n+1, m) has to be known
                !   - pi_nm(n, m) has to be known
                function set_tau_m_n(x, m, n, pi_nm, tau_nm, n_max_pi, n_max_tau) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max_pi, -n_max_pi:n_max_pi), intent(inout) :: pi_nm
                    real(kind=8), dimension(0:n_max_tau, -n_max_tau:n_max_tau), intent(inout) :: tau_nm
                    integer(kind=4) :: n_max_pi
                    integer(kind=4) :: n_max_tau

                    logical :: rv

                    tau_nm(n,m) = 0.0_8

                    if (pi_nm(n+1, m) .ne. 0.0_8) then
                        tau_nm(n,m) = (n - m + 1)/m * pi_nm(n+1, m)
                    end if

                    if (pi_nm(n, m) .ne. 0.0_8) then
                         tau_nm(n,m) = tau_nm(n,m) - (n + 1)/m * x * pi_nm(n, m)
                    end if

                    rv = .true.
                end function

                ! second line eq. A3
                !
                ! eq. A3: Electromagnetic scattering by an aggregate of spheres: errata 2001
                !
                ! Restriction
                ! ----
                !   - n .ne. 0
                !   - tau_nm(n, m) has to be known
                !   - tau_nm(n-1, m) has to be known
                function set_tau_0_n_plus_1(x, n, tau_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    tau_nm(n+1,0) = 0.0_8

                    if (tau_nm(n, 0) .ne. 0.0_8) then
                        tau_nm(n+1,0) = (2*n + 1)/n * x * tau_nm(n, 0)
                    end if

                    if (tau_nm(n-1, 0) .ne. 0.0_8) then
                        tau_nm(n+1,0) = tau_nm(n+1,0) - (n + 1)/n * tau_nm(n-1, 0)
                    end if

                    rv = .true.
                end function

                ! first line eq. A5
                !
                ! Restriction
                ! ----
                !   - tau_nm(n, m) has to be known
                function set_tau_minus_m_n(m, n, tau_nm, n_max) result (rv)
                    implicit none
                    ! dummy
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm
                    integer(kind=4) :: n_max

                    logical :: rv

                    if (tau_nm(n, m) .eq. 0.0_8) then
                        tau_nm(n, -m) = 0.0_8
                    else
                        tau_nm(n, -m) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m) * tau_nm(n,m)
                        ! check if m is even or odd
                        if (IAND(m, 1_4) .eq. 1) then
                            ! m is odd
                            tau_nm(n, -m) = -tau_nm(n, -m)
                        end if
                    end if

                    rv = .true.
                end function

        end subroutine lib_math_associated_legendre_polynomial_theta_xu

        ! calculates the associated Legendre polynomial
        !
        ! Formula
        ! ----
        !   pm = P_1l(cos(theta)) / sin(theta)
        !
        !   pd = derivation[P_1l(cos(theta)), theta ]
        !
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       input value
        !   n_max: integer
        !       maximum degree of the polynomial, > 0
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pi_mn: list_list_real
        !       result of the associated Legendre polynomial (cos(theta)) / sin(theta)
        !       1st-level: order m of the associated Legendre function
        !       2en-level: degree n of the associated Legendre function
        !   tau_mn: double precision, dimension(n, -n:n)
        !       result of the deriviative of the associated Legendre polynomial
        !       1st-level: degree n of the deriviative of the associated Legendre function
        !       2nd-level: order m of the deriviative of the associated Legendre function
        !
        !
        ! Formula
        ! ----
        !   pi_mn
        !       LaTeX: $$ \pi_{m n}(\cos \theta)=\frac{m}{\sin \theta} P_{n}^{m}(\cos \theta) $$
        !
        !   tau_mn
        !       LaTeX: $$ \tau_{m n}(\cos \theta)=\frac{d}{d \theta} P_{n}^{m}(\cos \theta) $$
        !
        ! Note
        ! ----
        !   Functions are generated by a Mathematica script.
        !
        !   todo: add script
        subroutine lib_math_associated_legendre_polynomial_theta_wa(theta, n_max, pi_nm, tau_nm, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            integer(kind=4), intent(in) :: n_max

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            logical, optional :: condon_shortley_phase

            ! auxiliary
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8) :: x

            call init_list(pi_nm, 0, n_max+1)
            call init_list(tau_nm, 0, n_max+1)

            x = cos(theta)

            pi_nm%item(0)%item(0) = 0.0_8
            tau_nm%item(0)%item(0) = 0

            if (x .eq. 1.0_8 .or.&
                x .eq. -1.0_8) then
!                !$OMP PARALLEL DO PRIVATE(n, m)
!                do n=1, n_max
!                    do m=-n, n
!                        pi_nm%item(n)%item(m) = get_associated_legendre_polynomial_limit(n, m)
!                    end do
!                end do
!                !$OMP END PARALLEL DO

                !$OMP PARALLEL DO PRIVATE(n)
                do n = 1, n_max
                    pi_nm%item(n)%item(:) = 0
                    pi_nm%item(n)%item(-1) = -0.5_8 !0.5_8

                    pi_nm%item(n)%item(1) = -n * (n + 1) / 2 !n * (n + 1) / 2
                end do
                !$OMP END PARALLEL DO

                if (x .eq. -1.0_8) then
                    !$OMP PARALLEL DO PRIVATE(n)
                    do n = 2, n_max, 2
                        pi_nm%item(n)%item(-1) = - pi_nm%item(n)%item(-1)
                        pi_nm%item(n)%item(1) = - pi_nm%item(n)%item(1)
                    end do
                    !$OMP END PARALLEL DO
                end if
            else
                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=1, n_max
                    do m=-n, n
                        pi_nm%item(n)%item(m) = get_associated_legendre_polynomial(n, m, theta)
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

            if (x .eq. 1.0_8 .or.&
                x .eq. -1.0_8) then
!                !$OMP PARALLEL DO PRIVATE(n, m)
!                do n=1, n_max
!                    do m=-n, n
!                        tau_nm%item(n)%item(m) = get_associated_legendre_polynomial_derivative_limit(n, m)
!                    end do
!                end do
!                !$OMP END PARALLEL DO

                !$OMP PARALLEL DO PRIVATE(n)
                do n = 1, n_max
                    tau_nm%item(n)%item(:) = 0
                    tau_nm%item(n)%item(-1) = 0.5_8 !-0.5_8

                    tau_nm%item(n)%item(1) = -n * (n + 1) / 2 !n * (n + 1) / 2
                end do
                !$OMP END PARALLEL DO

                if (x .eq. -1.0_8) then
                    !$OMP PARALLEL DO PRIVATE(n)
                    do n = 1, n_max, 2
                        tau_nm%item(n)%item(-1) = - pi_nm%item(n)%item(-1)
                        tau_nm%item(n)%item(1) = - pi_nm%item(n)%item(1)
                    end do
                    !$OMP END PARALLEL DO
                end if
            else
                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=1, n_max
                    do m=-n, n
                        tau_nm%item(n)%item(m) = get_associated_legendre_polynomial_derivative(n, m, theta)
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

            ! HINT: the calculation of the polynomial includes the condon shortley phase
            if (present(condon_shortley_phase)) then
                if (.not. condon_shortley_phase) then
                    ! remove phase
                    !$OMP PARALLEL DO PRIVATE(n, m)
                    do n=0, n_max
                        do m=-n, n
                            if (IAND(abs(m), 1) .eq. 1) then
                                pi_nm%item(n)%item(m) = -pi_nm%item(n)%item(m)
                                tau_nm%item(n)%item(m) = -tau_nm%item(n)%item(m)
                            end if
                        end do
                    end do
                    !$OMP END PARALLEL DO
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            else
                ! std: remove phase
                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=0, n_max
                    do m=-n, n
                        if (IAND(abs(m), 1) .eq. 1) then
                            pi_nm%item(n)%item(m) = -pi_nm%item(n)%item(m)
                            tau_nm%item(n)%item(m) = -tau_nm%item(n)%item(m)
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

        end subroutine

        ! calculates the Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   fnu: integer
        !       degree of initial function, fnu .GE. 0
        !   n: integer
        !       number of members of the sequence, n .GE. 1
        !
        ! Results
        ! ----
        !   pm: double precision, dimension(0:n)
        !       result of the Legendre polynomial
        !   pd: double precision, dimension(0:n)
        !       result of the deriviative of the Legendre polynomial
        !
        subroutine lib_math_legendre_polynomial(x, fnu, n, pn, pd)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: fnu
            integer(kind=4), intent(in) :: n

            double precision, dimension(n), intent(inout) :: pn
            double precision, dimension(n), intent(inout) :: pd

!#ifndef __GFORTRAN__
!            interface
!                subroutine lpn ( n, x, pn, pd )
!                    ! dummy
!                    integer ( kind = 4 ), intent(in) :: n
!                    real ( kind = 8 ), intent(in) :: x
!                    real ( kind = 8 ), dimension(n), intent(inout) :: pd
!                    real ( kind = 8 ), dimension(n), intent(inout) :: pn
!                end subroutine
!            end interface
!#endif            
            ! auxiliary
            double precision, dimension(0:fnu+n-1) :: buffer_pn
            double precision, dimension(0:fnu+n-1) :: buffer_pd
            
            
            call lpn ( fnu+n-1, x, buffer_pn, buffer_pd )

            pn = buffer_pn(fnu:fnu+n-1)
            pd = buffer_pd(fnu:fnu+n-1)

        end subroutine lib_math_legendre_polynomial

        function lib_math_legendre_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 5.0_8 * 10.0_8**(-5.0_8)

            rv = 0

            ! associated Legendre Polynomial
            if (.not. test_lib_math_associated_legendre_polynomial_m1()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m1_without_phase()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m2_without_phase()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m0()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m0_2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_range()) rv = rv + 1
            if (.not. test_lib_math_associated_legendre_polynomial_range_2()) rv = rv + 1
            if (.not. test_lib_math_associated_legendre_polynomial_with_negative_m1()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_with_negative_m2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_with_negative_m3()) then
                rv = rv + 1
            end if
!            if (.not. test_lib_math_associated_legendre_polynomial_theta_m1()) then
!                rv = rv + 1
!            end if
!            if (.not. test_lib_math_associated_legendre_polynomial_theta_xu_n0_3()) then
!                rv = rv + 1
!            end if
            if (.not. test_lib_math_associated_legendre_polynomial_theta_wa_n0_3()) then
                rv = rv + 1
            end if

            ! Legendre Polynomial
            if (.not. test_lib_math_legendre_polynomial()) then
                rv = rv + 1
            end if

            print *, "-------------lib_math_legendre_test_functions-------------"
            if (rv == 0) then
                print *, "lib_math_legendre_test_functions tests: OK"
            else
                print *, rv,"lib_math_legendre_test_functions test(s) FAILED"
            end if
            print *, "----------------------------------------------------------"


            contains

            function test_lib_math_associated_legendre_polynomial_m1() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 1

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 0.000000000000000_8
                ground_truth_pm(2) = -0.979795897113271_8
                ground_truth_pm(3) = -0.587877538267963_8
                ground_truth_pm(4) = 1.17575507653593_8
                ground_truth_pm(5) = 1.33252242007405_8
                ground_truth_pm(6) = -0.870058756636585_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 0.204124145231932_8
                ground_truth_pd(3) = -2.81691320420066_8
                ground_truth_pd(4) = -3.18433666561813_8
                ground_truth_pd(5) = 5.01328900689624_8
                ground_truth_pd(6) = 9.23457633029258_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m1:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m1

            function test_lib_math_associated_legendre_polynomial_m1_without_phase() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 1

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = -0.000000000000000_8
                ground_truth_pm(2) = 0.979795897113271_8
                ground_truth_pm(3) = 0.587877538267963_8
                ground_truth_pm(4) = -1.17575507653593_8
                ground_truth_pm(5) = -1.33252242007405_8
                ground_truth_pm(6) = 0.870058756636585_8

                ground_truth_pd(1) = 0.000000000000000_8
                ground_truth_pd(2) = -0.204124145231932_8
                ground_truth_pd(3) = 2.81691320420066_8
                ground_truth_pd(4) = 3.18433666561813_8
                ground_truth_pd(5) = -5.01328900689624_8
                ground_truth_pd(6) = -9.23457633029258_8

                x = 0.2

!                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .false.)
                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m1_without_phase:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m1_without_phase

            function test_lib_math_associated_legendre_polynomial_m2() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 2

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=2
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 0.000000000000000_8
                ground_truth_pm(2) = -1.06581410364015D-16
                ground_truth_pm(3) = 2.88000000000000_8
                ground_truth_pm(4) = 2.88000000000000_8
                ground_truth_pm(5) = -5.18400000000000_8
                ground_truth_pm(6) = -8.87040000000000_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 4.81867632215780D-16
                ground_truth_pd(3) = -1.20000000000000_8
                ground_truth_pd(4) = 13.2000000000000_8
                ground_truth_pd(5) = 22.3200000000000_8
                ground_truth_pd(6) = -28.5600000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m2:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m2

            function test_lib_math_associated_legendre_polynomial_m2_without_phase() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 2

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=2
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 0.000000000000000_8
                ground_truth_pm(2) = -1.06581410364015D-16
                ground_truth_pm(3) = 2.88000000000000_8
                ground_truth_pm(4) = 2.88000000000000_8
                ground_truth_pm(5) = -5.18400000000000_8
                ground_truth_pm(6) = -8.87040000000000_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 4.81867632215780D-16
                ground_truth_pd(3) = -1.20000000000000_8
                ground_truth_pd(4) = 13.2000000000000_8
                ground_truth_pd(5) = 22.3200000000000_8
                ground_truth_pd(6) = -28.5600000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, .false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m2_without_phase:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m2_without_phase

            function test_lib_math_associated_legendre_polynomial_m0() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 0

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=0
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 1.00000000000000_8
                ground_truth_pm(2) = 0.200000000000000_8
                ground_truth_pm(3) = -0.440000000000000_8
                ground_truth_pm(4) = -0.280000000000000_8
                ground_truth_pm(5) = 0.232000000000000_8
                ground_truth_pm(6) = 0.307520000000000_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 1.00000000000000_8
                ground_truth_pd(3) = 0.600000000000000_8
                ground_truth_pd(4) = -1.20000000000000_8
                ground_truth_pd(5) = -1.36000000000000_8
                ground_truth_pd(6) = 0.888000000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m0:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m0

            function test_lib_math_associated_legendre_polynomial_m0_2() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                integer(kind=4), parameter :: m = 0

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=0
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=math.cos(0.5)
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 1.00000000000000_8
                ground_truth_pm(2) = 0.877582561890373_8
                ground_truth_pm(3) = 0.655226729401105_8
                ground_truth_pm(4) = 0.373304211751205_8
                ground_truth_pm(5) = 0.0818891693470757_8
                ground_truth_pm(6) = -0.169287256752937_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 1.00000000000000_8
                ground_truth_pd(3) = 2.63274768567112_8
                ground_truth_pd(4) = 4.27613364700553_8
                ground_truth_pd(5) = 5.24587716792955_8
                ground_truth_pd(6) = 5.01313617112921_8

                x = cos(0.5_8)

                call lib_math_associated_legendre_polynomial(x, m, fnu, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m0_2:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m0_2

            function test_lib_math_associated_legendre_polynomial_with_negative_m1() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 1
                integer(kind=4), parameter :: n = 4
                integer(kind=4), parameter :: m = 1

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(:,:), allocatable :: pm
                double precision, dimension(:,:), allocatable :: pd

                double precision, dimension(2, n) :: ground_truth_pm
                double precision, dimension(2, n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with mathematica 8
                !
                ! source code:
                !  >>>  LPwa[l_, m_, x_] := (-1)^m LegendreP[l, m, x];
                !  >>>  makeMneg[l_, m_] := (-1)^m (l - m)!/(l + m)!;
                !  >>>
                !  >>>  lmax = 4 ;
                !  >>>  theta = 0.2 ;
                !  >>>
                !  >>>  For[m1 = 1, m1 <= 3, m1++ ,
                !  >>>   For[l1 = m1, l1 <= lmax, l1++ ,
                !  >>>     buffer = FullSimplify[LPwa[l1, m1, x]];
                !  >>>     Print["l = ", l1, ", m = ", -m1, ": ",
                !  >>>       N[makeMneg[l1, m1] buffer /. x -> theta, 16]]
                !  >>>      Print["l = ", l1, ", m = ", m1, ": ",
                !  >>>       N[buffer /. x -> theta, 16]]
                !  >>>
                !  >>>     ]
                !  >>>    Print["------------------"]
                !  >>>   ]

                ground_truth_pm(1, 1) = -0.489898_8
                ground_truth_pm(2, 1) =  0.979796_8
                ground_truth_pm(1, 2) = -0.0979796_8
                ground_truth_pm(2, 2) =  0.587878_8
                ground_truth_pm(1, 3) =  0.0979796_8
                ground_truth_pm(2, 3) = -1.17576_8
                ground_truth_pm(1, 4) =  0.0666261_8
                ground_truth_pm(2, 4) = -1.33252_8

                ground_truth_pd(1, 1) = 0.102062_8
                ground_truth_pd(2, 1) = -0.204124_8
                ground_truth_pd(1, 2) = -0.469486_8
                ground_truth_pd(2, 2) = 2.81691_8
                ground_truth_pd(1, 3) = -0.265361_8
                ground_truth_pd(2, 3) = 3.18434_8
                ground_truth_pd(1, 4) = 0.250664_8
                ground_truth_pd(2, 4) = -5.01329_8

                x = 0.2

                call lib_math_associated_legendre_polynomial_with_negative_m(x, m, fnu, n, pm, pd, &
                                                                             condon_shortley_phase=.false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_with_negative_m1:"
                do i=1, n
                    buffer = abs(pm(1,i) - ground_truth_pm(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(1,i) - ground_truth_pd(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_with_negative_m1

            function test_lib_math_associated_legendre_polynomial_with_negative_m2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 2
                integer(kind=4), parameter :: n = 3
                integer(kind=4), parameter :: m = 2

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(:,:), allocatable :: pm
                double precision, dimension(:,:), allocatable :: pd

                double precision, dimension(2, n) :: ground_truth_pm
                double precision, dimension(2, n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with mathematica 8
                !
                ! source code:
                !  >>>  LPwa[l_, m_, x_] := (-1)^m LegendreP[l, m, x];
                !  >>>  makeMneg[l_, m_] := (-1)^m (l - m)!/(l + m)!;
                !  >>>
                !  >>>  lmax = 4 ;
                !  >>>  theta = 0.2 ;
                !  >>>
                !  >>>  For[m1 = 1, m1 <= 3, m1++ ,
                !  >>>   For[l1 = m1, l1 <= lmax, l1++ ,
                !  >>>     buffer = FullSimplify[LPwa[l1, m1, x]];
                !  >>>     Print["l = ", l1, ", m = ", -m1, ": ",
                !  >>>       N[makeMneg[l1, m1] buffer /. x -> theta, 16]]
                !  >>>      Print["l = ", l1, ", m = ", m1, ": ",
                !  >>>       N[buffer /. x -> theta, 16]]
                !  >>>
                !  >>>     ]
                !  >>>    Print["------------------"]
                !  >>>   ]

                ground_truth_pm(1, 1) = 0.12_8
                ground_truth_pm(2, 1) = 2.88_8
                ground_truth_pm(1, 2) = 0.024_8
                ground_truth_pm(2, 2) = 2.88_8
                ground_truth_pm(1, 3) = -0.0144_8
                ground_truth_pm(2, 3) = -5.184_8

                ground_truth_pd(1, 1) = -0.05_8
                ground_truth_pd(2, 1) = -1.2_8
                ground_truth_pd(1, 2) = 0.11_8
                ground_truth_pd(2, 2) = 13.2_8
                ground_truth_pd(1, 3) = 0.062_8
                ground_truth_pd(2, 3) = 22.32_8

                x = 0.2

                call lib_math_associated_legendre_polynomial_with_negative_m(x, m, fnu, n, pm, pd, &
                                                                             condon_shortley_phase=.false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_with_negative_m2:"
                do i=1, n
                    buffer = abs(pm(1,fnu+i-1) - ground_truth_pm(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(1,fnu+i-1) - ground_truth_pd(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_with_negative_m2

            function test_lib_math_associated_legendre_polynomial_with_negative_m3() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 3
                integer(kind=4), parameter :: n = 2
                integer(kind=4), parameter :: m = 3

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                double precision, dimension(:,:), allocatable :: pm
                double precision, dimension(:,:), allocatable :: pd

                double precision, dimension(2, n) :: ground_truth_pm
                double precision, dimension(2, n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with mathematica 8
                !
                ! source code:
                !  >>>  LPwa[l_, m_, x_] := (-1)^m LegendreP[l, m, x];
                !  >>>  makeMneg[l_, m_] := (-1)^m (l - m)!/(l + m)!;
                !  >>>
                !  >>>  lmax = 4 ;
                !  >>>  theta = 0.2 ;
                !  >>>
                !  >>>  For[m1 = 1, m1 <= 3, m1++ ,
                !  >>>   For[l1 = m1, l1 <= lmax, l1++ ,
                !  >>>     buffer = FullSimplify[LPwa[l1, m1, x]];
                !  >>>     Print["l = ", l1, ", m = ", -m1, ": ",
                !  >>>       N[makeMneg[l1, m1] buffer /. x -> theta, 16]]
                !  >>>      Print["l = ", l1, ", m = ", m1, ": ",
                !  >>>       N[buffer /. x -> theta, 16]]
                !  >>>
                !  >>>     ]
                !  >>>    Print["------------------"]
                !  >>>   ]

                ground_truth_pm(1, 1) = -0.0195959_8
                ground_truth_pm(2, 1) = 14.1091_8
                ground_truth_pm(1, 2) = -0.00391918_8
                ground_truth_pm(2, 2) = 19.7527_8

                ground_truth_pd(1, 1) = 0.0122474_8
                ground_truth_pd(2, 1) = -8.81816_8
                ground_truth_pd(1, 2) = -0.0171464_8
                ground_truth_pd(2, 2) = 86.418_8

                x = 0.2

                call lib_math_associated_legendre_polynomial_with_negative_m(x, m, fnu, n, pm, pd, &
                                                                             condon_shortley_phase=.false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_with_negative_m3:"
                do i=1, n
                    buffer = abs(pm(1,fnu+i-1) - ground_truth_pm(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(1,fnu+i-1) - ground_truth_pd(1,i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_with_negative_m3

            function test_lib_math_associated_legendre_polynomial_theta_m1() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: n_max = 5

                ! auxiliary
                integer(kind=4) :: i
                double precision :: theta

                real(kind=8), dimension(0:n_max, -n_max:n_max) :: pi_nm
                real(kind=8), dimension(0:n_max, -n_max:n_max) :: tau_nm

                double precision, dimension(0:n_max) :: ground_truth_pi_n1
                double precision, dimension(0:n_max) :: ground_truth_tau_n1

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=cos(0.2)
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_tau_n1(0) = 0.000000000000000_8
                ground_truth_tau_n1(1) = -4.93315487558688_8
                ground_truth_tau_n1(2) = -13.9084526582466_8
                ground_truth_tau_n1(3) = -25.2179729025493_8
                ground_truth_tau_n1(4) = -36.4802655764097_8
                ground_truth_tau_n1(5) = -44.8436276630205_8

                theta = 0.2_8

                call genlgp(theta, ground_truth_pi_n1, n_max+1)

                call lib_math_associated_legendre_polynomial_theta_xu(theta, n_max, pi_nm, tau_nm)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_theta_m1:"
                do i=0, n_max
                    buffer = abs(pi_nm(i, 1) - ground_truth_pi_n1(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n_max
                    buffer = abs(tau_nm(i, 1) - ground_truth_tau_n1(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_theta_m1

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
                  real(kind=8) theta, pnmllg, costh, rn
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
          end subroutine genlgp

          function test_lib_math_associated_legendre_polynomial_theta_xu_n0_3() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: n_max = 3

                ! auxiliary
                integer(kind=4) :: i
                integer(kind=4) :: ii
                double precision :: theta

                real(kind=8), dimension(0:n_max, -n_max:n_max) :: pi_nm
                real(kind=8), dimension(0:n_max, -n_max:n_max) :: tau_nm

                double precision, dimension(0:n_max, -n_max:n_max) :: ground_truth_pi_nm
                double precision, dimension(0:n_max, -n_max:n_max) :: ground_truth_tau_nm

                double precision :: buffer

                ! Values were generated with WolframAplpha
                !
                ! source code:
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 0, m = 0, t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 1, m = [-1, 0, 1], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 2, m = [-2, -1, 0, 1, 2], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 3, m = [-3, -2, -1, 0, 1, 2, 3], t = 2

!                ground_truth_pi_nm(0, -3) = -2.077165092735346_8
!                ground_truth_pi_nm(0, -2) = 0_8
!                ground_truth_pi_nm(0, -1) = -1.712759410407380_8
                ground_truth_pi_nm(0,  0) = 0_8
!                ground_truth_pi_nm(0,  1) = 0_8
!                ground_truth_pi_nm(0,  2) = 0_8
!                ground_truth_pi_nm(0,  3) = 0_8

!                ground_truth_pi_nm(1, -3) = -1.341772398969518_8
!                ground_truth_pi_nm(1, -2) = 0.0_8
                ground_truth_pi_nm(1, -1) = -0.5000000000000000_8
                ground_truth_pi_nm(1,  0) = 0_8
                ground_truth_pi_nm(1,  1) = -1.000000000000000_8
!                ground_truth_pi_nm(1,  2) = 0_8
!                ground_truth_pi_nm(1,  3) = 0_8

!                ground_truth_pi_nm(2, -3) = -0.4958414335756772_8
                ground_truth_pi_nm(2, -2) = -0.2273243567064204_8
                ground_truth_pi_nm(2, -1) = 0.2080734182735712_8
                ground_truth_pi_nm(2,  0) = 0.0_8
                ground_truth_pi_nm(2,  1) = 1.248440509641427_8
                ground_truth_pi_nm(2,  2) = 5.455784560954090_8
!                ground_truth_pi_nm(2,  3) = 0_8

                ground_truth_pi_nm(3, -3) = -0.05167636315198787_8
                ground_truth_pi_nm(3, -2) = 0.09460031191349103_8
                ground_truth_pi_nm(3, -1) = 0.01676363151987872_8
                ground_truth_pi_nm(3,  0) = 0.0_8
                ground_truth_pi_nm(3,  1) = 0.2011635782385447_8
                ground_truth_pi_nm(3,  2) = -11.35203742961892_8
                ground_truth_pi_nm(3,  3) = -37.20698146943127_8

                ground_truth_tau_nm(0,  0) = 0_8

                ground_truth_tau_nm(1, -1) = 0.0_8
                ground_truth_tau_nm(1,  0) = 0.0_8
                ground_truth_tau_nm(1,  1) = 0.0_8

                ground_truth_tau_nm(2, -2) =  0.1040367091367856_8
                ground_truth_tau_nm(2, -1) =  0.4546487134128408_8
                ground_truth_tau_nm(2,  0) = 0_8
                ground_truth_tau_nm(2,  1) = 2.727892280477045_8
                ground_truth_tau_nm(2,  2) = -2.496881019282854_8

                ground_truth_tau_nm(3, -3) =  0.04730015595674552_8
                ground_truth_tau_nm(3, -2) =  0.1634109052159030_8
                ground_truth_tau_nm(3, -1) =  -0.4730015595674552_8
                ground_truth_tau_nm(3,  0) = 0_8
                ground_truth_tau_nm(3,  1) = -5.676018714809462_8
                ground_truth_tau_nm(3,  2) = -19.60930862590836_8
                ground_truth_tau_nm(3,  3) = 34.05611228885677_8

                theta = 2.0_8

                call lib_math_associated_legendre_polynomial_theta_xu(theta, n_max, pi_nm, tau_nm, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_theta_Xu_n0_3:"
                do i=0, n_max
                    !do ii=-n_max, n_max
                    do ii=-i, i
                        buffer = abs(pi_nm(i, ii) - ground_truth_pi_nm(i, ii))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i, " m: ", ii , "difference: ", buffer, " : FAILED"
                            print *, "      pi_nm: ", pi_nm(i, ii)
                            print *, "      truth: ", ground_truth_pi_nm(i, ii)
                            rv = .false.
                        else
                            print *, "  n: ", i, " m: ", ii, "OK"
                        end if
                    end do
                    print*, ""
                end do

                print*, "  deriviation:"
                do i=0, n_max
                    do ii=-i, i
                        buffer = abs(tau_nm(i, ii) - ground_truth_tau_nm(i, ii))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i, " m: ", ii , "difference: ", buffer, " : FAILED"
                            print *, "  ","                                    tau_nm: ", tau_nm(i, ii)
                            print *, "  ","                                    truth : ", ground_truth_tau_nm(i, ii)
                            rv = .false.
                        else
                            print *, "  n: ", i, " m: ", ii, "OK"
                        end if
                    end do
                    print*, ""
                end do

            end function test_lib_math_associated_legendre_polynomial_theta_xu_n0_3

            function test_lib_math_associated_legendre_polynomial_theta_wa_n0_3() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: n_max = 3

                ! auxiliary
                integer(kind=4) :: n
                integer(kind=4) :: m
                double precision :: theta

                type(list_list_real) :: pi_nm
                type(list_list_real) :: tau_nm

                type(list_list_real) :: ground_truth_pi_nm
                type(list_list_real) :: ground_truth_tau_nm

                double precision :: buffer

                call init_list(ground_truth_pi_nm, 0, n_max+1)
                call init_list(ground_truth_tau_nm, 0, n_max+1)

                ! Values were generated with WolframAplpha.
                !
                ! source code:
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 0, m = 0, t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 1, m = [-1, 0, 1], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 2, m = [-2, -1, 0, 1, 2], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 3, m = [-3, -2, -1, 0, 1, 2, 3], t = 2
                ground_truth_pi_nm%item(0)%item(0) = 0_8

                ground_truth_pi_nm%item(1)%item(-1) = -0.5000000000000000_8
                ground_truth_pi_nm%item(1)%item(0) = 0_8
                ground_truth_pi_nm%item(1)%item(1) = -1.000000000000000_8

                ground_truth_pi_nm%item(2)%item(-2) = -0.2273243567064204_8
                ground_truth_pi_nm%item(2)%item(-1) = 0.2080734182735712_8
                ground_truth_pi_nm%item(2)%item( 0) = 0.0_8
                ground_truth_pi_nm%item(2)%item( 1) = 1.248440509641427_8
                ground_truth_pi_nm%item(2)%item( 2) = 5.455784560954090_8

                ground_truth_pi_nm%item(3)%item(-3) = -0.05167636315198787_8
                ground_truth_pi_nm%item(3)%item(-2) = 0.09460031191349103_8
                ground_truth_pi_nm%item(3)%item(-1) = 0.01676363151987872_8
                ground_truth_pi_nm%item(3)%item( 0) = 0.0_8
                ground_truth_pi_nm%item(3)%item( 1) = 0.2011635782385447_8
                ground_truth_pi_nm%item(3)%item( 2) = -11.35203742961892_8
                ground_truth_pi_nm%item(3)%item( 3) = -37.20698146943127_8

                ! Values were generated with Mathematica
                ! -> libmath_generated_functions/src/lib/legendre_poly/generator/legendre2_onlyD.nb
                ground_truth_tau_nm%item(0)%item(0) = 0_8

                ground_truth_tau_nm%item(1)%item(-1) = -0.208073_8
                ground_truth_tau_nm%item(1)%item( 0) = -0.909297_8
                ground_truth_tau_nm%item(1)%item( 1) = 0.416147_8

                ground_truth_tau_nm%item(2)%item(-2) = -0.0946003_8
                ground_truth_tau_nm%item(2)%item(-1) = -0.326822_8
                ground_truth_tau_nm%item(2)%item( 0) = 1.1352_8
                ground_truth_tau_nm%item(2)%item( 1) = 1.96093_8
                ground_truth_tau_nm%item(2)%item( 2) = -2.27041_8

                ground_truth_tau_nm%item(3)%item(-3) = -0.021505_8
                ground_truth_tau_nm%item(3)%item(-2) = -0.0546107_8
                ground_truth_tau_nm%item(3)%item(-1) =  0.437075_8
                ground_truth_tau_nm%item(3)%item( 0) = 0.182918_8
                ground_truth_tau_nm%item(3)%item( 1) = -5.2449_8
                ground_truth_tau_nm%item(3)%item( 2) = -6.55329_8
                ground_truth_tau_nm%item(3)%item( 3) = 15.4836_8

                theta = 2.0_8

                call lib_math_associated_legendre_polynomial_theta_wa(theta, n_max, pi_nm, tau_nm, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_theta_wolfram_alpha_n0_3:"
                do n=0, n_max
                    !do ii=-n_max, n_max
                    do m=-n, n
                        buffer = abs(pi_nm%item(n)%item(m) - ground_truth_pi_nm%item(n)%item(m))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", n, " m: ", m , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", n, " m: ", m, "OK"
                        end if
                    end do
                    print*, ""
                end do

                print*, "  deriviation:"
                do n=0, n_max
                    !do ii=-n_max, n_max
                    do m=-n, n
                        buffer = abs(tau_nm%item(n)%item(m) - ground_truth_tau_nm%item(n)%item(m))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", n, " m: ", m , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", n, " m: ", m, "OK"
                        end if
                    end do
                    print*, ""
                end do

            end function test_lib_math_associated_legendre_polynomial_theta_wa_n0_3

            function test_lib_math_legendre_polynomial() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 6
                double precision, dimension(n) :: pm
                double precision, dimension(n) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>> var('l')
                !  >>> var('x')
                !  >>>
                !  >>> P(l,x) = legendre_P(l,x)
                !  >>> P_d(l,x) = derivative(legendre_P(l,x), x)
                !  >>>
                !  >>> x=0.4
                !  >>>
                !  >>> for i in range(0,6):
                !  >>>     value = numerical_approx(P(i,x))
                !  >>>     print("l = {}: {}".format(i, value))
                !  >>>
                !  >>> print("\nderivative:")
                !  >>> for i in range(0,6):
                !  >>>     value = numerical_approx(P_d(i,x))
                !  >>>     print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 1.00000000000000_8
                ground_truth_pm(2) = 0.400000000000000_8
                ground_truth_pm(3) = -0.260000000000000_8
                ground_truth_pm(4) = -0.440000000000000_8
                ground_truth_pm(5) = -0.113000000000000_8
                ground_truth_pm(6) = 0.270640000000000_8

                ground_truth_pd(1) = -0.000000000000000_8
                ground_truth_pd(2) = 1.00000000000000_8
                ground_truth_pd(3) = 1.20000000000000_8
                ground_truth_pd(4) = -0.300000000000000_8
                ground_truth_pd(5) = -1.88000000000000_8
                ground_truth_pd(6) = -1.31700000000000_8

                x = 0.4

                call lib_math_legendre_polynomial(x, fnu, n, pm, pd)


                rv = .true.
                print *, "test_lib_math_legendre_polynomial:"
                do i=1, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=1, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_legendre_polynomial

            function test_lib_math_associated_legendre_polynomial_range() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 1
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 1

                ! auxiliary
                integer :: i
                double precision :: x

                type(list_list_real) :: pm
                type(list_list_real) :: pd

                double precision, dimension(n) :: ground_truth_pm
                double precision, dimension(n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(1) = 0.979795897113271_8
                ground_truth_pm(2) = 0.587877538267963_8
                ground_truth_pm(3) = -1.17575507653593_8
                ground_truth_pm(4) = -1.33252242007405_8
                ground_truth_pm(5) = 0.870058756636585_8

                ground_truth_pd(1) = -0.204124145231932_8
                ground_truth_pd(2) = 2.81691320420066_8
                ground_truth_pd(3) = 3.18433666561813_8
                ground_truth_pd(4) = -5.01328900689624_8
                ground_truth_pd(5) = -9.23457633029258_8

                x = 0.2

                call lib_math_associated_legendre_polynomial_range(x, fnu, n, pm, pd)

                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_range:"
                do i=fnu, fnu+n-1
                    buffer = abs(pm%item(i)%item(m) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=fnu, fnu+n-1
                    buffer = abs(pd%item(i)%item(m) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do
            end function

            function test_lib_math_associated_legendre_polynomial_range_2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: fnu = 0
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 1

                ! auxiliary
                integer :: i
                integer :: ii
                double precision :: x

                type(list_list_real) :: pm
                type(list_list_real) :: pd

                type(list_list_real) :: ground_truth_pm
                type(list_list_real) :: ground_truth_pd

                double precision :: buffer

                call init_list(ground_truth_pm, fnu, n)
                call init_list(ground_truth_pd, fnu, n)

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>> LP[m_, n_, x_] := (-1)^m LegendreP[n, m, x];
                !  >>> LPDt = D[LP[md, nd, xd], xd];
                !  >>> LPD[m_, n_, x_] := LPDt /. {md -> m, nd -> n, xd -> x};
                !  >>>
                !  >>> x = 0.2
                !  >>>
                !  >>> nStart = 0;
                !  >>> nStop = 4;
                !  >>>
                !  >>> For[n = nStart, n <= nStop, n = n + 1,
                !  >>>  For[m = -n, m <= n, m = m + 1,
                !  >>>   val = LP[m, n, x];
                !  >>>   Print["%item(", n, ")%item(", m, ") = ", val];
                !  >>>   ]
                !  >>>  ]
                !  >>> Print[" --- derivative ---"]
                !  >>> For[n = nStart, n <= nStop, n = n + 1,
                !  >>>  For[m = -n, m <= n, m = m + 1,
                !  >>>   val = LPD[m, n, x];
                !  >>>   Print["%item(", n, ")%item(", m, ") = ", val];
                !  >>>   ]
                !  >>>  ]
                !  >>>
                ground_truth_pm%item(0)%item(0) = 1.
                ground_truth_pm%item(1)%item(-1) = -0.489898
                ground_truth_pm%item(1)%item(0) = 0.2
                ground_truth_pm%item(1)%item(1) = 0.979796
                ground_truth_pm%item(2)%item(-2) = 0.12
                ground_truth_pm%item(2)%item(-1) = -0.0979796
                ground_truth_pm%item(2)%item(0) = -0.44
                ground_truth_pm%item(2)%item(1) = 0.587878
                ground_truth_pm%item(2)%item(2) = 2.88
                ground_truth_pm%item(3)%item(-3) = -0.0195959
                ground_truth_pm%item(3)%item(-2) = 0.024
                ground_truth_pm%item(3)%item(-1) = 0.0979796
                ground_truth_pm%item(3)%item(0) = -0.28
                ground_truth_pm%item(3)%item(1) = -1.17576
                ground_truth_pm%item(3)%item(2) = 2.88
                ground_truth_pm%item(3)%item(3) = 14.1091
                ground_truth_pm%item(4)%item(-4) = 0.0024
                ground_truth_pm%item(4)%item(-3) = -0.00391918
                ground_truth_pm%item(4)%item(-2) = -0.0144
                ground_truth_pm%item(4)%item(-1) = 0.0666261
                ground_truth_pm%item(4)%item(0) = 0.232
                ground_truth_pm%item(4)%item(1) = -1.33252
                ground_truth_pm%item(4)%item(2) = -5.184
                ground_truth_pm%item(4)%item(3) = 19.7527
                ground_truth_pm%item(4)%item(4) = 96.768

                ground_truth_pd%item(0)%item(0) = 2.89121D-17
                ground_truth_pd%item(1)%item(-1) = 0.102062
                ground_truth_pd%item(1)%item(0) = 1.
                ground_truth_pd%item(1)%item(1) = -0.204124
                ground_truth_pd%item(2)%item(-2) = -0.05
                ground_truth_pd%item(2)%item(-1) = -0.469486
                ground_truth_pd%item(2)%item(0) = 0.6
                ground_truth_pd%item(2)%item(1) = 2.81691
                ground_truth_pd%item(2)%item(2) = -1.2
                ground_truth_pd%item(3)%item(-3) = 0.0122474
                ground_truth_pd%item(3)%item(-2) = 0.11
                ground_truth_pd%item(3)%item(-1) = -0.265361
                ground_truth_pd%item(3)%item(0) = -1.2
                ground_truth_pd%item(3)%item(1) = 3.18434
                ground_truth_pd%item(3)%item(2) = 13.2
                ground_truth_pd%item(3)%item(3) = -8.81816
                ground_truth_pd%item(4)%item(-4) = -0.002
                ground_truth_pd%item(4)%item(-3) = -0.0171464
                ground_truth_pd%item(4)%item(-2) = 0.062
                ground_truth_pd%item(4)%item(-1) = 0.250664
                ground_truth_pd%item(4)%item(0) = -1.36
                ground_truth_pd%item(4)%item(1) = -5.01329
                ground_truth_pd%item(4)%item(2) = 22.32
                ground_truth_pd%item(4)%item(3) = 86.418
                ground_truth_pd%item(4)%item(4) = -80.64

                x = 0.2

                call lib_math_associated_legendre_polynomial_range(x, fnu, n, pm, pd)

                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_range:"
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        buffer = abs(pm%item(i)%item(ii) - ground_truth_pm%item(i)%item(ii))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n:", i , " m:", ii, " difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n:", i , " m:", ii, ": OK"
                        end if
                    end do
                end do

                print*, "  deriviation:"
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        buffer = abs(pd%item(i)%item(ii) - ground_truth_pd%item(i)%item(ii))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n:", i , " m:", ii, " difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n:", i , " m:", ii, ": OK"
                        end if
                    end do
                end do
            end function test_lib_math_associated_legendre_polynomial_range_2

        end function lib_math_legendre_test_functions

end module lib_math_legendre
