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

!#define _PRINT_NOTE_

module lib_mie_ss_helper_functions
    use libmath
    use lib_mie_type
    use lib_mie_illumination
    implicit none

    private

    ! public parameter
    integer, parameter, public :: SINGLE_SPHERE_N_MAX = 45

    ! public functions
!    public :: lib_mie_ss_hf_contructor
    public :: lib_mie_ss_hf_destructor
    public :: lib_mie_ss_hf_get_n_c
!    public :: lib_mie_ss_hf_calc_triple_sum
    public :: lib_mie_ss_hf_get_coefficient_a_nm_b_nm
    public :: lib_mie_ss_hf_test_convergence_core

    public :: lib_mie_ss_helper_functions_test_functions
    public :: lib_mie_ss_helper_functions_benchmark

    ! public interfaces
    public :: lib_mie_ss_hf_init_coeff_a_n_b_n
    public :: lib_mie_ss_hf_get_coefficients_a_n_b_n

    ! --- interace ---
    interface lib_mie_ss_hf_init_coeff_a_n_b_n
        module procedure lib_mie_ss_hf_init_coeff_a_n_b_n_real
        module procedure lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx
    end interface

    interface lib_mie_ss_hf_get_coefficients_a_n_b_n
        module procedure get_coefficients_a_b_real_barberh
        module procedure get_coefficients_a_b_cmplx_barberh
    end interface

    interface lib_mie_ss_hf_get_An
        module procedure get_An_real
        module procedure get_An_cmplx
    end interface lib_mie_ss_hf_get_An

    ! ~~~ interface ~~~

    ! --- caching ---
    type cache_coefficients_a_b_real_barberh_type
        double precision :: x
        double precision :: m
        integer :: n_max
        complex(kind=8), dimension(:), allocatable :: a_n
        complex(kind=8), dimension(:), allocatable :: b_n
    end type

    type cache_coefficients_a_b_cmplx_barberh_type
        double precision :: x
        complex(kind=8) :: m
        integer :: n_max
        complex(kind=8), dimension(:), allocatable :: a_n
        complex(kind=8), dimension(:), allocatable :: b_n
    end type

    type(cache_coefficients_a_b_real_barberh_type), dimension(:), allocatable :: cache_coefficients_a_b_real_barberh
    logical :: cache_coefficients_a_b_real_barberh_enabled = .false.

    type(cache_coefficients_a_b_cmplx_barberh_type), dimension(:), allocatable :: cache_coefficients_a_b_cmplx_barberh
    logical :: cache_coefficients_a_b_cmplx_barberh_enabled = .false.
    ! ~~~ caching ~~~

    contains

!        ! Initializes the caching where the size of n_max_ab and n_max_pq represents the number of cached values.
!        !
!        ! Argument
!        ! ----
!        !   n_max_ab: integer, dimension(:)
!        !       max degree of the a b coefficients
!        !   x: double precision, dimension(size(n_max_ab))
!        !       size parameter: x = k*r = 2 * PI * N * r / lambda
!        !       k: wavenumber
!        !       r: distance
!        !       N: refractive index of the medium
!        !       lambda: wave length
!        !   m: double complex, dimension(size(n_max_ab))
!        !       relative refractive index: m = N_1 / N
!        !       N_1: refractive index of the particle [double precision, double complex]
!        !       N: refractive index of the medium [double precision]
!        !   n_max_pq: integer, dimension(:)
!        !       max degree of the p q coefficients
!        !   alpha: double precision, dimension(size(n_max_pq)), optional(std: 0)
!        !       incident angle with respect to the z-axis
!        !       codomain: 0..Pi
!        !   beta: double precision, dimension(size(n_max_pq)), optional(std: 0)
!        !       angle between the x axis and the projection of the wave vector on the x-y plane
!        !       codomain: 0..2Pi
!        subroutine lib_mie_ss_hf_contructor(n_max_ab, x, m, &
!                                         n_max_pq, alpha, beta)
!            implicit none
!            ! dummy
!            integer, dimension(:), intent(in) :: n_max_ab
!            double precision, dimension(size(n_max_ab)), intent(in) :: x
!            double complex, dimension(size(n_max_ab)), intent(in) :: m
!            integer, dimension(:), intent(in) :: n_max_pq
!            double precision, dimension(size(n_max_pq)), intent(in), optional :: alpha
!            double precision, dimension(size(n_max_pq)), intent(in), optional :: beta
!
!            ! auxiliary
!            integer :: i
!
!            integer :: counter_real
!            integer :: counter_cmplx
!
!            integer, dimension(size(n_max_ab)) :: n_max_ab_real
!            double precision, dimension(size(n_max_ab)) :: x_real
!            double precision, dimension(size(n_max_ab)) :: m_real
!
!            integer, dimension(size(n_max_ab)) :: n_max_ab_cmplx
!            double precision, dimension(size(n_max_ab)) :: x_cmplx
!            double complex, dimension(size(n_max_ab)) :: m_cmplx
!
!            double precision, dimension(size(n_max_pq)) :: m_alpha
!            double precision, dimension(size(n_max_pq)) :: m_beta
!
!            if (present(alpha)) then
!                m_alpha = alpha
!            else
!                m_alpha = 0
!            end if
!
!            if (present(beta)) then
!                m_beta = beta
!            else
!                m_beta = 0
!            end if
!
!            ! separate lists by pure real and complex "m"
!            counter_real = 0
!            counter_cmplx = 0
!            do i=1, size(n_max_ab)
!                if (aimag(m(i)) .eq. 0D0) then
!                    counter_real = counter_real + 1
!                    n_max_ab_real(counter_real) = n_max_ab(i)
!                    x_real(counter_real) = x(i)
!                    m_real(counter_real) = real(m(i))
!                else
!                    counter_cmplx = counter_cmplx + 1
!                    n_max_ab_cmplx(counter_cmplx) = n_max_ab(i)
!                    x_cmplx(counter_cmplx) = x(i)
!                    m_cmplx(counter_cmplx) = m(i)
!                end if
!            end do
!
!            if (counter_real .gt. 0) then
!                call lib_mie_ss_hf_init_coeff_a_n_b_n_real(x_real(1:counter_real), &
!                                                        m_real(1:counter_real), &
!                                                        n_max_ab_real(1:counter_real))
!            end if
!
!            if (counter_cmplx .gt. 0) then
!                call lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx(x_cmplx(1:counter_cmplx), &
!                                                         m_cmplx(1:counter_cmplx), &
!                                                         n_max_ab_cmplx(1:counter_cmplx))
!            end if
!
!            call lib_mie_ss_hf_init_coeff_p0_q0(m_alpha, m_beta, n_max_pq)
!
!        end subroutine lib_mie_ss_hf_contructor

        ! Argument
        ! ----
        !   x: double precision, dimension(:)
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: real, dimension(size(x))
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   n_max: integer, dimension(size(:))
        !       maximum degree of the polynomials
        subroutine lib_mie_ss_hf_init_coeff_a_n_b_n_real(x, m, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: x
            double precision, dimension(size(x)), intent(in) :: m
            integer, dimension(size(x)), intent(in) :: n_max

            ! auxiliary
            integer :: i

            complex(kind=8), dimension(:), allocatable :: a_n
            complex(kind=8), dimension(:), allocatable :: b_n

            ! --- init: cache_coefficients_a_b_real_barberh_x ---
            if (allocated(cache_coefficients_a_b_real_barberh)) then
                deallocate(cache_coefficients_a_b_real_barberh)
            end if

            allocate(cache_coefficients_a_b_real_barberh(size(x)))
            cache_coefficients_a_b_real_barberh_enabled = .false.

            do i=1, size(x)
                cache_coefficients_a_b_real_barberh%x = x(i)
                cache_coefficients_a_b_real_barberh%m = m(i)
                cache_coefficients_a_b_real_barberh%n_max = n_max(i)

                allocate(a_n(n_max(i)))
                allocate(b_n(n_max(i)))

                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x(i), m(i), (/ 1, n_max(i) /), a_n, b_n)

                cache_coefficients_a_b_real_barberh(i)%a_n = a_n
                cache_coefficients_a_b_real_barberh(i)%b_n = b_n

                deallocate(a_n)
                deallocate(b_n)
            end do

            cache_coefficients_a_b_real_barberh_enabled = .true.

        end subroutine lib_mie_ss_hf_init_coeff_a_n_b_n_real

        ! Argument
        ! ----
        !   x: double precision
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   n_max: integer, dimension(:)
        !       maximum degree of the polynomials
        subroutine lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx(x, m, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: x
            complex(kind=8), dimension(size(x)), intent(in) :: m
            integer, dimension(size(x)), intent(in) :: n_max

            ! auxiliary
            integer :: i

            complex(kind=8), dimension(:), allocatable :: a_n
            complex(kind=8), dimension(:), allocatable :: b_n

            ! --- init: cache_coefficients_a_b_cmplx_barberh_x ---
            if (allocated(cache_coefficients_a_b_cmplx_barberh)) then
                deallocate(cache_coefficients_a_b_cmplx_barberh)
            end if

            allocate(cache_coefficients_a_b_cmplx_barberh(size(x)))
            cache_coefficients_a_b_cmplx_barberh_enabled = .false.

            do i=1, size(x)
                cache_coefficients_a_b_cmplx_barberh%x = x(i)
                cache_coefficients_a_b_cmplx_barberh%m = m(i)
                cache_coefficients_a_b_cmplx_barberh%n_max = n_max(i)

                allocate(a_n(n_max(i)))
                allocate(b_n(n_max(i)))

                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x(i), m(i), (/ 1, n_max(i) /), a_n, b_n)

                cache_coefficients_a_b_cmplx_barberh(i)%a_n = a_n
                cache_coefficients_a_b_cmplx_barberh(i)%b_n = b_n

                deallocate(a_n)
                deallocate(b_n)
            end do

            cache_coefficients_a_b_cmplx_barberh_enabled = .true.

        end subroutine lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx



        subroutine lib_mie_ss_hf_destructor
            implicit none

            ! --- deallocate caching ---
            if (allocated(cache_coefficients_a_b_real_barberh)) then
                deallocate(cache_coefficients_a_b_real_barberh)
            end if
            cache_coefficients_a_b_real_barberh_enabled = .false.

            if (allocated(cache_coefficients_a_b_cmplx_barberh)) then
                deallocate(cache_coefficients_a_b_cmplx_barberh)
            end if
            cache_coefficients_a_b_cmplx_barberh_enabled = .false.
        end subroutine

        ! calculates the scattering coefficients
        !
        ! HINT: The refractive index has only a real value.
        !
        ! Arguments
        ! ----
        !   x: double precision
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: double precision
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   mu: double precision
        !       permeability of the medium
        !   mu1: double precision
        !       permeability of the particle
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Absorption and Scattering of Light by Small Particles,Bohren and Huffman,  eq. (4.53)
        subroutine get_coefficients_a_b_real_bohrenh(x, m, mu, mu1, n, a_n, b_n)
            implicit none
            ! dummy
            double precision :: x
            double precision :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            complex(kind=8), dimension(n(2)-n(1)+1) :: numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: denominator

            double precision, dimension(n(2)-n(1)+1) :: j_n_x
            double precision, dimension(n(2)-n(1)+1) :: j_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_n_x
            double precision, dimension(n(2)-n(1)+1) :: s_dn_x
            double precision, dimension(n(2)-n(1)+1) :: s_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_dn_mx

            complex(kind=8), dimension(n(2)-n(1)+1) :: h_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_dn_x

            double precision :: mx

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            s_dn_x = lib_math_riccati_s_derivative(x, n(1), number_of_n, s_n_x)
            j_n_x = s_n_x / x

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx) * m
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x

            numerator = cmplx(mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx, 0.0, kind=8)
            denominator = cmplx(mu * m*m * j_n_mx, 0.0, kind=8) * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = cmplx(mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx, 0.0, kind=8)
            denominator = cmplx(mu1 * j_n_mx, 0.0, kind=8) * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_real_bohrenh

        ! calculates the scattering coefficients
        !
        ! Arguments
        ! ----
        !   x: complex
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   mu: double precision
        !       permeability of the medium
        !   mu1: double precision
        !       permeability of the particle
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. (4.53)
        subroutine get_coefficients_a_b_cmplx_bohrenh(x, m, mu, mu1, n, a_n, b_n)
            implicit none
            ! dummy
            complex(kind=8) :: x
            complex(kind=8) :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            complex(kind=8), dimension(n(2)-n(1)+1) :: numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+1) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: j_n_mx
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_dn_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_n_mx
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_dn_mx

            complex(kind=8), dimension(n(2)-n(1)+1) :: h_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_dn_x

            complex(kind=8) :: mx

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            s_dn_x = lib_math_riccati_s_derivative(x, n(1), number_of_n, s_n_x)
            j_n_x = s_n_x / x

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx) * m
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x

            ! --- test ---

            numerator = mu * m*m * j_n_mx * s_dn_x
            numerator = - mu1 * j_n_x * s_dn_mx

            denominator = mu * m*m * j_n_mx * xi_dn_x
            denominator = - mu1 * h_n_x * s_dn_mx

            ! ~~~ test ~~~
            numerator = mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx
            denominator = mu * m*m * j_n_mx * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx
            denominator = mu1 * j_n_mx * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_cmplx_bohrenh

        ! calculates the scattering coefficients
        !
        ! Arguments
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: real
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !   caching: logical (std: .true.)
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n, caching)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: m
            integer(kind=4), dimension(2), intent(in) :: n

            complex(kind=8), dimension(n(2)-n(1)+1), intent(inout) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1), intent(inout) :: b_n

            logical, intent(in), optional :: caching

            ! auxiliary
            integer :: i
            complex(kind=8) :: numerator
            complex(kind=8) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+2) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+2) :: h_n_x

            real(kind=8), dimension(n(2)-n(1)+1) :: An

            real(kind=8) :: mx
            real(kind=8) :: n_div_x
            real(kind=8) :: mAn
            real(kind=8) :: An_div_m
            real(kind=8) :: buffer

            integer(kind=4) :: number_of_n

            integer :: cache_no

            ! 0: calculate
            ! >0: array number
            ! -1: search
            cache_no = 0

            if (cache_coefficients_a_b_real_barberh_enabled) then
                if (present(caching)) then
                    if (caching) then
                        cache_no = -1
                    end if
                else
                    cache_no = -1
                end if
            end if

            ! search
            if (cache_no .eq. -1) then
                do i=1, size(cache_coefficients_a_b_real_barberh)
                    if ((cache_coefficients_a_b_real_barberh(i)%x .eq. x) &
                        .and. (cache_coefficients_a_b_real_barberh(i)%m .eq. m) &
                        .and. (cache_coefficients_a_b_real_barberh(i)%n_max .ge. n(2)) ) then
                        cache_no = i
                        exit
                    end if
                end do
            end if

            if (cache_no .gt. 0) then
                a_n = cache_coefficients_a_b_real_barberh(cache_no)%a_n(n(1):n(2))
                b_n = cache_coefficients_a_b_real_barberh(cache_no)%b_n(n(1):n(2))
            else

                number_of_n = n(2) - n(1) + 1

                mx = m*x

                j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
                h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

                An = lib_mie_ss_hf_get_An(x, m, n(1), number_of_n)

                !$OMP PARALLEL DO PRIVATE(i, mAn, n_div_x, buffer, numerator, denominator)
                do i=1, number_of_n
                    n_div_x = real(n(1)+i-1, kind=8) / x
                    mAn = m * An(i)
                    An_div_m = An(i) / m

                    buffer = An_div_m + n_div_x
                    numerator = buffer * j_n_x(i+1) - j_n_x(i)
                    denominator = buffer * h_n_x(i+1) - h_n_x(i)

                    a_n(i) = numerator / denominator

                    buffer = mAn + n_div_x
                    numerator = buffer * j_n_x(i+1) - j_n_x(i)
                    denominator = buffer * h_n_x(i+1) - h_n_x(i)

                    b_n(i) = numerator / denominator
                end do
                !$OMP END PARALLEL DO

            end if

        end subroutine get_coefficients_a_b_real_barberh

        ! calculates the scattering coefficients
        !
        ! Arguments
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n, caching)
            implicit none
            ! dummy
            double precision :: x
            complex(kind=8) :: m
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            logical, intent(in), optional :: caching

            ! auxiliary
            integer :: i
            complex(kind=8) :: numerator
            complex(kind=8) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+2) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+2) :: h_n_x

            complex(kind=8), dimension(n(2)-n(1)+1) :: An

            complex(kind=8) :: mx
            complex(kind=8) :: n_div_x
            complex(kind=8) :: mAn
            complex(kind=8) :: An_div_m
            complex(kind=8) :: buffer

            integer(kind=4) :: number_of_n

            integer :: cache_no

            ! 0: calculate
            ! >0: array number
            ! -1: search
            cache_no = 0

            if (cache_coefficients_a_b_cmplx_barberh_enabled) then
                if (present(caching)) then
                    if (caching) then
                        cache_no = -1
                    end if
                else
                    cache_no = -1
                end if
            end if

            ! search
            if (cache_no .eq. -1) then
                do i=1, size(cache_coefficients_a_b_cmplx_barberh)
                    if ((cache_coefficients_a_b_cmplx_barberh(i)%x .eq. x) &
                        .and. (cache_coefficients_a_b_cmplx_barberh(i)%m .eq. m) &
                        .and. (cache_coefficients_a_b_cmplx_barberh(i)%n_max .ge. n(2)) ) then
                        cache_no = i
                        exit
                    end if
                end do
            end if

            if (cache_no .gt. 0) then
                a_n = cache_coefficients_a_b_cmplx_barberh(cache_no)%a_n(n(1):n(2))
                b_n = cache_coefficients_a_b_cmplx_barberh(cache_no)%b_n(n(1):n(2))
            else

                number_of_n = n(2) - n(1) + 1

                mx = m*x

                j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
                h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

                An = lib_mie_ss_hf_get_An(x, m, n(1), number_of_n)

                !$OMP PARALLEL DO PRIVATE(i, mAn, n_div_x, An_div_m, buffer, numerator, denominator)
                do i=1, number_of_n
                    n_div_x = real(n(1)+i-1, kind=8) / x
                    mAn = m * An(i)
                    An_div_m = An(i) / m

                    buffer = An_div_m + n_div_x
                    numerator = buffer * j_n_x(i+1) - j_n_x(i)
                    denominator = buffer * h_n_x(i+1) - h_n_x(i)

                    a_n(i) = numerator / denominator

                    buffer = mAn + n_div_x
                    numerator = buffer * j_n_x(i+1) - j_n_x(i)
                    denominator = buffer * h_n_x(i+1) - h_n_x(i)

                    b_n(i) = numerator / denominator
                end do
                !$OMP END PARALLEL DO
            end if

        end subroutine get_coefficients_a_b_cmplx_barberh

        ! Calculates the logarithmic derivative An
        !
        ! Formula: A_n = [mx j_n(mx)]' / (mx j_n(mx))
        !
        ! Argument
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: real
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   fnu: integer
        !       order of initial function, fnu.GE.0
        !   n: integer
        !       number of members of the sequence, n.GE.1
        !
        ! Returns
        ! ----
        !   An: complex(fnu:fnu + n - 1)
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_real(x, m, fnu, n) result (An)
            implicit none
            ! dummy
            real(kind=8) :: x
            real(kind=8) :: m
            integer :: fnu
            integer :: n

            real(kind=8), dimension(fnu:fnu+n-1) :: An

            ! parameter

            ! auxiliary
            integer :: i
            integer :: m_n_max
            integer :: m_n_mx
            real(kind=8) :: m_mx
            real(kind=8) :: m_buffer


            m_mx = m * x
            m_n_max = fnu + n - 1

            ! eq. (4.20)
            m_n_mx = max(lib_mie_ss_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

            m_buffer = 0
            do i=m_n_mx, fnu+1, -1
                m_buffer = get_An_minus_1_real(i, m_mx, m_buffer)

                if (i .le. m_n_max+1) then
                    An(i-1) = m_buffer
                end if
            end do

        end function get_An_real

        ! Calculates the logarithmic derivative An
        !
        ! Formula: A_n = [mx j_n(mx)]' / (mx j_n(mx))
        !
        ! Argument
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   fnu: integer
        !       order of initial function, fnu.GE.0
        !   n: integer
        !       number of members of the sequence, n.GE.1
        !
        ! Returns
        ! ----
        !   An: complex(fnu:fnu + n - 1)
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_cmplx(x, m, fnu, n) result (An)
            implicit none
            ! dummy
            double precision :: x
            complex(kind=8) :: m
            integer :: fnu
            integer :: n

            complex(kind=8), dimension(fnu:fnu+n-1) :: An

            ! parameter

            ! auxiliary
            integer :: i
            integer :: m_n_max
            integer :: m_n_mx
            complex(kind=8) :: m_mx
            complex(kind=8) :: m_buffer


            m_mx = m * x
            m_n_max = fnu + n - 1

            ! eq. (4.20)
            m_n_mx = max(lib_mie_ss_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

            m_buffer = cmplx(0.0, 0.0, kind=8)
            do i=m_n_mx, fnu+1, -1
                m_buffer = get_An_minus_1_cmplx(i, m_mx, m_buffer)

                if (i .le. m_n_max+1) then
                    An(i-1) = m_buffer
                end if
            end do

        end function get_An_cmplx

        ! Argument
        ! ----
        !   n: integer
        !       degree of the logarithmic derivative
        !   mx: complex
        !       x: real
        !           size parameter: x = k*r = 2 * PI * N * r / lambda
        !           k: wavenumber
        !           r: distance
        !           N: refractive index of the medium
        !           lambda: wave length
        !       m: real
        !           relative refractive index: m = N_1 / N
        !           N_1: refractive index of the particle
        !           N: refractive index of the medium
        !   An: complex
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_minus_1_real(n, mx, An) result (rv)
            implicit none
            ! dummy
            integer :: n
            real(kind=8) :: mx
            real(kind=8) :: An

            real(kind=8) :: rv

            ! auxiliary
            real(kind=8) :: n_div_mx

            n_div_mx = real(n, kind=8) / mx

            rv = n_div_mx - 1.0_8/(An + n_div_mx)

        end function get_An_minus_1_real

        ! Argument
        ! ----
        !   n: integer
        !       degree of the logarithmic derivative
        !   mx: complex
        !       x: complex
        !           size parameter: x = k*r = 2 * PI * N * r / lambda
        !           k: wavenumber
        !           r: distance
        !           N: refractive index of the medium
        !           lambda: wave length
        !       m: complex
        !           relative refractive index: m = N_1 / N
        !           N_1: refractive index of the particle
        !           N: refractive index of the medium
        !   An: complex
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_minus_1_cmplx(n, mx, An) result (rv)
            implicit none
            ! dummy
            integer :: n
            complex(kind=8) :: mx
            complex(kind=8) :: An

            complex(kind=8) :: rv

            ! auxiliary
            complex(kind=8) :: n_div_mx

            n_div_mx = real(n, kind=8) / mx

            rv = n_div_mx - 1.0_8/(An + n_div_mx)

        end function get_An_minus_1_cmplx

        ! calculates highest degree n for a convergent algorithm
        !
        ! Argument
        ! ----
        !   x: double precision
        !       size parameter
        !       x = r * 2 * Pi / lambda
        !       r: radius
        !       lambda: wavelength
        !
        ! Reference: Light Scattering by Particles: Computational Methods, PW Barber, S C Hill, eq. 4.16
        function lib_mie_ss_hf_get_n_c(x) result (rv)
            implicit none
            ! dummy
            double precision :: x

            integer :: rv

            ! auxiliary
            double precision :: dummy

            dummy = x + 4.05_8 * x**(1.0_8/3.0_8) + 2.0_8

            rv = int(ceiling(dummy))

        end function lib_mie_ss_hf_get_n_c

!        ! Arguments
!        ! ----
!        !   simulation_parameter: type(simulation_parameter_type)
!        !   sphere: type(sphere_type), dimension(:)
!        !       list of spheres
!        !   sphere_parameter: type(sphere_parameter_type), dimension(:)
!        !       list of shared sphere parameters
!        !   z_selector: integer
!        !       parameter of the spherical harmonics
!        !       values:
!        !           1: spherical Bessel function first kind   j_n
!        !           2: spherical Bessel function second kind  y_n
!        !           3: spherical Hankel function first kind   h^(1)_n
!        !           4: spherical Hankel function second kind  h^(2)_n
!        !   f: double precision, optional (std: 1)
!        !       numerical factor (0, 1]
!        !       "In our actual calculations, some multi-
!        !        sphere systems do not converge if f 5 1, but they do
!        !        converge when the value of f is reduced to, say, 0.7."[1]
!        !
!        ! Returns
!        ! ----
!        !   a_j_nm: type(list_list_cmplx)
!        !   b_j_nm: type(list_list_cmplx)
!        !
!        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30, 35
!        subroutine lib_mie_ss_hf_calc_triple_sum(simulation_parameter, &
!                                              sphere, sphere_parameter, sphere_j, &
!                                              z_selector, &
!                                              a_j_nm, b_j_nm, &
!                                              f)
!            implicit none
!            ! dummy
!            type(lib_mie_simulation_parameter_type) :: simulation_parameter
!            type(lib_mie_sphere_type), dimension(:), intent(in) :: sphere
!            type(lib_mie_sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
!            integer :: sphere_j
!            integer(kind=1) :: z_selector
!
!            type(list_list_cmplx), intent(inout) :: a_j_nm
!            type(list_list_cmplx), intent(inout) :: b_j_nm
!
!            double precision, intent(in), optional :: f
!
!            ! auxiliary
!            integer :: n
!            integer :: nu
!
!            integer :: l
!            integer :: i
!
!            double precision :: lambda
!            double precision :: n_medium
!            type(cartesian_coordinate_real_type) :: k
!
!            integer, dimension(2) :: n_range
!            integer, dimension(2) :: nu_range
!
!            type(list_list_cmplx) :: p
!            type(list_list_cmplx) :: q
!            type(list_cmplx) :: a_j_n
!            type(list_cmplx) :: b_j_n
!            type(list_cmplx) :: a_l_nu
!            type(list_cmplx) :: b_l_nu
!            type(list_list_cmplx) :: a_l_numu
!            type(list_list_cmplx) :: b_l_numu
!            type(list_4_cmplx) :: A_nmnumu
!            type(list_4_cmplx) :: B_nmnumu
!            type(cartesian_coordinate_real_type) :: d_0_j
!            type(cartesian_coordinate_real_type) :: d_0_l
!            type(cartesian_coordinate_real_type) :: x
!
!            type(list_list_cmplx) :: buffer_sum_a
!            type(list_list_cmplx) :: buffer_sum_b
!            type(list_list_cmplx), dimension(:), allocatable :: buffer_sum_a_l
!            type(list_list_cmplx), dimension(:), allocatable :: buffer_sum_b_l
!
!            type(list_list_cmplx) :: a_j_nm_old
!            type(list_list_cmplx) :: b_j_nm_old
!
!            i = sphere(sphere_j)%sphere_parameter_index
!            n_range = sphere_parameter(i)%n_range
!            d_0_j = sphere(sphere_j)%d_0_j
!            a_j_n = sphere_parameter(i)%a_n
!            b_j_n = sphere_parameter(i)%b_n
!
!            if (present(f)) then
!                a_j_nm_old = sphere(sphere_j)%a_nm
!                b_j_nm_old = sphere(sphere_j)%b_nm
!            end if
!
!            lambda = simulation_parameter%illumination%lambda_0
!            n_medium = simulation_parameter%refractive_index_medium
!            k = simulation_parameter%illumination%wave_vector_0 * n_medium
!
!
!            allocate(buffer_sum_a_l(lbound(sphere, 1):ubound(sphere, 1)))
!            allocate(buffer_sum_b_l(lbound(sphere, 1):ubound(sphere, 1)))
!
!            ! calculate the triple sum
!            !$OMP PARALLEL DO PRIVATE(l, i, nu_range, a_l_nu, b_l_nu, d_0_l, &
!            !$OMP&  a_nmnumu, b_nmnumu, x, p, q, nu)
!            do l=ubound(sphere, 1), lbound(sphere, 1)
!                if (l .ne. sphere_j) then
!                    i = sphere(l)%sphere_parameter_index
!                    nu_range = sphere_parameter(i)%n_range
!                    a_l_numu = sphere(l)%a_nm
!                    b_l_numu = sphere(l)%b_nm
!                    d_0_l = sphere(l)%d_0_j
!
!                    x = (d_0_l - d_0_j) * n_medium / lambda
!
!                    call lib_mie_vector_spherical_harmonics_translation_coefficient(x, &
!                                                                                    n_range, nu_range, z_selector,&
!                                                                                    a_nmnumu, b_nmnumu)
!
!
!                    call calc_inner_2sum(a_l_numu, b_l_numu, a_nmnumu, b_nmnumu, &
!                                         buffer_sum_a_l(l), buffer_sum_b_l(l))
!
!                end if
!            end do
!            !$OMP END PARALLEL DO
!
!            call init_list(buffer_sum_a, n_range(1), n_range(2) - n_range(1) + 1, cmplx(0, 0, kind=lib_math_type_kind))
!            call init_list(buffer_sum_b, n_range(1), n_range(2) - n_range(1) + 1, cmplx(0, 0, kind=lib_math_type_kind))
!            ! final summation of the triple sum
!            do l=lbound(sphere, 1), ubound(sphere, 1)
!                buffer_sum_a = buffer_sum_a + buffer_sum_a_l(l)
!                buffer_sum_b = buffer_sum_b + buffer_sum_b_l(l)
!            end do
!
!            deallocate(buffer_sum_a_l)
!            deallocate(buffer_sum_b_l)
!
!
!            call lib_mie_ss_hf_get_p_q_j_j(k, d_0_j, n_range, p, q)
!
!            call init_list(a_j_nm, n_range(1), n_range(2) - n_range(1) + 1)
!            call init_list(b_j_nm, n_range(1), n_range(2) - n_range(1) + 1)
!
!            a_j_nm = p - buffer_sum_a
!            b_j_nm = q - buffer_sum_b
!
!            !$OMP PARALLEL DO PRIVATE(n)
!            do n=n_range(1), n_range(2)
!                a_j_nm%item(n)%item(:) = a_j_n%item(n) * a_j_nm%item(n)%item(:)
!                b_j_nm%item(n)%item(:) = b_j_n%item(n) * b_j_nm%item(n)%item(:)
!            end do
!            !$OMP END PARALLEL DO
!
!            if (present(f)) then
!                a_j_nm = (1.0d0 - f) * a_j_nm_old + f * a_j_nm
!                b_j_nm = (1.0d0 - f) * b_j_nm_old + f * b_j_nm
!            end if
!
!        end subroutine
!
!        ! Calculates the two inner sums of eq. 30 (interactive scattering coefficients)
!        !
!        ! Argument
!        ! ----
!        !   a_l_munu: type(list_list_cmplx)
!        !       Mie coefficient of the l-th sphere
!        !       restriction: same dimension as b_l_munu
!        !   b_l_munu: type(list_list_cmplx)
!        !       Mie coefficient of the l-th sphere
!        !       restriction: same dimension as a_l_munu
!        !   a_munumn: type(list_4_cmplx)
!        !       translation coefficient
!        !       restriction: same dimension as b_munumn
!        !   b_munumn: type(list_4_cmplx)
!        !       translation coefficient
!        !       restriction: same dimension as a_munumn
!        !
!        ! Returns
!        ! ----
!        !   a_2sum: double complex
!        !
!        !
!        ! Reference:Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
!        subroutine calc_inner_2sum(a_l_munu, b_l_munu, a_munumn, b_munumn, &
!                                   a_2sum, b_2sum)
!            implicit none
!            ! dummy
!            type(list_list_cmplx), intent(in) :: a_l_munu
!            type(list_list_cmplx), intent(in) :: b_l_munu
!            type(list_4_cmplx), intent(in) :: a_munumn
!            type(list_4_cmplx), intent(in) :: b_munumn
!
!            type(list_list_cmplx), intent(inout) :: a_2sum
!            type(list_list_cmplx), intent(inout) :: b_2sum
!
!            ! auxiliary
!            integer :: n
!            integer :: m
!            integer :: nu
!            integer :: mu
!
!            integer, dimension(2) :: n_range
!            integer, dimension(2) :: nu_range
!
!            double complex :: buffer_cmplx
!
!            ! set boundaries
!            n_range(1) = lbound(a_munumn%item, 1)
!            n_range(2) = ubound(a_munumn%item, 1)
!            nu_range(1) = max(lbound(a_munumn%item(n_range(1))%item(0)%item, 1), 1)
!            nu_range(2) = ubound(a_munumn%item(n_range(1))%item(0)%item, 1)
!
!            call init_list(a_2sum, n_range(1), n_range(2)-n_range(1)+1)
!            call init_list(b_2sum, n_range(1), n_range(2)-n_range(1)+1)
!
!            ! calculate the summand
!            !$OMP PARALLEL DO PRIVATE(n, m)
!            do n=n_range(1), n_range(2)
!                do m=-n, n
!                    !$OMP PARALLEL DO PRIVATE(nu, mu, buffer_cmplx)
!                    do nu=nu_range(1), nu_range(2)
!                        do mu=-nu, nu
!                            ! summand of a_j_mn
!                            ! first line eq. 30
!                            buffer_cmplx = a_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
!                            buffer_cmplx = buffer_cmplx &
!                                           + b_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
!                            a_2sum%item(n)%item(m) = buffer_cmplx
!
!                            ! summand of b_j_mn
!                            ! second line eq. 30
!                            buffer_cmplx = a_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
!                            buffer_cmplx = buffer_cmplx &
!                                           + b_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
!                            b_2sum%item(n)%item(m) = buffer_cmplx
!                        end do
!                    end do
!                    !$OMP END PARALLEL DO
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!        end subroutine calc_inner_2sum

        ! Argument
        ! ----
        !   x: double precision
        !       size parameter
        !   m: double complex
        !       relative refractive index
        !       m = n_particle / n_medium
        !   p_nm: type(list_list_cmplx)
        !       illuminiation coefficient
        !   q_nm: type(list_list_cmplx)
        !       illuminiation coefficient
        !
        ! Returns
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       interactive scattering coefficient
        !   b_nm: type(list_list_cmplx)
        !       interactive scattering coefficient
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 12
        subroutine lib_mie_ss_hf_get_coefficient_a_nm_b_nm(x, n_particle, n_medium, p_nm, q_nm, a_nm, b_nm)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            double complex, intent(in) :: n_particle
            double precision, intent(in) :: n_medium
            type(list_list_cmplx), intent(in) :: p_nm
            type(list_list_cmplx), intent(in) :: q_nm

            type(list_list_cmplx) :: a_nm
            type(list_list_cmplx) :: b_nm

            ! auxiliary
            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n
            integer, dimension(2) :: n_range

            n_range(1) = lbound(p_nm%item, 1)
            n_range(2) = ubound(p_nm%item, 1)

            call init_list(a_n, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(b_n, n_range(1), n_range(2) - n_range(1) + 1)

            if (aimag(n_particle) .eq. 0) then
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x, real(n_particle) / n_medium, n_range, &
                                                         a_n%item, b_n%item)
            else
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x, n_particle / n_medium, n_range, &
                                                         a_n%item, b_n%item)
            end if


            call init_list(a_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(b_nm, n_range(1), n_range(2)-n_range(1)+1)

            a_nm = a_n * p_nm
            b_nm = b_n * q_nm
        end subroutine lib_mie_ss_hf_get_coefficient_a_nm_b_nm



        ! Argument
        ! ----
        !   p_nm: type(list_list_cmplx)
        !       illumination coefficient
        !   q_nm: type(list_list_cmplx)
        !       illumination coefficient
        !   a_nm: type(list_list_cmplx)
        !       scattering coefficient
        !   b_nm: type(list_list_cmplx)
        !       scatteirng coefficient
        !   n_c: integer
        !       calculated highest degree n for a convergent algorithm
        !
        ! Returns
        ! ----
        !   rv: integer, dimension(2)
        !       rv(1):
        !           -1: series may converge or diverge
        !            0: series diverges
        !           >0: series converges absolutely (degree n = rv(1))
        !       rv(2): relative error [1/1000]
        !           >=0: abs( cross_section(n=rv(1)) - cross_section(n=n_max) )**" / abs(cross_section(n=n_max))**2
        !            -1: rv(1) .le. 0
        function lib_mie_ss_hf_test_convergence_core(p_nm, q_nm, a_nm, b_nm, n_c) result(rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm

            type(list_list_cmplx), intent(in) :: p_nm
            type(list_list_cmplx), intent(in) :: q_nm

            integer, intent(in) :: n_c

            integer, dimension(2) :: rv

            ! auxiliary
            integer :: n
            integer :: m
            integer(kind=4), dimension(2) :: n_range

            type(list_list_real) :: summand_sca
            type(list_list_real) :: summand_ext

            double precision, dimension(:), allocatable :: c_sca_n
            double precision, dimension(:), allocatable :: c_ext_n
            double precision :: buffer_sca
            double precision :: buffer_ext

            integer :: r_sca
            integer :: r_ext

            n_range(1) = lbound(p_nm%item, 1)
            n_range(2) = ubound(p_nm%item, 1)

            allocate( c_sca_n(n_range(1):n_range(2)) )
            allocate( c_ext_n(n_range(1):n_range(2)) )

            call get_cross_section_core(p_nm, q_nm, a_nm, b_nm, &
                                        summand_sca, summand_ext)

            !$OMP PARALLEL DO PRIVATE(n, m, buffer_sca, buffer_ext)
            do n = n_range(1), n_range(2)
                buffer_sca = 0
                buffer_ext = 0
                do m = -n, n
                    buffer_sca = buffer_sca + summand_sca%item(n)%item(m)
                    buffer_ext = buffer_ext + summand_sca%item(n)%item(m)
                end do
                c_sca_n(n) = buffer_sca
                c_ext_n(n) = buffer_ext
            end do
            !$OMP END PARALLEL DO

            ! test with n_c
            r_sca = lib_math_convergence_root_test(c_sca_n(n_range(1):n_c))
            r_ext = lib_math_convergence_root_test(c_ext_n(n_range(1):n_c))

            if (r_sca .gt. 0 .and. r_ext .gt. 0) then
                rv(1) = n_c
            else
                ! test with "n_max" = 45
                r_sca = lib_math_convergence_root_test(c_sca_n)
                r_ext = lib_math_convergence_root_test(c_ext_n)

                if (r_sca .gt. 0 .and. r_ext .gt. 0) then
                    rv(1) = max(r_sca, r_ext, n_c)
                else
                    rv(1) = min(r_sca, r_ext)
                end if
            end if

            if (rv(1) .gt. 0) then
                buffer_sca = sum(c_sca_n(n_range(1):n_c))
                buffer_ext = sum(c_ext_n(n_range(1):n_c))

                buffer_sca = abs(buffer_sca - sum(c_sca_n))**2 / abs(sum(c_sca_n))**2
                buffer_ext = abs(buffer_ext - sum(c_ext_n))**2 / abs(sum(c_ext_n))**2

                rv(2) = int(ceiling(max(buffer_sca, buffer_ext) * 1000))
            else
                rv(2) = -1
            end if
        end function lib_mie_ss_hf_test_convergence_core

        ! Argument
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       scattering coefficient
        !   b_nm: type(list_list_cmplx)
        !       scattering coefficient
        !   k: double precision
        !       wave number: k = 2 PI / lambda * n_medium
        !       - lambda: wave length
        !       - n_medium: refractive index of the medium
        !
        ! Returns
        ! ----
        !   c: double precision, dimension(3)
        !       c(1): scattering cross section
        !       c(2): extinction cross section
        !       c(3): absorption cross section
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 57
        function get_cross_section(p_0_nm, q_0_nm, a_nm, b_nm, k) result (c)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: p_0_nm
            type(list_list_cmplx), intent(in) :: q_0_nm
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            double precision, intent(in) :: k

            double precision, dimension(3) :: c

            ! dummy
            integer :: n
            integer :: m

            double precision :: c_sca
            double precision :: c_ext
            double precision :: c_abs

            type(list_list_real) :: summand_sca
            type(list_list_real) :: summand_ext

            double precision :: buffer_sca
            double precision :: buffer_ext

            call get_cross_section_core(p_0_nm, q_0_nm, a_nm, b_nm, summand_sca, summand_ext)

            c_sca = 0
            c_ext = 0
            do n = lbound(a_nm%item, 1), ubound(a_nm%item, 1)
!                m=1
                buffer_sca = 0
                buffer_ext = 0
                do m = -n, n
                    buffer_sca = buffer_sca + summand_sca%item(n)%item(m)
                    buffer_ext = buffer_ext + summand_ext%item(n)%item(m)
                end do
                c_sca = c_sca + buffer_sca
                c_ext = c_ext + buffer_ext
            end do

            c_sca = c_sca * 4D0 * PI / (k**2)
            c_ext = c_ext * 4D0 * PI / (k**2)
            c_abs = c_ext - c_sca

            c(1) = c_sca
            c(2) = c_ext
            c(3) = c_abs

        end function get_cross_section

        subroutine get_cross_section_core(p_0_nm, q_0_nm, a_nm, b_nm, &
                                          summand_sca, summand_ext)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: p_0_nm
            type(list_list_cmplx), intent(in) :: q_0_nm
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
!            double precision, intent(in) :: k

            type(list_list_real) :: summand_sca
            type(list_list_real) :: summand_ext

            ! auxiliary
            integer :: n
            integer :: m
            double precision :: buffer_n
            double precision :: buffer_real

            call init_list(summand_sca, &
                           lbound(a_nm%item, 1), &
                           ubound(a_nm%item, 1) - lbound(a_nm%item, 1) + 1, &
                           0D0)
            call init_list(summand_ext, &
                           lbound(a_nm%item, 1), &
                           ubound(a_nm%item, 1) - lbound(a_nm%item, 1) + 1, &
                           0D0)

            !$OMP PARALLEL DO PRIVATE(n, m, buffer_n, buffer_real)
            do n = lbound(a_nm%item, 1), ubound(a_nm%item, 1)
                buffer_n = dble(n * (n + 1) * (2 * n + 1))
                do m = -n, n
                    buffer_real = buffer_n * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)
                    summand_sca%item(n)%item(m) = buffer_real &
                                                  * ( abs(a_nm%item(n)%item(m))**2 &
                                                      + abs(b_nm%item(n)%item(m))**2 )

                    summand_ext%item(n)%item(m) = buffer_real &
                                                  * real( conjg(p_0_nm%item(n)%item(m)) * a_nm%item(n)%item(m)  &
                                                          + conjg(q_0_nm%item(n)%item(m)) * b_nm%item(n)%item(m) )
                end do
            end do
            !$OMP END PARALLEL DO

        end subroutine get_cross_section_core

        function lib_mie_ss_helper_functions_test_functions() result(rv)
            use file_io
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_get_coefficients_a_b_real_bohrenh()) then
                rv = rv + 1
            end if
!            if (.not. test_get_coefficients_a_b_cmplx_bohrenh()) then
!                rv = rv + 1
!            end if
            if (.not. test_get_coefficients_a_b_real_barberh()) then
                rv = rv + 1
            end if
            if (.not. test_get_coefficients_a_b_cmplx_barberh()) then
                rv = rv + 1
            end if
            if (.not. test_get_cross_section()) then
                rv = rv + 1
            end if

            print *, ""
            print *, "------lib_mie_scattering_by_a_sphere_test_functions------"
            if (rv == 0) then
                print *, "lib_mie_scattering_by_a_sphere_test_functions tests: OK"
            else
                print *, rv,"lib_mie_scattering_by_a_sphere_test_functions test(s) FAILED"
            end if
            print *, "---------------------------------------------------------"
            print *, ""

            contains

            function test_get_coefficients_a_b_real_bohrenh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    real(kind=8) :: x
                    real(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    mu = 1
                    mu1 = 1

                    m = 1.5
                    x= 10

                    ground_truth_a_n(1) = cmplx(0.938111_8, +0.240954_8, kind=8)
                    ground_truth_a_n(2) = cmplx(0.962707_8, +0.189478_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.994996_8, +0.0705625_8, kind=8)
                    ground_truth_a_n(4) = cmplx(0.99737_8, -0.0512159_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.980385_8, -0.138672_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.805329_8, +0.395947_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.940931_8, -0.235754_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.999006_8, -0.0315159_8, kind=8)

                    call get_coefficients_a_b_real_bohrenh(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_real_bohrenh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_real_bohrenh

                function test_get_coefficients_a_b_cmplx_bohrenh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    complex(kind=8) :: x
                    complex(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    ! Reference: Electromagnetic scattering on spherical polydispersions,  D.Deirmendjian, p. 27
                    m = cmplx(1.28, -1.37, kind=8)
                    x= cmplx(20, 0, kind=8)
!                    ground_truth_a_n(1) = cmplx(-0.22686_8+0.5_8, -0.12863_8, kind=8)
!                    ground_truth_b_n(1) = cmplx(0.22864_8+0.5_8, 0.13377_8, kind=8)

                    ground_truth_a_n(1) = cmplx(-181.13_8, -327.306_8, kind=8)
                    ground_truth_a_n(2) = cmplx(81.2324_8, +94.3237_8, kind=8)
                    ground_truth_a_n(3) = cmplx(-51.6918_8, -32.7439, kind=8)
                    ground_truth_a_n(4) = cmplx(36.652_8, +5.67418_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.367587_8, -0.463775_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.722992_8, +0.427339_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.159304_8, -0.340386_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.947238_8, +0.177162_8, kind=8)


                    call get_coefficients_a_b_cmplx_bohrenh(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_cmplx_bohrenh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_cmplx_bohrenh

                function test_get_coefficients_a_b_real_barberh() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    real(kind=8) :: x
                    real(kind=8) :: m
!                    double precision :: mu
!                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n


                    m = 1.5_8
                    x = 10.0_8
                    ! Reference: Light Scattering by Particles: Computational Methods, Barber Hill, program S2
                    ! WARNING: same algorithm
                    ground_truth_a_n(1) = cmplx(0.825333297_8, 0.379681736_8, kind=8)
                    ground_truth_a_n(2) = cmplx(0.999948084_8, 0.0072028893_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.970794678_8, 0.168381661_8, kind=8)
                    ground_truth_a_n(4) = cmplx(0.995235085_8, -0.068863928_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.997406423_8,0.0508609638_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.885268927_8,0.318697095_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.995330393_8,-0.0681746677_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.998444855_8,-0.0394046269_8, kind=8)

                    call get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n)

                    open(unit=99, file="temp/test_get_coefficients_a_b_real_barberh.csv", status='unknown')
                    rv = write_csv(99, (/ a_n, b_n /), (/ "n_", "ab" /))
                    close(99)

                    rv = .true.
                     print *, "test_get_coefficients_a_b_real_barberh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) / ground_truth_a_n(i) - 1d0)
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) / ground_truth_b_n(i) - 1d0)
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_real_barberh

                function test_get_coefficients_a_b_cmplx_barberh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    double precision :: x
                    complex(kind=8) :: m
!                    double precision :: mu
!                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    m = cmplx(1.28, -1.37, kind=8)
                    x= 10.0_8

                    ! Reference: Light Scattering by Particles: Computational Methods, Barber Hill, program S2
                    ! WARNING: same algorithm
                    ground_truth_a_n(1) = cmplx(-0.333899468_8, 0.472839594_8, kind=8)
                    ground_truth_a_n(2) = cmplx(1.1025393_8, -0.765782595_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.419178337_8, 0.996963501_8, kind=8)
                    ground_truth_a_n(4) = cmplx(-0.174542636_8, -0.790816486_8, kind=8)

                    ground_truth_b_n(1) = cmplx(1.31455839_8, -0.476593971_8, kind=8)
                    ground_truth_b_n(2) = cmplx(-0.0436286479_8, 0.753336608_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.495494395_8, -0.906942308_8, kind=8)
                    ground_truth_b_n(4) = cmplx(1.16279888_8, 0.575285077_8, kind=8)

                    call get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n)

                    open(unit=99, file="temp/test_get_coefficients_a_b_cmplx_barberh.csv", status='unknown')
                    rv = write_csv(99, (/ a_n, b_n /), (/ "n_", "ab" /))
                    close(99)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_cmplx_barberh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_cmplx_barberh

                function test_get_cross_section() result (rv)
                    use lib_field_polarisation
                    use toolbox
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    character(len=*), parameter :: medium_str = "H2O"
                    character(len=*), parameter :: file_name_refractive_index_medium = &
                                                "refractive_index/H2O_Hale_Querry_1973_25_degree_Celsius.csv"


!                    character(len=*), parameter :: particle_str = "Ag"
!                    character(len=*), parameter :: file_name_refractive_index_particle = &
!                                                "refractive_index/Ag_Johnson_Christy_1972_thick_film.csv"

                    character(len=*), parameter :: particle_str = "Au"
                    character(len=*), parameter :: file_name_refractive_index_particle = &
                                                "refractive_index/Au_Johnson_Christy_1972_thick_film.csv"

                    ! auxiliary
                    integer :: i
                    integer :: r

                    type(list_cmplx) :: a_n
                    type(list_cmplx) :: b_n

                    type(list_list_cmplx) :: a_nm
                    type(list_list_cmplx) :: b_nm
                    double precision :: k
                    double precision, dimension(3) :: c
                    double precision, dimension(:), allocatable :: c_sca
                    double precision, dimension(:), allocatable :: c_ext
                    double precision, dimension(:), allocatable :: c_abs

                    double precision :: lambda ! wave length
                    double precision :: lambda_start ! wave length
                    double precision :: lambda_stop ! wave length
                    double precision :: lambda_step ! wave length
                    double precision, dimension(:), allocatable :: lambda_list ! wave length
                    integer :: no

                    double precision :: n_medium

                    double precision :: r_particle
                    double complex :: n_particle

                    type(cartesian_coordinate_real_type) :: d_0_j

                    double precision :: x

                    integer(kind=4), dimension(2) :: n_range

                    type(list_list_cmplx) :: p
                    type(list_list_cmplx) :: q

                    type(cartesian_coordinate_real_type) :: k_cartesian
                    type(spherical_coordinate_real_type) :: k_spherical

                    character(len=150) :: file_name_output
                    integer :: u
                    character(len=25), dimension(4) :: header
                    character(len=25) :: str

                    double precision, dimension(:,:), allocatable :: data_refractive_index_medium
                    double precision, dimension(:), allocatable :: data_refractive_index_medium_interpolation

                    double precision, dimension(:,:), allocatable :: data_refractive_index_particle
                    double precision, dimension(:), allocatable :: data_refractive_index_particle_interpolation

                    type(lib_mie_illumination_parameter) :: illumination

                    lambda_start = 400 * unit_nm
                    lambda_stop = 800 * unit_nm
                    lambda_step = 4 * unit_nm

!                    r_particle = 25 * unit_nm
!                    file_name_output = "temp/c_sca_ag_r_25nm.csv"
                    ! https://refractiveindex.info/?shelf=main&book=Ag&page=Johnson
!                    n_particle = cmplx(0.040000, 7.1155, kind=8)

                    no = int( (lambda_stop - lambda_start) / lambda_step )

                    allocate( c_sca(no) )
                    allocate( c_ext(no) )
                    allocate( c_abs(no) )
                    allocate( lambda_list(no) )

                    do i = 1, no
                        lambda_list(i) = lambda_start + (i-1) * lambda_step
                    end do


                    allocate( data_refractive_index_medium_interpolation(3) )
                    call read_csv(file_name_refractive_index_medium, 3, data_refractive_index_medium)

                    allocate( data_refractive_index_particle_interpolation(3) )
                    call read_csv(file_name_refractive_index_particle, 3, data_refractive_index_particle)

                    allocate(illumination%plane_wave(1))

                    illumination%plane_wave(1)%beam_parameter%polarisation = &
                                            lib_field_polarisation_jones_vector_get_linear_h()

                    print *, "test_get_cross_section"
                    do r=10, 100, 20
!                    do r=100, 400, 100
                        r_particle = r * unit_nm
                        if (r .lt. 100) then
                            write(str, '(A1, I2)') "0", r
                        else
                            write(str, '(I3)') r
                        end if
                        file_name_output = "temp/cross_section/cross_section_" &
                                            // trim(particle_str) // "_r_" // trim(str) // "nm_in_" &
                                            // trim(medium_str) // ".csv"

                        do i = 1, no
                            lambda = lambda_list(i)
                            k = 2D0 * PI * n_medium / lambda

                            k_spherical = make_spherical(k, 0D0, 0D0)
                            k_cartesian = k_spherical

                            illumination%e_field_0 = 1d0
                            illumination%lambda_0 = lambda

                            illumination%plane_wave(1)%g = 1d0
                            illumination%plane_wave(1)%beam_parameter%wave_length_0 = lambda
                            illumination%plane_wave(1)%beam_parameter%e_field_0 = 1d0

                            call data_interpolation(data_refractive_index_particle, 1, lambda / unit_mu, &
                                                    data_refractive_index_particle_interpolation)
                            n_particle = dcmplx(data_refractive_index_particle_interpolation(2), &
                                                data_refractive_index_particle_interpolation(3))

                            call data_interpolation(data_refractive_index_medium, 1, lambda / unit_mu, &
                                                    data_refractive_index_medium_interpolation)
                            n_medium = real(dcmplx(data_refractive_index_medium_interpolation(2), &
                                                data_refractive_index_medium_interpolation(3)))

                            illumination%plane_wave(1)%beam_parameter%refractive_index_medium = n_medium

                            d_0_j%x = 0
                            d_0_j%y = 0
                            d_0_j%z = 0

                            illumination%plane_wave(1)%d_0_i = d_0_j

                            x = abs(k * r_particle)

                            n_range(1) = 1
                            n_range(2) = lib_mie_ss_hf_get_n_c(x)
                            if (n_range(2) .gt. SINGLE_SPHERE_N_MAX) then
                               print *, "WARNING: max degree (", SINGLE_SPHERE_N_MAX, ") reached: ", n_range(2)
                               n_range(2) = SINGLE_SPHERE_N_MAX
#ifdef _PRINT_NOTE_
                            else
                                print *, "NOTE: max degree = ", n_range(2)
#endif
                            end if
!                            n_range(2) = SINGLE_SPHERE_N_MAX

                            call lib_mie_illumination_init_plane_wave((/ 0D0 /), (/ 0D0 /), (/ n_range(2) /))
                            call lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx((/ x /), (/ n_particle / n_medium /), (/ n_range(2) /))

!                            call lib_mie_illumination_get_p_q_j_j(k_cartesian, d_0_j, n_range, p, q)
                            call lib_mie_illumination_get_p_q_j_j(illumination, n_medium, d_0_j, n_range, p, q)

                            allocate (a_n%item(n_range(2)-n_range(1)+1))
                            allocate (b_n%item(n_range(2)-n_range(1)+1))
                            if (aimag(n_particle) .eq. 0) then
                                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x, real(n_particle) / n_medium, n_range, &
                                                                         a_n%item, b_n%item)
                            else
                                call lib_mie_ss_hf_get_coefficients_a_n_b_n(x, n_particle / n_medium, n_range, &
                                                                         a_n%item, b_n%item)
                            end if

                            a_nm = a_n * p
                            b_nm = b_n * q

                            c = get_cross_section(p, q, a_nm, b_nm, k)
                            c_sca(i) = c(1)
                            c_ext(i) = c(2)
                            c_abs(i) = c(3)

                            deallocate (a_n%item)
                            deallocate (b_n%item)
                        end do

                        ! write to csv
                        header(1) = "lambda / m"
                        header(2) = "c_sca"
                        header(3) = "c_ext"
                        header(4) = "c_abs"
                        u = 99
                        open(unit=u, file=file_name_output, status='unknown')
                        rv = write_csv(u, header, lambda_list, &
                                                  c_sca, c_ext, c_abs)
                        close(u)
                    end do

                    rv = .true.
                end function test_get_cross_section
        end function lib_mie_ss_helper_functions_test_functions

        subroutine lib_mie_ss_helper_functions_benchmark()
            implicit none

            integer :: iterations

            iterations = 10**7

            call benchmark_get_coefficients_a_b_real_barberh(iterations, .false.)
            call benchmark_get_coefficients_a_b_real_barberh(iterations, .true.)

            call benchmark_get_coefficients_a_b_cmplx_barberh(iterations, .false.)
            call benchmark_get_coefficients_a_b_cmplx_barberh(iterations, .true.)

            contains

                subroutine benchmark_get_coefficients_a_b_real_barberh(iterations, use_caching)
                    implicit none
                    ! dummy
                    integer :: iterations
                    logical :: use_caching

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    double precision :: x
                    double precision :: m
!                    double precision :: mu
!                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 32/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    ! CPU-time
                    real :: test_start, test_finish
                    ! WALL-time
                    INTEGER :: test_count_start, test_count_finish, test_count_rate

                    m = 1.5d0
                    x= 20.0_8

                    call system_clock(test_count_start, test_count_rate)
                    call cpu_time(test_start)

                    if (use_caching) then
                        call lib_mie_ss_hf_init_coeff_a_n_b_n_real((/ x /), (/ m /), n)
                    end if

                    do i=1, iterations
                        call get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n, use_caching)
                    end do

                    call cpu_time(test_finish)
                    call system_clock(test_count_finish, test_count_rate)

                    print *, ""
                    print *, "------benchmark_get_coefficients_a_b_real_barberh------"
                    print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
                    print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
                    print *, "  Caching: ", use_caching
                    print *, ""
                    print *, "------------------------------------------------------------"
                    print *, ""

                end subroutine

                subroutine benchmark_get_coefficients_a_b_cmplx_barberh(iterations, use_caching)
                    implicit none
                    ! dummy
                    integer :: iterations
                    logical :: use_caching

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    double precision :: x
                    complex(kind=8) :: m
!                    double precision :: mu
!                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/ 1, 32 /)


                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    ! CPU-time
                    real :: test_start, test_finish
                    ! WALL-time
                    INTEGER :: test_count_start, test_count_finish, test_count_rate

                    m = cmplx(1.28, -1.37, kind=8)
                    x= 20.0_8

                    call system_clock(test_count_start, test_count_rate)
                    call cpu_time(test_start)

                    if (use_caching) then
                        call lib_mie_ss_hf_init_coeff_a_n_b_n_cmplx( (/ x /), (/ m /), n)
                    end if

                    do i=1, iterations
                        call get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n, use_caching)
                    end do

                    call cpu_time(test_finish)
                    call system_clock(test_count_finish, test_count_rate)

                    print *, ""
                    print *, "------benchmark_get_coefficients_a_b_cmplx_barberh------"
                    print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
                    print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
                    print *, "  Caching: ", use_caching
                    print *, ""
                    print *, "------------------------------------------------------------"
                    print *, ""

                end subroutine

        end subroutine
end module lib_mie_ss_helper_functions
