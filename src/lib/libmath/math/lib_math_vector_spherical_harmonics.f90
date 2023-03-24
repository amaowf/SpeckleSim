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

#define _DEBUG_


module lib_math_vector_spherical_harmonics
    !$  use omp_lib
    use lib_math_type
    use lib_math_type_operator
    use lib_math_bessel
    use lib_math_legendre
    use lib_math_factorial
    use lib_math_wigner
    use lib_math_constants
    implicit none

    private

    ! --- public ---
    public :: lib_math_vector_spherical_harmonic_constructor
    public :: lib_math_vector_spherical_harmonics_components
    public :: lib_math_vector_spherical_harmonics_translation_coefficient
    public :: lib_math_vector_spherical_harmonics_translate_coefficient

    interface lib_math_vector_spherical_harmonics_components
        module procedure lib_math_vector_spherical_harmonics_components_real_xu
        module procedure lib_math_vector_spherical_harmonics_components_cmplx_xu
    end interface

    interface lib_math_vector_spherical_harmonics_translation_coefficient
        module procedure lib_math_vector_spherical_harmonics_translation_coeff_r
        module procedure lib_math_vector_spherical_harmonics_translation_coeff_spher_r
        module procedure lib_math_vector_spherical_harmonics_translation_coeff_carte_r
    end interface

    interface lib_math_vector_spherical_harmonics_translate_coefficient
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r_spher
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r_carte
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r_multi
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r_spher_m
        module procedure lib_math_vector_spherical_harmonics_translate_coeff_r_carte_m
    end interface

    public :: lib_math_vector_spherical_harmonics_test_functions

    ! --- parameter ---
    integer(kind=1), parameter :: VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND = 4

    contains

        ! Constructor
        !
        ! - initialise calculation of the translation coefficient
        !   - Wigner 3j-Symbol
        !
        ! Argument
        ! ----
        !   j_max: integer
        !
        subroutine lib_math_vector_spherical_harmonic_constructor(j_max)
            implicit none
            ! dummy
            integer(kind=4), intent(in) :: j_max

            ! Wigner 3j-Symbol
            call fwig_table_init(2 * j_max, 3_4)

            ! temp_init must be called per thread
            !$  if (.false.) then
            call fwig_temp_init(2 * j_max)      ! single threaded
            !$  end if

        end subroutine lib_math_vector_spherical_harmonic_constructor

        ! calculation of the components of the vector spherical harmonic
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   r: double precision
        !       distance [m]
        !   k: double precision
        !       wavenumber [1/m]
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Results
        ! ----
        !   M_nm: type(list_spherical_coordinate_cmplx_type)
        !       M component of the vector spherical harmonic
        !   N_nm: type(list_spherical_coordinate_cmplx_type)
        !       N component of the vector spherical harmonic
        !
        ! LaTeX: $$ \mathbf{M}_{m n}^{(J)}=\left[\mathbf{i}_{\theta} i \pi_{m n}(\cos \theta)-\mathbf{i}_{\phi} \tau_{m n}(\cos \theta)\right] z_{n}^{(J)} )(k r) \exp (i m \phi) $$
        !        $$ \mathbf{N}_{m n}^{(J)}=\mathbf{i}_r n\left(n+1 ) P_{n}^{m} (\cos \theta\right) \frac{z_{n}^{(J)}(k r)}{k r} \exp(i m \phi)
        !           + \left[\mathbf{i}_{\theta} \tau_{m n}(\cos \theta)+\mathbf{i}_{\phi} i \pi_{m n}(\cos \theta)\right]
        !           \times \frac{1}{k r} \frac{d}{d r}\left[r z_{n}^{(J)}(k r)\right] \exp (i m \phi) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 2
        subroutine lib_math_vector_spherical_harmonics_components_real_xu(theta, phi, r, k, n_range, z_selector, &
                                                      M_nm, N_nm)
            implicit none
            
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: k
            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_nm

            ! auxiliary
            integer :: n
            integer :: m
            integer :: i

            double precision :: rho

            integer(kind=4) :: number_of_members_n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            type(list_list_real) :: p_nm
            type(list_list_real) :: p_d_nm

            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_n
            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_d_n

            double precision, dimension(:,:), allocatable :: buffer_p_n_m_neg
            double precision, dimension(:,:), allocatable :: buffer_p_d_n_m_neg

            double precision, dimension(n_range(2)-n_range(1)+1) :: z_n_real
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_n_cmplx
            double precision, dimension(n_range(2)-n_range(1)+1) :: z_d_real ! deriviative
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_d_cmplx ! deriviative

            ! Riccati-Bessel
            double precision, dimension(n_range(2)-n_range(1)+1) :: r_real
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_cmplx
            double precision, dimension(n_range(2)-n_range(1)+1) :: r_d_real ! deriviative
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_d_cmplx ! deriviative

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            complex(kind=8), dimension(-n_range(2):n_range(2)) :: exp_i_m_phi
            double precision :: cos_theta

            number_of_members_n = n_range(2) - n_range(1) + 1

            ! --- init ---
            call init_list(p_nm, n_range(1), number_of_members_n)
            call init_list(p_d_nm, n_range(1), number_of_members_n)

            if (allocated(M_nm)) then
                deallocate(M_nm)
            end if

            allocate( M_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (M_nm(i)%coordinate(-i:i))
            end do

            if (allocated(N_nm)) then
                deallocate(N_nm)
            end if

            allocate( N_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (N_nm(i)%coordinate(-i:i))
            end do


            ! --- pre-calculation ---
            rho = k * r

            cos_theta = cos(theta)

            do i=-n_range(2), n_range(2)
                exp_i_m_phi(i)= cmplx(cos(i * phi), sin(i * phi), kind=8)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_real = lib_math_riccati_s_derivative(rho, n_range(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_real = lib_math_riccati_c_derivative(rho, n_range(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function second kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case default
                    z_n_real = 0
                    z_n_cmplx = cmplx(0,0)
                    z_d_real = 0
                    z_d_cmplx = cmplx(0,0)

                    r_real = 0
                    r_cmplx = cmplx(0,0)
                    r_d_real = 0
                    r_d_cmplx = cmplx(0,0)
                    print*, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select


            call lib_math_associated_legendre_polynomial_theta(theta, n_range(2), pi_nm, tau_nm, .false.)


            ! set p_nm(m .eq. 0)
            call lib_math_associated_legendre_polynomial(cos_theta, 0, n_range(1), number_of_members_n , &
                                                         buffer_p_n, buffer_p_d_n, .false.)
            !$OMP PARALLEL DO PRIVATE(n, i)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                p_nm%item(n)%item(0) = buffer_p_n(i)
            end do
            !$OMP END PARALLEL DO

            ! set p_nm(m .ne. 0)
            !       ->n
            !          +
            !         ++
            !    ^   +++  <-  m
            !   m|  ++++
            !        +++  <- -m
            !         ++
            !          +
            !$OMP PARALLEL DO PRIVATE(m, n, buffer_p_n_m_neg, buffer_p_d_n_m_neg)
            do m=1, n_range(2)
                allocate(buffer_p_n_m_neg(2, n_range(2)-m+1))
                allocate(buffer_p_d_n_m_neg(2, n_range(2)-m+1))
                call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, m, m, n_range(2) - m + 1, &
                                                                         buffer_p_n_m_neg, buffer_p_d_n_m_neg, .false.)
                do n=m, n_range(2)
                    p_nm%item(n)%item( m) = buffer_p_n_m_neg(2, n)
                    p_nm%item(n)%item(-m) = buffer_p_n_m_neg(1, n)
                end do

                deallocate(buffer_p_n_m_neg)
                deallocate(buffer_p_d_n_m_neg)
            end do
            !$OMP END PARALLEL DO

            ! --- calculations of the components M and N ---
            ! M_mn
            ! first line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = pi_nm%item(n)%item(m) * z_n_real(i)
                            M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * exp_i_m_phi(m)

                            buffer_real = -tau_nm%item(n)%item(m) * z_n_real(i)
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m)
                            M_nm(n)%coordinate(m)%phi = buffer_cmplx

                            M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
                    if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  M_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                    !$OMP END PARALLEL DO
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = pi_nm%item(n)%item(m)
                            M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)

                            buffer_real = -tau_nm%item(n)%item(m)
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)
                            M_nm(n)%coordinate(m)%phi = buffer_cmplx

                            M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
                    if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  M_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                    !$OMP END PARALLEL DO
            end select

            ! N_mn
            ! second line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m) * z_n_real(i) / rho
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m)
                            N_nm(n)%coordinate(m)%rho = buffer_cmplx

                            buffer_cmplx = cmplx(r_d_real(i) / rho, 0, kind=8) * exp_i_m_phi(m)

                            N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                            N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
                    if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  N_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                    !$OMP END PARALLEL DO
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m) / rho
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m) * z_n_cmplx(i)
                            N_nm(n)%coordinate(m)%rho = buffer_cmplx

                            buffer_cmplx = r_d_cmplx(i) * exp_i_m_phi(m) / rho

                            N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                            N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
                    if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  N_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                    !$OMP END PARALLEL DO
            end select

        end subroutine lib_math_vector_spherical_harmonics_components_real_xu

        ! calculation of the components of the vector spherical harmonic
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   k: complex
        !       wavenumber [1/m]
        !   r: double precision
        !       distance [m]
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Results
        ! ----
        !   rv: complex, dimension(3)
        !       values of the spherical coordinates (rho, theta, phi)
        !
        ! LaTeX: $$ \mathbf{M}_{m n}^{(J)}=\left[\mathbf{i}_{\theta} i \pi_{m n}(\cos \theta)-\mathbf{i}_{\phi} \tau_{m n}(\cos \theta)\right] z_{n}^{(J)} )(k r) \exp (i m \phi) $$
        !        $$ \mathbf{N}_{m n}^{(J)}=\mathbf{i}_r n\left(n+1 ) P_{n}^{m} (\cos \theta\right) \frac{z_{n}^{(J)}(k r)}{k r} \exp(i m \phi)
        !           + \left[\mathbf{i}_{\theta} \tau_{m n}(\cos \theta)+\mathbf{i}_{\phi} i \pi_{m n}(\cos \theta)\right]
        !           \times \frac{1}{k r} \frac{d}{d r}\left[r z_{n}^{(J)}(k r)\right] \exp (i m \phi) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 2
        subroutine lib_math_vector_spherical_harmonics_components_cmplx_xu(theta, phi, k, r, n_range, z_selector, &
                                                      M_nm, N_nm)
            implicit none

            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            complex(kind=8), intent(in) :: k
            double precision, intent(in) :: r
            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_nm

            ! auxiliary
            integer :: n
            integer :: m
            integer :: i

            complex(kind=8) :: rho

            integer(kind=4) :: number_of_members_n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            type(list_list_real) :: p_nm
            type(list_list_real) :: p_d_nm

            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_n
            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_d_n

            double precision, dimension(:,:), allocatable :: buffer_p_n_m_neg
            double precision, dimension(:,:), allocatable :: buffer_p_d_n_m_neg

            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_n_cmplx
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_d_cmplx ! deriviative

            ! Riccati-Bessel
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_cmplx
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_d_cmplx ! deriviative

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            complex(kind=8), dimension(-n_range(2):n_range(2)) :: exp_i_m_phi
            double precision :: cos_theta

            number_of_members_n = n_range(2) - n_range(1) + 1

            ! --- init ---
            call init_list(p_nm, n_range(1), number_of_members_n)
            call init_list(p_d_nm, n_range(1), number_of_members_n)

            allocate( M_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (M_nm(i)%coordinate(-i:i))
            end do

            allocate( N_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (N_nm(i)%coordinate(-i:i))
            end do


            ! --- pre-calculation ---
            rho = k * r

            cos_theta = cos(theta)

            do i=-n_range(2), n_range(2)
                exp_i_m_phi(i)= cmplx(cos(i * phi), sin(i * phi), kind=8)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_cmplx = lib_math_riccati_s_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_cmplx = lib_math_riccati_c_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function first kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case default
                    z_n_cmplx = cmplx(0,0)
                    z_d_cmplx = cmplx(0,0)

                    r_cmplx = cmplx(0,0)
                    r_d_cmplx = cmplx(0,0)
                    print*, "lib_math_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select


            call lib_math_associated_legendre_polynomial_theta(theta, n_range(2), pi_nm, tau_nm, .false.)


            ! set p_nm(m .eq. 0)
            call lib_math_associated_legendre_polynomial(cos_theta, 0, n_range(1), number_of_members_n , &
                                                         buffer_p_n, buffer_p_d_n, .false.)
            !$OMP PARALLEL DO PRIVATE(n, i)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                p_nm%item(n)%item(0) = buffer_p_n(i)
            end do
            !$OMP END PARALLEL DO

            ! set p_nm(m .ne. 0)
            !       ->n
            !          +
            !         ++
            !    ^   +++  <-  m
            !   m|  ++++
            !        +++  <- -m
            !         ++
            !          +
            !$OMP PARALLEL DO PRIVATE(m, n, buffer_p_n_m_neg, buffer_p_d_n_m_neg)
            do m=1, n_range(2)
                allocate(buffer_p_n_m_neg(2, m:n_range(2)-m+1))
                allocate(buffer_p_d_n_m_neg(2, m:n_range(2)-m+1))
                call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, m, m, n_range(2) - m + 1, &
                                                                         buffer_p_n_m_neg, buffer_p_d_n_m_neg, .false.)
                do n=m, n_range(2)
                    p_nm%item(n)%item( m) = buffer_p_n_m_neg(2, n)
                    p_nm%item(n)%item(-m) = buffer_p_n_m_neg(1, n)
                end do
                deallocate(buffer_p_n_m_neg)
                deallocate(buffer_p_d_n_m_neg)
            end do
            !$OMP END PARALLEL DO

            ! --- calculations of the components M and N ---
            ! M_mn
            ! first line eq. (2)
            !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                do m=-n, n
                    buffer_real = pi_nm%item(n)%item(m)
                    M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)

                    buffer_real = -tau_nm%item(n)%item(m)
                    buffer_cmplx = cmplx(buffer_real, 0, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)
                    M_nm(n)%coordinate(m)%phi = buffer_cmplx

                    M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
            if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                print *, "  M_nm(n)%coordinate(m) is NaN"
                print *, "  n = ", n
                print *, "  m = ", m
                print *, "  z_selector: ", z_selector
            end if
#endif
                end do
            end do
            !$OMP END PARALLEL DO

            ! N_mn
            ! second line eq. (2)
            !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                do m=-n, n
                    buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m)
                    buffer_cmplx = cmplx(buffer_real, 0, kind=8) / rho * exp_i_m_phi(m) * z_n_cmplx(i)
                    N_nm(n)%coordinate(m)%rho = buffer_cmplx

                    buffer_cmplx = r_d_cmplx(i) * exp_i_m_phi(m) / rho

                    N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                    N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
            if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                print *, "lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                print *, "  N_nm(n)%coordinate(m) is NaN"
                print *, "  n = ", n
                print *, "  m = ", m
                print *, "  z_selector: ", z_selector
            end if
#endif
                end do
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_math_vector_spherical_harmonics_components_cmplx_xu

!        ! calculation of the components of the vector spherical harmonic
!        !
!        ! Argument
!        ! ----
!        !   theta: double precision
!        !       polar angle
!        !   phi: double precision
!        !       azimuthal angle
!        !   rho: double precision
!        !       dimensionless varibale rho = k*r
!        !       k: wavenumber
!        !       r: distance
!        !   m: integer, dimension(2)
!        !       first element: start index
!        !       second element: last index
!        !       CONDITION: first element .le. second element
!        !   n: integer, dimension(2)
!        !       first element: start index
!        !       second element: last index
!        !       CONDITION: first element .le. second element
!        !   z_selector: integer
!        !       1: spherical Bessel function first kind   j_n
!        !       2: spherical Bessel function second kind  y_n
!        !       3: spherical Hankel function first kind   h^(1)_n
!        !       4: spherical Hankel function second kind  h^(2)_n
!        !
!        ! Results
!        ! ----
!        !   rv: complex, dimension(3)
!        !       values of the spherical coordinates (rho, theta, phi)
!        !
!        ! LaTeX: $$ \begin{aligned} \mathbf{M}_{e m n}=& \frac{-m}{\sin \theta} \sin m \phi P_{n}^{m}(\cos \theta) z_{n}(\rho) \hat{\mathbf{e}}_{\theta} \\ &-\cos m \phi \frac{d P_{n}^{m}(\cos \theta)}{d \theta} z_{n}(\rho) \hat{\mathbf{e}}_{\phi} \end{aligned} $$
!        !
!        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.17, 4.18, 4.19, 4.20
!        subroutine lib_math_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
!                                                      M_emn, M_omn, N_emn, N_omn, &
!                                                      not_calc_Memn, not_calc_Momn, not_calc_Nemn, not_calc_Nomn)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: theta
!            double precision, intent(in) :: phi
!            double precision, intent(in) :: rho
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
!            integer(kind=1) :: z_selector
!
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_omn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_omn
!
!            logical, optional :: not_calc_Memn
!            logical, optional :: not_calc_Momn
!            logical, optional :: not_calc_Nemn
!            logical, optional :: not_calc_Nomn
!
!            ! auxiliary
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n
!
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_dnm
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm_divided_by_sin_theta
!
!            double precision, dimension(n(2)-n(1)+1) :: z_n_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_n_cmplx
!            double precision, dimension(n(2)-n(1)+1) :: z_d_real ! deriviative
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative
!
!            double precision, dimension(n(2)-n(1)+1) :: z_divided_by_rho_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_divided_by_rho_cmplx
!
!            ! Riccati-Bessel
!            double precision, dimension(n(2)-n(1)+1) :: r_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: r_cmplx
!            double precision, dimension(n(2)-n(1)+1) :: r_d_real ! deriviative
!            complex(kind=8), dimension(n(2)-n(1)+1) :: r_d_cmplx ! deriviative
!
!            double precision :: cos_theta
!            double precision :: sin_theta
!            double precision :: minus_sin_theta
!            double precision, dimension(m(2)-m(1)+1) :: m_divided_by_sin_theta
!
!            double precision, dimension(n(2)-n(1)+1) :: cos_m_phi
!            double precision, dimension(n(2)-n(1)+1) :: sin_m_phi
!
!
!            logical :: m_not_calc_Memn
!            logical :: m_not_calc_Momn
!            logical :: m_not_calc_Nemn
!            logical :: m_not_calc_Nomn
!
!            double precision :: buffer_real
!            complex(kind=8) :: buffer_cmplx
!
!
!            number_of_members_m = m(2) - m(1) + 1
!            number_of_members_n = n(2) - n(1) + 1
!
!            ! --- standard value ---
!            m_not_calc_Memn = .false.
!            m_not_calc_Momn = .false.
!            m_not_calc_Nemn = .false.
!            m_not_calc_Nomn = .false.
!
!            if (present(not_calc_Memn)) then
!                m_not_calc_Memn = not_calc_Memn
!            end if
!
!            if (present(not_calc_Momn)) then
!                m_not_calc_Momn = not_calc_Momn
!            end if
!
!            if (present(not_calc_Nemn)) then
!                m_not_calc_Nemn = not_calc_Nemn
!            end if
!            if (present(not_calc_Nomn)) then
!                m_not_calc_Nomn = not_calc_Nomn
!            end if
!
!            ! --- init ---
!            if (.not. m_not_calc_Memn) then
!                allocate(M_emn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (M_emn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Momn) then
!                allocate(M_omn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (M_omn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Nemn) then
!                allocate(N_emn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (N_emn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Nemn) then
!                allocate(N_omn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (N_omn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!
!            ! --- pre-calculation ---
!            cos_theta = cos(theta)
!            sin_theta = sin(theta)
!            minus_sin_theta = -sin(theta)
!
!
!            do i=1, number_of_members_m
!                cos_m_phi(i) = cos(i*phi)
!                sin_m_phi(i) = sin(i*phi)
!                m_divided_by_sin_theta(i) = i / sin_theta
!            end do
!
!            select case (z_selector)
!                case(1)
!                    ! spherical Bessel function first kind j_n
!                    ! internal: calculation with Riccati-Bessel functions: S_n
!                    r_d_real = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, r_real)
!                    z_n_real = r_real / rho
!                case(2)
!                    ! spherical Bessel function second kind y_n
!                    ! internal: calculation with Riccati-Bessel functions: C_n
!                    r_d_real = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, r_real)
!                    z_n_real = r_real / rho
!                case(3)
!                    ! spherical Hankel function first kind   h^(1)_n
!                    ! internal: calculation with Riccati-Bessel functions: Xi_n
!                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n(1), number_of_members_n, r_cmplx)
!                    z_n_cmplx = r_cmplx / rho
!                case(4)
!                    ! spherical Hankel function first kind   h^(2)_n
!                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
!                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n(1), number_of_members_n, r_cmplx)
!                    z_n_cmplx = r_cmplx / rho
!                case default
!                    z_n_real = 0
!                    z_n_cmplx = cmplx(0,0)
!                    z_d_real = 0
!                    z_d_cmplx = cmplx(0,0)
!
!                    r_real = 0
!                    r_cmplx = cmplx(0,0)
!                    r_d_real = 0
!                    r_d_cmplx = cmplx(0,0)
!                    print*, "lib_math_vector_spherical_harmonics_M_emn: ERROR"
!                    print*, "  undefined z_selector value: ", z_selector
!                    return
!            end select
!
!
!            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
!            p_nm = p_nm_divided_by_sin_theta * sin_theta ! pontential of an error
!
!            ! --- calculations of the components M and N ---
!            ! M_emn
!            ! eq. (4.17)
!            if (.not. m_not_calc_Memn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = -cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_emn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)
!
!                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_emn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_emn(i)%coordinate(ii)%phi = buffer_cmplx
!
!                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                end select
!            end if
!
!            ! M_omn
!            ! eq. (4.18)
!            if (.not. m_not_calc_Momn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_omn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)
!
!                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_omn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_omn(i)%coordinate(ii)%phi = buffer_cmplx
!
!                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                end select
!            end if
!
!            ! N_emn
!            ! eq. (4.19)
!            if (.not. m_not_calc_Nemn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = z_divided_by_rho_real(ii) * cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                N_emn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho &
!                                              * r_d_real(ii)
!                                N_emn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
!                                N_emn(i)%coordinate(ii)%rho = buffer_cmplx
!
!                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_emn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_emn(i)%coordinate(ii)%phi = buffer_cmplx
!                            end do
!                        end do
!                end select
!            end if
!
!           ! N_omn
!           ! eq. (4.20)
!            if (.not. m_not_calc_Nomn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = z_divided_by_rho_real(ii) * sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                N_omn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_omn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
!                                N_omn(i)%coordinate(ii)%rho = buffer_cmplx
!
!                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_omn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_omn(i)%coordinate(ii)%phi = buffer_cmplx
!                            end do
!                        end do
!                end select
!            end if
!
!
!        end subroutine lib_math_vector_spherical_harmonics_components_real

        ! calculation of the components of the vector spherical harmonic
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   rho: double precision
        !       dimensionless varibale rho = k*r
        !       k: wavenumber
        !       r: distance
        !   m: intger, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Results
        ! ----
        !   rv: complex, dimension(3)
        !       values of the spherical coordinates (rho, theta, phi)
        !
        ! LaTeX: $$ \begin{aligned} \mathbf{M}_{e m n}=& \frac{-m}{\sin \theta} \sin m \phi P_{n}^{m}(\cos \theta) z_{n}(\rho) \hat{\mathbf{e}}_{\theta} \\ &-\cos m \phi \frac{d P_{n}^{m}(\cos \theta)}{d \theta} z_{n}(\rho) \hat{\mathbf{e}}_{\phi} \end{aligned} $$
        !
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.17
        subroutine lib_math_vector_spherical_harmonics_components_cmplx(theta, phi, rho, m, n, z_selector, &
                                                      M_emn, M_omn, N_emn, N_omn, &
                                                      not_calc_Memn, not_calc_Momn, not_calc_Nemn, not_calc_Nomn)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            complex(kind=8), intent(in) :: rho
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_omn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_omn

            logical, optional :: not_calc_Memn
            logical, optional :: not_calc_Momn
            logical, optional :: not_calc_Nemn
            logical, optional :: not_calc_Nomn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_dnm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm_divided_by_sin_theta

            complex(kind=8), dimension(n(2)-n(1)+1) :: z_n_cmplx
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative

            complex(kind=8), dimension(n(2)-n(1)+1) :: z_divided_by_rho_cmplx

            ! Riccati-Bessel
            complex(kind=8), dimension(n(2)-n(1)+1) :: r_cmplx
            complex(kind=8), dimension(n(2)-n(1)+1) :: r_d_cmplx ! deriviative

            double precision :: cos_theta
            double precision :: sin_theta
            double precision :: minus_sin_theta
            double precision, dimension(m(2)-m(1)+1) :: m_divided_by_sin_theta

            double precision, dimension(n(2)-n(1)+1) :: cos_m_phi
            double precision, dimension(n(2)-n(1)+1) :: sin_m_phi


            logical :: m_not_calc_Memn
            logical :: m_not_calc_Momn
            logical :: m_not_calc_Nemn
            logical :: m_not_calc_Nomn

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx


            number_of_members_m = m(2) - m(1) + 1
            number_of_members_n = n(2) - n(1) + 1

            ! --- standard value ---
            m_not_calc_Memn = .false.
            m_not_calc_Momn = .false.
            m_not_calc_Nemn = .false.
            m_not_calc_Nomn = .false.

            if (present(not_calc_Memn)) then
                m_not_calc_Memn = not_calc_Memn
            end if

            if (present(not_calc_Momn)) then
                m_not_calc_Momn = not_calc_Momn
            end if

            if (present(not_calc_Nemn)) then
                m_not_calc_Nemn = not_calc_Nemn
            end if
            if (present(not_calc_Nomn)) then
                m_not_calc_Nomn = not_calc_Nomn
            end if

            ! --- init ---
            if (.not. m_not_calc_Memn) then
                allocate(M_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Momn) then
                allocate(M_omn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_omn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (N_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_omn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (N_omn(i)%coordinate(number_of_members_n))
                end do
            end if


            ! --- pre-calculation ---
            cos_theta = cos(theta)
            sin_theta = sin(theta)
            minus_sin_theta = -sin(theta)


            do i=1, number_of_members_m
                cos_m_phi(i) = cos(i*phi)
                sin_m_phi(i) = sin(i*phi)
                m_divided_by_sin_theta(i) = i / sin_theta
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_cmplx = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_cmplx = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function first kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case default
                    z_n_cmplx = cmplx(0,0)
                    z_d_cmplx = cmplx(0,0)

                    r_cmplx = cmplx(0,0)
                    r_d_cmplx = cmplx(0,0)
                    print*, "lib_math_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select

!            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
            p_nm = p_nm_divided_by_sin_theta * sin_theta

            ! --- calculations of the components M and N ---
            ! M_emn
            ! eq. (4.17)
            if (.not. m_not_calc_Memn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_emn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_emn(i)%coordinate(ii)%phi = buffer_cmplx

                        M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                    end do
                end do
            end if

            ! M_omn
            ! eq. (4.18)
            if (.not. m_not_calc_Momn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_omn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_omn(i)%coordinate(ii)%phi = buffer_cmplx

                        M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                    end do
                end do
            end if

            ! N_emn
            ! eq. (4.19)
            if (.not. m_not_calc_Nemn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                        N_emn(i)%coordinate(ii)%rho = buffer_cmplx

                        buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_emn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_emn(i)%coordinate(ii)%phi = buffer_cmplx
                    end do
                end do
            end if

           ! N_omn
           ! eq. (4.20)
            if (.not. m_not_calc_Nomn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                        N_omn(i)%coordinate(ii)%rho = buffer_cmplx

                        buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_omn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_omn(i)%coordinate(ii)%phi = buffer_cmplx
                    end do
                end do
            end if


        end subroutine lib_math_vector_spherical_harmonics_components_cmplx

        ! Calculation of the translation transformation coefficients
        ! from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   x: double precision
        !       normalized distance: x = k * d_lj
        !       k: wave number
        !       d_lj: distance from origin l to origin j
        !   theta: double precision
        !       polar coordinate [rad]
        !   phi: double precision
        !       azimuthal coordinate [rad]
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   A_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !   B_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translation_coeff_r(x, theta, phi, &
                                                                           n_range, nu_range, &
                                                                           z_selector, &
                                                                           A_mnkl, B_mnkl)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: nu_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx), intent(inout) :: A_mnkl
            type(list_4_cmplx), intent(inout) :: B_mnkl

            ! auxiliary
            integer(kind=4) :: m
            integer(kind=4) :: n
            integer(kind=4) :: mu
            integer(kind=4) :: nu
            integer(kind=4) :: p
            integer(kind=4) :: q

            real(kind=8) :: a
            real(kind=8) :: b

            integer(kind=4) :: q_max
            integer(kind=4) :: q_max_a
            integer(kind=4) :: q_max_b

            double precision :: factorial

            double precision :: cos_theta

            type(list_list_real) :: p_k_minus_m_p
            type(list_list_real) :: dp_k_minus_m_p
            double precision, dimension(:, :), allocatable :: buffer_pm
            double precision, dimension(:, :), allocatable :: buffer_pd

            double precision, dimension(:), allocatable :: z_n_real
            complex(kind=8), dimension(:), allocatable :: z_n_cmplx

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx
            complex(kind=8) :: buffer_cmplx_a
            complex(kind=8) :: buffer_cmplx_b
            complex(kind=8) :: buffer_cmplx_e_m_mu_phi

            integer, dimension(2) :: n_max_range

            logical :: calc_a
            logical :: calc_b

            !$  integer(kind=4) :: j_max
            !$  logical :: thread_first_run

            ! --- init ---
            n_max_range(1) = min(n_range(1), nu_range(1))
            n_max_range(2) = max(n_range(2), nu_range(2))

            ! eq. 28
            ! p = n + nu - 2 * q
            ! worst case:
            !   q = 0
            !   n = nu = n_max_range(2)
            !   p = 2 * n_max_range(2) - 0
            !   b(p=p+1)
            !   => j_max = 2 * n_max_range(2) + 1

            !$  j_max = 2 * n_max_range(2) + 1
            !$  call fwig_thread_temp_init(2 * j_max)     ! multi threaded

            allocate (buffer_pm(2, 0:2*n_max_range(2)+1+1))
            allocate (buffer_pd(2, 0:2*n_max_range(2)+1+1))

            call init_list(p_k_minus_m_p, 0, 2*n_max_range(2)+1+1)

            call init_list(A_mnkl, n_range(1), n_range(2)-n_range(1)+1, &
                                   nu_range(1), nu_range(2)-nu_range(1)+1)
            call init_list(B_mnkl, n_range(1), n_range(2)-n_range(1)+1, &
                                   nu_range(1), nu_range(2)-nu_range(1)+1)

            ! --- pre-calc ---

            ! z function
            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    allocate(z_n_real(0:2*n_max_range(2)+1))
                    z_n_real = lib_math_bessel_spherical_first_kind(x, 0, 2*n_max_range(2)+1+1)
                case(2)
                    ! spherical Bessel function second kind y_n
                    allocate(z_n_real(0:2*n_max_range(2)+1))
                    z_n_real = lib_math_bessel_spherical_second_kind(x, 0, 2*n_max_range(2)+1+1)
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    allocate(z_n_cmplx(0:2*n_max_range(2)+1))
                    z_n_cmplx = lib_math_hankel_spherical_1(x, 0, 2*n_max_range(2)+1+1)
                case(4)
                    ! spherical Hankel function second kind   h^(2)_n
                    allocate(z_n_cmplx(0:2*n_max_range(2)+1))
                    z_n_cmplx = lib_math_hankel_spherical_2(x, 0, 2*n_max_range(2)+1+1)
                case default
                    z_n_real = 0
                    z_n_cmplx = cmplx(0,0)

                    print*, "lib_math_vector_spherical_harmonics_translation_coeff_r: ERROR"
                    print*, "  undefined z_selector value[1-4]: ", z_selector
                    return
            end select

            ! Legendre Polynomial
            cos_theta = cos(theta)
            call lib_math_associated_legendre_polynomial_range(cos_theta, 0, 2*n_max_range(2)+1, &
                                                               p_k_minus_m_p, dp_k_minus_m_p)
!            call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, 0, &
!                                                                         0, 2*n_max_range(2)+1+1, &
!                                                                         buffer_pm, buffer_pd)
!            do n=0, 2*n_max_range(2)+1
!                p_k_minus_m_p%item(n)%item(0) = buffer_pm(1,n)
!            end do
!            deallocate(buffer_pm)
!            deallocate(buffer_pd)
!
!            do m=1, 2*n_max_range(2)+1
!                call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, m, &
!                                                                             m, 2*n_max_range(2)-m+1+1, &
!                                                                             buffer_pm, buffer_pd)
!                do n=m, 2*n_max_range(2)+1
!                        p_k_minus_m_p%item(n)%item(-m) = buffer_pm(1,n)
!                        p_k_minus_m_p%item(n)%item(m) = buffer_pm(2,n)
!                end do
!                deallocate(buffer_pm)
!                deallocate(buffer_pd)
!            end do

            ! --- calculate A and B
            !$  thread_first_run = .true.
            !$OMP PARALLEL DO PRIVATE(n, m, nu, mu, q_max_a, q_max_b, q_max, q, p, &
            !$OMP  factorial, buffer_cmplx_a, buffer_cmplx_b, buffer_cmplx, &
            !$OMP  buffer_real, buffer_cmplx_e_m_mu_phi, a, b, calc_a, calc_b) &
            !$OMP  FIRSTPRIVATE(thread_first_run)
            do n=n_range(1), n_range(2)
                !$  if (thread_first_run) then
                !$    call fwig_thread_temp_init(2 * j_max)     ! multi threaded
                !$    thread_first_run = .false.
                !$  endif
                do m=-n, n
                    do nu=nu_range(1), nu_range(2)
                        do mu=-nu, nu
                            q_max_a = get_q_max_a(-m, n, mu, nu)
                            q_max_b = get_q_max_b(-m, n, mu, nu)
                            q_max = max(q_max_a, q_max_b)

                            factorial = (2 * nu + 1)
                            factorial = factorial * lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(n, m)
                            factorial = factorial * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(nu, mu)
                            factorial = factorial / real( 2 * n * ( n + 1 ) , kind=8)

                            buffer_cmplx_a = cmplx(0,0)
                            buffer_cmplx_b = cmplx(0,0)
                            do q=0, q_max
                                p = n + nu - 2 * q

                                if (q .ge. 0 &
                                    .and. q .le. q_max_a) then
                                    calc_a = .true.
                                else
                                    calc_a = .false.
                                end if

                                if (q .ge. 1 &
                                    .and. q .le. q_max_b) then
                                    calc_b = .true.
                                else
                                    calc_b = .false.
                                end if

                                call ab_xu_cruzan_eq34(-m, n, mu, nu, p, a, b, &
                                                       calc_a = calc_a, calc_b = calc_b)

                                if (calc_a) then

                                    buffer_real = (n*(n+1) + nu*(nu+1) - p*(p+1)) &
                                                  * a * p_k_minus_m_p%item(p)%item(mu-m)
                                    select case (z_selector)
                                        case(1,2)
                                            buffer_real = buffer_real * z_n_real(p)
                                            buffer_cmplx_a = buffer_cmplx_a &
                                                             + cmplx(0,1)**p * buffer_real
                                        case(3,4)
                                            buffer_cmplx = buffer_real * z_n_cmplx(p)
                                            buffer_cmplx_a = buffer_cmplx_a &
                                                             + cmplx(0,1)**p * buffer_cmplx
                                    end select
                                end if


                                if (calc_b) then

                                    buffer_real = sqrt( ( real(p+1, kind=8)**2 - real(n-nu, kind=8)**2 ) &
                                                            * ( real(n+nu+1, kind=8)**2 - real(p+1, kind=8)**2 ) ) &
                                                      * b *  p_k_minus_m_p%item(p+1)%item(mu-m)

                                    select case (z_selector)
                                        case(1,2)
                                            buffer_real = buffer_real * z_n_real(p+1)
                                            buffer_cmplx_b = buffer_cmplx_b &
                                                             + cmplx(0,1)**(p+1) * buffer_real
                                        case(3,4)
                                            buffer_cmplx = buffer_real * z_n_cmplx(p+1)
                                            buffer_cmplx_b = buffer_cmplx_b &
                                                             + cmplx(0,1)**(p+1) * buffer_cmplx
                                    end select
                                end if
                            end do
                            buffer_real = (mu-m) * phi
                            buffer_cmplx_e_m_mu_phi = cmplx(cos(buffer_real), sin(buffer_real), kind=8)

                            buffer_cmplx = factorial * buffer_cmplx_e_m_mu_phi

                            A_mnkl%item(n)%item(m)%item(nu)%item(mu) = (-1)**(-m) * buffer_cmplx * buffer_cmplx_a

                            B_mnkl%item(n)%item(m)%item(nu)%item(mu) = (-1)**(-m+1) * buffer_cmplx * buffer_cmplx_b
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_math_vector_spherical_harmonics_translation_coeff_r

        ! Calculation of the translation transformation coefficients
        ! from the l-th coordinate system to the j-th coordinate system
        !
        ! Argument
        ! ----
        !   x: type(spherical_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   A_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !   B_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !
        ! Reference: Experimental and theoretical results of light scattering by aggregates of spheres, Yu-lin Xu and Bo Å. S. Gustafson
        subroutine lib_math_vector_spherical_harmonics_translation_coeff_spher_r(x, &
                                                                                n_range, nu_range, z_selector,&
                                                                                A_mnkl, B_mnkl)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: nu_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx), intent(inout) :: A_mnkl
            type(list_4_cmplx), intent(inout) :: B_mnkl

            call lib_math_vector_spherical_harmonics_translation_coeff_r(x%rho, x%theta, x%phi, &
                                                                                 n_range, nu_range, z_selector,&
                                                                                 A_mnkl, B_mnkl)

        end subroutine lib_math_vector_spherical_harmonics_translation_coeff_spher_r

        ! Calculation of the translation transformation coefficients
        ! from the l-th coordinate system to the j-th coordinate system
        !
        ! Argument
        ! ----
        !   x: type(cartesian_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   A_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !   B_mnkl: type(list_4_cmplx)
        !       4 dimensional matrix
        !
        ! Reference: Experimental and theoretical results of light scattering by aggregates of spheres, Yu-lin Xu and Bo Å. S. Gustafson
        subroutine lib_math_vector_spherical_harmonics_translation_coeff_carte_r(x, &
                                                                                n_range, nu_range, z_selector,&
                                                                                A_mnkl, B_mnkl)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: nu_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx), intent(inout) :: A_mnkl
            type(list_4_cmplx), intent(inout) :: B_mnkl

            ! auxiliary
            type(spherical_coordinate_real_type) :: m_x

            m_x = x

            call lib_math_vector_spherical_harmonics_translation_coeff_r(m_x%rho, m_x%theta, m_x%phi, &
                                                                                 n_range, nu_range, z_selector,&
                                                                                 A_mnkl, B_mnkl)

        end subroutine lib_math_vector_spherical_harmonics_translation_coeff_carte_r


        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: double precision
        !       normalized distance: x = k * d_lj
        !       k: wave number
        !       d_lj: distance from origin l to origin j
        !   theta: double precision
        !       polar coordinate [rad]
        !   phi: double precision
        !       azimuthal coordinate [rad]
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r(a_nm, b_nm, &
                                                                         x, theta, phi, &
                                                                         n_range, nu_range, &
                                                                         z_selector, &
                                                                         a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            double precision, intent(in) :: x
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), intent(inout) :: a_nmt
            type(list_list_cmplx), intent(inout) :: b_nmt

            ! auxiliary
            integer :: n
            integer :: m

            type(list_4_cmplx) :: A_mnkl
            type(list_4_cmplx) :: B_mnkl

            call lib_math_vector_spherical_harmonics_translation_coefficient(x, theta, phi, n_range, nu_range, z_selector,&
                                                                             A_mnkl, B_mnkl)

            call init_list(a_nmt, nu_range(1), nu_range(2) - nu_range(1) + 1)
            call init_list(b_nmt, nu_range(1), nu_range(2) - nu_range(1) + 1)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    a_nmt%item(n)%item(m) = sum(a_nm * A_mnkl%item(n)%item(m) + b_nm * B_mnkl%item(n)%item(m))
                    b_nmt%item(n)%item(m) = sum(a_nm * B_mnkl%item(n)%item(m) + b_nm * A_mnkl%item(n)%item(m))
                end do
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r

        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: type(spherical_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_spher(a_nm, b_nm, &
                                                                               x, &
                                                                               n_range, nu_range, &
                                                                               z_selector, &
                                                                               a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            type(spherical_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), intent(inout) :: a_nmt
            type(list_list_cmplx), intent(inout) :: b_nmt

            ! auxiliary
            call lib_math_vector_spherical_harmonics_translate_coeff_r(a_nm, b_nm, &
                                                                       x%rho, x%theta, x%phi, &
                                                                       n_range, nu_range, &
                                                                       z_selector, &
                                                                       a_nmt, b_nmt)

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_spher

        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: type(cartesian_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_carte(a_nm, b_nm, &
                                                                               x, &
                                                                               n_range, nu_range, &
                                                                               z_selector, &
                                                                               a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            type(cartesian_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), intent(inout) :: a_nmt
            type(list_list_cmplx), intent(inout) :: b_nmt

            ! auxiliary
            type(spherical_coordinate_real_type) :: m_x

            m_x = x

            call lib_math_vector_spherical_harmonics_translate_coeff_r(a_nm, b_nm, &
                                                                       m_x%rho, m_x%theta, m_x%phi, &
                                                                       n_range, nu_range, &
                                                                       z_selector, &
                                                                       a_nmt, b_nmt)

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_carte

        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: double precision
        !       normalized distance: x = k * d_lj
        !       k: wave number
        !       d_lj: distance from origin l to origin j
        !   theta: double precision
        !       polar coordinate [rad]
        !   phi: double precision
        !       azimuthal coordinate [rad]
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_multi(a_nm, b_nm, &
                                                                               x, theta, phi, &
                                                                               n_range, nu_range, &
                                                                               z_selector, &
                                                                               a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), dimension(:), intent(in) :: a_nm
            type(list_list_cmplx), dimension(lbound(a_nm, 1):ubound(a_nm, 1)), intent(in) :: b_nm
            double precision, intent(in) :: x
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: a_nmt
            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: b_nmt

            ! auxiliary
            integer :: i
            integer :: n
            integer :: m

            type(list_4_cmplx) :: A_mnkl
            type(list_4_cmplx) :: B_mnkl

            if (allocated(a_nmt)) deallocate(a_nmt)
            allocate(a_nmt(lbound(a_nm, 1):ubound(a_nm, 1)))

            if (allocated(b_nmt)) deallocate(b_nmt)
            allocate(b_nmt(lbound(a_nm, 1):ubound(a_nm, 1)))

            call lib_math_vector_spherical_harmonics_translation_coefficient(x, theta, phi, n_range, nu_range, z_selector,&
                                                                             A_mnkl, B_mnkl)

            !$OMP PARALLEL DO PRIVATE(i)
            do i = lbound(a_nm, 1), ubound(a_nm, 1)
                call init_list(a_nmt(i), nu_range(1), nu_range(2) - nu_range(1) + 1)
                call init_list(b_nmt(i), nu_range(1), nu_range(2) - nu_range(1) + 1)

                !$OMP PARALLEL DO PRIVATE(n, m)
                do n = n_range(1), n_range(2)
                    do m = -n, n
                        a_nmt(i)%item(n)%item(m) = sum(a_nm(i) * A_mnkl%item(n)%item(m) + b_nm(i) * B_mnkl%item(n)%item(m))
                        b_nmt(i)%item(n)%item(m) = sum(a_nm(i) * B_mnkl%item(n)%item(m) + b_nm(i) * A_mnkl%item(n)%item(m))
                    end do
                end do
                !$OMP END PARALLEL DO
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_multi

        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: type(spherical_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_spher_m(a_nm, b_nm, &
                                                                               x, &
                                                                               n_range, nu_range, &
                                                                               z_selector, &
                                                                               a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), dimension(:), intent(in) :: a_nm
            type(list_list_cmplx), dimension(lbound(a_nm, 1):ubound(a_nm, 1)), intent(in) :: b_nm
            type(spherical_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: a_nmt
            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: b_nmt

            ! auxiliary
            call lib_math_vector_spherical_harmonics_translate_coeff_r_multi(a_nm, b_nm, &
                                                                             x%rho, x%theta, x%phi, &
                                                                             n_range, nu_range, &
                                                                             z_selector, &
                                                                             a_nmt, b_nmt)

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_spher_m

        ! Calculation of the coefficients from the l-th coordinate system to the j-th coordinate system
        !
        !
        !
        !            z ^
        !              |
        !          K_l |----> x
        !             /
        !            /
        !     theta / d_lj
        !      z ^ /
        !        |/
        !    K_j |-----> x
        !
        !   theta = angle(z_j, d_lj)
        !
        !
        !              z
        !          K_l o----> x
        !             /
        !            /
        !           / d_lj
        !          /
        !        z/ phi
        !    K_j o-----> x
        !
        !   phi = angle(d_lj, x_j)
        !
        !
        ! Argument
        ! ----
        !   a_nm, b_nm: type(list_list_cmplx), dimension(:)
        !       coefficients of system l
        !   x: type(cartesian_coordinate_real_type)
        !       coordinate of the l-th coordinate system respect to the j-th coordinate system
        !   n_range: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Returns
        ! ----
        !   a_nmt, b_nmt: type(list_list_cmplx)
        !       translated coefficients at system j
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_carte_m(a_nm, b_nm, &
                                                                               x, &
                                                                               n_range, nu_range, &
                                                                               z_selector, &
                                                                               a_nmt, b_nmt)
            implicit none
            ! dummy
            type(list_list_cmplx), dimension(:), intent(in) :: a_nm
            type(list_list_cmplx), dimension(lbound(a_nm, 1):ubound(a_nm, 1)), intent(in) :: b_nm
            type(cartesian_coordinate_real_type), intent(in) :: x
            integer(kind=4), dimension(2), intent(in) :: n_range
            integer(kind=4), dimension(2), intent(in) :: nu_range
            integer(kind=1), intent(in) :: z_selector

            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: a_nmt
            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: b_nmt

            ! auxiliary
            type(spherical_coordinate_real_type) :: m_x

            m_x = x

            call lib_math_vector_spherical_harmonics_translate_coeff_r_multi(a_nm, b_nm, &
                                                                             m_x%rho, m_x%theta, m_x%phi, &
                                                                             n_range, nu_range, &
                                                                             z_selector, &
                                                                             a_nmt, b_nmt)

        end subroutine lib_math_vector_spherical_harmonics_translate_coeff_r_carte_m

        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        !            eq. 25
        function get_q_max_a(m, n, k, l) result (rv)
            implicit none
            ! dummy
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            integer(kind=4), intent(in) :: k
            integer(kind=4), intent(in) :: l

            integer(kind=4) :: rv

            rv = min(n, l, int(floor( real(n + l - abs(m + k))/2.0 )))
        end function

        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        !            eq. 60
        function get_q_max_b(m, n, k, l) result (rv)
            implicit none
            ! dummy
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            integer(kind=4), intent(in) :: k
            integer(kind=4), intent(in) :: l

            integer(kind=4) :: rv

            rv = min(n, l, int(floor( real(n + l + 1 - abs(m + k))/2.0 )))
        end function

        ! Reference: [1] Experimental and theoretical results of light scattering by aggregates of spheres, Yu-lin Xu and Bo Å. S. Gustafson
        !                eq. 5
        !            [2] Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        !                eq. 28
        function get_p(n, l, q) result (rv)
            implicit none
            ! dummy
            integer(kind=4), intent(in) :: n
            integer(kind=4), intent(in) :: l
            integer(kind=4), intent(in) :: q

            integer(kind=4) :: rv

            rv = n + l - 2 * q
        end function

        ! todo: optimize for a(p=p) and b(p=p+1)
        !
        ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories, Yu-lin Xu
        !            eq. 5, eq. 61
        subroutine ab_xu_cruzan_eq34(m, n, mu, nu, p, a, b, calc_a, calc_b)
            implicit none
            ! dummy
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            integer(kind=4), intent(in) :: mu
            integer(kind=4), intent(in) :: nu
            integer(kind=4), intent(in) :: p

            real(kind=8), intent(inout) :: a
            real(kind=8), intent(inout) :: b

            logical, intent(in), optional :: calc_a
            logical, intent(in), optional :: calc_b

            ! auxiliary
            real(kind=8), dimension(3) :: buffer_factorial
            real(kind=8), dimension(2) :: buffer_wigner
            real(kind=8) :: minus_1_m_mu

            logical :: m_calc_a
            logical :: m_calc_b

            m_calc_a = .true.
            m_calc_b = .true.
            if ( present(calc_a) ) then
                m_calc_a = calc_a
            end if

            if ( present(calc_b) ) then
                m_calc_b = calc_b
            end if

            minus_1_m_mu = (-1d0)**(m + mu)

            ! --- a: eq. 5 ---
            buffer_factorial(1) = lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(n, m)
            buffer_factorial(2) = lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(nu, mu)

            buffer_wigner(2) = lib_math_wigner_3j(n, nu, p, 0, 0, 0)

            if ( m_calc_a ) then
                buffer_factorial(3) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(p, m + mu)

                buffer_wigner(1) = lib_math_wigner_3j(n, nu, p, m, mu, -m-mu)

                a = (2_4 * p + 1_4) * sqrt(buffer_factorial(1) * buffer_factorial(2) * buffer_factorial(3)) &
                     * buffer_wigner(1) * buffer_wigner(2)
                a = minus_1_m_mu * a
            else
                a = 0
            end if

            ! --- b: eq. 61 ---
            if (  m_calc_b ) then
                buffer_factorial(3) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(p + 1_4, m + mu)

                buffer_wigner(1) = lib_math_wigner_3j(n, nu, p+1_4, m, mu, -m-mu)

                b = (2_4 * p + 3_4) * sqrt(buffer_factorial(1) * buffer_factorial(2) * buffer_factorial(3)) &
                     * buffer_wigner(1) * buffer_wigner(2)
                b = minus_1_m_mu * b
            else
                b = 0
            end if

        end subroutine ab_xu_cruzan_eq34

        function lib_math_vector_spherical_harmonics_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)
            ! CPU-time
            real :: test_start, test_finish
            ! WALL-time
            INTEGER :: test_count_start, test_count_finish, test_count_rate

            rv = 0

            call system_clock(test_count_start, test_count_rate)
            call cpu_time(test_start)

            if (.not. test_lib_math_vector_spherical_harmonics_components_real_xu()) then
                rv = rv + 1
            end if
!            if (.not. test_lib_math_vector_spherical_harmonics_components_cmplx_xu()) then
!                rv = rv + 1
!            end if
            if (.not. test_ab_xu_cruzan_eq34()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_vector_spherical_harmonics_translation_coeff_r()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_vector_spherical_harmonics_translation_coeff_r_v2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_vector_spherical_harmonics_components_real_xu_2()) then
                rv = rv + 1
            end if


            call cpu_time(test_finish)
            call system_clock(test_count_finish, test_count_rate)

            print *, "----lib_math_vector_spherical_harmonics_test_functions----"
            print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
            print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
            print *, ""
            if (rv == 0) then
                print *, "lib_math_vector_spherical_harmonics_test_functions tests: OK"
            else
                print *, rv,"lib_math_vector_spherical_harmonics_test_functions test(s) FAILED"
            end if
            print *, "---------------------------------------------------------"
            print *, ""

            contains

!            function test_lib_math_vector_spherical_harmonics_components_real() result (rv)
!                implicit none
!                ! dummy
!                logical :: rv
!
!                ! auxiliary
!                double precision :: theta
!                double precision :: phi
!                double precision :: rho
!                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
!                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
!                integer(kind=1) :: z_selector
!
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_emn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_omn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_emn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_omn
!
!!                logical :: not_calc_Memn
!!                logical :: not_calc_Momn
!!                logical :: not_calc_Nemn
!!                logical :: not_calc_Nomn
!
!                theta = 0
!                phi = 0
!                rho = 100
!
!                z_selector = 3
!
!                m = (/ 1, 1 /)
!                n = (/ 1, 3 /)
!
!                call lib_math_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
!                                                                        M_emn, M_omn, N_emn, N_omn)
!
!
!            end function

            function test_lib_math_vector_spherical_harmonics_components_real_xu() result (rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                character(len=*), parameter :: file_name_M_mn = &
                     "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_M_mn.csv"
                 character(len=*), parameter :: file_name_N_mn = &
                     "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_N_mn.csv"

                ! auxiliary
                integer :: i
                integer :: ii
                integer :: n_value
                integer :: m_value

                double precision :: buffer
                complex(kind=8) :: buffer_cmplx

                double precision :: theta
                double precision :: phi
                double precision :: k
                double precision :: r
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_mn

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_N_mn

                double precision, dimension(:,:), allocatable :: csv_data
                integer :: csv_columns

                theta = 0.2_8
                phi = 0_8
                ! n = 1
                ! lam = 10**-6 m
                ! k = n * 2 Pi / lam
                k = 2.0_8 * PI !* 10.0_8**(6.0_8)
                ! r = 50 * 10**-6 m
                r = 50.0_8 !* 10.0_8**(-6.0_8)

                z_selector = 3

                n = (/ 1, 4 /)

                rv = .false.

                call init_list(ground_truth_M_mn, n(1), n(2)-n(1)+1)
                call init_list(ground_truth_N_mn, n(1), n(2)-n(1)+1)

                ! load ground truth M_mn
                if (file_exists(file_name_M_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_M_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_M_mn

                    return
                end if

                ! load ground truth M_mn
                if (file_exists(file_name_N_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_N_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_components_real_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_N_mn

                    return
                end if

                ! calculate M_mn, N_mn
                call lib_math_vector_spherical_harmonics_components_real_xu(theta, phi, r, k, n, z_selector, &
                                                                           M_mn, N_mn)

                ! evaluate
                rv = .true.
                print *, "test_lib_math_vector_spherical_harmonics_components_real_xu:"
                print *, "  M_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(M_mn(i)%coordinate(ii) - ground_truth_M_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "    n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "    n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

                print *, "  N_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(N_mn(i)%coordinate(ii) - ground_truth_N_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "    n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "    n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

            end function test_lib_math_vector_spherical_harmonics_components_real_xu

            function test_lib_math_vector_spherical_harmonics_components_cmplx_xu() result (rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                character(len=*), parameter :: file_name_M_mn = &
                     "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_M_mn_cmplx.csv"
                 character(len=*), parameter :: file_name_N_mn = &
                     "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_N_mn_cmplx.csv"

                ! auxiliary
                integer :: i
                integer :: ii
                integer :: n_value
                integer :: m_value

                double precision :: buffer
                complex(kind=8) :: buffer_cmplx

                double precision :: theta
                double precision :: phi
                complex(kind=8) :: k
                double precision :: r
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_mn

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_N_mn

                double precision, dimension(:,:), allocatable :: csv_data
                integer :: csv_columns

                theta = 0.2
                phi = 0
                ! Ag (Silver), lambda = 1 mu
                ! https://refractiveindex.info/?shelf=main&book=Ag&page=Johnson
                k = cmplx(0.04, 7.1155, kind=8) * 2.0 * PI / 10.0_8**(-6.0_8)
                r = 2.0_8 * 10.0_8**(-6.0_8)

                z_selector = 3

                n = (/ 1, 5 /)

                rv = .false.

                call init_list(ground_truth_M_mn, n(1), n(2)-n(1)+1)
                call init_list(ground_truth_N_mn, n(1), n(2)-n(1)+1)

                ! load ground truth M_mn
                if (file_exists(file_name_M_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_M_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_components_cmplx_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_M_mn

                    return
                end if

                ! load ground truth M_mn
                if (file_exists(file_name_N_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_N_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_components_cmplx_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_N_mn

                    return
                end if

                ! calculate M_mn, N_mn
                call lib_math_vector_spherical_harmonics_components_cmplx_xu(theta, phi, k, r, n, z_selector, &
                                                                           M_mn, N_mn)

                ! evaluate
                rv = .true.
                print *, "test_lib_math_vector_spherical_harmonics_components_cmplx_xu:"
                print *, "  M_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(M_mn(i)%coordinate(ii) - ground_truth_M_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "    n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "    n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

                print *, "  N_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(N_mn(i)%coordinate(ii) - ground_truth_N_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "    n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "    n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

            end function test_lib_math_vector_spherical_harmonics_components_cmplx_xu

            function test_lib_math_vector_spherical_harmonics_components_real_xu_2() result (rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                character(len=*), parameter :: file_name = "temp/spherical_harmonics.csv"
                integer, parameter :: no = 1000

                ! auxiliary
                integer :: i

                double precision :: theta
                double precision :: phi
                double precision :: k
                double precision :: r
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

                integer :: u
                character(len=25), dimension(4) :: header
                double precision, dimension(no) :: degree_list
                double complex, dimension(no) :: M_values
                double complex, dimension(no) :: N_values
                double complex, dimension(no) :: MN_sum

                theta = 0.2_8
                phi = 0_8
                ! n = 1
                ! lam = 10**-6 m
                ! k = n * 2 Pi / lam
                k = 2.0_8 * PI !* 10.0_8**(6.0_8)
                ! r = 50 * 10**-6 m
                r = 50.0_8 !* 10.0_8**(-6.0_8)

                z_selector = 3

                n = (/ 1, 2 /)

                rv = .false.

                do i=1, no
                    theta = (i-1) * PI / no
                    degree_list(i) = theta

                    ! calculate M_mn, N_mn
                    call lib_math_vector_spherical_harmonics_components_real_xu(theta, phi, r, k, n, z_selector, &
                                                                           M_nm, N_nm)

                    M_values(i) = M_nm(n(2))%coordinate(0)%rho
                    N_values(i) = N_nm(n(2))%coordinate(0)%rho
                    MN_sum(i) = M_values(i) + N_values(i)
                end do

                ! write to csv
                header(1) = "theta / rad"
                header(2) = "M"
                header(3) = "N"
                header(4) = "sum"
                u = 99
                open(unit=u, file=file_name, status='unknown')
                rv = write_csv(u, header, degree_list, M_values, N_values, MN_sum)
                close(u)

                ! evaluate
                rv = .true.
                print *, "test_lib_math_vector_spherical_harmonics_components_real_xu_2:"

            end function test_lib_math_vector_spherical_harmonics_components_real_xu_2

            function test_ab_xu_cruzan_eq34() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, dimension(2), parameter :: q = (/ 1, 5 /)

                ! auxiliary
                integer(kind=4) :: i

                integer(kind=4) :: m
                integer(kind=4) :: n
                integer(kind=4) :: mu
                integer(kind=4) :: nu
                integer(kind=4) :: p

                real(kind=8), dimension(q(1):q(2)) :: a
                real(kind=8), dimension(q(1):q(2)) :: b

                real(kind=8), dimension(q(1):q(2)) :: ground_truth_a
                real(kind=8), dimension(q(1):q(2)) :: ground_truth_b

                real(kind=8) :: buffer

                m = 15!-15
                n = 16
                mu = -15
                nu = 20

                ground_truth_a(1) = -2.289221071014754D-8
                ground_truth_b(1) = -4.708684376179547D-9

                ground_truth_a(2) = 2.691860941883802D-7
                ground_truth_b(2) = 8.05644974652361D-8

                ground_truth_a(3) = -1.989821339021589D-6
                ground_truth_b(3) = -7.524763818766572D-7

                ground_truth_a(4) = 1.033853835595824D-5
                ground_truth_b(4) = 4.673349949711971D-6

                ground_truth_a(5) = -3.995549233873365D-5
                ground_truth_b(5) = -2.099627385849716D-5

                do i=q(1), q(2)
                    ! eq. 5
                    p = n + nu - 2 * i
                    call ab_xu_cruzan_eq34(m, n, mu, nu, p, a(i), b(i))
                end do

                rv = .true.
                print *, "test_ab_xu_cruzan_eq34:"
                print *, "  a:"
                do i=q(1), q(2)
                    buffer = ground_truth_a(i) - a(i)
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    q: ", i, "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "    q: ", i , ": OK"
                    end if
                end do

                print *, "  b:"
                do i=q(1), q(2)
                    buffer = ground_truth_b(i) - b(i)
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    q: ", i, "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "    q: ", i , ": OK"
                    end if
                end do

            end function test_ab_xu_cruzan_eq34

            function test_lib_math_vector_spherical_harmonics_translation_coeff_r() result(rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 4
                character(len=*), parameter :: file_name = &
                        "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_AB_mnmunu.csv"
                character(len=*), parameter :: format = '(10(X)A, ES19.9)'
                ! auxiliary
                integer :: i
                double precision :: buffer
                double complex :: buffer_cmplx
                double precision :: x
                double precision :: theta
                double precision :: phi
                integer(kind=4), dimension(2) :: n_range
                integer(kind=4), dimension(2) :: nu_range
                integer(kind=1) :: z_selector

                type(list_4_cmplx) :: A_mnkl
                type(list_4_cmplx) :: B_mnkl


                integer :: m
                integer :: n
                integer :: mu
                integer :: nu
                type(list_4_cmplx) :: ground_truth_A
                type(list_4_cmplx) :: ground_truth_B

                double precision :: erg
                double precision :: ground_truth_erg
                integer :: erg_exponent
                double precision :: erg_mantissa

                double precision, dimension(:,:), allocatable :: csv_data
                double precision, dimension(:,:), allocatable :: csv_data_calc
                integer :: csv_columns

                character(len=50), dimension(8) :: header

                z_selector = 3
                x = 2.0
                theta = 0.5
                phi = 0.5

                n_range  = (/ 1, 5 /)
                nu_range = (/ 1, 5 /)


                call init_list(ground_truth_A, n_range(1), n_range(2)-n_range(1)+1, &
                                               nu_range(1), nu_range(2)-nu_range(1)+1)
                call init_list(ground_truth_B, n_range(1), n_range(2)-n_range(1)+1, &
                                               nu_range(1), nu_range(2)-nu_range(1)+1)
                n_range  = (/ 1, 5 /)
                nu_range = (/ 1, 5 /)

                call lib_math_factorial_initialise_caching(90)

                ! load ground truth
                if (file_exists(file_name)) then
                    csv_columns = 8
                    call read_csv(file_name, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        m = int(csv_data(i, 1))
                        n = int(csv_data(i, 2))
                        mu = int(csv_data(i, 3))
                        nu = int(csv_data(i, 4))

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_A%item(n)%item(m)%item(nu)%item(mu) = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_B%item(n)%item(m)%item(nu)%item(mu) = buffer_cmplx
                    end do
!                m =  (/ -2,  8,  0, -2 /)
!                n =  (/ 11, 10, 10,  6 /)
!                mu = (/ 3 , -9,  0, -2 /)
!                nu = (/ 9 , 12, 10, 10 /)
!                m = (/ -2, 8 /)
!                n = (/ 11, 10 /)
!                mu = (/ 3, -9 /)
!                nu = (/ 9, 12 /)
!                ground_truth_A(1) = cmplx(.7726121583d+12, .103425m5820d+13, kind=8)
!                ground_truth_B(1) = cmplx(.1222239141d+11, -0.9130398908d+10, kind=8)
!
!                ground_truth_A(2) = cmplx(.3663964990d+35, -.2762412192d+35, kind=8)
!                ground_truth_B(2) = cmplx(-.8370892023d+32, -.1110285257d+33, kind=8)
!
!                ground_truth_A(3) = cmplx(.2969682019d+00, -.1928601440d+18, kind=8)
!                ground_truth_B(3) = cmplx(0.0, 0.0, kind=8)
!
!                ground_truth_A(4) = cmplx(.1377011649d-01, .2385575934d+13, kind=8)
!                ground_truth_B(4) = cmplx(-.3282035237d+12, .1587043209d-02, kind=8)
!
!                n_range = (/ 1, maxval(n) /)
!                nu_range = (/ 1, maxval(nu) /)
                call lib_math_vector_spherical_harmonics_translation_coeff_r(x, theta, phi, &
                                                                             n_range, nu_range, z_selector, &
                                                                             A_mnkl, B_mnkl)

                allocate(csv_data_calc, mold=csv_data)
                i = lbound(csv_data_calc, 1) - 1
                do n = n_range(1), n_range(2)
                do m = -n, n
                do nu = nu_range(1), nu_range(2)
                do mu = -nu, nu
                    i = i + 1
                    csv_data_calc(i, 1) = m
                    csv_data_calc(i, 2) = n
                    csv_data_calc(i, 3) = mu
                    csv_data_calc(i, 4) = nu

                    buffer_cmplx = A_mnkl%item(n)%item(m)%item(nu)%item(mu)
                    csv_data_calc(i, 5) = real(buffer_cmplx)
                    csv_data_calc(i, 6) = aimag(buffer_cmplx)

                    buffer_cmplx = B_mnkl%item(n)%item(m)%item(nu)%item(mu)
                    csv_data_calc(i, 7) = real(buffer_cmplx)
                    csv_data_calc(i, 8) = aimag(buffer_cmplx)
                end do
                end do
                end do
                end do

!                m,n,mu,nu,Re[A],Im[A],Re[B],Im[B]
                header(1) = "m"
                header(2) = "n"
                header(3) = "mu"
                header(4) = "nu"
                header(5) = "Re[A]"
                header(6) = "Im[A]"
                header(7) = "Re[B]"
                header(8) = "Im[B]"

                open(unit=99, file= "temp/Translation_Coefficient.csv", status='unknown')
                rv = write_csv(99, csv_data_calc, header)
                close(99)

                rv = .true.
                print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r:"
                print *, "  A:"
                do n=n_range(1), n_range(2)
                do m=-n, n
                do nu=nu_range(1), nu_range(2)
                do mu=-nu, nu
                    print *, "  n=", n, "m=", m, "nu=", nu, "mu=", mu
                    erg = real(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = real(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Re) difference: ", buffer, " : FAILED"
                        print format, " A = ", real(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", real(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Re) OK"
                    end if

                    erg = aimag(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = aimag(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Im) difference: ", buffer, " : FAILED"
                        print format, " A = ", aimag(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", aimag(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Im) OK"
                    end if
                    print *, ""
                end do
                end do
                end do
                end do

                print *, "  B:"
                do n=n_range(1), n_range(2)
                do m=-n, n
                do nu=nu_range(1), nu_range(2)
                do mu=-nu, nu
                    print *, "  n=", n, "m=", m, "nu=", nu, "mu=", mu
                    erg = real(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = real(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg+1))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Re) difference: ", buffer, " : FAILED"
                        print format, " B = ", real(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", real(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Re) OK"
                    end if

                    erg = aimag(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = aimag(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg+1))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Im) difference: ", buffer, " : FAILED"
                        print format, " B = ", aimag(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", aimag(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Im) OK"
                    end if
                    print *, ""
                end do
                end do
                end do
                end do
!                rv = .true.
!                print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r:"
!                print *, "  A:"
!                do i=1, d
!                    erg = real(A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = real(ground_truth_A(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Re) difference: ", buffer, " : FAILED"
!                        print *, "    ", A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Re) OK"
!                    end if
!
!                    erg = aimag(A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = aimag(ground_truth_A(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Im) difference: ", buffer, " : FAILED"
!                        print *, "    ", A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Im) OK"
!                    end if
!                    print *, ""
!                end do
!
!                print *, "  B:"
!                do i=1, d
!                    erg = real(B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg+1))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = real(ground_truth_B(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Re) difference: ", buffer, " : FAILED"
!                        print *, "    ", B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Re) OK"
!                    end if
!
!                    erg = aimag(B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg+1))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = aimag(ground_truth_B(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Im) difference: ", buffer, " : FAILED"
!                        print *, "    ", B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Im) OK"
!                    end if
!                    print *, ""
!                end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name

                    rv = .false.
                    return
                end if

            end function test_lib_math_vector_spherical_harmonics_translation_coeff_r

            function test_lib_math_vector_spherical_harmonics_translation_coeff_r_v2() result(rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 4
                character(len=*), parameter :: file_name = &
                        "ground_truth/lib_math_vector_spherical_harmonics/ground_truth_AB_mnmunu_minus_phi.csv"
                character(len=*), parameter :: format = '(10(X)A, ES19.9)'
                ! auxiliary
                integer :: i
                double precision :: buffer
                double complex :: buffer_cmplx
                double precision :: x
                double precision :: theta
                double precision :: phi
                integer(kind=4), dimension(2) :: n_range
                integer(kind=4), dimension(2) :: nu_range
                integer(kind=1) :: z_selector

                type(list_4_cmplx) :: A_mnkl
                type(list_4_cmplx) :: B_mnkl


                integer :: m
                integer :: n
                integer :: mu
                integer :: nu
                type(list_4_cmplx) :: ground_truth_A
                type(list_4_cmplx) :: ground_truth_B

                double precision :: erg
                double precision :: ground_truth_erg
                integer :: erg_exponent
                double precision :: erg_mantissa

                double precision, dimension(:,:), allocatable :: csv_data
                integer :: csv_columns

                z_selector = 3
                x = 2.0
                theta = 0.5
                phi = 2 * PI - 0.5

                n_range  = (/ 1, 11 /)
                nu_range = (/ 1, 12 /)

                call init_list(ground_truth_A, n_range(1), n_range(2)-n_range(1)+1, &
                                               nu_range(1), nu_range(2)-nu_range(1)+1)
                call init_list(ground_truth_B, n_range(1), n_range(2)-n_range(1)+1, &
                                               nu_range(1), nu_range(2)-nu_range(1)+1)
                n_range  = (/ 1, 5 /)
                nu_range = (/ 1, 5 /)
                ! load ground truth
                if (file_exists(file_name)) then
                    csv_columns = 8
                    call read_csv(file_name, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        m = int(csv_data(i, 1))
                        n = int(csv_data(i, 2))
                        mu = int(csv_data(i, 3))
                        nu = int(csv_data(i, 4))

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_A%item(n)%item(m)%item(nu)%item(mu) = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_B%item(n)%item(m)%item(nu)%item(mu) = buffer_cmplx
                    end do
!                m =  (/ -2,  8,  0, -2 /)
!                n =  (/ 11, 10, 10,  6 /)
!                mu = (/ 3 , -9,  0, -2 /)
!                nu = (/ 9 , 12, 10, 10 /)
!                m = (/ -2, 8 /)
!                n = (/ 11, 10 /)
!                mu = (/ 3, -9 /)
!                nu = (/ 9, 12 /)
!                ground_truth_A(1) = cmplx(.7726121583d+12, .103425m5820d+13, kind=8)
!                ground_truth_B(1) = cmplx(.1222239141d+11, -0.9130398908d+10, kind=8)
!
!                ground_truth_A(2) = cmplx(.3663964990d+35, -.2762412192d+35, kind=8)
!                ground_truth_B(2) = cmplx(-.8370892023d+32, -.1110285257d+33, kind=8)
!
!                ground_truth_A(3) = cmplx(.2969682019d+00, -.1928601440d+18, kind=8)
!                ground_truth_B(3) = cmplx(0.0, 0.0, kind=8)
!
!                ground_truth_A(4) = cmplx(.1377011649d-01, .2385575934d+13, kind=8)
!                ground_truth_B(4) = cmplx(-.3282035237d+12, .1587043209d-02, kind=8)
!
!                n_range = (/ 1, maxval(n) /)
!                nu_range = (/ 1, maxval(nu) /)
                call lib_math_vector_spherical_harmonics_translation_coeff_r(x, theta, phi, &
                                                                             n_range, nu_range, z_selector, &
                                                                             A_mnkl, B_mnkl)

                rv = .true.
                print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r:"
                print *, "  A:"
                do n=n_range(1), n_range(2)
                do m=-n, n
                do nu=nu_range(1), nu_range(2)
                do mu=-nu, nu
                    print *, "  n=", n, "m=", m, "nu=", nu, "mu=", mu
                    erg = real(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = real(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Re) difference: ", buffer, " : FAILED"
                        print format, " A = ", real(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", real(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Re) OK"
                    end if

                    erg = aimag(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = aimag(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Im) difference: ", buffer, " : FAILED"
                        print format, " A = ", aimag(A_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", aimag(ground_truth_A%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Im) OK"
                    end if
                    print *, ""
                end do
                end do
                end do
                end do

                print *, "  B:"
                do n=n_range(1), n_range(2)
                do m=-n, n
                do nu=nu_range(1), nu_range(2)
                do mu=-nu, nu
                    print *, "  n=", n, "m=", m, "nu=", nu, "mu=", mu
                    erg = real(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = real(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg+1))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Re) difference: ", buffer, " : FAILED"
                        print format, " B = ", real(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", real(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Re) OK"
                    end if

                    erg = aimag(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                    ground_truth_erg = aimag(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))

                    erg_exponent = int(log(abs(erg+1))/log(10D0))
                    erg_mantissa = erg / 10d0**erg_exponent
                    buffer = ground_truth_erg / 10d0**erg_exponent - erg_mantissa
                    if (abs(buffer) .gt. ground_truth_e) then
                        print *, "    (Im) difference: ", buffer, " : FAILED"
                        print format, " B = ", aimag(B_mnkl%item(n)%item(m)%item(nu)%item(mu))
                        print format, "GT = ", aimag(ground_truth_B%item(n)%item(m)%item(nu)%item(mu))
                        rv = .false.
                    else
                        print *, "    (Im) OK"
                    end if
                    print *, ""
                end do
                end do
                end do
                end do
!                rv = .true.
!                print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r:"
!                print *, "  A:"
!                do i=1, d
!                    erg = real(A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = real(ground_truth_A(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Re) difference: ", buffer, " : FAILED"
!                        print *, "    ", A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Re) OK"
!                    end if
!
!                    erg = aimag(A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = aimag(ground_truth_A(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Im) difference: ", buffer, " : FAILED"
!                        print *, "    ", A_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Im) OK"
!                    end if
!                    print *, ""
!                end do
!
!                print *, "  B:"
!                do i=1, d
!                    erg = real(B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg+1))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = real(ground_truth_B(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Re) difference: ", buffer, " : FAILED"
!                        print *, "    ", B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Re) OK"
!                    end if
!
!                    erg = aimag(B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i)))
!                    erg_exponent = int(log(abs(erg+1))/log(10D0))
!                    erg_mantissa = erg / 10d0**erg_exponent
!                    buffer = aimag(ground_truth_B(i)) / 10d0**erg_exponent - erg_mantissa
!                    if (abs(buffer) .gt. ground_truth_e) then
!                        print *, "    ", i, " (Im) difference: ", buffer, " : FAILED"
!                        print *, "    ", B_mnkl%item(n(i))%item(m(i))%item(nu(i))%item(mu(i))
!                        rv = .false.
!                    else
!                        print *, "    ", i , ": (Im) OK"
!                    end if
!                    print *, ""
!                end do
                else
                    print *, "test_lib_math_vector_spherical_harmonics_translation_coeff_r: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name

                    rv = .false.
                    return
                end if

            end function test_lib_math_vector_spherical_harmonics_translation_coeff_r_v2

        end function lib_math_vector_spherical_harmonics_test_functions

end module lib_math_vector_spherical_harmonics
