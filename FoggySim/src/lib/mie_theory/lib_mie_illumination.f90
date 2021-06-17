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

module lib_mie_illumination
    use libmath

    use lib_field_polarisation
    use lib_field_plane_wave
    use lib_field_gaussian_beam
    implicit none

    private

    ! public functions
    public :: lib_mie_illumination_get_p_q_j_j
    public :: lib_mie_illumination_init_plane_wave
    public :: lib_mie_illumination_destructor

    ! public  types
    public :: lib_mie_illumination_parameter
    public :: lib_mie_illumination_plane_wave_parameter

!    interface lib_mie_illumination_get_p_q_j_j
!        module procedure lib_mie_illumination_get_p_q_j_j_single_plane_wave
!        module procedure lib_mie_illumination_get_p_q_j_j_multi_plane_wave
!    end interface

    type lib_mie_illumination_plane_wave_parameter
        ! ratio of the e-field of this plane wave to the "global" e-field:  e_x_field_0 / e_field_0
        double precision :: g = 1
        type(cartesian_coordinate_real_type) :: d_0_i ! [m]

        type(lib_field_plane_wave_type) :: beam_parameter
!         ! |k| = 2 Pi / lambda, wave_vector = [1/m]
!        type(cartesian_coordinate_real_type) :: wave_vector_0
!        ! true: transversal electric (TE) mode, false: tranversal magnetic (TM) mode
!        logical :: te_mode = .false.
!        ! additional phase: exp(I phase)
!        double precision :: phase = 0
    end type lib_mie_illumination_plane_wave_parameter

    type lib_mie_illumination_gaussian_beam_parameter
        ! ratio of the e-field of this Gaussian beam to the "global" e-field:  e_x_field_0 / e_field_0
        ! g = (0..1]
        double precision :: g = 1
        ! position of the origin of the beam coordinate system [m]
        type(cartesian_coordinate_real_type) :: d_0_i
!        ! popagation direction: polar angle [0, Pi] [rad]
!        double precision :: theta = 0
!        ! propagation direction: azimuthal angle [0, 2 Pi) [rad]
!        double precision :: phi = 0
        ! beam parameter
        type(lib_field_gaussian_beam_hermite_type) :: beam_parameter
        ! 1: Localised Approximation
        !       restictions:, TEM 00theta = 0, phi = 0
        !           - circular beam: w0 = wx0
        !           - mode: TEM 00
        !           - popagation direction: along the z axis
        !               - thete = 0, phi = 0
        ! 2: Localised Approximation + additions theorem
        !       restictions:, TEM 00theta = 0, phi = 0
        !           - circular beam: w0 = wx0
        !           - mode: TEM 00
        !           - popagation direction: along the z axis
        !               - thete = 0, phi = 0
        ! 3: plane wave approximation
        !       HINT:
        !           - beam_parameter%e_field_0 = 1, will be set automatically
        !             otherwise the scaling factor "g" doesn't work correctly
        !             -> e_field = e_field_0 * g * beam_parameter%e_field_0
        integer :: calculation_type
        ! necessary for calculation_type = 2
        integer(kind=1) :: z_selector_translation
    end type lib_mie_illumination_gaussian_beam_parameter

    ! Illumination Parameter
    ! ----
    !
    ! Plane Wave
    ! -----
    !             _________
    !             ___k^____
    !             ____|____
    !             _________ plane wave
    !                 z
    !                 ^
    !             K_i |
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_i
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_i: illumination coordinate system
    !
    !
    ! Gaussian Beam
    ! -----
    !          _______________
    !            ___________
    !             ___k,z___
    !               __^__   Gaussian beam
    !                _|_
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_i
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_i: illumination coordinate system
    !
    type lib_mie_illumination_parameter
        double precision :: e_field_0 ! [V/m]
        double precision :: lambda_0 ! wave_length_vaccum
        type(lib_mie_illumination_plane_wave_parameter), dimension(:), allocatable :: plane_wave
        type(lib_mie_illumination_gaussian_beam_parameter), dimension(:), allocatable :: gaussian_beam
    end type lib_mie_illumination_parameter

    ! --- caching ---
    type cache_coefficients_p_0_q_0_plane_wave_type
        double precision :: alpha
        double precision :: beta
        integer :: n_max
        type(list_list_cmplx) :: p_0
        type(list_list_cmplx) :: q_0
    end type

    type(cache_coefficients_p_0_q_0_plane_wave_type), dimension(:), allocatable :: cache_coefficients_p_0_q_0_plane_wave
    logical :: cache_coefficients_p_0_q_0_plane_wave_enabled = .false.
    ! ~~~ caching ~~~

    contains

        ! Argument
        ! ----
        !   alpha: double precision, dimension(:)
        !       polar angle [0, Pi)
        !   beta: double precision, dimension(:)
        !       azimuthal angle [0, 2 Pi)
        !   n_max: integer, dimension(:)
        !       maximum degree of the polynomials
        subroutine lib_mie_illumination_init_plane_wave(alpha, beta, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: alpha
            double precision, dimension(lbound(alpha, 1):ubound(alpha, 1)), intent(in) :: beta
            integer(kind=4), dimension(lbound(alpha, 1):ubound(alpha, 1)) :: n_max

            ! auxiliary
            integer :: i

            ! --- init: cache_coefficients_a_b_cmplx_barberh_x ---
            if (allocated(cache_coefficients_p_0_q_0_plane_wave)) then
                deallocate(cache_coefficients_p_0_q_0_plane_wave)
            end if

            allocate(cache_coefficients_p_0_q_0_plane_wave(size(alpha)))
            cache_coefficients_p_0_q_0_plane_wave_enabled = .false.

            do i=1, size(alpha)
                cache_coefficients_p_0_q_0_plane_wave(i)%alpha = alpha(i)
                cache_coefficients_p_0_q_0_plane_wave(i)%beta = beta(i)
                cache_coefficients_p_0_q_0_plane_wave(i)%n_max = n_max(i)

                call get_p_q_j_j_plane_wave_core(alpha(i), beta(i), (/1, n_max(i)/), &
                                                 cache_coefficients_p_0_q_0_plane_wave(i)%p_0, &
                                                 cache_coefficients_p_0_q_0_plane_wave(i)%q_0, &
                                                 caching=.false.)
            end do

            cache_coefficients_p_0_q_0_plane_wave_enabled = .true.

        end subroutine lib_mie_illumination_init_plane_wave

        subroutine lib_mie_illumination_destructor
            implicit none

            if (allocated(cache_coefficients_p_0_q_0_plane_wave)) then
                deallocate(cache_coefficients_p_0_q_0_plane_wave)
            end if
            cache_coefficients_p_0_q_0_plane_wave_enabled = .false.
        end subroutine
        subroutine lib_mie_illumination_get_p_q_j_j(illumination, n_medium, d_0_j, n_range, p_nm, q_nm, caching)

            implicit none
            ! dummy
            type(lib_mie_illumination_parameter), intent(in) :: illumination
            double precision :: n_medium
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            logical, intent(in), optional :: caching

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            ! auxiliary
            logical :: m_caching
            type(list_list_cmplx) :: buffer_p_nm
            type(list_list_cmplx) :: buffer_q_nm

            m_caching = .true.
            if (present(caching)) m_caching = caching

            call init_list(p_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))
            call init_list(q_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))

            if (allocated(illumination%plane_wave)) then
                call lib_mie_illumination_get_p_q_j_j_multi_plane_wave(illumination%plane_wave, n_medium, d_0_j, n_range, &
                                                                       buffer_p_nm, buffer_q_nm, &
                                                                       caching = m_caching)
                call move_alloc(buffer_p_nm%item, p_nm%item)
                call move_alloc(buffer_q_nm%item, q_nm%item)
            end if

            if (allocated(illumination%gaussian_beam)) then
                call lib_mie_illumination_get_p_q_gaussian_beam(illumination%gaussian_beam, n_medium, d_0_j, n_range,&
                                                                         buffer_p_nm, buffer_q_nm)

                call remove_zeros(buffer_p_nm, dcmplx(1d-290, 1d-290), .true.)
                call remove_zeros(buffer_q_nm, dcmplx(1d-290, 1d-290), .true.)

                p_nm = p_nm + buffer_p_nm
                q_nm = q_nm + buffer_q_nm
            end if

        end subroutine

        ! Calculates the coefficients (of vector spherical components) of a plane incident wave
        ! - transverse magnetic mode (TM) (std)
        ! - transverse electric mode (TE)
        !
        !                 z
        !                 ^
        !             K_j |
        !                 --> x
        !                ^
        !               /
        !           z  /d_j_0
        !           ^ /
        !       K_0 |/
        !           --> x
        !      ____________
        !      ____k^______
        !      _____|______
        !      ____________ plane wave
        !
        ! Argument
        ! ----
        !   k: type(cartesian_coordinate_real_type)
        !       wave vector
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       vector from the
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   te_mode: logical (std: .false.)
        !       false: TM mode referred to the incident plane (k-vector, z-axis)
        !       true: TE mode referred to the incident plane (k-vector, z-axis)
        !   phase: double precision (std: 0)
        !       add an additional phase
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !   q: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine lib_mie_illumination_get_p_q_j_j_single_plane_wave(k, d_0_j, n_range, p, q, &
                                                                      te_mode, phase, caching)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: k
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p
            type(list_list_cmplx), intent(inout) :: q

            logical, optional :: te_mode
            double precision, optional :: phase
            logical, optional :: caching

            ! auxiliary
            double precision :: alpha
            double precision :: beta

            type(spherical_coordinate_real_type) :: k_spherical

            double precision :: dot_k_d

            complex(kind=8) :: exp_k_d

            double precision :: m_phase

            logical :: m_te_mode
            type(list_list_cmplx) :: buffer

            m_phase = 0
            if (present(phase)) m_phase = phase

            m_te_mode = .false.
            if (present(te_mode)) m_te_mode = te_mode

            ! --- pre-calc ---
            k_spherical = k

            alpha = k_spherical%theta
            beta = k_spherical%phi

            if (abs(d_0_j) .gt. 0.0 .or. m_phase .ne. 0) then
                ! j-th coordinate system
                dot_k_d = dot_product(k, d_0_j) + m_phase
                exp_k_d = cmplx(cos(dot_k_d), sin(dot_k_d), kind=8)
                call get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, exp_k_d, caching=caching)
            else
                ! 0-th coordinate system
                call get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, caching=caching)
            end if

            if (m_te_mode) then
                call move_alloc(p%item, buffer%item)

                call move_alloc(q%item, p%item)
                call move_alloc(buffer%item, q%item)
            end if

        end subroutine lib_mie_illumination_get_p_q_j_j_single_plane_wave

        ! Calculates the coefficients (of vector spherical components) of multiple plane incident waves
        ! - transverse magnetic mode (TM)
        !
        !                 z
        !                 ^
        !             K_j |
        !                 --> x
        !                ^
        !               /
        !           z  /d_j_0
        !           ^ /
        !       K_0 |/
        !           --> x
        !      ____________
        !      ____k^______
        !      _____|______
        !      ____________ plane wave
        !
        ! Argument
        ! ----
        !   k: type(cartesian_coordinate_real_type)
        !       wave vector
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       vector from the
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !   q: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !
        subroutine lib_mie_illumination_get_p_q_j_j_multi_plane_wave(illumination, n_medium, d_0_j, n_range, p_nm, q_nm, caching)
            implicit none
            ! dummy
            type(lib_mie_illumination_plane_wave_parameter), dimension(:), intent(in) :: illumination
            double precision :: n_medium
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            logical, optional :: caching

            ! auxiliary
            integer :: i

            type(cartesian_coordinate_real_type) :: k
            type(spherical_coordinate_real_type) :: k_spherical
            type(cartesian_coordinate_real_type) :: d_i_j

            type(list_list_cmplx) :: buffer_a
            type(list_list_cmplx) :: buffer_b

            type(list_list_cmplx), dimension(size(illumination)) :: buffer_p_nm
            type(list_list_cmplx), dimension(size(illumination)) :: buffer_q_nm

            double precision :: buffer_abs
            double precision :: buffer_arg

            !$OMP PARALLEL DO PRIVATE(i, d_i_j, k, buffer_a, buffer_b, buffer_abs, buffer_arg)
            do i = 1, size(illumination)
                d_i_j = d_0_j - illumination(i)%d_0_i
                k_spherical = make_spherical(2d0 * PI * n_medium / illumination(i)%beam_parameter%wave_length_0, &
                                             illumination(i)%beam_parameter%theta, &
                                             illumination(i)%beam_parameter%phi)
                k = k_spherical

                ! calculate TM
                buffer_abs = abs(illumination(i)%beam_parameter%polarisation%x)
                buffer_arg = atan2(aimag(illumination(i)%beam_parameter%polarisation%x), &
                                   real(illumination(i)%beam_parameter%polarisation%x))
                buffer_arg = buffer_arg + illumination(i)%beam_parameter%phase

                if (buffer_abs .gt. 0) then
                    call lib_mie_illumination_get_p_q_j_j_single_plane_wave(k, d_i_j, n_range, &
                                                                            buffer_a, buffer_b,&
                                                                            phase = buffer_arg, &
                                                                            te_mode = .false., &
                                                                            caching=caching)

                    buffer_p_nm(i) = buffer_abs * buffer_a
                    buffer_q_nm(i) = buffer_abs * buffer_b
                else
                    call init_list(buffer_p_nm(i), 0, 1, dcmplx(0, 0))
                    call init_list(buffer_q_nm(i), 0, 1, dcmplx(0, 0))
                end if

                ! calculate TE
                buffer_abs = abs(illumination(i)%beam_parameter%polarisation%y)
                buffer_arg = atan2(aimag(illumination(i)%beam_parameter%polarisation%y), &
                                   real(illumination(i)%beam_parameter%polarisation%y))
                buffer_arg = buffer_arg + illumination(i)%beam_parameter%phase

                if (buffer_abs .gt. 0) then
                    call lib_mie_illumination_get_p_q_j_j_single_plane_wave(k, d_i_j, n_range, &
                                                                            buffer_a, buffer_b,&
                                                                            phase = buffer_arg, &
                                                                            te_mode = .true., &
                                                                            caching=caching)

                    buffer_p_nm(i) = buffer_p_nm(i) + buffer_abs * buffer_a
                    buffer_q_nm(i) = buffer_q_nm(i) + buffer_abs * buffer_b
                end if


            end do
            !$OMP END PARALLEL DO

            call init_list(p_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(q_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

            do i = 1, size(illumination)
                p_nm = p_nm + illumination(i)%g * buffer_p_nm(i)
                q_nm = q_nm + illumination(i)%g * buffer_q_nm(i)
            end do

        end subroutine lib_mie_illumination_get_p_q_j_j_multi_plane_wave

!        ! Argument
!        ! ----
!        !   illumination: type(lib_mie_illumination_parameter)
!        !       illumination parameter
!        !   d_0_j: type(cartesian_coordinate_real_type)
!        !       vector from the
!        !   n_range: integer, dimension(2)
!        !       first and last index (degree) of the sum to calculate the electical field
!        !       e.g. n = (/ 1, 45 /) <-- size parameter
!        !
!        ! Returns
!        ! ----
!        !   p: type(list_list_cmplx)
!        !       coefficient of vector spherical componets
!        !   q: type(list_list_cmplx)
!        !       coefficient of vector spherical componets
!        !
!        ! Reference: Expansion of an arbitrarily oriented, located, and shaped beam in spheroidal coordinates, Feng Xu
!        !            eq. 6, 7
!        subroutine lib_mie_illumination_get_p_q_j_j_arbitrary_field(illumination, n_medium, d_0_j, n_range, p_nm, q_nm)
!            implicit none
!            ! dummy
!            type(lib_mie_illumination_parameter), intent(in) :: illumination
!            double precision :: n_medium
!            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
!            integer(kind=4), dimension(2),intent(in) :: n_range
!
!            type(list_list_cmplx), intent(inout) :: p_nm
!            type(list_list_cmplx), intent(inout) :: q_nm
!
!            ! auxiliaray
!            integer :: n
!            integer :: m
!
!            integer :: z_selector
!
!            double precision :: prefactor
!            double complex :: prefactor_r
!
!            double precision :: kr
!
!
!            prefactor = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n, abs(m)) * kr / (4d0 * PI)
!
!            prefactor_r = dcmplx(0,1)**(n+1) * prefactor / get_z_function(x, selector, n_range)
!
!
!
!
!
!        end subroutine lib_mie_illumination_get_p_q_j_j_arbitrary_field
!
!        function get_z_function(x, selector, n_range) result (rv)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: x
!            integer, intent(in) :: selector
!            integer, dimension(2), intent(in) :: n_range
!
!            type(list_cmplx) :: rv
!
!
!            ! auxiliary
!            integer :: n
!            double precision, dimension(:), allocatable :: z_n_real
!            double complex, dimension(:), allocatable :: z_n_cmplx
!
!            ! z function
!            select case (z_selector)
!                case(1)
!                    ! spherical Bessel function first kind j_n
!                    allocate(z_n_real(0:n_range(2)+1))
!                    z_n_real = lib_math_bessel_spherical_first_kind(x, 0, n_range(2)+1)
!                case(2)
!                    ! spherical Bessel function second kind y_n
!                    allocate(z_n_real(0:n_range(2)+1))
!                    z_n_real = lib_math_bessel_spherical_second_kind(x, 0, n_range(2)+1)
!                case(3)
!                    ! spherical Hankel function first kind   h^(1)_n
!                    allocate(z_n_cmplx(0:n_range(2)+1))
!                    z_n_cmplx = lib_math_hankel_spherical_1(x, 0, n_range(2)+1)
!                case(4)
!                    ! spherical Hankel function second kind   h^(2)_n
!                    allocate(z_n_cmplx(0:n_range(2)+1))
!                    z_n_cmplx = lib_math_hankel_spherical_2(x, 0, n_range(2)+1)
!                case default
!                    z_n_real = 0
!                    z_n_cmplx = cmplx(0,0)
!
!                    print*, "get_z_function: ERROR"
!                    print*, "  undefined z_selector value[1-4]: ", z_selector
!                    return
!            end select
!
!            call init_list(rv, n_range(1), n_range(2)-n_range(1)+1)
!
!            select case (z_selector)
!                case(1, 2)
!                    do n = n_range(1), n_range(2)
!                        rv%item(n) = dcmplx(z_n_real(n), 0)
!                    end do
!
!                    deallocate(z_n_real)
!                case(3, 4)
!                    do n = n_range(1), n_range(2)
!                        rv%item(n) = z_n_cmplx(n)
!                    end do
!
!                    deallocate(z_n_cmplx)
!            end select
!        end function



        ! Calculates the coefficients (of vector spherical components) of a plane incident wave
        ! - transverse magnetic mode (TM)
        !
        !                 z
        !                 ^
        !             K_j |
        !                 --> x
        !                ^
        !               /
        !           z  /d_j_0
        !           ^ /
        !       K_0 |/
        !           --> x
        !      ____________
        !      ____k^______
        !      _____|______
        !      ____________ plane wave
        !
        ! Argument
        ! ----
        !   alpha: double precision
        !       polar angle [0..Pi)
        !       if 0: propagation direction along the z-axis
        !   beta: double precision
        !       azimuthal angle [0..2Pi)
        !       if 0 and alpha = 0: The E-field oscillates along the x-axis.
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   exp_k_d: complex, optional (std: 1.0, 0-th coordinate system)
        !       formula: exp(i dot(k, d_0_j))
        !           k: wave vector
        !           d_0_j: vector from 0-th to j-th coordinate system
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !
        !   q: type(list_list_cmplx)
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, exp_k_d, caching)
            implicit none
            ! dummy
            double precision, intent(in) :: alpha
            double precision, intent(in) :: beta
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p
            type(list_list_cmplx), intent(inout) :: q

            complex(kind=8), intent(in), optional :: exp_k_d

            logical, optional :: caching

            ! auxiliary
            integer :: i
            double precision :: buffer_real_n
            integer(kind=4) :: m
            integer(kind=4) :: n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            double precision :: sin_beta
            double precision :: cos_beta

            logical :: m_caching
            ! 0: p_0 and q_0 are not pre calculated
            ! >0: element number of the cache array
            integer :: cache_no

            if (cache_coefficients_p_0_q_0_plane_wave_enabled) then
                if (present(caching)) then
                    m_caching = caching
                else
                    m_caching = .true.
                end if
            else
                m_caching = .false.
            end if

            if (m_caching) then
                do i=1, size(cache_coefficients_p_0_q_0_plane_wave)
                    if ((cache_coefficients_p_0_q_0_plane_wave(i)%alpha .eq. alpha) &
                        .and. (cache_coefficients_p_0_q_0_plane_wave(i)%beta .eq. beta) &
                        .and. (cache_coefficients_p_0_q_0_plane_wave(i)%n_max .le. n_range(2))) then
                        cache_no = i
                    else
                        cache_no = 0
                    end if
                end do
            else
                cache_no = 0
            end if

            ! --- init ---
!            if (alpha .eq. 0.d0) then
!                allocate(p%item(n_range(1):n_range(2)))
!                allocate(q%item(n_range(1):n_range(2)))
!                do i=n_range(1), n_range(2)
!                    allocate(p%item(i)%item(-1:1))
!
!                    allocate(q%item(i)%item(-1:1))
!                end do
!            else
                call init_list(p, n_range(1), n_range(2)-n_range(1)+1)
                call init_list(q, n_range(1), n_range(2)-n_range(1)+1)
!            end if

            ! --- pre-calc ---
            if (cache_no .eq. 0) then
                sin_beta = sin(beta)
                cos_beta = cos(beta)

                call lib_math_associated_legendre_polynomial_theta(alpha, n_range(2), pi_nm, tau_nm)
            end if

            ! errata eq. (1) Equations (21) on p. 4577
            do n=n_range(1), n_range(2)
                buffer_real_n = 1.0d0 / real(n*(n+1), kind=8)
                if (alpha .eq. 0.d0) then
                    p%item(n)%item(:) = 0.d0
                    q%item(n)%item(:) = 0.d0
                    call get_value(-1_4, n, p%item(n)%item(-1), q%item(n)%item(-1))
                    call get_value(1_4, n, p%item(n)%item(1), q%item(n)%item(1))
                else
                    do m=-n, n
                        call get_value(m, n, p%item(n)%item(m), q%item(n)%item(m))
                    end do
                end if
            end do

            contains
                ! errata eq. (1) Equations (21) on p. 4577
                subroutine get_value(m, n, p, q)
                    implicit none
                    ! dummy
                    integer(kind=4), intent(in) :: m
                    integer(kind=4), intent(in) :: n

                    complex(kind=lib_math_type_kind), intent(inout) :: p
                    complex(kind=lib_math_type_kind), intent(inout) :: q

                    ! auxiliaray
                    double precision :: buffer_real
                    double complex :: buffer_cmplx

                    if (cache_no .gt. 0) then
                        p = cache_coefficients_p_0_q_0_plane_wave(cache_no)%p_0%item(n)%item(m)
                        q = cache_coefficients_p_0_q_0_plane_wave(cache_no)%q_0%item(n)%item(m)
                    else
                        buffer_real = -m * beta
                        buffer_cmplx = cmplx(cos(buffer_real), sin(buffer_real), kind=8)

                        buffer_cmplx = buffer_real_n * buffer_cmplx

                        p = tau_nm%item(n)%item(m) * buffer_cmplx
                        q = pi_nm%item(n)%item(m) * buffer_cmplx
                    end if

                    if (present(exp_k_d)) then
                        p = exp_k_d * p
                        q = exp_k_d * q
                    end if

                end subroutine

        end subroutine get_p_q_j_j_plane_wave_core

        ! Calculates the beam shape coefficients of a Gaussian beam
        !
        ! Argument
        ! ----
        !   illumination: type(lib_mie_illumination_parameter)
        !   n_medium: double precision
        !       refractive index of the medium
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       point of evaluation (e.g. location of the sphere)
        !   n_range: integer, dimension(2)
        !
        ! Returns
        ! ----
        !   p_nm, q_nm: type(list_list_cmplx)
        !       beam shape coefficients
        !
        subroutine lib_mie_illumination_get_p_q_gaussian_beam(illumination, n_medium, d_0_j, n_range, p_nm, q_nm)
            implicit none
            ! dummy
            type(lib_mie_illumination_gaussian_beam_parameter), dimension(:), intent(in) :: illumination
            double precision :: n_medium
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            ! auxiliary
            integer :: i
            integer :: n
            integer :: m
            integer :: first_beam
            integer :: last_beam
            integer :: calculation_type

            double precision :: wave_length
            double precision :: w0
            type(cartesian_coordinate_real_type) :: d

            type(list_list_cmplx) :: buffer_a_nm
            type(list_list_cmplx) :: buffer_b_nm

            type(list_list_cmplx), dimension(:), allocatable :: buffer_p_nm
            type(list_list_cmplx), dimension(:), allocatable :: buffer_q_nm

            integer(kind=1) :: z_selector

            first_beam = lbound(illumination, 1)
            last_beam = ubound(illumination, 1)

            allocate(buffer_p_nm(first_beam:last_beam))
            allocate(buffer_q_nm(first_beam:last_beam))


            do i = first_beam, last_beam
                wave_length = illumination(i)%beam_parameter%wave_length_0

                calculation_type = illumination(i)%calculation_type


                select case(illumination(i)%calculation_type)
                    case (1)
                        ! localised approximation
                        ! Reference: Generalized Lorenz-MieTheories, Gérard Gouesbet Gérard Gréhan, program GNMF
                        w0 = illumination(i)%beam_parameter%waist_x0
                        d = d_0_j - illumination(i)%d_0_i

                        call lib_mie_illumination_get_p_q_single_gaussian_beam_la(wave_length, n_medium, w0, d, n_range,&
                                                                                  buffer_a_nm, buffer_b_nm)

                        buffer_p_nm(i) = illumination(i)%g * buffer_a_nm
                        buffer_q_nm(i) = illumination(i)%g * buffer_b_nm
                    case (2)
                        ! localised approximation + addition theorem
                        ! Reference: Generalized Lorenz-MieTheories, Gérard Gouesbet Gérard Gréhan, program GNMF
                        w0 = illumination(i)%beam_parameter%waist_x0
                        d%x = 0d0
                        d%y = 0d0
                        d%z = 0d0

                        call lib_mie_illumination_get_p_q_single_gaussian_beam_la(wave_length, n_medium, w0, d, n_range,&
                                                                                  buffer_a_nm, buffer_b_nm)

                        d = (illumination(i)%d_0_i - d_0_j) * 2d0 * PI * n_medium / wave_length

                        z_selector = illumination(i)%z_selector_translation

                        call lib_math_vector_spherical_harmonics_translate_coefficient(buffer_a_nm, buffer_b_nm, &
                                                                                       d, &
                                                                                       n_range, n_range, &
                                                                                       z_selector, &
                                                                                       buffer_p_nm(i), buffer_q_nm(i))

                        buffer_p_nm(i) = illumination(i)%g * buffer_p_nm(i)
                        buffer_q_nm(i) = illumination(i)%g * buffer_q_nm(i)
                    case (3)
                        ! plane wave approximation

                        ! HINT: beam_parameter%e_field_0 = 1, will be set automatically
                        !       otherwise the scaling factor "g" doesn't work correctly
                        !       -> e_field = e_field_0 * g * beam_parameter%e_field_0

                        d = d_0_j - illumination(i)%d_0_i

                        call lib_mie_illumination_get_p_q_single_gaussian_beam_pw(illumination(i), n_medium, &
                                                                                  d, n_range, &
                                                                                  buffer_a_nm, buffer_b_nm)

                        buffer_p_nm(i) = illumination(i)%g * buffer_a_nm
                        buffer_q_nm(i) = illumination(i)%g * buffer_b_nm
                end select
            end do

            call init_list(p_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))
            call init_list(q_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))

            do i = first_beam, last_beam
                p_nm = p_nm + buffer_p_nm(i)
                q_nm = q_nm + buffer_q_nm(i)
            end do

            do n = n_range(1), n_range(2)
                do m = -n, n
                    if (isnan(real(p_nm%item(n)%item(m)))) then
                        print *, "lib_mie_illumination_get_p_q_gaussian_beam: ERROR"
                        print *, "  Re(p_nm) is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if

                    if (isnan(aimag(p_nm%item(n)%item(m)))) then
                        print *, "lib_mie_illumination_get_p_q_gaussian_beam: ERROR"
                        print *, "  Im(p_nm) is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if

                    if (isnan(real(q_nm%item(n)%item(m)))) then
                        print *, "lib_mie_illumination_get_p_q_gaussian_beam: ERROR"
                        print *, "  Re(q_nm) is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if

                    if (isnan(aimag(q_nm%item(n)%item(m)))) then
                        print *, "lib_mie_illumination_get_p_q_gaussian_beam: ERROR"
                        print *, "  Im(q_nm) is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if
                end do
            end do

        end subroutine lib_mie_illumination_get_p_q_gaussian_beam


        ! Calculates the beam shape coefficients of a circular Gaussian beam with the localised approximation method.
        !
        ! Argument
        ! ----
        !   illumination: type(lib_mie_illumination_parameter)
        !   n_medium: double precision
        !       refractive index of the medium
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       point of evaluation (e.g. location of the sphere)
        !   n_range: integer, dimension(2)
        !   esj: double precision, optional (std: 10**-30)
        !       "accuracy of g(n,m). The program stops when the term to add is smaller than esj"
        ! Returns
        ! ----
        !   p_nm, q_nm: type(list_list_cmplx)
        !       beam shape coefficients
        !
        ! Reference: Generalized Lorenz-MieTheories, Gérard Gouesbet Gérard Gréhan, program GNMF
        subroutine lib_mie_illumination_get_p_q_single_gaussian_beam_la(wave_length, n_medium, w0, d, n_range, &
                                                                        p_nm, q_nm, esj)
!            use gnm_mod
            implicit none
            ! dummy
            double precision, intent(in) :: wave_length
            double precision :: n_medium
            double precision :: w0
            type(cartesian_coordinate_real_type), intent(in) :: d
            integer, dimension(2),intent(in) :: n_range

            double precision, intent(in), optional :: esj

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            ! auxiliary
!            integer :: n
!            integer :: m
!            integer :: nmax
!            REAL(KIND=DBL)::s,l,x0ad,y0ad,z0ad,m_esj
!
!            m_esj=1.0d-30
!            if (present(esj)) m_esj = esj
!
!            nmax = n_range(2)
!
!
!            s=wave_length/(w0*2.0_DBL*PI * n_medium)
!            l=w0/s
!            x0ad=d%x/w0
!            y0ad=d%y/w0
!            z0ad=d%z/l
!
!            CALL gnmf(nmax,s,x0ad,y0ad,z0ad,m_esj)
!
!            call init_list(p_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))
!            call init_list(q_nm, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))
!
!            do n = n_range(1), n_range(2)
!                if (n .le. MAX_M) then
!                    m = n
!                else
!                    m = MAX_M
!                end if
!                p_nm%item(n)%item(-m:m) = gte(n, -m:m)
!                q_nm%item(n)%item(-m:m) = gtm(n, -m:m)
!            end do
!
!            do n = n_range(1), n_range(2)
!                do m = -n, n
!                    if (isnan(real(p_nm%item(n)%item(m)))) then
!                        print *, "lib_mie_illumination_get_p_q_single_gaussian_beam_la: ERROR"
!                        print *, "  Re(p_nm) is NaN"
!                        print * , "  n = ", n
!                        print * , "  m = ", m
!                    end if
!
!                    if (isnan(aimag(p_nm%item(n)%item(m)))) then
!                        print *, "lib_mie_illumination_get_p_q_single_gaussian_beam_la: ERROR"
!                        print *, "  Im(p_nm) is NaN"
!                        print * , "  n = ", n
!                        print * , "  m = ", m
!                    end if
!
!                    if (isnan(real(q_nm%item(n)%item(m)))) then
!                        print *, "lib_mie_illumination_get_p_q_single_gaussian_beam_la: ERROR"
!                        print *, "  Re(q_nm) is NaN"
!                        print * , "  n = ", n
!                        print * , "  m = ", m
!                    end if
!
!                    if (isnan(aimag(q_nm%item(n)%item(m)))) then
!                        print *, "lib_mie_illumination_get_p_q_single_gaussian_beam_la: ERROR"
!                        print *, "  Im(q_nm) is NaN"
!                        print * , "  n = ", n
!                        print * , "  m = ", m
!                    end if
!                end do
!            end do

        end subroutine lib_mie_illumination_get_p_q_single_gaussian_beam_la

        subroutine lib_mie_illumination_get_p_q_single_gaussian_beam_pw(beam_parameter, n_medium, d_0_j, n_range, p_nm, q_nm)
            implicit none
            ! dummy
            type(lib_mie_illumination_gaussian_beam_parameter), intent(in) :: beam_parameter
            double precision :: n_medium
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            ! auxiliary
            type(cartesian_coordinate_real_type) :: evaluation_point_x
            type(lib_mie_illumination_plane_wave_parameter) :: plane_wave


            plane_wave%g = beam_parameter%g
            plane_wave%d_0_i = beam_parameter%d_0_i

            evaluation_point_x = d_0_j - beam_parameter%d_0_i
            plane_wave%beam_parameter = lib_field_gaussian_beam_hermite_get_plane_wave_approximation( &
                                                                            beam_parameter%beam_parameter, &
                                                                            evaluation_point_x)

            call lib_mie_illumination_get_p_q_j_j_multi_plane_wave((/ plane_wave /), &
                                                                   n_medium, d_0_j, n_range, &
                                                                   p_nm, q_nm)

        end subroutine

!        !
!        subroutine lib_mie_illumination_get_p_q_elliptical_gaussian_beam(illumination, n_medium, d_0_j,  n_range, p_nm, q_nm)
!            implicit none
!            ! dummy
!            type(lib_mie_illumination_parameter), intent(in) :: illumination
!            double precision :: n_medium
!            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
!            integer(kind=4), dimension(2),intent(in) :: n_range
!
!            type(list_list_cmplx), intent(inout) :: p_nm
!            type(list_list_cmplx), intent(inout) :: q_nm
!
!            ! auxiliary
!            double precision :: X0
!            double precision :: Y0
!            double precision :: Z0
!
!            double precision :: s_x
!            double precision :: s_y
!            double precision :: s_z
!
!
!        end subroutine
!
!        ! Aurgument
!        ! ----
!        !   k: double precision
!        !       wave number k = 2 Pi / lambda [1/m]
!        !   n_medium: double precisio
!        !       refractive index of the medium
!        !   w0x: double precision
!        !       beam waist radius along the x-axis
!        !   w0y: double precision
!        !       beam waist radius along the y-axis
!        !   d_0_i: type(cartesian_coordinate_real_type)
!        subroutine lib_mie_illumination_get_p_q_single_elliptical_gaussian_beam(k, n_medium, w0x, w0y, d, n_range, &
!                                                                                p_nm, q_nm)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: k
!            double precision :: n_medium
!            double precision :: w0x
!            double precision :: w0y
!            type(cartesian_coordinate_real_type), intent(in) :: d
!            integer(kind=4), dimension(2),intent(in) :: n_range
!
!            type(list_list_cmplx), intent(inout) :: p_nm
!            type(list_list_cmplx), intent(inout) :: q_nm
!
!            ! auxiliary
!
!            double precision :: X_0
!            double precision :: Y_0
!            double precision :: Z_0
!
!            double precision :: s_x
!            double precision :: s_y
!            double precision :: s_z
!
!            double precision :: sqrt_Q
!            double precision :: Qsx
!            double precision :: Qsy
!
!            double precision :: QsxX_0
!            double precision :: QsyY_0
!
!            double precision :: buffer_real
!            double complex :: buffer_cmplx
!
!            double precision :: L
!            double precision :: psi ! Psi_00_sh
!            double precision :: F_tilde
!            double precision :: Xi
!            double precision :: D_tilde
!
!            X_0 = k * n_medium * d%x
!            Y_0 = k * n_medium * d%y
!            Z_0 = k * n_medium * d%z
!
!            s_x = 1 / (k * n_medium * w0x)
!            s_y = 1 / (k * n_medium * w0y)
!
!            L = gaussian_beam_get_parameter_L(n, m)
!
!            ! --- eq. 13 ---
!            sqrt_Q = sqrt(Q_x * Q_y)
!            Qsx = Q_x * s_x**2
!            Qsy = Q_y * s_y**2
!
!            QsxX0 = Qsx * X_0
!            QsyY0 = Qsy * Y_0
!
!            QsxX0X0 = QsxX0 * X_0
!            QsyY0Y0 = QsyY0 * Y_0
!
!            buffer_real = Z_0 - (Qsx + Qsy) * L / 2d0 - ( QsxX0X0 - QsyY0Y0 )
!
!            psi = cmplx(0, sqrt_Q) * exp(cmplx(0, buffer_real))
!            ! ~~~ eq. 13 ~~~
!
!            ! --- eq. 17 ---
!            F_tilde = 2 * sqrt(L * (QsxX0X0 + QyY0Y0) )
!            ! ~~~ eq. 17 ~~~
!
!            ! --- eq. 18 ---
!            if (Y_0 .eq. 0d0) then
!                Xi = PI / 2d0 ! acos(0)
!            else if (X_0 .eq. 0d0) then
!                Xi = acos(QsyY0 / abs(QsyY0))
!            else
!                Xi = acos(QsyY0 / sqrt(QsxX0**2 + QsyY0**2))
!            end if
!            ! ~~~ eq. 18 ~~~
!
!            ! --- eq. 30 ---
!
!            ! ~~~ eq. 30 ~~~
!
!
!        end subroutine
!
!        ! Argument
!        ! ----
!        !   n: integer
!        !       degree of the polynomial
!        !   m: integer
!        !       order of the polynomial
!        !   MLA: logical, optional (std: .false.)
!        !       false: Orginal Localised Approximation
!        !       true: Modified Localised Approximation
!        !
!        ! Returns
!        ! ----
!        !   L: double precision
!        !       parameter L
!        !       "The localization operator ˆG changes the radial parameter kr to L^1∕2"
!        !
!        ! refrence: Compact formulation of the beam shapecoefficients for elliptical Gaussian beambased on localized approximation
!        !           JIANQI SHEN,* XIAOWEI JIA,AND HAITAO YU, page 2
!        function gaussian_beam_get_parameter_L(n , m, MLA, square_root) result (L)
!            implicit none
!            ! dummy
!            integer, intent(in) :: n
!            integer, intent(in) :: m
!            logical, intent(in), optional :: MLA
!            logical, intent(in), optional :: square_root
!
!            double precision :: L
!
!            ! auxiliary
!            logical :: m_MLA
!            logical :: m_square_root
!
!            m_MLA = .false.
!            if (present(MLA)) m_MLA = MLA
!
!            if (m_MLA) then
!                L = ( dble(n) + 0.5d0 )**2 - / dble(abs(m)) + 0.5d0 )**2
!            else
!                L = ( dble(n) + 0.5d0 )**2
!            end if
!
!        end function gaussian_beam_get_parameter_L
!
!        ! Argument
!        ! ----
!        !   n: integer
!        !       degree of the polynomial
!        !   m: integer
!        !       order of the polynomial
!        !
!        ! Returns
!        ! ----
!        !   Z: double complex
!        !       normalisation factor Z
!        !
!        ! refrence: Compact formulation of the beam shapecoefficients for elliptical Gaussian beambased on localized approximation
!        !           JIANQI SHEN,* XIAOWEI JIA,AND HAITAO YU, eq. 5
!        function gaussian_beam_get_normalisation_factor_Z(n, m) returns(Z)
!            implicit none
!            ! dummy
!            integer, intent(in) :: n
!            integer, intent(in) :: m
!
!            double complex :: Z
!
!            ! auxiliary
!            double precision :: buffer
!
!            if (m .eq. 0) then
!                buffer = dble(n*(n+1)) / (dble(n) + 0.5d0)
!                Z = cmplx(0, buffer)
!            else
!                buffer = gaussian_beam_get_parameter_L(n,m)**(dble(abs(m)-1) / 2d0 )
!                Z = cmplx(0, -1)**(dble(abs(m)-1) / 2d0 )
!                Z = buffer * Z
!            end if
!        end function gaussian_beam_get_normalisation_factor_Z
!
!        ! refrence: Compact formulation of the beam shapecoefficients for elliptical Gaussian beambased on localized approximation
!        !           JIANQI SHEN,* XIAOWEI JIA,AND HAITAO YU, eq. 15
!        function gaussian_beam_get_Q_xy(Z_0, s_xy) result(Q_xy)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: Z_0
!            double precision, intent(in) :: s_xy
!
!            double precision :: Q_xy
!
!            Q_xy = 1 / (cmplx(-2 * Z_0 * s_xy**2 ,1) )
!
!        end function gaussian_beam_get_Q_xy


end module lib_mie_illumination
