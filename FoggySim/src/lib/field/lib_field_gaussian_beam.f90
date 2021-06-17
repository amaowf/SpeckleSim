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

module lib_field_gaussian_beam
    use libmath
    use lib_constants
    use lib_field_polarisation
    use lib_field_plane_wave
    implicit none

    private

    ! public function
    public :: lib_field_gaussian_beam_hermite_get_field
    public :: lib_field_gaussian_beam_hermite_get_propagation_direction
    public :: lib_field_gaussian_beam_hermite_get_plane_wave_approximation

    public :: lib_field_gaussian_beam_test_functions

    ! public types
    public :: lib_field_gaussian_beam_hermite_type

    type lib_field_gaussian_beam_hermite_type
        double precision :: e_field_0                   ! [V/m]
        double precision :: wave_length_0               ! [m]
        double precision :: refractive_index_medium     ! [1]
        double precision :: phase = 0
        double precision :: waist_x0                    ! beam waist radius along the x-axis [m]
        double precision :: waist_y0                    ! beam waist radius along the y-axis [m]
        integer :: tem_m = 0                            ! Hermite-Gaussian Mode along the x axis
        integer :: tem_n = 0                            ! Hermite-Gaussian Mode along the y axis
        type(jones_vector_type) :: polarisation
        double precision :: theta = 0                   ! propagation direction: polar angle [0, Pi] [rad]
        double precision :: phi = 0                     ! propagation direction: azimuthal angle [0, 2 Pi) [rad]
        ! 1: exp(-i(k*z - omega*t))
        ! 2: exp(i(k*z - omega*t))
        integer :: convention = 1
    end type

    interface assignment (=)
        module procedure lib_field_gaussian_beam_eq_plane_wave_parameter
    end interface

    contains

        subroutine lib_field_gaussian_beam_eq_plane_wave_parameter(lhs, rhs)
            implicit none
            ! dummy
            type(lib_field_plane_wave_type), intent(inout) :: lhs
            type(lib_field_gaussian_beam_hermite_type), intent(in) :: rhs


            lhs%e_field_0 = rhs%e_field_0
            lhs%wave_length_0 = rhs%wave_length_0
            lhs%refractive_index_medium = rhs%refractive_index_medium
            lhs%phase = rhs%phase
            lhs%polarisation = rhs%polarisation
            lhs%theta = rhs%theta
            lhs%phi = rhs%phi
            lhs%convention = rhs%convention

        end subroutine lib_field_gaussian_beam_eq_plane_wave_parameter

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermsubroutnite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !
        !                  z   k
        !                   ^ /
        !                   |/       o
        !               K_B --->
        !                       x
        !
        !          z
        !           ^
        !           |
        !       K_0 --->
        !               x
        !
        !   K_0: world coordinate system
        !   K_B: beam coordinate system
        !   o: evaluation_point_x (K_B)
        !
        ! Returns
        ! ----
        !   e_field: type(cartesian_coordinate_cmplx_type)
        !       electical field component at "evaluation_point_x"
        !   h_field: type(cartesian_coordinate_cmplx_type)
        !       magnetical field component at "evaluation_point_x"
        !
        subroutine lib_field_gaussian_beam_hermite_get_field(parameter, evaluation_point_x, &
                                                             e_field, h_field)
            implicit none
            ! dummy
            type(lib_field_gaussian_beam_hermite_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            type(cartesian_coordinate_cmplx_type) :: e_field
            type(cartesian_coordinate_cmplx_type) :: h_field

            ! auxiliary
            double complex :: field
            double precision :: wave_impedance

            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: point_x

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(PI / 2d0 + parameter%phi, parameter%theta, PI / 2d0 + parameter%phi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            field = lib_fiel_gaussian_beam_hermite_scalar(parameter%e_field_0, &
                                                          parameter%wave_length_0, &
                                                          parameter%refractive_index_medium, &
                                                          parameter%waist_x0, parameter%waist_y0, &
                                                          point_x, &
                                                          phase=parameter%phase, &
                                                          tem_m=parameter%tem_m, &
                                                          tem_n=parameter%tem_n, &
                                                          convention=parameter%convention)
            e_field = field * parameter%polarisation;

            ! Quantum Optics: An Introduction, Mark Fox
            ! eq. 2.25 but plane wave
            wave_impedance = const_z_0 / parameter%refractive_index_medium

            rot  = lib_math_get_matrix_rot(-PI / 2d0, 0d0, 0d0)
            h_field = rot * e_field / wave_impedance

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(-PI / 2d0 - parameter%phi, &
                                               -parameter%theta, &
                                               -PI / 2d0 - parameter%phi)
                e_field = rot * e_field
                h_field = rot * h_field
            end if

        end subroutine lib_field_gaussian_beam_hermite_get_field

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !
        !                  z   k
        !                   ^ /
        !                   |/
        !               K_B --->
        !                       x
        !
        !          z
        !           ^
        !           |
        !       K_0 --->
        !               x
        !
        !   K_0: world coordinate system
        !   K_B: beam coordinate system
        !   k: k_vector_n (K_B)
        !
        ! Returns
        ! ----
        !   k_vector_n: type(cartesian_coordinate_real_type)
        !       normalised k vector at x
        !
        function lib_field_gaussian_beam_hermite_get_propagation_direction(parameter, evaluation_point_x) result (k_vector_n)
            implicit none

            ! dummy
            type(lib_field_gaussian_beam_hermite_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            type(cartesian_coordinate_real_type) :: k_vector_n

            ! auxiliary
            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: point_x

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(PI / 2d0 + parameter%phi, parameter%theta, PI / 2d0 + parameter%phi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            k_vector_n = lib_fiel_gaussian_beam_hermite_scalar_propagation_direction(parameter%e_field_0, &
                                                                           parameter%wave_length_0, &
                                                                           parameter%refractive_index_medium, &
                                                                           parameter%waist_x0, parameter%waist_y0, &
                                                                           point_x, &
                                                                           tem_m=parameter%tem_m, &
                                                                           tem_n=parameter%tem_n, &
                                                                           convention=parameter%convention)

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(-PI / 2d0 - parameter%phi, -parameter%theta, -PI / 2d0 - parameter%phi)
                k_vector_n = rot * k_vector_n
            end if

        end function lib_field_gaussian_beam_hermite_get_propagation_direction

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !
        !
        !                  z   k     \
        !                   ^ /     \ \ plane wave approximation
        !                   |/       o \
        !               K_B --->      \ \
        !                       x
        !
        !          z
        !           ^
        !           |
        !       K_0 --->
        !               x
        !
        !   K_0: world coordinate system
        !   K_B: beam coordinate system
        !   o: evaluation_point_x (at K_B)
        !
        ! Returns
        ! ----
        !   plane_wave: type(lib_field_plane_wave_type)
        !       parameter of the plane wave of Gaussian beam at the point "evaluation_point_x"
        !
        function lib_field_gaussian_beam_hermite_get_plane_wave_approximation(parameter, evaluation_point_x) result(plane_wave)
            implicit none
            ! dummy
            type(lib_field_gaussian_beam_hermite_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            type(lib_field_plane_wave_type) :: plane_wave

            ! auxiliary
            type(cartesian_coordinate_real_type) :: k_vector_n
            type(spherical_coordinate_real_type) :: k_vector_n_spherical

            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: x_axis_pw
            type(cartesian_coordinate_real_type) :: y_axis_pw
            type(cartesian_coordinate_real_type) :: point_x

            double precision :: absolute
            double precision :: argument

            double precision :: buffer


            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(PI / 2d0 + parameter%phi, parameter%theta, PI / 2d0 + parameter%phi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            plane_wave = parameter

            ! get propagation direction
            k_vector_n = lib_field_gaussian_beam_hermite_get_propagation_direction(parameter, evaluation_point_x)

            k_vector_n_spherical = k_vector_n

            plane_wave%theta = k_vector_n_spherical%theta
            plane_wave%phi = k_vector_n_spherical%phi

            ! get absolut value und phase of the gaussian beam at the evaluation point
            call get_absolute_value_and_argument(parameter%e_field_0, &
                                                 parameter%wave_length_0, &
                                                 parameter%refractive_index_medium, &
                                                 parameter%waist_x0, &
                                                 parameter%waist_y0, &
                                                 point_x, &
                                                 absolute, argument, &
                                                 phase = parameter%phase, &
                                                 tem_m = parameter%tem_m, &
                                                 tem_n = parameter%tem_n, &
                                                 convention = parameter%convention)

            ! plane wave parameter at the origin
            select case(parameter%convention)
                case (1)
                    ! exp(-i(k*z - omega*t))
                    plane_wave%phase = argument + mod(2d0 * PI * point_x%z / parameter%wave_length_0, 2d0 * PI)

                case (2)
                    ! exp(i(k*z - omega*t))
!                    plane_wave%phase = argument - mod(2d0 * PI * point_x%z / parameter%wave_length_0, 2d0 * PI)
                    plane_wave%phase = argument + mod(2d0 * PI * point_x%z / parameter%wave_length_0, 2d0 * PI)

                case default
                    print *, "lib_field_gaussian_beam_hermite_get_plane_wave_approximation: ERROR"
                    print *, "  convention is not defined: ", parameter%convention
            end select

            ! adapt range: [0, 2 Pi)
            plane_wave%phase = mod(plane_wave%phase, 2d0 * PI)
            if (plane_wave%phase .lt. 0d0) then
                plane_wave%phase = 2d0 * PI + plane_wave%phase
            end if

            plane_wave%e_field_0 = absolute

            ! get the polarisation of the plane wave at the evaluation point
            x_axis_pw = make_cartesian(1d0, 0d0, 0d0)
            y_axis_pw = make_cartesian(0d0, 1d0, 0d0)

            if (plane_wave%phi .ne. 0 .or. plane_wave%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(PI / 2d0 + parameter%phi, &
                                               parameter%theta, &
                                               PI / 2d0 + parameter%phi)
                x_axis_pw = rot * x_axis_pw
                y_axis_pw = rot * y_axis_pw
            end if

            plane_wave%polarisation%x = parameter%polarisation%x * abs(x_axis_pw%x) &
                                        + parameter%polarisation%y * abs(y_axis_pw%x)

            plane_wave%polarisation%y = parameter%polarisation%y * abs(y_axis_pw%y) &
                                        + parameter%polarisation%x * abs(x_axis_pw%y)

            buffer = abs(abs(make_cartesian(plane_wave%polarisation%x, &
                                            plane_wave%polarisation%y, &
                                            dcmplx(0,0))))

            plane_wave%polarisation%x = plane_wave%polarisation%x / buffer
            plane_wave%polarisation%y = plane_wave%polarisation%y / buffer

        end function lib_field_gaussian_beam_hermite_get_plane_wave_approximation

        ! Argument
        ! ----
        !   e_field_0: double precision
        !       electical field [V/m]
        !   wave_length: double precision
        !       vacuum wave length [m]
        !   n_medium: double precision
        !       refractive index of the medium
        !   waist_x0: double precision
        !       beam waist at z = 0 in along the x axis
        !   waist_y0: double precision
        !       beam waist at z = 0 in along the y axis
        !   x: type(cartesian_coordinate_real)
        !       evaluation point [m^3]
        !   tem_m, tem_n: integer, optional(std: 0)
        !       order of the Hermite-Gaussian mode
        !   convention: integer, optional (std: 1)
        !       1: exp(-i(k*z - omega*t))
        !       2: exp(i(k*z - omega*t))
        !
        ! Returns
        ! ----
        !   field: double cmplx
        !       Gaussian field distribution with a propagation along the z-axis
        !
        !
        !          \     z       / <-- Gaussian beam
        !           \     ^     /
        !            \ K_B|    /
        !          -->|   --> |<-- 2 w_0x,y
        !            /     x,y \
        !           /           \
        !
        !       K_B: beam coordinate system
        !       z:   propagation along the z-axis
        !
        !
        ! TeX: $$ \begin{aligned} E_{m n}(x, y, z, t)=& \sqrt{\frac{w_{x, 0}}{w_{x}(z)}} H_{m}\left(\frac{x \sqrt{2}}{w_{x}(z)}\right) e^{-\frac{x^{2}}{w_{x}^{2}(z)}} e^{-i \frac{k}{2 R_{x}(z)}} \\ &\times \sqrt{\frac{w_{y, 0}}{w_{y}(z)} H_{n}\left(\frac{y \sqrt{2}}{w_{y}(z)}\right.}{w_{y}(z)}) e^{\frac{y^{2}}{w_{y}^{2}(z)}}{e^{w_{y}^{2}(z)}}-i \frac{k}{2 R_{y}(z)} \\ & \times E_{0} \sqrt{\frac{1}{2^{m+n} m ! n ! \pi}} e^{-i(k z-\omega t)} \times \\ & \times e^{i(m+1 / 2) \arctan \left(z / z_{R, x}\right)} e^{i(n+1 / 2) \arctan \left(z / z_{R, y}\right)} \end{aligned}$$
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.28 without the time-dependent term
        function lib_fiel_gaussian_beam_hermite_scalar(e_field_0, wave_length_0, n_medium, &
                                                  waist_x0, waist_y0, &
                                                  x, phase, tem_m, tem_n, &
                                                  convention) result(field)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_x0
            double precision, intent(in) :: waist_y0
            type(cartesian_coordinate_real_type), intent(in) :: x
            double precision, intent(in), optional :: phase
            integer, intent(in), optional :: tem_m
            integer, intent(in), optional :: tem_n
            integer, intent(in), optional :: convention

            double complex :: field

            ! auxiliaray
            double precision :: m_phase
            integer :: m_tem_m
            integer :: m_tem_n
            integer :: m_convention

            double precision :: w_x
            double precision :: zr_x
            double precision :: cr_x

            double precision :: w_y
            double precision :: zr_y
            double precision :: cr_y

            double precision :: k

            integer :: h_min_i
            integer :: h_max_i
            double precision, dimension(:,:), allocatable :: h

            double precision :: buffer_real

            double complex :: buffer_cmplx

            m_phase = 0
            if (present(phase)) m_phase = phase

            m_tem_m = 0
            if (present(tem_m)) m_tem_m = tem_m

            m_tem_n = 0
            if (present(tem_n)) m_tem_n = tem_n

            m_convention = 1
            if (present(convention)) m_convention = convention


            ! pre calc
            k = 2d0 * PI * n_medium / wave_length_0

            zr_x = lib_field_gaussian_beam_rayleigh_length(waist_x0, wave_length_0, n_medium)
            w_x = lib_field_gaussian_beam_waist(waist_x0, zr_x, x%z)
            cr_x = lib_field_gaussian_beam_curvature(zr_x, x%z)

            zr_y = lib_field_gaussian_beam_rayleigh_length(waist_y0, wave_length_0, n_medium)
            w_y = lib_field_gaussian_beam_waist(waist_y0, zr_y, x%z)
            cr_y = lib_field_gaussian_beam_curvature(zr_y, x%z)


            h_min_i = min(m_tem_m, m_tem_n)
            h_max_i = max(m_tem_m, m_tem_n)

            call lib_math_hermite_polynomial((/ x%x * sqrt(2d0) / w_x, &
                                               x%y * sqrt(2d0) / w_y /), &
                                            h_min_i, h_max_i - h_min_i + 1, h)

            buffer_real = sqrt(waist_x0 / w_x) * h(1, m_tem_m) * exp(- x%x**2d0 / w_x**2d0)
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%x**2 * cr_x / 2d0))

            buffer_real = sqrt(waist_y0 / w_y) * h(2, m_tem_n) * exp(- x%y**2d0 / w_y**2d0)
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%y**2 * cr_y / 2d0)) * buffer_cmplx

            buffer_real = e_field_0 * sqrt(1d0 &
                                           / (2d0**dble(m_tem_m + m_tem_n) &
                                              * lib_math_factorial_get_factorial(m_tem_m) &
                                              * lib_math_factorial_get_factorial(m_tem_n) * PI))

            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%z + phase)) * buffer_cmplx

            buffer_real = (dble(m_tem_m) + 0.5d0) * atan2(x%z, zr_x) &
                          + (dble(m_tem_n) + 0.5d0) * atan2(x%z, zr_y)
            buffer_cmplx = exp(dcmplx(0, buffer_real)) * buffer_cmplx

            select case(m_convention)
                case(1)
                    ! exp(-i(k*z - omega*t))
                    field = buffer_cmplx
                case(2)
                    ! exp(i(k*z - omega*t))
                    field = conjg(buffer_cmplx)

                case default
                    print *, "lib_fiel_gaussian_beam_hermite_scalar: ERROR"
                    print *, "  convention is not defined: ", m_convention
            end select


        end function lib_fiel_gaussian_beam_hermite_scalar

        ! Argument
        ! ----
        !   e_field_0: double precision
        !       electical field [V/m]
        !   wave_length: double precision
        !       vacuum wave length [m]
        !   n_medium: double precision
        !       refractive index of the medium
        !   waist_x0: double precision
        !       beam waist at z = 0 in along the x axis
        !   waist_y0: double precision
        !       beam waist at z = 0 in along the y axis
        !   x: type(cartesian_coordinate_real)
        !       evaluation point [m^3]
        !   tem_m, tem_n: integer, optional(std: 0)
        !       order of the Hermite-Gaussian mode
        !   convention: integer, optional (std: 1)
        !       1: exp(-i(k*z - omega*t))
        !       2: exp(i(k*z - omega*t))
        !
        ! Returns
        ! ----
        !   k_vector_n: type(cartesian_coordinate_real_type)
        !       normalised k vector at x
        !
        !
        !          \     z       / <-- Gaussian beam
        !           \     ^     /
        !            \ K_B|    /
        !          -->|   --> |<-- 2 w_0x,y
        !            /     x,y \
        !           /           \
        !
        !       K_B: beam coordinate system
        !       z:   propagation along the z-axis
        !
        ! TeX: $$ \begin{aligned} E_{m n}(x, y, z, t)=& \sqrt{\frac{w_{x, 0}}{w_{x}(z)}} H_{m}\left(\frac{x \sqrt{2}}{w_{x}(z)}\right) e^{-\frac{x^{2}}{w_{x}^{2}(z)}} e^{-i \frac{k}{2 R_{x}(z)}} \\ &\times \sqrt{\frac{w_{y, 0}}{w_{y}(z)} H_{n}\left(\frac{y \sqrt{2}}{w_{y}(z)}\right.}{w_{y}(z)}) e^{\frac{y^{2}}{w_{y}^{2}(z)}}{e^{w_{y}^{2}(z)}}-i \frac{k}{2 R_{y}(z)} \\ & \times E_{0} \sqrt{\frac{1}{2^{m+n} m ! n ! \pi}} e^{-i(k z-\omega t)} \times \\ & \times e^{i(m+1 / 2) \arctan \left(z / z_{R, x}\right)} e^{i(n+1 / 2) \arctan \left(z / z_{R, y}\right)} \end{aligned}$$
        function lib_fiel_gaussian_beam_hermite_scalar_propagation_direction(e_field_0, wave_length_0, n_medium, &
                                                                             waist_x0, waist_y0, &
                                                                             x, tem_m, tem_n, &
                                                                             convention) result(k_vector_n)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_x0
            double precision, intent(in) :: waist_y0
            type(cartesian_coordinate_real_type), intent(in) :: x
            integer, intent(in), optional :: tem_m
            integer, intent(in), optional :: tem_n
            integer, intent(in), optional :: convention

            type(cartesian_coordinate_real_type) :: k_vector_n

            k_vector_n%x = -get_derivative_xy(wave_length_0, n_medium, waist_x0, x%x, x%z)
            k_vector_n%y = -get_derivative_xy(wave_length_0, n_medium, waist_y0, x%y, x%z)
            k_vector_n%z = -get_derivative_z(wave_length_0, n_medium, waist_x0, waist_y0, tem_m, tem_n, x%x, x%y, x%z, &
                                             convention)

            k_vector_n = k_vector_n / abs(k_vector_n)

        end function

        function get_derivative_xy(wave_length_0, n_medium, waist_xy0, xy, z) result(rv)
            !dummy
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_xy0
            double precision, intent(in) :: xy
            double precision, intent(in) :: z

            double precision :: rv

            ! auxiliary
            double precision :: numerator
            double precision :: denominator

            if (xy .eq. 0 .or. z .eq. 0) then
                rv = 0d0
            else
                numerator = - 2d0 * PI * wave_length_0 * xy * z * n_medium
                denominator = (PI * n_medium)**2d0 * waist_xy0**4d0 + (wave_length_0 * z)**2d0

                rv = numerator / denominator
            end if
        end function get_derivative_xy

        function get_derivative_z(wave_length_0, n_medium, waist_x0, waist_y0, m, n, x, y, z, convention) result(rv)
            !dummy
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_x0
            double precision, intent(in) :: waist_y0
            integer, intent(in) :: m
            integer, intent(in) :: n
            double precision, intent(in) :: x
            double precision, intent(in) :: y
            double precision, intent(in) :: z
            integer, intent(in), optional :: convention

            double precision :: rv

            ! auxiliary
            double precision :: numerator
            double precision :: denominator

            double precision :: pi_n_medium_2
            double precision :: lambda_z_2

            double precision :: waist_x0_2
            double precision :: waist_x0_4

            double precision :: waist_y0_2
            double precision :: waist_y0_4

            integer :: m_convention

            m_convention = 1
            if (present(convention)) m_convention = convention

            ! pre-calc
            pi_n_medium_2 = (PI * n_medium)**2d0
            lambda_z_2 = (wave_length_0 * z)**2d0

            waist_x0_2 = waist_x0**2d0
            waist_x0_4 = waist_x0_2**2d0

            waist_y0_2 = waist_y0**2d0
            waist_y0_4 = waist_y0_2**2d0

            ! calc
!            numerator = (2d0 * m + 1d0) * waist_x0_2 + 2d0 * x**2d0
            numerator = 2d0 * m * waist_x0_2 + 2d0 * x**2d0
            denominator = pi_n_medium_2 * waist_x0_4 + lambda_z_2

            rv = numerator / denominator

!            numerator = (2d0 * n + 1d0) * waist_y0_2 + 2d0 * y**2d0
            numerator = 2d0 * n * waist_y0_2 + 2d0 * y**2d0
            denominator = pi_n_medium_2 * waist_y0_4 + lambda_z_2

            rv = rv + numerator / denominator

            numerator = - 4d0 * pi_n_medium_2 * x**2d0 * waist_x0_4
            denominator = (pi_n_medium_2 * waist_x0_4 + lambda_z_2)**2d0

            rv = rv + numerator / denominator

            numerator = - 4d0 * pi_n_medium_2 * y**2d0 * waist_y0_4
            denominator = (pi_n_medium_2 * waist_y0_4 + lambda_z_2)**2d0

            rv = rv + numerator / denominator

            rv = rv * wave_length_0**2d0 / 2d0 - 2d0

            rv = PI * n_medium * rv / wave_length_0

!            select case (m_convention)
!                case (1)
!                    ! exp(-i(k*z - omega*t))
!                case (2)
!                    ! exp(i(k*z - omega*t))
!                    rv = - PI * n_medium * rv / wave_length_0
!
!                case default
!                    print *, "lib_field_gaussian_beam__get_derivative_z: ERROR"
!                    print *, "  convention is not defined: ", m_convention
!            end select

        end function get_derivative_z

        subroutine get_absolute_value_and_argument(e_field_0, wave_length_0, n_medium, waist_x0, waist_y0, x, &
                                                   absolute, argument, &
                                                   phase, tem_m, tem_n, convention)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_x0
            double precision, intent(in) :: waist_y0
            type(cartesian_coordinate_real_type), intent(in) :: x

            double precision, intent(in), optional :: phase
            integer, intent(in), optional :: tem_m
            integer, intent(in), optional :: tem_n
            integer, intent(in), optional :: convention

            double precision :: absolute
            double precision :: argument

            ! auxiliary
            double precision :: zrx
            double precision :: zry

            double precision :: wx
            double precision :: wy

            double precision :: crx
            double precision :: cry

            double precision :: k
            double precision :: kz

            integer :: h_min_i
            integer :: h_max_i
            double precision, dimension(:,:), allocatable :: h

            double precision :: m_phase
            integer :: m_tem_m
            integer :: m_tem_n
            integer :: m_convention

            m_phase = 0
            if (present(phase)) m_phase = phase

            m_tem_m = 0
            if (present(tem_m)) m_tem_m = tem_m

            m_tem_n = 0
            if (present(tem_n)) m_tem_n = tem_n

            m_convention = 1
            if (present(convention)) m_convention = convention

            k = 2d0 * PI * n_medium / wave_length_0

            kz = mod(k * x%z, 2d0 * PI)

            zrx = lib_field_gaussian_beam_rayleigh_length(waist_x0, wave_length_0, n_medium)
            zry = lib_field_gaussian_beam_rayleigh_length(waist_y0, wave_length_0, n_medium)

            wx = lib_field_gaussian_beam_waist(waist_x0, zrx, x%z)
            wy = lib_field_gaussian_beam_waist(waist_y0, zry, x%z)

            crx = lib_field_gaussian_beam_curvature(zrx, x%z)
            cry = lib_field_gaussian_beam_curvature(zry, x%z)


            ! Calculation of the argument of the complex electrictical field
            argument = -mod(k * crx * x%x*x%x / 2d0, 2d0 * PI)
            argument = argument -mod(k * cry * x%y*x%y / 2d0, 2d0 * PI)
            argument = argument - kz + m_phase
            argument = argument + (dble(m_tem_m) + 0.5d0) * atan2(x%z, zrx)
            argument = argument + (dble(m_tem_n) + 0.5d0) * atan2(x%z, zry)

            h_min_i = min(m_tem_m, m_tem_n)
            h_max_i = max(m_tem_m, m_tem_n)

            call lib_math_hermite_polynomial((/ x%x * sqrt(2d0) / wx, &
                                               x%y * sqrt(2d0) / wy /), &
                                            h_min_i, h_max_i - h_min_i + 1, h)

            if (h(1, m_tem_m) * h(2, m_tem_n) .lt. 0d0) then
                argument = argument - PI
            end if

            argument = mod(argument, 2d0 * PI)


            select case(m_convention)
                case(1)
                    ! exp(-i(k*z - omega*t))
                    argument = argument

                case(2)
                    ! exp(i(k*z - omega*t))
                    argument = -argument

                case default
                    print *, "lib_field_gaussian_beam__get_phase: ERROR"
                    print *, "  convention is not defined: ", m_convention
            end select

            if (argument .lt. 0d0) then
                argument = 2d0 * PI + argument
            end if

            ! Calculation of the absolute value of the complex electrical field
            absolute = sqrt(waist_x0 / wx) * h(1, m_tem_m) * exp(-(x%x / wx)**2d0) &
                       * sqrt(waist_y0 / wy) * h(2, m_tem_n) * exp(-(x%y / wy)**2d0) &
                       * e_field_0 &
                       * sqrt(1d0 / (dble(2**(m_tem_m + m_tem_n)) &
                                    * lib_math_factorial_get_factorial(m_tem_m) &
                                    * lib_math_factorial_get_factorial(m_tem_n) &
                                     * PI))

            absolute = abs(absolute)

        end subroutine get_absolute_value_and_argument

        ! Argument
        ! ----
        !   w0: double precision
        !       beam waist at z = 0
        !   zr: double precision
        !       Rayleigh length
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   w: double precision
        !       beam waist at z
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.29
        function lib_field_gaussian_beam_waist(w0, zr, z) result(w)
            implicit none
            ! dummy
            double precision, intent(in) :: w0
            double precision, intent(in) :: zr
            double precision, intent(in) :: z

            double precision :: w

            w = w0 * sqrt(1 + (z/zr)**2)

        end function lib_field_gaussian_beam_waist

        ! Argument
        ! ----
        !   zr: double precision
        !       Rayleigh length
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   c: double precision
        !       radius of curvature of the beam's wavefront
        !
        ! Formula
        ! ----
        !   c = R^-1
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            inverse of eq. 2.35
        function lib_field_gaussian_beam_curvature(zr, z) result(c)
            implicit none
            ! dummy
            double precision, intent(in) :: zr
            double precision, intent(in) :: z

            double precision :: c

            ! R[zr_, z_] := zr(z/zr + zr/z);
            ! cr[zr_, z_] := R[zr, z]^-1;
            c = z / (z**2 + zr**2)

        end function lib_field_gaussian_beam_curvature

        ! Argument
        ! ----
        !   w0: double precision
        !       beam waist at z = 0
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   zr: double precision
        !       Rayleigh length
        !
        ! Formula
        ! ----
        !   c = R^-1
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.30
        function lib_field_gaussian_beam_rayleigh_length(w0, wave_length_0, n_medium) result(zr)
            implicit none
            ! dummy
            double precision, intent(in) :: w0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium

            double precision :: zr

            zr = PI * w0**2 * n_medium / wave_length_0

        end function lib_field_gaussian_beam_rayleigh_length

        function lib_field_gaussian_beam_test_functions() result(rv)
            use lib_field
            implicit none
            ! dummy
            integer :: rv

            ! parameter
            double precision, parameter :: ground_truth_e = 1d-7

            rv = 0

            if (.not. test_get_absolute_value_and_argument_convention_1()) rv = rv + 1
            if (.not. test_get_absolute_value_and_argument_convention_2()) rv = rv + 1
            if (.not. test_lib_field_gaussian_beam_hermite_get_plane_wave_approx()) rv = rv + 1
            if (.not. test_lib_field_gaussian_beam_hermite_get_field_approximation()) rv = rv + 1
            if (.not. test_lib_field_gaussian_beam_hermite_get_field_approx_full()) rv = rv + 1
            if (.not. test_lib_field_gaussian_beam_hermite_get_propagation_direction()) rv = rv + 1
            if (.not. test_lib_field_gaussian_beam_hermite_get_field()) rv = rv + 1
!            if (.not. test_lib_field_gaussian_beam_hermite_get_field_2()) rv = rv + 1


            if (rv == 0) then
                print *, "lib_field_gaussian_beam_test_functions tests: OK"
            else
                print *, rv,"lib_field_gaussian_beam_test_functions test(s) FAILED"
            end if

            contains

            function test_get_absolute_value_and_argument_convention_1() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                double complex :: field_scalar
                double precision :: field_absolute
                double precision :: field_argument

                double precision :: buffer

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 1.0_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                x = 0
                y = 1 * unit_mu
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 1
                gauss_parameter%tem_n = 1
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 2 * unit_mu
                gauss_parameter%waist_y0 = 2 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = 0!PI / 8d0
                gauss_parameter%phi = 0!PI

                gauss_parameter%convention = 1

                rv = .true.
                print *, "test_get_absolute_value_and_argument_convention_1"
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        field_scalar = lib_fiel_gaussian_beam_hermite_scalar(gauss_parameter%e_field_0, &
                                                                             gauss_parameter%wave_length_0, &
                                                                             gauss_parameter%refractive_index_medium, &
                                                                             gauss_parameter%waist_x0, &
                                                                             gauss_parameter%waist_y0, &
                                                                             point_cartesian, &
                                                                             gauss_parameter%phase, &
                                                                             gauss_parameter%tem_m, &
                                                                             gauss_parameter%tem_n, &
                                                                             gauss_parameter%convention)

                        call get_absolute_value_and_argument(gauss_parameter%e_field_0, &
                                                             gauss_parameter%wave_length_0, &
                                                             gauss_parameter%refractive_index_medium, &
                                                             gauss_parameter%waist_x0, &
                                                             gauss_parameter%waist_y0, &
                                                             point_cartesian, &
                                                             field_absolute, &
                                                             field_argument, &
                                                             gauss_parameter%phase, &
                                                             gauss_parameter%tem_m, &
                                                             gauss_parameter%tem_n, &
                                                             gauss_parameter%convention)

                        print *, "  (", i, ",", ii, ")"
!                        buffer = abs(field_scalar) - field_absolute
!                        if (abs(buffer) .lt. ground_truth_e) then
!                            print *, "  absolute: OK"
!                        else
!                            print *, "  absolute: FAILED"
!                            print *, "              diff = ", buffer
!                            print *, "       abs(scalar) = ", abs(field_scalar)
!                            print *, "    field_absolute = ", field_absolute
!                            rv = .false.
!                        end if
!
!                        buffer = atan2(aimag(field_scalar), real(field_scalar))
!                        if (buffer .lt. 0d0) buffer = 2d0 * PI + buffer
!                        buffer = buffer - field_argument
!                        if (abs(buffer) .lt. ground_truth_e) then
!                            print *, "  argument: OK"
!                        else
!                            print *, "  argument: FAILED"
!                            print *, "              diff = ", buffer
!                            print *, "       arg(scalar) = ", atan2(aimag(field_scalar), real(field_scalar))
!                            print *, "    field_argument = ", field_argument
!                            rv = .false.
!                        end if

                        buffer = real(field_scalar) - real(field_absolute * exp(dcmplx(0, field_argument)))
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  real: OK"
                        else
                            print *, "  real: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       abs(scalar) = ", real(field_scalar)
                            print *, "    field_absolute = ", real(field_absolute * exp(dcmplx(0, field_argument)))
                            rv = .false.
                        end if

                        buffer = aimag(field_scalar) - aimag(field_absolute * exp(dcmplx(0, field_argument)))
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  cmplx: OK"
                        else
                            print *, "  cmplx: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       abs(scalar) = ", aimag(field_scalar)
                            print *, "    field_absolute = ", aimag(field_absolute * exp(dcmplx(0, field_argument)))
                            rv = .false.
                        end if

                    end do
                end do

            end function test_get_absolute_value_and_argument_convention_1

            function test_get_absolute_value_and_argument_convention_2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                double complex :: field_scalar
                double precision :: field_absolute
                double precision :: field_argument

                double precision :: buffer

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 1.0_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 2 * unit_mu
                gauss_parameter%waist_y0 = 2 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = 0!PI / 8d0
                gauss_parameter%phi = 0!PI

                gauss_parameter%convention = 2

                rv = .true.
                print *, "test_get_absolute_value_and_argument_convention_2"
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        field_scalar = lib_fiel_gaussian_beam_hermite_scalar(gauss_parameter%e_field_0, &
                                                                             gauss_parameter%wave_length_0, &
                                                                             gauss_parameter%refractive_index_medium, &
                                                                             gauss_parameter%waist_x0, &
                                                                             gauss_parameter%waist_y0, &
                                                                             point_cartesian, &
                                                                             gauss_parameter%phase, &
                                                                             gauss_parameter%tem_m, &
                                                                             gauss_parameter%tem_n, &
                                                                             gauss_parameter%convention)

                        call get_absolute_value_and_argument(gauss_parameter%e_field_0, &
                                                             gauss_parameter%wave_length_0, &
                                                             gauss_parameter%refractive_index_medium, &
                                                             gauss_parameter%waist_x0, &
                                                             gauss_parameter%waist_y0, &
                                                             point_cartesian, &
                                                             field_absolute, &
                                                             field_argument, &
                                                             gauss_parameter%phase, &
                                                             gauss_parameter%tem_m, &
                                                             gauss_parameter%tem_n, &
                                                             gauss_parameter%convention)

                        print *, "  (", i, ",", ii, ")"
                        buffer = abs(field_scalar) - field_absolute
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  absolute: OK"
                        else
                            print *, "  absolute: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       abs(scalar) = ", abs(field_scalar)
                            print *, "    field_absolute = ", field_absolute
                            rv = .false.
                        end if

                        buffer = atan2(aimag(field_scalar), real(field_scalar))
                        if (buffer .lt. 0d0) buffer = 2d0 * PI + buffer
                        buffer = buffer - field_argument
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  argument: OK"
                        else
                            print *, "  argument: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       arg(scalar) = ", atan2(aimag(field_scalar), real(field_scalar))
                            print *, "    field_argument = ", field_argument
                            rv = .false.
                        end if

                        buffer = real(field_scalar) - real(field_absolute * exp(dcmplx(0, field_argument)))
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  real: OK"
                        else
                            print *, "  real: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       abs(scalar) = ", real(field_scalar)
                            print *, "    field_absolute = ", real(field_absolute * exp(dcmplx(0, field_argument)))
                            rv = .false.
                        end if

                        buffer = aimag(field_scalar) - aimag(field_absolute * exp(dcmplx(0, field_argument)))
                        if (abs(buffer) .lt. ground_truth_e) then
                            print *, "  cmplx: OK"
                        else
                            print *, "  cmplx: FAILED"
                            print *, "              diff = ", buffer
                            print *, "       abs(scalar) = ", aimag(field_scalar)
                            print *, "    field_absolute = ", aimag(field_absolute * exp(dcmplx(0, field_argument)))
                            rv = .false.
                        end if

                    end do
                end do

            end function test_get_absolute_value_and_argument_convention_2

            function test_lib_field_gaussian_beam_hermite_get_field() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.05_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 1.5 * unit_mu
                gauss_parameter%waist_y0 = 1.5 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = PI / 8d0
                gauss_parameter%phi = 0!PI

                gauss_parameter%convention = 1

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(1) + (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter, &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field)

                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_")

            end function test_lib_field_gaussian_beam_hermite_get_field

            function test_lib_field_gaussian_beam_hermite_get_field_2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type), dimension(2) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.075_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter(:)%e_field_0 = 1
                gauss_parameter(:)%refractive_index_medium = 1
                gauss_parameter(:)%tem_m = 0
                gauss_parameter(:)%tem_n = 0
!                gauss_paramete(:)r%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter(:)%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_paramete(:)r%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter(:)%waist_x0 = 2.5 * unit_mu
                gauss_parameter(:)%waist_y0 = 2.5 * unit_mu
                gauss_parameter(:)%wave_length_0 = 0.55 * unit_mu

                gauss_parameter(1)%theta = PI / 16d0
                gauss_parameter(1)%phi = 0

                gauss_parameter(2)%theta = PI / 16d0
                gauss_parameter(2)%phi = PI

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_y_values
                        y = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter(1), &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field)
                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter(2), &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field)
                        e_field(i,ii) = e_field(i,ii) + buffer_e_field
                        h_field(i,ii) = h_field(i,ii) + buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_2_")

            end function

            function test_lib_field_gaussian_beam_hermite_get_propagation_direction() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                double precision, parameter :: ground_truth_e = 1D-5

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_real_type) :: buffer_k_vector_n
!                type(cartesian_coordinate_real_type), dimension(:,:), allocatable :: k_vector_n
                type(cartesian_coordinate_real_type), dimension(:,:), allocatable :: ground_truth_k_vector_n

                double precision :: buffer

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 5d0 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

!                allocate(k_vector_n(no_x_values, no_y_values))
                allocate(ground_truth_k_vector_n(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 2.5 * unit_mu
                gauss_parameter%waist_y0 = 2.5 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = 0d0
                gauss_parameter%phi = 0d0

                ground_truth_k_vector_n(1,1)%x = -0.068322
                ground_truth_k_vector_n(1,1)%y = 0.068322
                ground_truth_k_vector_n(1,1)%z = 0.995321

                ground_truth_k_vector_n(1,2)%x = -0.0699719
                ground_truth_k_vector_n(1,2)%y = 0.034986
                ground_truth_k_vector_n(1,2)%z = 0.996935

                ground_truth_k_vector_n(1,3)%x = -0.070541
                ground_truth_k_vector_n(1,3)%y = 0.
                ground_truth_k_vector_n(1,3)%z = 0.997509

                ground_truth_k_vector_n(1,4)%x = -0.0699719
                ground_truth_k_vector_n(1,4)%y = -0.034986
                ground_truth_k_vector_n(1,4)%z = 0.996935

                ground_truth_k_vector_n(2,1)%x = -0.034986
                ground_truth_k_vector_n(2,1)%y = 0.0699719
                ground_truth_k_vector_n(2,1)%z = 0.996935

                ground_truth_k_vector_n(2,2)%x = -0.0358548
                ground_truth_k_vector_n(2,2)%y = 0.0358548
                ground_truth_k_vector_n(2,2)%z = 0.998714

                ground_truth_k_vector_n(2,3)%x = -0.0361547
                ground_truth_k_vector_n(2,3)%y = 0.
                ground_truth_k_vector_n(2,3)%z = 0.999346

                ground_truth_k_vector_n(2,4)%x = -0.0358548
                ground_truth_k_vector_n(2,4)%y = -0.0358548
                ground_truth_k_vector_n(2,4)%z = 0.998714

                ground_truth_k_vector_n(3,1)%x = 0.
                ground_truth_k_vector_n(3,1)%y = 0.070541
                ground_truth_k_vector_n(3,1)%z = 0.997509

                ground_truth_k_vector_n(3,2)%x = 0.
                ground_truth_k_vector_n(3,2)%y = 0.0361547
                ground_truth_k_vector_n(3,2)%z = 0.999346

                ground_truth_k_vector_n(3,3)%x = 0.
                ground_truth_k_vector_n(3,3)%y = 0.
                ground_truth_k_vector_n(3,3)%z = 1.

                ground_truth_k_vector_n(3,4)%x = 0.
                ground_truth_k_vector_n(3,4)%y = -0.0361547
                ground_truth_k_vector_n(3,4)%z = 0.999346

                ground_truth_k_vector_n(4,1)%x = 0.034986
                ground_truth_k_vector_n(4,1)%y = 0.0699719
                ground_truth_k_vector_n(4,1)%z = 0.996935

                ground_truth_k_vector_n(4,2)%x = 0.0358548
                ground_truth_k_vector_n(4,2)%y = 0.0358548
                ground_truth_k_vector_n(4,2)%z = 0.998714

                ground_truth_k_vector_n(4,3)%x = 0.0361547
                ground_truth_k_vector_n(4,3)%y = 0.
                ground_truth_k_vector_n(4,3)%z = 0.999346

                ground_truth_k_vector_n(4,4)%x = 0.0358548
                ground_truth_k_vector_n(4,4)%y = -0.0358548
                ground_truth_k_vector_n(4,4)%z = 0.998714


                print *, "test_lib_field_gaussian_beam_hermite_get_propagation_direction"
                rv = .true.

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        y = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        buffer_k_vector_n = lib_field_gaussian_beam_hermite_get_propagation_direction(gauss_parameter, &
                                                                                                      point_cartesian)

!                        k_vector_n(i,ii) = buffer_k_vector_n

                        buffer = ground_truth_k_vector_n(i,ii)%x - buffer_k_vector_n%x
                        if (abs(buffer) .lt. ground_truth_e ) then
                            print *, "  x(", i, ",", ii, "): OK"
                        else
                            print *, "  x(", i, ",", ii, "): FAILED"
                            print *, "     diff = ", buffer
                            print *, "       gt = ", ground_truth_k_vector_n(i,ii)%x
                            print *, "      k_x = ", buffer_k_vector_n%x
                            rv = .false.
                        end if

                        buffer = ground_truth_k_vector_n(i,ii)%y - buffer_k_vector_n%y
                        if (abs(buffer) .lt. ground_truth_e ) then
                            print *, "  y(" ,i, ",", ii, "): OK"
                        else
                            print *, "  y(" ,i, ",", ii, "): FAILED"
                            print * ,"     diff = ", buffer
                            print *, "       gt = ", ground_truth_k_vector_n(i,ii)%y
                            print *, "      k_y = ", buffer_k_vector_n%y
                            rv = .false.
                        end if

                        buffer = ground_truth_k_vector_n(i,ii)%z - buffer_k_vector_n%z
                        if (abs(buffer) .lt. ground_truth_e ) then
                            print *, "  z(", i, ",", ii, "): OK"
                        else
                            print *, "  z(", i, ",", ii, "): FAILED"
                            print * ,"     diff = ", buffer
                            print *, "       gt = ", ground_truth_k_vector_n(i,ii)%z
                            print *, "      k_z = ", buffer_k_vector_n%z
                            rv = .false.
                        end if

                    end do
                end do

            end function test_lib_field_gaussian_beam_hermite_get_propagation_direction

            function test_lib_field_gaussian_beam_hermite_get_plane_wave_approx() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                type(cartesian_coordinate_real_type) :: evaluation_point_x
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                type(lib_field_plane_wave_type) :: plane_wave_approximation

                type(cartesian_coordinate_cmplx_type) :: e_field_gauss
                type(cartesian_coordinate_cmplx_type) :: h_field_gauss
                type(cartesian_coordinate_cmplx_type) :: e_field_pw
                type(cartesian_coordinate_cmplx_type) :: h_field_pw

                double precision :: buffer

                evaluation_point_x = make_cartesian(0d0, 0d0, 5 * unit_mu)

                gauss_parameter%e_field_0 = 1
                gauss_parameter%wave_length_0 = 0.5 * unit_mu
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%waist_x0 = 2.5 * unit_mu
                gauss_parameter%waist_y0 = 2.5 * unit_mu
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
                gauss_parameter%theta = 0 !PI / 8d0
                gauss_parameter%phi = 0
                gauss_parameter%convention = 1

                plane_wave_approximation = lib_field_gaussian_beam_hermite_get_plane_wave_approximation(gauss_parameter, &
                                                                                                        evaluation_point_x)

                call lib_field_gaussian_beam_hermite_get_field(gauss_parameter, evaluation_point_x, &
                                                               e_field_gauss, h_field_gauss)
                call lib_field_plane_wave_get_field(plane_wave_approximation, evaluation_point_x, &
                                                    e_field_pw, h_field_pw)

                print *, " gauss: ", abs(e_field_gauss)
                print *, "approx: ", abs(e_field_pw)

                rv = .true.
                print *, "test_lib_field_gaussian_beam_hermite_get_plane_wave_approx: "

                buffer = abs(real(e_field_gauss%x) - real(e_field_pw%x))
                if (buffer .lt. ground_truth_e) then
                    print *, "  real(x): OK"
                else
                    print *, "  real(x): FAILED"
                    print *, "           diff = ", buffer
                    print *, "           gauss: ", real(e_field_gauss%x)
                    print *, "              pw: ", real(e_field_pw%x)
                    rv = .false.
                end if

                buffer = abs(aimag(e_field_gauss%x) - aimag(e_field_pw%x))
                if (buffer .lt. ground_truth_e) then
                    print *, "  aimag(x): OK"
                else
                    print *, "  aimag(x): FAILED"
                    print *, "            diff = ", buffer
                    print *, "            gauss: ", aimag(e_field_gauss%x)
                    print *, "               pw: ", aimag(e_field_pw%x)
                    rv = .false.
                end if


                buffer = abs(real(e_field_gauss%y) - real(e_field_pw%y))
                if (buffer .lt. ground_truth_e) then
                    print *, "  real(y): OK"
                else
                    print *, "  real(y): FAILED"
                    print *, "           diff = ", buffer
                    print *, "           gauss: ", real(e_field_gauss%y)
                    print *, "              pw: ", real(e_field_pw%y)
                    rv = .false.
                end if

                buffer = abs(aimag(e_field_gauss%y) - aimag(e_field_pw%y))
                if (buffer .lt. ground_truth_e) then
                    print *, "  aimag(y): OK"
                else
                    print *, "  aimag(y): FAILED"
                    print *, "            diff = ", buffer
                    print *, "            gauss: ", aimag(e_field_gauss%y)
                    print *, "               pw: ", aimag(e_field_pw%y)
                    rv = .false.
                end if


                buffer = abs(real(e_field_gauss%z) - real(e_field_pw%z))
                if (buffer .lt. ground_truth_e) then
                    print *, "  real(z): OK"
                else
                    print *, "  real(z): FAILED"
                    print *, "           diff = ", buffer
                    print *, "           gauss: ", real(e_field_gauss%z)
                    print *, "              pw: ", real(e_field_pw%z)
                    rv = .false.
                end if

                buffer = abs(aimag(e_field_gauss%z) - aimag(e_field_pw%z))
                if (buffer .lt. ground_truth_e) then
                    print *, "  aimag(z): OK"
                else
                    print *, "  aimag(z): FAILED"
                    print *, "            diff = ", buffer
                    print *, "            gauss: ", aimag(e_field_gauss%z)
                    print *, "               pw: ", aimag(e_field_pw%z)
                    rv = .false.
                end if

            end function test_lib_field_gaussian_beam_hermite_get_plane_wave_approx

            function test_lib_field_gaussian_beam_hermite_get_field_approximation() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_real_type) :: evaluation_point_x
                type(lib_field_plane_wave_type) :: plane_wave_approximation

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type) :: e_field_pw
                type(cartesian_coordinate_cmplx_type) :: h_field_pw

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.05_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 1.5 * unit_mu
                gauss_parameter%waist_y0 = 1.5 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = PI / 8d0
                gauss_parameter%phi = 0!PI

                gauss_parameter%convention = 1

                evaluation_point_x = make_cartesian(2 * unit_mu, 0d0, 5 * unit_mu)

                plane_wave_approximation = lib_field_gaussian_beam_hermite_get_plane_wave_approximation(gauss_parameter, &
                                                                                                        evaluation_point_x)

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(1) + (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter, &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field)

                        call lib_field_plane_wave_get_field(plane_wave_approximation, point_cartesian, &
                                                    e_field_pw, h_field_pw)

                        e_field(i,ii) = buffer_e_field - e_field_pw
                        h_field(i,ii) = buffer_h_field - h_field_pw
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_approx_")

            end function test_lib_field_gaussian_beam_hermite_get_field_approximation

            function test_lib_field_gaussian_beam_hermite_get_field_approx_full() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_real_type) :: evaluation_point_x
                type(lib_field_plane_wave_type) :: plane_wave_approximation

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type) :: e_field_pw
                type(cartesian_coordinate_cmplx_type) :: h_field_pw

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.05_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                gauss_parameter%waist_x0 = 2.5 * unit_mu
                gauss_parameter%waist_y0 = 2.5 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

                gauss_parameter%theta = 0!PI / 8d0
                gauss_parameter%phi = 0!PI

                gauss_parameter%convention = 1


!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        plane_wave_approximation = lib_field_gaussian_beam_hermite_get_plane_wave_approximation(&
                                                                                                        gauss_parameter, &
                                                                                                        point_cartesian)

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter, &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field)

                        call lib_field_plane_wave_get_field(plane_wave_approximation, point_cartesian, &
                                                    e_field_pw, h_field_pw)

!                        call lib_field_plane_wave_get_field(plane_wave_approximation, &
!                                                            make_cartesian(0d0, 0d0, 0d0), &
!                                                            e_field_pw, h_field_pw)

                        e_field(i,ii) = buffer_e_field - e_field_pw
                        h_field(i,ii) = buffer_h_field - h_field_pw
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_approx_full_")

            end function test_lib_field_gaussian_beam_hermite_get_field_approx_full

        end function lib_field_gaussian_beam_test_functions

end module lib_field_gaussian_beam
