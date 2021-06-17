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

module lib_field_plane_wave
    use libmath
    use lib_constants
    use lib_field_polarisation
    implicit none

    private

    ! public function
    public :: lib_field_plane_wave_get_field

    public :: lib_field_plane_wave_test_functions

    ! public types
    public :: lib_field_plane_wave_type

    type lib_field_plane_wave_type
        double precision :: e_field_0                   ! [V/m]
        double precision :: wave_length_0               ! [m]
        double precision :: refractive_index_medium     ! [1]
        double precision :: phase = 0
        type(jones_vector_type) :: polarisation         ! Jones vector
        double precision :: theta = 0                   ! propagation direction: polar angle [0, Pi] [rad]
        double precision :: phi = 0                     ! propagation direction: azimuthal angle [0, 2 Pi) [rad]
        ! 1: exp(-i(k*z - omega*t))
        ! 2: exp(i(k*z - omega*t))
        integer :: convention = 1
    end type

    contains

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermsubroutnite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !
        !          z            y
        !           ^ theta      ^    k
        !           |<->/ k      |   /
        !           |  /         |  /^
        !           | /          | / | phi
        !           |/           |/  v
        !           ------>      ------>
        !                  x            x
        !
        ! Returns
        ! ----
        !   e_field: type(cartesian_coordinate_cmplx_type)
        !       electical field component at "evaluation_point_x"
        !   h_field: type(cartesian_coordinate_cmplx_type)
        !       magnetical field component at "evaluation_point_x"
        !
        subroutine lib_field_plane_wave_get_field(parameter, evaluation_point_x, &
                                                  e_field, h_field)
            implicit none
            ! dummy
            type(lib_field_plane_wave_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            type(cartesian_coordinate_cmplx_type) :: e_field
            type(cartesian_coordinate_cmplx_type) :: h_field

            ! auxiliary
            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: point_x

            double precision :: wave_impedance
            double complex :: field

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(PI / 2d0 + parameter%phi, &
                                               parameter%theta, &
                                               PI / 2d0 + parameter%phi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            field = lib_field_plane_wave_scalar(parameter%e_field_0, &
                                                parameter%wave_length_0, &
                                                parameter%refractive_index_medium, &
                                                point_x%z, &
                                                phase=parameter%phase, &
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

        end subroutine lib_field_plane_wave_get_field

        ! Argument
        ! ----
        !   e_field_0: double precision
        !       electical field [V/m]
        !   wave_length: double precision
        !       vacuum wave length [m]
        !   n_medium: double precision
        !       refractive index of the medium
        !   phase: double precision
        !       additional phase of the plane wave
        !   convention: integer, optional (std: 1)
        !       1: exp(-i(k*z - omega*t))
        !       2: exp(i(k*z - omega*t))
        !
        ! Returns
        ! ----
        !   field: double cmplx
        !       field of a plane wave with a propagation along the z-axis at position "z"
        !
        !          z,k
        !           ^
        !        ___|___
        !        ___|___
        !           ---->
        !
        ! Reference: Experimentalphysik 2, Demtröder, eq. 7.9d
        function lib_field_plane_wave_scalar(e_field_0, wave_length_0, n_medium, z, phase, convention) &
                                             result (field)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: z
            double precision, intent(in), optional :: phase
            integer, intent(in), optional :: convention

            double complex :: field

            ! auxiliary
            double precision :: buffer_real

            double precision :: m_phase
            integer :: m_convention

            m_convention = 1
            if (present(convention)) m_convention = convention

            m_phase = 0
            if (present(phase)) m_phase = phase

            buffer_real = - 2d0 * PI * z * n_medium / wave_length_0 + m_phase

            select case(m_convention)
                case (1)
                    ! 1: exp(-i(k*z - omega*t))
                    field = e_field_0 * exp(dcmplx(0, buffer_real))
                case (2)
                    ! exp(i(k*z - omega*t))
                    field = e_field_0 * exp(dcmplx(0, -buffer_real))

                case default
                    print *, "lib_field_plane_wave_scalar: ERROR"
                    print *, "  convention is not defined: ", m_convention
             end select

        end function lib_field_plane_wave_scalar

        function lib_field_plane_wave_test_functions() result(rv)
            use lib_field
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_field_plane_wave_get_field()) rv = rv + 1
            if (.not. test_lib_field_plane_wave_get_field_c2()) rv = rv + 1

            if (rv == 0) then
                print *, "lib_field_plane_wave_test_functions tests: OK"
            else
                print *, rv,"lib_field_plane_wave_test_functions test(s) FAILED"
            end if

            contains

            function test_lib_field_plane_wave_get_field() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_plane_wave_type) :: plane_wave_parameter

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

                plane_wave_parameter%e_field_0 = 1
                plane_wave_parameter%refractive_index_medium = 1
!                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                plane_wave_parameter%wave_length_0 = 1 * unit_mu

                plane_wave_parameter%theta = PI / 8d0
                plane_wave_parameter%phi = 0

                plane_wave_parameter%convention = 1

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_plane_wave_get_field(plane_wave_parameter, &
                                                            point_cartesian, &
                                                            buffer_e_field, buffer_h_field)

                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/plane_wave_")

            end function test_lib_field_plane_wave_get_field

            function test_lib_field_plane_wave_get_field_c2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_plane_wave_type) :: plane_wave_parameter

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

                plane_wave_parameter%e_field_0 = 1
                plane_wave_parameter%refractive_index_medium = 1
!                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_rot(PI/4d0)
                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_linear_h()
!                plane_wave_parameter%polarisation = lib_field_polarisation_jones_vector_get_circular_plus()
                plane_wave_parameter%wave_length_0 = 1 * unit_mu

                plane_wave_parameter%theta = PI / 8d0
                plane_wave_parameter%phi = PI

                plane_wave_parameter%convention = 2

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= 1, no_y_values
                        z = y_range(2) - (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_plane_wave_get_field(plane_wave_parameter, &
                                                            point_cartesian, &
                                                            buffer_e_field, buffer_h_field)

                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/plane_wave_c2_")

            end function test_lib_field_plane_wave_get_field_c2
        end function
end module lib_field_plane_wave
