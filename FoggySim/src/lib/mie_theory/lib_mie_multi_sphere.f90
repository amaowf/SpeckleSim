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

module lib_mie_multi_sphere
    use libmath

    use lib_field
    use lib_field_plane_wave
    use lib_field_gaussian_beam
    use lib_field_polarisation

    use lib_mie_type
    use lib_mie_type_functions
    use lib_mie_ms_solver_interface
    use lib_mie_single_sphere

    use lib_mie_ms_data_container

    use lib_scene
    implicit none

    private


    public :: lib_mie_multi_sphere_test_functions
    public :: lib_mie_ms_benchmark

    logical :: m_use_ml_fmm
    logical :: m_init_with_single_sphere

    contains

!        ! Arguments
!        ! ----
!        !   k: type(spherical_coordinate_real_type)
!        !       ?
!        !   lambda: double precision
!        !       vacuum wave length
!        !   n_medium: double precision
!        !       refractive index of the medium
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
!        !   a_old, b_old: type(list_list_cmplx), optional
!        !       coefficient of the previous calculation, if f < 1
!        !   f: double precision, optional (std: 1)
!        !       numerical factor (0, 1]
!        !       "In our actual calculations, some multi-
!        !        sphere systems do not converge if f 5 1, but they do
!        !        converge when the value of f is reduced to, say, 0.7."[1]
!        ! Returns
!        ! ----
!        !
!        !
!        !
!        !
!        ! Reference: [1] Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
!        subroutine lib_mie_ms_get_interactive_scattering_coefficients(k, lambda, n_medium, &
!                                                                   sphere, sphere_parameter, sphere_j, &
!                                                                   z_selector, &
!                                                                   a, b)
!            implicit none
!            ! dummy
!            type(spherical_coordinate_real_type), intent(in) :: k
!            double precision, intent(in) :: lambda
!            double precision, intent(in) :: n_medium
!            type(lib_mie_sphere_type), dimension(:), intent(in) :: sphere
!            type(lib_mie_sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
!            integer :: sphere_j
!            integer(kind=1) :: z_selector
!
!            type(list_list_cmplx), intent(inout) :: a
!            type(list_list_cmplx), intent(inout) :: b
!
!            ! auxiliary
!            integer :: mu
!            integer :: nu
!
!
!        end subroutine lib_mie_ms_get_interactive_scattering_coefficients

        ! Argument
        ! ----
        !   use_ml_fmm: logical, optional (std: .false.)
        !       enables the ML FMM
        subroutine lib_mie_ms_constructor(n_range, use_ml_fmm, init_with_single_sphere)
            implicit none
            ! dummy
            integer, dimension(2), intent(in) :: n_range
            logical, intent(in), optional :: use_ml_fmm
            logical, intent(in), optional :: init_with_single_sphere

            ! auxiliary
            type(lib_mie_ms_solver_parameter_type) :: solver_parameter

            m_use_ml_fmm = .false.
            if (present(use_ml_fmm)) m_use_ml_fmm = use_ml_fmm

            m_init_with_single_sphere = .false.
            if (present(init_with_single_sphere)) m_init_with_single_sphere = init_with_single_sphere

            call fwig_table_init(4 * n_range(2), 3)
            call fwig_temp_init(4 * n_range(2))

            call lib_math_factorial_initialise_caching(2 * n_range(2))

            call lib_mie_ss_constructor(n_range)

            if (m_init_with_single_sphere) then
                call lib_mie_ss_calculate_scattering_coefficients_ab_nm()
            end if

            solver_parameter%use_initial_guess = m_init_with_single_sphere
            solver_parameter%convergence_tolerance = 1D-6
            solver_parameter%use_ml_fmm = use_ml_fmm

            call lib_mie_ms_solver_constructor(solver_parameter)

!            ! todo: init caching
!            !        - Tree: theta, phi
!            if (m_use_ml_fmm) then
!
!            end if

        end subroutine

        ! Argument
        ! ----
        !   update_initial_guess_with_n_max: integer, optional
        !       calculates the initial guess with n_max = *update_initial_guess_with_n_max*
        !   step_size: integer, optional
        !       if update_initial_guess_with_n_max was set, further iterations start with
        !           n_max = update_initial_guess_with_n_max + i * step_size, i = 1, ...
        !   dynamic_step_size: logical, optional (std: .false.)
        !       If the next step with *step_size* is too large, this is detected by a significant reduction in the loss,
        !       the *step* for the current iteration is reduced.
        subroutine lib_mie_ms_calculate_scattering_coefficients_ab_nm(update_initial_guess_with_n_max, &
                                                                      step_size, dynamic_step_size)
            use file_io
            implicit none
            ! dummy

            integer, intent(in), optional :: update_initial_guess_with_n_max
            integer, intent(in), optional :: step_size
            logical, intent(in), optional :: dynamic_step_size

            ! auxiliary
            integer, dimension(2) :: n_range
            integer, dimension(2) :: n_range_org
            double precision :: old_backward_error

            integer :: m_step_size
            logical :: m_dynamic_step_size

            ! initial values
!            call lib_mie_ss_calculate_scattering_coefficients_ab_nm()

            if ( present(update_initial_guess_with_n_max) ) then
                n_range_org = simulation_data%spherical_harmonics%n_range
                n_range = n_range_org

                n_range(2) = update_initial_guess_with_n_max

                if (present(step_size)) then

                    m_step_size = step_size

                    m_dynamic_step_size = .false.
                    if (present(dynamic_step_size)) then
                        m_dynamic_step_size = dynamic_step_size
                    end if

                    if (m_dynamic_step_size) then
                        call calc_with_dynmaic_step_size()
                    else
                        call calc_with_steps()
                    end if

                else

                    call lib_mie_ms_solver_set_n_range(n_range)

                    print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                    print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                    print *, "  "

                    call lib_mie_ms_solver_set_max_iterations(100)
                    call lib_mie_ms_solver_run()

                    n_range = n_range_org
                    call lib_mie_ms_solver_set_n_range(n_range)

                    print *, ""
                    print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                    print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                    print *, ""
                    call lib_mie_ms_solver_set_max_iterations(1000)
                    call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                    call lib_mie_ms_solver_run()

                end if

            else
                print *, ""
                print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                print *, ""
                call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                call lib_mie_ms_solver_run()
            end if

            contains

                subroutine calc_with_steps()
                    implicit none
                    ! auxiliary
                    integer :: n

                    do n=update_initial_guess_with_n_max, n_range_org(2) - m_step_size, step_size
                        n_range(2) = n
                        call lib_mie_ms_solver_set_n_range(n_range)

                        print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                        print *, "n = ", simulation_data%spherical_harmonics%n_range(2)
                        print *, ""

                        call lib_mie_ms_solver_set_max_iterations(100)
                        call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                        call lib_mie_ms_solver_run()
                    end do

                    n_range = n_range_org
                    call lib_mie_ms_solver_set_n_range(n_range)

                    print *, ""
                    print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                    print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                    print *, ""
                    call lib_mie_ms_solver_set_max_iterations(1000)
                    call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                    call lib_mie_ms_solver_run()
                end subroutine calc_with_steps

                subroutine calc_with_dynmaic_step_size()
                    implicit none
                    ! parameter
                    integer, parameter :: state_calculate = 1
                    integer, parameter :: state_test_next_step = 2
                    integer, parameter :: state_step_size_too_big = 3
                    integer, parameter :: state_finish = 4

                    ! auxiliaray
                    integer :: state

                    call lib_mie_ms_solver_set_n_range(n_range)

                    print *, ""
                    print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                    print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                    print *, ""

                    call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                    call lib_mie_ms_solver_run()

                    state = state_test_next_step

                    do
                        if (state .eq. state_test_next_step) then
                            if (n_range(2) + m_step_size .ge. n_range_org(2)) then
                                n_range = n_range_org
                                call lib_mie_ms_solver_set_n_range(n_range)
                            else
                                n_range(2) = n_range(2) + m_step_size
                                call lib_mie_ms_solver_set_n_range(n_range)
                            end if

                            if (m_step_size .gt. 1) then
                                print *, ""
                                print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                                print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                                print *, ""

                                call lib_mie_ms_solver_set_max_iterations(1)
                                call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                                call lib_mie_ms_solver_run()

                                if (lib_mie_ms_solver_get_backward_error() .lt. 1D-3) then
                                    if (n_range(2) + m_step_size .ge. n_range_org(2)) then
                                        n_range = n_range_org
                                    else
                                        n_range(2) = n_range(2) + m_step_size
                                    end if
                                    state = state_calculate
                                else
                                    if (m_step_size .gt. 1) then
                                        state = state_step_size_too_big
                                    else
                                        state = state_calculate
                                    end if
                                end if

                            else
                                if (n_range(2) + m_step_size .ge. n_range_org(2)) then
                                    n_range = n_range_org
                                else
                                    n_range(2) = n_range(2) + m_step_size
                                end if
                                state = state_calculate
                            end if

                        else if (state .eq. state_calculate) then
                            call lib_mie_ms_solver_set_n_range(n_range)

                            print *, ""
                            print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: NOTE"
                            print *, "  n = ", simulation_data%spherical_harmonics%n_range(2)
                            print *, ""

                            call lib_mie_ms_solver_set_max_iterations(1000)
                            call lib_mie_ms_solver_use_ml_fmm(m_use_ml_fmm)
                            call lib_mie_ms_solver_run()
                            old_backward_error = lib_mie_ms_solver_get_backward_error()
                            m_step_size = step_size

                            if (n_range(2) .eq. n_range_org(2)) then
                                state = state_finish
                            else
                                state = state_test_next_step
                            end if

                        else if (state .eq. state_step_size_too_big) then
                            if (int(m_step_size / 2D0) .gt. 1) then
                                m_step_size = int(m_step_size / 2D0)
                            else
                                m_step_size = 1
                            end if

                            state = state_test_next_step

                        else if (state .eq. state_finish) then
                            exit
                        end if
                    end do
                end subroutine calc_with_dynmaic_step_size

        end subroutine lib_mie_ms_calculate_scattering_coefficients_ab_nm

        ! Argument
        ! ----
        !   simulation: type(lib_mie_simulation_parameter_type)
        !       simulation data set
        !   x_0: cartesian_coordinate_real_type
        !
        function lib_mie_ms_get_field(x_0) result(field)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: x_0

            type(cartesian_coordinate_cmplx_type), dimension(2) :: field

            ! auxiliaray
            integer :: sphere_no
            integer :: parameter_no

            double precision :: e_field_0
            double precision :: k_0
            double precision :: n_medium

            type(cartesian_coordinate_real_type) :: x_j

            type(lib_mie_sphere_type) :: sphere
            type(lib_mie_sphere_parameter_type) :: sphere_parameter

            type(spherical_coordinate_cmplx_type), dimension(:, :), allocatable :: buffer_field


            e_field_0 = simulation_data%illumination%e_field_0
            k_0 = 2 * PI / simulation_data%illumination%lambda_0
            n_medium = simulation_data%refractive_index_medium

            allocate(buffer_field(lbound(simulation_data%sphere_list, 1):ubound(simulation_data%sphere_list, 1), 2))

            field(:)%x = dcmplx(0,0)
            field(:)%y = dcmplx(0,0)
            field(:)%z = dcmplx(0,0)

            do sphere_no = lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                sphere = simulation_data%sphere_list(sphere_no)
                parameter_no = sphere%sphere_parameter_index
                sphere_parameter = simulation_data%sphere_parameter_list(parameter_no)

                x_j = x_0 - sphere%d_0_j

                if (abs(x_j) .gt. sphere_parameter%radius) then
                    buffer_field(sphere_no, 1:2) = lib_mie_ss_get_field(x_0, &
                                                                     e_field_0, k_0, n_medium, &
                                                                     sphere, sphere_parameter)

                    field(1) = field(1) + make_cartesian(buffer_field(sphere_no, 1), x_j)
                    field(2) = field(2) + make_cartesian(buffer_field(sphere_no, 2), x_j)

                else
                    field(:)%x = dcmplx(0,0)
                    field(:)%y = dcmplx(0,0)
                    field(:)%z = dcmplx(0,0)
                    exit
                end if
            end do
        end function

        function lib_mie_multi_sphere_test_functions() result(rv)
            use file_io
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)
            ! CPU-time
            real :: test_start, test_finish
            ! WALL-time
            INTEGER :: test_count_start, test_count_finish, test_count_rate

            rv = 0

            call system_clock(test_count_start, test_count_rate)
            call cpu_time(test_start)

!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene()) rv = rv + 1
!            if (.not. test_lib_mie_ms_calculate_scattering_coeff_ab_nm_scene_eval()) rv = rv + 1
!            if (.not. test_lib_mie_ms_get_field_parallel_sphere_assemply()) rv = rv + 1
            if (.not. test_lib_mie_ms_get_field_sphere_grid_assemply()) rv = rv + 1
!            if (.not. test_lib_mie_ms_get_field_serial_sphere_assemply()) rv = rv + 1

            call cpu_time(test_finish)
            call system_clock(test_count_finish, test_count_rate)

            print *, ""
            print *, "------lib_mie_multi_sphere_test_functions------"
            print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
            print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
            print *, ""
            if (rv == 0) then
                print *, "lib_mie_multi_sphere_test_functions tests: OK"
            else
                print *, rv,"lib_mie_multi_sphere_test_functions test(s) FAILED"
            end if
            print *, "------------------------------------------------------------"
            print *, ""

            contains

            function test_lib_mie_ms_get_field_parallel_sphere_assemply() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_z_values

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 0.55 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 0
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 1
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                allocate(simulation_data%sphere_list(7))
                simulation_data%sphere_list(1)%sphere_parameter_index = 1
                simulation_data%sphere_list(2)%sphere_parameter_index = 1
                simulation_data%sphere_list(3)%sphere_parameter_index = 2
                simulation_data%sphere_list(4)%sphere_parameter_index = 1
                simulation_data%sphere_list(5)%sphere_parameter_index = 1
                simulation_data%sphere_list(6)%sphere_parameter_index = 2
                simulation_data%sphere_list(7)%sphere_parameter_index = 2

                sphere_d_0_j%x = -1.05 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 0
                simulation_data%sphere_list(1)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 1.05 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 0
                simulation_data%sphere_list(2)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 0
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 3 * unit_mu
                simulation_data%sphere_list(3)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = -1.5 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 3 * unit_mu
                simulation_data%sphere_list(4)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 1.5 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 3 * unit_mu
                simulation_data%sphere_list(5)%d_0_j = sphere_d_0_j
                sphere_d_0_j%x = -3 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 5 * unit_mu
                simulation_data%sphere_list(6)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 3 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 5 * unit_mu
                simulation_data%sphere_list(7)%d_0_j = sphere_d_0_j

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 0.5 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                n_range = simulation_data%spherical_harmonics%n_range
                n_range(2) = 6
                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)


                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(5), 6)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(n_range(2) / 2), 4)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(n_range(2) / 2), 2)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(, int(n_range(2) / 2))
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8)

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)

                print *, ""
                print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export
                x_range = (/ -5.0_8 * unit_mu, 5.0_8 * unit_mu /)
                z_range = (/ -5_8 * unit_mu, 10.0_8 * unit_mu /)
                step_size = 0.1_8 * unit_mu

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        call lib_field_plane_wave_get_field(simulation_data%illumination%plane_wave(1)%beam_parameter, &
                                                            x_0, &
                                                            field(1), field(2))

                        e_field_s(i,ii) = e_field_s(i,ii) + field(1)
                        h_field_s(i,ii) = h_field_s(i,ii) + field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/realTotal/")

                call fwig_temp_free()
                call fwig_table_free()

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_get_field_parallel_sphere_assemply

            function test_lib_mie_ms_get_field_sphere_grid_assemply() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_z_values

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 0.55 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 10 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 1
                no_spheres_z = 1
                distance_sphere = 5 * unit_mu
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do i = 1, no_spheres_x
                    do ii = 1, no_spheres_z
                        no = no + 1

                        if (mod(no, 2) .eq. 0) then
                            simulation_data%sphere_list(no)%sphere_parameter_index = 1
                        else
                            simulation_data%sphere_list(no)%sphere_parameter_index = 2
                        end if

                        sphere_d_0_j%x = i * distance_sphere
                        sphere_d_0_j%z = ii * distance_sphere

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                n_range = simulation_data%spherical_harmonics%n_range

                simulation_data%spherical_harmonics%n_range(2) = 15 !n_range(2) / 4

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)


                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(5), 6)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(n_range(2) / 3))
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(, int(n_range(2) / 2), 2)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(, int(n_range(2) / 2))
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)

                print *, ""
                print *, "lib_mie_ms_calculate_scattering_coefficients_ab_nm: "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export
                x_range = (/ 0_8 * unit_mu, (no_spheres_x+1) * distance_sphere /)
                z_range = (/ 0_8 * unit_mu, (no_spheres_z+1) * distance_sphere /)
                step_size = (x_range(2) - x_range(1)) / 201

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")

                call fwig_temp_free()
                call fwig_table_free()

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_get_field_sphere_grid_assemply

            function test_lib_mie_ms_get_field_serial_sphere_assemply() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_z_values

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 1
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 3

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 1 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 0
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                allocate(simulation_data%sphere_list(2))
                simulation_data%sphere_list(1)%sphere_parameter_index = 1
                simulation_data%sphere_list(2)%sphere_parameter_index = 1

                sphere_d_0_j%x = 0
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 0
                simulation_data%sphere_list(1)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 0
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 3 * unit_mu
                simulation_data%sphere_list(2)%d_0_j = sphere_d_0_j

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 0.5 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)
                n_range = simulation_data%spherical_harmonics%n_range
                call lib_mie_ms_constructor(n_range, use_ml_fmm=.false., init_with_single_sphere=.true.)


                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()

                ! evaluate and export
                x_range = (/ -5.0_8 * unit_mu, 5.0_8 * unit_mu /)
                z_range = (/ -5_8 * unit_mu, 10.0_8 * unit_mu /)
                step_size = 0.05_8 * unit_mu

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_get_field_serial_sphere_assemply
            ! calculates the scattering coefficients and compares the T-Matrix method against the analytical Mie solution
            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no

                integer :: n
                integer :: m

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                type(list_list_cmplx) :: a_nm_mie
                type(list_list_cmplx) :: b_nm_mie

                type(list_list_cmplx) :: a_nm_t_matrix
                type(list_list_cmplx) :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff
                double precision :: buffer

                double complex, dimension(:), allocatable :: buffer_list_mie
                double complex, dimension(:), allocatable :: buffer_list_t_matrix

                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 0.55 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 10 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 1
                no_spheres_z = 1
                distance_sphere = 4 * unit_mu
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do i = 1, no_spheres_x
                    do ii = 1, no_spheres_z
                        no = no + 1

!                        if (mod(no, 2) .eq. 0) then
                            simulation_data%sphere_list(no)%sphere_parameter_index = 1
!                        else
!                            simulation_data%sphere_list(no)%sphere_parameter_index = 2
!                        end if

                        sphere_d_0_j%x = -i * distance_sphere
                        sphere_d_0_j%z = ii * distance_sphere

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 1000 * unit_mu
                n_particle = dcmplx(1.5_8, 0)

!                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
!                if (n_range(1) .gt. 0) then
!                    n_range(2) = n_range(1)
!                    n_range(1) = 1
!                else
!                    print *, "test_lib_mie_ms_get_field: ERROR"
!                    print *, "  lib_mie_ss_test_convergence_plane_wave"
!                    print *, "  rv(1): ", n_range(1)
!                    print *, "  rv(2): ", n_range(2)
!                end if

                n_range(1) = 1
                n_range(2) = 45

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
!                r_particle = 5 * unit_mu
!                n_particle = dcmplx(1.33_8, 0)
!
!                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
!                if (n_range(1) .gt. 0) then
!                    n_range(2) = n_range(1)
!                    n_range(1) = 1
!                else
!                    print *, "test_lib_mie_ms_get_field: ERROR"
!                    print *, "  lib_mie_ss_test_convergence_plane_wave"
!                    print *, "  rv(1): ", n_range(1)
!                    print *, "  rv(2): ", n_range(2)
!                end if
!
!                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
!                                                                                          r_particle, n_particle, &
!                                                                                          n_range)
!
!                n_range = simulation_data%spherical_harmonics%n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(5), 6)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(n_range(2) / 2))
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(, int(n_range(2) / 2), 2)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(, int(n_range(2) / 2))
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()

!                a_nm_t_matrix = simulation_data%sphere_list(1)%a_nm
!                b_nm_t_matrix = simulation_data%sphere_list(1)%b_nm

                a_nm_mie = simulation_data%sphere_list(1)%a_nm
                b_nm_mie = simulation_data%sphere_list(1)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1: "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                call make_array((/a_nm_mie, b_nm_mie/), buffer_list_mie)
!                call make_array((/a_nm_t_matrix, b_nm_t_matrix/), buffer_list_t_matrix)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1_mie.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_mie)
                close(99)

!                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1_t_matrix.csv"), &
!                     status='unknown')
!                rv = write_csv_cmplx_array_1d(99, buffer_list_t_matrix)
!                close(99)
!
!                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1_diff.csv"), &
!                     status='unknown')
!                rv = write_csv_cmplx_array_1d(99, buffer_list_mie - buffer_list_t_matrix)
!                close(99)

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1:"
!                do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
!                    print *, "  i = ", i
!
!                    print *, "  a_nm"
!                    list_list_diff = a_nm_mie - a_nm_t_matrix
!
!                    do n = lbound(list_list_diff%item, 1), &
!                           ubound(list_list_diff%item, 1)
!                        print *, "  n = ", n
!                        do m = -n, n
!                            buffer = abs(list_list_diff%item(n)%item(m))
!                            if (buffer .gt. ground_truth_e) then
!                                print *, "    m: ", m , "difference: ", buffer, " : FAILED"
!!                                    print *, "            Mie: ", a_nm_mie%item(n)%item(m)
!!                                    print *, "       T-Matrix: ", a_nm_t_matrix%item(n)%item(m)
!                                rv = .false.
!                            else
!                                print *, "    m: ", m, ": OK"
!!                                    print *, "            Mie: ", a_nm_mie%item(n)%item(m)
!!                                    print *, "       T-Matrix: ", a_nm_t_matrix%item(n)%item(m)
!                            end if
!                        end do
!                    end do
!
!                    print *, "  b_nm"
!                    list_list_diff = b_nm_mie - b_nm_t_matrix
!
!                    do n = lbound(list_list_diff%item, 1), &
!                           ubound(list_list_diff%item, 1)
!                        print *, "  n = ", n
!                        do m = -n, n
!                            buffer = abs(list_list_diff%item(n)%item(m))
!                            if (buffer .gt. ground_truth_e) then
!                                print *, "    m: ", m , "difference: ", buffer, " : FAILED"
!!                                    print *, "            Mie: ", b_nm_mie%item(n)%item(m)
!!                                    print *, "       T-Matrix: ", b_nm_t_matrix%item(n)%item(m)
!                                rv = .false.
!                            else
!                                print *, "    m: ", m, ": OK"
!!                                    print *, "            Mie: ", a_nm_mie%item(n)%item(m)
!!                                    print *, "       T-Matrix: ", a_nm_t_matrix%item(n)%item(m)
!                            end if
!                        end do
!                    end do
!                end do

                x_range = (/ 0_8 * unit_mu, (no_spheres_x+1) * distance_sphere /)
                z_range = (/ 0_8 * unit_mu, (no_spheres_z+1) * distance_sphere /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real_v1/")

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v1

            ! calculates the scattering coefficients and compares the T-Matrix method against the ML FMM (1 level)
            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no

                integer :: n
                integer :: m

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff
                double precision :: buffer

                character(len=15) :: path

                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                double complex, dimension(:), allocatable :: buffer_list_t_matrix
                double complex, dimension(:), allocatable :: buffer_list_ml_fmm

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 1 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 10 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 3
                no_spheres_z = 5

                allocate (a_nm_t_matrix(no_spheres_x * no_spheres_z))
                allocate (b_nm_t_matrix(no_spheres_x * no_spheres_z))

                allocate (a_nm_ml_fmm(no_spheres_x * no_spheres_z))
                allocate (b_nm_ml_fmm(no_spheres_x * no_spheres_z))

                distance_sphere = 5 * unit_mu
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do i = 1, no_spheres_x
                    do ii = 1, no_spheres_z
                        no = no + 1

                        if (mod(no, 2) .eq. 0) then
                            simulation_data%sphere_list(no)%sphere_parameter_index = 1
                        else
                            simulation_data%sphere_list(no)%sphere_parameter_index = 2
                        end if

                        sphere_d_0_j%x = -i * distance_sphere
                        sphere_d_0_j%z = ii * distance_sphere

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 0.3 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 0.15 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                n_range = simulation_data%spherical_harmonics%n_range
                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(ceiling(n_range(2) * 0.75)), 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_t_matrix(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_t_matrix(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)

                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2 (T-Matrix): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(int(ceiling(n_range(2) * 0.75)), 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2 (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export

                call make_array((/a_nm_ml_fmm, b_nm_ml_fmm/), buffer_list_ml_fmm)
                call make_array((/a_nm_t_matrix, b_nm_t_matrix/), buffer_list_t_matrix)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2_ml_fmm.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2_t_matrix.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_t_matrix)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2_diff.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm - buffer_list_t_matrix)
                close(99)

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2:"
                do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                    print *, "  i = ", i

                    print *, "  a_nm"
                    list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)

                    do n = lbound(list_list_diff%item, 1), &
                           ubound(list_list_diff%item, 1)
                        print *, "  n = ", n
                        do m = -n, n
                            buffer = real(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (real): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, " (real): OK"
                            end if

                            buffer = aimag(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (imag): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, " (imag): OK"
                            end if
                        end do
                    end do

                    print *, "  b_nm"
                    list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)

                    do n = lbound(list_list_diff%item, 1), &
                           ubound(list_list_diff%item, 1)
                        print *, "  n = ", n
                        do m = -n, n
                            buffer = real(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (real): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, " (real): OK"
                            end if

                            buffer = aimag(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (imag): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, " (imag): OK"
                            end if
                        end do
                    end do
                end do

                x_range = (/ 0_8 * unit_mu, (no_spheres_x+1) * distance_sphere /)
                z_range = (/ 0_8 * unit_mu, (no_spheres_z+1) * distance_sphere /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = -x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real_v2/")

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v2

            ! calculates the scattering coefficients and compares the T-Matrix method against the ML FMM
            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no

                integer :: n
                integer :: m

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff
                double precision :: buffer

                character(len=15) :: path

                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                double complex, dimension(:), allocatable :: buffer_list_t_matrix
                double complex, dimension(:), allocatable :: buffer_list_ml_fmm

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 1 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 10 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 3
                no_spheres_z = 2

                allocate (a_nm_t_matrix(no_spheres_x * no_spheres_z))
                allocate (b_nm_t_matrix(no_spheres_x * no_spheres_z))

                allocate (a_nm_ml_fmm(no_spheres_x * no_spheres_z))
                allocate (b_nm_ml_fmm(no_spheres_x * no_spheres_z))

                distance_sphere = 10 * unit_mu
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do i = 1, no_spheres_x
                    do ii = 1, no_spheres_z
                        no = no + 1

                        if (mod(no, 2) .eq. 0) then
                            simulation_data%sphere_list(no)%sphere_parameter_index = 1
                        else
                            simulation_data%sphere_list(no)%sphere_parameter_index = 2
                        end if

                        sphere_d_0_j%x = -i * distance_sphere
                        sphere_d_0_j%z = ii * distance_sphere

                        sphere_d_0_j%y = 2*(ii-1) * unit_mu

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 3 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 2 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

!                n_range = simulation_data%spherical_harmonics%n_range
                n_range(2) = 12 !int(ceiling(real(n_range(2)) * 1.2))
                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(10, 1)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_t_matrix(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_t_matrix(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)

                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3 (T-Matrix): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3 (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export
                x_range = (/ 0_8 * unit_mu, (no_spheres_x+1) * distance_sphere /)
                z_range = (/ 0_8 * unit_mu, (no_spheres_z+1) * distance_sphere /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = -x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real_v3/")

                call make_array((/a_nm_ml_fmm, b_nm_ml_fmm/), buffer_list_ml_fmm)
                call make_array((/a_nm_t_matrix, b_nm_t_matrix/), buffer_list_t_matrix)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3_ml_fmm.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3_t_matrix.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_t_matrix)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3_diff.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm - buffer_list_t_matrix)
                close(99)

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3:"
                    do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                        print *, "  i = ", i

                        print *, "  a_nm"
                        list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do

                        print *, "  b_nm"
                        list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do
                    end do

                    call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v3

            ! calculates the scattering coefficients and compares the T-Matrix method against the ML FMM
            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no

                integer :: n
                integer :: m

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff
                double precision :: buffer

                character(len=15) :: path

                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                double complex, dimension(:), allocatable :: buffer_list_t_matrix
                double complex, dimension(:), allocatable :: buffer_list_ml_fmm

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 550 * unit_nm

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 10 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, n_medium, &
                                                                                             e_field_0, &
                                                                                             plane_wave_g, &
                                                                                             plane_wave_k, &
                                                                                             plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 10
                no_spheres_z = 10

                allocate (a_nm_t_matrix(no_spheres_x * no_spheres_z))
                allocate (b_nm_t_matrix(no_spheres_x * no_spheres_z))

                allocate (a_nm_ml_fmm(no_spheres_x * no_spheres_z))
                allocate (b_nm_ml_fmm(no_spheres_x * no_spheres_z))

                distance_sphere = 15.1 * unit_nm
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do ii = 1, no_spheres_z
                    do i = 1, no_spheres_x
                        no = no + 1
                        simulation_data%sphere_list(no)%sphere_parameter_index = 1

                        if (mod(ii, 2) .eq. 0) then
                            sphere_d_0_j%x = i * distance_sphere + distance_sphere / 4d0
                        else
                            sphere_d_0_j%x = i * distance_sphere - distance_sphere / 4d0
                        end if
                        sphere_d_0_j%z = ii * distance_sphere

!                        sphere_d_0_j%y = 2*(ii-1) * unit_mu

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 7.5 * unit_nm
                n_particle = dcmplx(2.5287_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                print *, " old: n_max = ", n_range(2)
                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
                print *, " new: n_max = ", n_range(2)

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)


!                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)
!
!                call system_clock(test_count_start_sub, test_count_rate_sub)
!                call cpu_time(test_start_sub)
!
!!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(10, 1)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
!                a_nm_t_matrix(:) = simulation_data%sphere_list(:)%a_nm
!                b_nm_t_matrix(:) = simulation_data%sphere_list(:)%b_nm
!
!                call cpu_time(test_finish_sub)
!                call system_clock(test_count_finish_sub, test_count_rate_sub)
!
!                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4 (T-Matrix): "
!                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
!                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
!                                                               / real(test_count_rate_sub)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
!                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4 (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export
                x_range = (/ 0_8 * unit_mu, (no_spheres_x) * distance_sphere * 1.2d0 /)
                z_range = (/ 0_8 * unit_mu, (no_spheres_z) * distance_sphere * 1.2d0 /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real_v4/")

                call make_array((/a_nm_ml_fmm, b_nm_ml_fmm/), buffer_list_ml_fmm)
                call make_array((/a_nm_t_matrix, b_nm_t_matrix/), buffer_list_t_matrix)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4_ml_fmm.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4_t_matrix.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_t_matrix)
                close(99)

                open(unit=99, file=trim("temp/test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4_diff.csv"), &
                     status='unknown')
                rv = write_csv_cmplx_array_1d(99, buffer_list_ml_fmm - buffer_list_t_matrix)
                close(99)

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4:"
                    do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                        print *, "  i = ", i

                        print *, "  a_nm"
                        list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do

                        print *, "  b_nm"
                        list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do
                    end do

                    call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_v4

            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                integer :: no

                integer :: n
                integer :: m

                integer :: no_spheres_x
                integer :: no_spheres_z
                double precision :: distance_sphere

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff
                double precision :: buffer

                character(len=15) :: path

                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                type(lib_mie_illumination_parameter) :: illumination

                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub


                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 550 * unit_nm

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                allocate(illumination%gaussian_beam(1))

                illumination%lambda_0 = lambda_0
                illumination%e_field_0 = 1

                illumination%gaussian_beam(:)%g = 1

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 0
                illumination%gaussian_beam(:)%d_0_i = buffer_car

                illumination%gaussian_beam(:)%calculation_type = 3

                illumination%gaussian_beam(:)%beam_parameter%refractive_index_medium = 1

                illumination%gaussian_beam(:)%beam_parameter%e_field_0 = 1
                illumination%gaussian_beam(:)%beam_parameter%phase = 0
                illumination%gaussian_beam(:)%beam_parameter%wave_length_0 = lambda_0
                illumination%gaussian_beam(:)%beam_parameter%polarisation = &
                                                lib_field_polarisation_jones_vector_get_linear_h()

                illumination%gaussian_beam(1)%beam_parameter%waist_x0 = 2 * unit_mu
                illumination%gaussian_beam(1)%beam_parameter%waist_y0 = 2 * unit_mu
                illumination%gaussian_beam(1)%beam_parameter%tem_m = 0
                illumination%gaussian_beam(1)%beam_parameter%tem_m = 0
                illumination%gaussian_beam(1)%beam_parameter%convention = 2

                simulation_data%illumination = illumination

                simulation_data%refractive_index_medium = n_medium

                ! set spheres
                no_spheres_x = 3
                no_spheres_z = 2

                allocate (a_nm_t_matrix(no_spheres_x * no_spheres_z))
                allocate (b_nm_t_matrix(no_spheres_x * no_spheres_z))

                allocate (a_nm_ml_fmm(no_spheres_x * no_spheres_z))
                allocate (b_nm_ml_fmm(no_spheres_x * no_spheres_z))

                distance_sphere = 20 * unit_mu
                allocate(simulation_data%sphere_list(no_spheres_x * no_spheres_z))

                sphere_d_0_j%y = 0

                no = 0
                do ii = 1, no_spheres_z
                    do i = 1, no_spheres_x
                        no = no + 1
                        simulation_data%sphere_list(no)%sphere_parameter_index = 1

                        if (mod(ii, 2) .eq. 0) then
                            sphere_d_0_j%x = i * distance_sphere + distance_sphere / 4d0 - distance_sphere * no_spheres_x / 2
                        else
                            sphere_d_0_j%x = i * distance_sphere - distance_sphere / 4d0 - distance_sphere * no_spheres_x / 2
                        end if
                        sphere_d_0_j%z = ii * distance_sphere

                        simulation_data%sphere_list(no)%d_0_j = sphere_d_0_j

                    end do
                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(2))

                ! set 1
                r_particle = 2 * unit_mu
                n_particle = dcmplx(2.5287_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)


                n_range(2) = 12 !int(ceiling(real(n_range(2)) * 1.2))
                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(10, 1)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_t_matrix(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_t_matrix(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)

                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss (T-Matrix): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)

!                call system_clock(test_count_start_sub, test_count_rate_sub)
!                call cpu_time(test_start_sub)
!
!                n_range(2) = 12 !int(ceiling(real(n_range(2)) * 1.2))
!                simulation_data%spherical_harmonics%n_range = n_range
!
!                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)
!
!!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
!                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
!                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm
!
!                call cpu_time(test_finish_sub)
!                call system_clock(test_count_finish_sub, test_count_rate_sub)
!
!
!                print *, ""
!                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss (ML-FMM): "
!                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
!                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
!                                                               / real(test_count_rate_sub)
!                print *, ""

                ! evaluate and export
                x_range = (/ -distance_sphere * 1.2d0 - distance_sphere * no_spheres_x / 2, &
                             no_spheres_x * distance_sphere * 1.2d0 - distance_sphere * no_spheres_x / 2 /)
                z_range = (/ -distance_sphere * 1.2d0, &
                             (no_spheres_z + 1) * distance_sphere * 1.2d0 /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        call lib_field_gaussian_beam_hermite_get_field(illumination%gaussian_beam(1)%beam_parameter,&
                                                                       x_0, &
                                                                       field(1), field(2))

                        e_field_s(i,ii) = e_field_s(i,ii) + field(1)
                        h_field_s(i,ii) = h_field_s(i,ii) + field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/realTotal/")

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss:"
                    do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                        print *, "  i = ", i

                        print *, "  a_nm"
                        list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do

                        print *, "  b_nm"
                        list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)

                        do n = lbound(list_list_diff%item, 1), &
                               ubound(list_list_diff%item, 1)
                            print *, "  n = ", n
                            do m = -n, n
                                buffer = abs(list_list_diff%item(n)%item(m))
                                if (buffer .gt. ground_truth_e) then
                                    print *, "    m: ", m , "difference: ", buffer, " : FAILED"
                                    print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                    print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                    rv = .false.
                                else
                                    print *, "    m: ", m, ": OK"
                                end if
                            end do
                        end do
                    end do

                    call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_gauss

            function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene() result(rv)
                use file_io
                use lib_mie_ms_ml_fmm_interface
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer, dimension(2) :: range_i
                integer :: ii

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix


                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_y_values
                integer :: no_z_values

                type(lib_scene_opject_type) :: scene
                type(cartesian_coordinate_real_type) :: d_o_j
                type(cartesian_coordinate_real_type) :: cap_plane_point
                type(cartesian_coordinate_real_type) :: cap_plane_normal
                double precision :: sphere_radius
                double precision :: lattice_sphere_radius

                type(cartesian_coordinate_real_type), dimension(:), allocatable :: sphere_coordinates
                type(cartesian_coordinate_real_type) :: buffer_coordinate
                double precision, dimension(:,:), allocatable :: data_list

                integer :: u
                character(len=50) :: filename


                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 1550 * unit_nm

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, &
                                                                                        n_medium, &
                                                                                        e_field_0, &
                                                                                        plane_wave_g, &
                                                                                        plane_wave_k, &
                                                                                        plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres

                scene%d_0_o = make_cartesian(0d0,0d0,0d0)
                d_o_j = make_cartesian(0d0,0d0,0d0)
                sphere_radius = 15 * unit_mu
                lattice_sphere_radius = 2.75 / 2.0 * unit_mu
                cap_plane_point = make_cartesian(0d0, 0d0, 0d0)
                cap_plane_normal = make_cartesian(0d0, 0d0, 1d0)

                allocate(scene%hcp_sphere(1))
                scene%hcp_sphere(1) = lib_scene_generator_hcp_lattice_fill_sphere(d_o_j, &
                                                sphere_radius, lattice_sphere_radius, &
                                                cap_plane_point=cap_plane_point, &
                                                cap_plane_normal=cap_plane_normal)

                sphere_coordinates = list_filter(scene%hcp_sphere(1)%hcp_cuboid%hcp_lattice_coordiantes, &
                                                 scene%hcp_sphere(1)%inside_sphere)
                sphere_coordinates = scene%d_0_o + scene%hcp_sphere(1)%d_o_j + sphere_coordinates

                range_i(1) = lbound(sphere_coordinates, 1)
                range_i(2) = ubound(sphere_coordinates, 1)

                allocate(data_list(range_i(1):range_i(2), 3))

                do i = range_i(1), range_i(2)
                    buffer_coordinate = lib_math_get_matrix_rot_z(PI/2) * sphere_coordinates(i)
                    data_list(i, 1) = buffer_coordinate%x
                    data_list(i, 2) = buffer_coordinate%y
                    data_list(i, 3) = buffer_coordinate%z
!                    data_list(i, 1) = sphere_coordinates(i)%x
!                    data_list(i, 2) = sphere_coordinates(i)%y
!                    data_list(i, 3) = sphere_coordinates(i)%z
                end do

                u = 99
                filename = "temp/hcp_sphere.csv"
                open(unit=u, file=trim(filename), status='unknown')
                rv = write_csv(u, data_list)
                close(u)

                allocate (a_nm_t_matrix(size(sphere_coordinates)))
                allocate (b_nm_t_matrix(size(sphere_coordinates)))

                allocate (a_nm_ml_fmm(size(sphere_coordinates)))
                allocate (b_nm_ml_fmm(size(sphere_coordinates)))

                allocate (simulation_data%sphere_list(size(sphere_coordinates)))

                simulation_data%sphere_list(:)%sphere_parameter_index = 1

                simulation_data%sphere_list(:)%d_0_j = sphere_coordinates(:)

!                do i = 1, size(sphere_coordinates)
!
!                    simulation_data%sphere_list(i)%d_0_j = sphere_coordinates(i)
!                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(1))

                ! set 1
!                r_particle = 0.8 * unit_mu
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

!                if (n_range(2) .gt. 12) then
!                    n_range(2) = 12
!                    print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene: NOTE"
!                    print *, "  n_range is limited to ", n_range(2)
!                end if

                n_range(2) = 5
                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set evaluation points
                x_range = (/ -sphere_radius * 1.5d0, sphere_radius * 1.5d0 /)
                y_range = (/ -sphere_radius * 1.5d0, sphere_radius * 1.5d0 /)
                z_range = (/ -sphere_radius * 1.5d0, sphere_radius * 2d0 /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

!                allocate(simulation_data%evaluation_points(no_x_values * no_z_values))
!
!                x = 0
!                y = 0
!                z = 0
!                do i=1, no_x_values
!                    x = x_range(1) + (i-1) * step_size
!                    do ii=1, no_z_values
!                        z = z_range(1) + (ii-1) * step_size
!
!                        x_0%x = x
!                        x_0%y = y
!                        x_0%z = z
!
!                        simulation_data%evaluation_points(no_z_values * (i - 1) + ii)%coordinate = x_0
!                    end do
!                 end do

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
!                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)

                call fwig_table_init(4 * (n_range(2)+1), 3)
!                call fwig_temp_init(4 * (n_range(2)+1))
                call fwig_thread_temp_init(4 * (n_range(2)+1))

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export

!                call lib_mie_ms_ml_fmm_calculate_evaluation_points()

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/xz_")

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        call lib_field_plane_wave_get_field(simulation_data%illumination%plane_wave(1)%beam_parameter, &
                                                            x_0, &
                                                            field(1), field(2))

                        e_field_s(i,ii) = e_field_s(i,ii) + field(1)
                        h_field_s(i,ii) = h_field_s(i,ii) + field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/realTotal/xz_")

                ! xy plane

                range_i(1) = lbound(sphere_coordinates, 1)
                range_i(2) = ubound(sphere_coordinates, 1)

                z = sphere_coordinates(range_i(1))%z
                do i = range_i(1)+1, range_i(2)
                    if (z .lt. sphere_coordinates(i)%z) then
                        z = sphere_coordinates(i)%z
                    end if
                end do
                z = z + lattice_sphere_radius

                if (allocated(e_field_s)) deallocate(e_field_s)
                allocate(e_field_s(no_x_values, no_y_values))

                if (allocated(h_field_s)) deallocate(h_field_s)
                allocate(h_field_s(no_x_values, no_y_values))

                x = 0
                y = 0
!                z = -lattice_sphere_radius

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_y_values
                        y = y_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/xy_")

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_y_values
                        y = y_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        call lib_field_plane_wave_get_field(simulation_data%illumination%plane_wave(1)%beam_parameter, &
                                                            x_0, &
                                                            field(1), field(2))

                        e_field_s(i,ii) = e_field_s(i,ii) + field(1)
                        h_field_s(i,ii) = h_field_s(i,ii) + field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/realTotal/xy_")

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene:"
!                do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
!                    print *, "  i = ", i
!
!                    print *, "  a_nm"
!                    list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)
!
!                    do n = lbound(list_list_diff%item, 1), &
!                           ubound(list_list_diff%item, 1)
!                        print *, "  n = ", n
!                        do m = -n, n
!                            buffer = abs(list_list_diff%item(n)%item(m))
!                            if (buffer .gt. ground_truth_e) then
!                                print *, "    m: ", m , "difference: ", buffer, " : FAILED"
!                                print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
!                                print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
!                                rv = .false.
!                            else
!                                print *, "    m: ", m, ": OK"
!                            end if
!                        end do
!                    end do
!
!                    print *, "  b_nm"
!                    list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)
!
!                    do n = lbound(list_list_diff%item, 1), &
!                           ubound(list_list_diff%item, 1)
!                        print *, "  n = ", n
!                        do m = -n, n
!                            buffer = abs(list_list_diff%item(n)%item(m))
!                            if (buffer .gt. ground_truth_e) then
!                                print *, "    m: ", m , "difference: ", buffer, " : FAILED"
!                                print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
!                                print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
!                                rv = .false.
!                            else
!                                print *, "    m: ", m, ": OK"
!                            end if
!                        end do
!                    end do
!                end do

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene

            function test_lib_mie_ms_calculate_scattering_coeff_ab_nm_scene_eval() result(rv)
                use file_io
                use lib_mie_ms_ml_fmm_interface
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                double precision :: buffer
                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer, dimension(2) :: range_i
                integer :: ii
                integer :: n
                integer :: m

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix

                type(list_list_cmplx) :: list_list_diff


                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                type(lib_scene_opject_type) :: scene
                type(cartesian_coordinate_real_type) :: d_o_j
                type(cartesian_coordinate_real_type) :: cap_plane_point
                type(cartesian_coordinate_real_type) :: cap_plane_normal
                double precision :: sphere_radius
                double precision :: lattice_sphere_radius

                type(cartesian_coordinate_real_type), dimension(:), allocatable :: sphere_coordinates
                double precision, dimension(:,:), allocatable :: data_list

                integer :: u
                character(len=50) :: filename
                character(len=50), dimension(:), allocatable :: header


                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 550 * unit_nm

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, &
                                                                                        n_medium, &
                                                                                        e_field_0, &
                                                                                        plane_wave_g, &
                                                                                        plane_wave_k, &
                                                                                        plane_wave_d_0_i)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres

                scene%d_0_o = make_cartesian(0d0,0d0,0d0)
                d_o_j = make_cartesian(0d0,0d0,0d0)
                sphere_radius = 50 * unit_nm
                lattice_sphere_radius = 8 * unit_nm
                cap_plane_point = make_cartesian(0d0, 0d0, 0d0)
                cap_plane_normal = make_cartesian(0d0, 0d0, 1d0)

                allocate(scene%hcp_sphere(1))
                scene%hcp_sphere(1) = lib_scene_generator_hcp_lattice_fill_sphere(d_o_j, &
                                                sphere_radius, lattice_sphere_radius, &
                                                cap_plane_point=cap_plane_point, &
                                                cap_plane_normal=cap_plane_normal)

                sphere_coordinates = list_filter(scene%hcp_sphere(1)%hcp_cuboid%hcp_lattice_coordiantes, &
                                                 scene%hcp_sphere(1)%inside_sphere)
                sphere_coordinates = scene%d_0_o + scene%hcp_sphere(1)%d_o_j + sphere_coordinates

                range_i(1) = lbound(sphere_coordinates, 1)
                range_i(2) = ubound(sphere_coordinates, 1)

                allocate(data_list(range_i(1):range_i(2), 3))

                do i = range_i(1), range_i(2)
                    data_list(i, 1) = sphere_coordinates(i)%x
                    data_list(i, 2) = sphere_coordinates(i)%y
                    data_list(i, 3) = sphere_coordinates(i)%z
                end do

                u = 99
                filename = "temp/hcp_sphere.csv"
                open(unit=u, file=trim(filename), status='unknown')
                rv = write_csv(u, data_list)
                close(u)

                allocate (a_nm_t_matrix(size(sphere_coordinates)))
                allocate (b_nm_t_matrix(size(sphere_coordinates)))

                allocate (a_nm_ml_fmm(size(sphere_coordinates)))
                allocate (b_nm_ml_fmm(size(sphere_coordinates)))

                allocate (simulation_data%sphere_list(size(sphere_coordinates)))

                simulation_data%sphere_list(:)%sphere_parameter_index = 1

                simulation_data%sphere_list(:)%d_0_j = sphere_coordinates(:)

!                do i = 1, size(sphere_coordinates)
!
!                    simulation_data%sphere_list(i)%d_0_j = sphere_coordinates(i)
!                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(1))

                ! set 1
!                r_particle = 0.8 * unit_mu
                r_particle = 7.5 * unit_nm
                n_particle = dcmplx(2.5287_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

!                if (n_range(2) .gt. 12) then
!                    n_range(2) = 12
!                    print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene: NOTE"
!                    print *, "  n_range is limited to ", n_range(2)
!                end if

                n_range(2) = 6
                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set evaluation points
                x_range = (/ -sphere_radius * 1.5d0, sphere_radius * 1.5d0 /)
                z_range = (/ -sphere_radius * 1.5d0, sphere_radius * 2d0 /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))


                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
!                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .false., init_with_single_sphere = .true.)

                call fwig_table_init(4 * (n_range(2)+1), 3)
!                call fwig_temp_init(4 * (n_range(2)+1))
                call fwig_thread_temp_init(4 * (n_range(2)+1))

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_t_matrix(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_t_matrix(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene (T-Matrix): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

!                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
!                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = .true., init_with_single_sphere = .true.)

                call fwig_table_init(4 * (n_range(2)+1), 3)
!                call fwig_temp_init(4 * (n_range(2)+1))
                call fwig_thread_temp_init(4 * (n_range(2)+1))

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                ! evaluate and export



!                call lib_mie_ms_ml_fmm_calculate_evaluation_points()

!                allocate(e_field_s(no_x_values, no_z_values))
!                allocate(h_field_s(no_x_values, no_z_values))
!
!                x = 0
!                y = 0
!                z = 0
!
!                do i=1, no_x_values
!                    x = x_range(1) + (i-1) * step_size
!                    do ii=1, no_z_values
!                        z = z_range(1) + (ii-1) * step_size
!
!                        x_0%x = x
!                        x_0%y = y
!                        x_0%z = z
!
!                        field = lib_mie_ms_get_field(x_0)
!                        e_field_s(i,ii) = field(1)
!                        h_field_s(i,ii) = field(2)
!                    end do
!                 end do

!                do i=1, no_x_values
!                    x = x_range(1) + (i-1) * step_size
!                    do ii=1, no_z_values
!                        z = z_range(1) + (ii-1) * step_size
!
!                        x_0%x = x
!                        x_0%y = y
!                        x_0%z = z
!
!                        field = lib_mie_ms_get_field(x_0)
!                        simulation_data%evaluation_points(no_z_values * (i - 1) + ii)%coordinate = x_0
!                        e_field_s(i,ii) = field(1)
!                        h_field_s(i,ii) = field(2)
!                    end do
!                 end do

!                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")
!
!                x = 0
!                y = 0
!                z = 0
!
!                do i=1, no_x_values
!                    x = x_range(1) + (i-1) * step_size
!                    do ii=1, no_z_values
!                        z = z_range(1) + (ii-1) * step_size
!
!                        x_0%x = x
!                        x_0%y = y
!                        x_0%z = z
!
!                        call lib_field_plane_wave_get_field(simulation_data%illumination%plane_wave(1)%beam_parameter, &
!                                                            x_0, &
!                                                            field(1), field(2))
!
!                        e_field_s(i,ii) = e_field_s(i,ii) + field(1)
!                        h_field_s(i,ii) = h_field_s(i,ii) + field(2)
!                    end do
!                 end do
!
!                rv = lib_field_export(e_field_s, h_field_s, "temp/realTotal/")

                !      1   , 2, 3,
                ! sphere_no, n, m,
                !     4           ,     5            ,     6           ,     7            ,
                ! a_nm_real_ml_fmm, a_nm_cmplx_ml_fmm, b_nm_real_ml_fmm, b_nm_cmplx_ml_fmm,
                !     8             ,       9            ,       10          ,       11           ,
                ! a_nm_real_t_matrix, a_nm_cmplx_t_matrix, b_nm_real_t_matrix, b_nm_cmplx_t_matrix,
                !     12        ,       13       ,     14        ,    15
                ! diff_a_nm_real, diff_a_nm_cmplx, diff_b_nm_real, diff_b_nm_cmplx
                if (allocated(header)) deallocate(header)
                allocate(header(15))
                header(1) = "sphere_no"
                header(2) = "n"
                header(3) = "m"
                header(4) = "a_nm_real_ml_fmm"
                header(5) = "a_nm_cmplx_ml_fmm"
                header(6) = "b_nm_real_ml_fmm"
                header(7) = "b_nm_cmplx_ml_fmm"
                header(8) =  "a_nm_real_t_matrix"
                header(9) =  "a_nm_cmplx_t_matrix"
                header(10) = "b_nm_real_t_matrix"
                header(11) = "b_nm_cmplx_t_matrix"
                header(12) = "diff_a_nm_real_ml_fmm"
                header(13) = "diff_a_nm_cmplx_ml_fmm"
                header(14) = "diff_b_nm_real_ml_fmm"
                header(15) = "diff_b_nm_cmplx_ml_fmm"

                range_i(1) = lbound(sphere_coordinates, 1)
                range_i(2) = ubound(sphere_coordinates, 1)
                i = (range_i(2) - range_i(1) + 1) * ( (n_range(2) + 1)**2 - n_range(1)**2 )

                if (allocated(data_list)) deallocate(data_list)
                allocate(data_list(i, size(header)))

                rv = .true.
                print *, "test_lib_mie_ms_calculate_scattering_coefficients_ab_nm_scene:"
                ii = 0
                do i=lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                    print *, "  i = ", i

                    print *, "  a_nm"
                    list_list_diff = a_nm_ml_fmm(i) - a_nm_t_matrix(i)

                    do n = lbound(list_list_diff%item, 1), &
                           ubound(list_list_diff%item, 1)
                        print *, "  n = ", n
                        do m = -n, n

                            ii = ii + 1
                            data_list(ii, 1) = i  ! sphere no
                            data_list(ii, 2) = n  ! n
                            data_list(ii, 3) = m  ! m

                            data_list(ii, 4) =   real(a_nm_ml_fmm(i)%item(n)%item(m))
                            data_list(ii, 5) =  aimag(a_nm_ml_fmm(i)%item(n)%item(m))
                            data_list(ii, 6) =   real(b_nm_ml_fmm(i)%item(n)%item(m))
                            data_list(ii, 7) =  aimag(b_nm_ml_fmm(i)%item(n)%item(m))

                            data_list(ii, 8) =   real(a_nm_t_matrix(i)%item(n)%item(m))
                            data_list(ii, 9) =  aimag(a_nm_t_matrix(i)%item(n)%item(m))
                            data_list(ii, 10) =  real(b_nm_t_matrix(i)%item(n)%item(m))
                            data_list(ii, 11) = aimag(b_nm_t_matrix(i)%item(n)%item(m))

                            data_list(ii, 12) =  real(list_list_diff%item(n)%item(m))
                            data_list(ii, 13) = aimag(list_list_diff%item(n)%item(m))
                            data_list(ii, 14) =  real(list_list_diff%item(n)%item(m))
                            data_list(ii, 15) = aimag(list_list_diff%item(n)%item(m))

                            buffer = real(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (real): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, ": (real) OK"
                            end if

                            buffer = aimag(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (imag): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", a_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", a_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, ": (imag) OK"
                            end if
                        end do
                    end do

                    print *, "  b_nm"
                    list_list_diff = b_nm_ml_fmm(i) - b_nm_t_matrix(i)

                    do n = lbound(list_list_diff%item, 1), &
                           ubound(list_list_diff%item, 1)
                        print *, "  n = ", n
                        do m = -n, n
                            buffer = real(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (real): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, ": (real) OK"
                            end if
                        end do

                        do m = -n, n
                            buffer = aimag(list_list_diff%item(n)%item(m))
                            if (buffer .gt. ground_truth_e) then
                                print *, "    m: ", m , "difference (imag): ", buffer, " : FAILED"
                                print *, "         ML-FMM: ", b_nm_ml_fmm(i)%item(n)%item(m)
                                print *, "       T-Matrix: ", b_nm_t_matrix(i)%item(n)%item(m)
                                rv = .false.
                            else
                                print *, "    m: ", m, ": (imag) OK"
                            end if
                        end do
                    end do
                end do

                u = 99
                filename = "temp/coefficient_error.csv"
                open(unit=u, file=trim(filename), status='unknown')
                rv = write_csv(u, data_list, header)
                close(u)

                call lib_mie_ms_data_container_destructor()

            end function test_lib_mie_ms_calculate_scattering_coeff_ab_nm_scene_eval

        end function

        subroutine lib_mie_ms_benchmark()
            implicit none
            ! dummy

            ! auxiliary

            call benchmark_sphere_scene

            contains

            subroutine benchmark_sphere_scene()
                use file_io
                implicit none
                ! auxiliry
                integer :: i

                double precision, dimension(5) :: res

                double precision :: sphere_radius_start
                double precision :: sphere_radius_stop

                double precision :: sphere_radius
                double precision :: lattice_sphere_radius
                double precision :: particle_radius

                double precision :: step_size
                integer :: no_of_steps

                logical :: use_ml_fmm
                logical :: plot_result

                character(len=50) :: filename
                character(len=150), dimension(5) :: header
                integer :: u
                logical :: rv_file
                double precision, dimension(:,:), allocatable :: csv_data


                use_ml_fmm = .true.
                plot_result = .true.

                sphere_radius_start = 20 * unit_nm
                sphere_radius_stop = 64 * unit_nm

                lattice_sphere_radius = 7.6 * unit_nm
                particle_radius = 7.5 * unit_nm

                no_of_steps = 12

                allocate(csv_data(no_of_steps, 5))
                csv_data = 0

                step_size = (sphere_radius_stop - sphere_radius_start) / no_of_steps

                do i = 1, no_of_steps

                    sphere_radius = sphere_radius_start + step_size * (i - 1)

                    res = calculate_sphere_scene(sphere_radius, lattice_sphere_radius, particle_radius, &
                                                 plot_result, use_ml_fmm)

                    csv_data(i, :) = res

                    u = 99
                    filename = "temp/hcp_sphere_benchmark_ml_fmm.csv"
                    header(1) = "no_spheres"
                    header(2) = "calc: CPU-Time"
                    header(3) = "calc: WALL-Time"
                    header(4) = "plot: CPU-Time"
                    header(5) = "plot: WALL-Time"

                    open(unit=u, file=trim(filename), status='unknown')
                    rv_file = write_csv(u, csv_data, header)
                    close(u)

                    call lib_mie_ms_data_container_destructor
                end do

                use_ml_fmm = .false.
                do i = 1, no_of_steps

                    sphere_radius = sphere_radius_start + step_size * (i - 1)

                    res = calculate_sphere_scene(sphere_radius, lattice_sphere_radius, particle_radius, &
                                                 plot_result, use_ml_fmm)

                    csv_data(i, :) = res

                    u = 99
                    filename = "temp/hcp_sphere_benchmark_t_matrix.csv"
                    header(1) = "no_spheres"
                    header(2) = "calc: CPU-Time"
                    header(3) = "calc: WALL-Time"
                    header(4) = "plot: CPU-Time"
                    header(5) = "plot: WALL-Time"

                    open(unit=u, file=trim(filename), status='unknown')
                    rv_file = write_csv(u, csv_data, header)
                    close(u)

                    call lib_mie_ms_data_container_destructor
                end do

            end subroutine

            !
            !
            !
            ! rv(1): numer of spheres
            ! rv(2): interactive coefficients: CPU-time
            ! rv(3): interactive coefficients: WALL-time
            !
            !
            function calculate_sphere_scene(sphere_radius, lattice_sphere_radius, r_particle, &
                                            plot_results, use_ml_fmm) result(rv)
                use file_io
                implicit none
                ! dummy
                double precision, intent(in) :: sphere_radius
                double precision, intent(in) :: lattice_sphere_radius
                double precision, intent(in) :: r_particle
                logical, intent(in) :: plot_results
                logical, intent(in) :: use_ml_fmm

                double precision, dimension(5) :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1

                ! auxiliaray
                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer, dimension(2) :: range_i
                integer :: ii

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_ml_fmm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_ml_fmm

                type(list_list_cmplx), dimension(:), allocatable :: a_nm_t_matrix
                type(list_list_cmplx), dimension(:), allocatable :: b_nm_t_matrix


                ! plot
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size
                integer :: no_x_values
                integer :: no_z_values

                type(lib_scene_opject_type) :: scene
                type(cartesian_coordinate_real_type) :: d_o_j
                type(cartesian_coordinate_real_type) :: cap_plane_point
                type(cartesian_coordinate_real_type) :: cap_plane_normal

                type(cartesian_coordinate_real_type), dimension(:), allocatable :: sphere_coordinates
                double precision, dimension(:,:), allocatable :: data_list

                integer :: u
                character(len=50) :: filename
                logical :: rv_file


                ! CPU-time
                real :: test_start_sub, test_finish_sub
                ! WALL-time
                INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

                rv(:) = 0


                simulation_data%spherical_harmonics%z_selector_incident_wave = 1
                simulation_data%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation_data%spherical_harmonics%z_selector_translation_gt_r = 3
                simulation_data%spherical_harmonics%z_selector_translation_le_r = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 550 * unit_nm

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation_data%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0 * unit_mu
                buffer_car%y = 0
                buffer_car%z = -10 * unit_mu
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation_data%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, &
                                                                                        n_medium, &
                                                                                        e_field_0, &
                                                                                        plane_wave_g, &
                                                                                        plane_wave_k, &
                                                                                        plane_wave_d_0_i)

                simulation_data%illumination%plane_wave(1)%beam_parameter%polarisation = &
                        lib_field_polarisation_jones_vector_get_linear_rot(PI / 4d0)

                simulation_data%refractive_index_medium = n_medium

                ! set spheres

                scene%d_0_o = make_cartesian(0d0,0d0,0d0)
                d_o_j = make_cartesian(0d0,0d0,0d0)
                cap_plane_point = make_cartesian(0d0, 0d0, 0d0)
                cap_plane_normal = make_cartesian(0d0, 0d0, 1d0)

                allocate(scene%hcp_sphere(1))
                scene%hcp_sphere(1) = lib_scene_generator_hcp_lattice_fill_sphere(d_o_j, &
                                                sphere_radius, lattice_sphere_radius, &
                                                cap_plane_point=cap_plane_point, &
                                                cap_plane_normal=cap_plane_normal)

                sphere_coordinates = list_filter(scene%hcp_sphere(1)%hcp_cuboid%hcp_lattice_coordiantes, &
                                                 scene%hcp_sphere(1)%inside_sphere)
                sphere_coordinates = scene%d_0_o + scene%hcp_sphere(1)%d_o_j + sphere_coordinates

                rv(1) = real(size(sphere_coordinates))

                range_i(1) = lbound(sphere_coordinates, 1)
                range_i(2) = ubound(sphere_coordinates, 1)

                allocate(data_list(range_i(1):range_i(2), 3))

                do i = range_i(1), range_i(2)
                    data_list(i, 1) = sphere_coordinates(i)%x
                    data_list(i, 2) = sphere_coordinates(i)%y
                    data_list(i, 3) = sphere_coordinates(i)%z
                end do

                u = 99
                filename = "temp/hcp_sphere.csv"
                open(unit=u, file=trim(filename), status='unknown')
                rv_file = write_csv(u, data_list)
                close(u)

                allocate (a_nm_t_matrix(size(sphere_coordinates)))
                allocate (b_nm_t_matrix(size(sphere_coordinates)))

                allocate (a_nm_ml_fmm(size(sphere_coordinates)))
                allocate (b_nm_ml_fmm(size(sphere_coordinates)))

                allocate (simulation_data%sphere_list(size(sphere_coordinates)))

                simulation_data%sphere_list(:)%sphere_parameter_index = 1

                simulation_data%sphere_list(:)%d_0_j = sphere_coordinates(:)

!                do i = 1, size(sphere_coordinates)
!
!                    simulation_data%sphere_list(i)%d_0_j = sphere_coordinates(i)
!                end do

                ! set sphere parameter
                allocate(simulation_data%sphere_parameter_list(1))

                ! set 1
!                r_particle = 1 * unit_mu
!                r_particle = 7.5 * unit_nm
                n_particle = dcmplx(2.5287_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

!                if (n_range(2) .gt. 8) then
!                    n_range(2) = 8
!                    print *, "calculate_sphere_scene: NOTE"
!                    print *, "  n_range is limited to 8"
!                end if

                simulation_data%spherical_harmonics%n_range = n_range

                simulation_data%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set evaluation points
                x_range = (/ -sphere_radius * 1.5d0, sphere_radius * 1.5d0 /)
                z_range = (/ -sphere_radius * 1.5d0, sphere_radius * 2d0 /)
                step_size = (x_range(2) - x_range(1)) / 200

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

!                allocate(simulation_data%evaluation_points(no_x_values * no_z_values))
!
!                x = 0
!                y = 0
!                z = 0
!                do i=1, no_x_values
!                    x = x_range(1) + (i-1) * step_size
!                    do ii=1, no_z_values
!                        z = z_range(1) + (ii-1) * step_size
!
!                        x_0%x = x
!                        x_0%y = y
!                        x_0%z = z
!
!                        simulation_data%evaluation_points(no_z_values * (i - 1) + ii)%coordinate = x_0
!                    end do
!                 end do

                call system_clock(test_count_start_sub, test_count_rate_sub)
                call cpu_time(test_start_sub)

                n_range(2) = int(ceiling(real(n_range(2)) * 1.2))
                simulation_data%spherical_harmonics%n_range = n_range

                call lib_mie_ms_constructor(n_range, use_ml_fmm = use_ml_fmm, init_with_single_sphere = .true.)

!                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(8, 2)
                call lib_mie_ms_calculate_scattering_coefficients_ab_nm()
                a_nm_ml_fmm(:) = simulation_data%sphere_list(:)%a_nm
                b_nm_ml_fmm(:) = simulation_data%sphere_list(:)%b_nm

                call cpu_time(test_finish_sub)
                call system_clock(test_count_finish_sub, test_count_rate_sub)


                print *, ""
                print *, "calculate_sphere_scene (ML-FMM): "
                print '("  CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
                print '("  WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                               / real(test_count_rate_sub)
                print *, ""

                rv(2) = test_finish_sub-test_start_sub
                rv(3) = (test_count_finish_sub-test_count_start_sub) / real(test_count_rate_sub)

                ! evaluate and export
                if (plot_results) then
!                call lib_mie_ms_ml_fmm_calculate_evaluation_points()
                    allocate(e_field_s(no_x_values, no_z_values))
                    allocate(h_field_s(no_x_values, no_z_values))

                    x = 0
                    y = 0
                    z = 0

                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            x_0%x = x
                            x_0%y = y
                            x_0%z = z

                            field = lib_mie_ms_get_field(x_0)
                            e_field_s(i,ii) = field(1)
                            h_field_s(i,ii) = field(2)
                        end do
                     end do

    !                do i=1, no_x_values
    !                    x = x_range(1) + (i-1) * step_size
    !                    do ii=1, no_z_values
    !                        z = z_range(1) + (ii-1) * step_size
    !
    !                        x_0%x = x
    !                        x_0%y = y
    !                        x_0%z = z
    !
    !                        field = lib_mie_ms_get_field(x_0)
    !                        simulation_data%evaluation_points(no_z_values * (i - 1) + ii)%coordinate = x_0
    !                        e_field_s(i,ii) = field(1)
    !                        h_field_s(i,ii) = field(2)
    !                    end do
    !                 end do

                    rv_file = lib_field_export(e_field_s, h_field_s, "temp/real/")

                    x = 0
                    y = 0
                    z = 0

                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            x_0%x = x
                            x_0%y = y
                            x_0%z = z

                            call lib_field_plane_wave_get_field(simulation_data%illumination%plane_wave(1)%beam_parameter, &
                                                                x_0, &
                                                                field(1), field(2))

                            e_field_s(i,ii) = e_field_s(i,ii) + field(1)
                            h_field_s(i,ii) = h_field_s(i,ii) + field(2)
                        end do
                     end do

                    rv_file = lib_field_export(e_field_s, h_field_s, "temp/realTotal/")
                end if

            end function

        end subroutine
end module lib_mie_multi_sphere
