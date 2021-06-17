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


program main
    !$  use omp_lib
    use file_io
    use libmath
    use lib_data_types
    use libtree

    use lib_hash_function
    use lib_sort

    use lib_scene_generator

    use lib_field

    use libmlfmm

    use lib_field_polarisation
    use lib_field_plane_wave
    use lib_field_gaussian_beam

    !use lib_mie_single_sphere
    !use lib_mie_ss_helper_functions

    !use lib_mie_ms_solver_interface_helper_functions
    !use lib_mie_multi_sphere
    !use lib_mie_ms_solver_interface

    use light_scattering
    implicit none

    integer :: error_counter

    integer(kind=4) :: l
    integer(kind=4) :: m
    double precision :: x
!    double precision :: erg

    ! CPU-time
    real :: test_start, test_finish
    ! WALL-time
    INTEGER :: test_count_start, test_count_finish, test_count_rate


    l = 3
    m = 1
    x = 0.01_8

!    polarization: parallel (1) perpendicular (2)
!    size parameter: x
!    index of refraction: real,imaginary (+ for absorption)
!    npnts: number of grid points (npnts by npnts)
!
!    e.g.:
!        ip = 1        ! polarisation
!        x = 2*pi * 10 ! sphere radius: a=10um;  wave length: 1um; a / wave length = 10
!        cmr = 1.33    ! water: 1.33
!        cmi = 0       ! water: 0
!        npnts = 10
!   call S2(1, 2*3.14159265358979*10, 1.33, 0.0, 10)
!    call S2(1, 20.0, 1.50, 0.0, 10)
    !call S2(1, 20.0, 1.28, 1.37, 10)

    !call test_file_io

!    call OMP_set_num_threads(16)
!    call OMP_set_nested(.true.)
!    call OMP_set_dynamic(.true.)

    call system_clock(test_count_start, test_count_rate)
    call cpu_time(test_start)

   error_counter = 0
!   error_counter = error_counter + lib_scene_generator_test_functions()
    error_counter = error_counter + test_lib_math()
!    error_counter = error_counter + lib_sort_test_functions()
!    error_counter = error_counter + lib_test_hash_function()
!    error_counter = error_counter + lib_tree_hf_test_functions()
!    error_counter = error_counter + lib_tree_test_functions()
!    error_counter = error_counter + lib_ml_fmm_type_operator_test_functions()
!    error_counter = error_counter + lib_ml_fmm_hf_test_functions()
!    error_counter = error_counter + lib_ml_fmm_test_functions()
!    error_counter = error_counter + lib_field_test_functions()
!    error_counter = error_counter + lib_field_polarisation_operator_test_functions()
!    error_counter = error_counter + lib_field_plane_wave_test_functions()
!    error_counter = error_counter + lib_field_gaussian_beam_test_functions()
!    error_counter = error_counter + lib_mie_ss_helper_functions_test_functions()
!    error_counter = error_counter + lib_mie_single_sphere_test_functions()
!    error_counter = error_counter + lib_mie_ms_solver_interface_hf_helper_functions()
!    error_counter = error_counter + lib_mie_ms_solver_interface_test_functions()
!    error_counter = error_counter + lib_mie_multi_sphere_test_functions()

!    call lib_tree_benchmark
!    call lib_mie_ss_helper_functions_benchmark
!    call lib_mie_ms_benchmark

    call cpu_time(test_finish)
    call system_clock(test_count_finish, test_count_rate)

    print *, ""
    print *, "-------------MAIN------------------"
    print '("    CPU-Time = ",f10.3," seconds.")',test_finish-test_start
    print '("    WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
    print *, ""
    if (error_counter == 0) then
        print *, "All tests: OK"
    else
        print *, error_counter,"test(s) FAILED"
    end if
    print *, "-----------------------------------"

!    call lib_tree_benchmark()
!     call lib_tree_hf_benchmark()

    call lib_tree_destructor()
!    call lib_tree_hf_destructor()

end program main
