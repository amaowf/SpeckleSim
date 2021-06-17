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

module lib_mie_type
    use libmath
    use lib_mie_illumination
    implicit none

    public

    !                 z
    !                 ^
    !             K_j |
    !                 o--> x
    !                ^
    !               /
    !           z  /d_0_j
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_j: sphere coordinate system
    !   o: j-th sphere
    !
    type lib_mie_sphere_type
        type(cartesian_coordinate_real_type) :: d_0_j ! [m]
        integer :: sphere_parameter_index
        type(list_list_cmplx) :: a_nm
        type(list_list_cmplx) :: b_nm
    end type lib_mie_sphere_type

    ! sphere parameter
    type lib_mie_sphere_parameter_type
        double complex :: refractive_index
        double precision :: radius  ! [m]
        double precision :: size_parameter ! = k * x = 2 Pi * n_medium / lambda * x
        integer, dimension(2) :: n_range
        type(list_cmplx) :: a_n
        type(list_cmplx) :: b_n
    end type lib_mie_sphere_parameter_type

    ! Depending on the initial definition of the electromagnetic wave,
    ! the spherical harmonics have to be calculated differently.
    !
    ! case: exp(i (k x + omega t) // actual case
    !   z_selector_incident_wave = 1
    !   z_selector_scatterd_wave = 3
    !   z_selector_translation_le_r = 3
    !   z_selector_translation_gt_r = 1
    !
    ! case: exp(-i (k x + omega t)
    !   z_selector_incident_wave = 2
    !   z_selector_scatterd_wave = 4
    !   z_selector_translation = 4
    !
    ! n_range
    !   Defines the minimum and maximum degree for the entire simulation.
    !
    type lib_mie_vector_spherical_harmonics_type
        integer(kind=1) :: z_selector_incident_wave
        integer(kind=1) :: z_selector_scatterd_wave
        integer(kind=1) :: z_selector_translation_le_r
        integer(kind=1) :: z_selector_translation_gt_r
        integer, dimension(2) :: n_range
    end type lib_mie_vector_spherical_harmonics_type

    type lib_mie_evaluation_point_type
        type(cartesian_coordinate_real_type) :: coordinate
        type(cartesian_coordinate_cmplx_type) :: e_field
        type(cartesian_coordinate_cmplx_type) :: h_field
    end type

    ! simulation parameter
    type lib_mie_simulation_parameter_type
        double precision :: refractive_index_medium
        type(lib_mie_illumination_parameter) :: illumination
        type(lib_mie_sphere_type), dimension(:), allocatable :: sphere_list
        type(lib_mie_sphere_parameter_type), dimension(:), allocatable :: sphere_parameter_list
        type(lib_mie_vector_spherical_harmonics_type) :: spherical_harmonics
        type(lib_mie_evaluation_point_type), dimension(:), allocatable :: evaluation_points
    end type lib_mie_simulation_parameter_type

end module lib_mie_type
