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

module lib_scene_type
    use libmath
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
    !       K_o |/
    !           --> x
    !
    ! K_o: object coordinate system
    ! K_j: sphere coordinate system
    !   o: j-th sphere
    !
    type lib_scene_object_sphere_type
        type(cartesian_coordinate_real_type) :: d_o_j ! [m]
        integer :: sphere_parameter_index
        double precision :: radius  ! [m]
    end type


    !
    !         . . . . z . . . .
    !          . . . .^  . . .
    !         . . K_j | . . . .
    !          . . . .-->.x. .
    !         . . . .^. . . . .
    !          . . ./. . . . .
    !         . z ./d_o_j . . .
    !           ^ /
    !       K_o |/
    !           --> x
    !
    !   K_o: object coordinate system
    !   o: spheres at the hcp_lattice_coordiantes
    !
    type lib_scene_object_hcp_cuboid_type
        type(cartesian_coordinate_real_type) :: d_o_j
        double precision :: sphere_radius
        type(cartesian_coordinate_real_type), dimension(:), allocatable :: hcp_lattice_coordiantes
    end type

    !
    !               . z .
    !            . . .^  . .
    !           . K_j | . . .
    !          . . . .-->.x. .
    !           . . .^. . . .
    !            . ./. . . .
    !           z ./d_o_j
    !           ^ /
    !       K_o |/
    !           --> x
    !
    !   K_o: object coordinate system
    !   o: spheres at the hcp_lattice_coordiantes
    !
    type lib_scene_object_hcp_sphere_type
        type(cartesian_coordinate_real_type) :: d_o_j
        type(lib_scene_object_hcp_cuboid_type) :: hcp_cuboid
        logical, dimension(:), allocatable :: inside_sphere
    end type

    !                 z
    !                 ^
    !             K_o |
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_o
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_o: object coordinate system
    !
    type lib_scene_opject_type
        type(cartesian_coordinate_real_type) :: d_0_o ! [m]
        type(lib_scene_object_sphere_type), dimension(:), allocatable :: sphere
        type(lib_scene_object_hcp_cuboid_type), dimension(:), allocatable :: hcp_cuboid
        type(lib_scene_object_hcp_sphere_type), dimension(:), allocatable :: hcp_sphere
    end type

end module lib_scene_type
