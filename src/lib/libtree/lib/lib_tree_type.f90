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

! spatial dimension, value = [2,3]
#define _FMM_DIMENSION_ 3

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
!
! Constraint
! ----
!          |  _FMM_DIMENSION_  |
!   value  |    2     |   3    |
!   -----------------------------
!   single | [4,8]    | [8]    |
!   double | [8,16]   | [8,16] |
!
#define _UINDEX_BYTES_ 8

module lib_tree_type
    use libmath
    implicit none

    ! parameter
    integer(kind=1), public, parameter :: TREE_DIMENSIONS = _FMM_DIMENSION_ ! dimensions

    integer(kind=1), public, parameter :: CORRESPONDENCE_VECTOR_KIND = 4    ! limited by the total number of elements -> 2**(8*4) = 4,294,967,296 elements

#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
    integer(kind=1), public, parameter :: COORDINATE_BINARY_BYTES = 8
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
    integer(kind=1), public, parameter :: COORDINATE_BINARY_BYTES = 4
#endif

    integer(kind=1), public, parameter :: UINDEX_BYTES = _UINDEX_BYTES_

    ! data element association: source(X), target(Y) or both(XY)
    ! Reference: Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, chapter 2. Data hierarchies
    integer(kind=1), public, parameter :: HIERARCHY_X = 1
    integer(kind=1), public, parameter :: HIERARCHY_Y = 2
    integer(kind=1), public, parameter :: HIERARCHY_XY = 3

    integer(kind=1), parameter, public :: IGNORE_ENTRY = -1
    ! ~ parameter ~

    ! type definitions
    type lib_tree_spatial_point
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, dimension(TREE_DIMENSIONS) :: x
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, dimension(TREE_DIMENSIONS) :: x
#endif
    end type lib_tree_spatial_point

    type lib_tree_universal_index
        integer(kind=UINDEX_BYTES) :: n
!        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), &
!            dimension(TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
!            :: interleaved_bits
        integer(kind=1) :: l
    end type lib_tree_universal_index

    type lib_tree_correspondece_vector_element
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: data_element_number
        integer(kind=2) :: number_of_hash_runs
    end type lib_tree_correspondece_vector_element
    ! ~ type definitions ~

    ! --- type definition ---
    type lib_tree_data_element
        type(lib_tree_spatial_point) :: point_x
        type(lib_tree_universal_index) :: uindex
        integer(kind=1) :: element_type
        integer(kind=1) :: hierarchy
        double precision :: expansion
    end type lib_tree_data_element

    type lib_tree_traversal_data
        ! --- input ---
        type(cartesian_coordinate_real_type) :: ray_origin
        type(cartesian_coordinate_real_type) :: ray_direction

        ! --- optional ---
        logical :: is_terminal = .false.

        ! --- payload ---
        integer, dimension(:), allocatable :: i
        double precision, dimension(:), allocatable :: d
        ! dimension(1:N)
        ! N: number of traversed boxes
        type(lib_tree_universal_index), dimension(:), allocatable :: traversed_box_list
        integer :: index = -1
    end type

end module lib_tree_type
