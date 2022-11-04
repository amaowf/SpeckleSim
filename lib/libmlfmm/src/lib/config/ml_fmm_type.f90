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

! LIB: Mulitlevel Fast Multipole Method - type definitions
!
module ml_fmm_type
    use libmath
    use libtree
    implicit none

    integer(kind=1), parameter, public :: LIB_ML_FMM_COEFFICIENT_KIND = 8

    integer(kind=1), parameter, public :: LIB_ML_FMM_COEFFICIENT_TYPE_C = 1
    integer(kind=1), parameter, public :: LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE = 2
    integer(kind=1), parameter, public :: LIB_ML_FMM_COEFFICIENT_TYPE_D = 3

    type lib_ml_fmm_data
        type(lib_tree_data_element), dimension(:), allocatable :: X
        type(lib_tree_data_element), dimension(:), allocatable :: Y
        type(lib_tree_data_element), dimension(:), allocatable :: XY
    end type

    type lib_ml_fmm_v
        type(list_list_cmplx) :: a_nm
        type(list_list_cmplx) :: b_nm
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: r
        complex(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: c
    end type

    ! e.g. A, B, C or D coefficient
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    !
    !                   -----
    !                   | 0 |                 l = 0
    !                   -----
    !         -----               -----
    !         | 0 |               | 1 |       l = 1
    !         -----               -----
    !    -----     -----     -----     -----
    !    | 0 |     | 1 |     | 2 |     | 3 |  l = 2
    !    -----     -----     -----     -----
    !      ^
    !
    ! coefficients of the box uindex(n,l)
    !   e.g. (0,2): C_i: type(lib_ml_coefficient)
    !
    !
    type lib_ml_fmm_coefficient
        type(list_list_cmplx) :: a_nm
        type(list_list_cmplx) :: b_nm
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: r
        complex(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: c
        type(lib_tree_universal_index) :: uindex
    end type

    type lib_ml_fmm_hashed_coeffcient_index
        integer(kind=UINDEX_BYTES) :: array_position
        integer(kind=2) :: number_of_hash_runs
    end type

    ! hierarchy of the data structure for the ML FMM
    ! Reference: Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, chapter 2. Data hierarchies
    !
    ! This type represents a dataset of a hierarchy level.
    !
    ! Example
    ! ----
    !
    !                            -----
    !                            | 0 |                 l = 0
    !                            -----
    !                  -----      XY       -----
    !                  | 0 |               | 1 |       l = 1
    !                  -----            v  -----
    !             -----  Y  -----     -----  X  -----
    !             | 0 |     | 1 |     | 2 |     | 3 |  l = 2    <---
    !             -----     -----     -----     -----
    !     type:     Y         Y         X         -
    !
    !     Variable                    Value
    !    -----------------------------------
    !     coefficient_type              C
    !     number_of_boxes               3         <-- these boxes are not empty
    !     type                       [X,Y,Y]      <-- this level contains boxes of both hierarchies (X and Y)
    !     size(coefficient_list)        3
    !
    !
    !   coefficient_list at level l=2
    !       -------------------
    !       |  2  |  0  |  1  |
    !       -------------------
    !        <-X-> <-Y-> <-Y->
    !
    !   uindex_list_X l=2
    !       -------
    !       |  2  |
    !       -------
    !   uindex_list_Y l=2
    !      -------------
    !      |  0  |  1  |
    !      -------------
    !   uindex_list_XY l=2
    !       -------
    !       |-----|  (not allocated
    !       -------
    !
    !   Get the coefficients of box(2,2)
    !
    !   Case: is_hashed .eq. .false.
    !   -----
    !       coefficient_list_index
    !          n|0             v       3|
    !           -------------------------
    !           |  2  |  3  |  1  | -1  |
    !           -------------------------
    !
    !       coefficient_list
    !           |1               3|
    !           -------------------
    !           |  *  |     |     |
    !           -------------------
    type lib_ml_fmm_hierarchy
        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: coefficient_list
        ! [X, Y, XY] hierarchy type
        integer(kind=1), dimension(:), allocatable :: hierarchy_type
        ! [C, D^~, D]
        integer(kind=1), dimension(:), allocatable :: coefficient_type
        type(lib_ml_fmm_hashed_coeffcient_index), dimension(:), allocatable :: hashed_coefficient_list_index
        ! is_hashed .eqv. .true.: list entry correspondence with the box index n
        ! is_hashed .eqv. .false: list entry correspondence with the coefficient_list_index index
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: coefficient_list_index
        ! true: access coefficient with hashed uindex%n, false: access coefficient with uindex%n
        logical :: is_hashed
        ! number of not empty boxes
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=2) :: maximum_number_of_hash_runs
    end type
end module ml_fmm_type
