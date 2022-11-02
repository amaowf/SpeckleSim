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

module lib_math_type
    implicit none

    integer(kind=1), parameter :: lib_math_type_kind = 8

    ! cartesian coordinates
    type cartesian_coordinate_real_type
        real(kind=lib_math_type_kind) :: x
        real(kind=lib_math_type_kind) :: y
        real(kind=lib_math_type_kind) :: z
    end type cartesian_coordinate_real_type

    type list_cartesian_coordinate_real_type
        type(cartesian_coordinate_real_type), dimension(:), allocatable :: coordinate
    end type list_cartesian_coordinate_real_type

    type cartesian_coordinate_cmplx_type
        complex(kind=lib_math_type_kind) :: x
        complex(kind=lib_math_type_kind) :: y
        complex(kind=lib_math_type_kind) :: z
    end type cartesian_coordinate_cmplx_type

    type list_cartesian_coordinate_cmplx_type
        type(cartesian_coordinate_cmplx_type), dimension(:), allocatable :: coordinate
    end type list_cartesian_coordinate_cmplx_type

    type cartresian_coordinate_rot_matrix_type
        real(kind=lib_math_type_kind) :: r_11
        real(kind=lib_math_type_kind) :: r_12
        real(kind=lib_math_type_kind) :: r_13
        real(kind=lib_math_type_kind) :: r_21
        real(kind=lib_math_type_kind) :: r_22
        real(kind=lib_math_type_kind) :: r_23
        real(kind=lib_math_type_kind) :: r_31
        real(kind=lib_math_type_kind) :: r_32
        real(kind=lib_math_type_kind) :: r_33
    end type

    ! spherical coordinates
    type spherical_coordinate_real_type
        real(kind=lib_math_type_kind) :: rho
        real(kind=lib_math_type_kind) :: theta
        real(kind=lib_math_type_kind) :: phi
    end type spherical_coordinate_real_type

    type list_spherical_coordinate_real_type
        type(spherical_coordinate_real_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_real_type

    type spherical_coordinate_cmplx_type
        complex(kind=lib_math_type_kind) :: rho
        complex(kind=lib_math_type_kind) :: theta
        complex(kind=lib_math_type_kind) :: phi
    end type spherical_coordinate_cmplx_type

    type list_spherical_coordinate_cmplx_type
        type(spherical_coordinate_cmplx_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_cmplx_type

    ! lists
    type list_logical
        logical, dimension(:), allocatable :: item
    end type

    type list_list_logical
        type(list_logical), dimension(:), allocatable :: item
    end type

    type list_integer_sys
        integer, dimension(:), allocatable :: item
    end type

    type list_list_integer_sys
        type(list_integer_sys), dimension(:), allocatable :: item
    end type

    type list_integer
        integer(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_integer
        type(list_integer), dimension(:), allocatable :: item
    end type

    type list_real
        real(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_real
        type(list_real), dimension(:), allocatable :: item
    end type

    type list_list_list_real
        type(list_list_real), dimension(:), allocatable :: item
    end type

    type list_4_real
        type(list_list_list_real), dimension(:), allocatable :: item
    end type

    type list_cmplx
        complex(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_cmplx
        type(list_cmplx), dimension(:), allocatable :: item
    end type

    type list_list_list_cmplx
        type(list_list_cmplx), dimension(:), allocatable :: item
    end type

    type list_4_cmplx
        type(list_list_list_cmplx), dimension(:), allocatable :: item
    end type

    ! probability distribution
    ! Deirmendjian modified gamma distribution
    !
    ! Reference: Comparison of laser beam propagation at 785 nm and 1550 nm in fog and haze for optical wireless communications,
    !            Kim, Isaac I. and {McArthur}, Bruce and Korevaar, Eric J.,
    !            eq. (8)
    type lib_math_gamma_distribution_deirmedjian_type
        real(kind=lib_math_type_kind) :: a
        real(kind=lib_math_type_kind) :: alpha
        real(kind=lib_math_type_kind) :: gamma
        real(kind=lib_math_type_kind) :: b
    end type

    type triangle_type
        type(cartesian_coordinate_real_type), dimension(3) :: vertex
        type(cartesian_coordinate_real_type) :: normal
    end type

end module lib_math_type
