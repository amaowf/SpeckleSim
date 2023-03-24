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

module lib_math_type_operator
    use lib_math_type
    use lib_math_constants
    implicit none

    private

    ! --- public ---
    public :: operator (+)
    public :: lib_math_real_array_add
    public :: lib_math_cmplx_array_add

    public :: operator (-)
    public :: lib_math_real_array_sub
    public :: lib_math_cmplx_array_sub

    public :: operator (*)
    public :: operator (/)
    public :: assignment (=)
    public :: sum
    public :: abs
    public :: spherical_abs
    public :: cartesian_abs
    public :: unit_vector

    public :: init_list
    public :: deallocate_list
    public :: make_list
    public :: make_array
    public :: get_structure
    public :: get_range
    public :: remove_zeros

    public :: make_cartesian
    public :: make_spherical

    public :: cross_product
    public :: dot_product

    public :: lib_math_type_operator_test_functions

    public :: lib_math_get_matrix_rot_x
    public :: lib_math_get_matrix_rot_y
    public :: lib_math_get_matrix_rot_z
    public :: lib_math_get_matrix_rot

    public :: list_filter

    ! ---- operator ----
    interface operator (+)
        ! real coordinates
        module procedure lib_math_cartesian_operator_real_add
        module procedure lib_math_cartesian_operator_real_0d_add_array
        module procedure lib_math_cartesian_operator_real_array_add_0d

        ! cmplx coordinates
        module procedure lib_math_cartesian_operator_add
        module procedure lib_math_cartesian_operator_add_array
        module procedure lib_math_cartesian_operator_0d_add_array
        module procedure lib_math_cartesian_operator_array_add_0d
        module procedure lib_math_list_cartesian_operator_add_array_cmplx
        module procedure lib_math_list_cartesian_operator_0d_add_array_cmplx
        module procedure lib_math_list_cartesian_operator_array_add_0d_cmplx

        module procedure lib_math_spherical_operator_add
        module procedure lib_math_spherical_operator_add_array
        module procedure lib_math_spherical_operator_0d_add_array
        module procedure lib_math_spherical_operator_array_add_0d
        module procedure lib_math_list_spherical_operator_add_array_cmplx
        module procedure lib_math_list_spherical_operator_0d_add_array_cmplx
        module procedure lib_math_list_spherical_operator_array_add_0d_cmplx

        ! list
        module procedure lib_math_list_real_add
        module procedure lib_math_list_cmplx_add

        ! list list
        module procedure lib_math_list_list_real_add
        module procedure lib_math_list_list_cmplx_add
    end interface

    interface operator (-)
        ! real coordinates
        module procedure lib_math_cartesian_operator_real_sub
        module procedure lib_math_cartesian_operator_real_array_sub_0d

        ! cmplx coordinates
        module procedure lib_math_cartesian_operator_sub
        module procedure lib_math_cartesian_operator_sub_array
        module procedure lib_math_cartesian_operator_array_sub_0d
        module procedure lib_math_list_cartesian_operator_sub_array_cmplx
        module procedure lib_math_list_cartesian_operator_array_sub_0d_cmplx

        module procedure lib_math_spherical_operator_sub
        module procedure lib_math_spherical_operator_sub_array
        module procedure lib_math_spherical_operator_array_sub_0d
        module procedure lib_math_list_spherical_operator_sub_array_cmplx
        module procedure lib_math_list_spherical_operator_array_sub_0d_cmplx

        ! list
        module procedure lib_math_list_real_sub
        module procedure lib_math_list_cmplx_sub

        ! list list
        module procedure lib_math_list_list_real_sub
        module procedure lib_math_list_list_cmplx_sub
    end interface

    interface operator (*)
        ! real coordinates
        module procedure lib_math_cartesian_operator_real_scalar_real_mul
        module procedure lib_math_cartesian_operator_real_scalar_cmplx_mul
        module procedure lib_math_cartesian_operator_real_mul_scalar_real
        module procedure lib_math_cartesian_operator_real_mul_scalar_cmplx

        module procedure lib_math_cartesian_operator_real_array_real_mul_array
        module procedure lib_math_cartesian_operator_real_array_cmplx_mul_array
        module procedure lib_math_cartesian_operator_real_0d_real_mul_array
        module procedure lib_math_cartesian_operator_real_0d_cmplx_mul_array

        module procedure lib_math_cartesian_operator_matrix_mul_real

        ! cmplx coordinates
        module procedure lib_math_cartesian_operator_scalar_real_mul
        module procedure lib_math_cartesian_operator_scalar_cmplx_mul
        module procedure lib_math_cartesian_operator_mul_scalar_real
        module procedure lib_math_cartesian_operator_mul_scalar_cmplx
        module procedure lib_math_cartesian_operator_array_real_mul_array
        module procedure lib_math_cartesian_operator_array_cmplx_mul_array
        module procedure lib_math_cartesian_operator_0d_real_mul_array
        module procedure lib_math_cartesian_operator_0d_cmplx_mul_array
        module procedure lib_math_list_cartesian_operator_array_real_mul_array_c
        module procedure lib_math_list_cartesian_operator_array_c_mul_array_cmplx
        module procedure lib_math_list_cartesian_operator_real_mul_array_cmplx
        module procedure lib_math_list_cartesian_operator_cmplx_mul_array_cmplx

        module procedure lib_math_spherical_operator_scalar_real_mul
        module procedure lib_math_spherical_operator_scalar_cmplx_mul
        module procedure lib_math_spherical_operator_mul_scalar_real
        module procedure lib_math_spherical_operator_mul_scalar_cmplx
        module procedure lib_math_spherical_operator_array_real_mul_array
        module procedure lib_math_spherical_operator_array_cmplx_mul_array
        module procedure lib_math_spherical_operator_0d_real_mul_array
        module procedure lib_math_spherical_operator_0d_cmplx_mul_array
        module procedure lib_math_list_spherical_operator_array_real_mul_array_c
        module procedure lib_math_list_spherical_operator_array_c_mul_array_cmplx
        module procedure lib_math_list_spherical_operator_real_mul_array_cmplx
        module procedure lib_math_list_spherical_operator_cmplx_mul_array_cmplx

        module procedure lib_math_cartesian_operator_matrix_mul_cmplx

        module procedure lib_math_cartesian_operator_matrix_mul_matrix

        ! list
        module procedure lib_math_list_real_mul_real
        module procedure lib_math_list_real_mul_cmplx

        module procedure lib_math_list_real_mul_list_real
        module procedure lib_math_list_cmplx_mul_list_cmplx

        ! list list
        module procedure lib_math_list_list_real_mul_real
        module procedure lib_math_real_mul_list_list_cmplx
        module procedure lib_math_list_list_cmplx_mul_real

        module procedure lib_math_list_list_cmplx_mul_cmplx
        module procedure lib_math_cmplx_mul_list_list_cmplx

        module procedure lib_math_list_cmplx_mul_list_list_cmplx
        module procedure lib_math_list_list_cmplx_mul_list_list_cmplx
    end interface

    interface operator (/)
        ! real coordinates
        module procedure lib_math_cartesian_operator_real_divide_by_scalar_real
        module procedure lib_math_cartesian_operator_real_divide_by_scalar_cmplx

        ! cmplx coordinates
        module procedure lib_math_cartesian_operator_divide_by_scalar_real
        module procedure lib_math_cartesian_operator_divide_by_scalar_cmplx
        module procedure lib_math_cartesian_operator_array_divide_by_scalar_array_real
        module procedure lib_math_cartesian_operator_array_divide_by_scalar_array_cmplx
        module procedure lib_math_cartesian_operator_array_divide_by_0d_real
        module procedure lib_math_cartesian_operator_array_divide_by_0d_cmplx
        module procedure lib_math_list_cartesian_operator_array_divide_by_real_array
        module procedure lib_math_list_cartesian_operator_c_array_divide_by_c_array
        module procedure lib_math_list_cartesian_operator_array_divide_by_real
        module procedure lib_math_list_cartesian_operator_array_divide_by_cmplx

        module procedure lib_math_spherical_operator_divide_by_scalar_real
        module procedure lib_math_spherical_operator_divide_by_scalar_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_real
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_0d_real
        module procedure lib_math_spherical_operator_array_divide_by_0d_cmplx
        module procedure lib_math_list_spherical_operator_array_divide_by_real_array
        module procedure lib_math_list_spherical_operator_c_array_divide_by_c_array
        module procedure lib_math_list_spherical_operator_array_divide_by_real
        module procedure lib_math_list_spherical_operator_array_divide_by_cmplx
    end interface

    interface assignment (=)
        module procedure lib_math_list_real_assignment
        module procedure lib_math_list_cmplx_assignment

        module procedure lib_math_list_list_real_assignment
        module procedure lib_math_list_list_cmplx_assignment

        module procedure lib_math_spherical_point_to_cartesian_point
        module procedure lib_math_cartesian_point_to_spherical_point

        module procedure lib_math_array_to_cartesian_point
        module procedure lib_math_cartesian_point_to_array
    end interface

    interface spherical_abs
        module procedure lib_math_spherical_operator_abs_cmplx
    end interface

    interface cartesian_abs
        module procedure lib_math_cartesian_operator_abs_cmplx
        module procedure lib_math_cartesian_operator_abs_real
    end interface

    interface abs
        module procedure lib_math_cartesian_operator_abs_cmplx
        module procedure lib_math_cartesian_operator_abs_real

        module procedure lib_math_spherical_operator_abs_cmplx
    end interface

    interface unit_vector
        module procedure lib_math_cartesian_operator_unit_vector_real
    end interface

    interface sum
        module procedure lib_math_list_list_cmplx_sum
        module procedure lib_math_list_list_real_sum
    end interface

    interface init_list
        module procedure lib_math_list_list_logical_init
        module procedure lib_math_list_list_integer_sys_init
        module procedure lib_math_list_list_integer_init
        module procedure lib_math_list_list_real_init
        module procedure lib_math_list_4_real_init
        module procedure lib_math_list_list_cmplx_init
        module procedure lib_math_list_list_cmplx_init_with_shape
        module procedure lib_math_list_4_cmplx_init
        module procedure lib_math_list_real_init
        module procedure lib_math_list_cmplx_init
        module procedure lib_math_list_spherical_coordinate_cmplx_type_init
    end interface

    interface deallocate_list
        module procedure lib_math_list_real_deallocate
        module procedure lib_math_list_cmplx_deallocate

        module procedure lib_math_list_list_real_deallocate
        module procedure lib_math_list_list_cmplx_deallocate
    end interface

    interface make_list
        module procedure lib_math_array_make_list_list_cmplx
        module procedure lib_math_array_make_list_of_list_list_cmplx
    end interface

    interface make_array
        module procedure lib_math_list_list_make_array_cmplx
        module procedure lib_math_list_of_list_list_make_array_cmplx
    end interface

    interface get_structure
        module procedure lib_math_get_list_list_structure_cmplx
        module procedure lib_math_get_list_of_list_list_structure_cmplx
    end interface

    interface get_range
        module procedure lib_math_get_list_list_cmplx_with_range
    end interface

    interface remove_zeros
        module procedure lib_math_list_cmplx_remove_zeros
        module procedure lib_math_list_list_cmplx_remove_zeros
    end interface

    interface make_cartesian
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_a
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_c
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_s

        module procedure lib_math_make_cartesian_coordiante_real
        module procedure lib_math_make_cartesian_coordiante_cmplx
    end interface

    interface make_spherical
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_a
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_c
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_s
        module procedure lib_math_make_spherical_coordiante
    end interface

    interface cross_product
        module procedure lib_math_cartesian_cross_product_real
        module procedure lib_math_cartesian_cross_product_real_list
        module procedure lib_math_cartesian_cross_product_cmplx
        module procedure lib_math_cartesian_cross_product_cmplx_list
    end interface

    interface dot_product
        module procedure lib_math_cartesian_dot_product_real
        module procedure lib_math_cartesian_dot_product_real_list
        module procedure lib_math_cartesian_dot_product_cmplx
        module procedure lib_math_cartesian_dot_product_cmplx_list
    end interface

    interface list_filter
        module procedure list_filter_cartesian_coordiante_real
        module procedure list_filter_cartesian_coordiante_cmplx
        module procedure list_filter_spherical_coordiante_real
        module procedure list_filter_spherical_coordiante_cmplx
    end interface

    contains

! ---- single spherical coordinate ----
        function lib_math_spherical_operator_add(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho + rhs%rho
            rv%phi = lhs%phi + rhs%phi
            rv%theta = lhs%theta + rhs%theta

        end function

        function lib_math_spherical_operator_sub(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho - rhs%rho
            rv%phi = lhs%phi - rhs%phi
            rv%theta = lhs%theta - rhs%theta

        end function

        function lib_math_spherical_operator_scalar_real_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs * rhs%rho
            rv%phi = lhs * rhs%phi
            rv%theta = lhs * rhs%theta

        end function

        function lib_math_spherical_operator_scalar_cmplx_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=8), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs * rhs%rho
            rv%phi = lhs * rhs%phi
            rv%theta = lhs * rhs%theta

        end function

        function lib_math_spherical_operator_mul_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho * rhs
            rv%phi = lhs%phi * rhs
            rv%theta = lhs%theta * rhs

        end function

        function lib_math_spherical_operator_mul_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho * rhs
            rv%phi = lhs%phi * rhs
            rv%theta = lhs%theta * rhs

        end function

        function lib_math_spherical_operator_divide_by_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho / rhs
            rv%phi = lhs%phi / rhs
            rv%theta = lhs%theta / rhs

        end function

        function lib_math_spherical_operator_divide_by_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho / rhs
            rv%phi = lhs%phi / rhs
            rv%theta = lhs%theta / rhs

        end function

        function lib_math_spherical_operator_abs_cmplx(rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=8) :: rv

            rv = sqrt(rhs%rho*rhs%rho + rhs%phi*rhs%phi + rhs%theta*rhs%theta)

        end function

! ---- array spherical coordinate ----
        function lib_math_spherical_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_0d_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- list_spherical_coordinate ----
        function lib_math_list_spherical_operator_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) + rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if

        end function

        function lib_math_list_spherical_operator_sub_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) - rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if
        end function

        function lib_math_list_spherical_operator_0d_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs + rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_add_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_sub_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_real_mul_array_c(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(lhs,1):ubound(lhs,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_spherical_operator_array_real_mul_array_c: ERROR"
            end if

        end function

        function lib_math_list_spherical_operator_array_c_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_spherical_operator_real_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_cmplx_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_divide_by_real_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_spherical_operator_c_array_divide_by_c_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_spherical_operator_array_divide_by_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_divide_by_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- single cartesian coordinate ----
        function lib_math_cartesian_operator_add(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x + rhs%x
            rv%y = lhs%y + rhs%y
            rv%z = lhs%z + rhs%z

        end function

        function lib_math_cartesian_operator_sub(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x - rhs%x
            rv%y = lhs%y - rhs%y
            rv%z = lhs%z - rhs%z

        end function

        function lib_math_cartesian_operator_scalar_real_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_scalar_cmplx_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_mul_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_mul_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_divide_by_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

        function lib_math_cartesian_operator_divide_by_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

                function lib_math_cartesian_operator_real_add(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_real_type) :: rv

            rv%x = lhs%x + rhs%x
            rv%y = lhs%y + rhs%y
            rv%z = lhs%z + rhs%z

        end function

                function lib_math_cartesian_operator_real_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            type (cartesian_coordinate_real_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_real_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_sub(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_real_type) :: rv

            rv%x = lhs%x - rhs%x
            rv%y = lhs%y - rhs%y
            rv%z = lhs%z - rhs%z

        end function

        function lib_math_cartesian_operator_real_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_scalar_real_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_real_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_real_scalar_cmplx_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_real_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_real_mul_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_real_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_real_mul_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_real_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_real_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_real_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_real_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_real_divide_by_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_real_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

        function lib_math_cartesian_operator_real_divide_by_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

        function lib_math_cartesian_operator_abs_real(rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: rhs

            real(kind=lib_math_type_kind) :: rv

            rv = sqrt(rhs%x*rhs%x + rhs%y*rhs%y + rhs%z*rhs%z)

        end function

        function lib_math_cartesian_operator_unit_vector_real(rhs) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type) :: rv

            ! auxiliary
            real(kind=lib_math_type_kind) :: buffer

            buffer = abs(rhs)

            if (buffer .gt. 0) then
                rv = rhs / buffer
            else
                rv = make_cartesian(0d0, 0d0, 0d0)
                print *, "lib_math_cartesian_operator_unit_vector_real: ERROR"
                print *, "  division by zero"
            end if

        end function

        function lib_math_cartesian_operator_abs_cmplx(rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=8) :: rv

            rv = sqrt(rhs%x*rhs%x + rhs%y*rhs%y + rhs%z*rhs%z)

        end function

! ---- array cartesian coordinate ----
        function lib_math_cartesian_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_scalar_array_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_scalar_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_0d_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_matrix_mul_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(cartresian_coordinate_rot_matrix_type), intent(in) :: lhs
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type) :: rv

            rv%x = lhs%r_11 * rhs%x + lhs%r_12 * rhs%y + lhs%r_13 * rhs%z
            rv%y = lhs%r_21 * rhs%x + lhs%r_22 * rhs%y + lhs%r_23 * rhs%z
            rv%z = lhs%r_31 * rhs%x + lhs%r_32 * rhs%y + lhs%r_33 * rhs%z

        end function

        function lib_math_cartesian_operator_matrix_mul_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(cartresian_coordinate_rot_matrix_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%r_11 * rhs%x + lhs%r_12 * rhs%y + lhs%r_13 * rhs%z
            rv%y = lhs%r_21 * rhs%x + lhs%r_22 * rhs%y + lhs%r_23 * rhs%z
            rv%z = lhs%r_31 * rhs%x + lhs%r_32 * rhs%y + lhs%r_33 * rhs%z

        end function

        ! Formula: c = a * b
        !
        !
        !                 b_11 b_12 b_13
        !                 b_21 b_22 b_23
        !                 b_31 b_32 b_33
        !
        ! a_11 a_12 a_13  c_11 c_12 c_13
        ! a_21 a_22 a_23  c_21 c_22 c_23
        ! a_31 a_32 a_33  c_31 c_32 c_33
        !
        !
        ! c_11 = a_11 * b_11 + a_12 * b_21 + a_13 * b_31
        ! c_12 = a_11 * b_12 + a_12 * b_22 + a_13 * b_32
        ! c_13 = a_11 * b_13 + a_12 * b_23 + a_13 * b_33
        !
        ! c_21 = a_21 * b_11 + a_22 * b_21 + a_23 * b_31
        ! c_22 = a_21 * b_12 + a_22 * b_22 + a_23 * b_32
        ! c_23 = a_21 * b_13 + a_22 * b_23 + a_23 * b_33
        !
        ! c_31 = a_31 * b_11 + a_32 * b_21 + a_33 * b_31
        ! c_32 = a_31 * b_12 + a_32 * b_22 + a_33 * b_32
        ! c_33 = a_31 * b_13 + a_32 * b_23 + a_33 * b_33
        !
        function lib_math_cartesian_operator_matrix_mul_matrix(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(cartresian_coordinate_rot_matrix_type), intent(in) :: lhs
            type(cartresian_coordinate_rot_matrix_type), intent(in) :: rhs

            type(cartresian_coordinate_rot_matrix_type) :: rv

            rv%r_11 = lhs%r_11 * rhs%r_11 + lhs%r_12 * rhs%r_21 + lhs%r_13 * rhs%r_31
            rv%r_12 = lhs%r_11 * rhs%r_12 + lhs%r_12 * rhs%r_22 + lhs%r_13 * rhs%r_32
            rv%r_13 = lhs%r_11 * rhs%r_13 + lhs%r_12 * rhs%r_23 + lhs%r_13 * rhs%r_33

            rv%r_21 = lhs%r_21 * rhs%r_11 + lhs%r_22 * rhs%r_21 + lhs%r_23 * rhs%r_31
            rv%r_22 = lhs%r_21 * rhs%r_12 + lhs%r_22 * rhs%r_22 + lhs%r_23 * rhs%r_32
            rv%r_23 = lhs%r_21 * rhs%r_13 + lhs%r_22 * rhs%r_23 + lhs%r_23 * rhs%r_33

            rv%r_31 = lhs%r_31 * rhs%r_11 + lhs%r_32 * rhs%r_21 + lhs%r_33 * rhs%r_31
            rv%r_32 = lhs%r_31 * rhs%r_12 + lhs%r_32 * rhs%r_22 + lhs%r_33 * rhs%r_32
            rv%r_33 = lhs%r_31 * rhs%r_13 + lhs%r_32 * rhs%r_23 + lhs%r_33 * rhs%r_33

        end function

        ! "coordinate system rotations of the x-, y-, and z-axes in a counterclockwise direction 
        !  when looking towards the origin give the matrices"
        !
        ! Argument
        ! ----
        !   alpha: double precision
        !       angle [rad]
        !
        ! Returns
        ! ----
        !   rv: type(cartresian_coordinate_rot_matrix_type)
        !       rotation matrix
        !
        ! Reference: http://mathworld.wolfram.com/RotationMatrix.html
        function lib_math_get_matrix_rot_x(alpha) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: alpha

            type(cartresian_coordinate_rot_matrix_type) :: rv

            ! auxiliary
            double precision :: cos_alpha
            double precision :: sin_alpha

            cos_alpha = cos(alpha)
            sin_alpha = sin(alpha)

            rv%r_11 = 1
            rv%r_12 = 0
            rv%r_13 = 0

            rv%r_21 = 0
            rv%r_22 = cos_alpha
            rv%r_23 = sin_alpha

            rv%r_31 = 0
            rv%r_32 = -sin_alpha
            rv%r_33 = cos_alpha

        end function lib_math_get_matrix_rot_x

        ! "coordinate system rotations of the x-, y-, and z-axes in a counterclockwise direction 
        !  when looking towards the origin give the matrices"
        !
        ! Argument
        ! ----
        !   alpha: double precision
        !       angle [rad]
        !
        ! Returns
        ! ----
        !   rv: type(cartresian_coordinate_rot_matrix_type)
        !       rotation matrix
        !
        ! Reference: http://mathworld.wolfram.com/RotationMatrix.html
        function lib_math_get_matrix_rot_y(beta) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: beta

            type(cartresian_coordinate_rot_matrix_type) :: rv

            ! auxiliary
            double precision :: cos_beta
            double precision :: sin_beta

            cos_beta = cos(beta)
            sin_beta = sin(beta)

            rv%r_11 = cos_beta
            rv%r_12 = 0
            rv%r_13 = -sin_beta

            rv%r_21 = 0
            rv%r_22 = 1
            rv%r_23 = 0

            rv%r_31 = sin_beta
            rv%r_32 = 0
            rv%r_33 = cos_beta

        end function lib_math_get_matrix_rot_y

        ! "coordinate system rotations of the x-, y-, and z-axes in a counterclockwise direction 
        !  when looking towards the origin give the matrices"
        !
        ! Argument
        ! ----
        !   alpha: double precision
        !       angle [rad]
        !
        ! Returns
        ! ----
        !   rv: type(cartresian_coordinate_rot_matrix_type)
        !       rotation matrix
        !
        ! Reference: http://mathworld.wolfram.com/RotationMatrix.html
        function lib_math_get_matrix_rot_z(gamma) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: gamma

            type(cartresian_coordinate_rot_matrix_type) :: rv

            ! auxiliary
            double precision :: cos_gamma
            double precision :: sin_gamma

            cos_gamma = cos(gamma)
            sin_gamma = sin(gamma)

            rv%r_11 = cos_gamma
            rv%r_12 = sin_gamma
            rv%r_13 = 0

            rv%r_21 = -sin_gamma
            rv%r_22 = cos_gamma
            rv%r_23 = 0

            rv%r_31 = 0
            rv%r_32 = 0
            rv%r_33 = 1

        end function lib_math_get_matrix_rot_z

        ! Argument
        ! ----
        !   phi: double precsision
        !       rotation about the z-axis [rad]
        !   theta: double precsision
        !       rotation about the x'-axis [0, Pi] [rad]
        !   psi: double precision
        !       rotation about the z'-axis [rad]
        !
        ! Returns
        ! ----
        !   rv: type(cartresian_coordinate_rot_matrix_type)
        !       rotation matrix
        !
        ! Convention: x-convention
        !
        ! Refrence: http://mathworld.wolfram.com/EulerAngles.html
        !           eq. 6..14
        !
        function lib_math_get_matrix_rot(phi, theta, psi) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: psi
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi

            type(cartresian_coordinate_rot_matrix_type) :: rv

            ! auxiliary
            double precision :: cos_theta
            double precision :: sin_theta

            double precision :: cos_phi
            double precision :: sin_phi

            double precision :: cos_psi
            double precision :: sin_psi

            cos_phi = cos(phi)
            sin_phi = sin(phi)

            cos_theta = cos(theta)
            sin_theta = sin(theta)

            cos_psi = cos(psi)
            sin_psi = sin(psi)

            ! Refrence: http://mathworld.wolfram.com/EulerAngles.html
            !           eq. 6..14
            rv%r_11 = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi
            rv%r_12 = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi
            rv%r_13 = sin_psi * sin_theta

            rv%r_21 = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi
            rv%r_22 = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi
            rv%r_23 = cos_psi * sin_theta

            rv%r_31 = sin_theta * sin_phi
            rv%r_32 = -sin_theta * cos_phi
            rv%r_33 = cos_theta

        end function lib_math_get_matrix_rot

! ---- list_cartesian_coordinate ----
        function lib_math_list_cartesian_operator_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) + rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if

        end function

        function lib_math_list_cartesian_operator_sub_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) - rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if
        end function

        function lib_math_list_cartesian_operator_0d_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs + rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_add_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_sub_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_real_mul_array_c(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(lhs,1):ubound(lhs,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_cartesian_operator_array_real_mul_array_c: ERROR"
            end if

        end function

        function lib_math_list_cartesian_operator_array_c_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_real_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_cmplx_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_divide_by_real_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_c_array_divide_by_c_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_array_divide_by_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_divide_by_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- list_list ----

        subroutine lib_math_list_real_deallocate(list)
            ! dummy
            type(list_real), intent(inout) :: list

            if (allocated(list%item)) then
                deallocate(list%item)
            end if
        end subroutine

        subroutine lib_math_list_cmplx_deallocate(list)
            ! dummy
            type(list_cmplx), intent(inout) :: list

            if (allocated(list%item)) then
                deallocate(list%item)
            end if
        end subroutine

        subroutine lib_math_list_real_assignment(lhs, rhs)
            implicit none
            ! dummy
            type(list_real), intent(in) :: rhs

            type(list_real), intent(inout) :: lhs

            call lib_math_list_real_deallocate(lhs)

            if (allocated(rhs%item)) then
                lhs%item = rhs%item
            end if
        end subroutine

        subroutine lib_math_list_cmplx_assignment(lhs, rhs)
            implicit none
            ! dummy
            type(list_cmplx), intent(in) :: rhs

            type(list_cmplx), intent(inout) :: lhs

            call lib_math_list_cmplx_deallocate(lhs)

            if (allocated(rhs%item)) then
                lhs%item = rhs%item
            end if
        end subroutine

        subroutine lib_math_list_list_real_deallocate(list)
            ! dummy
            type(list_list_real), intent(inout) :: list

            ! auxiliary
            integer :: i

            if (allocated(list%item)) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = lbound(list%item, 1), ubound(list%item, 1)
                    if (allocated(list%item(i)%item)) then
                        deallocate(list%item(i)%item)
                    end if
                end do
                !$OMP END PARALLEL DO
                deallocate(list%item)
            end if
        end subroutine lib_math_list_list_real_deallocate

        subroutine lib_math_list_list_cmplx_deallocate(list)
            ! dummy
            type(list_list_cmplx), intent(inout) :: list

            ! auxiliary
            integer :: i

            if (allocated(list%item)) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = lbound(list%item, 1), ubound(list%item, 1)
                    if (allocated(list%item(i)%item)) then
                        deallocate(list%item(i)%item)
                    end if
                end do
                !$OMP END PARALLEL DO
                deallocate(list%item)
            end if
        end subroutine lib_math_list_list_cmplx_deallocate

        subroutine lib_math_list_list_real_assignment(lhs, rhs)
            implicit none
            ! dummy
            type(list_list_real), intent(in) :: rhs

            type(list_list_real), intent(inout) :: lhs

            call lib_math_list_list_real_deallocate(lhs)

            if (allocated(rhs%item)) then
                lhs%item = rhs%item
            end if
        end subroutine

        subroutine lib_math_list_list_cmplx_assignment(lhs, rhs)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx), intent(inout) :: lhs

            call lib_math_list_list_cmplx_deallocate(lhs)

            if (allocated(rhs%item)) then
                lhs%item = rhs%item
            end if
        end subroutine

        ! Arguments
        ! ----
        !   list: type (list_list_logical)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_logical_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_list_logical), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            logical, intent(in), optional :: init_value

            ! auxiliary
            integer :: i
            integer :: ii

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        list%item(i)%item(ii) = init_value
                    end do
                end do
            end if

        end subroutine lib_math_list_list_logical_init

        ! Arguments
        ! ----
        !   list: type (list_list_integer_sys)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_integer_sys_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_list_integer_sys), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            integer, intent(in), optional :: init_value

            ! auxiliary
            integer :: i
            integer :: ii

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        list%item(i)%item(ii) = init_value
                    end do
                end do
            end if

        end subroutine lib_math_list_list_integer_sys_init

        ! Arguments
        ! ----
        !   list: type (list_list_integer)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_integer_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_list_integer), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            integer(kind=lib_math_type_kind), intent(in), optional :: init_value

            ! auxiliary
            integer :: i
            integer :: ii

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        list%item(i)%item(ii) = init_value
                    end do
                end do
            end if

        end subroutine lib_math_list_list_integer_init

        ! Arguments
        ! ----
        !   list: type (list_list_real)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_real_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_list_real), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            real(kind=lib_math_type_kind), intent(in), optional :: init_value

            ! auxiliary
            integer :: i
            integer :: ii

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    do ii=-i, i
                        list%item(i)%item(ii) = init_value
                    end do
                end do
            end if

        end subroutine lib_math_list_list_real_init

        ! Arguments
        ! ----
        !   list: type (list_list_cmplx)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_cmplx_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_list_cmplx), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            complex(kind=lib_math_type_kind), intent(in), optional :: init_value

            ! auxiliary
            integer :: i

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    list%item(i)%item = init_value
                end do
            end if

        end subroutine lib_math_list_list_cmplx_init

        ! Arguments
        ! ----
        !   list: type (list_list_cmplx)
        !       derived data type to initialize
        !   shape: type (list_list_cmplx)
        !       shape with which the list is to be initialized
        subroutine lib_math_list_list_cmplx_init_with_shape(list, shape)
            implicit none
            ! dummy
            type (list_list_cmplx), intent(inout) :: list
            type (list_list_cmplx), intent(in) :: shape

            ! auxiliary
            integer :: i

            call deallocate_list(list)

            allocate(list%item(lbound(shape%item, 1):ubound(shape%item, 1)))

            do i = lbound(list%item, 1), ubound(list%item, 1)
                allocate(list%item(i)%item(lbound(shape%item(i)%item, 1):ubound(shape%item(i)%item, 1)))
            end do

        end subroutine lib_math_list_list_cmplx_init_with_shape

        ! Arguments
        ! ----
        !   list: type (list_real)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_real_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_real), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            real(kind=lib_math_type_kind), intent(in), optional :: init_value

            ! auxiliary
            integer :: i

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    list%item(i) = init_value
                end do
            end if

        end subroutine lib_math_list_real_init

        ! Arguments
        ! ----
        !   list: type (list_cmplx)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_cmplx_init(list, fnu, n, init_value)
            implicit none
            ! dummy
            type (list_cmplx), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n
            complex(kind=lib_math_type_kind), intent(in), optional :: init_value

            ! auxiliary
            integer :: i

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )

            if (present(init_value)) then
                do i=fnu, fnu+n-1
                    list%item(i) = init_value
                end do
            end if

        end subroutine lib_math_list_cmplx_init

        ! Arguments
        ! ----
        !   list: type (list_4_real)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_4_real_init(list, fnu_1, n_1, fnu_2, n_2)
            implicit none
            ! dummy
            type (list_4_real), intent(inout) :: list
            integer, intent(in) :: fnu_1
            integer, intent(in) :: n_1
            integer, intent(in) :: fnu_2
            integer, intent(in) :: n_2

            ! auxiliary
            integer :: n
            integer :: m
            integer :: nu

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu_1:fnu_1+n_1-1) )
            do n=fnu_1, fnu_1+n_1-1
                allocate (list%item(n)%item(-n:n) )
                do m=-n, n
                    allocate (list%item(n)%item(m)%item(fnu_2:fnu_2+n_2-1))
                    do nu=fnu_2, fnu_2+n_2-1
                        allocate (list%item(n)%item(m)%item(nu)%item(-nu:nu))
                    end do
                end do
            end do

        end subroutine lib_math_list_4_real_init

        ! Arguments
        ! ----
        !   list: type (list_4_cmplx)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_4_cmplx_init(list, fnu_1, n_1, fnu_2, n_2)
            implicit none
            ! dummy
            type (list_4_cmplx), intent(inout) :: list
            integer, intent(in) :: fnu_1
            integer, intent(in) :: n_1
            integer, intent(in) :: fnu_2
            integer, intent(in) :: n_2

            ! auxiliary
            integer :: n
            integer :: m
            integer :: nu

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu_1:fnu_1+n_1-1) )
            do n=fnu_1, fnu_1+n_1-1
                allocate (list%item(n)%item(-n:n) )
                do m=-n, n
                    allocate (list%item(n)%item(m)%item(fnu_2:fnu_2+n_2-1))
                    do nu=fnu_2, fnu_2+n_2-1
                        allocate (list%item(n)%item(m)%item(nu)%item(-nu:nu))
                    end do
                end do
            end do

        end subroutine lib_math_list_4_cmplx_init

        ! Argument
        ! ----
        !   rhs: type(list_list_cmplx)
        !       list of list with complex numbers
        !
        ! Returns
        ! ----
        !   rv: complex(kind=lib_math_type_kind)
        !       the sum of all elements of *rhs*
        function lib_math_list_list_real_sum(rhs) result(rv)
            implicit none
            ! dummy
            type(list_list_real), intent(in) :: rhs

            real(kind=lib_math_type_kind) :: rv

            ! auxiliry
            integer :: n
            real(kind=lib_math_type_kind), dimension(lbound(rhs%item, 1):ubound(rhs%item, 1)) :: buffer

            !$OMP PARALLEL DO PRIVATE(n)
            do n = lbound(rhs%item, 1), ubound(rhs%item, 1)
                buffer(n) = sum(rhs%item(n)%item)
            end do
            !$OMP END PARALLEL DO

            rv = sum(buffer)

        end function lib_math_list_list_real_sum

        ! Argument
        ! ----
        !   rhs: type(list_list_cmplx)
        !       list of list with complex numbers
        !
        ! Returns
        ! ----
        !   rv: complex(kind=lib_math_type_kind)
        !       the sum of all elements of *rhs*
        function lib_math_list_list_cmplx_sum(rhs) result(rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: rhs

            complex(kind=lib_math_type_kind) :: rv

            ! auxiliry
            integer :: n
            complex(kind=lib_math_type_kind), dimension(lbound(rhs%item, 1):ubound(rhs%item, 1)) :: buffer

            if (allocated(rhs%item)) then
                !$OMP PARALLEL DO PRIVATE(n)
                do n = lbound(rhs%item, 1), ubound(rhs%item, 1)
                    if (allocated(rhs%item(n)%item)) then
                        buffer(n) = sum(rhs%item(n)%item)
                    else
                        buffer(n) = 0
                    end if
                end do
                !$OMP END PARALLEL DO

                rv = sum(buffer)
            else
                rv = 0
            end if

        end function lib_math_list_list_cmplx_sum

        ! Argument
        ! ----
        !   list: type(list_list_cmplx)
        !       list of lists like: list%item(2)%item(-2:2)
        !
        ! Returns
        ! ----
        !   array: complex(kind=lib_math_type_kind), dimension(:)
        !       all elements of list as an array
        !
        !   list%item(n)%item(m): element (n,m) -> array element
        !
        !   array | (1,-1) | (1,0) | (1,1) | (2,-2) | (2,-1) | (2,0) | ...
        !
        subroutine lib_math_list_list_make_array_cmplx(list, array, fnu, n)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: list
            integer, intent(in), optional :: fnu
            integer, intent(in), optional :: n

            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(inout) :: array

            ! auxiliary
            integer :: i
            integer :: n_min
            integer :: n_max
            integer :: first
            integer :: last
            integer, dimension(2) :: n_range

            if (present(fnu) .and. present(n)) then
                n_range(1) = fnu
                n_range(2) = fnu + n - 1
                n_min = max(n_range(1), lbound(list%item, 1))
                n_max = min(n_range(2), ubound(list%item, 1))
            else
                n_range(1) = lbound(list%item, 1)
                n_range(2) = ubound(list%item, 1)
                n_min = n_range(1)
                n_max = n_range(2)
            end if

            if (allocated(array)) then
                deallocate(array)
            end if
            i = (1 + n_range(2))**2 - n_range(1)**2
            allocate(array(i))
            array = dcmplx(0,0)

            if ( n_min .gt. n_range(1) ) then
                ! (1 + (n_min-1))**2 - n_range(1)**2
                first = n_min**2 - n_range(1)**2
            else
                first = 1
            end if
            do i = n_min, n_max
                last = first + 2 * i ! +1 -1
                array(first:last) = list%item(i)%item
                first = last + 1
            end do

        end subroutine lib_math_list_list_make_array_cmplx

        ! Argument
        ! ----
        !   array: complex(kind=lib_math_type_kind), dimension(:)
        !       all elements of list as an array
        !   fnu: integer
        !       start index of the list
        !   n: integer
        !       number of elements, n .GE. 1
        !
        !
        ! Returns
        ! ----
        !   list: type(list_list_cmplx)
        !       list of lists like: list%item(2)%item(-2:2)
        !
        !   list%item(n)%item(m): element (n,m) -> array element
        !
        !   array | (1,-1) | (1,0) | (1,1) | (2,-2) | (2,-1) | (2,0) | ...
        !
        subroutine lib_math_array_make_list_list_cmplx(array, fnu, n, list)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: array
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            type(list_list_cmplx), intent(inout) :: list

            ! auxiliary
            integer :: m_i
            integer :: m_no_of_elements
            integer :: m_n
            integer :: m_m
            integer :: m_first
            integer :: m_last

            m_no_of_elements = (1 + fnu+n-1)**2 - fnu**2
            if (size(array, 1) .ge. m_no_of_elements) then
                call init_list(list, fnu, n)

!                m_i = lbound(array, 1)
                m_first = lbound(array, 1)
                do m_n = fnu, fnu + n - 1
                    m_last = m_first + 2 * m_n ! + 1 - 1
                    list%item(m_n)%item(:) = array(m_first:m_last)
                    m_first = m_last + 1
!                    do m_m = -m_n, m_n
!                        list%item(m_n)%item(m_m) = array(m_i)
!                        m_i = m_i + 1
!                    end do
                end do
            else
                call init_list(list, fnu, n, cmplx(0, 0, kind=lib_math_type_kind))

                m_i = ubound(array, 1)
                do m_n = fnu, fnu + n - 1
                    do m_m = -m_n, m_n
                        list%item(m_n)%item(m_m) = array(m_i)
                        m_i = m_i + 1
                        if (m_i .ge. m_no_of_elements) then
                            return
                        end if
                    end do
                end do
            end if
        end subroutine lib_math_array_make_list_list_cmplx

        ! Argument
        ! ----
        !   list: type(list_list_cmplx), dimension(:)
        !       list of lists like: list%item(2)%item(-2:2)
        !
        ! Returns
        ! ----
        !   array: complex(kind=lib_math_type_kind), dimension(:)
        !       all elements of list as an array
        !
        !   list%item(n)%item(m): element (n,m) -> array element
        !
        !   array | (1,-1) | (1,0) | (1,1) | (2,-2) | (2,-1) | (2,0) | ...
        !
        subroutine lib_math_list_of_list_list_make_array_cmplx(list, array)
            implicit none
            type(list_list_cmplx), dimension(:), intent(in) :: list

            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(inout) :: array

            ! auxiliary
            integer :: i
            integer, dimension(lbound(list, 1):ubound(list, 1)) :: no_of_elements
            integer :: start
            integer :: last
            integer, dimension(lbound(list, 1):ubound(list, 1), 2) :: n_range
            complex(kind=lib_math_type_kind), dimension(:), allocatable :: buffer_array

            do i = lbound(list, 1), ubound(list, 1)
                n_range(i, 1) = lbound(list(i)%item, 1)
                n_range(i, 2) = ubound(list(i)%item, 1)

                no_of_elements(i) = (1 + n_range(i, 2))**2 - (n_range(i, 1))**2
            end do

            if (allocated(array)) then
                deallocate(array)
            end if
            allocate(array(sum(no_of_elements)))

            !$OMP PARALLEL DO PRIVATE(i, buffer_array, start, last)
            do i = lbound(list, 1), ubound(list, 1)
                call lib_math_list_list_make_array_cmplx(list(i), buffer_array)

                start = sum(no_of_elements(lbound(list, 1):i)) - no_of_elements(i) + 1
                last = start + no_of_elements(i) - 1

                array(start:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_math_list_of_list_list_make_array_cmplx

        subroutine lib_math_array_make_list_of_list_list_cmplx(array, fnu, n, list)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: array

            integer, dimension(:) :: fnu
            integer, dimension(lbound(fnu, 1):ubound(fnu, 1)) :: n
            type(list_list_cmplx), dimension(:), allocatable, intent(inout) :: list

            ! auxiliary
            integer :: i
            integer :: start
            integer :: last
            integer, dimension(lbound(fnu, 1):ubound(fnu, 1)) :: no_of_elements
            integer, dimension(lbound(fnu, 1):ubound(fnu, 1), 2) :: n_range
            type(list_list_cmplx) :: buffer_list


            if (allocated(list)) then
                deallocate(list)
            end if
            allocate(list(lbound(fnu, 1):ubound(fnu, 1)))

            no_of_elements = 0
            do i = lbound(fnu, 1), ubound(fnu, 1)
                n_range(i, 1) = fnu(i)
                n_range(i, 2) = fnu(i) + n(i) - 1

                no_of_elements(i) = (1 + n_range(i, 2))**2 - n_range(i, 1)**2

                call init_list(list(i), n_range(i, 1), n_range(i, 2) - n_range(i, 1) + 1)
            end do

            if (sum(no_of_elements) .eq. size(array, 1)) then

                !$OMP PARALLEL DO PRIVATE(i, buffer_list, start, last)
                do i = lbound(fnu, 1), ubound(fnu, 1)

                    start = sum(no_of_elements(lbound(fnu, 1):i)) - no_of_elements(i) + 1
                    last = start + no_of_elements(i) - 1
                    call lib_math_array_make_list_list_cmplx(array(start:last), fnu(i), n(i), buffer_list)

                    list(i) = buffer_list
                end do
                !$OMP END PARALLEL DO
            end if

        end subroutine lib_math_array_make_list_of_list_list_cmplx

        subroutine lib_math_get_list_list_structure_cmplx(list, fnu, n)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: list
            integer, intent(inout) :: fnu
            integer, intent(inout) :: n

            fnu = lbound(list%item, 1)
            n = ubound(list%item, 1) - lbound(list%item, 1) + 1
        end subroutine lib_math_get_list_list_structure_cmplx

        subroutine lib_math_get_list_of_list_list_structure_cmplx(list, fnu, n)
            implicit none
            ! dummy
            type(list_list_cmplx), dimension(:), intent(in) :: list
            integer, dimension(:), allocatable, intent(inout) :: fnu
            integer, dimension(:), allocatable, intent(inout) :: n

            ! auxiliary
            integer :: i

            if (allocated(fnu)) deallocate(fnu)
            allocate(fnu(lbound(list, 1):ubound(list, 1)))

            if (allocated(n)) deallocate(n)
            allocate(n(lbound(list, 1):ubound(list, 1)))

            do i = lbound(list, 1), ubound(list, 1)
                call lib_math_get_list_list_structure_cmplx(list(i), fnu(i), n(i))
            end do

        end subroutine lib_math_get_list_of_list_list_structure_cmplx

        function lib_math_get_list_list_cmplx_with_range(list, fnu, n) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: i
            integer :: n_min
            integer :: n_max

            call init_list(rv, fnu, n, dcmplx(0,0))

            n_min = max(fnu, lbound(list%item, 1))
            n_max = min(fnu + n - 1, ubound(list%item, 1))

            !$OMP PARALLEL DO PRIVATE(i)
            do i = n_min, n_max
                rv%item(i) = list%item(i)
            end do
            !$OMP END PARALLEL DO

        end function lib_math_get_list_list_cmplx_with_range

        ! removes zeros of a list
        subroutine lib_math_list_list_cmplx_remove_zeros(list, threshold, set_zero)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(inout) :: list
            complex(kind=lib_math_type_kind), intent(in), optional :: threshold
            logical, intent(in), optional :: set_zero

            ! dummy
            integer :: n
            integer :: m
            integer :: n_min
            integer :: n_max
            type(list_list_cmplx) :: buffer_list

            logical :: m_set_zero

            m_set_zero = .false.
            if (present(threshold)) then
                if (present(set_zero)) m_set_zero = set_zero
            end if

            n_min = -1
            n_max = -1

            ! get boundaries
            do n = lbound(list%item, 1), ubound(list%item, 1)
                do m = -n, n
                    if (present(threshold)) then
                        if( abs(list%item(n)%item(m)) .gt. abs(threshold)) then
                             n_min = n
                             exit
                        end if
                    else
                        if( list%item(n)%item(m) .ne. cmplx(0, 0, kind=lib_math_type_kind) ) then
                             n_min = n
                             exit
                        end if
                    end if
                end do

                if (n_min .ne. -1) then
                    exit
                end if
            end do

            do n = ubound(list%item, 1), lbound(list%item, 1), -1
                do m = -n, n
                    if (present(threshold)) then
                        if( abs(list%item(n)%item(m)) .gt. abs(threshold)) then
                             n_max = n
                             exit
                        end if
                    else
                        if( list%item(n)%item(m) .ne. cmplx(0, 0, kind=lib_math_type_kind) ) then
                             n_max = n
                             exit
                        end if
                    end if
                end do

                if (n_max .ne. -1) then
                    exit
                end if
            end do

            ! create list with removed zeros
            if (n_min .ne. -1 &
                .and. n_max .ne. -1) then

                if (n_min .gt. lbound(list%item, 1) &
                    .or. n_max .lt. ubound(list%item, 1)) then

                    call init_list(buffer_list, n_min, n_max - n_min + 1)

                    do n = n_min, n_max
                        buffer_list%item(n) = list%item(n)
                    end do

                    call move_alloc(buffer_list%item, list%item)
                end if
            else
                ! list contains only zeros
                call init_list(list, 1, 1, cmplx(0, 0, kind=lib_math_type_kind))
            end if

            if (m_set_zero) then
                do n = lbound(list%item, 1), ubound(list%item, 1)
                    do m = -n, n
                        if (abs(real(list%item(n)%item(m))) .le. real(threshold)) then
                            list%item(n)%item(m) = cmplx(0 , aimag(list%item(n)%item(m)), kind=lib_math_type_kind)
                        end if
                        if (abs(aimag(list%item(n)%item(m))) .le. aimag(threshold)) then
                            list%item(n)%item(m) = cmplx(real(list%item(n)%item(m)), 0, kind=lib_math_type_kind)
                        end if
                    end do
                end do
            end if

        end subroutine lib_math_list_list_cmplx_remove_zeros

        ! removes zeros of a list
        subroutine lib_math_list_cmplx_remove_zeros(list, threshold, set_zero)
            implicit none
            ! dummy
            type(list_cmplx), intent(inout) :: list
            complex(kind=lib_math_type_kind), intent(in), optional :: threshold
            logical, intent(in), optional :: set_zero

            ! dummy
            integer :: n
            integer :: n_min
            integer :: n_max
            type(list_cmplx) :: buffer_list
            logical :: m_set_zero

            m_set_zero = .false.
            if (present(threshold)) then
                if (present(set_zero)) m_set_zero = set_zero
            end if

            n_min = -1
            n_max = -1

            ! get boundaries
            do n = lbound(list%item, 1), ubound(list%item, 1)
                if (present(threshold)) then
                    if( abs(list%item(n)) .gt. abs(threshold)) then
                         n_min = n
                         exit
                    end if
                else
                    if( list%item(n) .ne. cmplx(0, 0, kind=lib_math_type_kind) ) then
                         n_min = n
                         exit
                    end if
                end if
            end do

            do n = ubound(list%item, 1), lbound(list%item, 1)
                if (present(threshold)) then
                    if( abs(list%item(n)) .gt. abs(threshold)) then
                         n_max = n
                         exit
                    end if
                else
                    if( list%item(n) .ne. cmplx(0, 0, kind=lib_math_type_kind) ) then
                         n_max = n
                         exit
                    end if
                end if
            end do

            ! create list with removed zeros
            if (n_min .ne. -1 &
                .and. n_max .ne. -1) then

                if (n_min .gt. lbound(list%item, 1) &
                    .or. n_max .lt. ubound(list%item, 1)) then

                    call init_list(buffer_list, n_min, n_max - n_min + 1)

                    do n = n_min, n_max
                        buffer_list%item(n) = list%item(n)
                    end do
                end if
            else
                ! list contains only zeros
                call init_list(list, 1, 1, cmplx(0, 0, kind=lib_math_type_kind))
            end if

            if (m_set_zero) then
                do n = lbound(list%item, 1), lbound(list%item, 1)
                    if (abs(list%item(n)) .le. abs(threshold)) then
                        list%item(n) = cmplx(0 ,0, kind=lib_math_type_kind)
                    end if
                end do
            end if

        end subroutine lib_math_list_cmplx_remove_zeros

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: type(list_list_real)
        !   rhs: type(list_list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_real)
        function lib_math_list_list_real_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_real), intent(in) :: lhs
            type(list_list_real), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: n
            integer :: m

            if (lbound(lhs%item, 1) .eq. lbound(rhs%item, 1) .and. &
                ubound(lhs%item, 1) .eq. ubound(rhs%item, 1) ) then
                call init_list(rv, lbound(lhs%item, 1), ubound(lhs%item, 1) - lbound(lhs%item, 1) + 1)

                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=lbound(lhs%item, 1), ubound(lhs%item, 1)
                    do m=-n, n
                        rv%item(n)%item(m) = lhs%item(n)%item(m) + rhs%item(n)%item(m)
                    end do
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_list_cmplx_add: ERROR"
                print *, "  size of lhs and rhs are not equal"
            end if

        end function lib_math_list_list_real_add

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: real, dimension(:), allocatable
        !   rhs: real, dimension(:), allocatable
        !
        ! Retruns
        ! ----
        !   rv: real, dimension(:), allocatable
        function lib_math_real_array_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: rhs

            real(kind=lib_math_type_kind), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs) .and. allocated(rhs)) then
                n_range_mutual(1) = max( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_mutual(2) = min( ubound(lhs, 1), ubound(rhs, 1) )

                n_range_full(1) = min( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_full(2) = max( ubound(lhs, 1), ubound(rhs, 1) )

                allocate (rv(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv(i) = lhs(i) + rhs(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs, 1) .lt. lbound(rhs, 1)) then
                    do i = lbound(lhs, 1), lbound(rhs, 1) - 1
                        rv(i) = lhs(i)
                    end do
                end if

                if (lbound(rhs, 1) .lt. lbound(lhs, 1)) then
                    do i = lbound(rhs, 1), lbound(lhs, 1) - 1
                        rv(i) = rhs(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs, 1) .gt. ubound(rhs, 1)) then
                    do i = ubound(rhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = lhs(i)
                    end do
                end if

                if (ubound(rhs, 1) .gt. ubound(lhs, 1)) then
                    do i = ubound(lhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = rhs(i)
                    end do
                end if

            else if (allocated(lhs)) then
                rv = lhs
            else if (allocated(rhs)) then
                rv = rhs
            else
                print *, "lib_math_real_array_add: ERROR"
            end if

        end function lib_math_real_array_add

        ! Elementwise subtraction
        !
        ! Arguments
        ! ----
        !   lhs: real, dimension(:), allocatable
        !   rhs: real, dimension(:), allocatable
        !
        ! Retruns
        ! ----
        !   rv: real, dimension(:), allocatable
        function lib_math_real_array_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: rhs

            real(kind=lib_math_type_kind), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs) .and. allocated(rhs)) then
                n_range_mutual(1) = max( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_mutual(2) = min( ubound(lhs, 1), ubound(rhs, 1) )

                n_range_full(1) = min( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_full(2) = max( ubound(lhs, 1), ubound(rhs, 1) )

                allocate (rv(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv(i) = lhs(i) - rhs(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs, 1) .lt. lbound(rhs, 1)) then
                    do i = lbound(lhs, 1), lbound(rhs, 1) - 1
                        rv(i) = lhs(i)
                    end do
                end if

                if (lbound(rhs, 1) .lt. lbound(lhs, 1)) then
                    do i = lbound(rhs, 1), lbound(lhs, 1) - 1
                        rv(i) = -1 * rhs(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs, 1) .gt. ubound(rhs, 1)) then
                    do i = ubound(rhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = lhs(i)
                    end do
                end if

                if (ubound(rhs, 1) .gt. ubound(lhs, 1)) then
                    do i = ubound(lhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = -1 * rhs(i)
                    end do
                end if

            else if (allocated(lhs)) then
                rv = lhs
            else if (allocated(rhs)) then
                rv = -1 * rhs
            else
                print *, "lib_math_real_array_sub: ERROR"
            end if

        end function lib_math_real_array_sub

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: double complex, dimension(:), allocatable
        !   rhs: double complex, dimension(:), allocatable
        !
        ! Retruns
        ! ----
        !   rv: double complex, dimension(:), allocatable
        function lib_math_cmplx_array_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: rhs

            complex(kind=lib_math_type_kind), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs) .and. allocated(rhs)) then
                n_range_mutual(1) = max( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_mutual(2) = min( ubound(lhs, 1), ubound(rhs, 1) )

                n_range_full(1) = min( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_full(2) = max( ubound(lhs, 1), ubound(rhs, 1) )

                allocate (rv(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv(i) = lhs(i) + rhs(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs, 1) .lt. lbound(rhs, 1)) then
                    do i = lbound(lhs, 1), lbound(rhs, 1) - 1
                        rv(i) = lhs(i)
                    end do
                end if

                if (lbound(rhs, 1) .lt. lbound(lhs, 1)) then
                    do i = lbound(rhs, 1), lbound(lhs, 1) - 1
                        rv(i) = rhs(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs, 1) .gt. ubound(rhs, 1)) then
                    do i = ubound(rhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = lhs(i)
                    end do
                end if

                if (ubound(rhs, 1) .gt. ubound(lhs, 1)) then
                    do i = ubound(lhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = rhs(i)
                    end do
                end if

            else if (allocated(lhs)) then
                rv = lhs
            else if (allocated(rhs)) then
                rv = rhs
            else
                print *, "lib_math_cmplx_array_add: ERROR"
            end if

        end function lib_math_cmplx_array_add

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: double complex, dimension(:), allocatable
        !   rhs: double complex, dimension(:), allocatable
        !
        ! Retruns
        ! ----
        !   rv: double complex, dimension(:), allocatable
        function lib_math_cmplx_array_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), allocatable, intent(in) :: rhs

            complex(kind=lib_math_type_kind), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs) .and. allocated(rhs)) then
                n_range_mutual(1) = max( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_mutual(2) = min( ubound(lhs, 1), ubound(rhs, 1) )

                n_range_full(1) = min( lbound(lhs, 1), lbound(rhs, 1) )
                n_range_full(2) = max( ubound(lhs, 1), ubound(rhs, 1) )

                allocate (rv(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv(i) = lhs(i) - rhs(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs, 1) .lt. lbound(rhs, 1)) then
                    do i = lbound(lhs, 1), lbound(rhs, 1) - 1
                        rv(i) = lhs(i)
                    end do
                end if

                if (lbound(rhs, 1) .lt. lbound(lhs, 1)) then
                    do i = lbound(rhs, 1), lbound(lhs, 1) - 1
                        rv(i) = -1 * rhs(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs, 1) .gt. ubound(rhs, 1)) then
                    do i = ubound(rhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = lhs(i)
                    end do
                end if

                if (ubound(rhs, 1) .gt. ubound(lhs, 1)) then
                    do i = ubound(lhs, 1) + 1, ubound(rhs, 1)
                        rv(i) = -1 * rhs(i)
                    end do
                end if

            else if (allocated(lhs)) then
                rv = lhs
            else if (allocated(rhs)) then
                rv = -1 * rhs
            else
                print *, "lib_math_cmplx_array_sub: ERROR"
            end if

        end function lib_math_cmplx_array_sub

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: type(list_cmplx)
        !   rhs: type(list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_cmplx)
        function lib_math_list_cmplx_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_cmplx), intent(in) :: lhs
            type(list_cmplx), intent(in) :: rhs

            type(list_cmplx) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                allocate (rv%item(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) + rhs%item(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = rhs%item(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    do i = ubound(rhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = rhs%item(i)
                    end do
                end if

            else if (allocated(lhs%item)) then
                rv%item = lhs%item
            else if (allocated(rhs%item)) then
                rv%item = rhs%item
            else
                print *, "lib_math_list_cmplx_add: ERROR"
            end if

        end function lib_math_list_cmplx_add

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: type(list_real)
        !   rhs: type(list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_real)
        function lib_math_list_real_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_real), intent(in) :: lhs
            type(list_real), intent(in) :: rhs

            type(list_real) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                allocate (rv%item(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) + rhs%item(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = rhs%item(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    do i = ubound(rhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = rhs%item(i)
                    end do
                end if

            else if (allocated(lhs%item)) then
                rv%item = lhs%item
            else if (allocated(rhs%item)) then
                rv%item = rhs%item
            else
                print *, "lib_math_list_real_add: ERROR"
            end if

        end function lib_math_list_real_add

        ! Elementwise subtraction
        !
        ! Arguments
        ! ----
        !   lhs: type(list_cmplx)
        !   rhs: type(list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_cmplx)
        function lib_math_list_cmplx_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_cmplx), intent(in) :: lhs
            type(list_cmplx), intent(in) :: rhs

            type(list_cmplx) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                allocate (rv%item(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) - rhs%item(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = -1 * rhs%item(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    do i = ubound(rhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = -1 * rhs%item(i)
                    end do
                end if

            else if (allocated(lhs%item)) then
                rv%item = lhs%item
            else if (allocated(rhs%item)) then
                rv%item = -1 * rhs%item
            else
                print *, "lib_math_list_cmplx_sub: ERROR"
            end if

        end function lib_math_list_cmplx_sub

        ! Elementwise subtraction
        !
        ! Arguments
        ! ----
        !   lhs: type(list_real)
        !   rhs: type(list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_real)
        function lib_math_list_real_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_real), intent(in) :: lhs
            type(list_real), intent(in) :: rhs

            type(list_real) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                allocate (rv%item(n_range_full(1):n_range_full(2)))

                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) - rhs%item(i)
                end do

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = -1 * rhs%item(i)
                    end do
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    do i = ubound(rhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = -1 * rhs%item(i)
                    end do
                end if

            else if (allocated(lhs%item)) then
                rv%item = lhs%item
            else if (allocated(rhs%item)) then
                rv%item = -1 * rhs%item
            else
                print *, "lib_math_list_real_sub: ERROR"
            end if

        end function lib_math_list_real_sub

        ! Elementwise addition
        !
        ! Arguments
        ! ----
        !   lhs: type(list_list_cmplx)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_list_list_cmplx_add(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                !call init_list(rv, n_range_full(1), n_range_full(2) - n_range_full(1) + 1)
                allocate(rv%item(n_range_full(1):n_range_full(2)))
                !$OMP PARALLEL DO PRIVATE(i)
                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) + rhs%item(i)
                end do
                !$OMP END PARALLEL DO

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = rhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = ubound(rhs%item, 1) + 1, ubound(lhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = rhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

            else if (allocated(lhs%item)) then
                rv = lhs
            else if (allocated(rhs%item)) then
                rv = rhs
            else
                print*, "lib_math_list_list_cmplx_add: ERROR"
            end if

        end function lib_math_list_list_cmplx_add

        ! Elementwise subtraction
        !
        ! Arguments
        ! ----
        !   lhs: type(list_list_real)
        !   rhs: type(list_list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_real)
        function lib_math_list_list_real_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_real), intent(in) :: lhs
            type(list_list_real), intent(in) :: rhs

            type(list_list_real) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2) :: n_range_mutual
            integer, dimension(2) :: n_range_full

            real(kind=lib_math_type_kind) :: minus_1

            minus_1 = real(-1, kind=lib_math_type_kind)

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                n_range_mutual(1) = max( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_mutual(2) = min( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                n_range_full(1) = min( lbound(lhs%item, 1), lbound(rhs%item, 1) )
                n_range_full(2) = max( ubound(lhs%item, 1), ubound(rhs%item, 1) )

                !call init_list(rv, n_range_full(1), n_range_full(2) - n_range_full(1) + 1)
                allocate(rv%item(n_range_full(1):n_range_full(2)))
                !$OMP PARALLEL DO PRIVATE(i)
                do i = n_range_mutual(1), n_range_mutual(2)
                    rv%item(i) = lhs%item(i) - rhs%item(i)
                end do
                !$OMP END PARALLEL DO

                ! fill entries with i < n_range_mutual(1)
                if (lbound(lhs%item, 1) .lt. lbound(rhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = lbound(lhs%item, 1), lbound(rhs%item, 1) - 1
                        rv%item(i) = lhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                if (lbound(rhs%item, 1) .lt. lbound(lhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = lbound(rhs%item, 1), lbound(lhs%item, 1) - 1
                        rv%item(i) = minus_1 * rhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                ! fill entries with i > n_range_mutual(2)
                if (ubound(lhs%item, 1) .gt. ubound(rhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = ubound(rhs%item, 1) + 1, ubound(lhs%item, 1)
                        rv%item(i) = lhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

                if (ubound(rhs%item, 1) .gt. ubound(lhs%item, 1)) then
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = ubound(lhs%item, 1) + 1, ubound(rhs%item, 1)
                        rv%item(i) = minus_1 * rhs%item(i)
                    end do
                    !$OMP END PARALLEL DO
                end if

            else if (allocated(lhs%item)) then
                rv = lhs
            else if (allocated(rhs%item)) then
                rv = minus_1 * rhs
            else
                print*, "lib_math_list_list_cmplx_sub: ERROR"
            end if

        end function lib_math_list_list_real_sub

        ! Elementwise subtraction
        !
        ! Arguments
        ! ----
        !   lhs: type(list_list_cmplx)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_list_list_cmplx_sub(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
!            integer :: n
!            integer :: m

            rv = lib_math_list_list_cmplx_add(lhs, lib_math_list_list_cmplx_mul_real(real(-1, kind=lib_math_type_kind), rhs))

!            if (lbound(lhs%item, 1) .eq. lbound(rhs%item, 1) .and. &
!                ubound(lhs%item, 1) .eq. ubound(rhs%item, 1) ) then
!
!                call init_list(rv, lbound(lhs%item, 1), ubound(lhs%item, 1) - lbound(lhs%item, 1) + 1)
!
!                !$OMP PARALLEL DO PRIVATE(n, m)
!                do n=lbound(lhs%item, 1), ubound(lhs%item, 1)
!                    do m=-n, n
!                        rv%item(n)%item(m) = lhs%item(n)%item(m) - rhs%item(n)%item(m)
!                    end do
!                end do
!                !$OMP END PARALLEL DO
!            else
!                print *, "lib_math_list_list_cmplx_sub: ERROR"
!                print *, "  size of lhs and rhs are not equal"
!            end if
        end function lib_math_list_list_cmplx_sub

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: real(kind=lib_math_type_kind)
        !   rhs: type(list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_real)
        function lib_math_list_real_mul_real(lhs, rhs) result (rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type(list_real), intent(in) :: rhs

            type(list_real) :: rv

            ! auxiliary
            integer :: i

            call init_list(rv, lbound(rhs%item, 1), ubound(rhs%item, 1) - lbound(rhs%item, 1) + 1)

            do i=lbound(rhs%item, 1), ubound(rhs%item, 1)
                rv%item(i) = lhs * rhs%item(i)
            end do
        end function lib_math_list_real_mul_real

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: cmplx(kind=lib_math_type_kind)
        !   rhs: type(list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_cmplx)
        function lib_math_list_real_mul_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type(list_real), intent(in) :: rhs

            type(list_cmplx) :: rv

            ! auxiliary
            integer :: i

            call init_list(rv, lbound(rhs%item, 1), ubound(rhs%item, 1) - lbound(rhs%item, 1) + 1)

            do i=lbound(rhs%item, 1), ubound(rhs%item, 1)
                rv%item(i) = lhs * rhs%item(i)
            end do
        end function lib_math_list_real_mul_cmplx

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: real(kind=lib_math_type_kind)
        !   rhs: type(list_list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_real)
        function lib_math_list_list_real_mul_real(lhs, rhs) result (rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type(list_list_real), intent(in) :: rhs

            type(list_list_real) :: rv

            ! auxiliary
            integer :: i

            call init_list(rv, lbound(rhs%item, 1), ubound(rhs%item, 1) - lbound(rhs%item, 1) + 1)

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%item, 1), ubound(rhs%item, 1)
                rv%item(i)%item = lhs * rhs%item(i)%item
            end do
            !$OMP END PARALLEL DO
        end function lib_math_list_list_real_mul_real

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: real(kind=lib_math_type_kind)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_real_mul_list_list_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            rv = rhs * lhs

        end function lib_math_real_mul_list_list_cmplx

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: real(kind=lib_math_type_kind)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_list_list_cmplx_mul_real(lhs, rhs) result (rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: n
            integer :: m

            call init_list(rv, rhs)

            if (lhs .eq. real(1, kind = lib_math_type_kind)) then
                rv = rhs
            else if (lhs .eq. real(0, kind = lib_math_type_kind)) then
                call deallocate_list(rv)
            else
                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=lbound(rhs%item, 1), ubound(rhs%item, 1)
                    do m=lbound(rhs%item(n)%item, 1), ubound(rhs%item(n)%item, 1)
                        rv%item(n)%item(m) = lhs * rhs%item(n)%item(m)
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

        end function lib_math_list_list_cmplx_mul_real

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: complex(kind=lib_math_type_kind)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_cmplx_mul_list_list_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: n
            integer :: m

            call init_list(rv, rhs)

            if (lhs .eq. cmplx(1, 0, kind = lib_math_type_kind)) then
                rv = rhs
            else if (lhs .eq. cmplx(0, 0, kind = lib_math_type_kind)) then
                call deallocate_list(rv)
            else
                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=lbound(rhs%item, 1), ubound(rhs%item, 1)
                    do m=lbound(rhs%item(n)%item, 1), ubound(rhs%item(n)%item, 1)
                        rv%item(n)%item(m) = lhs * rhs%item(n)%item(m)
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

        end function lib_math_cmplx_mul_list_list_cmplx

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: type(list_list_cmplx)
        !   rhs: complex(kind=lib_math_type_kind)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_list_list_cmplx_mul_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            rv = lib_math_cmplx_mul_list_list_cmplx(rhs, lhs)

        end function lib_math_list_list_cmplx_mul_cmplx

        ! Elementwise multiplication
        !
        ! HINT
        ! ----
        !   An uallocated list means that there are only zeros.
        !
        ! Arguments
        ! ----
        !   lhs: type(list_cmplx)
        !   rhs: type(list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_cmplx)
        !
        ! Note
        ! ----
        !   lhs:   l-----u
        !   rhs:       l-------u
        !
        !   lhs:       l-------u
        !   rhs:   l-----u
        !
        !
        !   lhs:   l-----u
        !   rhs:     l-u
        !
        !
        !   lhs:     l-u
        !   rhs:   l-----u
        !
        !
        !   overlap: l.u <= r.l  or  l.l >= r.u
        function lib_math_list_list_cmplx_mul_list_list_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: n

            integer :: n_min
            integer :: n_max

            logical :: calculate

            calculate = .false.
            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                if (lbound(lhs%item, 1) .ge. ubound(rhs%item, 1) &
                    .or. ubound(lhs%item, 1) .le. lbound(rhs%item, 1)) then
                    calculate = .true.
                else if (lbound(lhs%item, 1) .ge. lbound(rhs%item, 1) &
                         .and. ubound(lhs%item, 1) .le. ubound(rhs%item, 1)) then
                    calculate = .true.
                else if (lbound(lhs%item, 1) .le. lbound(rhs%item, 1) &
                         .and. ubound(lhs%item, 1) .ge. ubound(rhs%item, 1)) then
                    calculate = .true.
                end if

                if (calculate) then
                    n_min = max(lbound(lhs%item, 1), lbound(rhs%item, 1))
                    n_max = min(ubound(lhs%item, 1), ubound(rhs%item, 1))

                    allocate(rv%item(n_min:n_max))

                    !$OMP PARALLEL DO PRIVATE(n)
                    do n=lbound(rv%item, 1), ubound(rv%item, 1)
                        rv%item(n) = lhs%item(n) * rhs%item(n)
                    end do
                    !$OMP END PARALLEL DO
                end if
            end if

        end function lib_math_list_list_cmplx_mul_list_list_cmplx

        ! Elementwise multiplication
        !
        ! Arguments
        ! ----
        !   lhs: type(list_cmplx)
        !   rhs: type(list_list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_list_cmplx)
        function lib_math_list_cmplx_mul_list_list_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_cmplx), intent(in) :: lhs
            type(list_list_cmplx), intent(in) :: rhs

            type(list_list_cmplx) :: rv

            ! auxiliary
            integer :: n
            integer :: m

            integer :: n_min
            integer :: n_max

            if (lbound(lhs%item, 1) .eq. lbound(rhs%item, 1) &
                .and. ubound(lhs%item, 1) .eq. ubound(rhs%item, 1)) then

                n_min = lbound(lhs%item, 1)
                n_max = ubound(lhs%item, 1)

                call init_list(rv, n_min, n_max - n_min + 1)


                !$OMP PARALLEL DO PRIVATE(n, m)
                do n=n_min, n_max
                    do m=lbound(rhs%item(n)%item, 1), ubound(rhs%item(n)%item, 1)
                        !rv%item(n)%item(m) = lhs%item(n) * rhs%item(n)%item(m)
                        rv%item(n)%item = lhs%item(n) * rhs%item(n)%item
                    end do
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_cmplx_mul_list_list_cmplx: ERROR"
            end if

        end function lib_math_list_cmplx_mul_list_list_cmplx

        ! Elementwise multiplication
        !
        ! HINT
        ! ----
        !   An uallocated list means that there are only zeros.
        !
        ! Arguments
        ! ----
        !   lhs: type(list_real)
        !   rhs: type(list_real)
        !
        ! Retruns
        ! ----
        !   rv: type(list_real)
        !
        ! Note
        ! ----
        !   lhs:   l-----u
        !   rhs:       l-------u
        !
        !   lhs:       l-------u
        !   rhs:   l-----u
        !
        !
        !   lhs:   l-----u
        !   rhs:     l-u
        !
        !
        !   lhs:     l-u
        !   rhs:   l-----u
        !
        !
        !   overlap: l.u <= r.l  or  l.l >= r.u
        function lib_math_list_real_mul_list_real(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_real), intent(in) :: lhs
            type(list_real), intent(in) :: rhs

            type(list_real) :: rv

            ! auxiliary
            integer :: n

            integer :: n_min
            integer :: n_max

            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                if (lbound(lhs%item, 1) .ge. ubound(rhs%item, 1) &
                    .or. ubound(lhs%item, 1) .le. lbound(rhs%item, 1)) then

                    n_min = max(lbound(lhs%item, 1), lbound(rhs%item, 1))
                    n_max = min(ubound(lhs%item, 1), ubound(rhs%item, 1))

                    allocate(rv%item(n_min:n_max))

                    !$OMP PARALLEL DO PRIVATE(n)
                    do n=lbound(rv%item, 1), ubound(rv%item, 1)
                        rv%item(n) = lhs%item(n) * rhs%item(n)
                    end do
                    !$OMP END PARALLEL DO
                end if
            end if

        end function lib_math_list_real_mul_list_real

        ! Elementwise multiplication
        !
        ! HINT
        ! ----
        !   An uallocated list means that there are only zeros.
        !
        ! Arguments
        ! ----
        !   lhs: type(list_cmplx)
        !   rhs: type(list_cmplx)
        !
        ! Retruns
        ! ----
        !   rv: type(list_cmplx)
        !
        ! Note
        ! ----
        !   lhs:   l-----u
        !   rhs:       l-------u
        !
        !   lhs:       l-------u
        !   rhs:   l-----u
        !
        !
        !   lhs:   l-----u
        !   rhs:     l-u
        !
        !
        !   lhs:     l-u
        !   rhs:   l-----u
        !
        !
        !   overlap: l.u <= r.l  or  l.l >= r.u
        function lib_math_list_cmplx_mul_list_cmplx(lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(list_cmplx), intent(in) :: lhs
            type(list_cmplx), intent(in) :: rhs

            type(list_cmplx) :: rv

            ! auxiliary
            integer :: n

            integer :: n_min
            integer :: n_max

            logical :: calculate

            calculate = .false.
            if (allocated(lhs%item) .and. allocated(rhs%item)) then
                if (lbound(lhs%item, 1) .ge. ubound(rhs%item, 1) &
                    .or. ubound(lhs%item, 1) .le. lbound(rhs%item, 1)) then
                    calculate = .true.
                else if (lbound(lhs%item, 1) .ge. lbound(rhs%item, 1) &
                         .and. ubound(lhs%item, 1) .le. ubound(rhs%item, 1)) then
                    calculate = .true.
                else if (lbound(lhs%item, 1) .le. lbound(rhs%item, 1) &
                         .and. ubound(lhs%item, 1) .ge. ubound(rhs%item, 1)) then
                    calculate = .true.
                end if

                if (calculate) then
                    n_min = max(lbound(lhs%item, 1), lbound(rhs%item, 1))
                    n_max = min(ubound(lhs%item, 1), ubound(rhs%item, 1))

                    allocate(rv%item(n_min:n_max))

                    !$OMP PARALLEL DO PRIVATE(n)
                    do n=lbound(rv%item, 1), ubound(rv%item, 1)
                        rv%item(n) = lhs%item(n) * rhs%item(n)
                    end do
                    !$OMP END PARALLEL DO
                end if
            end if

        end function lib_math_list_cmplx_mul_list_cmplx

        ! Arguments
        ! ----
        !   list: type (list_spherical_coordinate_cmplx_type)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_spherical_coordinate_cmplx_type_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            ! auxiliary
            integer :: i

            if( allocated(list) ) then
                deallocate( list )
            end if

            allocate( list(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list(i)%coordinate(-i:i))
            end do

        end subroutine lib_math_list_spherical_coordinate_cmplx_type_init

        ! reference: http://mathworld.wolfram.com/SphericalCoordinates.html
        ! ISO 31-11
        subroutine lib_math_spherical_point_to_cartesian_point(lhs, rhs)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type), intent(inout) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: theta
            real(kind=lib_math_type_kind) :: phi

            real(kind=lib_math_type_kind) :: cos_theta
            real(kind=lib_math_type_kind) :: sin_theta
            real(kind=lib_math_type_kind) :: cos_phi
            real(kind=lib_math_type_kind) :: sin_phi

            real(kind=lib_math_type_kind) :: r_sin_theta

            theta = rhs%theta
            phi = rhs%phi

            if (abs(theta) .gt. 2d0 * PI) then
                theta = modulo(theta, 2 * PI)
            end if

            if (theta .lt. 0d0) then
                theta = 2 * PI + theta
            end if

            if (theta .gt. PI) then
                theta = 2 * PI - theta
                phi = PI + phi
            end if

            if (phi .ge. 2 * PI) then
                phi = modulo(phi, 2 * PI)
            end if

            cos_theta = cos(theta)
            sin_theta = sin(theta)
            cos_phi = cos(phi)
            sin_phi = sin(phi)

            r_sin_theta = rhs%rho * sin_theta

            lhs%x = r_sin_theta * cos_phi
            lhs%y = r_sin_theta * sin_phi
            lhs%z = rhs%rho * cos_theta

        end subroutine lib_math_spherical_point_to_cartesian_point

        ! reference: http://mathworld.wolfram.com/SphericalCoordinates.html
        ! ISO 31-11
        subroutine lib_math_cartesian_point_to_spherical_point(lhs, rhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(spherical_coordinate_real_type), intent(inout) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: r

            r = sqrt(rhs%x * rhs%x + rhs%y * rhs%y + rhs%z * rhs%z)

            lhs%rho = r

            if (r .eq. 0) then
                lhs%theta = 0
            else
                lhs%theta = acos(rhs%z / r)
            end if

            if (rhs%x .eq. 0 .and. rhs%y .eq. 0) then
                lhs%phi = 0
            else
                lhs%phi = atan(rhs%y, rhs%x)
            end if

        end subroutine lib_math_cartesian_point_to_spherical_point

        subroutine lib_math_array_to_cartesian_point(lhs, rhs)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(3), intent(in) :: rhs
            type(cartesian_coordinate_real_type), intent(inout) :: lhs

            lhs%x = rhs(1)
            lhs%y = rhs(2)
            lhs%z = rhs(3)

        end subroutine

        subroutine lib_math_cartesian_point_to_array(lhs, rhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: rhs
            real(kind=lib_math_type_kind), dimension(3), intent(inout) :: lhs

            lhs(1) = rhs%x
            lhs(2) = rhs%y
            lhs(3) = rhs%z

        end subroutine

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   theta: real
        !       polar coordinate of the complex vector (rhs)
        !   phi: real
        !       azimuthal coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, theta, phi) result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            real(kind=lib_math_type_kind), intent(in) :: theta
            real(kind=lib_math_type_kind), intent(in) :: phi

            type(cartesian_coordinate_cmplx_type) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: cos_theta
            real(kind=lib_math_type_kind) :: sin_theta
            real(kind=lib_math_type_kind) :: cos_phi
            real(kind=lib_math_type_kind) :: sin_phi

            cos_theta = cos(theta)
            sin_theta = sin(theta)
            cos_phi = cos(phi)
            sin_phi = sin(phi)

            lhs%x = rhs%rho * (sin_theta * cos_phi) &
                    + rhs%theta * (cos_theta * cos_phi) &
                    - rhs%phi * sin_phi
            lhs%y = rhs%rho * (sin_theta * sin_phi) &
                      + rhs%theta * (cos_theta * sin_phi) &
                      + rhs%phi * cos_phi
            lhs%z = rhs%rho * cos_theta &
                    - rhs%theta * sin_theta

        end function lib_math_spherical_components_to_cartesian_components_cmplx_a

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   coordinate: type(cartesian_coordinate_real_type)
        !       cartesian coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_c(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            type(cartesian_coordinate_real_type), intent(in) :: coordinate

            type(cartesian_coordinate_cmplx_type) :: lhs

            ! auxiliary
            type(spherical_coordinate_real_type) :: coordinate_spherical

            coordinate_spherical = coordinate

            lhs = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                coordinate_spherical%theta, &
                                                                                coordinate_spherical%phi)
        end function lib_math_spherical_components_to_cartesian_components_cmplx_c

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   coordinate: type(spherical_coordinate_real_type)
        !       spherical coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_s(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            type(spherical_coordinate_real_type), intent(in) :: coordinate

            type(cartesian_coordinate_cmplx_type) :: lhs

            lhs = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                    coordinate%theta, &
                                                                                    coordinate%phi)
        end function lib_math_spherical_components_to_cartesian_components_cmplx_s

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   theta: real
        !       polar coordinate of the complex vector (rhs)
        !   phi: real
        !       azimuthal coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, theta, phi) result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            real(kind=lib_math_type_kind), intent(in) :: theta
            real(kind=lib_math_type_kind), intent(in) :: phi

            type(spherical_coordinate_cmplx_type) :: lhs

            ! auxiliary
            complex(kind=lib_math_type_kind) :: cos_theta
            complex(kind=lib_math_type_kind) :: sin_theta
            complex(kind=lib_math_type_kind) :: cos_phi
            complex(kind=lib_math_type_kind) :: sin_phi

            cos_theta = cos(theta)
            sin_theta = sin(theta)
            cos_phi = cos(phi)
            sin_phi = sin(phi)

            lhs%rho = rhs%x * (sin_theta *cos_phi) &
                    + rhs%y * (sin_theta * sin_phi) &
                    + rhs%z * cos_theta
            lhs%theta = rhs%x * (cos_theta * cos_phi) &
                      + rhs%y * (cos_theta * sin_phi) &
                      - rhs%z * sin_theta
            lhs%phi = - rhs%x * sin_phi &
                      + rhs%y * cos_phi

        end function lib_math_cartesian_components_to_spherical_components_cmplx_a

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   coordinate: type(cartesian_coordinate_real_type)
        !       coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_c(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type(cartesian_coordinate_real_type), intent(in) :: coordinate

            type(spherical_coordinate_cmplx_type) :: lhs

            ! auxiliary
            type(spherical_coordinate_real_type) :: coordinate_spherical

            coordinate_spherical = coordinate

            lhs = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                coordinate_spherical%theta, &
                                                                                coordinate_spherical%phi)
        end function lib_math_cartesian_components_to_spherical_components_cmplx_c

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   coordinate: type(spherical_coordinate_real_type)
        !       coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_s(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type(spherical_coordinate_real_type), intent(in) :: coordinate

            type(spherical_coordinate_cmplx_type) :: lhs

            lhs = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                coordinate%theta, &
                                                                                coordinate%phi)
        end function lib_math_cartesian_components_to_spherical_components_cmplx_s

        function lib_math_make_cartesian_coordiante_real(x, y, z) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind) :: x
            real(kind=lib_math_type_kind) :: y
            real(kind=lib_math_type_kind) :: z

            type(cartesian_coordinate_real_type) :: rv

            rv%x = x
            rv%y = y
            rv%z = z

        end function lib_math_make_cartesian_coordiante_real

        function lib_math_make_cartesian_coordiante_cmplx(x, y, z) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind) :: x
            complex(kind=lib_math_type_kind) :: y
            complex(kind=lib_math_type_kind) :: z

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = x
            rv%y = y
            rv%z = z

        end function lib_math_make_cartesian_coordiante_cmplx

        ! Argument
        ! ----
        !   rho: real(kind=lib_math_type_kind)
        !       vector length
        !   theta: real(kind=lib_math_type_kind)
        !       polar angle [rad]
        !   phi: real(kind=lib_math_type_kind)
        !       azimuthal angle [rad]
        function lib_math_make_spherical_coordiante(rho, theta, phi) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind) :: rho
            real(kind=lib_math_type_kind) :: theta
            real(kind=lib_math_type_kind) :: phi

            type(spherical_coordinate_real_type) :: rv

            rv%theta = theta
            rv%phi = phi

            if (abs(rv%theta) .gt. 2d0 * PI) then
                rv%theta = modulo(rv%theta, 2 * PI)
            end if

            if (rv%theta .lt. 0d0) then
                rv%theta = 2 * PI + rv%theta
            end if

            if (rv%theta .gt. PI) then
                rv%theta = 2 * PI - rv%theta
                rv%phi = PI + rv%phi
            end if

            if (rv%phi .ge. 2 * PI) then
                rv%phi = modulo(rv%phi, 2 * PI)
            end if

            rv%rho = rho
!            rv%theta = theta
!            rv%phi = phi

        end function lib_math_make_spherical_coordiante

        function lib_math_cartesian_cross_product_real (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: lhs
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type) :: rv

            rv%x = lhs%y * rhs%z - lhs%z * rhs%y
            rv%y = lhs%z * rhs%x - lhs%x * rhs%z
            rv%z = lhs%x * rhs%y - lhs%y * rhs%x

        end function

        function lib_math_cartesian_cross_product_real_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_real_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            type(cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_cross_product_real(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_cross_product_cmplx (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%y * rhs%z - lhs%z * rhs%y
            rv%y = lhs%z * rhs%x - lhs%x * rhs%z
            rv%z = lhs%x * rhs%y - lhs%y * rhs%x

        end function

        function lib_math_cartesian_cross_product_cmplx_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_cross_product_cmplx(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_dot_product_real (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: lhs
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            real(kind=lib_math_type_kind) :: rv

            rv = lhs%x * rhs%x + lhs%y * rhs%y + lhs%z * rhs%z

        end function

        function lib_math_cartesian_dot_product_real_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_real_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_dot_product_real(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_dot_product_cmplx (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=lib_math_type_kind) :: rv

            rv = lhs%x * rhs%x + lhs%y * rhs%y + lhs%z * rhs%z

        end function

        function lib_math_cartesian_dot_product_cmplx_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_dot_product_cmplx(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function list_filter_cartesian_coordiante_real(list, use_point) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: list
            logical, dimension(lbound(list, 1):ubound(list, 1)), intent(in) :: use_point

            type(cartesian_coordinate_real_type), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i
            integer :: x

            integer, dimension(2):: i_range

            i_range(1) = lbound(list, 1)
            i_range(2) = ubound(list, 1)

            x = count(use_point)
            allocate(rv(x))

            x = 0
            do i = i_range(1), i_range(2)
                if (use_point(i)) then
                    x = x + 1
                    rv(x) = list(i)
                end if
            end do

        end function list_filter_cartesian_coordiante_real

        function list_filter_cartesian_coordiante_cmplx(list, use_point) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: list
            logical, dimension(lbound(list, 1):ubound(list, 1)), intent(in) :: use_point

            type(cartesian_coordinate_cmplx_type), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i
            integer :: x

            integer, dimension(2):: i_range

            i_range(1) = lbound(list, 1)
            i_range(2) = ubound(list, 1)

            x = count(use_point)
            allocate(rv(x))

            x = 0
            do i = i_range(1), i_range(2)
                if (use_point(i)) then
                    x = x + 1
                    rv(x) = list(i)
                end if
            end do

        end function list_filter_cartesian_coordiante_cmplx

        function list_filter_spherical_coordiante_real(list, use_point) result(rv)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), dimension(:), intent(in) :: list
            logical, dimension(lbound(list, 1):ubound(list, 1)), intent(in) :: use_point

            type(spherical_coordinate_real_type), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i
            integer :: x

            integer, dimension(2):: i_range

            i_range(1) = lbound(list, 1)
            i_range(2) = ubound(list, 1)

            x = count(use_point)
            allocate(rv(x))

            x = 0
            do i = i_range(1), i_range(2)
                if (use_point(i)) then
                    x = x + 1
                    rv(x) = list(i)
                end if
            end do

        end function list_filter_spherical_coordiante_real

        function list_filter_spherical_coordiante_cmplx(list, use_point) result(rv)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), dimension(:), intent(in) :: list
            logical, dimension(lbound(list, 1):ubound(list, 1)), intent(in) :: use_point

            type(spherical_coordinate_cmplx_type), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i
            integer :: x

            integer, dimension(2):: i_range

            i_range(1) = lbound(list, 1)
            i_range(2) = ubound(list, 1)

            x = count(use_point)
            allocate(rv(x))

            x = 0
            do i = i_range(1), i_range(2)
                if (use_point(i)) then
                    x = x + 1
                    rv(x) = list(i)
                end if
            end do

        end function list_filter_spherical_coordiante_cmplx

        function lib_math_type_operator_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! parameter
            double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

            rv = 0

            if (.not. test_lib_math_list_spherical_operator_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_sub_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_0d_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_add_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_sub_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_real_mul_array_c()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_c_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_real_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_c_array_divide_by_r_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_c_array_divide_by_c_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_real()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_cmplx()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_list_cartesian_operator_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_sub_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_0d_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_add_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_sub_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_real_mul_array_c()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_real_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_c_array_divide_by_r_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_c_array_divide_by_c_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_divide_by_real()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_divide_by_cmplx()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_cartesian_point_to_spherical_point()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_spherical_point_to_cartesian_point()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_spherical_point_to_cartesian_point_2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_spherical_components_to_cartesian_components_c_a()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_components_to_spherical_components_c_a()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_cross_product_real_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_cross_product_cmplx_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_dot_product_real_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_dot_product_cmplx_list()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_list_cmplx_add()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cmplx_mul()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_list_list_cmplx_add()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_list_cmplx_mul()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_list_of_list_list_make_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_array_make_list_of_list_list_cmplx()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_get_matrix_rot()) rv = rv + 1

            print *, "--------lib_math_type_operator_test_functions--------"
            if (rv == 0) then
                print *, "lib_math_type_operator_test_functions tests: OK"
            else
                print *, rv,"lib_math_type_operator_test_functions test(s) FAILED"
            end if
            print *, "-----------------------------------------------------"

            contains

                function get_test_values_cartesian_1() result(rv)
                    implicit none
                    ! dummy
                    type (cartesian_coordinate_cmplx_type) :: rv

                    rv%y = cmplx(1,2, kind=8)
                    rv%x = cmplx(3,4, kind=8)
                    rv%z = cmplx(5,6, kind=8)
                end function

                function get_test_values_cartesian_2() result(rv)
                    implicit none
                    ! dummy
                    type (cartesian_coordinate_cmplx_type) :: rv

                    rv%y = cmplx(6,5, kind=8)
                    rv%x = cmplx(4,3, kind=8)
                    rv%z = cmplx(2,1, kind=8)
                end function

                function get_test_values_spherical_1() result(rv)
                    implicit none
                    ! dummy
                    type (spherical_coordinate_cmplx_type) :: rv

                    rv%phi = cmplx(1,2, kind=8)
                    rv%rho = cmplx(3,4, kind=8)
                    rv%theta = cmplx(5,6, kind=8)
                end function

                function get_test_values_spherical_2() result(rv)
                    implicit none
                    ! dummy
                    type (spherical_coordinate_cmplx_type) :: rv

                    rv%phi = cmplx(6,5, kind=8)
                    rv%rho = cmplx(4,3, kind=8)
                    rv%theta = cmplx(2,1, kind=8)
                end function

                function evaluate(error, counter, str) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: error
                    integer(kind=4) :: counter
                    character(len=*), optional :: str

                    logical :: rv

                    ! auxiliary
                    double precision, parameter :: ground_truth_e = 10.0_8**(-15.0_8)

                    if (error .gt. ground_truth_e) then
                        if (present(str)) then
                            print *, "  ", str, " ", counter , "error: ", error, " : FAILED"
                        else
                            print *, "  ", counter, "error: ", error, " : FAILED"
                        end if
                        rv = .false.
                    else
                        if (present(str)) then
                            print *, "  ", str, " ", counter , "OK"
                        else
                            print *, "  ", counter, "OK"
                        end if
                        rv = .true.
                    end if
                end function evaluate

                function test_lib_math_list_spherical_operator_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_add_array_cmplx

                function test_lib_math_list_spherical_operator_sub_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-5,-3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-1,1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(3,5, kind=8)

                    value = lib_math_list_spherical_operator_sub_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_sub_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_sub_array_cmplx

                function test_lib_math_list_spherical_operator_0d_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_0d_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_0d_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_0d_add_array_cmplx

                function test_lib_math_list_spherical_operator_array_add_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = get_test_values_spherical_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_array_add_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_add_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_add_0d_cmplx

                function test_lib_math_list_spherical_operator_array_sub_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = get_test_values_spherical_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(5,3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,-1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-3,-5, kind=8)

                    value = lib_math_list_spherical_operator_array_sub_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_sub_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_sub_0d_cmplx

                function test_lib_math_list_spherical_operator_array_real_mul_array_c() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_spherical_coordinate_cmplx_type)  :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = 2.0
                    lhs(2) = 5.0

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(30,25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(20,15, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(10,5, kind=8)

                    value = lib_math_list_spherical_operator_array_real_mul_array_c(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_real_mul_array_c:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_real_mul_array_c

                function test_lib_math_list_spherical_operator_array_c_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = cmplx(0, 1, kind=8)
                    lhs(2) = cmplx(2, 0, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-2,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-4,3, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-6,5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(12,10, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(8,6, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(4,2, kind=8)

                    value = lib_math_list_spherical_operator_array_c_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_c_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_c_mul_array_cmplx

                function test_lib_math_list_spherical_operator_real_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = 4

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(4,8, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(12,16, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(20,24, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(24,20, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(16,12, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(8,4, kind=8)

                    value = lib_math_list_spherical_operator_real_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_real_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_real_mul_array_cmplx

                function test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = cmplx(0, 4, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-8,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-16,12, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-24,20, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-20,24, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-12,16, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-4,8, kind=8)

                    value = lib_math_list_spherical_operator_cmplx_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx

                function test_lib_math_list_spherical_operator_c_array_divide_by_r_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs(1) = 2.0
                    rhs(2) = 4.0

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0.5,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(1.5,2, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(2.5,3, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.5,1.25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,0.75, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5,0.25, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_real_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_c_array_divide_by_r_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_c_array_divide_by_r_array

                function test_lib_math_list_spherical_operator_c_array_divide_by_c_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs(1) = cmplx(0, 2, kind=8)
                    rhs(2) = cmplx(0, 4, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.25, -1.5, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(0.75, -1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.25, -0.5, kind=8)

                    value = lib_math_list_spherical_operator_c_array_divide_by_c_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_c_array_divide_by_c_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_c_array_divide_by_c_array

                function test_lib_math_list_spherical_operator_array_divide_by_real() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = 2

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0.5, 1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(1.5, 2, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(2.5, 3, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(3, 2.5, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(2, 1.5, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(1, 0.5, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_real(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_real:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_real

                function test_lib_math_list_spherical_operator_array_divide_by_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = cmplx(0, 2, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(2.5, -3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1.5, -2, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5, -1, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, "phi  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_cmplx

                function test_lib_math_list_cartesian_operator_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_add_array_cmplx

                function test_lib_math_list_cartesian_operator_sub_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(-5,-3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(-1,1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(3,5, kind=8)

                    value = lib_math_list_cartesian_operator_sub_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_sub_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_sub_array_cmplx

                function test_lib_math_list_cartesian_operator_0d_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_0d_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_0d_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_0d_add_array_cmplx

                function test_lib_math_list_cartesian_operator_array_add_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = get_test_values_cartesian_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_array_add_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_add_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_add_0d_cmplx

                function test_lib_math_list_cartesian_operator_array_sub_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = get_test_values_cartesian_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(5,3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1,-1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(-3,-5, kind=8)

                    value = lib_math_list_cartesian_operator_array_sub_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_sub_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_sub_0d_cmplx

                function test_lib_math_list_cartesian_operator_array_real_mul_array_c() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_cartesian_coordinate_cmplx_type)  :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = 2.0
                    lhs(2) = 5.0

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(30,25, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(20,15, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(10,5, kind=8)

                    value = lib_math_list_cartesian_operator_array_real_mul_array_c(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_real_mul_array_c:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_real_mul_array_c

                function test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = cmplx(0, 1, kind=8)
                    lhs(2) = cmplx(2, 0, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(-2,1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(-4,3, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(-6,5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(12,10, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(8,6, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(4,2, kind=8)

                    value = lib_math_list_cartesian_operator_array_c_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_real_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = 4

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(4,8, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(12,16, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(20,24, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(24,20, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(16,12, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(8,4, kind=8)

                    value = lib_math_list_cartesian_operator_real_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_real_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_real_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = cmplx(0, 4, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(-8,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(-16,12, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(-24,20, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(-20,24, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(-12,16, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(-4,8, kind=8)

                    value = lib_math_list_cartesian_operator_cmplx_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_c_array_divide_by_r_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs(1) = 2.0
                    rhs(2) = 4.0

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0.5,1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(1.5,2, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(2.5,3, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(1.5,1.25, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1,0.75, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.5,0.25, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_real_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_c_array_divide_by_r_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_c_array_divide_by_r_array

                function test_lib_math_list_cartesian_operator_c_array_divide_by_c_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs(1) = cmplx(0, 2, kind=8)
                    rhs(2) = cmplx(0, 4, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(1.25, -1.5, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(0.75, -1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.25, -0.5, kind=8)

                    value = lib_math_list_cartesian_operator_c_array_divide_by_c_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_c_array_divide_by_c_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_c_array_divide_by_c_array

                function test_lib_math_list_cartesian_operator_array_divide_by_real() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = 2

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0.5, 1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(1.5, 2, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(2.5, 3, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(3, 2.5, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(2, 1.5, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(1, 0.5, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_real(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_divide_by_real:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_divide_by_real

                function test_lib_math_list_cartesian_operator_array_divide_by_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = cmplx(0, 2, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(2.5, -3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1.5, -2, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.5, -1, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_divide_by_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, "y  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_divide_by_cmplx

                function test_lib_math_cartesian_point_to_spherical_point() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 1

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type) :: rhs
                    type(spherical_coordinate_real_type), dimension(d) :: lhs
                    type(spherical_coordinate_real_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%x = 1
                    rhs%y = 1
                    rhs%z = 1

                    ground_truth_lhs(1)%rho = sqrt(3.0_8)
                    ground_truth_lhs(1)%theta = 0.9553166181245092_8
                    ground_truth_lhs(1)%phi = 0.25 * PI

                    i=1

                    call lib_math_cartesian_point_to_spherical_point(lhs(i), rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_point_to_spherical_point:"

                    buffer = abs(lhs(i)%rho - ground_truth_lhs(i)%rho)
                    if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                    buffer = abs(lhs(i)%theta - ground_truth_lhs(i)%theta)
                    if (.not. evaluate(buffer, i, "theta")) rv = .false.

                    buffer = abs(lhs(i)%phi - ground_truth_lhs(i)%phi)
                    if (.not. evaluate(buffer, i, "phi  ")) rv = .false.

                end function test_lib_math_cartesian_point_to_spherical_point

                function test_lib_math_spherical_point_to_cartesian_point() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 1

                    ! auxiliary
                    integer :: i
                    type(spherical_coordinate_real_type) :: rhs
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%rho = 1
                    rhs%theta = 0.25 * PI
                    rhs%phi = 0.25 * PI

                    ground_truth_lhs(1)%x = 0.5_8
                    ground_truth_lhs(1)%y = 0.5_8
                    ground_truth_lhs(1)%z = 1.0_8 / sqrt(2.0_8)

                    i=1

                    call lib_math_spherical_point_to_cartesian_point(lhs(i), rhs)

                    rv = .true.
                    print *, "test_lib_math_spherical_point_to_cartesian_point:"

                    buffer = abs(lhs(i)%x - ground_truth_lhs(i)%x)
                    if (.not. evaluate(buffer, i, "x")) rv = .false.

                    buffer = abs(lhs(i)%y - ground_truth_lhs(i)%y)
                    if (.not. evaluate(buffer, i, "y")) rv = .false.

                    buffer = abs(lhs(i)%z - ground_truth_lhs(i)%z)
                    if (.not. evaluate(buffer, i, "z")) rv = .false.

                end function test_lib_math_spherical_point_to_cartesian_point

                function test_lib_math_spherical_point_to_cartesian_point_2() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 1

                    ! auxiliary
                    integer :: i
                    type(spherical_coordinate_real_type) :: rhs
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%rho = 1
                    rhs%theta = -0.25 * PI - 4 * PI
                    rhs%phi = 1.25 * PI - 4 * PI

                    ground_truth_lhs(1)%x = 0.5_8
                    ground_truth_lhs(1)%y = 0.5_8
                    ground_truth_lhs(1)%z = 1.0_8 / sqrt(2.0_8)

                    i=1

                    call lib_math_spherical_point_to_cartesian_point(lhs(i), rhs)

                    rv = .true.
                    print *, "test_lib_math_spherical_point_to_cartesian_point:"

                    buffer = abs(lhs(i)%x - ground_truth_lhs(i)%x)
                    if (.not. evaluate(buffer, i, "x")) rv = .false.

                    buffer = abs(lhs(i)%y - ground_truth_lhs(i)%y)
                    if (.not. evaluate(buffer, i, "y")) rv = .false.

                    buffer = abs(lhs(i)%z - ground_truth_lhs(i)%z)
                    if (.not. evaluate(buffer, i, "z")) rv = .false.

                end function test_lib_math_spherical_point_to_cartesian_point_2

                function test_lib_math_spherical_components_to_cartesian_components_c_a() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 6

                    ! auxiliary
                    integer :: i
                    type(spherical_coordinate_cmplx_type) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: theta
                    real(kind=lib_math_type_kind), dimension(d) :: phi

                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%rho = cmplx(1,0)
                    rhs%theta = cmplx(0,0)
                    rhs%phi = cmplx(0,0)

                    theta(1) = 0.25 * PI
                    phi(1) = 0

                    theta(2) = 0.25 * PI
                    phi(2) = 0.25 * PI

                    theta(3) = 0.25 * PI
                    phi(3) = 0.75 * PI

                    theta(4) = 0.25 * PI
                    phi(4) = 1.25 * PI

                    theta(5) = 0.25 * PI
                    phi(5) = 1.5 * PI

                    theta(6) = 0
                    phi(6) = 0

                    ground_truth_lhs(1)%x = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%y = cmplx(0,0)
                    ground_truth_lhs(1)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(2)%x = 0.5_8
                    ground_truth_lhs(2)%y = 0.5_8
                    ground_truth_lhs(2)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(3)%x = -0.5_8
                    ground_truth_lhs(3)%y = 0.5_8
                    ground_truth_lhs(3)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(4)%x = - 0.5_8
                    ground_truth_lhs(4)%y = - 0.5_8
                    ground_truth_lhs(4)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(5)%x = cmplx(0,0)
                    ground_truth_lhs(5)%y = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(6)%x = cmplx(0,0)
                    ground_truth_lhs(6)%y = cmplx(0,0)
                    ground_truth_lhs(6)%z = 1.0_8

                    do i=1, d
                        lhs(i) = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                               theta(i), &
                                                                                               phi(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_spherical_components_to_cartesian_components_c_a:"
                    do i=1, d
                        buffer = abs(lhs(i)%x - ground_truth_lhs(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(lhs(i)%y - ground_truth_lhs(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(lhs(i)%z - ground_truth_lhs(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_spherical_components_to_cartesian_components_c_a

                function test_lib_math_cartesian_components_to_spherical_components_c_a() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 6

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: theta
                    real(kind=lib_math_type_kind), dimension(d) :: phi

                    type(spherical_coordinate_cmplx_type), dimension(d) :: lhs
                    type(spherical_coordinate_cmplx_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%x = cmplx(0,0)
                    rhs%y = cmplx(0,0)
                    rhs%z = cmplx(1,0)

                    theta(1) = 0.25 * PI
                    phi(1) = 0

                    theta(2) = 0.25 * PI
                    phi(2) = 0.25 * PI

                    theta(3) = 0.25 * PI
                    phi(3) = 0.75 * PI

                    theta(4) = 0.25 * PI
                    phi(4) = 1.25 * PI

                    theta(5) = 0.25 * PI
                    phi(5) = 1.5 * PI

                    theta(6) = 0
                    phi(6) = 0

                    ground_truth_lhs(1)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%phi = cmplx(0,0)

                    ground_truth_lhs(2)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(2)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(2)%phi = cmplx(0,0)

                    ground_truth_lhs(3)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(3)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(3)%phi = cmplx(0,0)

                    ground_truth_lhs(4)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(4)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(4)%phi = cmplx(0,0)

                    ground_truth_lhs(5)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%phi = cmplx(0,0)

                    ground_truth_lhs(6)%rho = 1.0_8
                    ground_truth_lhs(6)%theta = cmplx(0,0)
                    ground_truth_lhs(6)%phi = cmplx(0,0)

                    do i=1, d
                        lhs(i) = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                               theta(i), &
                                                                                               phi(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_cartesian_components_to_spherical_components_c_a:"
                    do i=1, d
                        buffer = abs(lhs(i)%rho - ground_truth_lhs(i)%rho)
                        if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                        buffer = abs(lhs(i)%theta - ground_truth_lhs(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.

                        buffer = abs(lhs(i)%phi - ground_truth_lhs(i)%phi)
                        if (.not. evaluate(buffer, i, "phi  ")) rv = .false.
                    end do

                end function test_lib_math_cartesian_components_to_spherical_components_c_a

                function test_lib_math_cartesian_cross_product_real_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimensioN(d) :: rhs
                    type(cartesian_coordinate_real_type), dimension(d) :: res

                    type(cartesian_coordinate_real_type), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1)%x = 0
                    ground_truth_res(1)%y = 0
                    ground_truth_res(1)%z = 1
                    ground_truth_res(2)%x = 1
                    ground_truth_res(2)%y = 0
                    ground_truth_res(2)%z = 0
                    ground_truth_res(3)%x = 0
                    ground_truth_res(3)%y = 1
                    ground_truth_res(3)%z = 0
                    ground_truth_res(4)%x = -50
                    ground_truth_res(4)%y = 250
                    ground_truth_res(4)%z = -223.5

                    res = lib_math_cartesian_cross_product_real_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_cross_product_real_list:"
                    do i=1, d
                        buffer = abs(res(i)%x - ground_truth_res(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(res(i)%y - ground_truth_res(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(res(i)%z - ground_truth_res(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_cartesian_cross_product_real_list

                function test_lib_math_cartesian_cross_product_cmplx_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimensioN(d) :: rhs
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: res

                    type(cartesian_coordinate_cmplx_type), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = (1, 1)
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = (1, 1)
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = (1, 1)
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = (1, 1)
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = (1, 1)
                    rhs(3)%x = (1, 1)
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1)%x = 0
                    ground_truth_res(1)%y = 0
                    ground_truth_res(1)%z = (0, 2)
                    ground_truth_res(2)%x = (0, 2)
                    ground_truth_res(2)%y = 0
                    ground_truth_res(2)%z = 0
                    ground_truth_res(3)%x = 0
                    ground_truth_res(3)%y = (0, 2)
                    ground_truth_res(3)%z = 0
                    ground_truth_res(4)%x = -50
                    ground_truth_res(4)%y = 250
                    ground_truth_res(4)%z = -223.5

                    res = lib_math_cartesian_cross_product_cmplx_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_cross_product_cmplx_list:"
                    do i=1, d
                        buffer = abs(res(i)%x - ground_truth_res(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(res(i)%y - ground_truth_res(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(res(i)%z - ground_truth_res(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_cartesian_cross_product_cmplx_list

                function test_lib_math_cartesian_dot_product_real_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimensioN(d) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: res

                    real(kind=lib_math_type_kind), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1) = 0
                    ground_truth_res(2) = 0
                    ground_truth_res(3) = 0
                    ground_truth_res(4) = 65.5

                    res = lib_math_cartesian_dot_product_real_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_dot_product_real_list:"
                    do i=1, d
                        buffer = abs(res(i) - ground_truth_res(i))
                        if (.not. evaluate(buffer, i)) rv = .false.
                    end do

                end function test_lib_math_cartesian_dot_product_real_list

                function test_lib_math_cartesian_dot_product_cmplx_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimensioN(d) :: rhs
                    complex(kind=lib_math_type_kind), dimension(d) :: res

                    complex(kind=lib_math_type_kind), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1) = 0
                    ground_truth_res(2) = 0
                    ground_truth_res(3) = 0
                    ground_truth_res(4) = 65.5

                    res = lib_math_cartesian_dot_product_cmplx_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_dot_product_cmplx_list:"
                    do i=1, d
                        buffer = abs(res(i) - ground_truth_res(i))
                        if (.not. evaluate(buffer, i)) rv = .false.
                    end do

                end function test_lib_math_cartesian_dot_product_cmplx_list

                function test_lib_math_list_of_list_list_make_array_cmplx() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    type(list_list_cmplx), dimension(:), allocatable :: list
                    complex(kind=lib_math_type_kind), dimension(:), allocatable :: array
                    complex(kind=lib_math_type_kind), dimension(:), allocatable :: ground_truth_array

                    double precision :: buffer

                    allocate(list(2))

                    call init_list(list(1), 1, 1)
                    call init_list(list(2), 1, 2)

                    list(1)%item(1)%item = (/ 1, 2, 3 /)

                    list(2)%item(1)%item = (/ 1, 2, 3 /)
                    list(2)%item(2)%item = (/ 4, 5, 6, 7, 8 /)

                    i = (1+1)**2-1**2 + (1+2)**2 - 1**2
                    allocate(ground_truth_array(i))
                    ground_truth_array = (/ 1, 2, 3, 1, 2, 3, 4, 5, 6, 7, 8 /)

                    call lib_math_list_of_list_list_make_array_cmplx(list, array)

                    rv = .true.
                    print *, "test_lib_math_list_of_list_list_make_array_cmplx:"
                    do i=1, size(ground_truth_array, 1)
                        buffer = abs(array(i) - ground_truth_array(i))
                        if (.not. evaluate(buffer, i)) rv = .false.
                    end do

                end function test_lib_math_list_of_list_list_make_array_cmplx

                function test_lib_math_array_make_list_of_list_list_cmplx() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    integer :: n
                    integer :: m
                    complex(kind=lib_math_type_kind), dimension(:), allocatable :: array
                    type(list_list_cmplx), dimension(:), allocatable :: list
                    type(list_list_cmplx), dimension(:), allocatable :: ground_truth_list
                    integer, dimension(:), allocatable :: list_fnu
                    integer, dimension(:), allocatable :: list_n

                    double precision :: buffer
                    character(len=30) :: str

                    allocate(ground_truth_list(2))

                    call init_list(ground_truth_list(1), 1, 1)
                    call init_list(ground_truth_list(2), 1, 2)

                    ground_truth_list(1)%item(1)%item = (/ 1, 2, 3 /)

                    ground_truth_list(2)%item(1)%item = (/ 1, 2, 3 /)
                    ground_truth_list(2)%item(2)%item = (/ 4, 5, 6, 7, 8 /)

                    i = (1+1)**2-1**2 + (1+2)**2 - 1**2
                    allocate(array(i))
                    array = (/ 1, 2, 3, 1, 2, 3, 4, 5, 6, 7, 8 /)

!                    list_fnu(:) = 1
!                    list_n(1) = 1
!                    list_n(2) = 2

                    call lib_math_get_list_of_list_list_structure_cmplx(ground_truth_list, list_fnu, list_n)

                    call lib_math_array_make_list_of_list_list_cmplx(array, list_fnu, list_n, list)

                    rv = .true.
                    print *, "test_lib_math_array_make_list_of_list_list_cmplx:"
                    do i=1, size(ground_truth_list, 1)
                        do n = list_fnu(i), list_fnu(i) + list_n(i) - 1
                            do m = -n, n
                                buffer = abs(list(i)%item(n)%item(m) - ground_truth_list(i)%item(n)%item(m))
                                write(str, '(2X, A, 1X, I2, 3X, A, 1X, I2)') "n =", n, "m =", m
                                if (.not. evaluate(buffer, i, trim(str))) rv = .false.
                            end do
                        end do
                    end do

                end function test_lib_math_array_make_list_of_list_list_cmplx

                function test_lib_math_list_cmplx_add() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    type(list_cmplx) :: lhs
                    type(list_cmplx) :: rhs
                    type(list_cmplx) :: res

                    allocate(lhs%item(-1:-1))
                    lhs%item(-1) = cmplx(0, 0, kind = lib_math_type_kind)

                    allocate(rhs%item(-1:1))
                    rhs%item(-1) = cmplx(1, 0, kind = lib_math_type_kind)
                    rhs%item( 0) = cmplx(2, 0, kind = lib_math_type_kind)
                    rhs%item( 1) = cmplx(4, 0, kind = lib_math_type_kind)

                    res = lhs + rhs


                    rv = .true.
                    print *, "test_lib_math_list_cmplx_add:"
                    do i = lbound(rhs%item, 1), ubound(rhs%item, 1)
                        if (res%item(i) .ne. rhs%item(i)) then
                            write(*, '(4X, A, 1X, I2, 2X, A)') "i =", i, "FAILED"
                            rv = .false.
                        else
                            write(*, '(4X, A, 1X, I2, 2X, A)') "i =", i, "OK"
                        end if
                    end do

                end function test_lib_math_list_cmplx_add

                function test_lib_math_list_cmplx_mul() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    type(list_cmplx) :: lhs
                    type(list_cmplx) :: rhs
                    type(list_cmplx) :: res
                    type(list_cmplx) :: ground_truth_res

                    allocate(lhs%item(2:2))
                    lhs%item(2) = cmplx(2, 0, kind=lib_math_type_kind)

                    allocate(rhs%item(3))
                    rhs%item(1) = cmplx(1, 0, kind=lib_math_type_kind)
                    rhs%item(2) = cmplx(2, 0, kind=lib_math_type_kind)
                    rhs%item(3) = cmplx(4, 0, kind=lib_math_type_kind)

                    allocate(ground_truth_res%item(2:2))
                    ground_truth_res%item(2) = cmplx(4, 0, kind=lib_math_type_kind)

                    res = lhs * rhs


                    rv = .true.
                    print *, "test_lib_math_list_cmplx_mul:"
                    do i = lbound(ground_truth_res%item, 1), ubound(ground_truth_res%item, 1)
                        if (res%item(i) .ne. ground_truth_res%item(i)) then
                            write(*, '(4X, A, 1X, I2, 2X, A)') "i =", i, "FAILED"
                            rv = .false.
                        else
                            write(*, '(4X, A, 1X, I2, 2X, A)') "i =", i, "OK"
                        end if
                    end do

                end function test_lib_math_list_cmplx_mul

                function test_lib_math_list_list_cmplx_add() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    integer :: ii
                    type(list_list_cmplx) :: lhs
                    type(list_list_cmplx) :: rhs
                    type(list_list_cmplx) :: res

                    allocate(lhs%item(-1:-1))
                    allocate(lhs%item(-1)%item, source = (/dcmplx(0,0)/))

                    allocate(rhs%item(-1:1))
                    allocate(rhs%item(-1)%item, source = (/dcmplx(1, 0), dcmplx(2,0)/))
                    allocate(rhs%item(0)%item, source = (/dcmplx(2, 0), dcmplx(3,0)/))
                    allocate(rhs%item(1)%item, source = (/dcmplx(4, 0), dcmplx(5,0)/))

                    res = lhs + rhs


                    rv = .true.
                    print *, "test_lib_math_list_list_cmplx_add:"
                    do i = lbound(rhs%item, 1), ubound(rhs%item, 1)
                        do ii = lbound(rhs%item(i)%item, 1), ubound(rhs%item(i)%item, 1)
                            if (res%item(i)%item(ii) .ne. rhs%item(i)%item(ii)) then
                                write(*, '(4X, A, 1X, I2, 2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "ii = ", ii, "FAILED"
                                rv = .false.
                            else
                                write(*, '(4X, A, 1X, I2, 2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "ii = ", ii, "OK"
                            end if
                        end do
                    end do

                end function test_lib_math_list_list_cmplx_add

                function test_lib_math_list_list_cmplx_mul() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    integer :: ii
                    type(list_list_cmplx) :: lhs
                    type(list_list_cmplx) :: rhs
                    type(list_list_cmplx) :: res
                    type(list_list_cmplx) :: ground_truth_res

                    allocate(lhs%item(2:2))
                    allocate(lhs%item(2)%item, source = (/dcmplx(2,0)/))

                    allocate(rhs%item(3))
                    allocate(rhs%item(1)%item, source = (/dcmplx(1, 0), dcmplx(2,0)/))
                    allocate(rhs%item(2)%item, source = (/dcmplx(2, 0), dcmplx(3,0)/))
                    allocate(rhs%item(3)%item, source = (/dcmplx(4, 0), dcmplx(5,0)/))

                    allocate(ground_truth_res%item(2:2))
                    !allocate(ground_truth_res%item(1)%item, source = (/dcmplx(1, 0), dcmplx(2,0)/))
                    allocate(ground_truth_res%item(2)%item, source = (/dcmplx(4, 0)/))
                    !allocate(ground_truth_res%item(3)%item, source = (/dcmplx(4, 0), dcmplx(5,0)/))

                    res = lhs * rhs


                    rv = .true.
                    print *, "test_lib_math_list_list_cmplx_mul:"
                    do i = lbound(ground_truth_res%item, 1), ubound(ground_truth_res%item, 1)
                        do ii = lbound(ground_truth_res%item(i)%item, 1), ubound(ground_truth_res%item(i)%item, 1)
                            if (res%item(i)%item(ii) .ne. ground_truth_res%item(i)%item(ii)) then
                                write(*, '(4X, A, 1X, I2, 2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "ii = ", ii, "FAILED"
                                rv = .false.
                            else
                                write(*, '(4X, A, 1X, I2, 2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "ii = ", ii, "OK"
                            end if
                        end do
                    end do

                end function test_lib_math_list_list_cmplx_mul

                function test_lib_math_get_matrix_rot() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliaray
                    integer :: i
                    double precision :: buffer
                    character(len=30) :: str

                    double precision :: phi
                    double precision :: theta
                    double precision :: psi
                    type(cartresian_coordinate_rot_matrix_type) :: rot

                    type(cartesian_coordinate_real_type), dimension(3) :: point
                    type(cartesian_coordinate_real_type), dimension(3) :: point_rot
                    type(cartesian_coordinate_real_type), dimension(3) :: ground_truth_point_rot


                    phi = PI / 2d0
                    theta = PI / 2d0
                    psi = PI / 2d0
                    rot = lib_math_get_matrix_rot(phi, theta, psi)

                    point(1)%x = 1
                    point(1)%y = 0
                    point(1)%z = 0

                    point_rot(1) = rot * point(1)

                    ground_truth_point_rot(1)%x = 0
                    ground_truth_point_rot(1)%y = 0
                    ground_truth_point_rot(1)%z = 1

                    ! ---

                    point(2)%x = 0
                    point(2)%y = 1
                    point(2)%z = 0

                    point_rot(2) = rot * point(2)

                    ground_truth_point_rot(2)%x = 0
                    ground_truth_point_rot(2)%y = -1
                    ground_truth_point_rot(2)%z = 0

                    ! ---

                    point(3)%x = 0
                    point(3)%y = 0
                    point(3)%z = 1

                    point_rot(3) = rot * point(3)

                    ground_truth_point_rot(3)%x = 1
                    ground_truth_point_rot(3)%y = 0
                    ground_truth_point_rot(3)%z = 0


                    rv = .true.
                    print *, "test_lib_math_get_matrix_rot:"
                    do i=1, size(ground_truth_point_rot, 1)
                        ! x
                        buffer = ground_truth_point_rot(i)%x - point_rot(i)%x
                        write(str, '(2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "x:"
                        if (.not. evaluate(buffer, i, trim(str))) rv = .false.

                        ! y
                        buffer = ground_truth_point_rot(i)%y - point_rot(i)%y
                        write(str, '(2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "y:"
                        if (.not. evaluate(buffer, i, trim(str))) rv = .false.

                        ! z
                        buffer = ground_truth_point_rot(i)%z - point_rot(i)%z
                        write(str, '(2X, A, 1X, I2, 3X, A, 1X)') "i =", i, "z:"
                        if (.not. evaluate(buffer, i, trim(str))) rv = .false.

                        print *, ""
                    end do

                end function test_lib_math_get_matrix_rot

        end function lib_math_type_operator_test_functions

end module lib_math_type_operator
