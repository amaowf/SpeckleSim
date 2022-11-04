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

module lib_ml_fmm_type_operator
    use ml_fmm_type
    implicit none

    private

    ! ----- public functions -----
    public :: lib_ml_fmm_type_operator_constructor

    !public :: operator (*)
    public :: operator (+)
!    public :: operator (-)
!
    public :: operator (.eq.)
!    public :: operator (.ne.)

    public :: assignment (=)

    public :: lib_ml_fmm_type_operator_set_coefficient_zero
!    public :: lib_ml_fmm_type_operator_allocate_coefficient_list
!    public :: lib_ml_fmm_type_operator_deallocate_coefficient_list
!    public :: lib_ml_fmm_type_operator_set_coefficient
!    public :: lib_ml_fmm_type_operator_get_coefficient

    public :: lib_ml_fmm_get_C
    public :: lib_ml_fmm_get_v_y_j

    public :: lib_ml_fmm_get_u_B_i
    public :: lib_ml_fmm_get_u_phi_i_j

    public :: lib_ml_fmm_translation_RR
    public :: lib_ml_fmm_translation_SR
    public :: lib_ml_fmm_translation_SS

    ! ---- public type definitions -----
    public :: ml_fmm_type_operator_procedures
    public :: lib_ml_fmm_procedure_handles
!    public :: ml_fmm_coefficient_add_operator

    ! ----- operator -----
    interface operator (+)
        module procedure lib_ml_fmm_type_operator_coefficient_add
!        module procedure lib_ml_fmm_type_operator_v_add
        module procedure lib_ml_fmm_type_operator_v_add_0D
    end interface

!    interface operator (-)
!!        module procedure lib_ml_fmm_type_operator_v_sub
!!        module procedure lib_ml_fmm_type_operator_v_sub_0D
!    end interface

!    interface operator (*)
!        module procedure lib_ml_fmm_u_dot_coefficient_operator
!    end interface

    interface operator (.eq.)
        module procedure lib_ml_fmm_coefficient_eq
    end interface

!    interface operator (.ne.)
!        module procedure lib_ml_fmm_coefficient_ne
!    end interface

    ! Sum_q=0_p-1[C_q R_q(y_j − x_∗)] = C o R(y_j − x_∗)
    !                                     ^ cor-operator
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
!    interface operator (.DoR.)
!        module procedure lib_ml_fmm_dor_operator
!    end interface

    interface assignment (=)
        module procedure lib_ml_fmm_coefficient_assignment
    end interface

    ! ----- interfaces -----
    interface
        ! Argument
        ! ----
        !   lhs: type (lib_ml_fmm_coefficient)
        !   rhs: type (lib_ml_fmm_coefficient)
        !
        ! Returns
        ! ----
        !   rv: type (lib_ml_fmm_coefficient)
        !       summation of both coefficients
        function ml_fmm_coefficient_add_operator(lhs,rhs) result (rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_coefficient_add_operator

!        function ml_fmm_u_dot_coefficient_operator(lhs,rhs) result (rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), intent(in) :: lhs
!            type (lib_ml_fmm_coefficient), intent(in) :: rhs
!            type (lib_ml_fmm_coefficient) :: rv
!        end function ml_fmm_u_dot_coefficient_operator

        ! Argument
        ! ----
        !   D: type(lib_ml_fmm_coefficient)
        !       D coefficient of the box (n, l_max)
        !   x_c: type(lib_tree_spatial_point)
        !       centre of the box (n, l_max)
        !       HINT: unscale with lib_tree_get_unscaled_point
        !   y_j: type(lib_tree_spatial_point)
        !       scaled point of the j-th element
        !   element_number_j: integer
        !       number of the j-th element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !
        ! Returns
        ! ----
        !   rv: type(lib_ml_fmm_v)
        !       R expansion around y_j
        !
        ! Reference:  Data_Structures_Optimal_Choice_of_Parameters_and_C, eq. 38
        function lib_ml_fmm_dor_operator(D, x_c, data_element_j, element_number_j) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: D
            type(lib_tree_spatial_point), intent(in) :: x_c
            type(lib_tree_data_element), intent(in) :: data_element_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv
        end function lib_ml_fmm_dor_operator

        function ml_fmm_coefficient_eq(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_coefficient), intent(in) :: rhs
            logical :: rv
        end function

!        function ml_fmm_coefficient_ne(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            type(lib_ml_fmm_coefficient), intent(in) :: lhs
!            type(lib_ml_fmm_coefficient), intent(in) :: rhs
!            logical :: rv
!        end function

        subroutine ml_fmm_coefficient_set_zero(coefficient)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
        end subroutine

!        subroutine ml_fmm_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!        end subroutine
!
!        subroutine ml_fmm_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            type(lib_ml_fmm_coefficient_list_list), intent(inout) :: coefficient_list
!        end subroutine

!        subroutine ml_fmm_set_coefficient(coefficient, uindex, hierarchy)
!            use ml_fmm_type
!            use lib_tree_public
!            implicit none
!            ! dummy
!            type(lib_tree_universal_index), intent(in) :: uindex
!            type(lib_ml_fmm_coefficient), intent(in) :: coefficient
!            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
!        end subroutine
!
!        function ml_fmm_get_coefficient(uindex, hierarchy) result(coefficient)
!            use ml_fmm_type
!            use lib_tree_public
!            implicit none
!            ! dummy
!            type(lib_tree_universal_index), intent(in) :: uindex
!            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
!            type(lib_ml_fmm_coefficient) :: coefficient
!        end function

!        function ml_fmm_v_add_operator(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
!            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
!            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv
!        end function

        function ml_fmm_v_add_0D_operator(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_v), intent(in) :: rhs
            type (lib_ml_fmm_v) :: rv
        end function

!        function ml_fmm_v_sub_operator(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
!            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
!            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv
!        end function
!
!        function ml_fmm_v_sub_0D_operator(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), intent(in) :: lhs
!            type (lib_ml_fmm_v), intent(in) :: rhs
!            type (lib_ml_fmm_v) :: rv
!        end function

!        ! Basis function: A
!        function lib_ml_fmm_get_A(x) result(A)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_A) :: A
!        end function

!        ! Basis function: A_i
!        function lib_ml_fmm_get_A_i(x, data_element) result(A_i)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_data_element) :: data_element
!            type(lib_ml_fmm_coefficient) :: A_i
!        end function

!        ! Basis function: B
!        function lib_ml_fmm_get_B(x) result(B)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_B) :: B
!        end function

        ! Argument
        ! ----
        !   x: type(lib_tree_spatial_point)
        !       centre of the box (n, l_max)
        !       HINT: unscale with lib_tree_get_unscaled_point
        !   data_element: type(lib_tree_data_element)
        !       data element
        !   element_number: integer
        !       number of the data element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !
        ! Returns
        ! ----
        !   B_i: type(lib_ml_fmm_coefficient)
        !       S expansion around x of the data element
        !
        ! Reference:  Data_Structures_Optimal_Choice_of_Parameters_and_C, eq. 32
        function lib_ml_fmm_get_u_B_i(x, data_element, element_number) result(B_i)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element), intent(in) :: data_element
            integer(kind=4), intent(in) :: element_number

            type(lib_ml_fmm_coefficient) :: B_i
        end function


        ! Argument
        ! ----
        !   x: type(lib_tree_spatial_point)
        !       centre of the box
        !   data_element: type(lib_tree_data_element), dimension(:)
        !       data element list of the box
        !   element_number: integer(kind=4), dimension(:)
        !       number of the element
        !       HINT: X, Y, and XY lists are concatenated
        !
        ! Returns
        ! ----
        !   C: type(lib_ml_fmm_coefficient)
        !       coefficient of the box
        !
        !
        !
        function lib_ml_fmm_get_C(x, data_element, element_number) result(C)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element), dimension(:), intent(in) :: data_element
            integer(kind=4), dimension(:), intent(in) :: element_number

            type(lib_ml_fmm_coefficient) :: C
        end function

!        ! Basis function: S
!        function lib_ml_fmm_get_S(x) result(S)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_S) :: S
!        end function
!
!        ! Basis function: R
!        function lib_ml_fmm_get_R(x) result(R)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_R) :: R
!        end function

!        ! Local expansion (inner or Regular expansion)
!        !
!        ! Restriction
!        ! ----
!        ! Example calculation:
!        !       phi_i(y) = A_i (x_∗) o R(y − x_∗)     (8)
!        !
!        !   "Here the series is valid in the domain |y − x_∗| .le. r_c |x_i − x_∗| (see Fig. 2),
!        !   where 0 < r c < 1 is some real number."
!        !
!        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   i
!        !   x_*
!        !   y
!        !
!        ! Returns
!        ! ----
!        !   phi_i
!        !
!        !
!        function lib_ml_fmm_expansion_R(i, x, y) result(phi_i)
!            use lib_tree_type
!            implicit none
!            ! dummy
!            integer, intent(in) :: i
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_spatial_point), intent(in) :: y
!            integer :: phi_i                   ! todo: define type
!        end function lib_ml_fmm_expansion_R

!        ! Far field expansion (outer, Singular or multipole expansion)
!        !
!        ! Restriction
!        ! ----
!        ! Example calculation:
!        !   "Any function phi_i(y) has a complementary expansion valid outside
!        !   a d-dimensional sphere centered at y = x_∗ with radius R_c |x_i − x_∗| :
!        !       phi_i(y) = B_i (x_∗) o S(y − x_∗), |y - x_*| .ge. R_c |x_i - x_*|,     (10)
!        !
!        !   where R_c > 1 is a real number similar to r_c".
!        !   "Even though for many physical fields, such as the Green’s function
!        !   for Laplace’s equation, the function S(y − x_∗) is singular at y = x_∗ ,
!        !   this condition is not necessary. In particular we can have S = R."
!        !
!        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   i
!        !   x_*
!        !   y
!        !
!        ! Returns
!        ! ----
!        !   phi_i
!        !
!        !
!        function lib_ml_fmm_expansion_S(i, x, y) result(phi_i_j)
!            use lib_tree_type
!            implicit none
!            ! dummy
!            integer, intent(in) :: i
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_spatial_point), intent(in) :: y
!            integer :: phi_i_j                   ! todo: define type
!        end function lib_ml_fmm_expansion_S

        ! Argument
        ! ----
        !   data_element_i: type(lib_tree_data_element)
        !       i-th data element
        !   element_number_i: integer
        !       number of the i-th data element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !   y_j: type(lib_tree_spatial_point)
        !       scaled point of the j-th data element
        !       HINT: unscale with lib_tree_get_unscaled_point
        !   element_number_j: integer
        !       number of the j-th data element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !
        ! Returns
        ! ----
        !   rv: type(lib_ml_fmm_v)
        !       the result of calculation of u_i * phi_i(y_j)
        !
        ! Reference:  Data_Structures_Optimal_Choice_of_Parameters_and_C, eq. 38
        function lib_ml_fmm_get_u_phi_i_j(data_element_i, element_number_i, data_element_j, element_number_j) result(rv)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_data_element), intent(in) :: data_element_i
            integer(kind=4), intent(in) :: element_number_i
            type(lib_tree_data_element), intent(in) :: data_element_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv

        end function

        ! Argument
        ! ----
        !   data_element: type(lib_tree_data_element), dimension(:)
        !       data element
        !   element_number_i: integer
        !       number of the i-th data element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !   y_j: type(lib_tree_spatial_point), dimension(:)
        !       scaled point of the j-th data element
        !       HINT: unscale with lib_tree_get_unscaled_point
        !   element_number_j: integer
        !       number of the j-th data element at the concatenated element data list.
        !       CONVENTION:
        !           Internal representation of the element data list
        !           from the 1-st to the N-th element.
        !           ------------------------------------
        !           |1    X    |     Y    |     XY    N|
        !           ------------------------------------
        !           X-, Y-, XY- hierarchy
        !
        ! Returns
        ! ----
        !   rv: type(lib_ml_fmm_v)
        !       the result of calculation of u_i * phi_i(y_j)
        !
        ! Reference:  Data_Structures_Optimal_Choice_of_Parameters_and_C, eq. 38

        function lib_ml_fmm_get_v_y_j(data_element_y_j, element_number_j, data_element_e2, element_number_e2, D) result(rv)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_data_element), intent(in) :: data_element_y_j
            integer(kind=CORRESPONDENCE_VECTOR_KIND), intent(in) :: element_number_j
            type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_element_e2
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(in) :: element_number_e2
            type(lib_ml_fmm_coefficient), intent(in) :: D

            type(lib_ml_fmm_v) :: rv

        end function

        ! Translation: local-to-local (Regular-to-Regular)
        !
        ! Arguments
        ! ----
        !   A_i_1
        !   x_1
        !   x_2
        !
        ! Returns
        ! ----
        !   A_i_2
        !
        function lib_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
        end function

        ! Translation: far-to-local (Singular-to-Regular)
        !
        ! Arguments
        ! ----
        !   B_i_1
        !       set of expansion coefficients
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_1: spatial point
        !       origin of coordinate system 2
        !
        ! Returns
        ! ----
        !   A_i_2
        !       set of expansion coefficients
        !
        function lib_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
        end function

        ! Translation: far-to-far (Singular-to-Singular)
        !
        ! Arguments
        ! ----
        !   B_i_1
        !       set of expansion coefficients
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_2: spatial point
        !       origin of coordinate system 2
        !
        ! Returns
        ! ----
        !   B_i_2
        !       set of expansion coefficients
        !
        function lib_ml_fmm_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: B_i_2
        end function
    end interface

    type lib_ml_fmm_procedure_handles
!        procedure(lib_ml_fmm_get_A_i), pointer, nopass :: get_A_i => null()
        procedure(lib_ml_fmm_get_u_B_i), pointer, nopass :: get_u_B_i => null()
!        procedure(lib_ml_fmm_get_S), pointer, nopass :: get_S => null()
!        procedure(lib_ml_fmm_get_R), pointer, nopass :: get_R => null()
!        procedure(lib_ml_fmm_expansion_R), pointer, nopass :: expansion_R => null()
!        procedure(lib_ml_fmm_expansion_S), pointer, nopass :: expansion_S => null()
        procedure(lib_ml_fmm_dor_operator), pointer, nopass :: dor => null()
        procedure(lib_ml_fmm_get_u_phi_i_j), pointer, nopass :: get_u_phi_i_j => null()
        procedure(lib_ml_fmm_translation_RR), pointer, nopass :: get_translation_RR => null()
        procedure(lib_ml_fmm_translation_SR), pointer, nopass :: get_translation_SR => null()
        procedure(lib_ml_fmm_translation_SS), pointer, nopass :: get_translation_SS => null()

        procedure(lib_ml_fmm_get_C), pointer, nopass :: get_c => null()
        procedure(lib_ml_fmm_get_v_y_j), pointer, nopass :: get_v_y_j => null()
    end type lib_ml_fmm_procedure_handles

    type ml_fmm_type_operator_procedures
        procedure(ml_fmm_coefficient_add_operator), pointer, nopass :: coefficient_add => null()
!        procedure(ml_fmm_u_dot_coefficient_operator), pointer, nopass :: u_dot_coefficient => null()
        procedure(ml_fmm_coefficient_set_zero), pointer, nopass :: coefficient_set_zero => null()
!        procedure(ml_fmm_allocate_coefficient_list), pointer, nopass :: allocate_coefficient_list => null()
!        procedure(ml_fmm_set_coefficient), pointer, nopass :: set_coefficient => null()
!        procedure(ml_fmm_get_coefficient), pointer, nopass :: get_coefficient => null()

!        procedure(ml_fmm_v_add_operator), pointer, nopass :: v_add => null()
        procedure(ml_fmm_v_add_0D_operator), pointer, nopass :: v_add_0D => null()
!        procedure(ml_fmm_v_sub_operator), pointer, nopass :: v_sub => null()
!        procedure(ml_fmm_v_sub_0D_operator), pointer, nopass :: v_sub_0D => null()

!        procedure(ml_fmm_deallocate_coefficient_list), pointer, nopass :: deallocate_coefficient_list => null()
        procedure(ml_fmm_coefficient_eq), pointer, nopass :: coefficient_eq => null()
!        procedure(ml_fmm_coefficient_eq), pointer, nopass :: coefficient_ne => null()
    end type

    ! ----- member procedures -----
#ifdef __GFORTRAN__
    procedure(ml_fmm_coefficient_add_operator), pointer :: m_coefficient_add => null()
!    procedure(ml_fmm_u_dot_coefficient_operator), pointer :: m_u_dot_coefficient => null()
!    procedure(lib_ml_fmm_dor_operator), pointer :: m_dor => null()
    procedure(ml_fmm_coefficient_set_zero), pointer :: m_coefficient_set_zero => null()
!    procedure(ml_fmm_allocate_coefficient_list), pointer :: m_allocate_coefficient_list => null()
!    procedure(ml_fmm_set_coefficient), pointer :: m_set_coefficient => null()
!    procedure(ml_fmm_get_coefficient), pointer :: m_get_coefficient => null()

!    procedure(ml_fmm_v_add_operator), pointer :: m_v_add => null()
    procedure(ml_fmm_v_add_0D_operator), pointer :: m_v_add_0D => null()
!    procedure(ml_fmm_v_sub_operator), pointer :: m_v_sub => null()
!    procedure(ml_fmm_v_sub_0D_operator), pointer :: m_v_sub_0D => null()
	
	!    procedure(ml_fmm_deallocate_coefficient_list), pointer :: m_deallocate_coefficient_list => null()
    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_eq => null()
!    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_ne => null()
#else
    procedure(ml_fmm_coefficient_add_operator), pointer :: m_coefficient_add
!    procedure(ml_fmm_u_dot_coefficient_operator), pointer :: m_u_dot_coefficient
!    procedure(ml_fmm_dor_operator), pointer :: m_dor
    procedure(ml_fmm_coefficient_set_zero), pointer :: m_coefficient_set_zero
!    procedure(ml_fmm_allocate_coefficient_list), pointer :: m_allocate_coefficient_list
!    procedure(ml_fmm_set_coefficient), pointer :: m_set_coefficient
!    procedure(ml_fmm_get_coefficient), pointer :: m_get_coefficient
	
!	procedure(ml_fmm_v_add_operator), pointer :: m_v_add
    procedure(ml_fmm_v_add_0D_operator), pointer :: m_v_add_0D
!    procedure(ml_fmm_v_sub_operator), pointer :: m_v_sub
!    procedure(ml_fmm_v_sub_0D_operator), pointer :: m_v_sub_0D
	
	!    procedure(ml_fmm_deallocate_coefficient_list), pointer :: m_deallocate_coefficient_list => null()
    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_eq
!    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_ne
#endif

    contains

    subroutine lib_ml_fmm_type_operator_constructor(operator_procedures)
!    subroutine lib_ml_fmm_type_operator_constructor(coefficient_add, u_dot_coefficient, cor, &
!                                                    set_coefficient_zero)!, &
!                                                    allocate_coefficient_list, &
!                                                    set_coefficient, get_coefficient, &
!                                                    deallocate_coefficient_list)
        use ml_fmm_type
        implicit none
        ! dummy
!        procedure(ml_fmm_coefficient_add_operator) :: coefficient_add
!        procedure(ml_fmm_u_dot_coefficient_operator) :: u_dot_coefficient
!        procedure(ml_fmm_cor_operator) :: cor
!        procedure(ml_fmm_coefficient_set_zero) :: set_coefficient_zero
!        procedure(ml_fmm_allocate_coefficient_list) :: allocate_coefficient_list
!        procedure(ml_fmm_set_coefficient) :: set_coefficient
!        procedure(ml_fmm_get_coefficient) :: get_coefficient
!        procedure(ml_fmm_deallocate_coefficient_list) :: deallocate_coefficient_list
!
!        m_coefficient_add => coefficient_add
!        m_u_dot_coefficient => u_dot_coefficient
!        m_cor => cor
!        m_coefficient_set_zero => set_coefficient_zero
!        m_allocate_coefficient_list => allocate_coefficient_list
!        m_set_coefficient => set_coefficient
!        m_get_coefficient => get_coefficient
!        m_deallocate_coefficient_list => deallocate_coefficient_list


        type(ml_fmm_type_operator_procedures) :: operator_procedures
        m_coefficient_add => operator_procedures%coefficient_add
!        m_u_dot_coefficient => operator_procedures%u_dot_coefficient
!        m_dor => operator_procedures%dor
        m_coefficient_set_zero => operator_procedures%coefficient_set_zero
!        m_allocate_coefficient_list => operator_procedures%allocate_coefficient_list
!        m_set_coefficient => operator_procedures%set_coefficient
!        m_get_coefficient => operator_procedures%get_coefficient

        !m_v_add => operator_procedures%v_add
        m_v_add_0D => operator_procedures%v_add_0D
        !m_v_sub => operator_procedures%v_sub
        !m_v_sub_0D => operator_procedures%v_sub_0D
!        m_deallocate_coefficient_list => operator_procedures%deallocate_coefficient_list

        m_coefficient_eq => operator_procedures%coefficient_eq
!        m_coefficient_ne => operator_procedures%coefficient_ne

    end subroutine lib_ml_fmm_type_operator_constructor

    subroutine lib_ml_fmm_coefficient_assignment(lhs, rhs)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: rhs

        type(lib_ml_fmm_coefficient), intent(inout) :: lhs

        lhs%a_nm = rhs%a_nm
        lhs%b_nm = rhs%b_nm

        if (allocated(lhs%r)) then
            deallocate(lhs%r)
        end if

        if (allocated(rhs%r)) then
            lhs%r = rhs%r
        end if

        if (allocated(lhs%c)) then
            deallocate(lhs%c)
        end if

        if (allocated(rhs%c)) then
            lhs%c = rhs%c
        end if

        lhs%uindex = rhs%uindex

    end subroutine

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (associated (m_coefficient_set_zero)) then
            call m_coefficient_set_zero(coefficient)
        else
            print *, "lib_ml_fmm_type_operator_set_coefficient_zero:  ERROR"
            print *, "  m_coefficient_set_zero is not associated"
        end if
    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero

!    subroutine lib_ml_fmm_type_operator_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_allocate_coefficient_list) ) then
!                call m_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            else
!                print *, "lib_ml_fmm_type_operator_allocate_coefficient_list:  ERROR"
!                print *, "  m_allocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_allocate_coefficient_list
!
!    subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_deallocate_coefficient_list) ) then
!                call m_deallocate_coefficient_list(coefficient_list)
!            else
!                print *, "lib_ml_fmm_type_operator_deallocate_coefficient_list:  ERROR"
!                print *, "  m_deallocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list

!        subroutine lib_ml_fmm_type_operator_set_coefficient(coefficient, uindex, hierarchy)
!            use ml_fmm_type
!            use lib_tree_type
!            implicit none
!            ! dummy
!            type(lib_tree_universal_index), intent(inout) :: uindex
!            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
!            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
!
!            if ( associated(m_set_coefficient)) then
!                call m_set_coefficient(coefficient, uindex, hierarchy)
!            else
!                print *, "lib_ml_fmm_type_operator_set_coefficient:  ERROR"
!                print *, "  m_set_coefficient is not associated"
!            end if
!        end subroutine
!
!        function lib_ml_fmm_type_operator_get_coefficient(uindex, hierarchy) result(coefficient)
!            use ml_fmm_type
!            use lib_tree_type
!            implicit none
!            ! dummy
!            type(lib_tree_universal_index), intent(inout) :: uindex
!            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
!            type(lib_ml_fmm_coefficient) :: coefficient
!
!            if ( associated(m_set_coefficient)) then
!                coefficient = m_get_coefficient(uindex, hierarchy)
!            else
!                print *, "lib_ml_fmm_type_operator_get_coefficient:  ERROR"
!                print *, "  m_get_coefficient is not associated"
!            end if
!        end function
        
        function lib_ml_fmm_type_operator_coefficient_add(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv
            
            if ( associated(m_coefficient_add)) then
                rv = m_coefficient_add(lhs, rhs)
            else
                print *, "lib_ml_fmm_type_operator_coefficient_add:  ERROR"
                print *, "  m_coefficient_add is not associated"
            end if
            
        end function

!        function lib_ml_fmm_type_operator_v_add(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
!            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
!            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv
!
!            if ( associated(m_v_add)) then
!                rv = m_v_add(lhs, rhs)
!            else
!                print *, "lib_ml_fmm_type_operator_v_add:  ERROR"
!                print *, "  m_v_add is not associated"
!            end if
!        end function

        function lib_ml_fmm_type_operator_v_add_0d(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_v), intent(in) :: rhs
            type (lib_ml_fmm_v) :: rv

            if ( associated(m_v_add_0d)) then
                rv = m_v_add_0d(lhs, rhs)
            else
                print *, "lib_ml_fmm_type_operator_v_add_0d:  ERROR"
                print *, "  m_v_add_0d is not associated"
            end if
        end function

!        function lib_ml_fmm_type_operator_v_sub(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
!            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
!            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv
!
!            if ( associated(m_v_sub)) then
!                rv = m_v_sub(lhs, rhs)
!            else
!                print *, "lib_ml_fmm_type_operator_v_sub:  ERROR"
!                print *, "  m_v_sub is not associated"
!            end if
!        end function

!        function lib_ml_fmm_type_operator_v_sub_0d(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), intent(in) :: lhs
!            type (lib_ml_fmm_v), intent(in) :: rhs
!            type (lib_ml_fmm_v) :: rv
!
!            if ( associated(m_v_sub_0d)) then
!                rv = m_v_sub_0D(lhs, rhs)
!            else
!                print *, "lib_ml_fmm_type_operator_v_sub_0d:  ERROR"
!                print *, "  m_v_sub_0d is not associated"
!            end if
!        end function
        
!        function lib_ml_fmm_u_dot_coefficient_operator(lhs,rhs) result (rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type (lib_ml_fmm_v), intent(in) :: lhs
!            type (lib_ml_fmm_coefficient), intent(in) :: rhs
!            type (lib_ml_fmm_coefficient) :: rv
!
!            if (associated(m_u_dot_coefficient)) then
!                rv = m_u_dot_coefficient(lhs, rhs)
!            else
!                print *, "lib_ml_fmm_u_dot_coefficient_operator:  ERROR"
!                print *, "  m_u_dot_coefficient is not associated"
!            end if
!
!        end function lib_ml_fmm_u_dot_coefficient_operator
        
        function lib_ml_fmm_coefficient_eq(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_coefficient), intent(in) :: rhs
            logical :: rv

            rv = m_coefficient_eq(lhs, rhs)

        end function lib_ml_fmm_coefficient_eq
        
!        function lib_ml_fmm_coefficient_ne(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            type(lib_ml_fmm_coefficient), intent(in) :: lhs
!            type(lib_ml_fmm_coefficient), intent(in) :: rhs
!            logical :: rv
!
!            rv = m_coefficient_ne(lhs, rhs)
!
!        end function lib_ml_fmm_coefficient_ne
        
!        function lib_ml_fmm_dor_operator(lhs, rhs) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient), intent(in) :: lhs
!            type(lib_tree_spatial_point), intent(in) :: rhs
!            type(lib_ml_fmm_v) :: rv
!
!            rv = m_dor(lhs, rhs)
!
!        end function lib_ml_fmm_dor_operator

!        function lib_ml_fmm_dor_operator(D, x_c, y_j, element_number_j) result(rv)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient), intent(in) :: D
!            type(lib_tree_spatial_point), intent(in) :: x_c
!            type(lib_tree_spatial_point), intent(in) :: y_j
!            integer(kind=4), intent(in) :: element_number_j
!            type(lib_ml_fmm_v) :: rv
!
!            rv = m_dor(D, x_c, y_j, element_number_j)
!
!        end function lib_ml_fmm_dor_operator

end module lib_ml_fmm_type_operator
