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

module ml_fmm_math
    use ml_fmm_type
    use lib_ml_fmm_type_operator
    implicit none

    private

    ! --- public functions ---
    public :: ml_fmm_type_operator_get_procedures

    public :: ml_fmm_get_test_procedures
    public :: lib_ml_fmm_type_operator_test_functions

    contains

    ! Argument
    ! ----
    !   m_procedure_type: integer, optional (std=0)
    !       operators manipulating the vector v, u and the coefficients C, D
    !       0: r
    !       1: a_nm and b_nm: type(list_list_cmplx)
    function ml_fmm_type_operator_get_procedures(procedure_type) result (operator_procedures)
        implicit none
        ! dummy
        integer, intent(in), optional :: procedure_type
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! auxiliary
        integer :: m_procedure_type

        m_procedure_type = 0
        if (present(procedure_type)) m_procedure_type = procedure_type

        if (m_procedure_type .eq. 0) then
            ! load test procedure functions
            operator_procedures%coefficient_add => ml_fmm_coefficient_add_operator_list_2_cmplx
            operator_procedures%coefficient_set_zero => test_set_coefficient_zero
!            operator_procedures%u_dot_coefficient => test_u_dot_coefficient
            operator_procedures%coefficient_eq => test_coefficient_eq
!            operator_procedures%coefficient_ne => test_coefficient_ne
!            operator_procedures%v_add => test_v_operator_add
            operator_procedures%v_add_0D => test_v_operator_add_0D
!            operator_procedures%v_sub => test_v_operator_sub
!            operator_procedures%v_sub_0D => test_v_operator_sub_0D
        else if (m_procedure_type .eq. 1) then
            ! load FoggySim specialised procedure functions
            operator_procedures%coefficient_add => ml_fmm_coefficient_add_operator_list_2_cmplx
            operator_procedures%v_add_0D => ml_fmm_type_operator_v_add_0d_list_2_cmplx
            operator_procedures%coefficient_set_zero => lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx
        else if (m_procedure_type .eq. 2) then
            ! load SpeckleSim specialised procedure functions
            operator_procedures%coefficient_add => ml_fmm_coefficient_add_operator_list_2_cmplx
            operator_procedures%v_add_0D => ml_fmm_type_operator_v_add_0d_list_2_cmplx
            operator_procedures%coefficient_set_zero => lib_ml_fmm_type_operator_set_coefficient_zero_list_2_c_alt
        end if



    end function ml_fmm_type_operator_get_procedures

    function ml_fmm_get_test_procedures() result (handle)
        implicit none
        ! dummy
        type(lib_ml_fmm_procedure_handles) :: handle

        handle%get_u_B_i => test_get_u_B_i
        handle%get_u_phi_i_j => test_get_u_phi_i_j
        handle%get_translation_RR  => test_translation_RR
        handle%get_translation_SR  => test_translation_SR
        handle%get_translation_SS  => test_translation_SS
        handle%dor => test_dor
    end function

    function ml_fmm_coefficient_add_operator_list_2_cmplx(lhs,rhs) result (rv)
        use libmath
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
        type (lib_ml_fmm_coefficient) :: rv

        rv%a_nm = lhs%a_nm + rhs%a_nm
        rv%b_nm = lhs%b_nm + rhs%b_nm

        if (allocated(lhs%c) .or. allocated(rhs%c)) rv%c = lib_math_cmplx_array_add(lhs%c, rhs%c)
        if (allocated(lhs%r) .or. allocated(rhs%r)) rv%r = lib_math_real_array_add(lhs%r, rhs%r)

    end function ml_fmm_coefficient_add_operator_list_2_cmplx

    function ml_fmm_type_operator_v_add_0d_list_2_cmplx(lhs, rhs) result(rv)
        use libmath
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        rv%a_nm = lhs%a_nm + rhs%a_nm
        rv%b_nm = lhs%b_nm + rhs%b_nm

        if (allocated(lhs%c) .or. allocated(rhs%c)) rv%c = lib_math_cmplx_array_add(lhs%c, rhs%c)
        if (allocated(lhs%r) .or. allocated(rhs%r)) rv%r = lib_math_real_array_add(lhs%r, rhs%r)

    end function ml_fmm_type_operator_v_add_0d_list_2_cmplx

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        call init_list(coefficient%a_nm, 1 ,1, dcmplx(0,0))
        call init_list(coefficient%b_nm, 1 ,1, dcmplx(0,0))

    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_c_alt(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        call deallocate_list(coefficient%a_nm)
        call deallocate_list(coefficient%b_nm)

        allocate(coefficient%a_nm%item(1))
        allocate(coefficient%a_nm%item(1)%item(1))
        coefficient%a_nm%item(1)%item(1) = 0d0

        allocate(coefficient%b_nm%item(1))
        allocate(coefficient%b_nm%item(1)%item(1))
        coefficient%b_nm%item(1)%item(1) = 0d0

        if (allocated(coefficient%c)) then
            deallocate(coefficient%c)
        end if
        allocate(coefficient%c, source = (/cmplx(0, 0, kind=lib_math_type_kind)/))

        if (allocated(coefficient%r)) then
            deallocate(coefficient%r)
        end if
        allocate(coefficient%r, source = (/real(0, kind=lib_math_type_kind)/))

    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_c_alt

!    function lib_ml_fmm_type_operator_coefficient_eq_list_2_cmplx(lhs, rhs) result (rv)
!        implicit none
!        ! dummy
!        type(lib_ml_fmm_coefficient), intent(in) :: lhs
!        type(lib_ml_fmm_coefficient), intent(in) :: rhs
!        logical :: rv
!
!        ! auxiliary
!        integer :: n
!        integer :: m
!
!        rv = .true.
!        if ( (lbound(lhs%a_nm%item, 1) .eq. lbound(rhs%a_nm%item)) &
!             .and. (ubound(lhs%a_nm%item, 1) .eq. ubound(rhs%a_nm%item)) ) then
!
!            do n = lbound(lhs%a_nm%item, 1), ubound(lhs%a_nm%item, 1)
!                do m = -n, n
!                    if (lhs%a_nm%item(n)%item(m) .ne. rhs%a_nm%item(n)%item(m)) then
!                        rv = .false.
!                        return
!                    end if
!
!                    if (lhs%b_nm%item(n)%item(m) .ne. rhs%b_nm%item(n)%item(m)) then
!                        rv = .false.
!                        return
!                    end if
!                end do
!            end do
!
!        else
!            rv = .false.
!        end if
!
!    end function lib_ml_fmm_type_operator_coefficient_eq_list_2_cmplx

    ! ----- test functions ----
    function lib_ml_fmm_type_operator_test_functions() result (error_counter)
        implicit none
        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_ml_fmm_type_operator_constructor()) then
            error_counter = error_counter + 1
        end if

        print *, "-------------lib_ml_fmm_type_operator_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_ml_fmm_type_operator_test_functions tests: OK"
        else
            print *, error_counter,"lib_ml_fmm_type_operator_test_functions test(s) FAILED"
        end if
        print *, "--------------------------------------------------------------------"

        contains

        function test_lib_ml_fmm_type_operator_constructor() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=1), parameter :: length = 3
            type(lib_ml_fmm_v) :: lhs_u
            type(lib_ml_fmm_coefficient) :: lhs_coeff
            type(lib_ml_fmm_coefficient) :: rhs_coeff
            type(lib_ml_fmm_coefficient) :: res_coeff
            type(lib_tree_spatial_point) :: rhs_R
!            type(lib_ml_fmm_v) :: res_v

            type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures

            ml_fmm_operator_procedures%coefficient_add => ml_fmm_coefficient_add_operator_list_2_cmplx
            ml_fmm_operator_procedures%coefficient_set_zero => test_set_coefficient_zero
!            ml_fmm_operator_procedures%dor => test_dor
!            ml_fmm_operator_procedures%u_dot_coefficient => test_u_dot_coefficient
!            ml_fmm_operator_procedures%v_add => test_v_operator_add
            ml_fmm_operator_procedures%v_add_0D => test_v_operator_add_0D
!            ml_fmm_operator_procedures%v_sub => test_v_operator_sub
!            ml_fmm_operator_procedures%v_sub_0D => test_v_operator_sub_0D

            call lib_ml_fmm_type_operator_constructor(ml_fmm_operator_procedures)
!            call lib_ml_fmm_type_operator_constructor(test_c_add, test_u_dot_coefficient, test_cor, &
!                                                      test_set_coefficient_zero)!, &
!!                                                      test_allocate_coefficient_list, &
!!                                                      test_set_coefficient, test_get_coefficient, &
!!                                                      test_deallocate_coefficient_list)

            allocate (lhs_u%r(length))
            allocate (lhs_coeff%r(length))
            allocate (rhs_coeff%r(length))

            lhs_coeff%r = (/5, 4, 3/)
            call lib_ml_fmm_type_operator_set_coefficient_zero(lhs_coeff)
            lhs_coeff%r = (/5, 4, 3/)
            rhs_coeff%r = (/2, 1, 3/)

            lhs_u%r = (/1, 2, 3/)

            res_coeff = lhs_coeff + rhs_coeff

!            res_coeff = lhs_u * res_coeff

            rhs_R%x(:2) = (/1.5, 7.0/)
!            res_v = lhs_coeff .dor. rhs_R

            ! todo: add condition
            rv = .true.

        end function test_lib_ml_fmm_type_operator_constructor

    end function lib_ml_fmm_type_operator_test_functions

    function test_c_add(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
        type (lib_ml_fmm_coefficient) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        if (size(lhs%r) .eq. size(rhs%r)) then
            length = size(lhs%r)
            allocate (rv%r(length))
            do i=1, length
                rv%r(i) = lhs%r(i) + rhs%r(i)
            end do
        end if

        if (allocated(lhs%c) .and. allocated(rhs%c)) then
            if (size(lhs%c) .eq. size(rhs%c)) then
                length = size(lhs%c)
                allocate (rv%c(length))
                do i=1, length
                    rv%c(i) = lhs%c(i) + rhs%c(i)
                end do
            end if
        end if

        rv%a_nm = lhs%a_nm + rhs%a_nm
        rv%b_nm = lhs%b_nm + rhs%b_nm

    end function

    function test_u_dot_coefficient(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_coefficient), intent(in) :: rhs
        type (lib_ml_fmm_coefficient) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        if (size(lhs%r) .eq. size(rhs%r)) then
            length = size(lhs%r)
            allocate (rv%r(length))
            do i=1, length
                rv%r(i) = lhs%r(i) * rhs%r(i)
            end do
        end if
    end function test_u_dot_coefficient

    function test_dor(D, x_c, data_element_j, element_number_j)  result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: D
        type(lib_tree_spatial_point), intent(in) :: x_c
        type(lib_tree_data_element), intent(in) :: data_element_j
        integer(kind=4), intent(in) :: element_number_j
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        allocate (rv%r(1))
        rv%r(1) = 0
        length = size(D%r)
        do i=1, length
            rv%r(1) = rv%r(1) + D%r(i) * 2!(y_j%x(1) - x_c%x(1))
        end do

        if (allocated(D%c)) then
            allocate (rv%c(1))
            rv%c(1) = 0
            length = size(D%c)
            do i=1, length
                rv%c(1) = rv%c(1) + D%c(i) * 2!(y_j%x(1) - x_c%x(1))
            end do
        end if

    end function test_dor

    function test_coefficient_eq(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: lhs
        type(lib_ml_fmm_coefficient), intent(in) :: rhs
        logical :: rv

        rv = .true.

        if (lhs%r(1) .ne. rhs%r(1)) then
            rv = .false.
        end if

        if (allocated(lhs%c) .or. allocated(rhs%c)) then
            if (lhs%c(1) .ne. rhs%c(1)) then
                rv = .false.
            end if
        end if

        if (allocated(lhs%a_nm%item) .or. allocated(rhs%a_nm%item)) then
            if (allocated(lhs%a_nm%item) .and. allocated(rhs%a_nm%item)) then
                if (lhs%a_nm%item(1)%item(1) .ne. rhs%a_nm%item(1)%item(1)) then
                    rv = .false.
                end if
            else
                if (sum(lhs%a_nm) .ne. 0d0 .or. sum(rhs%a_nm) .ne. 0) then
                    rv = .false.
                end if
            end if
        end if

        if (allocated(lhs%b_nm%item) .or. allocated(rhs%b_nm%item)) then
            if (allocated(lhs%b_nm%item) .and. allocated(rhs%b_nm%item)) then
                if (lhs%b_nm%item(1)%item(1) .ne. rhs%b_nm%item(1)%item(1)) then
                    rv = .false.
                end if
            else
                if (sum(lhs%b_nm) .ne. 0d0 .or. sum(rhs%b_nm) .ne. 0) then
                    rv = .false.
                end if
            end if
        end if

    end function test_coefficient_eq

    function test_coefficient_ne(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: lhs
        type(lib_ml_fmm_coefficient), intent(in) :: rhs
        logical :: rv

        if (test_coefficient_eq(lhs, rhs)) then
            rv = .false.
        else
            rv = .true.
        end if

    end function test_coefficient_ne

    subroutine test_set_coefficient_zero(coefficient)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (allocated(coefficient%r)) then
            coefficient%r(:) = 0
        else
            allocate(coefficient%r(1))
            coefficient%r(1) = 0
        end if

        if (allocated(coefficient%c)) then
            coefficient%c(:) = 0
        else
            allocate(coefficient%c(1))
            coefficient%c(1) = 0
        end if

    end subroutine

    function test_v_operator_add(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
        type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
        type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

        ! auxilary
        integer :: i

        do i=1, size(lhs)
            rv(i)%r = lhs(i)%r + rhs(i)%r
        end do
    end function

    function test_v_operator_add_0d(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        allocate(rv%r(1))
        rv%r(1) = lhs%r(1) + rhs%r(1)
    end function

    function test_v_operator_sub(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
        type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
        type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

        ! auxilary
        integer :: i

        do i=1, size(lhs)
            allocate(rv(i)%r, source = lhs(i)%r - rhs(i)%r)
!            rv(i)%r = lhs(i)%r - rhs(i)%r
        end do
    end function

    function test_v_operator_sub_0d(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        allocate(rv%r, source = lhs%r - rhs%r)
!        rv%r = lhs%r - rhs%r
    end function

    ! Argument
    ! ----
    !   x: type(lib_tree_spatial_point)
    !       normalised point of the box centre
    !   data_element: type(lib_tree_data_element)
    !       internal "Tree" data element
    !   element_number: integer
    !       number of the element: m_ml_fmm_u(element_number)
    !
    ! dummyesult
    ! ----
    !   u_B_i: type(lib_ml_fmm_coefficient)
    !       calculation of one summand of the sum eq. 32
    !
    ! Visualisation
    ! ----
    !
    !  own data structure with informaion of each element:
    !       - position x
    !       - e.g. radius r
    !       - ...
    !
    !   -------------------------------------------------
    !   | x, r, ... | x, r, ... | x, r, ... | x, r, ... |
    !   -------------------------------------------------
    !         1           2          ...          N
    !
    !   tree data structure
    !       - normalised position ^x
    !       - hierarchy type h
    !
    !   -------------------------------------------------
    !   | ^x, h     | ^x, h     | ^x, h     | ^x, h     |
    !   -------------------------------------------------
    !      ^  1           2          ...          N     <-- element number
    !      |
    !      "data_element" with the number "element_number" (aka no)
    !
    !   Matrix Vector representation:
    !       asdf
    !         -       -   -  -     -  -
    !        |         | |    |   |    |
    !        | _______ | |    |   |    |
    !        | __XXO__ | | no | = |    | ;
    !        |         | |    |   |    |
    !        |         | |    |   |    |
    !         -       -   -  -     -  -
    !
    !   Alternative representation:
    !
    !       transformation of element "no"
    !         ---------
    !         |     no|
    !         |  ^x   |
    !         |       |
    !         ---------
    !
    ! HINT
    ! ----
    !   unscale position with: x_unscaled = lib_tree_get_unscaled_point(x)
    !
    !
    !
    ! dummyeference: Data_Structures_Optimal_Choice_of_Parameters_and_C, Gumerov, Duraiswami, Borovikov
    !            e.q. 32
    function test_get_u_B_i(x, data_element, element_number) result(u_B_i)
        use lib_tree_public
        use ml_fmm_type
        use lib_ml_fmm_data_container
        implicit none
        ! dummy
        type(lib_tree_spatial_point), intent(in) :: x
        type(lib_tree_data_element), intent(in) :: data_element
        integer(kind=4), intent(in) :: element_number

        type(lib_ml_fmm_coefficient) :: u_B_i

!        allocate(B_i%r, source = (/data_element%uindex%n + abs(x)/))
        allocate(u_B_i%r, source = (/real(data_element%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))

        allocate(u_B_i%c, source =(/cmplx(data_element%uindex%n, 0, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))

        allocate(u_B_i%a_nm%item(1))
        allocate(u_B_i%a_nm%item(1)%item, source =(/cmplx(0, data_element%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))

        allocate(u_B_i%b_nm%item(1))
        allocate(u_B_i%b_nm%item(1)%item, source =(/cmplx(data_element%uindex%n, 0, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))

    end function

    function test_get_u_phi_i_j(data_element_i, element_number_i, data_element_j, element_number_j) result(rv)
        use lib_tree_public
        use ml_fmm_type
        use lib_ml_fmm_data_container
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(in) :: data_element_i
        integer(kind=4), intent(in) :: element_number_i
        type(lib_tree_data_element), intent(in) :: data_element_j
        integer(kind=4), intent(in) :: element_number_j
        type(lib_ml_fmm_v) :: rv

!        allocate(rv%r, source = (/data_element_i%uindex%n + abs(y_j)/))
        allocate(rv%r(1))
        rv%r(1) = real(data_element_i%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)

        allocate(rv%c(1))
        rv%c(1) = cmplx(data_element_i%uindex%n, 0, kind=LIB_ML_FMM_COEFFICIENT_KIND)

        allocate(rv%a_nm%item(1))
        allocate(rv%a_nm%item(1)%item(1))
        rv%a_nm%item(1)%item(1) = cmplx(0, data_element_i%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)

        allocate(rv%b_nm%item(1))
        allocate(rv%b_nm%item(1)%item(1))
        rv%b_nm%item(1)%item(1) = cmplx(data_element_i%uindex%n, 0, kind=LIB_ML_FMM_COEFFICIENT_KIND)

    end function

    function test_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: A_i_2


!        allocate(A_i_2%r, source = (/A_i_1%r + abs(x_2 - x_1)/))
        allocate(A_i_2%r, source = (/A_i_1%r(1)/))

        allocate(A_i_2%c, source = A_i_1%c)

        allocate(A_i_2%a_nm%item, source = A_i_1%a_nm%item)

        allocate(A_i_2%b_nm%item, source = A_i_1%b_nm%item)

    end function

    function test_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: A_i_2

!        allocate(A_i_2%r, source = (/B_i_1%r + abs(x_2 - x_1)/))
        allocate(A_i_2%r, source = (/B_i_1%r(1)/))

        allocate(A_i_2%c, source = B_i_1%c)

        allocate(A_i_2%a_nm%item, source = B_i_1%a_nm%item)

        allocate(A_i_2%b_nm%item, source = B_i_1%b_nm%item)

    end function

    function test_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: B_i_2

!        allocate(B_i_2%r, source = (/B_i_1%r(1) + abs(x_2 - x_1)/))
        allocate(B_i_2%r, source = (/B_i_1%r(1)/))

        allocate(B_i_2%c, source = B_i_1%c)

        allocate(B_i_2%a_nm%item, source = B_i_1%a_nm%item)

        allocate(B_i_2%b_nm%item, source = B_i_1%b_nm%item)

    end function
end module ml_fmm_math
