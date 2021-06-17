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

module lib_sie_ml_fmm_interface
    use lib_tree
    use lib_tree_type
    use lib_ml_fmm
    use lib_ml_fmm_type_operator
    use lib_ml_fmm_helper_functions
    use ml_fmm_type
    use ml_fmm_math

    implicit none

    private

    public :: lib_sie_constructor
    public :: lib_sie_destructor

    public :: lib_sie_ml_fmm_calculate_vector_b

    contains

        ! Argument
        ! ----
        !   data_elements: type(lib_ml_fmm_data)
        !       positions of the data elements of the X, Y, and XY hierarchy
        !   s_opt: integer
        !       maximum number of elements at box (n, l_max)
        !
        subroutine lib_sie_constructor(data_elements, s_opt)
            implicit none
            ! dummy
            type(lib_ml_fmm_data), intent(inout) :: data_elements
            integer, intent(in) :: s_opt

            ! auxiliary

            type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures


            ml_fmm_operator_procedures = ml_fmm_type_operator_get_procedures(1)
            ml_fmm_procedures = lib_sie_ml_fmm_get_procedures()

            call lib_ml_fmm_constructor(data_elements, &
                                        ml_fmm_operator_procedures, &
                                        ml_fmm_procedures, &
                                        tree_s_opt=s_opt, &
                                        final_sum_calc_y_hierarchy = .false., &
                                        final_sum_calc_xy_hierarchy = .true. ,&
                                        use_own_sum=.true.)

        end subroutine

        subroutine lib_sie_destructor()
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine

        ! Argument
        ! ----
        !   vector_x: double complex, dimension(:)
        !       vector x: M x = b
        !   vector_b: double complex, dimension(:)
        !       vector x:
        subroutine lib_sie_ml_fmm_calculate_vector_b(vector_x, vector_b)
!            use lib_sie_data_container
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: b
            type(lib_ml_fmm_v), dimension(:), allocatable :: x

            ! todo: pre-process: x

            call lib_ml_fmm_run(x, b)

            ! todo: post-process: b

        end subroutine


        function lib_sie_ml_fmm_get_procedures() result(handle)
            implicit none
            ! dummy
            type(lib_ml_fmm_procedure_handles) :: handle

            ! Upward Pass
            ! Step 1
            handle%get_c => lib_sie_ml_fmm_get_C ! todo: replace null() with own function
            ! Step 2
            handle%get_translation_SS => null() ! todo: replace null() with own function

            ! Downward Pass
            ! Step 1
            handle%get_translation_SR => null() ! todo: replace null() with own function
            ! Step 2
            handle%get_translation_RR => null() ! todo: replace null() with own function

            ! Final Summation
            handle%get_v_y_j => null() ! todo: replace null() with own function

        end function

        ! todo: define procedures
        ! - get_c
        ! - get_translation_SS
        ! - get_translation_SR
        ! - get_translation_RR
        ! - get_v_y_j

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
        function lib_sie_ml_fmm_get_C(x, data_element, element_number) result(C)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element), dimension(:), intent(in) :: data_element
            integer(kind=4), dimension(:), intent(in) :: element_number

            type(lib_ml_fmm_coefficient) :: C
        end function
end module lib_sie_ml_fmm_interface
