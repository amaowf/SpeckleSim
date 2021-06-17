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

module lib_mie_ms_ml_fmm_interface
    use libmath

    use libtree
    use libmlfmm

    use lib_mie_type
    use lib_mie_type_functions
    use lib_mie_ms_solver_interface_helper_functions

    use lib_mie

    implicit none

    private

    public :: lib_mie_ms_ml_fmm_constructor
    public :: lib_mie_ms_ml_fmm_destructor
    public :: lib_mie_ms_ml_fmm_calculate_vector_b
    public :: lib_mie_ms_ml_fmm_calculate_evaluation_points

    contains

        ! HINT: simulation_data @ lib_mie_ms_data_container has to be initialised
        !       - sphere_list
        subroutine lib_mie_ms_ml_fmm_constructor()
            use lib_ml_fmm_type_operator
            use ml_fmm_math
            use lib_mie_ms_data_container
            implicit none
            ! dummy

            ! auxiliary
            integer :: i
            integer :: no

            type(ml_fmm_type_operator_procedures) :: operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures
            type(lib_ml_fmm_data) :: data_elements

            if (allocated(simulation_data%evaluation_points)) then
                allocate(data_elements%Y(size(simulation_data%evaluation_points)))

                do i = lbound(simulation_data%evaluation_points, 1), ubound(simulation_data%evaluation_points, 1)
                    no = i - lbound(simulation_data%evaluation_points, 1) + 1

                    data_elements%Y(no)%element_type = 1 ! value .ne. -1
                    data_elements%Y(no)%hierarchy = HIERARCHY_Y

                    data_elements%Y(no)%point_x%x(1) = simulation_data%evaluation_points(i)%coordinate%x
                    data_elements%Y(no)%point_x%x(2) = simulation_data%evaluation_points(i)%coordinate%y
                    data_elements%Y(no)%point_x%x(3) = simulation_data%evaluation_points(i)%coordinate%z
                end do
            end if

            allocate(data_elements%XY(size(simulation_data%sphere_list)))

            do i = lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                no = i - lbound(simulation_data%sphere_list, 1) + 1

                data_elements%XY(no)%element_type = 1 ! value .ne. -1
                data_elements%XY(no)%hierarchy = HIERARCHY_XY

                data_elements%XY(no)%point_x%x(1) = simulation_data%sphere_list(i)%d_0_j%x
                data_elements%XY(no)%point_x%x(2) = simulation_data%sphere_list(i)%d_0_j%y
                data_elements%XY(no)%point_x%x(3) = simulation_data%sphere_list(i)%d_0_j%z
            end do

            operator_procedures = ml_fmm_type_operator_get_procedures(1)
            ml_fmm_procedures = lib_mie_ms_ml_fmm_get_procedures()

            call lib_ml_fmm_constructor(data_elements, operator_procedures, ml_fmm_procedures, tree_s_opt = 10, &
                                        final_sum_calc_y_hierarchy = .false., final_sum_calc_xy_hierarchy = .true.)

        end subroutine lib_mie_ms_ml_fmm_constructor

        subroutine lib_mie_ms_ml_fmm_destructor
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine lib_mie_ms_ml_fmm_destructor

        ! Argument
        ! ----
        !   vector_x: double complex, dimension(:)
        !       vector x: M x = b
        !   vector_b: double complex, dimension(:)
        !       vector x:
        subroutine lib_mie_ms_ml_fmm_calculate_vector_b(vector_x, vector_b)
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i
            integer :: j

            integer :: counter

            integer :: first_sphere
            integer :: last_sphere

            integer(kind=4), dimension(2) :: n_range

            type(list_list_cmplx) :: a_nm
            type(list_list_cmplx) :: b_nm

            type(lib_ml_fmm_v), dimension(:), allocatable :: b
            type(lib_ml_fmm_v), dimension(:), allocatable :: x

            integer :: offset

            first_sphere = lbound(simulation_data%sphere_list, 1)
            last_sphere = ubound(simulation_data%sphere_list, 1)

            n_range = simulation_data%spherical_harmonics%n_range

            if (allocated(simulation_data%evaluation_points)) then
                offset = size(simulation_data%evaluation_points)
            else
                offset = 0
            end if

            allocate ( x(offset + size(simulation_data%sphere_list)) )

            !$OMP PARALLEL DO PRIVATE(j, i, a_nm, b_nm)
            do j = first_sphere, last_sphere
                i = j - first_sphere + 1

                call lib_mie_ms_solver_hf_get_list_list_cmplx_from_array(vector_x, i, n_range, &
                                                                         a_nm, b_nm)

                call move_alloc(a_nm%item, x(offset + i)%a_nm%item)
                call move_alloc(b_nm%item, x(offset + i)%b_nm%item)
            end do
            !$OMP END PARALLEL DO

            call lib_ml_fmm_run(x, b)

            counter = (1 + n_range(2))**2 - n_range(1)**2
            if (allocated(vector_b)) deallocate(vector_b)
            allocate(vector_b(2 * counter * size(simulation_data%sphere_list)))

            !$OMP PARALLEL DO PRIVATE(i, a_nm, b_nm)
            do i = 1, size(simulation_data%sphere_list)
                a_nm = b(offset + i)%a_nm
                b_nm = b(offset + i)%b_nm
                call lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array(a_nm, b_nm, i, n_range, vector_b)
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_ml_fmm_calculate_vector_b

        function lib_mie_ms_ml_fmm_get_procedures() result (handle)
!            use ml_fmm_type
            use lib_ml_fmm_type_operator
            implicit none
            ! dummy
            type(lib_ml_fmm_procedure_handles) :: handle

            handle%get_u_phi_i_j => lib_mie_ms_ml_fmm_get_u_phi_i_j
            handle%get_u_B_i => lib_mie_ms_ml_fmm_get_u_B_i
            handle%get_translation_SS => lib_mie_ms_ml_fmm_translation_SS
            handle%get_translation_SR => lib_mie_ms_ml_fmm_translation_SR
            handle%get_translation_RR => lib_mie_ms_ml_fmm_translation_RR
            handle%dor => lib_mie_ms_ml_fmm_dor

        end function lib_mie_ms_ml_fmm_get_procedures

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
        ! dummyreference: Data_Structures_Optimal_Choice_of_Parameters_and_C, Gumerov, Duraiswami, Borovikov
        !            e.q. 32
        function lib_mie_ms_ml_fmm_get_u_B_i(x, data_element, element_number) result(u_B_i)
            use libmath
            use lib_tree_public
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element), intent(in) :: data_element
            integer(kind=4), intent(in) :: element_number

            type(lib_ml_fmm_coefficient) :: u_B_i

            ! auxiliary
            integer :: i
            integer :: no

            integer :: n
            integer :: m

            integer(kind=1) :: z_selector

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            type(cartesian_coordinate_real_type) :: k_d_j_l
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_x_1_nm
            type(list_list_cmplx) :: buffer_x_2_nm

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            i = element_number - lbound(simulation_data%sphere_list, 1) + 1
            d_0_j = simulation_data%sphere_list(i)%d_0_j
            no = simulation_data%sphere_list(i)%sphere_parameter_index
            n_range_j = simulation_data%sphere_parameter_list(no)%n_range

            n_range = simulation_data%spherical_harmonics%n_range

            n_range_j(1) = max(n_range_j(1), n_range(1))
            n_range_j(2) = min(n_range_j(2), n_range(2))

            x_unscaled = lib_tree_get_unscaled_point(x)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                       / simulation_data%illumination%lambda_0 &
                       * d_j_l

            if (abs(d_j_l) .eq. 0d0) then
!                u_B_i%a_nm = simulation_data%sphere_list(i)%a_nm
!                u_B_i%b_nm = simulation_data%sphere_list(i)%b_nm
                u_B_i%a_nm = m_ml_fmm_u(element_number)%a_nm
                u_B_i%b_nm = m_ml_fmm_u(element_number)%b_nm

                return

            ! Xu - 1998 - Efficient Evaluation of Vector Translation Coeffic
            else if (abs(d_j_l) .gt. simulation_data%sphere_parameter_list(no)%radius) then
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r
            else
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r
            end if

            ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
            z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave

            call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                    n_range, n_range_j, z_selector, &
                    a_nmnumu, b_nmnumu)

!            buffer_x_1_nm = simulation_data%sphere_list(i)%a_nm
!            buffer_x_2_nm = simulation_data%sphere_list(i)%b_nm

            ! vector u aka vector x (solver_ Mx = b)
            buffer_x_1_nm = m_ml_fmm_u(element_number)%a_nm
            buffer_x_2_nm = m_ml_fmm_u(element_number)%b_nm

            call init_list(buffer_1_nm, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(buffer_2_nm, n_range(1), n_range(2) - n_range(1) + 1)

            !$OMP PARALLEL DO PRIVATE(n)
            do n = n_range(1), n_range(2)
                !$OMP PARALLEL DO PRIVATE(m)
                do m = -n, n
                    if (n_range(2) .le. n_range_j(2) &
                        .and. n_range(1) .ge. n_range_j(1)) then
                        buffer_1_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + b_nmnumu%item(n)%item(m) * buffer_x_2_nm)

                        buffer_2_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + a_nmnumu%item(n)%item(m) * buffer_x_2_nm)
                    else
                        buffer_1_nm%item(n)%item(m) = dcmplx(0,0)
                        buffer_2_nm%item(n)%item(m) = dcmplx(0,0)
                    end if
                end do
                !$OMP END PARALLEL DO
            end do
            !$OMP END PARALLEL DO

            call move_alloc(buffer_1_nm%item, u_B_i%a_nm%item)
            call move_alloc(buffer_2_nm%item, u_B_i%b_nm%item)

            if (allocated(buffer_x_1_nm%item)) deallocate(buffer_x_1_nm%item)
            if (allocated(buffer_x_2_nm%item)) deallocate(buffer_x_2_nm%item)

            if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
            if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)

        end function lib_mie_ms_ml_fmm_get_u_B_i

        function lib_mie_ms_ml_fmm_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: B_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            type(cartesian_coordinate_real_type) :: k_d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

!            z_selector = simulation_data%spherical_harmonics%z_selector_scatterd_wave

            ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
            z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                       / simulation_data%illumination%lambda_0 &
                       * d_j_l
            call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            call init_list(B_i_2%a_nm, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(B_i_2%b_nm, n_range(1), n_range(2) - n_range(1) + 1)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    B_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * B_i_1%b_nm)

                    B_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * B_i_1%b_nm)
                end do
            end do
            !$OMP END PARALLEL DO

            if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
            if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)

        end function lib_mie_ms_ml_fmm_translation_SS

        function lib_mie_ms_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            type(cartesian_coordinate_real_type) :: k_d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            !z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave

            ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
            z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                       / simulation_data%illumination%lambda_0 &
                       * d_j_l
            call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            call init_list(A_i_2%a_nm, n_range(1), n_range(2) - n_range(1) +1)
            call init_list(A_i_2%b_nm, n_range(1), n_range(2) - n_range(1) +1)

            !$OMP PARALLEL DO PRIVATE(n)
            do n = n_range(1), n_range(2)
                !$OMP PARALLEL DO PRIVATE(m)
                do m = -n, n
                    A_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * B_i_1%b_nm)

                    A_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * B_i_1%b_nm)
                end do
                !$OMP END PARALLEL DO
            end do
            !$OMP END PARALLEL DO

            if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
            if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)

        end function lib_mie_ms_ml_fmm_translation_SR

        function lib_mie_ms_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            type(cartesian_coordinate_real_type) :: k_d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
            z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                       / simulation_data%illumination%lambda_0 &
                       * d_j_l
            call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            call init_list(A_i_2%a_nm, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(A_i_2%b_nm, n_range(1), n_range(2) - n_range(1) + 1)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    A_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * A_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * A_i_1%b_nm)

                    A_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * A_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * A_i_1%b_nm)
                end do
            end do
            !$OMP END PARALLEL DO

            if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
            if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)

        end function lib_mie_ms_ml_fmm_translation_RR

        ! Argument
        ! ----
        !   data_element_i: type(lib_tree_data_element)
        !       tree specific data of the i-th element
        !   element_number_i: integer
        !       element number according to the tree list (1:N)
        !   data_element_j: type(lib_tree_data_element)
        !       tree specific data of the j-th element
        !   element_number_j: integer
        !       element number according to the tree list (1:N)
        !
        ! Returns
        ! ----
        !   rv: type(lib_ml_fmm_v)
        !       one summand of the sum, eq. 38 [1]
        !
        ! Hint
        ! ----
        !   Coordinates of the tree specific data variabels are scaled.
        !   Get the uncalled coordinated with the function call lib_tree_get_unscaled_point().
        !
        !
        ! Reference: [1] Data Structures, Optimal Choice of Parameters, and Complexity Results for
        !                Generalized Multilevel Fast Multipole Methods in d Dimensions,
        !                Nail Gumerov, Ramani Duraiswami, Eugene Borovikov
        function lib_mie_ms_ml_fmm_get_u_phi_i_j(data_element_i, element_number_i, &
                                                 data_element_j, element_number_j) result(rv)
            use lib_tree_public
            use ml_fmm_type
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_tree_data_element), intent(in) :: data_element_i
            integer(kind=4), intent(in) :: element_number_i
            type(lib_tree_data_element), intent(in) :: data_element_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv

            ! auxiliary
            integer :: j
            integer :: l

            integer :: n
            integer :: m

            integer :: no

!            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            type(spherical_coordinate_real_type) :: d_j_l_spherical
            type(cartesian_coordinate_real_type) :: k_d_j_l

            integer, dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l
            integer(kind=1) :: z_selector

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_b_1_nm
            type(list_list_cmplx) :: buffer_b_2_nm

            type(list_list_cmplx) :: buffer_x_1_nm
            type(list_list_cmplx) :: buffer_x_2_nm

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            type(spherical_coordinate_cmplx_type), dimension(2) :: field
            type(cartesian_coordinate_cmplx_type), dimension(2) :: field_cart



            if (data_element_j%hierarchy .eq. HIERARCHY_Y) then

                if (j .eq. l) then
                    print *, "lib_mie_ms_ml_fmm_get_u_phi_i_j: ERROR"
                    print *, "  HIERARCHY_Y: j .eq. l: ", j
                else

                    j = element_number_j - lbound(simulation_data%evaluation_points, 1) + 1
                    d_0_j = simulation_data%evaluation_points(j)%coordinate

                    l = element_number_i - lbound(simulation_data%sphere_list, 1) + 1

                    buffer_x_1_nm = m_ml_fmm_u(element_number_i)%a_nm
                    buffer_x_2_nm = m_ml_fmm_u(element_number_i)%b_nm

                    d_0_l = simulation_data%sphere_list(l)%d_0_j

                    !z_selector = simulation_data%spherical_harmonics%z_selector_scatterd_wave

                    ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
                    z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave

                    d_j_l = d_0_l - d_0_j
                    d_j_l_spherical = d_j_l

                    field = lib_mie_get_field(d_j_l_spherical%theta, &
                                              d_j_l_spherical%phi, &
                                              d_j_l_spherical%rho, &
                                              simulation_data%illumination%e_field_0, &
                                              simulation_data%illumination%lambda_0, &
                                              simulation_data%refractive_index_medium, &
                                              buffer_x_1_nm, buffer_x_2_nm, &
                                              z_selector=z_selector)

                    field_cart(1) = make_cartesian(field(1), d_j_l_spherical)
                    field_cart(2) = make_cartesian(field(2), d_j_l_spherical)

                    allocate(rv%c(6))
                    rv%c(1) = field_cart(1)%x
                    rv%c(2) = field_cart(1)%y
                    rv%c(3) = field_cart(1)%z

                    rv%c(4) = field_cart(2)%x
                    rv%c(5) = field_cart(2)%y
                    rv%c(6) = field_cart(2)%z

                end if

            else if (data_element_j%hierarchy .eq. HIERARCHY_XY) then

                j = element_number_j - lbound(simulation_data%sphere_list, 1) + 1

                l = element_number_i - lbound(simulation_data%sphere_list, 1) + 1
                buffer_x_1_nm = m_ml_fmm_u(element_number_i)%a_nm
                buffer_x_2_nm = m_ml_fmm_u(element_number_i)%b_nm

                if (j .eq. l) then
                    buffer_b_1_nm = buffer_x_1_nm
                    buffer_b_2_nm = buffer_x_2_nm
                else
                    d_0_j = simulation_data%sphere_list(j)%d_0_j
                    no = simulation_data%sphere_list(j)%sphere_parameter_index
                    n_range_j = simulation_data%sphere_parameter_list(no)%n_range

                    a_n = simulation_data%sphere_parameter_list(no)%a_n
                    b_n = simulation_data%sphere_parameter_list(no)%b_n

                    n_range = simulation_data%spherical_harmonics%n_range
                    z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave


                    d_0_l = simulation_data%sphere_list(l)%d_0_j
                    no = simulation_data%sphere_list(l)%sphere_parameter_index
                    n_range_l = simulation_data%sphere_parameter_list(no)%n_range

                    n_range_j(1) = max(n_range_j(1), n_range(1))
                    n_range_j(2) = min(n_range_j(2), n_range(2))

                    n_range_l(1) = max(n_range_l(1), n_range(1))
                    n_range_l(2) = min(n_range_l(2), n_range(2))

                    d_j_l = d_0_l - d_0_j
                    k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                           / simulation_data%illumination%lambda_0 &
                           * d_j_l
                    call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                            n_range_j, n_range_l, z_selector, &
                            a_nmnumu, b_nmnumu)

                    call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                    call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                    call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                    call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                    ! calculate element j of vector v aka vector b (solver: Mx=b)

                    !$OMP PARALLEL DO PRIVATE(n, m) &
                    !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                    do n = n_range(1), n_range(2)
                        !$OMP PARALLEL DO PRIVATE(m) &
                        !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                        do m = -n, n
                            if (n_range(2) .le. n_range_j(2) &
                                .and. n_range(1) .ge. n_range_j(1)) then
                                buffer_1_nm = a_nmnumu%item(n)%item(m)
                                buffer_2_nm = b_nmnumu%item(n)%item(m)

                                buffer_b_1_nm%item(n)%item(m) = a_n%item(n) * sum(buffer_1_nm * buffer_x_1_nm &
                                                                    + buffer_2_nm * buffer_x_2_nm)

                                buffer_b_2_nm%item(n)%item(m) = b_n%item(n) * sum(buffer_2_nm * buffer_x_1_nm &
                                                                    + buffer_1_nm * buffer_x_2_nm)
                            else
                                buffer_b_1_nm%item(n)%item(m) = dcmplx(0,0)
                                buffer_b_2_nm%item(n)%item(m) = dcmplx(0,0)
                            end if
                        end do
                        !$OMP END PARALLEL DO
                    end do
                    !$OMP END PARALLEL DO
                end if

                call move_alloc(buffer_b_1_nm%item, rv%a_nm%item)
                call move_alloc(buffer_b_2_nm%item, rv%b_nm%item)

                if (allocated(buffer_1_nm%item)) deallocate(buffer_1_nm%item)
                if (allocated(buffer_2_nm%item)) deallocate(buffer_2_nm%item)

                if (allocated(a_n%item)) deallocate(a_n%item)
                if (allocated(b_n%item)) deallocate(b_n%item)

                if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
                if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)
            end if

        end function lib_mie_ms_ml_fmm_get_u_phi_i_j

        ! Argument
        ! ----
        !   D: type(lib_ml_fmm_coefficient)
        !
        !
        ! Hint
        ! ----
        !   Coordinates of the tree specific data variabels are scaled.
        !   Get the uncalled coordinated with the function call lib_tree_get_unscaled_point().
        !
        function lib_mie_ms_ml_fmm_dor(D, x_c, data_element_j, element_number_j)  result (rv)
            use lib_tree_public
            use ml_fmm_type
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: D
            type(lib_tree_spatial_point), intent(in) :: x_c
            type(lib_tree_data_element), intent(in) :: data_element_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv

            ! auxiliary
            integer :: j

            integer :: n
            integer :: m

            integer :: no

            type(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_j_l
            type(spherical_coordinate_real_type) :: d_j_l_spherical
            type(cartesian_coordinate_real_type) :: k_d_j_l

            integer, dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l
            integer(kind=1) :: z_selector

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_b_1_nm
            type(list_list_cmplx) :: buffer_b_2_nm

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            type(spherical_coordinate_cmplx_type), dimension(2) :: field
            type(cartesian_coordinate_cmplx_type), dimension(2) :: field_cart

            x_unscaled = lib_tree_get_unscaled_point(data_element_j%point_x - x_c)

            d_j_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            if (data_element_j%hierarchy .eq. HIERARCHY_Y) then

                j = element_number_j - lbound(simulation_data%evaluation_points, 1) + 1

                !z_selector = simulation_data%spherical_harmonics%z_selector_scatterd_wave

                ! Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
                z_selector = simulation_data%spherical_harmonics%z_selector_incident_wave

                d_j_l_spherical = d_j_l

                field = lib_mie_get_field(d_j_l_spherical%theta, &
                                          d_j_l_spherical%phi, &
                                          d_j_l_spherical%rho, &
                                          simulation_data%illumination%e_field_0, &
                                          simulation_data%illumination%lambda_0, &
                                          simulation_data%refractive_index_medium, &
                                          D%a_nm, D%b_nm, &
                                          z_selector=z_selector)

                field_cart(1) = make_cartesian(field(1), d_j_l_spherical)
                field_cart(2) = make_cartesian(field(2), d_j_l_spherical)

                allocate(rv%c(6))
                rv%c(1) = field_cart(1)%x
                rv%c(2) = field_cart(1)%y
                rv%c(3) = field_cart(1)%z

                rv%c(4) = field_cart(2)%x
                rv%c(5) = field_cart(2)%y
                rv%c(6) = field_cart(2)%z

            else if (data_element_j%hierarchy .eq. HIERARCHY_XY) then

                k_d_j_l = 2 * PI * simulation_data%refractive_index_medium &
                           / simulation_data%illumination%lambda_0 &
                           * d_j_l

                j = element_number_j - lbound(simulation_data%sphere_list, 1) + 1
                no = simulation_data%sphere_list(j)%sphere_parameter_index
                n_range_j = simulation_data%sphere_parameter_list(no)%n_range

                a_n = simulation_data%sphere_parameter_list(no)%a_n
                b_n = simulation_data%sphere_parameter_list(no)%b_n

                n_range_l(1) = lbound(d%a_nm%item, 1)
                n_range_l(2) = ubound(d%a_nm%item, 1)

                n_range = simulation_data%spherical_harmonics%n_range

                z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r

                n_range_j(1) = max(n_range_j(1), n_range(1))
                n_range_j(2) = min(n_range_j(2), n_range(2))

                n_range_l(1) = max(n_range_l(1), n_range(1))
                n_range_l(2) = min(n_range_l(2), n_range(2))

                if (abs(d_j_l) .eq. 0d0) then
                    buffer_b_1_nm = D%a_nm
                    buffer_b_2_nm = D%b_nm
                else
                    call lib_math_vector_spherical_harmonics_translation_coefficient(k_d_j_l, &
                            n_range_j, n_range_l, z_selector, &
                            a_nmnumu, b_nmnumu)

                    call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                    call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                    call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                    call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                    ! calculate element j of vector v aka vector b (solver: Mx=b)

                    !$OMP PARALLEL DO PRIVATE(n, m) &
                    !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                    do n = n_range(1), n_range(2)
                        do m = -n, n
                            if (n_range(2) .le. n_range_j(2) &
                                .and. n_range(1) .ge. n_range_j(1)) then
                                buffer_1_nm = a_nmnumu%item(n)%item(m)
                                buffer_2_nm = b_nmnumu%item(n)%item(m)

                                buffer_b_1_nm%item(n)%item(m) = a_n%item(n) * sum(buffer_1_nm * D%a_nm&
                                                                                  + buffer_2_nm * D%b_nm)

                                buffer_b_2_nm%item(n)%item(m) = b_n%item(n) * sum(buffer_2_nm * D%a_nm&
                                                                                  + buffer_1_nm * D%b_nm)
                            else
                                buffer_b_1_nm%item(n)%item(m) = dcmplx(0,0)
                                buffer_b_2_nm%item(n)%item(m) = dcmplx(0,0)
                            end if
                        end do
                    end do
                    !$OMP END PARALLEL DO
                end if

                rv%a_nm = buffer_b_1_nm
                rv%b_nm = buffer_b_2_nm

                if (allocated(a_nmnumu%item)) deallocate(a_nmnumu%item)
                if (allocated(b_nmnumu%item)) deallocate(b_nmnumu%item)
            end if

        end function lib_mie_ms_ml_fmm_dor

        subroutine lib_mie_ms_ml_fmm_calculate_evaluation_points()
            use libmath
            use lib_ml_fmm
            use lib_mie_ms_data_container
            implicit none
            ! dummy

            ! auxiliary
            integer :: i
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v
            type(cartesian_coordinate_cmplx_type), dimension(2) :: field

            if (allocated(simulation_data%evaluation_points)) then
                call lib_ml_fmm_final_summation(calculate_y_hierarchy = .true., &
                                                calculate_xy_hierarchy = .false.)

                vector_v = lib_ml_fmm_get_vector_v()

                do i = 1, size(simulation_data%evaluation_points)
                    field(1)%x = vector_v(i)%c(1)
                    field(1)%y = vector_v(i)%c(2)
                    field(1)%z = vector_v(i)%c(3)

                    field(2)%x = vector_v(i)%c(4)
                    field(2)%y = vector_v(i)%c(5)
                    field(2)%z = vector_v(i)%c(6)

                    simulation_data%evaluation_points(i)%e_field = field(1)
                    simulation_data%evaluation_points(i)%h_field = field(2)
                end do
            end if

        end subroutine lib_mie_ms_ml_fmm_calculate_evaluation_points

end module lib_mie_ms_ml_fmm_interface
