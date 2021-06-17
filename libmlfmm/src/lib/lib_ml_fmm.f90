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

!#define _TRACE_

#define _FMM_DIMENSION_ 3

! LIB: Mulitlevel Fast Multipole Method
!
module lib_ml_fmm
    use libtree
    use ml_fmm_type
    use ml_fmm_math
    use lib_ml_fmm_type_operator
    use lib_ml_fmm_helper_functions
    use lib_ml_fmm_data_container
    implicit none

    private

    ! --- public functions ---
    public :: lib_ml_fmm_constructor
    public :: lib_ml_fmm_destructor
    public :: lib_ml_fmm_run
    public :: lib_ml_fmm_final_summation
    public :: lib_ml_fmm_get_vector_v

    public :: lib_ml_fmm_get_neighbourhood_size_k

    public :: lib_ml_fmm_test_functions

    public :: lib_ml_fmm_hf_test_functions

    ! --- member ---
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles

    type(lib_ml_fmm_hierarchy), dimension(:), allocatable, public :: m_ml_fmm_hierarchy

    ! Tree parameters
    integer(kind=UINDEX_BYTES) :: m_tree_neighbourhood_size_k
    integer(kind=4) :: m_tree_s_opt
    integer(kind=1) :: m_tree_l_min
    integer(kind=1) :: m_tree_l_max

    integer(kind=UINDEX_BYTES), dimension(3) :: m_data_element_hierarchy_type_length

    logical :: m_final_sum_calc_y_hierarchy
    logical :: m_final_sum_calc_xy_hierarchy
    logical :: m_use_own_sum

    contains

    ! Argument
    ! ----
    !   tree_s_opt: integer, optional (std: 1)
    !       Maximum number of elements in a box at the lowest level (l_max). *tree_l_max* is
    !       calculated based on this parameter.
    !   tree_bounding_box: lib_tree_spatial_point, dimension(2), optional
    !       Defines two opposite coreners of the box at level 0. Without this argument,
    !       the box at level 0 is automatically calculated.
    !       HINT:
    !           - All elements has to be inside of the box.
    !           - bounding_box(1)%x(i) < bounding_box(1)%x(i); i = 1..TREE_DIMENSION
    !   tree_l_max: integer, optional
    !       Specifies the lowest level.
    !       Can be used as an alternaive to the *s* parameter.
    !
    ! Procedure
    ! ----
    !   1. call lib_tree_constructor
    !       - determine s_opt
    !       - determine l_min and l_max
    !   2. create ml fmm data set
    !       - C per box and per level (l_max up to l_min)
    !       - D per box and per level (l_min down to l_max)
    subroutine lib_ml_fmm_constructor(data_elements, operator_procedures, ml_fmm_procedures, &
                                      tree_s_opt, tree_bounding_box, tree_l_max, &
                                      final_sum_calc_y_hierarchy, final_sum_calc_xy_hierarchy, &
                                      use_own_sum)
        implicit none
        ! dummy
        type(lib_ml_fmm_data), intent(inout) :: data_elements
        type(ml_fmm_type_operator_procedures), intent(in) :: operator_procedures
        type(lib_ml_fmm_procedure_handles), intent(in) :: ml_fmm_procedures

        integer, intent(in), optional :: tree_s_opt
        type(lib_tree_spatial_point), dimension(2), intent(in), optional :: tree_bounding_box
        integer(kind=1), intent(in), optional :: tree_l_max

        logical, intent(in), optional :: final_sum_calc_y_hierarchy
        logical, intent(in), optional :: final_sum_calc_xy_hierarchy
        logical, intent(in), optional :: use_own_sum

        ! auxiliaray
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated
        type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: correspondence_vector
        integer(kind=UINDEX_BYTES), dimension(3) :: length

        type(lib_tree_spatial_point), dimension(2) :: m_tree_bounding_box

        m_final_sum_calc_y_hierarchy = .true.
        if (present(final_sum_calc_y_hierarchy)) m_final_sum_calc_y_hierarchy = final_sum_calc_y_hierarchy

        m_final_sum_calc_xy_hierarchy = .true.
        if (present(final_sum_calc_xy_hierarchy)) m_final_sum_calc_xy_hierarchy = final_sum_calc_xy_hierarchy

        m_use_own_sum = .false.
        if (present(use_own_sum)) m_use_own_sum = use_own_sum

        ! concatenate X-, Y- and XY-hierarchy data
        !   length(1) = size(data_elements%X)
        !   length(2) = size(data_elements%Y)
        !   length(3) = size(data_elements%XY)
        allocate(data_concatenated, source=lib_ml_fmm_concatenate_data_array(data_elements, length))

        m_data_element_hierarchy_type_length = length

        m_tree_s_opt = 1
        if (present(tree_s_opt)) m_tree_s_opt = int(tree_s_opt, 4)

        if (present(tree_bounding_box)) then
            m_tree_bounding_box = tree_bounding_box
        else
            m_tree_bounding_box(1)%x(:) = 0
            m_tree_bounding_box(2)%x(:) = 0
        end if

        ! initiate the Tree
        if (present(tree_l_max)) then
            call lib_tree_constructor(data_concatenated, &
                                      l_max = tree_l_max, &
                                      bounding_box = m_tree_bounding_box)
            m_tree_l_max = tree_l_max
        else
            call lib_tree_constructor(data_concatenated, &
                                      s = m_tree_s_opt, &
                                      bounding_box = m_tree_bounding_box)
            m_tree_l_max = lib_tree_get_level_max(m_tree_s_opt)
        end if


        data_concatenated = lib_tree_get_element_list()

        ! initiate the ml fmm type operators
!        operator_procedures = ml_fmm_type_operator_get_procedures()
        call lib_ml_fmm_type_operator_constructor(operator_procedures)

        ! initiate the ml fmm procedures
        m_ml_fmm_handles = ml_fmm_procedures

        ! adjust Tree parameters
        m_tree_neighbourhood_size_k = 1!lib_ml_fmm_hf_get_neighbourhood_size(R_c1, r_c2)
        m_tree_l_min = lib_tree_get_level_min(m_tree_neighbourhood_size_k)

        if (m_tree_l_max .lt. m_tree_l_min) m_tree_l_min = m_tree_l_max

        ! initiate the X- and Y-hierarchy
        correspondence_vector = lib_tree_get_correspondence_vector()
        call lib_ml_fmm_hf_create_hierarchy(data_concatenated, &
                                            length, m_tree_l_min, m_tree_l_max,&
                                            m_ml_fmm_hierarchy)
        print *, "lib_ml_fmm_constructor: INFO"
        print *, "  m_tree_l_min: ", m_tree_l_min
        print *, "  m_tree_l_max: ", m_tree_l_max
        print *, "  m_tree_s_opt: ", m_tree_s_opt

        if (allocated(m_ml_fmm_u)) then
            deallocate(m_ml_fmm_u)
        end if
        if (allocated(m_ml_fmm_v)) then
            deallocate(m_ml_fmm_v)
        end if

        allocate (m_ml_fmm_u(sum(length)))
        allocate (m_ml_fmm_v(sum(length)))
    end subroutine lib_ml_fmm_constructor

    ! clean up
    ! - coefficents
    ! - Tree
    subroutine lib_ml_fmm_destructor()
        implicit none

        ! auxilary
        ! example for an indiviual deallocation
!        integer(kind=1) :: i
!        integer(kind=4) :: ii

        if (allocated(m_ml_fmm_hierarchy)) then
            deallocate(m_ml_fmm_hierarchy)
        end if

        if (allocated(m_ml_fmm_u)) then
            deallocate(m_ml_fmm_u)
        end if

        if (allocated(m_ml_fmm_v)) then
            deallocate(m_ml_fmm_v)
        end if

!        if ( allocated(m_ml_fmm_C) ) then
!            ! HINT: m_ml_fmm_C has a deep structure. If the automatic deallocation fails,
!            !       each subelement should be deallocated individually.
!            deallocate (m_ml_fmm_C)
!            ! example for an indiviual deallocation
!!            do i=1, size(m_ml_fmm_C)
!!                do ii=1, size(m_ml_fmm_C(i)%r)
!!                    if (allocated (m_ml_fmm_C(i)%r))
!!                        deallocate (m_ml_fmm_C(i)%r)
!!                    end if
!!                end do
!!            end do
!        end if
!
!        if ( allocated(m_ml_fmm_D) ) then
!            ! HINT: m_ml_fmm_D has a deep structure. If the automatic deallocation fails,
!            !       each subelement should be deallocated individually.
!            !       (as like m_ml_fmm_C)
!            deallocate(m_ml_fmm_D)
!        end if

        call lib_tree_destructor()

    end subroutine lib_ml_fmm_destructor

    ! Argument
    ! ----
    !   vector_u: type(lib_ml_fmm_v), dimension(:)
    !       vector u of the matrix vector product PHI * u = v
    !       HINT: length(vector_u) = no of all data elements (X-, Y- and XY-, hierarchy)
    !       CONVENTION:
    !           Internal representation of the element data list
    !           from the 1-st to the N-th element.
    !           ------------------------------------
    !           |1    X    |     Y    |     XY    N|
    !           ------------------------------------
    !           X-, Y-, XY- hierarchy
    !   use_move_alloc_vector_u: logical, optional (std: .true.)
    !
    ! Returns
    ! ----
    !   vector_v: type(lib_ml_fmm_v), dimension(:)
    !       vector v  of the matrix vector product PHI * u = v
    !       HINT: length(vector_u) = no of all data elements (X-, Y- and XY-, hierarchy)
    !       CONVENTION:
    !           Internal representation of the element data list
    !           from the 1-st to the N-th element.
    !           ------------------------------------
    !           |1    X    |     Y    |     XY    N|
    !           ------------------------------------
    !           X-, Y-, XY- hierarchy
    subroutine lib_ml_fmm_run(vector_u, vector_v, use_move_alloc_vector_u)
        implicit none
        ! dummy
        type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_u
        type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_v

        logical, intent(in), optional :: use_move_alloc_vector_u

        ! auxiliary
        logical :: m_use_move_alloc_vector_u

        m_use_move_alloc_vector_u = .true.
        if (present(use_move_alloc_vector_u)) m_use_move_alloc_vector_u = use_move_alloc_vector_u

        if (m_use_move_alloc_vector_u) then
            call move_alloc(vector_u, m_ml_fmm_u)
        else
            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if
        end if

        call lib_ml_fmm_calculate_upward_pass()
        call lib_ml_fmm_calculate_downward_pass()
        call lib_ml_fmm_final_summation(m_final_sum_calc_y_hierarchy, m_final_sum_calc_xy_hierarchy)

        allocate(vector_v, source=m_ml_fmm_v)
!        call move_alloc(m_ml_fmm_v, vector_v)

    end subroutine lib_ml_fmm_run

    function lib_ml_fmm_get_vector_v() result(rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_v), dimension(:), allocatable :: rv

        allocate(rv, source=m_ml_fmm_v)

    end function

    ! Concatenates the arrays at the data_elements to one array of the type lib_tree_data_element
    !
    ! Argument
    ! ----
    !   data_elements: lib_ml_fmm_data
    !       lists of lib_tree_data_element arrays
    !   length: integer array, out
    !       size of each lib_tree_data_element array
    !
    ! Returns
    ! ----
    !   data_concatenated: lib_tree_data_element array
    !       concatenation of the lib_tree_element arrays of data_elements
    !
    function lib_ml_fmm_concatenate_data_array(data_elements, length) result(data_concatenated)
        implicit none
        type(lib_ml_fmm_data), intent(inout) :: data_elements
        integer(kind=UINDEX_BYTES), dimension(3), intent(out) :: length
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated

        ! auxiliaray
        integer(kind=1) :: i
        integer(kind=UINDEX_BYTES) :: start, last

        length(:) = 0
        if( allocated(data_elements%X) ) then
            length(HIERARCHY_X) = size(data_elements%X)
        end if

        if( allocated(data_elements%Y) ) then
            length(HIERARCHY_Y) = size(data_elements%Y)
        end if

        if( allocated(data_elements%XY) ) then
            length(HIERARCHY_XY) = size(data_elements%XY)
        end if

        ! concatenate all three datasets
        allocate( data_concatenated(sum(length)))

        do i=1, 3
            if (length(i) .gt. 0) then
                if (i .eq. 1) then
                    start = 1
                else
                    start = sum(length(:i-1)) + 1
                end if
                last = start + length(i) - 1
                if (i .eq. HIERARCHY_X) then
                    data_concatenated(start:last) = data_elements%X
                    data_concatenated(start:last)%hierarchy = HIERARCHY_X
                else if (i .eq. HIERARCHY_Y) then
                    data_concatenated(start:last) = data_elements%Y
                    data_concatenated(start:last)%hierarchy = HIERARCHY_Y
                else if (i .eq. HIERARCHY_XY) then
                    data_concatenated(start:last) = data_elements%XY
                    data_concatenated(start:last)%hierarchy = HIERARCHY_XY
                end if
            end if
        end do
    end function



    subroutine lib_ml_fmm_calculate_upward_pass()
        implicit none
#ifdef _TRACE_
        ! CPU-time
        real :: test_start_sub, test_finish_sub
        ! WALL-time
        INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

        call system_clock(test_count_start_sub, test_count_rate_sub)
        call cpu_time(test_start_sub)
#endif

        call lib_ml_fmm_calculate_upward_pass_step_1()
#ifdef _TRACE_
        call cpu_time(test_finish_sub)
        call system_clock(test_count_finish_sub, test_count_rate_sub)

        print *, "TRACE: lib_ml_fmm_calculate_upward_pass_step_1 ..done"
        print '("    CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
        print '("    WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)

        call system_clock(test_count_start_sub, test_count_rate_sub)
        call cpu_time(test_start_sub)
#endif

        call lib_ml_fmm_calculate_upward_pass_step_2()
#ifdef _TRACE_
        call cpu_time(test_finish_sub)
        call system_clock(test_count_finish_sub, test_count_rate_sub)

        print *, "TRACE: lib_ml_fmm_calculate_upward_pass_step_2 ..done"
        print '("    CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
        print '("    WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)
#endif
    end subroutine

    ! Upward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_max
    !       1. get all elements of the X-hierarchy
    !       2. calc B coefficients
    !       3. multiply with u_i
    !       4. calculate the sum of all elements of a box
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(32)
    !
    subroutine lib_ml_fmm_calculate_upward_pass_step_1()
        use lib_tree
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient) :: C

        ! auxiliaray
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        type(lib_tree_universal_index) :: uindex
        logical :: ignore_box

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        number_of_boxes = size(m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index)

        uindex%l = m_tree_l_max
        if (m_ml_fmm_hierarchy(m_tree_l_max)%is_hashed) then
            !$OMP PARALLEL DO PRIVATE(i, hierarchy_type, C, ignore_box) &
            !$OMP  FIRSTPRIVATE(uindex)
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((uindex%n .ge. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                    C%uindex = uindex
                    if (.not. ignore_box) then
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(i) = C
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end if
            end do
            !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, C, ignore_box) &
            !$OMP  FIRSTPRIVATE(uindex)
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index(i)
                if (list_index .gt. 0) then
                    hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)

                    if ((hierarchy_type .eq. HIERARCHY_X) .or. &
                        (hierarchy_type .eq. HIERARCHY_XY)) then

                        C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                        C%uindex = uindex
								!print*, 'n of uindex in step 1', uindex%n
                        if (.not. ignore_box) then
                            m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(list_index) = C
                            m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                        end if
                    end if
                end if
            end do
            !$OMP END PARALLEL DO
        end if

    end subroutine lib_ml_fmm_calculate_upward_pass_step_1

    !
    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       universal index of the a Tree-box
    !
    ! Returns
    ! ----
    !   C_i: type(lib_ml_fmm_coefficient)
    !       C coefficient of the Tree-box
    function lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box) result(C_i)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        logical, intent(out) :: ignore_box
        type(lib_ml_fmm_coefficient) :: C_i

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: buffer_C_i

        type(lib_tree_data_element), dimension(:), allocatable :: data_element
        integer(kind=4), dimension(:), allocatable :: element_number
        integer(kind=4) :: m_element_number

        type(lib_tree_spatial_point) :: x_c

        integer(kind=UINDEX_BYTES) :: i

        call lib_ml_fmm_type_operator_set_coefficient_zero(C_i)

        data_element = lib_tree_get_domain_e1(uindex, element_number)
        if ((allocated (data_element)) &
            .and. (size(data_element) .gt. 0)) then
            ignore_box = .false.
            x_c = lib_tree_get_centre_of_box(uindex)
            if (m_use_own_sum) then
                C_i = m_ml_fmm_handles%get_c(x_c, data_element, element_number)
            else
                do i=1, size(data_element)

                    if (data_element(i)%hierarchy .eq. HIERARCHY_X) then
                        m_element_number = element_number(i)
                    else if (data_element(i)%hierarchy .eq. HIERARCHY_XY) then
                        m_element_number = element_number(i) &
                                           - int(sum(m_data_element_hierarchy_type_length(HIERARCHY_X:HIERARCHY_Y)), &
                                                 CORRESPONDENCE_VECTOR_KIND)
                    else
                        exit
    !                    print *, "lib_ml_fmm_get_C_i_from_elements_at_box: ERROR"
    !                    print *, "  Wrong hierarchy of data element", i, " :", data_element(i)%hierarchy
                    end if

                    buffer_C_i = m_ml_fmm_handles%get_u_B_i(x_c, data_element(i), m_element_number)
                    C_i = C_i + buffer_C_i
                end do
            end if
        else
            ignore_box = .true.
        end if
    end function lib_ml_fmm_get_C_i_from_elements_at_box

    ! Upward pass - step 2
    !
    ! Procedure
    ! ----
    !   for each box at l_max -1, ..., l_min
    !       1. "reexpanding v^(1)_Children(X;n,l),l+1 (y) near the center of box (n, l) and
    !          summing up the contribution of all the child boxes"
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(33)
    !
    subroutine lib_ml_fmm_calculate_upward_pass_step_2()
        implicit none
        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C
        type(lib_tree_universal_index) :: uindex

        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=1) :: l

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type


        do l=m_tree_l_max-int(1, 1), m_tree_l_min, -int(1, 1)
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                !$OMP PARALLEL DO PRIVATE(i, hierarchy_type, C) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            C = lib_ml_fmm_get_C_of_box(uindex)
                            C%uindex = uindex
                            m_ml_fmm_hierarchy(l)%coefficient_list(i) = C
                            m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, C) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            C = lib_ml_fmm_get_C_of_box(uindex)
                            C%uindex = uindex
                            m_ml_fmm_hierarchy(l)%coefficient_list(list_index) = C
                            m_ml_fmm_hierarchy(l)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            end if
        end do

    end subroutine lib_ml_fmm_calculate_upward_pass_step_2

    ! Argument
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       universal index of a box
    !
    ! Returns
    ! ----
    !   C: type(lib_ml_fmm_coefficient)
    !       expansion coefficients of the box uindex
    function lib_ml_fmm_get_C_of_box(uindex) result(C)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_ml_fmm_coefficient) :: C

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C_child
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_children
        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_child

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: buffer_C

        x_c = lib_tree_get_centre_of_box(uindex)

        ! get child boxes
        uindex_children = lib_tree_get_children(uindex)

        allocate(buffer_C(size(uindex_children)))

        !$OMP PARALLEL DO PRIVATE(i, x_c_child, list_index, hierarchy_type)
        do i=1, size(uindex_children)
            x_c_child = lib_tree_get_centre_of_box(uindex_children(i))
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex_children(i)%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C_child = m_ml_fmm_hierarchy(uindex_children(i)%l)%coefficient_list(list_index)
                    buffer_C(i) = m_ml_fmm_handles%get_translation_SS(C_child, x_c_child, x_c)
                end if
            end if
        end do
        !$OMP END PARALLEL DO

        call lib_ml_fmm_type_operator_set_coefficient_zero(C)

        do i=1, size(uindex_children)
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex_children(i)%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = C + buffer_C(i)
                end if
            end if
        end do

    end function lib_ml_fmm_get_C_of_box

    subroutine lib_ml_fmm_calculate_downward_pass
        implicit none
#ifdef _TRACE_
        ! CPU-time
        real :: test_start_sub, test_finish_sub
        ! WALL-time
        INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

        call system_clock(test_count_start_sub, test_count_rate_sub)
        call cpu_time(test_start_sub)
#endif

        call lib_ml_fmm_calculate_downward_pass_step_1()
#ifdef _TRACE_
        call cpu_time(test_finish_sub)
        call system_clock(test_count_finish_sub, test_count_rate_sub)

        print *, "TRACE: lib_ml_fmm_calculate_downward_pass_step_1 ..done"
        print '("    CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
        print '("    WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)

        call system_clock(test_count_start_sub, test_count_rate_sub)
        call cpu_time(test_start_sub)

#endif

        call lib_ml_fmm_calculate_downward_pass_step_2()
#ifdef _TRACE_
        call cpu_time(test_finish_sub)
        call system_clock(test_count_finish_sub, test_count_rate_sub)

        print *, "TRACE: lib_ml_fmm_calculate_downward_pass_step_2 ..done"
        print '("    CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
        print '("    WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)
#endif

    end subroutine

    ! Downward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_min, ..., l_max
    !       - SR-transformation of the C coefficients into the middle of box I_4
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(34)
    !
    subroutine lib_ml_fmm_calculate_downward_pass_step_1()
        implicit none
        ! auxilary
        type(lib_tree_universal_index) :: uindex
        type(lib_ml_fmm_coefficient) :: D_tilde

        integer(kind=1) :: l
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: coefficient_list
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: list_index_list

        call lib_ml_fmm_type_operator_set_coefficient_zero(D_tilde)
        do l=m_tree_l_min, m_tree_l_max
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            ! create coefficient buffer
            allocate(coefficient_list(number_of_boxes))
            allocate(list_index_list(number_of_boxes))
            list_index_list(:) = -1

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                !$OMP PARALLEL DO PRIVATE(i, hierarchy_type, D_tilde) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                            D_tilde%uindex = uindex
                            list_index_list(i) = i
                            coefficient_list(i) = D_tilde
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, D_tilde) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                            D_tilde%uindex = uindex
                            list_index_list(i) = list_index
                            coefficient_list(i) = D_tilde
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            end if

            ! wirte to hierarchy
!            where (list_index_list .gt. -1)
!                m_ml_fmm_hierarchy(l)%coefficient_list(list_index_list) = coefficient_list
!                m_ml_fmm_hierarchy(l)%coefficient_type(list_index_list) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
!            end where
            !$OMP PARALLEL DO PRIVATE(i)
            do i=1, number_of_boxes
                if (list_index_list(i) .gt. -1) then
                    m_ml_fmm_hierarchy(l)%coefficient_list(list_index_list(i)) = coefficient_list(i)
                    m_ml_fmm_hierarchy(l)%coefficient_type(list_index_list(i)) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
                end if
            end do
            !$OMP END PARALLEL DO


            ! clean up
            if (allocated(coefficient_list)) then
                deallocate(coefficient_list)
            end if
            if (allocated(list_index_list)) then
                deallocate(list_index_list)
            end if
        end do


    end subroutine lib_ml_fmm_calculate_downward_pass_step_1

    ! Translates the C coefficient of a box from the regular to the local expansion.
    !
    ! Argument
    ! ----
    !   uindex: lib_tree_universal_index
    !       index of the box
    !
    ! Returns
    ! ----
    !   D_tilde: lib_ml_fmm_coefficient
    !       C coefficient of the box *uindex* S|R-translated
    !
    function lib_ml_fmm_get_D_tilde_of_box(uindex) result(D_tilde)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex
        type(lib_ml_fmm_coefficient) :: D_tilde

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C
        type(lib_tree_universal_index), dimension(:), allocatable :: boxes_i4

        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_neighbour

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: buffer_D_tilde

        allocate(boxes_i4, source=lib_tree_get_domain_i4(uindex, m_tree_neighbourhood_size_k))


        x_c = lib_tree_get_centre_of_box(uindex)

        allocate(buffer_D_tilde(size(boxes_i4)))

        !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, C, x_c_neighbour)
        do i=1, size(boxes_i4)
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, boxes_i4(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = m_ml_fmm_hierarchy(uindex%l)%coefficient_list(list_index)
                    x_c_neighbour = lib_tree_get_centre_of_box(boxes_i4(i))

                    buffer_D_tilde(i) = m_ml_fmm_handles%get_translation_SR(C, x_c, x_c_neighbour)
                end if
            end if
        end do
        !$OMP END PARALLEL DO

        call lib_ml_fmm_type_operator_set_coefficient_zero(D_tilde)

        do i=1, size(boxes_i4)
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, boxes_i4(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    D_tilde = D_tilde + buffer_D_tilde(i)
                end if
            end if
        end do

    end function lib_ml_fmm_get_D_tilde_of_box

    ! Downward pass - step 2
    !
    ! Procedure
    ! ----
    !   for boxes at l_min
    !       - D = D_tilde
    !
    !   for boxes at l_min + 1, ..., l_max
    !       - adds the local expansion of a box with the RR translated local expansion of the parent box
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(36)
    !
    subroutine lib_ml_fmm_calculate_downward_pass_step_2()
        implicit none
        ! auxiliaray
        type(lib_ml_fmm_coefficient) :: D
        type(lib_tree_universal_index) :: uindex
        integer(kind=1) :: l
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

!        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: coefficient_list
!        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: list_index_list

        m_ml_fmm_hierarchy(m_tree_l_min)%coefficient_type(:) = LIB_ML_FMM_COEFFICIENT_TYPE_D

        do l=m_tree_l_min+int(1,1), m_tree_l_max
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            ! create coefficient buffer
!            allocate(coefficient_list(number_of_boxes))
!            allocate(list_index_list(number_of_boxes))
!            list_index_list(:) = -1

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                !$OMP PARALLEL DO PRIVATE(i, hierarchy_type, D) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D = lib_ml_fmm_get_D_of_box(uindex)
                            D%uindex = uindex
                            m_ml_fmm_hierarchy(l)%coefficient_list(i) = D
                            m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, D) &
                !$OMP  FIRSTPRIVATE(uindex)
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D = lib_ml_fmm_get_D_of_box(uindex)
                            D%uindex = uindex
                            m_ml_fmm_hierarchy(l)%coefficient_list(list_index) = D
                            m_ml_fmm_hierarchy(l)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                        end if
                    end if
                end do
                !$OMP END PARALLEL DO
            end if

            ! wirte to hierarchy
!            !$OMP PARALLEL DO PRIVATE(i)
!            do i=1, number_of_boxes
!                if (list_index_list(i) .gt. -1) then
!                    m_ml_fmm_hierarchy(l)%coefficient_list(list_index_list(i)) = coefficient_list(i)
!                    m_ml_fmm_hierarchy(l)%coefficient_type(list_index_list(i)) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
!                end if
!            end do
!            !$OMP END PARALLEL DO

            ! clean up
!            if (allocated(coefficient_list)) then
!                deallocate(coefficient_list)
!            end if
!            if (allocated(list_index_list)) then
!                deallocate(list_index_list)
!            end if
        end do

    end subroutine lib_ml_fmm_calculate_downward_pass_step_2

    ! Argument
    ! ----
    !   uindex: lib_tree_universal_index
    !       index of the box
    !
    ! Returns
    ! -----
    !   D: lib_ml_fmm_coefficient
    !       summation of the local exapnsion of the box and of the parent box
    !
    function lib_ml_fmm_get_D_of_box(uindex) result(D)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex
        type(lib_ml_fmm_coefficient) :: D

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: D_parent
        type(lib_ml_fmm_coefficient) :: D_tilde
        type(lib_tree_universal_index) :: uindex_parent
        integer(kind=1) :: coefficient_type

        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_parent

        uindex_parent = lib_tree_get_parent(uindex)

        if (uindex_parent%n .ge. 0) then
            x_c = lib_tree_get_centre_of_box(uindex)
            x_c_parent = lib_tree_get_centre_of_box(uindex_parent)
            D_tilde = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
            D_parent = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex_parent, coefficient_type)
            D = D_tilde + m_ml_fmm_handles%get_translation_RR(D_parent, x_c_parent, x_c)
        else
            print *, "lib_ml_fmm_get_D_of_box: ERROR"
            print *, "  parent box not found"
        end if

    end function lib_ml_fmm_get_D_of_box

    ! Final Summation
    !
    ! Procedure
    ! ----
    !   for boxes of the Y-hierarchy at l_max
    !       1. calculate DoR
    !       2. calculate the contribution of the source elements (X-hierarchy)
    !          at the E2 domain of the evaluation element (Y-hierarchy)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(38)
    subroutine lib_ml_fmm_final_summation(calculate_y_hierarchy, calculate_xy_hierarchy)
        implicit none
        logical, intent(in), optional :: calculate_y_hierarchy
        logical, intent(in), optional :: calculate_xy_hierarchy

        ! auxiliary
        type(lib_tree_universal_index) :: uindex
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        logical :: m_calculate_y_hierarchy
        logical :: m_calculate_xy_hierarchy
#ifdef _TRACE_
        ! CPU-time
        real :: test_start_sub, test_finish_sub
        ! WALL-time
        INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub

        call system_clock(test_count_start_sub, test_count_rate_sub)
        call cpu_time(test_start_sub)
#endif

        m_calculate_y_hierarchy = .true.
        if (present(calculate_y_hierarchy)) m_calculate_y_hierarchy = calculate_y_hierarchy

        m_calculate_xy_hierarchy = .true.
        if (present(calculate_xy_hierarchy)) m_calculate_xy_hierarchy = calculate_xy_hierarchy

        uindex%l = m_tree_l_max
        number_of_boxes = size(m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index)

        if (m_ml_fmm_hierarchy(uindex%l)%is_hashed) then
            !$OMP PARALLEL DO PRIVATE(i, hierarchy_type) &
            !$OMP  FIRSTPRIVATE(uindex)
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                if (uindex%n .ge. 0) then
                    hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                    if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then

                        call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex, &
                                    m_calculate_y_hierarchy, m_calculate_xy_hierarchy)
                    end if
                end if
            end do
            !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type) &
            !$OMP  FIRSTPRIVATE(uindex)
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                if (list_index .gt. 0) then
                    hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
                    if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then

                        call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex, &
                                    m_calculate_y_hierarchy, m_calculate_xy_hierarchy)
                    end if
                end if
            end do
            !$OMP END PARALLEL DO
        end if
#ifdef _TRACE_
        call cpu_time(test_finish_sub)
        call system_clock(test_count_finish_sub, test_count_rate_sub)

        print *, "TRACE: lib_ml_fmm_final_summation ..done"
        print '("    CPU-Time = ",f10.3," seconds.")', test_finish_sub-test_start_sub
        print '("    WALL-Time = ",f10.3," seconds.")', (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)
#endif

    end subroutine lib_ml_fmm_final_summation

    ! Argument
    ! ----
    !   uindex: lib_tree_universal_index
    !       uiniversal index of a box of the Y-hierarchy
    !
    !Procedure
    ! ----
    !   for elements of the box *uindex*
    !       1. calculate DoR
    !       2. calculate the contribution of the source elements (X-hierarchy)
    !          at the E2 domain of the evaluation element (Y-hierarchy)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(38)
    subroutine lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex, calculate_y_hierarchy, calculate_xy_hierarchy)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex
        logical, intent(in), optional :: calculate_y_hierarchy
        logical, intent(in), optional :: calculate_xy_hierarchy

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: D
        integer(kind=1) :: coefficient_type
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e1
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e1
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2

        integer(kind=UINDEX_BYTES) :: v_counter

        logical :: m_calculate_y_hierarchy
        logical :: m_calculate_xy_hierarchy

        m_calculate_y_hierarchy = .true.
        if (present(calculate_y_hierarchy)) m_calculate_y_hierarchy = calculate_y_hierarchy

        m_calculate_xy_hierarchy = .true.
        if (present(calculate_xy_hierarchy)) m_calculate_xy_hierarchy = calculate_xy_hierarchy


        allocate(data_element_e1, source = lib_tree_get_domain_e1(uindex, &
                                                                  element_number_e1))
        allocate(data_element_e2, source = lib_tree_get_domain_e2(m_tree_neighbourhood_size_k, &
                                                                  uindex, &
                                                                  element_number_e2))

        D = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
        !$OMP PARALLEL DO PRIVATE(i, v_counter)
        do i=1, size(data_element_e1)
            if ((data_element_e1(i)%hierarchy .eq. HIERARCHY_Y .and. m_calculate_y_hierarchy) .or. &
                (data_element_e1(i)%hierarchy .eq. HIERARCHY_XY .and. m_calculate_xy_hierarchy)) then
                v_counter = element_number_e1(i)

                if (m_use_own_sum) then
                    m_ml_fmm_v(v_counter) = m_ml_fmm_handles%get_v_y_j(data_element_e1(i), element_number_e1(i), &
                                                                       data_element_e2, element_number_e2, D)
                else
                    m_ml_fmm_v(v_counter) = lib_ml_fmm_calculate_v_y_j(data_element_e1(i), element_number_e1(i), &
                                                                       data_element_e2, element_number_e2, D)
                end if
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine lib_ml_fmm_calculate_all_v_y_j_at_uindex

    ! Argument
    ! ----
    !   data_element_y_j: lib_tree_data_element
    !       evaluation element
    !   data_element_e2: lib_tree_data_element array
    !       elements of the E2 domain of the evaluation element *data_element_y_j*
    !   D: lib_ml_fmm_coefficient
    !       expansion coefficient of the box containing the element *data_element_y_j*
    !
    !Procedure
    ! ----
    !   for a evaluation element
    !       1. calculate DoR
    !       2. calculate the contribution of the source elements (X-hierarchy)
    !          at the E2 domain of the evaluation element (Y-hierarchy)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(38)
    !
    ! todo: re-calculate the element number: should be correspondence with the x-, y-, or, xy-hierarchy lsit
    function lib_ml_fmm_calculate_v_y_j(data_element_y_j, element_number_j, data_element_e2, element_number_e2, D) result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(inout) :: data_element_y_j
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_number_j
        type(lib_tree_data_element), dimension(:), allocatable, intent(inout) :: data_element_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2
        type(lib_ml_fmm_coefficient), intent(inout) :: D
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_spatial_point) :: x_c

        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: m_element_number_j
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: m_element_number_i

        x_c = lib_tree_get_centre_of_box(data_element_y_j%uindex)

        if (data_element_y_j%hierarchy .eq. HIERARCHY_Y) then
            m_element_number_j = element_number_j - int(m_data_element_hierarchy_type_length(HIERARCHY_X), &
                                                        CORRESPONDENCE_VECTOR_KIND)
        else if (data_element_y_j%hierarchy .eq. HIERARCHY_XY) then
            m_element_number_j = element_number_j - int(sum(m_data_element_hierarchy_type_length(HIERARCHY_X:HIERARCHY_Y)), &
                                                        CORRESPONDENCE_VECTOR_KIND)
        end if

        rv = m_ml_fmm_handles%dor(D, x_c, data_element_y_j, m_element_number_j)

        do i=1, size(data_element_e2)
            if (data_element_e2(i)%hierarchy .eq. HIERARCHY_X &
                .or. data_element_e2(i)%hierarchy .eq. HIERARCHY_XY) then

                if (data_element_e2(i)%hierarchy .eq. HIERARCHY_X) then
                    m_element_number_i = element_number_e2(i)
                else if (data_element_e2(i)%hierarchy .eq. HIERARCHY_XY) then
                    m_element_number_i = element_number_e2(i) &
                                         - int(sum(m_data_element_hierarchy_type_length(HIERARCHY_X:HIERARCHY_Y)), &
                                               CORRESPONDENCE_VECTOR_KIND)
                end if

                rv = rv + m_ml_fmm_handles%get_u_phi_i_j(data_element_e2(i), m_element_number_i, &
                                                         data_element_y_j, m_element_number_j)
            end if
        end do

    end function lib_ml_fmm_calculate_v_y_j

    ! Returns
    ! ----
    !   rv: integer
    !       neighbourhood size
    function lib_ml_fmm_get_neighbourhood_size_k() result(rv)
        implicit none
        ! dummy
        integer :: rv

        rv = m_tree_neighbourhood_size_k

    end function lib_ml_fmm_get_neighbourhood_size_k

    ! ----- test functions -----
    function lib_ml_fmm_test_functions() result(error_counter)
        implicit none
        ! dummy
        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_ml_fmm_constructor()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_upward_pass()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_downward_pass_1()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_downward_pass_2()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_final_summation()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_upward_pass_v2()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_downward_pass_1_v2()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_calculate_downward_pass_2_v2()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        if (.not. test_lib_ml_fmm_final_summation_v2()) error_counter = error_counter + 1
        call lib_ml_fmm_destructor()

        print *, "-------------lib_ml_fmm_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_ml_fmm_test_functions tests: OK"
        else
            print *, error_counter,"lib_ml_fmm_test_functions test(s) FAILED"
        end if
        print *, "------------------------------------------------------"

        contains

        ! Test values
        !
        !  Application
        !  ----
        !
        !  level 2:
        !      C: 1    6   -   28  -   -   97  54  124
        !      n: 0    1   2   3   8   10  12  13  15
        !   type: X    X   Y   X   Y   Y   X   X   X
        !
        !
        !  level 3:
        !      C: 0  13  49  6   63  48  61  1   54  15  -   -   -   -
        !      n: 0  13  49  6   63  48  61  1   54  15  42  10  34  2
        !   type: <--             X                 -->  <--   Y   -->
        !
        !
        !
        !  Manual caluclation
        !  ----
        !
        !  Upward pass
        !  -----
        !
        !  level 2
        !          C:     1       6   -     28    -   -      97   54     124
        !          n:     0       1   2     3     8   10     12   13     15
        !       type:     XY      X   Y     X     Y   Y      X    X      X
        !             ----|----   |   |   --|---  |   |   ---|--  |    --|---
        !  level 3:
        !          C: 0   1   -   6   -   13  15  -   -   48  49  54   61  63
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   X   Y   X   Y   <-X ->  <-Y -> <--     X        -->
        !
        !
        !  Downward pass
        !  -----
        !  level 2
        !          C:     1       6   -     28    -   -      97   54     124
        !         D~:  97+54+124  -   275   -  275+7 275+35  -    -      -
        !          D:     275     -   275   -     282 310    -    -      -    <-- ! D=D~ ! reason: l = l_min
        !          n:     0       1   2     3     8   10     12   13     15
        !       type:     XY      X   Y     X     Y   Y      X    X      X
        !             ----|----   |   |   --|---  |   |   ---|--  |    --|---
        !  level 3:
        !         D~: -   - 6+28  -  7+28 -   -   28  0   -   -   -    -   -
        !          D: -   - 275+34- 275+35-   -282+28 310 -   -   -    -   -
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   X   Y   X   Y   <-X ->  <-Y -> <--     X        -->
        !
        !  Final Summation:
        !        2*D: -   -  618  -  620  -   -  620  620 -   -   -    -   -
        !sum(phi_e2): -   -   1   -   0   -   -   0   0   -   -   -    -   -
        !      final: -   -  619  -  620  -   -  620  620 -   -   -    -   -
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   X   Y   X   Y   <-X ->  <-Y -> <--     X        -->
        !
        !  --------------------------------- ---------------------------------
        !  |  21   |  23   | 29    | 31    | |  53   |  55   | 61    | 63    |
        !  |       |       |       |       | |       |       |    X  |    X  |
        !  |-------5---------------7-------| |-------13--------------15------|
        !  |  20   |  22   | 28    | 30    | |  52   |  54   | 60    | 62    |
        !  |       |       |       |       | |       |    X  |       |       |
        !  |---------------1---------------| |---------------3---------------|
        !  |  17   |  19   | 25    | 27    | |  49   |  51   | 57    | 59    |
        !  |       |       |       |       | |    X  |       |       |       |
        !  |-------4---------------6-------| |-------12--------------14------|
        !  |  16   |  18   | 24    | 26    | |  48   |  50   | 56    | 58    |
        !  |       |       |       |       | |    X  |       |       |       |
        !  --------------------------------- ---------------------------------
        !  --------------------------------- ---------------------------------
        !  |  5    |  7    | 13    | 15    | |  37   |  39   | 45    | 47    |
        !  |       |       |    X  |    X  | |       |       |       |       |
        !  |-------1---------------3-------| |-------9---------------11------|
        !  |  4    |  6    | 12    | 14    | |  36   |  38   | 44    | 46    |
        !  |       |    X  |       |       | |       |       |       |       |
        !  |---------------0---------------| |---------------2---------------|
        !  |  1    |  3    | 9     | 11    | |  33   |  35   | 41    | 43    |
        !  |    X  |       |       |       | |       |       |       |       |
        !  |-------0---------------2-------| |-------8---------------10------|
        !  |  0    |  2    | 8     | 10    | |  32   |  34   | 40    | 42    |
        !  |    X  |    Y  |       |     Y | |       |     Y |       |    Y  |
        !  --------------------------------- ---------------------------------
        subroutine setup_hierarchy_with_test_data_2D()
            implicit none

            ! auxiliary
            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            type(ml_fmm_type_operator_procedures) :: operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures
            type(lib_ml_fmm_data) :: data_elements

             ! --- generate test data ---
            allocate(data_elements%X, source=lib_tree_get_diagonal_test_dataset(list_length, element_type, HIERARCHY_X))
!            allocate(data_elements%Y, source=lib_tree_get_diagonal_test_dataset(4_8, element_type, HIERARCHY_Y, .true.))

            allocate(data_elements%Y(4))
            data_elements%Y(:)%hierarchy = HIERARCHY_Y
            data_elements%Y(:)%element_type = element_type

            data_elements%Y(1)%point_x%x(1) = 0.24975000321865082D0
            data_elements%Y(1)%point_x%x(2) = 0.00024999678134918213D0
            data_elements%Y(2)%point_x%x(1) = 0.49950000643730164D0
            data_elements%Y(2)%point_x%x(2) = 0.00049999356269836426D0
            data_elements%Y(3)%point_x%x(1) = 0.74924999475479126D0
            data_elements%Y(3)%point_x%x(2) = 0.00074999034404754639D0
            data_elements%Y(4)%point_x%x(1) = 0.99900001287460327D0
            data_elements%Y(4)%point_x%x(2) = 0.00099998712539672852D0


            operator_procedures = ml_fmm_type_operator_get_procedures()
            ml_fmm_procedures = ml_fmm_get_test_procedures()
            call lib_ml_fmm_constructor(data_elements, operator_procedures, ml_fmm_procedures)

        end subroutine setup_hierarchy_with_test_data_2D

        subroutine setup_ground_truth_uindex_list_2D(ground_truth_uindex_list_X_l_2, &
                                                     ground_truth_uindex_list_Y_l_2, &
                                                     ground_truth_uindex_list_XY_l_2, &
                                                     ground_truth_uindex_list_X_l_3, &
                                                     ground_truth_uindex_list_Y_l_3, &
                                                     ground_truth_uindex_list_XY_l_3)
            implicit none
            ! dummy
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_X_l_2
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_Y_l_2
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_XY_l_2

            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_X_l_3
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_Y_l_3
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_XY_l_3

            ! --- setup ground truth data (2D) ---
            allocate(ground_truth_uindex_list_X_l_2(5))
            allocate(ground_truth_uindex_list_Y_l_2(3))
            allocate(ground_truth_uindex_list_XY_l_2(1))

            ground_truth_uindex_list_X_l_2(:)%l = 2
            ground_truth_uindex_list_X_l_2(1)%n = 1
            ground_truth_uindex_list_X_l_2(2)%n = 3
            ground_truth_uindex_list_X_l_2(3)%n = 12
            ground_truth_uindex_list_X_l_2(4)%n = 13
            ground_truth_uindex_list_X_l_2(5)%n = 15

            ground_truth_uindex_list_Y_l_2(:)%l = 2
            ground_truth_uindex_list_Y_l_2(1)%n = 2
            ground_truth_uindex_list_Y_l_2(2)%n = 8
            ground_truth_uindex_list_Y_l_2(3)%n = 10

            ground_truth_uindex_list_XY_l_2(:)%l = 2
            ground_truth_uindex_list_XY_l_2(1)%n = 0

            allocate(ground_truth_uindex_list_X_l_3(10))
            allocate(ground_truth_uindex_list_Y_l_3(4))
            allocate(ground_truth_uindex_list_XY_l_3(0))

            ground_truth_uindex_list_X_l_3(:)%l = 3
            ground_truth_uindex_list_X_l_3(1)%n = 0
            ground_truth_uindex_list_X_l_3(2)%n = 13
            ground_truth_uindex_list_X_l_3(3)%n = 49
            ground_truth_uindex_list_X_l_3(4)%n = 6
            ground_truth_uindex_list_X_l_3(5)%n = 63
            ground_truth_uindex_list_X_l_3(6)%n = 48
            ground_truth_uindex_list_X_l_3(7)%n = 61
            ground_truth_uindex_list_X_l_3(8)%n = 1
            ground_truth_uindex_list_X_l_3(9)%n = 54
            ground_truth_uindex_list_X_l_3(10)%n = 15

            ground_truth_uindex_list_Y_l_3(:)%l = 3
            ground_truth_uindex_list_Y_l_3(1)%n = 42
            ground_truth_uindex_list_Y_l_3(2)%n = 10
            ground_truth_uindex_list_Y_l_3(3)%n = 34
            ground_truth_uindex_list_Y_l_3(4)%n = 2

            ground_truth_uindex_list_XY_l_3(:)%l = 3

        end subroutine setup_ground_truth_uindex_list_2D

        subroutine setup_ground_truth_C_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                           ground_truth_uindex_list_l_2, &
                                                           ground_truth_coefficient_list_l_3, &
                                                           ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(10))
            allocate(ground_truth_uindex_list_l_3(10))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/0.0D0/))
            allocate(ground_truth_coefficient_list_l_3(1)%c, source = (/dcmplx(0.0D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1)%item, source = (/dcmplx(0, 0D0), dcmplx(0, 0D0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1)%item, source = (/dcmplx(0D0,0), dcmplx(0D0, 0)/))
            ground_truth_uindex_list_l_3(1)%n = 0

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/13D0/))
            allocate(ground_truth_coefficient_list_l_3(2)%c, source = (/dcmplx(13D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1)%item, source = (/dcmplx(0, 13D0), dcmplx(0, 13D0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1)%item, source = (/dcmplx(13D0,0), dcmplx(13D0, 0)/))
            ground_truth_uindex_list_l_3(2)%n = 13

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/49D0/))
            allocate(ground_truth_coefficient_list_l_3(3)%c, source = (/dcmplx(49D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1)%item, source = (/dcmplx(0, 49D0), dcmplx(0, 49D0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1)%item, source = (/dcmplx(49D0,0), dcmplx(49D0, 0)/))
            ground_truth_uindex_list_l_3(3)%n = 49

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/6D0/))
            allocate(ground_truth_coefficient_list_l_3(4)%c, source = (/dcmplx(6D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1)%item, source =  (/dcmplx(0, 6D0), dcmplx(0, 6D0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1)%item, source =  (/dcmplx(6D0,0), dcmplx(6D0, 0)/))
            ground_truth_uindex_list_l_3(4)%n = 6

            allocate(ground_truth_coefficient_list_l_3(5)%r, source = (/63D0/))
            allocate(ground_truth_coefficient_list_l_3(5)%c, source = (/dcmplx(63D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(5)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(5)%a_nm%item(1)%item, source = (/dcmplx(0, 63D0), dcmplx(0, 63D0)/))
            allocate(ground_truth_coefficient_list_l_3(5)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(5)%b_nm%item(1)%item, source = (/dcmplx(63D0,0), dcmplx(63D0, 0)/))
            ground_truth_uindex_list_l_3(5)%n = 63

            allocate(ground_truth_coefficient_list_l_3(6)%r, source = (/48D0/))
            allocate(ground_truth_coefficient_list_l_3(6)%c, source = (/dcmplx(48D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(6)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(6)%a_nm%item(1)%item, source = (/dcmplx(0, 48D0), dcmplx(0, 48D0)/))
            allocate(ground_truth_coefficient_list_l_3(6)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(6)%b_nm%item(1)%item, source = (/dcmplx(48D0,0), dcmplx(48D0, 0)/))
            ground_truth_uindex_list_l_3(6)%n = 48

            allocate(ground_truth_coefficient_list_l_3(7)%r, source = (/61D0/))
            allocate(ground_truth_coefficient_list_l_3(7)%c, source = (/dcmplx(61D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(7)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(7)%a_nm%item(1)%item, source = (/dcmplx(0, 61D0), dcmplx(0, 61D0)/))
            allocate(ground_truth_coefficient_list_l_3(7)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(7)%b_nm%item(1)%item, source = (/dcmplx(61D0,0), dcmplx(61D0, 0)/))
            ground_truth_uindex_list_l_3(7)%n = 61

            allocate(ground_truth_coefficient_list_l_3(8)%r, source = (/1D0/))
            allocate(ground_truth_coefficient_list_l_3(8)%c, source = (/dcmplx(1D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(8)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(8)%a_nm%item(1)%item, source =  (/dcmplx(0, 1D0), dcmplx(0, 1D0)/))
            allocate(ground_truth_coefficient_list_l_3(8)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(8)%b_nm%item(1)%item, source =  (/dcmplx(1D0,0), dcmplx(1D0, 0)/))
            ground_truth_uindex_list_l_3(8)%n = 1

            allocate(ground_truth_coefficient_list_l_3(9)%r, source = (/54D0/))
            allocate(ground_truth_coefficient_list_l_3(9)%c, source = (/dcmplx(54D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(9)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(9)%a_nm%item(1)%item, source = (/dcmplx(0, 54D0), dcmplx(0, 54D0)/))
            allocate(ground_truth_coefficient_list_l_3(9)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(9)%b_nm%item(1)%item, source = (/dcmplx(54D0,0), dcmplx(54D0, 0)/))
            ground_truth_uindex_list_l_3(9)%n = 54

            allocate(ground_truth_coefficient_list_l_3(10)%r, source = (/15D0/))
            allocate(ground_truth_coefficient_list_l_3(10)%c, source = (/dcmplx(15D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(10)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(10)%a_nm%item(1)%item, source = (/dcmplx(0, 15D0), dcmplx(0, 15D0)/))
            allocate(ground_truth_coefficient_list_l_3(10)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(10)%b_nm%item(1)%item, source = (/dcmplx(15D0,0), dcmplx(15D0, 0)/))
            ground_truth_uindex_list_l_3(10)%n = 15


            allocate(ground_truth_coefficient_list_l_2(6))
            allocate(ground_truth_uindex_list_l_2(6))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/1D0/))
            allocate(ground_truth_coefficient_list_l_2(1)%c, source = (/dcmplx(1D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1)%item, source = (/dcmplx(0, 1D0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1)%item, source = (/dcmplx(1D0, 0)/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/6D0/))
            allocate(ground_truth_coefficient_list_l_2(2)%c, source = (/dcmplx(6D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1)%item, source = (/dcmplx(0, 6D0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1)%item, source = (/dcmplx(6D0, 0)/))
            ground_truth_uindex_list_l_2(2)%n = 1

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/28D0/))
            allocate(ground_truth_coefficient_list_l_2(3)%c, source = (/dcmplx(28D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1)%item, source = (/dcmplx(0, 28D0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1)%item, source = (/dcmplx(28D0, 0)/))
            ground_truth_uindex_list_l_2(3)%n = 3

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/97D0/))
            allocate(ground_truth_coefficient_list_l_2(4)%c, source = (/dcmplx(97D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1)%item, source = (/dcmplx(0, 97D0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1)%item, source = (/dcmplx(97D0, 0)/))
            ground_truth_uindex_list_l_2(4)%n = 12

            allocate(ground_truth_coefficient_list_l_2(5)%r, source = (/54D0/))
            allocate(ground_truth_coefficient_list_l_2(5)%c, source = (/dcmplx(54D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(5)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(5)%a_nm%item(1)%item, source = (/dcmplx(0, 54D0)/))
            allocate(ground_truth_coefficient_list_l_2(5)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(5)%b_nm%item(1)%item, source = (/dcmplx(54D0, 0)/))
            ground_truth_uindex_list_l_2(5)%n = 13

            allocate(ground_truth_coefficient_list_l_2(6)%r, source = (/124D0/))
            allocate(ground_truth_coefficient_list_l_2(6)%c, source = (/dcmplx(124D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(6)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(6)%a_nm%item(1)%item, source = (/dcmplx(0, 124D0)/))
            allocate(ground_truth_coefficient_list_l_2(6)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(6)%b_nm%item(1)%item, source = (/dcmplx(124D0, 0)/))
            ground_truth_uindex_list_l_2(6)%n = 15

        end subroutine setup_ground_truth_C_coefficient_list_2D

        subroutine setup_ground_truth_D_tilde_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                                  ground_truth_uindex_list_l_2, &
                                                                  ground_truth_coefficient_list_l_3, &
                                                                  ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(4))
            allocate(ground_truth_uindex_list_l_3(4))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/34D0/))
            allocate(ground_truth_coefficient_list_l_3(1)%c, source = (/dcmplx(34D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1)%item, source = (/dcmplx(0, 34D0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1)%item, source = (/dcmplx(34D0, 0)/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/35D0/))
            allocate(ground_truth_coefficient_list_l_3(2)%c, source = (/dcmplx(35D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1)%item, source = (/dcmplx(0, 35D0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1)%item, source = (/dcmplx(35D0, 0)/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/28D0/))
            allocate(ground_truth_coefficient_list_l_3(3)%c, source = (/dcmplx(28D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1)%item, source = (/dcmplx(0, 28D0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1)%item, source = (/dcmplx(28D0, 0)/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/0D0/))
            allocate(ground_truth_coefficient_list_l_3(4)%c, source = (/dcmplx(0D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1)%item, source = (/dcmplx(0, 0D0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1)%item, source = (/dcmplx(0D0, 0)/))
            ground_truth_uindex_list_l_3(4)%n = 42


            allocate(ground_truth_coefficient_list_l_2(4))
            allocate(ground_truth_uindex_list_l_2(4))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/275D0/))
            allocate(ground_truth_coefficient_list_l_2(1)%c, source = (/dcmplx(275D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1)%item, source = (/dcmplx(0, 275D0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1)%item, source = (/dcmplx(275D0, 0)/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/275D0/))
            allocate(ground_truth_coefficient_list_l_2(2)%c, source = (/dcmplx(275D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1)%item, source = (/dcmplx(0, 275D0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1)%item, source = (/dcmplx(275D0, 0)/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/282D0/))
            allocate(ground_truth_coefficient_list_l_2(3)%c, source = (/dcmplx(282D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1)%item, source = (/dcmplx(0, 282D0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1)%item, source = (/dcmplx(282D0, 0)/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/310D0/))
            allocate(ground_truth_coefficient_list_l_2(4)%c, source = (/dcmplx(310D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1)%item, source = (/dcmplx(0, 310D0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1)%item, source = (/dcmplx(310D0, 0)/))
            ground_truth_uindex_list_l_2(4)%n = 10

        end subroutine setup_ground_truth_D_tilde_coefficient_list_2D

        subroutine setup_ground_truth_D_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                            ground_truth_uindex_list_l_2, &
                                                            ground_truth_coefficient_list_l_3, &
                                                            ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(4))
            allocate(ground_truth_uindex_list_l_3(4))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/309D0/))
            allocate(ground_truth_coefficient_list_l_3(1)%c, source = (/dcmplx(309D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%a_nm%item(1)%item, source = (/dcmplx(0, 309D0)/))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(1)%b_nm%item(1)%item, source = (/dcmplx(309D0, 0)/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/310D0/))
            allocate(ground_truth_coefficient_list_l_3(2)%c, source = (/dcmplx(310D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%a_nm%item(1)%item, source = (/dcmplx(0, 310D0)/))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(2)%b_nm%item(1)%item, source = (/dcmplx(310D0, 0)/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/310D0/))
            allocate(ground_truth_coefficient_list_l_3(3)%c, source = (/dcmplx(310D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%a_nm%item(1)%item, source = (/dcmplx(0, 310D0)/))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(3)%b_nm%item(1)%item, source = (/dcmplx(310D0, 0)/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/310D0/))
            allocate(ground_truth_coefficient_list_l_3(4)%c, source = (/dcmplx(310D0, 0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%a_nm%item(1)%item, source = (/dcmplx(0, 310D0)/))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_3(4)%b_nm%item(1)%item, source = (/dcmplx(310D0, 0)/))
            ground_truth_uindex_list_l_3(4)%n = 42


            allocate(ground_truth_coefficient_list_l_2(4))
            allocate(ground_truth_uindex_list_l_2(4))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/275D0/))
            allocate(ground_truth_coefficient_list_l_2(1)%c, source = (/dcmplx(275D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%a_nm%item(1)%item, source = (/dcmplx(0, 275D0)/))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(1)%b_nm%item(1)%item, source = (/dcmplx(275D0, 0)/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/275D0/))
            allocate(ground_truth_coefficient_list_l_2(2)%c, source = (/dcmplx(275D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%a_nm%item(1)%item, source = (/dcmplx(0, 275D0)/))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(2)%b_nm%item(1)%item, source = (/dcmplx(275D0, 0)/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/282D0/))
            allocate(ground_truth_coefficient_list_l_2(3)%c, source = (/dcmplx(282D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%a_nm%item(1)%item, source = (/dcmplx(0, 282D0)/))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(3)%b_nm%item(1)%item, source = (/dcmplx(282D0, 0)/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/310D0/))
            allocate(ground_truth_coefficient_list_l_2(4)%c, source = (/dcmplx(310D0, 0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%a_nm%item(1)%item, source = (/dcmplx(0, 310D0)/))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1))
            allocate(ground_truth_coefficient_list_l_2(4)%b_nm%item(1)%item, source = (/dcmplx(310D0, 0)/))
            ground_truth_uindex_list_l_2(4)%n = 10

        end subroutine setup_ground_truth_D_coefficient_list_2D

        subroutine setup_ground_truth_final_summation_2D(vector_v)
            implicit none
            ! dummy
            type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_v

            ! auxiliary
            integer :: i

            allocate(vector_v(14))
            do i=11, 14
                allocate(vector_v(i)%r(1))
                allocate(vector_v(i)%c(1))
                allocate(vector_v(i)%a_nm%item(1))
                allocate(vector_v(i)%a_nm%item(1)%item(1))
                allocate(vector_v(i)%b_nm%item(1))
                allocate(vector_v(i)%b_nm%item(1)%item(1))
            end do

            vector_v(11)%r(1) = 619
            vector_v(11)%c(1) = dcmplx(619, 0)
            vector_v(11)%a_nm%item(1)%item(1) = dcmplx(0, 619)
            vector_v(11)%b_nm%item(1)%item(1) = dcmplx(619, 0)
            vector_v(12)%r(1) = 620
            vector_v(12)%c(1) = dcmplx(620, 0)
            vector_v(12)%a_nm%item(1)%item(1) = dcmplx(0, 620)
            vector_v(12)%b_nm%item(1)%item(1) = dcmplx(620, 0)
            vector_v(13)%r(1) = 620
            vector_v(13)%c(1) = dcmplx(620, 0)
            vector_v(13)%a_nm%item(1)%item(1) = dcmplx(0, 620)
            vector_v(13)%b_nm%item(1)%item(1) = dcmplx(620, 0)
            vector_v(14)%r(1) = 620
            vector_v(14)%c(1) = dcmplx(620, 0)
            vector_v(14)%a_nm%item(1)%item(1) = dcmplx(0, 620)
            vector_v(14)%b_nm%item(1)%item(1) = dcmplx(620, 0)


        end subroutine setup_ground_truth_final_summation_2D

        ! Test values_2
        !
        !  Application
        !  ----
        !
        !  level 2:
        !      C: 1    6   -   28  -   42  97  54  124
        !      n: 0    1   2   3   8   10  12  13  15
        !   type: XY   X   Y   X   Y   XY  X   XY  X
        !
        !
        !  level 3:
        !      C: 0  13  49  6   63  48  61  1   54  15  42  -   -   -
        !      n: 0  13  49  6   63  48  61  1   54  15  42  10  34  2
        !   type: <--        X          -->  XY  XY  X   XY  <-- Y -->
        !
        !
        !
        !  Manual caluclation
        !  ----
        !
        !  Upward pass
        !  -----
        !
        !  level 2
        !          C:     1       6   -     28    -   42     97   54     124
        !          n:     0       1   2     3     8   10     12   13     15
        !       type:     XY      X   Y     X     Y   XY     X    XY     X
        !             ----|----   |   |   --|---  |   |   ---|--  |    --|---
        !  level 3:
        !          C: 0   1   -   6   -   13  15  -   42  48  49  54   61  63
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   XY  Y   X   Y   <-X ->  Y   XY  < X  >  XY   X   X
        !
        !
        !  Downward pass
        !  -----
        !  level 2
        !          C:     1       6   -     28    -   42     97   54     124
        !         D~:97+54+124+42 -   317   -  317+7 317+35  - 1+6+28+42 -
        !          D:     317     -   317   -     324 352    -    77     -    <-- ! D=D~ ! reason: l = l_min
        !          n:     0       1   2     3     8   10     12   13     15
        !       type:     XY      X   Y     X     Y   XY     X    XY     X
        !             ----|----   |   |   --|---  |   |   ---|--  |    --|---
        !  level 3:
        !         D~: -   28 6+28 -  7+28 -   - 28+42 0   -   - 224+48 -   -
        !          D: - 345 317+34- 317+35-   -324+70 352 -   - 272+77 -   -
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   XY  Y   X   Y   <-X ->  Y   XY  < X  >  XY   X   X
        !
        !  Final Summation:
        !        2*D: -  690 702  -  704  -   -  788  704 -   -   698  -   -
        !sum(phi_e2): -   6   1   -   0   -   -   0   0   -   -   110  -   -
        !      final: -  696 703  -  704  -   -  788  704 -   -   808  -   -
        !          n: 0   1   2   6   10  13  15  34  42  48  49  54   61  63
        !       type: X   XY  Y   X   Y   <-X ->  Y   XY  < X  >  XY   X   X
        !
        !  --------------------------------- ---------------------------------
        !  |  21   |  23   | 29    | 31    | |  53   |  55   | 61    | 63    |
        !  |       |       |       |       | |       |       |    X  |    X  |
        !  |-------5---------------7-------| |-------13--------------15------|
        !  |  20   |  22   | 28    | 30    | |  52   |  54   | 60    | 62    |
        !  |       |       |       |       | |       |    XY |       |       |
        !  |---------------1---------------| |---------------3---------------|
        !  |  17   |  19   | 25    | 27    | |  49   |  51   | 57    | 59    |
        !  |       |       |       |       | |    X  |       |       |       |
        !  |-------4---------------6-------| |-------12--------------14------|
        !  |  16   |  18   | 24    | 26    | |  48   |  50   | 56    | 58    |
        !  |       |       |       |       | |    X  |       |       |       |
        !  --------------------------------- ---------------------------------
        !  --------------------------------- ---------------------------------
        !  |  5    |  7    | 13    | 15    | |  37   |  39   | 45    | 47    |
        !  |       |       |    X  |    X  | |       |       |       |       |
        !  |-------1---------------3-------| |-------9---------------11------|
        !  |  4    |  6    | 12    | 14    | |  36   |  38   | 44    | 46    |
        !  |       |    X  |       |       | |       |       |       |       |
        !  |---------------0---------------| |---------------2---------------|
        !  |  1    |  3    | 9     | 11    | |  33   |  35   | 41    | 43    |
        !  |    XY |       |       |       | |       |       |       |       |
        !  |-------0---------------2-------| |-------8---------------10------|
        !  |  0    |  2    | 8     | 10    | |  32   |  34   | 40    | 42    |
        !  |    X  |    Y  |       |     Y | |       |     Y |       |    XY |
        !  --------------------------------- ---------------------------------
        subroutine setup_hierarchy_with_test_data_2D_2()
            implicit none

            ! auxiliary
            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            type(ml_fmm_type_operator_procedures) :: operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures
            type(lib_ml_fmm_data) :: data_elements

             ! --- generate test data ---
            allocate(data_elements%X, source=lib_tree_get_diagonal_test_dataset(list_length, element_type, HIERARCHY_X))
!            allocate(data_elements%Y, source=lib_tree_get_diagonal_test_dataset(4_8, element_type, HIERARCHY_Y, .true.))

            allocate(data_elements%Y(4))
            data_elements%Y(:)%hierarchy = HIERARCHY_Y
            data_elements%Y(:)%element_type = element_type

            data_elements%Y(1)%point_x%x(1) = 0.24975000321865082D0
            data_elements%Y(1)%point_x%x(2) = 0.00024999678134918213D0
            data_elements%Y(2)%point_x%x(1) = 0.49950000643730164D0
            data_elements%Y(2)%point_x%x(2) = 0.00049999356269836426D0
            data_elements%Y(3)%point_x%x(1) = 0.74924999475479126D0
            data_elements%Y(3)%point_x%x(2) = 0.00074999034404754639D0
            data_elements%Y(4)%point_x%x(1) = 0.99900001287460327D0
            data_elements%Y(4)%point_x%x(2) = 0.00099998712539672852D0


            operator_procedures = ml_fmm_type_operator_get_procedures()
            ml_fmm_procedures = ml_fmm_get_test_procedures()
            call lib_ml_fmm_constructor(data_elements, operator_procedures, ml_fmm_procedures)

        end subroutine setup_hierarchy_with_test_data_2D_2

        subroutine setup_ground_truth_uindex_list_2D_2(ground_truth_uindex_list_X_l_2, &
                                                       ground_truth_uindex_list_Y_l_2, &
                                                       ground_truth_uindex_list_XY_l_2, &
                                                       ground_truth_uindex_list_X_l_3, &
                                                       ground_truth_uindex_list_Y_l_3, &
                                                       ground_truth_uindex_list_XY_l_3)
            implicit none
            ! dummy
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_X_l_2
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_Y_l_2
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_XY_l_2

            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_X_l_3
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_Y_l_3
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: ground_truth_uindex_list_XY_l_3

            ! --- setup ground truth data (2D) ---
            allocate(ground_truth_uindex_list_X_l_2(4))
            allocate(ground_truth_uindex_list_Y_l_2(3))
            allocate(ground_truth_uindex_list_XY_l_2(3))

            ground_truth_uindex_list_X_l_2(:)%l = 2
            ground_truth_uindex_list_X_l_2(1)%n = 1
            ground_truth_uindex_list_X_l_2(2)%n = 3
            ground_truth_uindex_list_X_l_2(3)%n = 12
            ground_truth_uindex_list_X_l_2(4)%n = 15

            ground_truth_uindex_list_Y_l_2(:)%l = 2
            ground_truth_uindex_list_Y_l_2(1)%n = 2
            ground_truth_uindex_list_Y_l_2(2)%n = 8

            ground_truth_uindex_list_XY_l_2(:)%l = 2
            ground_truth_uindex_list_XY_l_2(1)%n = 0
            ground_truth_uindex_list_XY_l_2(2)%n = 10
            ground_truth_uindex_list_XY_l_2(3)%n = 13

            allocate(ground_truth_uindex_list_X_l_3(8))
            allocate(ground_truth_uindex_list_Y_l_3(3))
            allocate(ground_truth_uindex_list_XY_l_3(3))

            ground_truth_uindex_list_X_l_3(:)%l = 3
            ground_truth_uindex_list_X_l_3(1)%n = 0
            ground_truth_uindex_list_X_l_3(2)%n = 13
            ground_truth_uindex_list_X_l_3(3)%n = 49
            ground_truth_uindex_list_X_l_3(4)%n = 6
            ground_truth_uindex_list_X_l_3(5)%n = 63
            ground_truth_uindex_list_X_l_3(6)%n = 48
            ground_truth_uindex_list_X_l_3(7)%n = 61
            ground_truth_uindex_list_X_l_3(8)%n = 15

            ground_truth_uindex_list_Y_l_3(:)%l = 3
            ground_truth_uindex_list_Y_l_3(1)%n = 10
            ground_truth_uindex_list_Y_l_3(2)%n = 34
            ground_truth_uindex_list_Y_l_3(3)%n = 2

            ground_truth_uindex_list_XY_l_3(:)%l = 3
            ground_truth_uindex_list_XY_l_3(1)%n = 1
            ground_truth_uindex_list_XY_l_3(2)%n = 54
            ground_truth_uindex_list_XY_l_3(3)%n = 42

        end subroutine setup_ground_truth_uindex_list_2D_2

        subroutine setup_ground_truth_C_coefficient_list_2D_2(ground_truth_coefficient_list_l_2, &
                                                              ground_truth_uindex_list_l_2, &
                                                              ground_truth_coefficient_list_l_3, &
                                                              ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(11))
            allocate(ground_truth_uindex_list_l_3(11))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/0.0D0/))
            ground_truth_uindex_list_l_3(1)%n = 0

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/13D0/))
            ground_truth_uindex_list_l_3(2)%n = 13

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/49D0/))
            ground_truth_uindex_list_l_3(3)%n = 49

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/6D0/))
            ground_truth_uindex_list_l_3(4)%n = 6

            allocate(ground_truth_coefficient_list_l_3(5)%r, source = (/63D0/))
            ground_truth_uindex_list_l_3(5)%n = 63

            allocate(ground_truth_coefficient_list_l_3(6)%r, source = (/48D0/))
            ground_truth_uindex_list_l_3(6)%n = 48

            allocate(ground_truth_coefficient_list_l_3(7)%r, source = (/61D0/))
            ground_truth_uindex_list_l_3(7)%n = 61

            allocate(ground_truth_coefficient_list_l_3(8)%r, source = (/1D0/))
            ground_truth_uindex_list_l_3(8)%n = 1

            allocate(ground_truth_coefficient_list_l_3(9)%r, source = (/54D0/))
            ground_truth_uindex_list_l_3(9)%n = 54

            allocate(ground_truth_coefficient_list_l_3(10)%r, source = (/15D0/))
            ground_truth_uindex_list_l_3(10)%n = 15

            allocate(ground_truth_coefficient_list_l_3(11)%r, source = (/42D0/))
            ground_truth_uindex_list_l_3(11)%n = 42


            allocate(ground_truth_coefficient_list_l_2(7))
            allocate(ground_truth_uindex_list_l_2(7))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/1D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/6D0/))
            ground_truth_uindex_list_l_2(2)%n = 1

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/28D0/))
            ground_truth_uindex_list_l_2(3)%n = 3

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/97D0/))
            ground_truth_uindex_list_l_2(4)%n = 12

            allocate(ground_truth_coefficient_list_l_2(5)%r, source = (/54D0/))
            ground_truth_uindex_list_l_2(5)%n = 13

            allocate(ground_truth_coefficient_list_l_2(6)%r, source = (/124D0/))
            ground_truth_uindex_list_l_2(6)%n = 15

            allocate(ground_truth_coefficient_list_l_2(7)%r, source = (/42D0/))
            ground_truth_uindex_list_l_2(7)%n = 10

        end subroutine setup_ground_truth_C_coefficient_list_2D_2

        subroutine setup_ground_truth_D_tilde_coefficient_list_2D_2(ground_truth_coefficient_list_l_2, &
                                                                  ground_truth_uindex_list_l_2, &
                                                                  ground_truth_coefficient_list_l_3, &
                                                                  ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(6))
            allocate(ground_truth_uindex_list_l_3(6))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/34D0/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/35D0/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/70D0/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/0D0/))
            ground_truth_uindex_list_l_3(4)%n = 42

            allocate(ground_truth_coefficient_list_l_3(5)%r, source = (/28D0/))
            ground_truth_uindex_list_l_3(5)%n = 1

            allocate(ground_truth_coefficient_list_l_3(6)%r, source = (/272D0/))
            ground_truth_uindex_list_l_3(6)%n = 54


            allocate(ground_truth_coefficient_list_l_2(5))
            allocate(ground_truth_uindex_list_l_2(5))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/317D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/317D0/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/324D0/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/352D0/))
            ground_truth_uindex_list_l_2(4)%n = 10

            allocate(ground_truth_coefficient_list_l_2(5)%r, source = (/77D0/))
            ground_truth_uindex_list_l_2(5)%n = 13

        end subroutine setup_ground_truth_D_tilde_coefficient_list_2D_2

        subroutine setup_ground_truth_D_coefficient_list_2D_2(ground_truth_coefficient_list_l_2, &
                                                            ground_truth_uindex_list_l_2, &
                                                            ground_truth_coefficient_list_l_3, &
                                                            ground_truth_uindex_list_l_3)
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

            ! auxiliary

            allocate(ground_truth_coefficient_list_l_3(4))
            allocate(ground_truth_uindex_list_l_3(4))

            ground_truth_uindex_list_l_3(:)%l = 3

            allocate(ground_truth_coefficient_list_l_3(1)%r, source = (/351D0/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%r, source = (/352D0/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%r, source = (/394D0/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%r, source = (/352D0/))
            ground_truth_uindex_list_l_3(4)%n = 42

            allocate(ground_truth_coefficient_list_l_3(5)%r, source = (/345D0/))
            ground_truth_uindex_list_l_3(5)%n = 1

            allocate(ground_truth_coefficient_list_l_3(6)%r, source = (/349D0/))
            ground_truth_uindex_list_l_3(6)%n = 54


            allocate(ground_truth_coefficient_list_l_2(5))
            allocate(ground_truth_uindex_list_l_2(5))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%r, source = (/317D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%r, source = (/317D0/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%r, source = (/324D0/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%r, source = (/352D0/))
            ground_truth_uindex_list_l_2(4)%n = 10

            allocate(ground_truth_coefficient_list_l_2(5)%r, source = (/77D0/))
            ground_truth_uindex_list_l_2(5)%n = 13

        end subroutine setup_ground_truth_D_coefficient_list_2D_2

        subroutine setup_ground_truth_final_summation_2D_2(vector_v)
            implicit none
            ! dummy
            type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_v

            ! auxiliary
            integer :: i

            allocate(vector_v(14))
            do i=11, 14
                allocate(vector_v(i)%r(1))
            end do

            vector_v(11)%r(1) = 703
            vector_v(12)%r(1) = 704
            vector_v(13)%r(1) = 788
            vector_v(14)%r(1) = 696
            vector_v(15)%r(1) = 704
            vector_v(16)%r(1) = 808


        end subroutine setup_ground_truth_final_summation_2D_2

        function test_lib_ml_fmm_constructor() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            integer(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: i
            integer(kind=UINDEX_BYTES) :: list_index

            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_X_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_Y_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_XY_l_2

            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_X_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_Y_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_XY_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()


            ! --- setup ground truth data (2D) ---
            call setup_ground_truth_uindex_list_2D(ground_truth_uindex_list_X_l_2, &
                                                   ground_truth_uindex_list_Y_l_2, &
                                                   ground_truth_uindex_list_XY_l_2, &
                                                   ground_truth_uindex_list_X_l_3, &
                                                   ground_truth_uindex_list_Y_l_3, &
                                                   ground_truth_uindex_list_XY_l_3)

            ! --- test ---
            rv = .true.

            ! number of levels
            if (size(m_ml_fmm_hierarchy) .eq. 2) then
                ! level 3
                do i=1, size(ground_truth_uindex_list_X_l_3)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_X_l_3(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(3)%hierarchy_type(list_index) .ne. HIERARCHY_X) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do

                do i=1, size(ground_truth_uindex_list_Y_l_3)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_Y_l_3(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(3)%hierarchy_type(list_index) .ne. HIERARCHY_Y) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do

                do i=1, size(ground_truth_uindex_list_XY_l_3)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_XY_l_3(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(3)%hierarchy_type(list_index) .ne. HIERARCHY_XY) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do

                ! level 2
                do i=1, size(ground_truth_uindex_list_X_l_2)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_X_l_2(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(2)%hierarchy_type(list_index) .ne. HIERARCHY_X) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do

                do i=1, size(ground_truth_uindex_list_Y_l_2)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_Y_l_2(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(2)%hierarchy_type(list_index) .ne. HIERARCHY_Y) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do

                do i=1, size(ground_truth_uindex_list_XY_l_2)
                    list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, &
                                                                   ground_truth_uindex_list_XY_l_2(i))
                    if (list_index .gt. 0) then
                        if (m_ml_fmm_hierarchy(2)%hierarchy_type(list_index) .ne. HIERARCHY_XY) then
                            rv = .false.
                        else
                            ! correct
                            continue
                        end if
                    else
                        rv = .false.
                    end if
                end do
            else
                rv = .false.
            end if

            ! ~~~ test ~~~

            if (rv) then
                print *, "test_lib_ml_fmm_constructor (2D): OK"
            else
                print *, "test_lib_ml_fmm_constructor (2D): FAILED"
            end if
#else
            print *, "test_lib_ml_fmm_constructor (3D): NOT DEFINED"
            rv = .true.
#endif

!            allocate(vector_u(list_length+4_8))
!            allocate(dummy(1))
!            dummy(1) = 1.0
!            do i=1, size(vector_u)
!                vector_u(i)%r = dummy
!            end do
!
!            allocate(vector_v, source = lib_ml_fmm_run(vector_u))

        end function test_lib_ml_fmm_constructor

        function test_lib_ml_fmm_calculate_upward_pass() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliaray
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            ! --- setup ground truth data ---
            call setup_ground_truth_C_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                        ground_truth_uindex_list_l_2, &
                                                        ground_truth_coefficient_list_l_3, &
                                                        ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_upward_pass()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_upward_pass"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_C)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_C)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_upward_pass (3D): NOT DEFINED"
            rv = .true.
#endif


        end function test_lib_ml_fmm_calculate_upward_pass

        function test_lib_ml_fmm_calculate_downward_pass_1() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()

            ! --- setup ground truth data ---
            call setup_ground_truth_D_tilde_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                        ground_truth_uindex_list_l_2, &
                                                        ground_truth_coefficient_list_l_3, &
                                                        ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_downward_pass_step_1()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_downward_pass_1"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_downward_pass_1 (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_1

        function test_lib_ml_fmm_calculate_downward_pass_2() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()
            call lib_ml_fmm_calculate_downward_pass_step_1()

            ! --- setup ground truth data ---
            call setup_ground_truth_D_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                          ground_truth_uindex_list_l_2, &
                                                          ground_truth_coefficient_list_l_3, &
                                                          ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_downward_pass_step_2()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_downward_pass_2"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_downward_pass_2 (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_2

        function test_lib_ml_fmm_final_summation() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            integer(kind=UINDEX_BYTES) :: i

            type(lib_ml_fmm_v), dimension(:), allocatable :: ground_truth_vector_v
            double precision :: diff

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()
            call lib_ml_fmm_calculate_downward_pass_step_1()
            call lib_ml_fmm_calculate_downward_pass_step_2()

            ! --- setup ground truth data ---
            call setup_ground_truth_final_summation_2D(ground_truth_vector_v)

            ! --- function to test ---
            call lib_ml_fmm_final_summation()


            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_final_summation"
            do i=1, size(ground_truth_vector_v)
                if (allocated(m_ml_fmm_v(i)%r)) then
                    diff = m_ml_fmm_v(i)%r(1) - ground_truth_vector_v(i)%r(1)
                    if ( diff .lt. 1D-14 ) then
                        print *, "  m_ml_fmm_v(i=", i, "):  OK"
                    else
                        rv = .false.
                        print *, "  m_ml_fmm_v(i=", i, "):  FAILED"
                        print *, "    m_ml_fmm_v: ", m_ml_fmm_v(i)%r(1)
                        print *, "    GT        : ", ground_truth_vector_v(i)%r(1)
                    end if
                end if
            end do

#else
            print *, "test_lib_ml_fmm_final_summation (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_final_summation

        function test_lib_ml_fmm_calculate_upward_pass_v2() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliaray
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            ! --- setup ground truth data ---
            call setup_ground_truth_C_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                        ground_truth_uindex_list_l_2, &
                                                        ground_truth_coefficient_list_l_3, &
                                                        ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_upward_pass()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_upward_pass_v2"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_C)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_C)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_upward_pass_v2 (3D): NOT DEFINED"
            rv = .true.
#endif


        end function test_lib_ml_fmm_calculate_upward_pass_v2

        function test_lib_ml_fmm_calculate_downward_pass_1_v2() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()

            ! --- setup ground truth data ---
            call setup_ground_truth_D_tilde_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                        ground_truth_uindex_list_l_2, &
                                                        ground_truth_coefficient_list_l_3, &
                                                        ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_downward_pass_step_1()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_downward_pass_1_v2"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_downward_pass_1_v2 (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_1_v2

        function test_lib_ml_fmm_calculate_downward_pass_2_v2() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i

            type(lib_tree_universal_index) :: uindex
            integer(kind=1) :: coefficient_type

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_2

            type(lib_ml_fmm_coefficient), dimension(:), allocatable :: ground_truth_coefficient_list_l_3
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_list_l_3

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()
            call lib_ml_fmm_calculate_downward_pass_step_1()

            ! --- setup ground truth data ---
            call setup_ground_truth_D_coefficient_list_2D(ground_truth_coefficient_list_l_2, &
                                                          ground_truth_uindex_list_l_2, &
                                                          ground_truth_coefficient_list_l_3, &
                                                          ground_truth_uindex_list_l_3)

            ! --- function to test ---
            call lib_ml_fmm_calculate_downward_pass_step_2()

            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_calculate_downward_pass_v2"
            do i=1, size(ground_truth_coefficient_list_l_2)
                uindex = ground_truth_uindex_list_l_2(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_2(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D)) then
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 2) n:", uindex%n, ":  FAILED"
                end if
            end do

            do i=1, size(ground_truth_coefficient_list_l_3)
                uindex = ground_truth_uindex_list_l_3(i)
                coefficient = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
                if ((coefficient .eq. ground_truth_coefficient_list_l_3(i)) .and. &
                    (coefficient_type .eq. LIB_ML_FMM_COEFFICIENT_TYPE_D)) then
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  OK"
                else
                    rv = .false.
                    print *, "  uinedx(l = 3) n:", uindex%n, ":  FAILED"
                end if
            end do

#else
            print *, "test_lib_ml_fmm_calculate_downward_pass_v2 (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_2_v2

        function test_lib_ml_fmm_final_summation_v2() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            integer(kind=UINDEX_BYTES) :: i

            type(lib_ml_fmm_v), dimension(:), allocatable :: ground_truth_vector_v
            double precision :: diff

#if (_FMM_DIMENSION_ == 2)
            ! --- generate test data & setup the heriarchy ---
            call setup_hierarchy_with_test_data_2D()

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%r = dummy
            end do

            if (allocated(m_ml_fmm_u)) then
                m_ml_fmm_u = vector_u
            else
                allocate(m_ml_fmm_u, source=vector_u)
            end if

            call lib_ml_fmm_calculate_upward_pass()
            call lib_ml_fmm_calculate_downward_pass_step_1()
            call lib_ml_fmm_calculate_downward_pass_step_2()

            ! --- setup ground truth data ---
            call setup_ground_truth_final_summation_2D(ground_truth_vector_v)

            ! --- function to test ---
            call lib_ml_fmm_final_summation()


            ! --- test ---
            rv = .true.

            print *, "test_lib_ml_fmm_final_summation_v2"
            do i=1, size(ground_truth_vector_v)
                if (allocated(m_ml_fmm_v(i)%r)) then
                    diff = m_ml_fmm_v(i)%r(1) - ground_truth_vector_v(i)%r(1)
                    if ( diff .lt. 1D-14 ) then
                        print *, "  m_ml_fmm_v(i=", i, "):  OK"
                    else
                        rv = .false.
                        print *, "  m_ml_fmm_v(i=", i, "):  FAILED"
                        print *, "    m_ml_fmm_v: ", m_ml_fmm_v(i)%r(1)
                        print *, "    GT        : ", ground_truth_vector_v(i)%r(1)
                    end if
                end if
            end do

#else
            print *, "test_lib_ml_fmm_final_summation_v2 (3D): NOT DEFINED"
            rv = .true.
#endif
        end function test_lib_ml_fmm_final_summation_v2

    end function lib_ml_fmm_test_functions

end module lib_ml_fmm
