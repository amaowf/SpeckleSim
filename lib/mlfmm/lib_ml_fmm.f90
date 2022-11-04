#define _FMM_DIMENSION_ 2

! LIB: Mulitlevel Fast Multipole Method
!
module lib_ml_fmm
    use lib_tree
    use lib_tree_type
    use ml_fmm_type
    use ml_fmm_math
    use lib_ml_fmm_type_operator
    use lib_ml_fmm_helper_functions
    implicit none

    private

    ! --- public functions ---
    public :: lib_ml_fmm_test_functions

    public :: lib_ml_fmm_hf_test_functions

    ! --- member ---
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles
    integer(kind=2) :: m_ml_fmm_p_truncation

    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the X-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the Y-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_v

    type(lib_ml_fmm_hierarchy), dimension(:), allocatable :: m_ml_fmm_hierarchy

    ! Tree parameters
    integer(kind=UINDEX_BYTES) :: m_tree_neighbourhood_size_k
    integer(kind=4) :: m_tree_s_opt
    integer(kind=1) :: m_tree_l_min
    integer(kind=1) :: m_tree_l_max


    contains

    ! Procedure
    ! 1. call lib_tree_constructor
    !       - determine s_opt
    !       - determine l_min and l_max
    ! 2. create ml fmm data set
    !       - C per box and per level (l_max up to l_min)
    !       - D per box and per level (l_min down to l_max)
    subroutine lib_ml_fmm_constructor(data_elements)
        implicit none
        ! dummy
        type(lib_ml_fmm_data), intent(inout) :: data_elements

        ! auxiliaray
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated
        type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: correspondence_vector
        integer(kind=UINDEX_BYTES), dimension(3) :: length
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! concatenate X-, Y- and XY-hierarchy data
        !   length(1) = size(data_elements%X)
        !   length(2) = size(data_elements%Y)
        !   length(3) = size(data_elements%XY)
        allocate(data_concatenated, source=lib_ml_fmm_concatenate_data_array(data_elements, length))

        m_tree_s_opt = 1 ! todo: calculate s

        ! initiate the Tree
        call lib_tree_constructor(data_concatenated, m_tree_s_opt)
        data_concatenated = lib_tree_get_element_list()

        ! initiate the ml fmm type operators
        operator_procedures = ml_fmm_type_operator_get_procedures()
        call lib_ml_fmm_type_operator_constructor(operator_procedures)

        ! initiate the ml fmm procedures
        m_ml_fmm_handles = ml_fmm_get_procedures()

        ! adjust Tree parameters
        m_tree_neighbourhood_size_k = 1!lib_ml_fmm_hf_get_neighbourhood_size(R_c1, r_c2)
        m_tree_l_min = lib_tree_get_level_min(m_tree_neighbourhood_size_k)
        m_tree_l_max = lib_tree_get_level_max(m_tree_s_opt)

        ! initiate the X- and Y-hierarchy
        correspondence_vector = lib_tree_get_correspondence_vector()
        call lib_ml_fmm_hf_create_hierarchy(data_concatenated, correspondence_vector, &
                                            length, m_tree_l_min, m_tree_l_max,&
                                            m_ml_fmm_hierarchy)

        if (allocated(m_ml_fmm_u)) then
            deallocate(m_ml_fmm_u)
        end if
        if (allocated(m_ml_fmm_v)) then
            deallocate(m_ml_fmm_v)
        end if

        allocate (m_ml_fmm_u(length(HIERARCHY_X) + length(HIERARCHY_XY)))
        allocate (m_ml_fmm_v(length(HIERARCHY_Y) + length(HIERARCHY_XY)))
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
!!                do ii=1, size(m_ml_fmm_C(i)%dummy)
!!                    if (allocated (m_ml_fmm_C(i)%dummy))
!!                        deallocate (m_ml_fmm_C(i)%dummy)
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

    function lib_ml_fmm_run(vector_u) result(vector_v)
        implicit none
        ! dummy
        type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_u
        type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

        if (allocated(m_ml_fmm_u)) then
            m_ml_fmm_u = vector_u
        else
            allocate(m_ml_fmm_u, source=vector_u)
        end if

        call lib_ml_fmm_calculate_upward_pass()
        call lib_ml_fmm_calculate_downward_pass()
        call lib_ml_fmm_final_summation()

        allocate(vector_v, source=m_ml_fmm_v)

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

        call lib_ml_fmm_calculate_upward_pass_step_1()
        call lib_ml_fmm_calculate_upward_pass_step_2()
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
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((uindex%n .ge. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                    if (.not. ignore_box) then
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(i) = C
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end if
            end do
        else
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((list_index .gt. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                    if (.not. ignore_box) then
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(list_index) = C
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end if
            end do
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

        type(lib_tree_spatial_point) :: x_c

        integer(kind=UINDEX_BYTES) :: i

        call lib_ml_fmm_type_operator_set_coefficient_zero(C_i)

        data_element = lib_tree_get_domain_e1(uindex, element_number)
        if ((allocated (data_element)) &
            .and. (size(data_element) .gt. 0)) then
            ignore_box = .false.
            x_c = lib_tree_get_centre_of_box(uindex)
            do i=1, size(data_element)
                buffer_C_i = m_ml_fmm_u(element_number(i)) * m_ml_fmm_handles%get_B_i(x_c, data_element(i))
                C_i = C_i + buffer_C_i
            end do
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
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            C = lib_ml_fmm_get_C_of_box(uindex)
                            m_ml_fmm_hierarchy(l)%coefficient_list(i) = C
                            m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                        end if
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            C = lib_ml_fmm_get_C_of_box(uindex)
                            m_ml_fmm_hierarchy(l)%coefficient_list(list_index) = C
                            m_ml_fmm_hierarchy(l)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                        end if
                    end if
                end do
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

        x_c = lib_tree_get_centre_of_box(uindex)

        ! get child boxes
        uindex_children = lib_tree_get_children(uindex)

        call lib_ml_fmm_type_operator_set_coefficient_zero(C)
        do i=1, size(uindex_children)
            x_c_child = lib_tree_get_centre_of_box(uindex_children(i))
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex_children(i)%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C_child = m_ml_fmm_hierarchy(uindex_children(i)%l)%coefficient_list(list_index)
                    C = C + m_ml_fmm_handles%get_translation_SS(C_child, x_c_child, x_c)
                end if
            end if
        end do
    end function lib_ml_fmm_get_C_of_box

    subroutine lib_ml_fmm_calculate_downward_pass
        implicit none

        call lib_ml_fmm_calculate_downward_pass_step_1()
        call lib_ml_fmm_calculate_downward_pass_step_2()


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
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                            list_index_list(i) = i
                            coefficient_list(i) = D_tilde
                        end if
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                            list_index_list(i) = list_index
                            coefficient_list(i) = D_tilde
                        end if
                    end if
                end do
            end if

            ! wirte to hierarchy
            do i=1, number_of_boxes
                if (list_index_list(i) .gt. -1) then
                    m_ml_fmm_hierarchy(l)%coefficient_list(list_index_list(i)) = coefficient_list(i)
                    m_ml_fmm_hierarchy(l)%coefficient_type(list_index_list(i)) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
                end if
            end do

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

        allocate(boxes_i4, source=lib_tree_get_domain_i4(uindex, m_tree_neighbourhood_size_k))

        call lib_ml_fmm_type_operator_set_coefficient_zero(D_tilde)

        x_c = lib_tree_get_centre_of_box(uindex)

        do i=1, size(boxes_i4)
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, boxes_i4(i))
            if (list_index .gt. 0) then
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = m_ml_fmm_hierarchy(uindex%l)%coefficient_list(list_index)
                    x_c_neighbour = lib_tree_get_centre_of_box(boxes_i4(i))

                    D_tilde = D_tilde + m_ml_fmm_handles%get_translation_SR(C, x_c, x_c_neighbour)
                end if
            end if
        end do

    end function

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

        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: coefficient_list
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: list_index_list

        m_ml_fmm_hierarchy(m_tree_l_min)%coefficient_type(:) = LIB_ML_FMM_COEFFICIENT_TYPE_D

        do l=m_tree_l_min+int(1,1), m_tree_l_max
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            ! create coefficient buffer
            allocate(coefficient_list(number_of_boxes))
            allocate(list_index_list(number_of_boxes))
            list_index_list(:) = -1

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (uindex%n .ge. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D = lib_ml_fmm_get_D_of_box(uindex)
                            m_ml_fmm_hierarchy(l)%coefficient_list(i) = D
                            m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                        end if
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    if (list_index .gt. 0) then
                        hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(list_index)
                        if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                             (hierarchy_type .eq. HIERARCHY_XY))) then
                            D = lib_ml_fmm_get_D_of_box(uindex)
                            m_ml_fmm_hierarchy(l)%coefficient_list(i) = D
                            m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                        end if
                    end if
                end do
            end if

            ! wirte to hierarchy
            do i=1, number_of_boxes
                if (list_index_list(i) .gt. -1) then
                    m_ml_fmm_hierarchy(l)%coefficient_list(list_index_list(i)) = coefficient_list(i)
                    m_ml_fmm_hierarchy(l)%coefficient_type(list_index_list(i)) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
                end if
            end do

            ! clean up
            if (allocated(coefficient_list)) then
                deallocate(coefficient_list)
            end if
            if (allocated(list_index_list)) then
                deallocate(list_index_list)
            end if
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
    subroutine lib_ml_fmm_final_summation()
        implicit none

        ! auxiliary
        type(lib_tree_universal_index) :: uindex
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        uindex%l = m_tree_l_max
        number_of_boxes = size(m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index)

        if (m_ml_fmm_hierarchy(uindex%l)%is_hashed) then
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                if (uindex%n .ge. 0) then
                    hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                    if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then

                        call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
                    end if
                end if
            end do
        else
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                if (list_index .gt. 0) then
                    hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
                    if (((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then

                        call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
                    end if
                end if
            end do
        end if

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
    subroutine lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: D
        integer(kind=1) :: coefficient_type
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e1
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e1
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2

        integer(kind=UINDEX_BYTES) :: v_counter


        allocate(data_element_e1, source = lib_tree_get_domain_e1(uindex, &
                                                                  element_number_e1))
        allocate(data_element_e2, source = lib_tree_get_domain_e2(m_tree_neighbourhood_size_k, &
                                                                  uindex, &
                                                                  element_number_e2))

        D = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
        v_counter = 0
        do i=1, size(data_element_e1)
            if ((data_element_e1(i)%hierarchy .eq. HIERARCHY_Y) .or. &
                (data_element_e1(i)%hierarchy .eq. HIERARCHY_XY)) then
                v_counter = v_counter + 1
                m_ml_fmm_v(v_counter) = lib_ml_fmm_calculate_v_y_j(data_element_e1(i), data_element_e2, D)
            end if
        end do
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
    function lib_ml_fmm_calculate_v_y_j(data_element_y_j, data_element_e2, D) result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(inout) :: data_element_y_j
        type(lib_tree_data_element), dimension(:), allocatable, intent(inout) :: data_element_e2
        type(lib_ml_fmm_coefficient), intent(inout) :: D
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: y_j

        x_c = lib_tree_get_centre_of_box(data_element_y_j%uindex)
        y_j = data_element_y_j%point_x

        rv = D .cor. (y_j - x_c)

        do i=1, size(data_element_e2)
            rv = rv + m_ml_fmm_handles%get_phi_i_j(data_element_e2(i), y_j)
        end do

    end function lib_ml_fmm_calculate_v_y_j

    ! ----- test functions -----
    function lib_ml_fmm_test_functions() result(error_counter)
        implicit none
        ! dummy
        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_ml_fmm_constructor()) then
            error_counter = error_counter + 1
        end if
        call lib_ml_fmm_destructor()
        if (.not. test_lib_ml_fmm_calculate_upward_pass()) then
            error_counter = error_counter + 1
        end if
        call lib_ml_fmm_destructor()
        if (.not. test_lib_ml_fmm_calculate_downward_pass_1()) then
            error_counter = error_counter + 1
        end if
        call lib_ml_fmm_destructor()
        if (.not. test_lib_ml_fmm_calculate_downward_pass_2()) then
            error_counter = error_counter + 1
        end if
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

            call lib_ml_fmm_constructor(data_elements)

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

            allocate(ground_truth_coefficient_list_l_3(1)%dummy, source = (/0.0D0/))
            ground_truth_uindex_list_l_3(1)%n = 0

            allocate(ground_truth_coefficient_list_l_3(2)%dummy, source = (/13D0/))
            ground_truth_uindex_list_l_3(2)%n = 13

            allocate(ground_truth_coefficient_list_l_3(3)%dummy, source = (/49D0/))
            ground_truth_uindex_list_l_3(3)%n = 49

            allocate(ground_truth_coefficient_list_l_3(4)%dummy, source = (/6D0/))
            ground_truth_uindex_list_l_3(4)%n = 6

            allocate(ground_truth_coefficient_list_l_3(5)%dummy, source = (/63D0/))
            ground_truth_uindex_list_l_3(5)%n = 63

            allocate(ground_truth_coefficient_list_l_3(6)%dummy, source = (/48D0/))
            ground_truth_uindex_list_l_3(6)%n = 48

            allocate(ground_truth_coefficient_list_l_3(7)%dummy, source = (/61D0/))
            ground_truth_uindex_list_l_3(7)%n = 61

            allocate(ground_truth_coefficient_list_l_3(8)%dummy, source = (/1D0/))
            ground_truth_uindex_list_l_3(8)%n = 1

            allocate(ground_truth_coefficient_list_l_3(9)%dummy, source = (/54D0/))
            ground_truth_uindex_list_l_3(9)%n = 54

            allocate(ground_truth_coefficient_list_l_3(10)%dummy, source = (/15D0/))
            ground_truth_uindex_list_l_3(10)%n = 15


            allocate(ground_truth_coefficient_list_l_2(6))
            allocate(ground_truth_uindex_list_l_2(6))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%dummy, source = (/1D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%dummy, source = (/6D0/))
            ground_truth_uindex_list_l_2(2)%n = 1

            allocate(ground_truth_coefficient_list_l_2(3)%dummy, source = (/28D0/))
            ground_truth_uindex_list_l_2(3)%n = 3

            allocate(ground_truth_coefficient_list_l_2(4)%dummy, source = (/97D0/))
            ground_truth_uindex_list_l_2(4)%n = 12

            allocate(ground_truth_coefficient_list_l_2(5)%dummy, source = (/54D0/))
            ground_truth_uindex_list_l_2(5)%n = 13

            allocate(ground_truth_coefficient_list_l_2(6)%dummy, source = (/124D0/))
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

            allocate(ground_truth_coefficient_list_l_3(1)%dummy, source = (/34D0/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%dummy, source = (/35D0/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%dummy, source = (/28D0/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%dummy, source = (/0D0/))
            ground_truth_uindex_list_l_3(4)%n = 42


            allocate(ground_truth_coefficient_list_l_2(4))
            allocate(ground_truth_uindex_list_l_2(4))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%dummy, source = (/275D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%dummy, source = (/275D0/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%dummy, source = (/282D0/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%dummy, source = (/310D0/))
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

            allocate(ground_truth_coefficient_list_l_3(1)%dummy, source = (/309D0/))
            ground_truth_uindex_list_l_3(1)%n = 2

            allocate(ground_truth_coefficient_list_l_3(2)%dummy, source = (/310D0/))
            ground_truth_uindex_list_l_3(2)%n = 10

            allocate(ground_truth_coefficient_list_l_3(3)%dummy, source = (/310D0/))
            ground_truth_uindex_list_l_3(3)%n = 34

            allocate(ground_truth_coefficient_list_l_3(4)%dummy, source = (/310D0/))
            ground_truth_uindex_list_l_3(4)%n = 42


            allocate(ground_truth_coefficient_list_l_2(4))
            allocate(ground_truth_uindex_list_l_2(4))
            ground_truth_uindex_list_l_2(:)%l = 2

            allocate(ground_truth_coefficient_list_l_2(1)%dummy, source = (/275D0/))
            ground_truth_uindex_list_l_2(1)%n = 0

            allocate(ground_truth_coefficient_list_l_2(2)%dummy, source = (/275D0/))
            ground_truth_uindex_list_l_2(2)%n = 2

            allocate(ground_truth_coefficient_list_l_2(3)%dummy, source = (/282D0/))
            ground_truth_uindex_list_l_2(3)%n = 8

            allocate(ground_truth_coefficient_list_l_2(4)%dummy, source = (/310D0/))
            ground_truth_uindex_list_l_2(4)%n = 10

        end subroutine setup_ground_truth_D_coefficient_list_2D

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
#endif

!            allocate(vector_u(list_length+4_8))
!            allocate(dummy(1))
!            dummy(1) = 1.0
!            do i=1, size(vector_u)
!                vector_u(i)%dummy = dummy
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
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1
            type(lib_ml_fmm_data) :: data_elements

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=UINDEX_BYTES) :: list_index

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
                vector_u(i)%dummy = dummy
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
#endif


        end function test_lib_ml_fmm_calculate_upward_pass

        function test_lib_ml_fmm_calculate_downward_pass_1() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1
            type(lib_ml_fmm_data) :: data_elements

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=UINDEX_BYTES) :: list_index

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
                vector_u(i)%dummy = dummy
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
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_1

        function test_lib_ml_fmm_calculate_downward_pass_2() result(rv)
            implicit none
            ! dummy
            logical rv

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1
            type(lib_ml_fmm_data) :: data_elements

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            type(lib_ml_fmm_coefficient) :: coefficient
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=UINDEX_BYTES) :: list_index

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
                vector_u(i)%dummy = dummy
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
#endif
        end function test_lib_ml_fmm_calculate_downward_pass_2

    end function lib_ml_fmm_test_functions

end module lib_ml_fmm
