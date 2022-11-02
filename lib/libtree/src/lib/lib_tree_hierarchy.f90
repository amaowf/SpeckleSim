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
! Created on Fri Jun 12 11:42:51 2020
!
! @author: Max Daiber-Huppert

module lib_tree_hierarchy
    use libmath
    use lib_tree_type
    use lib_tree
    implicit none

    private

    ! --- public ---
    public :: lib_tree_hierarchy_class
    public :: lib_tree_hierarchy_test_functions

    public :: lib_tree_hierarchy_box


    ! --- parameter ---
    integer(kind=2), parameter :: MAXIMUM_NUMBER_OF_HASH_RUNS = 200
    integer(kind=2), public, parameter :: LIB_TREE_HIERARCHY_MARGIN = 220

    ! --- type definitions ----

    type lib_tree_hierarchy_box
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: element_number
        type(lib_tree_universal_index) :: uindex
        logical :: child_exist
        logical :: parent_exist
    end type

    type lib_tree_hierarchy_element_hashed
        integer(kind=UINDEX_BYTES) :: array_position
        integer(kind=2) :: number_of_hash_runs
    end type

    ! hierarchy of the tree data structure
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
    type lib_tree_hierarchy_level_type
        type(lib_tree_hierarchy_box), dimension(:), allocatable :: box_list
        ! [X, Y, XY] hierarchy type
        integer(kind=1), dimension(:), allocatable :: hierarchy_type
        ! [C, D^~, D]
!        integer(kind=1), dimension(:), allocatable :: element_type
        type(lib_tree_hierarchy_element_hashed), dimension(:), allocatable :: hashed_element_list_index
        ! is_hashed .eqv. .true.: list entry correspondence with the box index n
        ! is_hashed .eqv. .false: list entry correspondence with the coefficient_list_index index
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: element_list_index
        ! true: access coefficient with hashed uindex%n, false: access coefficient with uindex%n
        logical :: is_hashed
        ! number of not empty boxes
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=2) :: maximum_number_of_hash_runs
    end type lib_tree_hierarchy_level_type

!   !!! ATTENTION !!!
!
!   Only one instance of this class can be used !
!   *lib_tree* is NOT ready
!
    type lib_tree_hierarchy_class
        private
        type(lib_tree_hierarchy_level_type), dimension(:), allocatable :: hierarchy_level
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: tree_s_opt
        integer(kind=1) :: tree_l_min
        integer(kind=1) :: tree_l_max
        integer(kind=COORDINATE_BINARY_BYTES) :: tree_neighbourhood_size_k
        integer(kind=UINDEX_BYTES), dimension(3) :: length
    contains
        procedure :: constructor => lib_tree_hierarchy_constructor
        final :: lib_tree_hierarchy_destructor
        procedure :: get_box => lib_tree_hierarchy_get_hierarchy_box
        procedure :: get_l_max => lib_tree_hierarchy_get_l_max
    end type lib_tree_hierarchy_class

    ! --- type definitions ---
    type hash_list_type
        integer(kind=UINDEX_BYTES) :: value
        integer(kind=2) :: hash_runs
    end type

contains
    ! Argument
    ! ----
    !   self: lib_tree_hierarchy_class
    !       object of the class
    !   data_elements: lib_tree_data_element, dimension(:)
    !       list of data elements to create a tree hierarchy
    !   length: integer, dimension(3), optional (std: all data_elements are located at the X-hierarchy)
    !       number of elements per hierarchy (X, Y, XY)
    !       HINT: the elements have to be passed in the order X, Y, XY
    !   tree_s_opt: integer, optional (std: 1)
    !       number of elements located in one box at level l_max
    !   tree_l_min: integer, optional (std: neighbourhood size k dependent)
    !       highes level of the tree hierarchy
    subroutine lib_tree_hierarchy_constructor(self, data_elements, &
        length, tree_s_opt, tree_l_min, bounding_box)
        implicit none
        ! dummy
        class(lib_tree_hierarchy_class), intent(inout) :: self
        type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_elements
        integer(kind=UINDEX_BYTES), dimension(3), intent(in), optional :: length
        integer, intent(in), optional :: tree_s_opt
        integer(kind=1), intent(in), optional :: tree_l_min
        type(lib_tree_spatial_point), dimension(2), intent(in), optional :: bounding_box

        ! auxiliaray
!        type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: correspondence_vector

        if (present(tree_s_opt)) then
            self%tree_s_opt = int(tree_s_opt)
        else
            self%tree_s_opt = 1
        end if

        if (present(length)) then
            self%length = length
        else
            self%length(1) = size(data_elements, 1)
            self%length(2:3) = 0
        end if

        ! initiate the Tree
        if (present(bounding_box)) then
            call lib_tree_constructor(data_elements, self%tree_s_opt, bounding_box)
        else
            call lib_tree_constructor(data_elements, self%tree_s_opt)
        end if

        ! adjust Tree parameters
        self%tree_neighbourhood_size_k = 1!lib_tree_hierarchy_hf_get_neighbourhood_size(R_c1, r_c2)
        if (present(tree_l_min)) then
            self%tree_l_min = tree_l_min
        else
            self%tree_l_min = lib_tree_get_level_min(self%tree_neighbourhood_size_k)
        end if
        self%tree_l_max = lib_tree_get_level_max(self%tree_s_opt)

        if (self%tree_l_max .lt. self%tree_l_min) then
            self%tree_l_min = self%tree_l_max
        end if

        ! initiate the X- and Y-hierarchy
!        call lib_tree_hierarchy_create_hierarchy(self, lib_tree_get_element_list(), data_element_expansion)
        call lib_tree_hierarchy_create_hierarchy(self, lib_tree_get_element_list())

        print *, "lib_tree_hierarchy_constructor: INFO"
        print *, "  tree_l_min: ", self%tree_l_min
        print *, "  tree_l_max: ", self%tree_l_max
        print *, "  tree_s_opt: ", self%tree_s_opt

    end subroutine lib_tree_hierarchy_constructor

    ! clean up
    ! - coefficents
    ! - Tree
    subroutine lib_tree_hierarchy_destructor(self)
        implicit none
        type(lib_tree_hierarchy_class), intent(inout) :: self

        ! auxilary
        ! example for an indiviual deallocation
        !        integer(kind=1) :: i
        !        integer(kind=4) :: ii

        !        if (allocated(m_ml_fmm_hierarchy)) then
        !            deallocate(m_ml_fmm_hierarchy)
        !        end if

        call lib_tree_destructor()

    end subroutine lib_tree_hierarchy_destructor

!    ! Argument
!    ! ----
!    !   data: array<lib_tree_data_element>
!    !       list of data points
!    !   correspondence_vector: type(lib_tree_correspondece_vector_element)
!    !
!    !
!    ! Result
!    ! ----
!    !   hierarchy
!    !       hierarchy of the X- and Y-hierarchy (sources and targets)
!    subroutine lib_tree_hierarchy_create_hierarchy_with_possible_empty_boxes(self, data_elements, data_element_expansion)
!        implicit none
!        ! dummy
!        type(lib_tree_hierarchy_class), intent(inout) :: self
!        type(lib_tree_data_element), dimension(:), intent(in) :: data_elements
!        double precision, dimension(:), allocatable, intent(in) :: data_element_expansion
!
!        ! auxiliary
!        integer(kind=UINDEX_BYTES) :: i
!        type(lib_tree_universal_index) :: buffer_uindex
!        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_X
!        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_Y
!        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_XY
!        integer(kind=1) :: l_th
!        integer(kind=1) :: hierarchy_type
!        integer, dimension(3) :: uindex_list_counter
!
!        if (allocated(self%hierarchy_level)) then
!            deallocate(self%hierarchy_level)
!        end if
!        allocate( self%hierarchy_level(self%tree_l_min:self%tree_l_max))
!
!        allocate( uindex_list_X(self%length(HIERARCHY_X)) )
!        allocate( uindex_list_Y(self%length(HIERARCHY_Y)) )
!        allocate( uindex_list_XY(self%length(HIERARCHY_XY)) )
!        uindex_list_counter(:) = 0
!
!        ! iterate over all data elements (l = l_th)
!        do i=1, size(data_elements)
!            hierarchy_type = data_elements(i)%hierarchy
!            buffer_uindex = data_elements(i)%uindex
!            if (hierarchy_type .eq. HIERARCHY_X) then
!                uindex_list_counter(HIERARCHY_X) = uindex_list_counter(HIERARCHY_X) + 1
!                uindex_list_X(uindex_list_counter(HIERARCHY_X)) = buffer_uindex
!            else if (hierarchy_type .eq. HIERARCHY_Y) then
!                uindex_list_counter(HIERARCHY_Y) = uindex_list_counter(HIERARCHY_Y) + 1
!                uindex_list_Y(uindex_list_counter(HIERARCHY_Y)) = buffer_uindex
!            else if (hierarchy_type .eq. HIERARCHY_XY) then
!                uindex_list_counter(HIERARCHY_XY) = uindex_list_counter(HIERARCHY_XY) + 1
!                uindex_list_XY(uindex_list_counter(HIERARCHY_XY)) = buffer_uindex
!            else
!                print *, "lib_tree_hierarchy_create_hierarchy: ERROR"
!                print *, "  hierarchy not defined"
!            end if
!        end do
!
!        l_th = buffer_uindex%l
!
!        ! if necessary, get indices at level l_max
!        if (l_th .gt. self%tree_l_max) then
!            call lib_tree_make_parent_uindex_list(uindex_list_X, l_th - self%tree_l_max)
!            call lib_tree_make_parent_uindex_list(uindex_list_Y, l_th - self%tree_l_max)
!            call lib_tree_make_parent_uindex_list(uindex_list_XY, l_th - self%tree_l_max)
!
!            call lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique(uindex_list_X, uindex_list_Y, uindex_list_XY)
!        end if
!
!        ! setup hierarcy at level l_max
!        call lib_tree_hierarchy_add_uindex_to_hierarchy(self%hierarchy_level, self%tree_l_max, &
!            uindex_list_X, uindex_list_Y, uindex_list_XY)
!
!        ! setup hierary up to the level l_min
!        if (self%tree_l_max .gt. self%tree_l_min) then
!            do i=self%tree_l_max-1, self%tree_l_min, -1
!                call lib_tree_make_parent_uindex_list(uindex_list_X)
!                call lib_tree_make_parent_uindex_list(uindex_list_Y)
!                call lib_tree_make_parent_uindex_list(uindex_list_XY)
!
!                call lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique(uindex_list_X, uindex_list_Y, uindex_list_XY)
!
!                call lib_tree_hierarchy_add_uindex_to_hierarchy(self%hierarchy_level, int(i, 1), &
!                    uindex_list_X, uindex_list_Y, uindex_list_XY)
!            end do
!        else
!            print *, "lib_tree_hierarchy_create_hierarchy: ERROR"
!            print *, "  l_max .le. l_min"
!            print *, "  l_max = ", self%tree_l_max
!            print *, "  l_min = ", self%tree_l_min
!        end if
!
!        ! add data elements
!        do i = lbound(data_elements, 1), ubound(data_elements, 1)
!            call lib_tree_hierarchy_add_element(self%hierarchy_level, data_elements(i), i)
!        end do
!
!    end subroutine lib_tree_hierarchy_create_hierarchy_with_possible_empty_boxes

    ! Argument
    ! ----
    !   data: array<lib_tree_data_element>
    !       list of data points
    !   correspondence_vector: type(lib_tree_correspondece_vector_element)
    !
    !
    ! Result
    ! ----
    !   hierarchy
    !       hierarchy of the X- and Y-hierarchy (sources and targets)
    !
    ! Todo
    ! ----
    ! - [ ] use less memory: uindex_list_X, ...
    subroutine lib_tree_hierarchy_create_hierarchy(self, data_elements)
        implicit none
        ! dummy
        type(lib_tree_hierarchy_class), intent(inout) :: self
        type(lib_tree_data_element), dimension(:), intent(in) :: data_elements

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii
        type(lib_tree_universal_index) :: buffer_uindex
        type(lib_tree_universal_index), dimension(:, :), allocatable :: uindex_list_X
        type(lib_tree_universal_index), dimension(:, :), allocatable :: uindex_list_Y
        type(lib_tree_universal_index), dimension(:, :), allocatable :: uindex_list_XY

        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_X_per_level
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_Y_per_level
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_XY_per_level

        integer(kind=1), dimension(3) :: l_min
        integer(kind=1), dimension(3) :: l_max
        integer(kind=1) :: hierarchy_type
        integer, dimension(:,:), allocatable :: uindex_list_counter

        allocate(uindex_list_counter(3, 0:self%tree_l_max))

        allocate( uindex_list_X(self%length(HIERARCHY_X), 0:self%tree_l_max) )
        allocate( uindex_list_Y(self%length(HIERARCHY_Y), 0:self%tree_l_max) )
        allocate( uindex_list_XY(self%length(HIERARCHY_XY), 0:self%tree_l_max) )
        uindex_list_counter(:, :) = 0

        ! iterate over all data elements to get the expansion dependent universal index
        do i=1, size(data_elements)
            hierarchy_type = data_elements(i)%hierarchy
            if (hierarchy_type .eq. HIERARCHY_X) then
                if (data_elements(i)%expansion .gt. 0d0) then
                    buffer_uindex = lib_tree_hierarchy_get_uindex_expansion(data_elements(i), self%tree_l_min, self%tree_l_max)
                    if (buffer_uindex%l .lt. self%tree_l_min) then
                        self%tree_l_min = buffer_uindex%l
                    end if
                else
                    buffer_uindex = data_elements(i)%uindex
                end if
                uindex_list_counter(HIERARCHY_X, buffer_uindex%l) = &
                    uindex_list_counter(HIERARCHY_X, buffer_uindex%l) + 1
                uindex_list_X(uindex_list_counter(HIERARCHY_X, buffer_uindex%l), buffer_uindex%l) = buffer_uindex

            else if (hierarchy_type .eq. HIERARCHY_Y) then
                if (data_elements(i)%expansion .gt. 0d0) then
                    buffer_uindex = lib_tree_hierarchy_get_uindex_expansion(data_elements(i), self%tree_l_min, self%tree_l_max)
                    if (buffer_uindex%l .lt. self%tree_l_min) then
                        self%tree_l_min = buffer_uindex%l
                    end if
                else
                    buffer_uindex = data_elements(i)%uindex
                end if
                uindex_list_counter(HIERARCHY_Y, buffer_uindex%l) = &
                    uindex_list_counter(HIERARCHY_Y, buffer_uindex%l) + 1
                uindex_list_Y(uindex_list_counter(HIERARCHY_Y, buffer_uindex%l), buffer_uindex%l) = buffer_uindex
            else if (hierarchy_type .eq. HIERARCHY_XY) then
                if (data_elements(i)%expansion .gt. 0d0) then
                    buffer_uindex = lib_tree_hierarchy_get_uindex_expansion(data_elements(i), self%tree_l_min, self%tree_l_max)
                    if (buffer_uindex%l .lt. self%tree_l_min) then
                        self%tree_l_min = buffer_uindex%l
                    end if
                else
                    buffer_uindex = data_elements(i)%uindex
                end if
                uindex_list_counter(HIERARCHY_XY, buffer_uindex%l) = &
                    uindex_list_counter(HIERARCHY_XY, buffer_uindex%l) + 1
                uindex_list_XY(uindex_list_counter(HIERARCHY_XY, buffer_uindex%l), buffer_uindex%l) = buffer_uindex
            else
                print *, "lib_tree_hierarchy_create_hierarchy: ERROR"
                print *, "  hierarchy not defined"
            end if
        end do

        ! get tree limits: l_min, l_max
        l_min(:) = IGNORE_ENTRY
        l_max(:) = self%tree_l_max
        do ii = HIERARCHY_X, HIERARCHY_XY
            do i = 0, self%tree_l_max
                if (uindex_list_counter(ii, i) .gt. 0) then
                    l_min(ii) = int(i, kind=1)
                    exit
                end if
            end do

            do i = self%tree_l_max, 0, -1
                if (uindex_list_counter(ii, i) .gt. 0) then
                    l_max(ii) = int(i, kind=1)
                    exit
                end if
            end do

            if (l_min(ii) .lt. self%tree_l_min &
                .and. l_min(ii) .ne. IGNORE_ENTRY) then
                self%tree_l_min = l_min(ii)
            end if

            if (l_max(ii) .lt. self%tree_l_max &
                .and. l_max(ii) .ne. IGNORE_ENTRY) then
                self%tree_l_max = l_max(ii)
            end if
        end do

        ! add parent boxes
        do i=self%tree_l_max, self%tree_l_min+1, -1

            ! --- X hierarchy ---
            if (uindex_list_counter(HIERARCHY_X, i) .gt. 0) then
                allocate(uindex_list_X_per_level(uindex_list_counter(HIERARCHY_X, i)))
                allocate(uindex_list_Y_per_level(uindex_list_counter(HIERARCHY_X, i-1)))

                ! current boxes
                uindex_list_X_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_X, i), i)
                ! current parent boxes, maybe without the parents of the current boxes
                uindex_list_Y_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_X, i-1), i-1)

                call lib_tree_hierarchy_add_parent_box(uindex_list_X_per_level, &
                                                       uindex_list_Y_per_level)

                uindex_list_counter(HIERARCHY_X, i-1) = size(uindex_list_Y_per_level)

                uindex_list_X(1:uindex_list_counter(HIERARCHY_X, i-1), i-1) = uindex_list_Y_per_level

                deallocate(uindex_list_X_per_level)
                deallocate(uindex_list_Y_per_level)
            end if

            ! --- Y hierarchy ---
            if (uindex_list_counter(HIERARCHY_Y, i) .gt. 0) then
                allocate(uindex_list_X_per_level(uindex_list_counter(HIERARCHY_Y, i)))
                allocate(uindex_list_Y_per_level(uindex_list_counter(HIERARCHY_Y, i-1)))

                ! current boxes
                uindex_list_X_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_Y, i), i)
                ! current parent boxes, maybe without the parents of the current boxes
                uindex_list_Y_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_Y, i-1), i-1)

                call lib_tree_hierarchy_add_parent_box(uindex_list_X_per_level, &
                                                       uindex_list_Y_per_level)

                uindex_list_counter(HIERARCHY_Y, i-1) = size(uindex_list_Y_per_level)

                uindex_list_Y(1:uindex_list_counter(HIERARCHY_Y, i-1), i-1) = uindex_list_Y_per_level

                deallocate(uindex_list_X_per_level)
                deallocate(uindex_list_Y_per_level)
            end if

            ! --- XY hierarchy ---
            if (uindex_list_counter(HIERARCHY_XY, i) .gt. 0) then
                allocate(uindex_list_X_per_level(uindex_list_counter(HIERARCHY_XY, i)))
                allocate(uindex_list_Y_per_level(uindex_list_counter(HIERARCHY_XY, i-1)))

                ! current boxes
                uindex_list_X_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_XY, i), i)
                ! current parent boxes, maybe without the parents of the current boxes
                uindex_list_Y_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_XY, i-1), i-1)

                call lib_tree_hierarchy_add_parent_box(uindex_list_X_per_level, &
                                                       uindex_list_Y_per_level)

                uindex_list_counter(HIERARCHY_XY, i-1) = size(uindex_list_XY_per_level)

                uindex_list_XY(1:uindex_list_counter(HIERARCHY_XY, i-1), i-1) = uindex_list_Y_per_level

                deallocate(uindex_list_X_per_level)
                deallocate(uindex_list_Y_per_level)
            end if
        end do


        ! create hierarchy
        if (allocated(self%hierarchy_level)) then
            deallocate(self%hierarchy_level)
        end if
        allocate( self%hierarchy_level(self%tree_l_min:self%tree_l_max))

        do i=self%tree_l_min, self%tree_l_max

            allocate(uindex_list_X_per_level(uindex_list_counter(HIERARCHY_X, i)))
            allocate(uindex_list_Y_per_level(uindex_list_counter(HIERARCHY_Y, i)))
            allocate(uindex_list_XY_per_level(uindex_list_counter(HIERARCHY_XY, i)))

            uindex_list_X_per_level = uindex_list_X(1:uindex_list_counter(HIERARCHY_X, i), i)
            uindex_list_Y_per_level = uindex_list_Y(1:uindex_list_counter(HIERARCHY_Y, i), i)
            uindex_list_XY_per_level = uindex_list_XY(1:uindex_list_counter(HIERARCHY_XY, i), i)



            call lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique(uindex_list_X_per_level, &
                                                                   uindex_list_Y_per_level, &
                                                                   uindex_list_XY_per_level)

            call lib_tree_hierarchy_add_uindex_to_hierarchy(self%hierarchy_level, int(i, 1), &
                                                            uindex_list_X_per_level, &
                                                            uindex_list_Y_per_level, &
                                                            uindex_list_XY_per_level)
            deallocate(uindex_list_X_per_level)
            deallocate(uindex_list_Y_per_level)
            deallocate(uindex_list_XY_per_level)
        end do

        ! add data elements
        do i = lbound(data_elements, 1), ubound(data_elements, 1)
            call lib_tree_hierarchy_add_element(self%hierarchy_level, data_elements(i), i)
        end do

        ! set *parent_exist* and *child_exist*
        call lib_tree_hierarchy_set_existence_of_parent_and_child(self%hierarchy_level)

    end subroutine lib_tree_hierarchy_create_hierarchy

    function lib_tree_hierarchy_get_l_max(self) result(rv)
        implicit none
        ! dummy
        class(lib_tree_hierarchy_class), intent(inout) :: self
        integer(kind=1) :: rv

        rv = self%tree_l_max

    end function lib_tree_hierarchy_get_l_max

    subroutine lib_tree_hierarchy_add_parent_box(uindex_child, uindex_parent)
        implicit none
        ! dummy
        type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_child
        type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_parent

        ! auxiliary
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_uindex
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: n
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: n_1
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: n_2

        if (allocated(uindex_child)) then

            buffer_uindex = lib_tree_get_parent_uindex_list(uindex_child)

            allocate(n_1(size(uindex_parent, 1)))
            n_1 = uindex_parent(:)%n

            allocate(n_2(size(buffer_uindex, 1)))
            n_2 = buffer_uindex(:)%n

            n = get_union(n_1, n_2)

            deallocate(n_1)
            deallocate(n_2)

            deallocate(uindex_parent)
            allocate(uindex_parent(size(n)))

            uindex_parent(:)%l = uindex_child(lbound(uindex_child, 1))%l - int(1, kind = 1)
            uindex_parent(:)%n = n
        end if

    end subroutine lib_tree_hierarchy_add_parent_box

    subroutine lib_tree_hierarchy_set_existence_of_parent_and_child(hierarchy)
        implicit none
        ! dummy
        type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy

        ! auxiliary
        integer :: l
        integer :: i
        integer(kind = UINDEX_BYTES) :: box_index
        type(lib_tree_universal_index), dimension(:), allocatable :: a

        hierarchy(ubound(hierarchy, 1))%box_list(:)%child_exist = .false.
        hierarchy(ubound(hierarchy, 1))%box_list(:)%parent_exist = .false.

        hierarchy(lbound(hierarchy, 1))%box_list(:)%child_exist = .false.
        hierarchy(lbound(hierarchy, 1))%box_list(:)%parent_exist = .false.

        do l = ubound(hierarchy, 1), lbound(hierarchy, 1) + 1, -1
            hierarchy(l)%box_list(:)%parent_exist = .false.
            hierarchy(l-1)%box_list(:)%child_exist = .false.

            allocate(a(size(hierarchy(l)%box_list)))
            a = hierarchy(l)%box_list(:)%uindex
            call lib_tree_make_parent_uindex_list(a)

            do i = lbound(a, 1), ubound(a, 1)
                box_index = lib_tree_hierarchy_get_hierarchy_index(hierarchy, a(i))

                if (box_index .eq. IGNORE_ENTRY) then
                    ! box (level l) has no parent box
                    print *, "lib_tree_hierarchy_set_existence_of_parent_and_child: ERROR"
                    print *, "  box without a parent box"
                    print *, "  l: ", a(i)%l
                    print *, "  n: ", a(i)%n
                else
                    hierarchy(l)%box_list(i)%parent_exist = .true.
                    hierarchy(l-1)%box_list(box_index)%child_exist = .true.
                end if
            end do

            deallocate(a)
        end do

    end subroutine lib_tree_hierarchy_set_existence_of_parent_and_child

    subroutine lib_tree_hierarchy_get_uindex_sets(a, b, &
                                                  intersection, &
                                                  relative_complement_of_b_in_a, &
                                                  relative_complement_of_a_in_b)
        implicit none
        ! dummy
        type(lib_tree_universal_index), dimension(:), allocatable, intent(in) :: a
        type(lib_tree_universal_index), dimension(:), allocatable, intent(in) :: b
        type(lib_tree_universal_index), dimension(:), allocatable, intent(out) :: intersection
        type(lib_tree_universal_index), dimension(:), allocatable, intent(out) :: relative_complement_of_b_in_a
        type(lib_tree_universal_index), dimension(:), allocatable, intent(out) :: relative_complement_of_a_in_b

        ! auxiliary
        integer(kind = UINDEX_BYTES), dimension(:), allocatable :: a_n
        integer(kind = UINDEX_BYTES), dimension(:), allocatable :: b_n
        integer(kind = UINDEX_BYTES), dimension(:), allocatable :: intersection_n
        integer(kind = UINDEX_BYTES), dimension(:), allocatable :: b_in_a_n
        integer(kind = UINDEX_BYTES), dimension(:), allocatable :: a_in_b_n

        allocate(a_n(size(a)))
        allocate(b_n(size(b)))

        a_n = a(:)%n
        b_n = b(:)%n

        call get_sets(a_n, b_n, intersection_n, b_in_a_n, a_in_b_n)

        deallocate(a_n)
        deallocate(b_n)

        allocate(intersection(size(intersection_n)))
        allocate(relative_complement_of_b_in_a(size(b_in_a_n)))
        allocate(relative_complement_of_a_in_b(size(a_in_b_n)))

        intersection(:)%l = a(lbound(a, 1))%l
        intersection(:)%n = intersection_n

        relative_complement_of_b_in_a(:)%l = a(lbound(a, 1))%l
        relative_complement_of_b_in_a(:)%n = b_in_a_n

        relative_complement_of_a_in_b(:)%l = a(lbound(a, 1))%l
        relative_complement_of_a_in_b(:)%n = a_in_b_n

    end subroutine lib_tree_hierarchy_get_uindex_sets


    function lib_tree_hierarchy_get_uindex_expansion(data_element, l_min, l_max) result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(in) :: data_element
        integer(kind = 1), intent(in) :: l_min
        integer(kind = 1), intent(in) :: l_max
        type(lib_tree_universal_index) :: rv

        ! auxiliary
        rv%l = lib_tree_get_level_from_edge_length(data_element%expansion)
        if (rv%l .gt. l_max) then
            if (data_element%uindex%l .gt. l_max) then
                rv = lib_tree_get_parent(data_element%uindex, &
                                         data_element%uindex%l - l_max)
            else
                rv = data_element%uindex
            end if
        else if (rv%l .lt. l_min) then
            rv = lib_tree_get_parent(data_element%uindex, &
                                     abs(rv%l - data_element%uindex%l))
        else
            rv = lib_tree_get_parent(data_element%uindex, &
                                     abs(rv%l - data_element%uindex%l))
        end if
    end function lib_tree_hierarchy_get_uindex_expansion

    ! Adds lists of universal indices of different hierachical type to the hierarchy.
    !
    ! Argument
    ! ----
    !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
    !       reference to the hierarchy object to be processed
    !   level: integer, in
    !       level of the hierarchy to be processed
    !   uindex_list_[X,Y,XY]: array<lib_tree_universal_index>, inout (read-only)
    !       list of uindex' of the type [X,Y,XY] to add to the hierarchy
    !
    subroutine lib_tree_hierarchy_add_uindex_to_hierarchy(hierarchy, level, uindex_list_X, uindex_list_Y, uindex_list_XY)
        implicit none
        ! dummy
        type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy
        integer(kind=1), intent(in) :: level
        type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_X
        type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_Y
        type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_XY

        !auxiliary
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: number_of_entries_log_2
        integer(kind=UINDEX_BYTES) :: counter
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: max_value
        logical :: no_hash

        number_of_boxes = size(uindex_list_X) &
            + size(uindex_list_Y) &
            + size(uindex_list_XY)
        hierarchy(level)%number_of_boxes = number_of_boxes

        allocate (hierarchy(level)%box_list(number_of_boxes))
        allocate (hierarchy(level)%hierarchy_type(number_of_boxes))
!        allocate (hierarchy(level)%element_type(number_of_boxes))

        hierarchy(level)%hierarchy_type(:) = IGNORE_ENTRY

        ! 2**(dl) < 256
        if (log(real(level)) / log(2.0) .lt. 8) then
            no_hash = .true.
        else
            no_hash = .false.
        end if
        !            no_hash = .false.

        ! calculates whether hash access or direct access via the universal index requires less memory
        !
        ! number_of_boxes * margin/100.0 < 2**(dimension * level)
        number_of_entries_log_2 = ceiling(log(real(number_of_boxes, 8) * LIB_TREE_HIERARCHY_MARGIN/100.0D0) / log(2.0D0))
        if ((number_of_entries_log_2 .lt. TREE_DIMENSIONS * level) .and. &
            .not. no_hash) then
            hierarchy(level)%is_hashed = .true.
            max_value = ceiling(real(number_of_boxes, 8) * LIB_TREE_HIERARCHY_MARGIN/100.0D0)
            allocate (hierarchy(level)%hashed_element_list_index(max_value))
            hierarchy(level)%hashed_element_list_index(:)%number_of_hash_runs = IGNORE_ENTRY
            hierarchy(level)%maximum_number_of_hash_runs = 0
            allocate (hierarchy(level)%element_list_index(number_of_boxes))
            hierarchy(level)%element_list_index(:) = IGNORE_ENTRY

            counter = 0
            call lib_tree_hierarchy_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_X, HIERARCHY_X, counter)
            call lib_tree_hierarchy_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_Y, HIERARCHY_Y, counter)
            call lib_tree_hierarchy_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_XY, HIERARCHY_XY, counter)
        else
            hierarchy(level)%is_hashed = .false.
            i = 2**(TREE_DIMENSIONS * level)
            allocate (hierarchy(level)%element_list_index(i))

            hierarchy(level)%element_list_index(:) = IGNORE_ENTRY

            ! setup the arrays coeeficient_list_index and hierarchy_type
            counter = 0

            call lib_tree_hierarchy_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_X, HIERARCHY_X, counter)
            call lib_tree_hierarchy_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_Y, HIERARCHY_Y, counter)
            call lib_tree_hierarchy_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_XY, HIERARCHY_XY, counter)
        end if

    end subroutine lib_tree_hierarchy_add_uindex_to_hierarchy

    ! Add the uindex to the hierarchy by using a hash function
    !
    ! Arguments
    ! ----
    !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
    !       reference to the hierarchy object to be processed
    !   level: integer, in
    !       level of the hierarchy to be processed
    !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
    !       list of uindex' to add to the hierarchy
    !   hierarchy_type: integer, in
    !       type of the uindex_list [HIERARCHY_X, HIERARCHY_Y, HIERARCHY_XY]
    !   counter: integer, inout
    !       counter of added uindex' to the hierarchy()
    !
    subroutine lib_tree_hierarchy_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list, hierarchy_type, counter)
        implicit none
        ! dummy
        type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy
        integer(kind=1), intent(in) :: level
        type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list
        integer(kind=1), intent(in) :: hierarchy_type
        integer(kind=UINDEX_BYTES) :: counter

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=2) :: ii
        integer(kind=UINDEX_BYTES) :: hash
        integer(kind=UINDEX_BYTES) :: max_value

        max_value = size(hierarchy(level)%hashed_element_list_index)

        ! interate over all indices
        do i=1, size(uindex_list)
            ! determine a valid hash
            hash = hash_fnv1a(uindex_list(i)%n, max_value)
            do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                if (hierarchy(level)%hashed_element_list_index(hash)%number_of_hash_runs .eq. IGNORE_ENTRY) then
                    ! saving the uindex in the hierarchy
                    counter = counter + 1
                    hierarchy(level)%hashed_element_list_index(hash)%number_of_hash_runs = ii
                    hierarchy(level)%hashed_element_list_index(hash)%array_position = counter
                    hierarchy(level)%element_list_index(counter) = uindex_list(i)%n
                    hierarchy(level)%hierarchy_type(counter) = hierarchy_type

                    if (hierarchy(level)%maximum_number_of_hash_runs .lt. ii) then
                        hierarchy(level)%maximum_number_of_hash_runs = ii
                    end if

                    exit
                else if (hierarchy(level)%element_list_index(counter+1) .ne. uindex_list(i)%n) then
                    ! hash collision -> re-hash
                    hash = IEOR(hash, int(ii, UINDEX_BYTES))
                    hash = hash_fnv1a(hash, max_value)
                else
                    ! uindex exists already in the hierarchy
                    if (hierarchy(level)%hierarchy_type(counter+1) .ne. hierarchy_type) then
                        ! hierarchy type differ -> merge hierarchy types
                        hierarchy(level)%hierarchy_type(counter+1) = HIERARCHY_XY
                    end if
#ifdef _DEBUG_
                    print *, "lib_tree_hierarchy_add_uindex_to_hierarchy_hashed: NOTE"
                    print *, "  uindex bypassed at level: ", level
                    print *, "  hierarchy_type", hierarchy_type
                    print *, "  uindex_list(i)%n: ", uindex_list(i)%n
#endif
                end if
            end do
            if (ii .gt. MAXIMUM_NUMBER_OF_HASH_RUNS) then
                print *, "lib_tree_hierarchy_add_uindex_to_hierarchy_hashed: ERROR"
                print *, "  uindex bypassed at level: ", level
                print *, "  hierarchy_type", hierarchy_type
                print *, "  uindex_list(i)%n: ", uindex_list(i)%n
            end if
        end do

    end subroutine lib_tree_hierarchy_add_uindex_to_hierarchy_hashed

    ! Add the uindex to the hierarchy by using a lookup table
    !
    ! Arguments
    ! ----
    !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
    !       reference to the hierarchy object to be processed
    !   level: integer, in
    !       level of the hierarchy to be processed
    !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
    !       list of uindex' to add to the hierarchy
    !   hierarchy_type: integer, in
    !       type of the uindex_list [HIERARCHY_X, HIERARCHY_Y, HIERARCHY_XY]
    !   counter: integer, inout
    !       counter of added uindex' to the hierarchy()
    !
    subroutine lib_tree_hierarchy_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list, hierarchy_type, counter)
        implicit none
        ! dummy
        type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy
        integer(kind=1), intent(in) :: level
        type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list
        integer(kind=1), intent(in) :: hierarchy_type
        integer(kind=UINDEX_BYTES) :: counter

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: n

        do i=1, size(uindex_list)
            counter = counter + 1
            n = uindex_list(i)%n + 1
            if (hierarchy(level)%element_list_index(n) .eq. IGNORE_ENTRY) then
                hierarchy(level)%element_list_index(n) = counter
                hierarchy(level)%hierarchy_type(counter) = hierarchy_type
                hierarchy(level)%box_list(counter)%uindex = uindex_list(i)
            else
                if (hierarchy(level)%hierarchy_type(counter) .ne. hierarchy_type) then
                    hierarchy(level)%hierarchy_type(counter) = HIERARCHY_XY
                end if
#ifdef _DEBUG_
                print *, "lib_tree_hierarchy_add_uindex_to_hierarchy: NOTE"
                print *, "  uindex bypassed at level: ", level
                print *, "  hierarchy_type", hierarchy_type
                print *, "  uindex_list(i)%n: ", n
#endif
            end if
        end do

    end subroutine lib_tree_hierarchy_add_uindex_to_hierarchy_lookup_table

    ! Argument
    ! ----
    !   hierarchy: lib_tree_hierarchy_level_type, dimension(:)
    !       levels of the tree hierarchy
    !   element: lib_tree_data_element
    !       element to be added
    !   element_number: integer
    !       element index of the *lib_tree_data_element_list* list
    !   element_expansion: double precision
    !       spherical expansion of the data point (element)
    !       0.0: Element is saved at level l_max.
    !       >0.0: Element is saved on the level <= l_max with the restriction that the element only occupies one box.
    subroutine lib_tree_hierarchy_add_element(hierarchy, element, element_number)
            implicit none
            ! dummy
            type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_tree_data_element), intent(in) :: element
            integer(kind=UINDEX_BYTES), intent(in) :: element_number

            ! auxiliary
            type(lib_tree_universal_index) :: uindex

            if (element%expansion .gt. 0d0) then
                uindex%l = lib_tree_get_level_from_edge_length(element%expansion)
                uindex = lib_tree_get_parent(element%uindex, abs(uindex%l - element%uindex%l))
            else
                uindex = element%uindex
            end if

            call lib_tree_hierarchy_add_element_number(hierarchy, uindex, element_number)

    end subroutine lib_tree_hierarchy_add_element

    ! Argument
    ! ----
    !   uindex: lib_tree_universal_index
    !       universal index of the element
    !   element_number: integer
    !       number of the element of the *data_elements* list at the constructor
    !
    ! Returns
    ! ----
    !   hierarchy: lib_tree_hierarchy_level_type, dimension(:)
    !       Where the element numbers are added
    subroutine lib_tree_hierarchy_add_element_number(hierarchy, uindex, element_number)
            implicit none
            ! dummy
            type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_tree_universal_index), intent(in) :: uindex
            integer(kind=UINDEX_BYTES), intent(in) :: element_number

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: list_index
            integer(kind=UINDEX_BYTES) :: buffer
            integer(kind=UINDEX_BYTES), dimension(:), allocatable :: buffer_list

            if (uindex%n .ge. 0) then
                list_index = lib_tree_hierarchy_get_hierarchy_index(hierarchy, uindex)

                if (allocated(hierarchy(uindex%l)%box_list(list_index)%element_number)) then
                    call move_alloc(hierarchy(uindex%l)%box_list(list_index)%element_number, buffer_list)
                    buffer = size(buffer_list, 1)
                    allocate(hierarchy(uindex%l)%box_list(list_index)%element_number(buffer + 1))
                    hierarchy(uindex%l)%box_list(list_index)%element_number(1:buffer) = buffer_list
                    hierarchy(uindex%l)%box_list(list_index)%element_number(buffer + 1) = element_number
                else
                    allocate(hierarchy(uindex%l)%box_list(list_index)%element_number(1))
                    hierarchy(uindex%l)%box_list(list_index)%element_number(1) = element_number
                end if
            end if
        end subroutine lib_tree_hierarchy_add_element_number

        ! Argument
        ! ----
        !   hierarchy: lib_tree_hierarchy_level_type, dimension(:)
        !       levels of the hierarchy with the element lists
        !   uindex: lib_tree_universal_index
        !       universal index of the box
        !
        ! Returns
        ! ----
        !   box: lib_tree_hierarchy_box
        !       list of elements in the box, if available
        function lib_tree_hierarchy_get_hierarchy_box(self, uindex) result(box)
            implicit none
            ! dummy
            class(lib_tree_hierarchy_class), intent(in) :: self
!            type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(in) :: hierarchy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_tree_hierarchy_box) :: box

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: list_index

            if (uindex%n .ge. 0 &
                .and. uindex%l .ge. lbound(self%hierarchy_level, 1) &
                .and. uindex%l .le. ubound(self%hierarchy_level, 1)) then
                list_index = lib_tree_hierarchy_get_hierarchy_index(self%hierarchy_level, uindex)
                if (list_index .ne. IGNORE_ENTRY) then
                    box = self%hierarchy_level(uindex%l)%box_list(list_index)
                else
                    box%uindex%l = IGNORE_ENTRY
                    box%uindex%n = IGNORE_ENTRY
                end if
            else
                box%uindex%l = IGNORE_ENTRY
                box%uindex%n = IGNORE_ENTRY
!                print *, "lib_tree_hierarchy_get_hierarchy_box: ERROR"
!                print *, "  index is not valid"
            end if
        end function lib_tree_hierarchy_get_hierarchy_box

        ! Argument
        ! ----
        !   hierarchy: lib_tree_hierarchy_level_type, dimension(:)
        !       levels of the hierarchy
        !   uindex: lib_tree_universal_index
        !       univesal index of the box
        !
        ! Returns
        ! ----
        !   rv: integer
        !       index of the *lib_tree_hierarchy_element* list
        function lib_tree_hierarchy_get_hierarchy_index(hierarchy, uindex) result(rv)
            implicit none
            ! dummy
            type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(in) :: hierarchy
            type(lib_tree_universal_index), intent(in) :: uindex
            integer(UINDEX_BYTES) :: rv

            if (hierarchy(uindex%l)%is_hashed) then
                rv = lib_tree_hierarchy_get_index_hierarchy_hashed(hierarchy, uindex)
            else
                rv = hierarchy(uindex%l)%element_list_index(uindex%n + int(1,1))
            end if

        end function lib_tree_hierarchy_get_hierarchy_index

        ! Argument
        ! ----
        !   hierarchy: lib_tree_hierarchy_level_type, dimension(:)
        !       levels of the hierarchy
        !   uindex: lib_tree_universal_index
        !       univesal index of the box
        !
        ! Returns
        ! ----
        !   rv: integer
        !       index of the *lib_tree_hierarchy_element_hashed* list
        function lib_tree_hierarchy_get_index_hierarchy_hashed(hierarchy, uindex) result(rv)
            implicit none
            ! dummy
            type(lib_tree_hierarchy_level_type), dimension(:), allocatable, intent(in) :: hierarchy
            type(lib_tree_universal_index), intent(in) :: uindex
            integer(UINDEX_BYTES) :: rv

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=UINDEX_BYTES) :: n
            integer(kind=UINDEX_BYTES) :: hash
            integer(kind=UINDEX_BYTES) :: max_value
            integer(kind=UINDEX_BYTES) :: element_list_index

            max_value = size(hierarchy(uindex%l)%hashed_element_list_index)
            hash = hash_fnv1a(uindex%n, max_value)
            rv = IGNORE_ENTRY
            do i=1, hierarchy(uindex%l)%maximum_number_of_hash_runs
                if (hierarchy(uindex%l)%hashed_element_list_index(hash)%number_of_hash_runs .gt. 0) then
                    element_list_index = hierarchy(uindex%l)%hashed_element_list_index(hash)%array_position
                    n = hierarchy(uindex%l)%element_list_index(element_list_index)
                    if ((uindex%n .eq. n) .and. &
                        (hierarchy(uindex%l)%hashed_element_list_index(hash)%number_of_hash_runs .eq. i)) then
                        rv = element_list_index
                        exit
                    else
                        hash = ieor(hash, i)
                        hash = hash_fnv1a(hash, max_value)
                    end if
                else
                    rv = IGNORE_ENTRY
                    exit
                end if
            end do

        end function lib_tree_hierarchy_get_index_hierarchy_hashed

        ! Returns a list of unique univeral indices
        !
        ! Argument
        ! ----
        !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
        !       list of universal indices
        !
        ! Returns
        ! ----
        !   rv: array<lib_tree_universal_index>, inout (read-only)
        !       list of unique universal indices
        !
        function lib_tree_hierarchy_make_uindex_list_unique(uindex_list) result(rv)
            implicit none
            ! dummy
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list
            type(lib_tree_universal_index), dimension(:), allocatable :: rv

            ! auxiliary
            type(hash_list_type), dimension(size(uindex_list) * LIB_TREE_HIERARCHY_MARGIN) :: hash_list
            type(lib_tree_universal_index) :: buffer_uindex
            integer(kind=4) :: max_value
            integer(kind=4) :: hash
            integer(kind=UINDEX_BYTES) :: counter
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=2) :: ii

            if (size(uindex_list) .gt. 0) then

                hash_list(:)%value = -1
                hash_list(:)%hash_runs = 0

                max_value = size(hash_list)

                counter = 0
                do i=1, size(uindex_list)
                    ! hash runs
                    do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                        hash = hash_fnv1a(uindex_list(i)%n, max_value)
                        if (hash_list(hash)%hash_runs .eq. 0) then
                            hash_list(hash)%value = uindex_list(i)%n
                            hash_list(hash)%hash_runs = ii
                            counter = counter + 1
                            exit
                        else if (hash_list(hash)%value .ne. uindex_list(i)%n) then
                            hash = IEOR(hash, int(ii, 4))
                            hash = hash_fnv1a(hash, max_value)
                        end if
                    end do
                end do

                allocate (rv(counter))

                counter = 0
                buffer_uindex%l = uindex_list(1)%l
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        buffer_uindex%n = hash_list(i)%value
                        rv(counter) = buffer_uindex
                    end if
                end do

            else
                rv = uindex_list
            end if
        end function lib_tree_hierarchy_make_uindex_list_unique

        ! Returns a list of unique univeral indices
        !
        ! Argument
        ! ----
        !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
        !       list of universal indices
        !
        ! Returns
        ! ----
        !   rv: array<lib_tree_universal_index>, inout (read-only)
        !       list of unique universal indices
        !
        subroutine lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique(uindex_list_X, uindex_list_Y, uindex_list_XY)
            implicit none
            ! dummy
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list_X
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list_Y
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list_XY

            ! auxiliary
            type(hash_list_type), dimension(:), allocatable :: hash_list
            integer(kind=1), dimension(:), allocatable :: hash_type
            type(lib_tree_universal_index) :: buffer_uindex
            integer(kind=4) :: max_value
            integer(kind=4) :: hash
            integer(kind=UINDEX_BYTES) :: counter_X
            integer(kind=UINDEX_BYTES) :: counter_Y
            integer(kind=UINDEX_BYTES) :: counter_XY
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=2) :: ii

            integer(kind=UINDEX_BYTES) :: length

            length = (size(uindex_list_X) + size(uindex_list_Y) + size(uindex_list_XY)) * &
                                                    LIB_TREE_HIERARCHY_MARGIN
            if (length .gt. 0) then
                allocate(hash_list(length))
                allocate(hash_type(length))

                hash_list(:)%value = -1
                hash_list(:)%hash_runs = 0

                max_value = size(hash_list)
                counter_X = 0
                counter_Y = 0
                counter_XY = 0

                ! X-hierarchy
                if (size(uindex_list_X) .gt. 0) then
                    buffer_uindex%l = uindex_list_X(1)%l
                    do i=1, size(uindex_list_X)
                        ! hash runs
                        do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                            hash = hash_fnv1a(uindex_list_X(i)%n, max_value)
                            if (hash_list(hash)%hash_runs .eq. 0) then
                                ! new entry at hash_list
                                hash_list(hash)%value = uindex_list_X(i)%n
                                hash_list(hash)%hash_runs = ii
                                counter_X = counter_X + 1
                                hash_type(hash) = HIERARCHY_X
                                exit
                            else if (hash_list(hash)%value .ne. uindex_list_X(i)%n) then
                                ! hash collision -> re-hash
                                hash = IEOR(hash, int(ii, 4))
                                hash = hash_fnv1a(hash, max_value)
                            else if (hash_list(hash)%value .eq. uindex_list_X(i)%n) then
                                ! uindex duplicate found
                                exit
                            else
                                print *, "lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique: ERROR"
                                print *, "  hash: uindex_list_X"
                                exit
                            end if
                        end do
                    end do
                end if

                ! Y-hierarchy
                if (size(uindex_list_Y) .gt. 0) then
                    buffer_uindex%l = uindex_list_Y(1)%l
                    do i=1, size(uindex_list_Y)
                        ! hash runs
                        do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                            hash = hash_fnv1a(uindex_list_Y(i)%n, max_value)
                            if (hash_list(hash)%hash_runs .eq. 0) then
                                ! new entry at hash_list
                                hash_list(hash)%value = uindex_list_Y(i)%n
                                hash_list(hash)%hash_runs = ii
                                counter_Y = counter_Y + 1
                                hash_type(hash) = HIERARCHY_Y
                                exit
                            else if (hash_list(hash)%value .ne. uindex_list_Y(i)%n) then
                                ! hash colision -> re-hash
                                hash = IEOR(hash, int(ii, 4))
                                hash = hash_fnv1a(hash, max_value)
                            else if (hash_list(hash)%value .eq. uindex_list_Y(i)%n) then
                                if (hash_type(hash) .eq. HIERARCHY_Y) then
                                    ! uindex duplicate found
                                    exit
                                else if (hash_type(hash) .eq. HIERARCHY_X) then
                                    ! uindex was previously of type x-hierarchy
                                    counter_X = counter_X - 1
                                    counter_XY = counter_XY + 1
                                    hash_type(hash) = HIERARCHY_XY
                                    exit
                                else
                                    ! should never be reached
                                    print *, "lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique: ERROR"
                                    print *, "  hash: uindex_list_Y: should never be reached"
                                    exit
                                end if
                            else
                                print *, "lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique: ERROR"
                                print *, "  hash: uindex_list_Y"
                                exit
                            end if
                        end do
                    end do
                end if

                ! XY-hierarchy
                if (size(uindex_list_XY) .gt. 0) then
                    buffer_uindex%l = uindex_list_XY(1)%l
                    do i=1, size(uindex_list_XY)
                        ! hash runs
                        do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                            hash = hash_fnv1a(uindex_list_XY(i)%n, max_value)
                            if (hash_list(hash)%hash_runs .eq. 0) then
                                ! new entry at hash_list
                                hash_list(hash)%value = uindex_list_XY(i)%n
                                hash_list(hash)%hash_runs = ii
                                counter_XY = counter_XY + 1
                                hash_type(hash) = HIERARCHY_XY
                                exit
                            else if (hash_list(hash)%value .ne. uindex_list_XY(i)%n) then
                                ! hash colision -> re-hash
                                hash = IEOR(hash, int(ii, 4))
                                hash = hash_fnv1a(hash, max_value)
                            else if (hash_list(hash)%value .eq. uindex_list_XY(i)%n) then
                                ! entry exists -> set a more general type (HIERARCHY_XY)
                                if (hash_type(hash) .eq. HIERARCHY_XY) then
                                    ! duplicate found
                                    exit
                                else if (hash_type(hash) .eq. HIERARCHY_X) then
                                    ! uindex was previously of type x-hierarchy
                                    counter_X = counter_X - 1
                                    counter_XY = counter_XY + 1
                                    hash_type(hash) = HIERARCHY_XY
                                    exit
                                else if (hash_type(hash) .eq. HIERARCHY_Y) then
                                    ! uindex was previously of type y-hierarchy
                                    counter_Y = counter_Y - 1
                                    counter_XY = counter_XY + 1
                                    hash_type(hash) = HIERARCHY_XY
                                    exit
                                else
                                    ! should never be reached
                                    print *, "lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique: ERROR"
                                    print *, "  hash: uindex_list_XY: should never be reached"
                                    exit
                                end if
                            else
                                print *, "lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique: ERROR"
                                print *, "  hash: uindex_list_XY"
                                exit
                            end if
                        end do
                    end do
                end if

                if ((counter_X + counter_Y + counter_XY) .gt. 0) then
                    if (allocated(uindex_list_X)) then
                        deallocate(uindex_list_X)
                        allocate(uindex_list_X(counter_X))
                    end if

                    if (allocated(uindex_list_Y)) then
                        deallocate(uindex_list_Y)
                        allocate(uindex_list_Y(counter_Y))
                    end if

                    if (allocated(uindex_list_XY)) then
                        deallocate(uindex_list_XY)
                        allocate(uindex_list_XY(counter_XY))
                    end if

                    counter_X = 0
                    counter_Y = 0
                    counter_XY = 0
                    do i=1, size(hash_list)
                        if (hash_list(i)%hash_runs .ne. 0) then
                            if (hash_type(i) .eq. HIERARCHY_X) then
                                counter_X = counter_X + 1
                                buffer_uindex%n = hash_list(i)%value
                                uindex_list_X(counter_X) = buffer_uindex
                            else if (hash_type(i) .eq. HIERARCHY_Y) then
                                counter_Y = counter_Y + 1
                                buffer_uindex%n = hash_list(i)%value
                                uindex_list_Y(counter_Y) = buffer_uindex
                            else if (hash_type(i) .eq. HIERARCHY_XY) then
                                counter_XY = counter_XY + 1
                                buffer_uindex%n = hash_list(i)%value
                                uindex_list_XY(counter_XY) = buffer_uindex
                            end if
                        end if
                    end do
                end if
            end if
        end subroutine lib_tree_hierarchy_make_uindex_list_X_Y_XY_unique

!    ! Argument
!    ! ----
!    !   uindex: type(lib_tree_universal_index)
!    !       universal index of a box
!    !
!    ! Returns
!    ! ----
!    !   C: type(lib_tree_hierarchy_coefficient)
!    !       expansion coefficients of the box uindex
!    function lib_tree_hierarchy_get_C_of_box(uindex) result(C)
!        implicit none
!        ! dummy
!        type(lib_tree_universal_index), intent (in) :: uindex
!        type(lib_tree_hierarchy_coefficient) :: C
!
!        ! auxiliary
!        type(lib_tree_hierarchy_coefficient) :: C_child
!        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_children
!        type(lib_tree_spatial_point) :: x_c
!        type(lib_tree_spatial_point) :: x_c_child
!
!        integer(kind=UINDEX_BYTES) :: i
!        integer(kind=UINDEX_BYTES) :: list_index
!        integer(kind=1) :: hierarchy_type
!
!        type(lib_tree_hierarchy_coefficient), dimension(:), allocatable :: buffer_C
!
!        x_c = lib_tree_get_centre_of_box(uindex)
!
!        ! get child boxes
!        uindex_children = lib_tree_get_children(uindex)
!
!        allocate(buffer_C(size(uindex_children)))
!
!        !$OMP PARALLEL DO PRIVATE(i, x_c_child, list_index, hierarchy_type)
!        do i=1, size(uindex_children)
!            x_c_child = lib_tree_get_centre_of_box(uindex_children(i))
!            list_index = lib_tree_hierarchy_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
!            if (list_index .gt. 0) then
!                hierarchy_type = m_ml_fmm_hierarchy(uindex_children(i)%l)%hierarchy_type(list_index)
!                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
!                    (hierarchy_type .eq. HIERARCHY_XY))) then
!                    C_child = m_ml_fmm_hierarchy(uindex_children(i)%l)%coefficient_list(list_index)
!                    buffer_C(i) = m_ml_fmm_handles%get_translation_SS(C_child, x_c_child, x_c)
!                end if
!            end if
!        end do
!        !$OMP END PARALLEL DO
!
!        call lib_tree_hierarchy_type_operator_set_coefficient_zero(C)
!
!        do i=1, size(uindex_children)
!            list_index = lib_tree_hierarchy_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
!            if (list_index .gt. 0) then
!                hierarchy_type = m_ml_fmm_hierarchy(uindex_children(i)%l)%hierarchy_type(list_index)
!                if (((hierarchy_type .eq. HIERARCHY_X) .or. &
!                    (hierarchy_type .eq. HIERARCHY_XY))) then
!                    C = C + buffer_C(i)
!                end if
!            end if
!        end do
!
!    end function lib_tree_hierarchy_get_C_of_box
!
!    ! Returns
!    ! ----
!    !   rv: integer
!    !       neighbourhood size
!    function lib_tree_hierarchy_get_neighbourhood_size_k() result(rv)
!        implicit none
!        ! dummy
!        integer :: rv
!
!        rv = m_tree_neighbourhood_size_k
!
!    end function lib_tree_hierarchy_get_neighbourhood_size_k

    function lib_tree_hierarchy_test_functions() result(rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        if (.not. test_lib_tree_hierarchy_constructor_without_length()) rv = rv + 1


    contains

        function test_lib_tree_hierarchy_constructor_without_length() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            integer :: n
            type(cartesian_coordinate_real_type) :: point
            type(lib_tree_hierarchy_class) :: h
            type(lib_tree_data_element), dimension(:), allocatable :: data_elements
!            integer(kind=UINDEX_BYTES), dimension(3) :: length ! optional
!            integer, intent(in), optional :: tree_s_opt ! optional
            type(lib_tree_universal_index) :: uindex
            type(lib_tree_hierarchy_box) :: box

            n = 10
            allocate(data_elements(n))

            if (TREE_DIMENSIONS .eq. 3) then
                do i = 1, n
                    point = dble(i - 1) * make_cartesian(1d0, 1d0, 1d0)
                    data_elements(i)%point_x%x = (/point%x, point%y, point%z/)
                    data_elements(i)%expansion = 1d0/i
                end do
                data_elements(n)%expansion = 0
                data_elements(:)%hierarchy = HIERARCHY_X
            else
                print *, "test_lib_tree_hierarchy_constructor_without_length: WARNING"
                print *, "  TREE_DIMENSION .ne. 3: ", TREE_DIMENSIONS
                rv = .true.
                return
            end if

            call h%constructor(data_elements)

            uindex%l = 0
            uindex%n = 0
            box = h%get_box(uindex)

        end function test_lib_tree_hierarchy_constructor_without_length

    end function lib_tree_hierarchy_test_functions

end module lib_tree_hierarchy
