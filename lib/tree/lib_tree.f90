#define _DEBUG_

#define _FMM_DIMENSION_ 2

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
!
! Constraint
! ----
!          |  _FMM_DIMENSION_  |
!   value  |    2     |   3    |
!   -----------------------------
!   single | [4,8]    | [8]    |
!   double | [8,16]   | [8,16] |
!
#define _UINDEX_BYTES_ 8

module lib_tree
    !$  use omp_lib
    use lib_tree_type
    use lib_tree_type_operator
    use lib_tree_helper_functions
    use lib_hash_function
    use lib_sort

    use toolbox
    implicit none
    ! Data Structures, Optimal Choice of Parameters, and Complexity Results for Generalized Multilevel Fast Multipole Methods in d Dimensions

    private

    ! --- public parameter ---
    public :: TREE_DIMENSIONS

    ! --- public functions ---
    public :: lib_tree_constructor
    public :: lib_tree_destructor

    public :: lib_tree_get_element_list
    public :: lib_tree_get_correspondence_vector

    public :: lib_tree_get_parent
    public :: lib_tree_get_children

    public :: lib_tree_get_domain_e1
    public :: lib_tree_get_domain_e2
    public :: lib_tree_get_domain_e3
    public :: lib_tree_get_domain_e4
    public :: lib_tree_get_domain_i2
    public :: lib_tree_get_domain_i3
    public :: lib_tree_get_domain_i4

    public :: lib_tree_get_level_min
    public :: lib_tree_get_level_max

    public :: lib_tree_get_number_of_boxes
    public :: lib_tree_get_centre_of_box

    public :: lib_tree_test_functions
    public :: lib_tree_benchmark

    public :: lib_tree_get_diagonal_test_dataset

    public :: lib_tree_hf_test_functions
    public :: lib_tree_hf_benchmark
!    public :: lib_tree_hf_destructor       ! only for debugging purpose

    ! --- public member variables ---
    public :: lib_tree_scaling_D
    public :: lib_tree_scaling_x_min
    public :: lib_tree_scaling_x_max

    ! --- public type definitions ---
    public :: lib_tree_spatial_point

    ! --- public operator definitions ---
    public :: operator (+)
    public :: operator (-)
    public :: operator (*)

    ! --- member---
    integer(kind=4), parameter :: LIB_TREE_MAX_HASH_RUNS = 400

    integer(kind=1), parameter :: LIB_TREE_ELEMENT_TYPE_EMPTY = -1


    ! --- module global ---
    type(lib_tree_data_element), dimension (:), allocatable :: lib_tree_data_element_list
    ! List of references of data elements at the lib_tree_data_element_list array.
    ! The access is granted by the hashed universal index instead by the universal index directly.
    type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: lib_tree_correspondence_vector
    ! List of data_element positions at the lib_tree_data_element_list..
    ! These references (positions) are sorted into ascending numerical order of the universal index.
    integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: lib_tree_correspondence_vector_sorted_data_elements

    ! threshold level: At this level of the tree the correspondence vector and the data element list are working.
    integer(kind=1) :: lib_tree_l_th
    ! maximum hash value (correspondence with the length of the lib_tree_data_element_list * margin)
    integer(kind=4) :: lib_tree_hash_max
    ! Maximum number of hash runs to access a data element from the lib_tree_correspondence_vector array.
    ! If the hash runs exceed this number, there is no data element with the wanted universal index.
    integer(kind=2) :: lib_tree_max_number_of_hash_runs

    ! scaling
    type(lib_tree_spatial_point) :: lib_tree_scaling_D
    type(lib_tree_spatial_point) :: lib_tree_scaling_x_min
    type(lib_tree_spatial_point) :: lib_tree_scaling_x_max

    contains

    ! Arguments
    ! ----
    !   element_list: lib_tree_data_element
    !
    !   s: integer
    !       maximum number of elements at a box
    !
    subroutine lib_tree_constructor(element_list, s)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), intent(inout) :: element_list
        integer(kind=CORRESPONDENCE_VECTOR_KIND), optional :: s    ! "optimisation value"

        ! auxiliary
!        type(lib_tree_data_element), dimension(size(element_list)) :: element_list_scaled

        ! auxiliary: create_correspondece_vector
        integer(kind=1) :: threshold_level
        integer(kind=2) :: margin

        ! auxiliary: get_level_max
        integer(kind=1) :: l_max

#ifdef _DEBUG_
        ! debug
        integer(kind=1) :: i
        real :: start, finish
        integer(kind=1) :: length
        double precision, dimension(6) :: cpu_time_delta
        character(len=75), dimension(6) :: note_string
#endif

        ! copy data to the module globle variable
        lib_tree_data_element_list = element_list
        ! scaling
#ifdef _DEBUG_
        call cpu_time(start)
#endif
        call lib_tree_get_scaled_element_list(lib_tree_data_element_list)
#ifdef _DEBUG_
        call cpu_time(finish)
        cpu_time_delta(1) = finish - start
        note_string(1) = "lib_tree_get_scaled_element_list"
#endif

        margin = 200
#if (_UINDEX_BYTES_ == 16)
        threshold_level = int(int(UINDEX_BYTES,2) * int(NUMBER_OF_BITS_PER_BYTE,2) * 1.0 / TREE_DIMENSIONS, 1) - 1
#else
        threshold_level = int(UINDEX_BYTES * NUMBER_OF_BITS_PER_BYTE * 1.0 / TREE_DIMENSIONS, 1)
#endif
#ifdef _DEBUG_
        print *, "threshold_level = ", threshold_level
        call cpu_time(start)
#endif
        call lib_tree_create_correspondence_vector(threshold_level, margin)
#ifdef _DEBUG_
        cpu_time_delta(2) = finish - start
        call cpu_time(finish)
        note_string(2) = "lib_tree_create_correspondence_vector"
#endif
        length = 2

        ! optimise threshold level
        if (present(s)) then
#ifdef _DEBUG_
            call cpu_time(start)
#endif
            call lib_tree_create_correspondece_vector_sorted_data_elements()
#ifdef _DEBUG_
            call cpu_time(finish)
            cpu_time_delta(3) = finish - start
            note_string(3) = "lib_tree_create_correspondece_vector_sorted_data_elements"
#endif
#ifdef _DEBUG_
            call cpu_time(start)
#endif
            length = 3
            l_max = lib_tree_get_level_max(s)
            if ( l_max .lt. threshold_level ) then
#ifdef _DEBUG_
                call cpu_time(finish)
                cpu_time_delta(4) = finish - start
                note_string(4) = "lib_tree_get_level_max"
#endif

                ! clean up for optimisation
                if (allocated(lib_tree_correspondence_vector_sorted_data_elements)) then
                    deallocate (lib_tree_correspondence_vector_sorted_data_elements)
                end if
                if (allocated(lib_tree_correspondence_vector)) then
                    deallocate (lib_tree_correspondence_vector)
                end if
#ifdef _DEBUG_
                print *, "threshold_level = ", l_max
                call cpu_time(start)
#endif
                call lib_tree_create_correspondence_vector(l_max, margin)
#ifdef _DEBUG_
                call cpu_time(finish)
                cpu_time_delta(5) = finish - start
                note_string(5) = "lib_tree_create_correspondence_vector"
#endif
#ifdef _DEBUG_
                call cpu_time(start)
#endif
                call lib_tree_create_correspondece_vector_sorted_data_elements()
#ifdef _DEBUG_
                call cpu_time(finish)
                cpu_time_delta(6) = finish - start
                note_string(6) = "lib_tree_create_correspondece_vector_sorted_data_elements"
#endif
                ! todo: optimise correspondence vector (if size(elements) ~~ size(correspondence vector), then without hash  algorithm)
                ! ~~~ end: optimization ~~~
                length = 6
            else
                print *, "lib_tree_constructor: NOTE"
                print *, "  number of elements per box could NOT be set with *s*"
            end if
        end if

#ifdef _DEBUG_
        print *, " ------- DEBUG: lib_tree_constructor (CPU time) -------"
        do i=1, length
            print *, note_string(i), cpu_time_delta(i), " seconds"
        end do
        print *, " ------------------------------------------------------"
#endif

    end subroutine lib_tree_constructor

    ! cleans up the memory
    subroutine lib_tree_destructor()
        !$OMP SINGLE
        if (allocated (lib_tree_correspondence_vector)) then
            deallocate (lib_tree_correspondence_vector)
        end if

        if (allocated (lib_tree_data_element_list)) then
            deallocate (lib_tree_data_element_list)
        end if

        if (allocated (lib_tree_correspondence_vector_sorted_data_elements)) then
            deallocate (lib_tree_correspondence_vector_sorted_data_elements)
        end if

        call lib_tree_hf_destructor()
        !$OMP END SINGLE
    end subroutine lib_tree_destructor

    function lib_tree_get_element_list() result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension (:), allocatable :: rv

        rv = lib_tree_data_element_list
    end function

    function lib_tree_get_correspondence_vector() result(rv)
        implicit none
        ! dummy
        type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: rv

        rv = lib_tree_correspondence_vector
    end function

    !
    ! D = max_d D_d                        (59)
    !
    ! \overline{x} = (x − x_min) / D       (61)
    !
    ! Modification of equation (59) to avoid ones as a result for the
    ! normalised coordinates \overline{x}:
    !
    ! \overline{x} = ((x − x_min) / D) * (1 + 2**(-BITS_MANTISSA))
    !
    ! BITS_MANTISSA
    !   number of bits of the Mantissa (floating point number)
    !
    subroutine lib_tree_get_scaled_element_list(element_list)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), intent(inout) :: element_list
!        type(lib_tree_data_element), dimension(size(element_list)) :: rv

        ! auxiliary
        type(lib_tree_spatial_point) :: D
        type(lib_tree_spatial_point) :: x_min
        type(lib_tree_spatial_point) :: x_max

        integer(kind=COORDINATE_BINARY_BYTES) :: i
        integer(kind=1) :: ii

        x_max%x(:) = 0
        x_min%x(:) = 1

        do ii=1, TREE_DIMENSIONS
            x_max%x(ii) = maxval(element_list(:)%point_x%x(ii))
            x_min%x(ii) = minval(element_list(:)%point_x%x(ii))
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
            D%x(ii) = (x_max%x(ii) - x_min%x(ii)) * (1 + 2D0**(-52))    ! double: BITS_MANTISSA = 52
#else
            D%x(ii) = (x_max%x(ii) - x_min%x(ii)) * (1 + 2.0**(-23))    ! single: BITS_MANTISSA = 23
#endif
!            D%x(ii) = (x_max%x(ii) - x_min%x(ii)) * nearest(1.0, -1.0)
        end do

        ! use for all dimension the same scaling factor
        D%x(:) = maxval(D%x)

        !$OMP PARALLEL DO PRIVATE(i, ii)
        do i=1, size(element_list)
            do ii=1, TREE_DIMENSIONS
!                rv(i)%point_x%x(ii) = (element_list(i)%point_x%x(ii) - x_min%x(ii)) / D%x(ii)
                element_list(i)%point_x%x(ii) = (element_list(i)%point_x%x(ii) - x_min%x(ii)) / D%x(ii)
            end do
        end do
        !$OMP END PARALLEL DO

        ! store for rescaling
        lib_tree_scaling_D = D
        lib_tree_scaling_x_min = x_min
        lib_tree_scaling_x_max = x_max

!        rv(:)%element_type = element_list(:)%element_type

    end subroutine lib_tree_get_scaled_element_list

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !   step: integer, optional (standard value: 1)
    !       value = 1 : function returns the universal index of the parent box
    !       value >= 1: function returns the universal index of the (grand*step)parent box
    !
    ! Returns
    ! ----
    !   the universal index of the parent box
    function lib_tree_get_parent(uindex, step) result (rv)
        implicit none

        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=1), intent (in), optional :: step
        type(lib_tree_universal_index) :: rv

        ! auxiliary
        integer(kind=1) :: m_step

        m_step = 1
        if(present(step))m_step=step

        ! check
        if (uindex%l .lt. m_step) then
            print *, "lib_tree_get_parent: ERROR"
            print *, "  (grand..)parent overflow: step > uindex%l"
            print *, "  step: ", step
            print *, "  uindex%l: ", uindex%l
        end if

        ! calc
        rv%n = lib_tree_hf_get_parent(uindex%n, m_step)
        rv%l = uindex%l - int(m_step,1)

    end function lib_tree_get_parent

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !
    ! Returns
    ! ----
    !   all the univerval indeces of the children's boxes.
    function lib_tree_get_children(uindex) result (rv)
        implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_tree_universal_index), dimension(2**TREE_DIMENSIONS) :: rv

        rv(:)%n = lib_tree_hf_get_children_all(uindex%n)
        rv(:)%l = uindex%l + int(1,1)

    end function lib_tree_get_children

    ! Arguments
    ! ----
    !   k
    !       k-neighbour of the box (n,l)
    !   uindex
    !       universal index of the box
    !
    ! Returns
    ! ----
    !   all the universal indeces of the k-neighours.
    function lib_tree_get_neighbours(k,uindex) result (rv)
        implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=UINDEX_BYTES) :: k
        type(lib_tree_universal_index), dimension(3**TREE_DIMENSIONS-1) :: rv

        rv(:)%n = lib_tree_hf_get_neighbour_all_xD(k, uindex%n, uindex%l)
        rv(:)%l = uindex%l

    end function lib_tree_get_neighbours

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !   element_number: array<integer>
    !       location of the element in the lib_tree_data_element_list
    !
    ! Returns
    ! ----
    !   the spatial points inside the box (n,l)
    !
    !   rv: type(lib_tree_data_element), dimension(:), ALLOCATABLE
    function lib_tree_get_domain_e1(uindex, element_number) result (rv)
        implicit none
        ! dummy arguments
        type(lib_tree_universal_index) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(out) :: element_number
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliar
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_buffer_list
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_1_uindex
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_2_uindex
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_data_element
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_rv
        integer(kind=1) :: l_diff
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: buffer_element_number
        integer(kind=UINDEX_BYTES) :: number_of_data_elements
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii
        integer(kind=UINDEX_BYTES) :: index_start
        integer(kind=UINDEX_BYTES) :: index_end

        ! ---- OMP ----
        ! simple implementation of a semaphore
        ! multiple read, single write
        !$  logical :: semaphore_write_number_of_data_elements
        ! read and write is initially possilble
        !$  semaphore_write_number_of_data_elements = .true.
        ! ~~~~ OMP ~~~~

        if (uindex%n .ne. TREE_BOX_IGNORE_ENTRY) then
            l_diff = lib_tree_l_th - uindex%l

!            allocate (buffer_rv(2**(l_diff*TREE_DIMENSIONS)))
            allocate (buffer_1_uindex(2**(l_diff*TREE_DIMENSIONS)))
            allocate (buffer_2_uindex(2**(l_diff*TREE_DIMENSIONS)))

            ! get all universal indices at the threshold level
            buffer_1_uindex(1) = uindex
            do i=1, l_diff
                !$OMP PARALLEL DO PRIVATE(index_start, index_end, i)
                do ii=1, 2**((i-1)*TREE_DIMENSIONS)
                    index_start = (ii-1)*2**TREE_DIMENSIONS+1
                    index_end = index_start + 2**TREE_DIMENSIONS-1
                    buffer_2_uindex(index_start:index_end) = lib_tree_get_children(buffer_1_uindex(ii))
                end do
                !$OMP END PARALLEL DO
                buffer_1_uindex = buffer_2_uindex
            end do
            deallocate (buffer_2_uindex)

            ! get all elements
            number_of_data_elements = 0
!            allocate(element_number_buffer_list(size(buffer_1_uindex)))
            !$OMP PARALLEL DO PRIVATE(buffer_data_element, buffer_element_number, i)
            do i=1, size(buffer_1_uindex)
                buffer_data_element = lib_tree_get_element_from_correspondence_vector(buffer_1_uindex(i), &
                                                                               buffer_element_number)
                ! ---- OMP semaphore ----
                !$  if (semaphore_write_number_of_data_elements .eqv. .false.) then
                !$      ! wait until write process is finished
                !$      do
                !$          if (semaphore_write_number_of_data_elements .eqv. .true.) then
                !$              exit
                !$          end if
                !$      end do
                !$  end if
                ! ~~~~ OMP semaphore ~~~~

                !$  semaphore_write_number_of_data_elements = .false.
                do ii=1, size(buffer_data_element)
                    if (buffer_data_element(ii)%element_type .ne. LIB_TREE_ELEMENT_TYPE_EMPTY) then
                        call concatenate(buffer_rv, buffer_data_element(ii))                            ! todo: optimize
                        call concatenate(element_number_buffer_list, buffer_element_number(ii))         ! todo: optimize
                        number_of_data_elements = number_of_data_elements + 1
                    end if
                end do
                !$  semaphore_write_number_of_data_elements = .true.

                ! clean up
                if (allocated(buffer_element_number)) then
                    deallocate(buffer_element_number)
                end if
                if (allocated(buffer_data_element)) then
                    deallocate(buffer_data_element)
                end if
            end do
            !$OMP END PARALLEL DO

            deallocate (buffer_1_uindex)

            if (number_of_data_elements .gt. 0) then
                allocate (rv(number_of_data_elements))
                allocate (element_number(number_of_data_elements))

                ii=1
                do i=1, size(buffer_rv)
                    if (buffer_rv(i)%element_type .ne. LIB_TREE_ELEMENT_TYPE_EMPTY) then
                        element_number(ii) = element_number_buffer_list(i)
                        rv(ii) = buffer_rv(i)
                        ii = ii + 1
                    end if
                end do
            else
                allocate(rv(0))
            end if

            if ( allocated(buffer_rv) ) then
                deallocate (buffer_rv)
            end if
        else
            allocate(rv(0))
        end if

    end function lib_tree_get_domain_e1

    ! Arguments
    ! ----
    !   k
    !       k-neighbour of the box (n,l)
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !   element_number: array<integer>
    !       location of the element in the lib_tree_data_element_list
    !
    ! Returns
    ! ----
    !   the spatial points in the k-neighorhood of box (n,l)
    function lib_tree_get_domain_e2(k,uindex, element_number) result (rv)
        implicit none
        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(out) :: element_number
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliar
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        type(lib_tree_universal_index), dimension(:, :), allocatable :: buffer_uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: buffer_element_number
!        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_temp
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_uindex_reduced
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_uindex_reduced_temp
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_data_element_list
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_rv

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii
!        integer(kind=UINDEX_BYTES) :: index_number

        ! get all neighbours
        allocate(buffer_uindex(k, 3**TREE_DIMENSIONS-1))
        !$OMP PARALLEL DO PRIVATE(i)
        do i=1, k
            buffer_uindex(i, :) = lib_tree_get_neighbours(i, uindex)
        end do
        !$OMP END PARALLEL DO


        ! reduce neigbours list (remove invalide universal indices)
        allocate(buffer_uindex_reduced_temp(k*(3**TREE_DIMENSIONS-1) + 1))
        number_of_boxes = 1
        buffer_uindex_reduced_temp(number_of_boxes) = uindex
        do i=1, k
            do ii=1, 3**TREE_DIMENSIONS-1
                if (buffer_uindex(i,ii)%n .ne. TREE_BOX_IGNORE_ENTRY) then
                    number_of_boxes = number_of_boxes + 1
                    buffer_uindex_reduced_temp(number_of_boxes) = buffer_uindex(i,ii)
                end if
            end do
        end do

        ! clean up
        if (allocated(buffer_uindex)) then
            deallocate (buffer_uindex)
        end if

        allocate(buffer_uindex_reduced(number_of_boxes))
        buffer_uindex_reduced = buffer_uindex_reduced_temp(:number_of_boxes)

        ! clean up
        if (allocated(buffer_uindex_reduced_temp)) then
            deallocate (buffer_uindex_reduced_temp)
        end if

        ! get data elements
        do i=1, number_of_boxes
            buffer_data_element_list = lib_tree_get_domain_e1(buffer_uindex_reduced(i), buffer_element_number)
!            if (allocated(buffer_rv)) then
                if ((allocated(buffer_data_element_list)) .and. (size(buffer_data_element_list).gt. 0)) then

                    call concatenate(buffer_rv, buffer_data_element_list)
!                    index_number = size(buffer_rv) + 1
!                    call reallocate(buffer_rv, size(buffer_data_element_list))  ! check replacement
!                    buffer_rv(index_number:) = buffer_data_element_list
                end if

                ! copy buffer_element_number into element_number
                if ((allocated(buffer_element_number)) .and. (size(buffer_element_number).gt. 0)) then

                    call concatenate(element_number, buffer_element_number)
!                    allocate(element_number_temp(size(element_number)))
!                    element_number_temp = element_number
!
!                    ii = size(buffer_element_number) + size(element_number)
!                    deallocate(element_number)
!                    allocate(element_number(ii))
!                    ii = size(element_number_temp)
!                    element_number(:ii) = element_number_temp
!                    element_number(ii+1:) = buffer_element_number
                end if

                if (allocated(buffer_element_number)) then
                    deallocate (buffer_element_number)
                end if
                if (allocated(buffer_data_element_list)) then
                    deallocate (buffer_data_element_list)
                end if
!            else
!                call move_alloc(buffer_element_number, element_number)
!                call move_alloc(buffer_data_element_list, buffer_rv)
!            end if
        end do

        ! clean up
        if (allocated(buffer_uindex_reduced)) then
            deallocate (buffer_uindex_reduced)
        end if

        call move_alloc(buffer_rv, rv)
    end function lib_tree_get_domain_e2

    ! Arguments
    ! ----
    !   k
    !       k-neighbour of the box (n,l)
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !   element_number: array<integer>
    !       location of the element in the lib_tree_data_element_list
    !
    ! Returns
    ! ----
    !   the spatial points outside the k-neighorhood of tree element (n,l)
    function lib_tree_get_domain_e3(k,uindex, element_number) result (rv)
        implicit none
        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(out) :: element_number
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliar
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2
        logical, dimension(:), allocatable :: data_element_list_regard

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii

        ! gets all element numbers that are to be skipped
        buffer_e2 = lib_tree_get_domain_e2(k, uindex, element_number_e2)

        if (allocated(element_number_e2)) then
            ii = size(element_number_e2)
            allocate (element_number(ii))

            call move_alloc(element_number_e2, element_number)

            ! clean up
            deallocate (buffer_e2)
        else
            ii = 0
        end if

        if (ii .gt. 0) then
            ! create a data element list without elements listed in *element_number*
            allocate(data_element_list_regard(size(lib_tree_data_element_list)))
            data_element_list_regard(:) = .true.
            do i=1, size(element_number)
                ii = element_number(i)
                data_element_list_regard(ii) = .false.
            end do

            i = size(lib_tree_data_element_list) - size(element_number)
            allocate (rv(i))
            ii = 1
            do i=1, size(lib_tree_data_element_list)
                if (data_element_list_regard(i)) then
                    rv(ii) = lib_tree_data_element_list(i)
                    ii = ii + 1
                end if
            end do

            ! clean up
            if (allocated(data_element_list_regard)) then
                deallocate (data_element_list_regard)
            end if
        else if (ii .eq. 0) then
            rv = lib_tree_data_element_list
        end if
    end function lib_tree_get_domain_e3

    ! Arguments
    ! ----
    !   k
    !       k-neighbour of the box (n,l)
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !   element_number: array<integer>
    !       location of the element in the lib_tree_data_element_list
    !
    ! Returns
    ! ----
    !   the spatial points in the k-neighorhood of parent box of box (n,l), which do not belong
    !   the k-neighbourhood of the box itself.
    function lib_tree_get_domain_e4(k, uindex, element_number) result (rv)
        implicit none
        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(out) :: element_number
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliary
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_e2
        type(lib_tree_data_element), dimension(:), allocatable :: buffer_e2_parent
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_list_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_list_e2_parent
        logical, dimension(:), allocatable :: data_element_list_regard

        integer(kind=UINDEX_BYTES) :: number_of_elements
        integer(kind=UINDEX_BYTES) :: number_of_elements_parent

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii

        allocate(buffer_e2_parent, source=lib_tree_get_domain_e2(k, lib_tree_get_parent(uindex), element_number_list_e2_parent))
        buffer_e2 = lib_tree_get_domain_e2(k, uindex, element_number_list_e2)

        if (allocated(buffer_e2_parent)) then
            number_of_elements_parent = size(buffer_e2_parent)
            deallocate(buffer_e2_parent)
        else
            number_of_elements_parent = 0
        end if

        if (allocated(buffer_e2)) then
            number_of_elements = size(buffer_e2)
            deallocate(buffer_e2)
        else
            number_of_elements = 0
        end if

        ! set up regarded element list
        allocate (data_element_list_regard(size(lib_tree_data_element_list)))
        data_element_list_regard(:) = .false.

        do i=1, number_of_elements_parent
            ii = element_number_list_e2_parent(i)
            data_element_list_regard(ii) = .true.
        end do

        do i=1, number_of_elements
            ii = element_number_list_e2(i)
            data_element_list_regard(ii) = .false.
        end do

        ! clean up
        if (allocated(element_number_list_e2)) then
            deallocate (element_number_list_e2)
        end if
        if (allocated(element_number_list_e2_parent)) then
            deallocate (element_number_list_e2_parent)
        end if

        ! create a data element list
        i = number_of_elements_parent - number_of_elements
        allocate (rv(i))
        ii = 1
        do i=1, size(lib_tree_data_element_list)
            if (data_element_list_regard(i)) then
                rv(ii) = lib_tree_data_element_list(i)
                ii = ii + 1
            end if
        end do

        ! clean up
        if (allocated(data_element_list_regard)) then
            deallocate (data_element_list_regard)
        end if
    end function lib_tree_get_domain_e4

    function lib_tree_get_domain_i2(uindex, k) result(uindex_i2)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i2

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        type(lib_tree_universal_index), dimension(:, :), allocatable :: buffer_uindex
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_uindex_reduced_temp

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii

        ! get all neighbours
        allocate(buffer_uindex(k, 3**TREE_DIMENSIONS-1))
        !$OMP PARALLEL DO PRIVATE(i)
        do i=1, k
            buffer_uindex(i, :) = lib_tree_get_neighbours(i, uindex)
        end do
        !$OMP END PARALLEL DO


        ! reduce neigbours list (remove invalide universal indices)
        allocate(buffer_uindex_reduced_temp(k*(3**TREE_DIMENSIONS-1) + 1))
        number_of_boxes = 1
        buffer_uindex_reduced_temp(number_of_boxes) = uindex
        do i=1, k
            do ii=1, 3**TREE_DIMENSIONS-1
                if (buffer_uindex(i,ii)%n .ne. TREE_BOX_IGNORE_ENTRY) then
                    number_of_boxes = number_of_boxes + 1
                    buffer_uindex_reduced_temp(number_of_boxes) = buffer_uindex(i,ii)
                end if
            end do
        end do

        ! clean up
        if (allocated(buffer_uindex)) then
            deallocate (buffer_uindex)
        end if

        allocate(uindex_i2(number_of_boxes))
        uindex_i2 = buffer_uindex_reduced_temp(:number_of_boxes)

        ! clean up
        if (allocated(buffer_uindex_reduced_temp)) then
            deallocate (buffer_uindex_reduced_temp)
        end if
    end function lib_tree_get_domain_i2

    function lib_tree_get_domain_i3(uindex, k) result(uindex_i3)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i3

        ! auxiliary
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i2
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: buffer_uindex

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes

        uindex_i2 = lib_tree_get_domain_i2(uindex, k)

        allocate(buffer_uindex(2**(TREE_DIMENSIONS * uindex%l)))

        ! initialize
        do i=1, size(buffer_uindex)
            buffer_uindex(i) = i - int(1, UINDEX_BYTES)
        end do

        ! ignore boxes of the domain i2
        do i=1, size(uindex_i2)
            buffer_uindex(uindex_i2(i)%n + 1) = TREE_BOX_IGNORE_ENTRY
        end do
        ! ignore uindex box
        buffer_uindex(uindex%n + 1) = TREE_BOX_IGNORE_ENTRY

        ! reduce the array buffer_uindex
        number_of_boxes = size(buffer_uindex) - size(uindex_i2) - 1
        if (number_of_boxes .gt. 0) then
            allocate (uindex_i3(number_of_boxes))
            number_of_boxes = 0
            do i=1, size(buffer_uindex)
                if (buffer_uindex(i) .ne. TREE_BOX_IGNORE_ENTRY) then
                    number_of_boxes = number_of_boxes + 1
                    uindex_i3(number_of_boxes)%n = buffer_uindex(i)
                    uindex_i3(number_of_boxes)%l = uindex%l
                end if
            end do
        end if
    end function lib_tree_get_domain_i3

    function lib_tree_get_domain_i4(uindex, k) result(uindex_i4)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i4

        ! parameter
        integer(kind=1), parameter :: HASH_MARGIN = 2 ! = 200%
        integer(kind=2), parameter :: HASH_MAX_RUNS = 200

        ! auxiliary
        type(lib_tree_universal_index) :: uindex_parent
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i2
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_parent_i2
        type(lib_tree_universal_index), dimension(:,:), allocatable :: uindex_parent_i2_children

        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: list_n
        integer(kind=2), dimension(:), allocatable :: number_of_hash_runs
        integer(kind=2) :: hash

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii
        integer(kind=2) :: hash_counter
        integer(kind=UINDEX_BYTES) :: buffer


        uindex_parent = lib_tree_get_parent(uindex)

        allocate(uindex_i2, source=lib_tree_get_domain_i2(uindex, k))
        allocate(uindex_parent_i2, source=lib_tree_get_domain_i2(uindex_parent, k))
        allocate(uindex_parent_i2_children(size(uindex_parent_i2), 2**TREE_DIMENSIONS))


        do i=1, size(uindex_parent_i2)
            uindex_parent_i2_children(i,:) = lib_tree_get_children(uindex_parent_i2(i))
        end do

        allocate(uindex_i4(size(uindex_parent_i2_children) - size(uindex_i2)))

        ! calculate the relative complement
        if (TREE_DIMENSIONS * uindex%l .gt. (log(HASH_MARGIN * real(size(uindex_parent_i2_children))) / log(2.0))) then
            ! use hash table
            allocate(list_n(HASH_MARGIN*size(uindex_parent_i2_children)), source=TREE_BOX_IGNORE_ENTRY)
            allocate(number_of_hash_runs(size(list_n)), source=0_2)

            do i=1, size(uindex_parent_i2)
                do ii=1, 2**TREE_DIMENSIONS
                    hash = int(hash_fnv1a(uindex_parent_i2_children(i,ii)%n, size(list_n)), 2)
                    ! get unique hash
                    do hash_counter=1, HASH_MAX_RUNS
                        if (number_of_hash_runs(hash) .eq. 0) then
                            ! add uindex to the list_n
                            list_n(hash) = uindex_parent_i2_children(i,ii)%n
                            number_of_hash_runs(hash) = hash_counter
                            exit
                        else if (list_n(hash) .ne. uindex_parent_i2_children(i,ii)%n) then
                            ! hash collision -> re-hash
                            hash = ieor(hash, hash_counter)
                            hash = int(hash_fnv1a(hash, size(list_n)), 2)
                        else
                            ! this should never happen
                            print *, "lib_tree_get_domain_i4: ERROR"
                            print *, "  uindex_parent_i2_children: the same uindex occurs several times. n = ", &
                                      uindex_parent_i2_children(i,ii)%n
                            exit
                        end if
                    end do
                end do
            end do

            do i=1, size(uindex_i2)
                hash = int(hash_fnv1a(uindex_i2(i)%n, size(list_n)), 2)
                ! get unique hash
                do hash_counter=1, HASH_MAX_RUNS
                    if (number_of_hash_runs(hash) .eq. 0) then
                        ! add uindex to the list_n: list_n(hash) = -n-1
                        list_n(hash) = -uindex_i2(i)%n - 1
                        number_of_hash_runs(hash) = hash_counter
                        exit
                    else if ((list_n(hash) .eq. uindex_i2(i)%n) .and. &
                             (number_of_hash_runs(hash) .eq. hash_counter)) then
                        ! has to be removed from the list_n -> list_n(hash) = -list_n(hash)-1
                        list_n(hash) = -list_n(hash)-1
                        exit
                    else if ((list_n(hash) .eq. (-uindex_i2(i)%n-1)) .and. &
                             (number_of_hash_runs(hash) .eq. hash_counter)) then
                        print *, "lib_tree_get_domain_i4: ERROR"
                        print *, "  uindex_i2: the same uindex occurs several times. n = ", uindex_i2(i)%n
                    else
                        ! hash collision
                        hash = ieor(hash, hash_counter)
                        hash = int(hash_fnv1a(hash, size(list_n)), 2)
                    end if
                end do
            end do
        else
            ! use lookup table
            allocate(list_n(2**(TREE_DIMENSIONS * uindex%l)))
            list_n(:) = TREE_BOX_IGNORE_ENTRY

            do i=1, size(uindex_parent_i2)
                do ii=1, 2**TREE_DIMENSIONS
                    buffer = uindex_parent_i2_children(i, ii)%n
                    list_n(buffer + 1) = buffer
                end do
            end do

            do i=1, size(uindex_i2)
                buffer = uindex_i2(i)%n
                list_n(buffer + 1) = TREE_BOX_IGNORE_ENTRY
            end do
!
!            ! create uindex_i4
!            buffer = 0
!            do i=1, size(list)
!                if (list(i) .ne. TREE_BOX_IGNORE_ENTRY) then
!                    buffer = buffer + 1
!                    uindex_i4(buffer)%n = list(i)
!                    uindex_i4(buffer)%l = uindex%l
!                end if
!            end do
!
!            ! clean up
!            if allocated(list) then
!                deallocate(list)
!            end if
        end if

        ! create uindex_i4
        buffer = 0
        do i=1, size(list_n)
            if (list_n(i) .ge. 0) then
                buffer = buffer + 1
                uindex_i4(buffer)%n = list_n(i)
                uindex_i4(buffer)%l = uindex%l
            end if
        end do

        ! clean up
        if (allocated(list_n)) then
            deallocate(list_n)
        end if

        if (allocated(number_of_hash_runs)) then
            deallocate(number_of_hash_runs)
        end if
    end function lib_tree_get_domain_i4

    ! Calculates the minimum tree level
    !
    ! Equation
    !   l_min = [1 + log_2(k + 1)]    (21)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Arguments
    ! ----
    !   k: integer
    !       number of neighbours
    !
    ! Returns
    ! ----
    !   l_min: integer
    !       minimum tree level
    !
    function lib_tree_get_level_min(k) result(l_min)
        implicit none
        ! dummy
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        integer(kind=1) :: l_min

        l_min = int(1 + floor(log(real(k+1)) / log(2.0), UINDEX_BYTES), 1)

    end function

    ! Calculates the maximum tree level, where a box contains a maximum of s elements
    !
    ! Arguments
    ! ----
    !   s: integer
    !       number of data elements with the same universal index at *l_max*
    !
    ! Returns
    ! ----
    !   l_max: integer
    !       maximum tree level
    !
    ! Dependency
    ! ----
    !    this function has to be called first: lib_tree_create_correspondece_vector_sorted_data_elements()
    !
    function lib_tree_get_level_max(s) result(l_max)
        implicit none
        ! dummy
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: s
        integer(kind=1) :: l_max

        ! auxiliary
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: i
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: m
        integer(kind=1) :: Bit_max, j
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: N

        type(lib_tree_universal_index) :: a
        type(lib_tree_universal_index) :: b

        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: buffer_index

        if (allocated(lib_tree_correspondence_vector_sorted_data_elements)) then
            N = size(lib_tree_data_element_list)
            l_max = 1
            Bit_max = lib_tree_l_th

            i = 0
            m = s
            do
                i = i + 1
                m = m + 1

                if (m > N) then
                    exit
                end if

                !            a = Interleaved(v(ind(i))
                !            b = Interleaved(v(ind(m))
                buffer_index = lib_tree_correspondence_vector_sorted_data_elements(i)
                a = lib_tree_data_element_list(buffer_index)%uindex

                buffer_index = lib_tree_correspondence_vector_sorted_data_elements(m)
                b = lib_tree_data_element_list(buffer_index)%uindex

                j = Bit_max + int(1,1)

                do
                    j = j - int(1,1)
                    a = lib_tree_get_parent(a)
                    b = lib_tree_get_parent(b)
                    if (a%n .eq. b%n) then
                        l_max = max(l_max , j)
                        exit
                    end if
                end do
            end do
        end if

    end function lib_tree_get_level_max

    ! Returns the threshold level of the Tree.
    ! Each datapoint has a unique universal index at this level.
    function lib_tree_get_level_threshold() result(rv)
        implicit none
        integer(kind=1) :: rv

        rv = lib_tree_l_th

    end function lib_tree_get_level_threshold

    ! Returns the number of boxes at level *l*.
    !
    ! Arguments
    ! ----
    !   l: integer
    !       number of the level [0..l_th]
    !
    ! Returns
    ! ----
    !   rv: integer
    !       number of boxes a t level *l*
    !
    function lib_tree_get_number_of_boxes(l) result (rv)
        implicit none
        ! dummy
        integer(kind=1) :: l
        integer(kind=UINDEX_BYTES) :: rv

        rv = 2**(l*TREE_DIMENSIONS)
    end function lib_tree_get_number_of_boxes

    ! Calculates the centre of a box(n,l).
    !
    !   x_m(n, l) = 2^(−l) ( n + 2^(−1))       (87)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Arguments
    ! ----
    !   n: integer
    !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !   l: integer
    !       number of the level
    ! Returns
    ! ----
    !   x_m: [float, double]
    !       counting system undependet value
    !
    function lib_tree_get_centre_of_box(uindex) result (rv)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        type(lib_tree_spatial_point) :: rv

        rv = lib_tree_hf_get_centre_of_box(uindex%n, uindex%l)

    end function lib_tree_get_centre_of_box

    ! This routine stores references (list entries) of the lib_tree_data_element_list array.
    ! These references are sorted into ascending numerical order of the universal index.
    subroutine lib_tree_create_correspondece_vector_sorted_data_elements()
        implicit none

        ! auxiliray
        integer(kind=UINDEX_BYTES), dimension(:), allocatable :: uindex_list
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: uindex_old_position_list
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: i

        if (allocated(lib_tree_data_element_list) .and. allocated(lib_tree_correspondence_vector)) then
            i = size(lib_tree_data_element_list)
            if (allocated(lib_tree_correspondence_vector_sorted_data_elements)) then
                deallocate(lib_tree_correspondence_vector_sorted_data_elements)
            end if
            allocate (lib_tree_correspondence_vector_sorted_data_elements(i))

            allocate (uindex_list(i))
            allocate (uindex_old_position_list(i))

            ! get uindex list
            uindex_list(:) = lib_tree_data_element_list(:)%uindex%n

            call lib_sort_hpsort_integer(i, uindex_list, uindex_old_position_list)

            call move_alloc(uindex_old_position_list, lib_tree_correspondence_vector_sorted_data_elements)
            ! clean up
            deallocate (uindex_list)
        end if

    end subroutine lib_tree_create_correspondece_vector_sorted_data_elements

    ! Sorting data and saving in the module globally variable *lib_tree_data*
    !
    ! Arguments
    ! ----
    !   element_list
    !       list of data elements
    !
    !   threshold_level
    !       Threshold value (l_max) at which two adjacent elements can still be distinguished.
    !
    !   margin
    !       length of the correspondece vector compared with the length of the element_list
    !       recommented values: 110-200
    !
    ! Returns
    ! ----
    !       the reduced correspondece vector
    !
    !
    ! Correspondence vector
    ! ----
    !   The i-th elements of this vector represents the data point of the i-th box.
    !
    !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C: equation (97) and (98)
    !
    !   element_list
    !       d=2, l_th=2
    !       -----------------
    !       | 3 | 5 | 7 | 1 |   // number = universal index
    !       -----------------
    !
    !   complete correspondence vector (can be sparse)
    !       -> dimension = 2**(d*l_th) = 2**(2*2) = 16
    !       -----------------------------------------------------------------
    !       |   | 3 |   | 0 |   | 1 |   | 2 |   |   |   |   |   |   |   |   |   // number = element_list index
    !       -----------------------------------------------------------------
    !         ^                           ^                               ^
    !  index: 0                           7                              16
    !
    !        |  1. calculate the hash of the universal index
    !        |  2. save elements at the position corresponding to the hash value
    !        V
    !   correspondence vector
    !       -------------------------
    !       |   |   |   |   |   |   |
    !       -------------------------
    !         ^                   ^
    !  index: 0                   5     // hash(universal index)
    !
    subroutine lib_tree_create_correspondence_vector(threshold_level, margin)
        implicit none
        ! dummy
        integer(kind=1), intent(in) :: threshold_level
        integer(kind=2), intent(in) :: margin

        ! auxiliary
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: correspondence_vector_dimension
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: correspondence_vector_dimension_log_2
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hashed_uindex
#ifdef _DEBUG_
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hashed_uindex_old = -1
#endif
        type(lib_tree_universal_index) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: i
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: buffer_element_number
        logical :: element_with_uindex_exists
        integer(kind=2) :: ii
#if (_UINDEX_BYTES_ == 16)
        integer(kind=2) :: number_of_bits
#else
        integer(kind=1) :: number_of_bits
#endif
        integer(kind=4) :: hash_overflow_counter
        logical :: element_saved

        ! ---- OMP ----
        ! simple implementation of a semaphore
        ! multiple read, single write
        !$  logical :: semaphore_write_correspondence_vector

        ! read and write is initially possilble
        !$  semaphore_write_correspondence_vector = .true.
        ! ~~~~ OMP ~~~~


        ! copy element list to the module global *lib_tree_data_element_list* list.

        lib_tree_l_th = threshold_level

        ! check if the significant bits doesn't exceed the number of bits of the universal index
        number_of_bits = threshold_level * TREE_DIMENSIONS
#if (_UINDEX_BYTES_ == 16)
        if (number_of_bits .le. (int(UINDEX_BYTES,2) * int(NUMBER_OF_BITS_PER_BYTE,2))) then
#else
        if (number_of_bits .le. (UINDEX_BYTES * NUMBER_OF_BITS_PER_BYTE)) then
#endif
            correspondence_vector_dimension = ceiling(size(lib_tree_data_element_list) * margin / 100.0)
            correspondence_vector_dimension_log_2 = int(log(real(correspondence_vector_dimension))/log(2.0), &
                                                              CORRESPONDENCE_VECTOR_KIND)
            ! if the number of boxes on the level l_th is smaler than the number of elements, the
            ! size of the correspondece_vector can be set to this number.
            !
            !   correspondence_vector_dimension > 2**(dimension * level) * margin/100.0
            if (correspondence_vector_dimension_log_2 .gt. (TREE_DIMENSIONS*threshold_level + log(margin/100.0) / log(2.0))) then
                correspondence_vector_dimension = ceiling(2**(TREE_DIMENSIONS*threshold_level) * margin/100.0)
            end if

            lib_tree_hash_max = correspondence_vector_dimension

            ! initiate lib_tree_correspondence_vector
            if (.not. allocated (lib_tree_correspondence_vector)) then
                allocate( lib_tree_correspondence_vector(correspondence_vector_dimension) )
            else
                if (size(lib_tree_correspondence_vector) .ne. correspondence_vector_dimension) then
                    deallocate (lib_tree_correspondence_vector)
                    allocate( lib_tree_correspondence_vector(correspondence_vector_dimension) )
                end if
            end if
            lib_tree_correspondence_vector(:)%number_of_hash_runs = 0

            ! iterate element_list
            lib_tree_max_number_of_hash_runs = 0
            hash_overflow_counter = 0
            !$OMP PARALLEL DO PRIVATE(uindex, hashed_uindex, element_saved, i, ii)
            do i=1, size(lib_tree_data_element_list)
                uindex = lib_tree_hf_get_universal_index(lib_tree_data_element_list(i)%point_x, threshold_level)
                lib_tree_data_element_list(i)%uindex = uindex
                ! find unique hashed universal index
                hashed_uindex = hash_fnv1a(uindex%n, lib_tree_hash_max)

                element_saved = .false.
                do ii=1, LIB_TREE_MAX_HASH_RUNS !huge(lib_tree_correspondence_vector(1)%number_of_hash_runs)-1
                    ! ---- OMP semaphore ----
                    !$  if (semaphore_write_correspondence_vector .eqv. .false.) then
                    !$      ! wait until write process is finished
                    !$      do
                    !$          if (semaphore_write_correspondence_vector .eqv. .true.) then
                    !$              exit
                    !$          end if
                    !$      end do
                    !$  end if
                    ! ~~~~ OMP semaphore ~~~~
                    ! save
                    element_with_uindex_exists = .false.
                    if ( allocated(lib_tree_correspondence_vector(hashed_uindex)%data_element_number) ) then
                        buffer_element_number = lib_tree_correspondence_vector(hashed_uindex)%data_element_number(1)
                        if ( lib_tree_data_element_list(buffer_element_number)%uindex%n .eq. uindex%n ) then
                            element_with_uindex_exists = .true.
                        end if
                    end if
                    if ((lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. 0) .or. &
                        element_with_uindex_exists) then
                        !$  semaphore_write_correspondence_vector = .false.
                        if ((hashed_uindex .gt. 0) .and. (hashed_uindex .le. lib_tree_hash_max)) then
                            call concatenate(lib_tree_correspondence_vector(hashed_uindex)%data_element_number, &
                                                                            i)
!                            lib_tree_correspondence_vector(hashed_uindex)%data_element_number = i
                            lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs = ii

                            if (lib_tree_max_number_of_hash_runs .lt. ii) then
                                lib_tree_max_number_of_hash_runs = ii
                            end if
                            ! unique hash found -> terminate the inner do loop immediatly
                            !$  semaphore_write_correspondence_vector = .true.
                            element_saved = .true.
                            exit
                        else
                            print *, "lib_tree_create_correspondence_vector ..ERROR"
                            print *, "  hashed_uindex out of range"
                        end if
                        !$  semaphore_write_correspondence_vector = .true.
                    end if

#ifdef _DEBUG_
                    if (hashed_uindex_old .eq. hashed_uindex) then
                        print *, "lib_tree_create_correspondence_vector .. WARNING"
                        print *, "  hashed uindex doesn't changed"
                        print *, "  run: ", ii, " data element: ", i, " hashed_uindex", hashed_uindex
                    end if

                    hashed_uindex_old = hashed_uindex
#endif

                    hashed_uindex =  IEOR(hashed_uindex, int(ii, CORRESPONDENCE_VECTOR_KIND))
                    hashed_uindex = hash_fnv1a(hashed_uindex, lib_tree_hash_max)
                end do
                if (.not. element_saved) then
                    hash_overflow_counter = hash_overflow_counter + 1
                end if
            end do
            !$OMP END PARALLEL DO
            if (hash_overflow_counter .ne. 0) then
                print *, "lib_tree_create_correspondence_vector  ..ERROR"
                print *, "    number of hash overflows : ", hash_overflow_counter
            end if
        else
            print *, "lib_tree_create_correspondence_vector: tree is too deep"
            print *, "   the biggest universal index would exceed (tree_dimension * threshold level * number_of_bits_per_byte)"
            print *, "   number_of_bits = ", number_of_bits
        end if

    end subroutine lib_tree_create_correspondence_vector

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level, has to be equal to the threshold level (*lib_tree_l_th*)
    !
    !   element_number: integer
    !       location of the element in the lib_tree_data_element_list
    !
    function lib_tree_get_element_from_correspondence_vector(uindex, element_number) result(rv)
        implicit none
        !dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(inout) :: element_number
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliary
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hashed_uindex
        integer(kind=2) :: i
        integer(kind=UINDEX_BYTES) :: ii
        logical :: element_found

        !$  logical :: opm_end_do_loop

        if (uindex%l .ne. lib_tree_l_th) then
            print *, "lib_tree_get_element_from_correspondence_vector:"
            print *, "    level is NOT equal to the threshold level"
            print *, "    l: ", uindex%l
            print *, "    l_th: ", lib_tree_l_th

            element_number = LIB_TREE_ELEMENT_TYPE_EMPTY
            rv%element_type = LIB_TREE_ELEMENT_TYPE_EMPTY
            return
        end if

        if (allocated(lib_tree_correspondence_vector)) then

            hashed_uindex = hash_fnv1a(uindex%n, lib_tree_hash_max)

            element_found = .false.
!            !$  opm_end_do_loop = .false.
!            !  $   OMP PARALLEL DO PRIVATE(i, hash_idum)
            do i=1, lib_tree_max_number_of_hash_runs
!                !$  if (.not. opm_end_do_loop) then
                if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. i) then
!                    !$  opm_end_do_loop = .true.
!                    print *, "element found"

                    element_number = lib_tree_correspondence_vector(hashed_uindex)%data_element_number
                    if (lib_tree_data_element_list(element_number(1))%uindex%n .eq. uindex%n) then
                        element_found = .true.
                        allocate( rv(size(element_number)) )

                        do ii=1, size(element_number)
                            rv(ii) = lib_tree_data_element_list(element_number(ii))
                        end do

                        exit
                    else
                        hashed_uindex =  IEOR(hashed_uindex, int(i, CORRESPONDENCE_VECTOR_KIND))
                        hashed_uindex = hash_fnv1a(hashed_uindex, lib_tree_hash_max)
                    end if
                else
                        hashed_uindex =  IEOR(hashed_uindex, int(i, CORRESPONDENCE_VECTOR_KIND))
                        hashed_uindex = hash_fnv1a(hashed_uindex, lib_tree_hash_max)
                end if
!                !$  end if

            end do
!            !   $   OMP END PARALLEL DO
            if (.not. element_found) then
                if (.not. allocated(rv)) then
                    allocate(rv(1))
                end if
                rv%element_type = LIB_TREE_ELEMENT_TYPE_EMPTY
#ifdef DEBUG
                print *, "Element could not be found: n=", uindex%n," ..note"
#endif
            end if

        end if


    end function lib_tree_get_element_from_correspondence_vector

    ! ----------------- toolbox -----------------

    ! reallocates the data element list with additional n elements
    !
    ! Arguments
    ! ----
    !   a: 1-dimensional lib_tree_data_element array
    !       original array, will be replaced with the resized array
    !   n: pos. integer
    !       number of additional array elements
    !
    ! copy of toolbox.reallocate_1d
    subroutine lib_tree_reallocate_1d_data_element_list(a,n)
        implicit none
        type(lib_tree_data_element),dimension(:),allocatable,intent(inout) :: a
        type(lib_tree_data_element),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        ni_old = size(a)

        allocate(temp(ni_old+n))

        temp(1:ni_old) = a

        call move_alloc(temp,a)

    end subroutine lib_tree_reallocate_1d_data_element_list

    function lib_tree_get_diagonal_test_dataset(list_length, element_type, hierarchy_type, mirror) result(element_list)
        implicit none
        ! dummy
        integer(kind=UINDEX_BYTES), intent(in) :: list_length
        integer(kind=1), intent(in) :: element_type
        integer(kind=1), intent(in) :: hierarchy_type
!        integer(kind=2), optional :: margin
        logical, optional :: mirror
        type(lib_tree_data_element), dimension(list_length) :: element_list

        ! auxiliary
        logical :: m_mirror
!        integer(kind=2) :: m_margin
        integer(kind=UINDEX_BYTES) :: i

        m_mirror = .false.
        if (present(mirror)) m_mirror = mirror

!        m_margin = 300
!        if (present(margin)) m_margin = margin

#if (_FMM_DIMENSION_ == 2)
        if (m_mirror) then
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = ((1.0 - 0.999) * i)/(1.0*list_length)
            end do
        else
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
        end if
#elif (_FMM_DIMENSION_ == 3)
        if (m_mirror) then
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = ((1.0 - 0.9) * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = ((1.0 - 0.9) * i)/(1.0*list_length)
            end do
        else
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
        end if
#endif
        element_list(:)%element_type = element_type
        element_list(:)%hierarchy = hierarchy_type
    end function

    ! ----------------- test functions -----------------
    function lib_tree_test_functions() result(error_counter)
        implicit none

        integer :: error_counter

        !$  print *, "number of threads: ", omp_get_num_threads()
        !$  print *, "max number of threads: ", omp_get_max_threads()

        error_counter = 0

        if (.not. test_lib_tree_create_correspondence_vector()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_element_from_correspondence_vector()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_create_correspondece_vector_sorted_data_elements()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_domain_e1()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_domain_e2()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_domain_e3()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_domain_e4()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_domain_i2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_domain_i4()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_level_min()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_level_max()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_number_of_boxes()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_get_scaled_element_list()) then
            error_counter = error_counter + 1
        end if
        call lib_tree_destructor()
        if (.not. test_lib_tree_constructor()) then
            error_counter = error_counter + 1
        end if

        print *, "-------------lib_tree_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_tree_test_functions tests: OK"
        else
            print *, error_counter,"lib_tree_test_functions test(s) FAILED"
        end if
        print *, "----------------------------------------------------"

    contains

        function test_lib_tree_constructor() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=8), parameter :: list_length = 10**6
            type(lib_tree_data_element), dimension(list_length) :: element_list

            integer(kind=8) :: i

            ! set up the environment
            call lib_tree_destructor()
#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0D0 - (0.999D0 * i)
                element_list(i)%point_x%x(2) = 1.0D0 - (0.999D0 * i)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0D0 - (0.999D0 * i)
                element_list(i)%point_x%x(2) = 1.0D0 - (0.999D0 * i)
                element_list(i)%point_x%x(3) = 1.0D0 - (0.999D0 * i)
            end do
#endif

            call lib_tree_constructor(element_list)

            rv = .true.
            print *, "test_lib_tree_constructor: OK, but there are no tests"
        end function test_lib_tree_constructor

        function test_lib_tree_get_scaled_element_list() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=4), parameter :: list_length = 10**5

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
!            type(lib_tree_data_element), dimension(list_length) :: element_list_rv

            integer(kind=4) :: i

            ! create test data
#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = i
                element_list(i)%point_x%x(2) = i
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = i
                element_list(i)%point_x%x(2) = i
                element_list(i)%point_x%x(3) = i
            end do
#endif
            element_list(:)%element_type = 1

            ! normalise element data point

            call lib_tree_get_scaled_element_list(element_list)

            print *, "test_lib_tree_get_scaled_element_list"
            print *, "  last x-value: ", element_list(list_length)%point_x%x(1)

            ! todo: evaluate element_list_rv

            rv = .true.

        end function test_lib_tree_get_scaled_element_list

        function test_lib_tree_get_domain_e1() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10

            integer(kind=1), parameter :: l_th = 5 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            type(lib_tree_data_element), dimension(:), allocatable :: ground_truth_element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            type(lib_tree_universal_index) :: uindex
            type(lib_tree_data_element), dimension(:), allocatable :: domain_element_list
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number

            integer :: ground_truth_number_of_elements
            integer :: number_of_elements

            ! ---- create correspondece vector ----

            margin = 300


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
            allocate(ground_truth_element_list(2))
            ground_truth_element_list(1) = element_list(1)
            ground_truth_element_list(2) = element_list(2)
            ground_truth_number_of_elements = 2
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
            allocate(ground_truth_element_list(2))
            ground_truth_element_list(1) = element_list(1)
            ground_truth_element_list(2) = element_list(2)
            ground_truth_number_of_elements = 2
#endif
            element_list(:)%element_type = 1

            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .ne. list_length) then
                rv = .false.
                print *, "test_lib_tree_get_domain_e1_create_correspondence_vector: FAILD"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"

                return
            end if

            ! ---- test get_domain_e1 ----

            uindex%n = 0
            uindex%l = 2

            domain_element_list = lib_tree_get_domain_e1(uindex, element_number)

            print *, "test_lib_tree_get_domain_e1:"
            rv = .true.
            if (allocated(domain_element_list)) then
                number_of_elements = size(domain_element_list)

                if (number_of_elements .eq. ground_truth_number_of_elements) then
                    do i=1, size(domain_element_list)
                        if ((ground_truth_element_list(i)%point_x%x(1) .eq. domain_element_list(i)%point_x%x(1)) .and. &
                            (ground_truth_element_list(i)%point_x%x(2) .eq. domain_element_list(i)%point_x%x(2))) then
                            rv = .true.
                            print *, "   ", i, ": OK"
                        else
                            rv = .false.
                            print *, "   ", i, ": FAILED"
                        end if
                    end do
                else
                    rv = .false.
                    print *, "   number of domain elements is NOT equal with the ground truth"
                    print *, "   number of domain elements: ", size(domain_element_list)
                    print *, "   number of ground truth elements: ", size(ground_truth_element_list)
                end if
            else
                rv = .false.
                print *, "   number of domain elements is NOT equal with the ground truth"
                print *, "   number of domain elements: ", 0
                print *, "   number of ground truth elements: ", size(ground_truth_element_list)
            end if

            ! clean up
            if (allocated(ground_truth_element_list)) then
                deallocate (ground_truth_element_list)
            end if
            if (allocated(domain_element_list)) then
                deallocate (domain_element_list)
            end if
            if (allocated(element_number)) then
                deallocate (element_number)
            end if

        end function test_lib_tree_get_domain_e1

        function test_lib_tree_get_domain_e2() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**1

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            type(lib_tree_universal_index) :: uindex
            integer(kind=UINDEX_BYTES) :: k
            type(lib_tree_data_element), dimension(:), allocatable :: domain_element_list
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number

            type(lib_tree_data_element), dimension(:), allocatable :: ground_truth_element_list
            integer :: ground_truth_number_of_elements
            integer :: number_of_elements
            ! ---- create correspondece vector ----

            margin = 300


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, 5
                ground_truth_element_list(i) = element_list(i)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, 5
                ground_truth_element_list(i) = element_list(i)
            end do
#endif
            element_list(:)%element_type = 1

            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .ne. list_length) then
                rv = .false.
                print *, "test_lib_tree_get_domain_e2_create_correspondence_vector: FAILD"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"

                return
            end if

            ! ---- test get_domain_e2 ----
            uindex%n = 0
            uindex%l = 2
            k = 1

            domain_element_list = lib_tree_get_domain_e2(k,uindex, element_number)

            print *, "test_lib_tree_get_domain_e2:"
            rv = .true.
            if (allocated(domain_element_list)) then
                number_of_elements = size(domain_element_list)

                if (number_of_elements .eq. ground_truth_number_of_elements) then
                    do i=1, size(domain_element_list)
                        if ((ground_truth_element_list(i)%point_x%x(1) .eq. domain_element_list(i)%point_x%x(1)) .and. &
                            (ground_truth_element_list(i)%point_x%x(2) .eq. domain_element_list(i)%point_x%x(2))) then
                            rv = .true.
                            print *, "   ", i, ": OK"
                        else
                            rv = .false.
                            print *, "   ", i, ": FAILED"
                        end if
                    end do
                else
                    rv = .false.
                    print *, "   number of domain elements is NOT equal with the ground truth"
                    print *, "   number of domain elements: ", size(domain_element_list)
                    print *, "   number of ground truth elements: ", size(ground_truth_element_list)
                end if
            else
                rv = .false.
                print *, "   number of domain elements is NOT equal with the ground truth"
                print *, "   number of domain elements: ", 0
                print *, "   number of ground truth elements: ", size(ground_truth_element_list)
            end if

            ! clean up
            if (allocated(ground_truth_element_list)) then
                deallocate (ground_truth_element_list)
            end if
            if (allocated(domain_element_list)) then
                deallocate (domain_element_list)
            end if
            if (allocated(element_number)) then
                deallocate (element_number)
            end if
        end function test_lib_tree_get_domain_e2

        function test_lib_tree_get_domain_e3() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**1

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            type(lib_tree_universal_index) :: uindex
            integer(kind=UINDEX_BYTES) :: k
            type(lib_tree_data_element), dimension(:), allocatable :: domain_element_list
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number

            type(lib_tree_data_element), dimension(:), allocatable :: ground_truth_element_list
            integer :: ground_truth_number_of_elements
            integer :: number_of_elements
            ! ---- create correspondece vector ----

            margin = 300


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, ground_truth_number_of_elements
                ground_truth_element_list(i) = element_list(i+ground_truth_number_of_elements)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, ground_truth_number_of_elements
                ground_truth_element_list(i) = element_list(i+ground_truth_number_of_elements)
            end do
#endif
            element_list(:)%element_type = 1

            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .ne. list_length) then
                rv = .false.
                print *, "test_lib_tree_get_domain_e3_create_correspondence_vector: FAILD"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"

                return
            end if

            ! ---- test get_domain_e3 ----
            uindex%n = 0
            uindex%l = 2
            k = 1

            domain_element_list = lib_tree_get_domain_e3(k,uindex, element_number)

            print *, "test_lib_tree_get_domain_e3:"
            rv = .true.
            if (allocated(domain_element_list)) then
                number_of_elements = size(domain_element_list)

                if (number_of_elements .eq. ground_truth_number_of_elements) then
                    do i=1, size(domain_element_list)
                        if ((ground_truth_element_list(i)%point_x%x(1) .eq. domain_element_list(i)%point_x%x(1)) .and. &
                            (ground_truth_element_list(i)%point_x%x(2) .eq. domain_element_list(i)%point_x%x(2))) then
                            rv = .true.
                            print *, "   ", i, ": OK"
                        else
                            rv = .false.
                            print *, "   ", i, ": FAILED"
                        end if
                    end do
                else
                    rv = .false.
                    print *, "   number of domain elements is NOT equal with the ground truth"
                    print *, "   number of domain elements: ", size(domain_element_list)
                    print *, "   number of ground truth elements: ", size(ground_truth_element_list)
                end if
            else
                rv = .false.
                print *, "   number of domain elements is NOT equal with the ground truth"
                print *, "   number of domain elements: ", 0
                print *, "   number of ground truth elements: ", size(ground_truth_element_list)
            end if

            ! clean up
            if (allocated(ground_truth_element_list)) then
                deallocate (ground_truth_element_list)
            end if
            if (allocated(domain_element_list)) then
                deallocate (domain_element_list)
            end if
            if (allocated(element_number)) then
                deallocate (element_number)
            end if
        end function test_lib_tree_get_domain_e3

        function test_lib_tree_get_domain_e4() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**1

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            type(lib_tree_universal_index) :: uindex
            integer(kind=UINDEX_BYTES) :: k
            type(lib_tree_data_element), dimension(:), allocatable :: domain_element_list
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number

            type(lib_tree_data_element), dimension(:), allocatable :: ground_truth_element_list
            integer :: ground_truth_number_of_elements
            integer :: number_of_elements
            ! ---- create correspondece vector ----

            margin = 300


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, ground_truth_number_of_elements
                ground_truth_element_list(i) = element_list(i+ground_truth_number_of_elements)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
            ground_truth_number_of_elements = 5
            allocate(ground_truth_element_list(ground_truth_number_of_elements))
            do i=1, ground_truth_number_of_elements
                ground_truth_element_list(i) = element_list(i+ground_truth_number_of_elements)
            end do
#endif
            element_list(:)%element_type = 1

            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .ne. list_length) then
                rv = .false.
                print *, "test_lib_tree_get_domain_e4_create_correspondence_vector: FAILD"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"

                return
            end if

            ! ---- test get_domain_e4 ----
            uindex%n = 0
            uindex%l = 2
            k = 1

            domain_element_list = lib_tree_get_domain_e4(k,uindex, element_number)

            print *, "test_lib_tree_get_domain_e4:"
            rv = .true.
            if (allocated(domain_element_list)) then
                number_of_elements = size(domain_element_list)

                if (number_of_elements .eq. ground_truth_number_of_elements) then
                    do i=1, size(domain_element_list)
                        if ((ground_truth_element_list(i)%point_x%x(1) .eq. domain_element_list(i)%point_x%x(1)) .and. &
                            (ground_truth_element_list(i)%point_x%x(2) .eq. domain_element_list(i)%point_x%x(2))) then
                            rv = .true.
                            print *, "   ", i, ": OK"
                        else
                            rv = .false.
                            print *, "   ", i, ": FAILED"
                        end if
                    end do
                else
                    rv = .false.
                    print *, "   number of domain elements is NOT equal with the ground truth"
                    print *, "   number of domain elements: ", size(domain_element_list)
                    print *, "   number of ground truth elements: ", size(ground_truth_element_list)
                end if
            else
                rv = .false.
                print *, "   number of domain elements is NOT equal with the ground truth"
                print *, "   number of domain elements: ", 0
                print *, "   number of ground truth elements: ", size(ground_truth_element_list)
            end if

            ! clean up
            if (allocated(ground_truth_element_list)) then
                deallocate (ground_truth_element_list)
            end if
            if (allocated(domain_element_list)) then
                deallocate (domain_element_list)
            end if
            if (allocated(element_number)) then
                deallocate (element_number)
            end if
        end function test_lib_tree_get_domain_e4

        function test_lib_tree_get_domain_i2() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_tree_universal_index) :: uindex
            integer(kind=UINDEX_BYTES) :: k
            type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i2
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_i2

            integer(kind=1) :: i
            k = 1

#if (_FMM_DIMENSION_ == 2)
            ! ---------------------------------
            ! |  5    |  7    | 13    | 15    |
            ! |       |       |       |       |
            ! |-------1---------------3-------|
            ! |  4    |  6    | 12    | 14    |
            ! |    *  |    *  |   *   |       |
            ! |---------------0---------------|
            ! |  1    |  3    | 9     | 11    |
            ! |    *  |    X* |   *   |       |
            ! |-------0---------------2-------|
            ! |  0    |  2    | 8     | 10    |
            ! |    *  |    *  |   *   |       |  X: uindex
            ! ---------------------------------  *: ground_truth_uindex_i2
            uindex%n = 3
            uindex%l = 2
            allocate(uindex_i2, source=lib_tree_get_domain_i2(uindex, k))

            allocate(ground_truth_uindex_i2(3**TREE_DIMENSIONS))
            ground_truth_uindex_i2(:)%l = uindex%l
            ground_truth_uindex_i2(1)%n = 0
            ground_truth_uindex_i2(2)%n = 1
            ground_truth_uindex_i2(3)%n = 2
            ground_truth_uindex_i2(4)%n = 3
            ground_truth_uindex_i2(5)%n = 4
            ground_truth_uindex_i2(6)%n = 6
            ground_truth_uindex_i2(7)%n = 8
            ground_truth_uindex_i2(8)%n = 9
            ground_truth_uindex_i2(9)%n = 12

            if (sum(uindex_i2(:)%n) .eq. sum(ground_truth_uindex_i2(:)%n)) then
                rv = .true.
                print *, "test_lib_tree_get_domain_i2: OK"
            else
                rv = .false.
                print *, "test_lib_tree_get_domain_i2: FAILED"
            end if
#elif (_FMM_DIMENSION_ == 3)
            !  /                               /
            ! ---------------------------------
            ! |       |       |       |       |
            ! |       |       |       |       |
            ! |-------1---------------3-------|
            ! |       |       |       |       |
            ! |       |       |       |       |
            ! |---------------0---------------|
            ! |  1    |  3    |       |       |
            ! |    *  |    *  |       |       |
            ! |-------0---------------2-------|
            ! |  0    |  2    |       |       |
            ! |    *X |    *  |       |       |  X: uindex
            ! ---------------------------------/ *: ground_truth_uindex_i2
            uindex%n = 0
            uindex%l = 2
            allocate(uindex_i2, source=lib_tree_get_domain_i2(uindex, k))

            allocate(ground_truth_uindex_i2(8))
            ground_truth_uindex_i2(:)%l = uindex%l
            ground_truth_uindex_i2(1)%n = 0
            ground_truth_uindex_i2(2)%n = 1
            ground_truth_uindex_i2(3)%n = 2
            ground_truth_uindex_i2(4)%n = 3
            ground_truth_uindex_i2(5)%n = 4
            ground_truth_uindex_i2(6)%n = 5
            ground_truth_uindex_i2(7)%n = 6
            ground_truth_uindex_i2(8)%n = 7

            if (sum(uindex_i2(:)%n) .eq. sum(ground_truth_uindex_i2(:)%n)) then
                rv = .true.
                print *, "test_lib_tree_get_domain_i2: OK"
            else
                rv = .false.
                print *, "test_lib_tree_get_domain_i2: FAILED"
            end if
#endif

        end function test_lib_tree_get_domain_i2

        function test_lib_tree_get_domain_i4() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_tree_universal_index) :: uindex
            integer(kind=UINDEX_BYTES) :: k
            type(lib_tree_universal_index), dimension(:), allocatable :: uindex_i4
            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_i4

            ! ---------------------------------
            ! |  5    |  7    | 13    | 15    |
            ! |    *  |    *  |    *  |    *  |
            ! |-------1---------------3-------|
            ! |  4    |  6    | 12    | 14    |
            ! |       |       |       |    *  |
            ! |---------------0---------------|
            ! |  1    |  3    | 9     | 11    |
            ! |       |    X  |       |    *  |
            ! |-------0---------------2-------|
            ! |  0    |  2    | 8     | 10    |
            ! |       |       |       |    *  |  X: uindex
            ! ---------------------------------  *: ground_truth_uindex_i4
            uindex%n = 3
            uindex%l = 3
            k = 1
            allocate(uindex_i4, source=lib_tree_get_domain_i4(uindex, k))

            allocate(ground_truth_uindex_i4(7))
            ground_truth_uindex_i4(:)%l = uindex%l
            ground_truth_uindex_i4(1)%n = 5
            ground_truth_uindex_i4(2)%n = 7
            ground_truth_uindex_i4(3)%n = 10
            ground_truth_uindex_i4(4)%n = 11
            ground_truth_uindex_i4(5)%n = 13
            ground_truth_uindex_i4(6)%n = 14
            ground_truth_uindex_i4(7)%n = 15

            if (sum(uindex_i4(:)%n) .eq. sum(ground_truth_uindex_i4(:)%n)) then
                rv = .true.
                print *, "test_lib_tree_get_domain_i4: OK"
            else
                rv = .false.
                print *, "test_lib_tree_get_domain_i4: FAILED"
            end if

        end function test_lib_tree_get_domain_i4

        function test_lib_tree_get_level_min() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=COORDINATE_BINARY_BYTES) :: k
            integer(kind=UINDEX_BYTES) :: l_min
            integer(kind=UINDEX_BYTES) :: ground_truth_l_min

            k = 1
            ground_truth_l_min = 2

            l_min = lib_tree_get_level_min(k)

            if (l_min .eq. ground_truth_l_min) then
                print *, "test_lib_tree_get_level_min: ", "OK"
                rv = .true.
            else
                print *, "test_lib_tree_get_level_min: ", "FAILED"
                rv = .false.
            end if
        end function test_lib_tree_get_level_min

        function test_lib_tree_get_level_max() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 5

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i

            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(list_length) &
                :: gt_correspondence_vector_sorted_data_elements ! gt: ground truth

            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(3) :: s
            integer(kind=1), dimension(3) :: ground_truth_l_max
            integer(kind=1) :: l_max
            ! set up the environment
            margin = 200

            call lib_tree_destructor()

#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 1.0 - (0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = 1.0 - (0.999 * i)/(1.0*list_length)
            end do
#endif
            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            gt_correspondence_vector_sorted_data_elements = (/5,4,3,2,1/)

            call lib_tree_create_correspondece_vector_sorted_data_elements()

            rv = .true.
            do i=1, list_length
                if (lib_tree_correspondence_vector_sorted_data_elements(i) .ne. &
                    gt_correspondence_vector_sorted_data_elements(i)) then
                    rv = .false.
                end if
            end do
            if (rv .eqv. .false.) then
                print *, "test_lib_tree_get_level_max: FAILED"
                print *, "  test_lib_tree_create_correspondece_vector: FAILED"
                return
            end if

            ! begin of the test of the test_lib_tree_get_level_max function

            !
            ! |-------------------0 <- level
            ! |         1         |
            ! |         |      *  |
            ! |         |         |
            ! |         |  *      |
            ! |         |         |
            ! |---------|---------|
            ! |    2   *|         |
            ! |    |    |         |
            ! |----|----|         |
            ! |   *|    |         |
            ! |*   |    |         |
            ! --------------------|
            print *, "test_lib_tree_get_level_max:"
            rv = .true.
            s = (/1, 2, 3/)
            ground_truth_l_max = (/int(3,1), int(2,1), int(1,1)/)
            do i=1, 3
                l_max = lib_tree_get_level_max(s(i))
                if (l_max .eq. ground_truth_l_max(i)) then
                    print *, "  s = ",s(i) , " : OK"
                else
                    print *, "  s = ",s(i) , " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_get_level_max

        function test_lib_tree_get_number_of_boxes() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=1) :: l
            integer(kind=UINDEX_BYTES) :: number_of_boxes
            integer(kind=UINDEX_BYTES) :: ground_truth_number_of_boxes

            l = 2
#if (_FMM_DIMENSION_ == 2)
            ground_truth_number_of_boxes = 16! 2**(2*2)
#elif (_FMM_DIMENSION_ == 3)
            ground_truth_number_of_boxes = 64! 2**(3*2)
#endif
            number_of_boxes = lib_tree_get_number_of_boxes(l)

            if (ground_truth_number_of_boxes .eq. number_of_boxes) then
                print *, "test_lib_tree_get_number_of_boxes: OK"
                rv = .true.
            else
                print *, "test_lib_tree_get_number_of_boxes: FAILED"
                rv = .false.
            end if

        end function test_lib_tree_get_number_of_boxes

        function test_lib_tree_create_correspondece_vector_sorted_data_elements() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 5

            integer(kind=1), parameter :: l_th = 8 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i

            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(list_length) &
                :: gt_correspondence_vector_sorted_data_elements ! gt: ground truth

            ! set up the environment
            margin = 200

#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 1.0 - (0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 1.0 - (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = 1.0 - (0.999 * i)/(1.0*list_length)
            end do
#endif
            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            ! begin of the test of the lib_tree_create_correspondece_vector_sorted_data_elements function
            gt_correspondence_vector_sorted_data_elements = (/5,4,3,2,1/)

            call lib_tree_create_correspondece_vector_sorted_data_elements()


            print *, "test_lib_tree_create_correspondece_vector_sorted_data_elements:"
            rv = .true.
            do i=1, list_length
                if (lib_tree_correspondence_vector_sorted_data_elements(i) .eq. &
                    gt_correspondence_vector_sorted_data_elements(i)) then
                    print *, "  ", i, ": OK"
                else
                    rv = .false.
                    print *, "  ", i, ": FAILED"
                end if
            end do
        end function test_lib_tree_create_correspondece_vector_sorted_data_elements

        function test_lib_tree_create_correspondence_vector() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 5

            integer(kind=1), parameter :: l_th = 2 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            margin = 200


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.999 * i)/(1.0*list_length)
            end do
#endif
            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + size(lib_tree_correspondence_vector(i)%data_element_number)
                end if
            end do


            if (number .eq. list_length) then
                rv = .true.
                print *, "test_lib_tree_create_correspondece_vector: OK"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
            else
                rv = .false.
                print *, "test_lib_tree_create_correspondece_vector: FAILED"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"
            end if

        end function test_lib_tree_create_correspondence_vector

        function test_lib_tree_get_element_from_correspondence_vector() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_tree_universal_index) :: uindex
            type(lib_tree_data_element), dimension(:), allocatable :: data_element
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: buffer_element_number

            ! generate dataset
            integer(kind=4), parameter :: list_length = 5

            integer(kind=1), parameter :: l_th = 8 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: wrong

            margin = 200


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.999 * i)/(1.0*list_length)
            end do
#endif
            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)


            ! test dataset
            wrong = 0
            do i=1, list_length
                uindex = lib_tree_hf_get_universal_index(element_list(i)%point_x, lib_tree_l_th)
                data_element = lib_tree_get_element_from_correspondence_vector(uindex, buffer_element_number)
                if ((element_list(i)%point_x%x(1) .ne. data_element(1)%point_x%x(1)) .or. &
                    (element_list(i)%point_x%x(2) .ne. data_element(1)%point_x%x(2))) then
                   wrong = wrong + 1
                end if
            end do

            if (wrong .eq. 0) then
                rv = .true.
                print *, "test_lib_tree_get_element_from_correspondence_vector: ", "OK"
            else
                print *, "test_lib_tree_get_element_from_correspondence_vector: ", "FAILED"
                print *, "  wrong elements: ", wrong
                print *, "  wrong elements [%]: ", 100.0 * wrong / list_length
                rv = .false.
            end if

        end function test_lib_tree_get_element_from_correspondence_vector

    end function lib_tree_test_functions

    subroutine lib_tree_benchmark
        implicit none

        call benchmark_lib_tree_create_correspondece_vector()
        call benchmark_lib_tree_get_element_from_correspondence_vector()

    contains

        subroutine benchmark_lib_tree_create_correspondece_vector()
            implicit none

            real :: start, finish
            double precision :: delta

            integer(kind=4), parameter :: list_length = 10**5
            integer(kind=1), parameter :: number_of_runs = 100

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i

            margin = 400


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.999 + (-0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
                element_list(i)%point_x%x(3) = 0.9 + (-0.9 * i)/list_length
            end do
#endif

            lib_tree_data_element_list = element_list

            call cpu_time(start)
            do i=1, number_of_runs
                call lib_tree_create_correspondence_vector(l_th, margin)
            end do
            call cpu_time(finish)
            print *, "benchmark_lib_tree_create_correspondece_vector:"
            delta = (finish-start)/number_of_runs
            print *, " create correspondence_vector time: ", delta, " seconds."
            print *, "   total run time: ", finish-start, " seconds"

        end subroutine benchmark_lib_tree_create_correspondece_vector

        subroutine benchmark_lib_tree_get_element_from_correspondence_vector()
            implicit none
            real :: start, finish
            double precision :: delta

            integer(kind=4), parameter :: list_length = 10**5
            integer(kind=4), parameter :: number_of_runs = 10**7

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: buffer_element_number
            type(lib_tree_data_element), dimension(:), allocatable :: data_element
            type(lib_tree_universal_index) :: uindex
            integer(kind=2) :: margin

            integer(kind=4) :: i

            margin = 400


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.999 + (-0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
                element_list(i)%point_x%x(3) = 0.9 + (-0.9 * i)/list_length
            end do
#endif
            lib_tree_data_element_list = element_list
            call lib_tree_create_correspondence_vector(l_th, margin)

            uindex = lib_tree_hf_get_universal_index(element_list(1)%point_x, lib_tree_l_th)

            call cpu_time(start)
            do i=1, number_of_runs
                data_element = lib_tree_get_element_from_correspondence_vector(uindex, buffer_element_number)
            end do
            call cpu_time(finish)
            print *, "benchmark_lib_tree_get_element_from_correspondence_vector:"
            delta = (finish-start)/number_of_runs
            print *, " get element from correspondence vector time: ", delta, " seconds."
            print *, "   total run time: ", finish-start, " seconds"

        end subroutine benchmark_lib_tree_get_element_from_correspondence_vector
    end subroutine lib_tree_benchmark
end module lib_tree
