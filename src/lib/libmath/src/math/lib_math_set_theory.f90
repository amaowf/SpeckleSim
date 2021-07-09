module lib_math_set_theory
    use lib_hash_function
    implicit none

    private

    public :: get_sets
    public :: get_union
    public :: get_intersection
    public :: get_relative_complement_of_b_in_a
    public :: make_unique

    public :: lib_math_set_theory_test_functions

    ! --- parameter ---
    integer, parameter :: MAXIMUM_NUMBER_OF_HASH_RUNS = 200
    double precision, parameter :: MARGIN = 2.2d0

    integer, parameter :: SET_INTERSECTION = 1
    integer, parameter :: SET_RELATIVE_COMPLEMENT_OF_B_IN_A = 2
    integer, parameter :: SET_RELATIVE_COMPLEMENT_OF_A_IN_B = 3

    ! --- type definitions ---
    type hash_list_type_4_byte
        integer(kind = 4) :: value
        integer :: hash_runs
        integer :: algebra_of_sets
    end type

    type hash_list_type_8_byte
        integer(kind = 8) :: value
        integer :: hash_runs
        integer :: algebra_of_sets
    end type

    ! --- interface ---
    interface get_union
        module procedure lib_math_set_theory_get_union_4_byte
        module procedure lib_math_set_theory_get_union_8_byte
    end interface

    interface get_sets
        module procedure lib_math_set_theory_get_sets_4_byte
        module procedure lib_math_set_theory_get_sets_8_byte
    end interface

    interface get_intersection
        module procedure lib_math_set_theory_get_intersection_4_byte
        module procedure lib_math_set_theory_get_intersection_8_byte
    end interface

    interface get_relative_complement_of_b_in_a
        module procedure lib_math_set_theory_get_relative_complement_of_b_in_a_4_byte
        module procedure lib_math_set_theory_get_relative_complement_of_b_in_a_8_byte
    end interface

    interface make_unique
        module procedure lib_math_set_theory_make_unique_4_byte
        module procedure lib_math_set_theory_make_unique_8_byte
    end interface

    contains

    ! Argument
    ! ----
    !   a: integer(kind = 4), dimension(:)
    !       set a
    !   b: integer(kind = 4), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 4), dimension(:)
    !       the union of a and b
    function lib_math_set_theory_get_union_4_byte(a, b) result(res)
        ! dummy
        integer(kind = 4), dimension(:), allocatable, intent(in) :: a
        integer(kind = 4), dimension(:), allocatable, intent(in) :: b
        integer(kind = 4), dimension(:), allocatable :: res

        ! auxiliary
        type(hash_list_type_4_byte), dimension(:), allocatable :: hash_list
        integer(kind = 4) :: max_value
        integer(kind = 4) :: hash
        integer :: counter
        integer :: i
        integer :: ii

        i = size(a, 1) + size(b, 1)

        if (i .gt. 0) then

            allocate(hash_list(int(i * MARGIN)))

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0

            max_value = size(hash_list)

            counter = 0

            ! has set a
            do i = lbound(a, 1), ubound(a, 1)
                ! hash runs
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    hash = hash_fnv1a(a(i), max_value)
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = a(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. a(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. a(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: a"
                        exit
                    end if
                end do
            end do

            ! hash set b
            do i = lbound(b, 1), ubound(b, 1)
                ! hash runs
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    hash = hash_fnv1a(b(i), max_value)
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = b(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. b(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. b(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: b"
                        exit
                    end if
                end do
            end do

            if (counter .gt. 0) then
                allocate(res(counter))

                counter = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        res(counter) = hash_list(i)%value
                    end if
                end do
            end if

        end if

    end function lib_math_set_theory_get_union_4_byte

    ! Argument
    ! ----
    !   a: integer(kind = 8), dimension(:)
    !       set a
    !   b: integer(kind = 8), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 8), dimension(:)
    !       the union of a and b
    function lib_math_set_theory_get_union_8_byte(a, b) result(res)
        ! dummy
        integer(kind = 8), dimension(:), allocatable, intent(in) :: a
        integer(kind = 8), dimension(:), allocatable, intent(in) :: b
        integer(kind = 8), dimension(:), allocatable :: res

        ! auxiliary
        type(hash_list_type_8_byte), dimension(:), allocatable :: hash_list
        integer(kind = 8) :: max_value
        integer(kind = 8) :: hash
        integer :: counter
        integer :: i
        integer :: ii

        i = size(a, 1) + size(b, 1)

        if (i .gt. 0) then

            allocate(hash_list(int(i * MARGIN)))

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0

            max_value = size(hash_list)

            counter = 0

            ! has set a
            do i = lbound(a, 1), ubound(a, 1)
                ! hash runs
                hash = hash_fnv1a(a(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = a(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. a(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. a(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: a"
                        exit
                    end if
                end do
            end do

            ! hash set b
            do i = lbound(b, 1), ubound(b, 1)
                ! hash runs
                hash = hash_fnv1a(b(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = b(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. b(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. b(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: b"
                        exit
                    end if
                end do
            end do

            if (counter .gt. 0) then
                allocate(res(counter))

                counter = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        res(counter) = hash_list(i)%value
                    end if
                end do
            end if

        end if

    end function lib_math_set_theory_get_union_8_byte

    ! Argument
    ! ----
    !   a: integer(kind = 4), dimension(:)
    !       set a
    !   b: integer(kind = 4), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 4), dimension(:)
    !       the union of a and b
    function lib_math_set_theory_get_intersection_4_byte(a, b) result(res)
        ! dummy
        integer(kind = 4), dimension(:), allocatable, intent(in) :: a
        integer(kind = 4), dimension(:), allocatable, intent(in) :: b
        integer(kind = 4), dimension(:), allocatable :: res

        ! auxiliary
        integer(kind = 4), dimension(:), allocatable :: buffer_b_in_a
        integer(kind = 4), dimension(:), allocatable :: buffer_a_in_b

        call get_sets(a, b, res, buffer_b_in_a, buffer_a_in_b)

    end function lib_math_set_theory_get_intersection_4_byte

    ! Argument
    ! ----
    !   a: integer(kind = 8), dimension(:)
    !       set a
    !   b: integer(kind = 8), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 8), dimension(:)
    !       the union of a and b
    function lib_math_set_theory_get_intersection_8_byte(a, b) result(res)
        ! dummy
        integer(kind = 8), dimension(:), allocatable, intent(in) :: a
        integer(kind = 8), dimension(:), allocatable, intent(in) :: b
        integer(kind = 8), dimension(:), allocatable :: res

        ! auxiliary
        integer(kind = 8), dimension(:), allocatable :: buffer_b_in_a
        integer(kind = 8), dimension(:), allocatable :: buffer_a_in_b

        call get_sets(a, b, res, buffer_b_in_a, buffer_a_in_b)

    end function lib_math_set_theory_get_intersection_8_byte

    ! Argument
    ! ----
    !   a: integer(kind = 4), dimension(:)
    !       set a
    !   b: integer(kind = 4), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 4), dimension(:)
    !       the realative complement of b in a
    function lib_math_set_theory_get_relative_complement_of_b_in_a_4_byte(a, b) result(res)
        ! dummy
        integer(kind = 4), dimension(:), allocatable, intent(in) :: a
        integer(kind = 4), dimension(:), allocatable, intent(in) :: b
        integer(kind = 4), dimension(:), allocatable :: res

        ! auxiliary
        integer(kind = 4), dimension(:), allocatable :: buffer_intersection
        integer(kind = 4), dimension(:), allocatable :: buffer_a_in_b

        call get_sets(a, b, buffer_intersection, res, buffer_a_in_b)

    end function lib_math_set_theory_get_relative_complement_of_b_in_a_4_byte

    ! Argument
    ! ----
    !   a: integer(kind = 8), dimension(:)
    !       set a
    !   b: integer(kind = 8), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   res: integer(kind = 8), dimension(:)
    !       the realative complement of b in a
    function lib_math_set_theory_get_relative_complement_of_b_in_a_8_byte(a, b) result(res)
        ! dummy
        integer(kind = 8), dimension(:), allocatable, intent(in) :: a
        integer(kind = 8), dimension(:), allocatable, intent(in) :: b
        integer(kind = 8), dimension(:), allocatable :: res

        ! auxiliary
        integer(kind = 8), dimension(:), allocatable :: buffer_intersection
        integer(kind = 8), dimension(:), allocatable :: buffer_a_in_b

        call get_sets(a, b, buffer_intersection, res, buffer_a_in_b)

    end function lib_math_set_theory_get_relative_complement_of_b_in_a_8_byte

    ! Argument
    ! ----
    !   a: integer(kind = 4), dimension(:)
    !       set a
    !   b: integer(kind = 4), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   intersection: integer, dimension(:)
    !       the intersection of a and b
    !   relative_complement_of_b_in_a: integer, dimension(:)
    !       the realative complement of b in a
    !   relative_complement_of_a_in_b: integer, dimension(:)
    !       the realative complement of a in b
    subroutine lib_math_set_theory_get_sets_4_byte(a, b, &
                                                   intersection, &
                                                   relative_complement_of_b_in_a, &
                                                   relative_complement_of_a_in_b)
        ! dummy
        integer(kind = 4), dimension(:), allocatable, intent(in) :: a
        integer(kind = 4), dimension(:), allocatable, intent(in) :: b
        integer(kind = 4), dimension(:), allocatable, intent(out) :: intersection
        integer(kind = 4), dimension(:), allocatable, intent(out) :: relative_complement_of_b_in_a
        integer(kind = 4), dimension(:), allocatable, intent(out) :: relative_complement_of_a_in_b

        ! auxiliary
        integer(kind = 4), dimension(:), allocatable :: a_unique
        integer(kind = 4), dimension(:), allocatable :: b_unique
        type(hash_list_type_4_byte), dimension(:), allocatable :: hash_list
        integer(kind = 4) :: max_value
        integer(kind = 4) :: hash
        integer, dimension(3) :: counter
        integer :: i
        integer :: ii

        i = size(a, 1) + size(b, 1)

        if (i .gt. 0) then

            allocate(a_unique(size(a)))
            allocate(b_unique(size(b)))

            a_unique = a
            b_unique = b

            call make_unique(a_unique)
            call make_unique(b_unique)

            i = size(a_unique, 1) + size(b_unique, 1)

            if (int(i * MARGIN) .lt. 64) then
                allocate(hash_list(64))
            else
                allocate(hash_list(int(i * MARGIN)))
            end if

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0
            hash_list(:)%algebra_of_sets = -1

            max_value = int(size(hash_list), kind = 4)

            ! hash set a
            do i = lbound(a_unique, 1), ubound(a_unique, 1)
                ! hash runs
                hash = hash_fnv1a(a_unique(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = a_unique(i)
                        hash_list(hash)%hash_runs = ii
                        hash_list(hash)%algebra_of_sets = SET_RELATIVE_COMPLEMENT_OF_B_IN_A
                        exit
                    else if (hash_list(hash)%value .ne. a_unique(i)) then
                        ! hash collision -> re-hash
!                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. a_unique(i)) then
                        ! uindex duplicate found
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash a_unique: duplicate found"
                        exit
                    else
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash: a_unique"
                        exit
                    end if
                end do
            end do

            counter(:) = 0
            counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) = size(a_unique)
            ! hash set b
            do i = lbound(b_unique, 1), ubound(b_unique, 1)
                ! hash runs
                hash = hash_fnv1a(b_unique(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = b_unique(i)
                        hash_list(hash)%hash_runs = ii
                        hash_list(hash)%algebra_of_sets = SET_RELATIVE_COMPLEMENT_OF_A_IN_B
                        counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) = counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) + 1
                        exit
                    else if (hash_list(hash)%value .ne. b_unique(i)) then
                        ! hash collision -> re-hash
!                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. b_unique(i)) then
                        ! uindex duplicate found
                        hash_list(hash)%algebra_of_sets = SET_INTERSECTION
                        counter(SET_INTERSECTION) = counter(SET_INTERSECTION) + 1
                        counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) = counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) - 1
                        exit
                    else
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash: b"
                        exit
                    end if
                end do
            end do

            if (counter(SET_INTERSECTION) .gt. 0) then
                allocate(intersection(counter(SET_INTERSECTION)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_INTERSECTION) then
                        ii = ii + 1
                        intersection(ii) = hash_list(i)%value
                    end if
                end do
            end if

            if (counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) .gt. 0) then
                allocate(relative_complement_of_b_in_a(counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_RELATIVE_COMPLEMENT_OF_B_IN_A) then
                        ii = ii + 1
                        relative_complement_of_b_in_a(ii) = hash_list(i)%value
                    end if
                end do
            end if

            if (counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) .gt. 0) then
                allocate(relative_complement_of_a_in_b(counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_RELATIVE_COMPLEMENT_OF_A_IN_B) then
                        ii = ii + 1
                        relative_complement_of_a_in_b(ii) = hash_list(i)%value
                    end if
                end do
            end if
        end if

    end subroutine lib_math_set_theory_get_sets_4_byte

    ! Argument
    ! ----
    !   a: integer(kind = 8), dimension(:)
    !       set a
    !   b: integer(kind = 8), dimension(:)
    !       set b
    !
    ! Results
    ! ----
    !   intersection: integer, dimension(:)
    !       the intersection of a and b
    !   relative_complement_of_b_in_a: integer, dimension(:)
    !       the realative complement of b in a
    !   relative_complement_of_a_in_b: integer, dimension(:)
    !       the realative complement of a in b
    subroutine lib_math_set_theory_get_sets_8_byte(a, b, &
                                                   intersection, &
                                                   relative_complement_of_b_in_a, &
                                                   relative_complement_of_a_in_b)
        ! dummy
        integer(kind = 8), dimension(:), allocatable, intent(in) :: a
        integer(kind = 8), dimension(:), allocatable, intent(in) :: b
        integer(kind = 8), dimension(:), allocatable, intent(out) :: intersection
        integer(kind = 8), dimension(:), allocatable, intent(out) :: relative_complement_of_b_in_a
        integer(kind = 8), dimension(:), allocatable, intent(out) :: relative_complement_of_a_in_b

        ! auxiliary
        integer(kind = 8), dimension(:), allocatable :: a_unique
        integer(kind = 8), dimension(:), allocatable :: b_unique
        type(hash_list_type_8_byte), dimension(:), allocatable :: hash_list
        integer(kind = 8) :: max_value
        integer(kind = 8) :: hash
        integer, dimension(3) :: counter
        integer :: i
        integer :: ii

        i = size(a, 1) + size(b, 1)

        if (i .gt. 0) then

            allocate(a_unique(size(a)))
            allocate(b_unique(size(b)))

            a_unique = a
            b_unique = b

            call make_unique(a_unique)
            call make_unique(b_unique)

            i = size(a_unique, 1) + size(b_unique, 1)

            if (int(i * MARGIN) .lt. 64) then
                allocate(hash_list(64))
            else
                allocate(hash_list(int(i * MARGIN)))
            end if

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0
            hash_list(:)%algebra_of_sets = -1

            max_value = int(size(hash_list), kind = 8)

            ! hash set a
            do i = lbound(a_unique, 1), ubound(a_unique, 1)
                ! hash runs
                hash = hash_fnv1a(a_unique(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = a_unique(i)
                        hash_list(hash)%hash_runs = ii
                        hash_list(hash)%algebra_of_sets = SET_RELATIVE_COMPLEMENT_OF_B_IN_A
                        exit
                    else if (hash_list(hash)%value .ne. a_unique(i)) then
                        ! hash collision -> re-hash
!                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. a_unique(i)) then
                        ! uindex duplicate found
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash a_unique: duplicate found"
                        exit
                    else
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash: a_unique"
                        exit
                    end if
                end do
            end do

            counter(:) = 0
            counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) = size(a_unique)
            ! hash set b
            do i = lbound(b_unique, 1), ubound(b_unique, 1)
                ! hash runs
                hash = hash_fnv1a(b_unique(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = b_unique(i)
                        hash_list(hash)%hash_runs = ii
                        hash_list(hash)%algebra_of_sets = SET_RELATIVE_COMPLEMENT_OF_A_IN_B
                        counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) = counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) + 1
                        exit
                    else if (hash_list(hash)%value .ne. b_unique(i)) then
                        ! hash collision -> re-hash
!                        hash = IEOR(hash, int(ii, 4))
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. b_unique(i)) then
                        ! uindex duplicate found
                        hash_list(hash)%algebra_of_sets = SET_INTERSECTION
                        counter(SET_INTERSECTION) = counter(SET_INTERSECTION) + 1
                        counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) = counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) - 1
                        exit
                    else
                        print *, "lib_math_set_theory_get_sets: ERROR"
                        print *, "  hash: b"
                        exit
                    end if
                end do
            end do

            if (counter(SET_INTERSECTION) .gt. 0) then
                allocate(intersection(counter(SET_INTERSECTION)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_INTERSECTION) then
                        ii = ii + 1
                        intersection(ii) = hash_list(i)%value
                    end if
                end do
            end if

            if (counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A) .gt. 0) then
                allocate(relative_complement_of_b_in_a(counter(SET_RELATIVE_COMPLEMENT_OF_B_IN_A)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_RELATIVE_COMPLEMENT_OF_B_IN_A) then
                        ii = ii + 1
                        relative_complement_of_b_in_a(ii) = hash_list(i)%value
                    end if
                end do
            end if

            if (counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B) .gt. 0) then
                allocate(relative_complement_of_a_in_b(counter(SET_RELATIVE_COMPLEMENT_OF_A_IN_B)))

                ii = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0 &
                        .and. hash_list(i)%algebra_of_sets .eq. SET_RELATIVE_COMPLEMENT_OF_A_IN_B) then
                        ii = ii + 1
                        relative_complement_of_a_in_b(ii) = hash_list(i)%value
                    end if
                end do
            end if
        end if

    end subroutine lib_math_set_theory_get_sets_8_byte

    ! Argument
    ! ----
    !   set: integer(kind = 4), dimension(:)
    !       Set to be made unique
    subroutine lib_math_set_theory_make_unique_4_byte(set)
        implicit none
        ! dummy
        integer(kind = 4), dimension(:), allocatable, intent(inout) :: set

        ! auxiliary
        type(hash_list_type_4_byte), dimension(:), allocatable :: hash_list
        integer(kind = 4) :: max_value
        integer(kind = 4) :: hash
        integer :: counter
        integer :: i
        integer :: ii

        i = size(set, 1)

        if (i .gt. 0) then
            if (int(i * MARGIN) .lt. 64) then
                allocate(hash_list(64))
            else
                allocate(hash_list(int(i * MARGIN)))
            end if

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0

            max_value = size(hash_list)

            counter = 0

            ! hash the set
            do i = lbound(set, 1), ubound(set, 1)
                ! hash runs
                hash = hash_fnv1a(set(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = set(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. set(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, ii)
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. set(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: set"
                        exit
                    end if
                end do

                if (ii .gt. MAXIMUM_NUMBER_OF_HASH_RUNS) then
                    print *, "lib_math_set_theory_make_unique: ERROR"
                    print *, "  Maximum number of hash runs exceeded."
                    print *, "  i: ", i
                    print *, "  set: ", set(i)
                end if
            end do

            if (counter .gt. 0) then
                deallocate(set)
                allocate(set(counter))

                counter = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        set(counter) = hash_list(i)%value
                    end if
                end do
            end if

        end if
    end subroutine lib_math_set_theory_make_unique_4_byte

    ! Argument
    ! ----
    !   set: integer(kind = 8), dimension(:)
    !       Set to be made unique
    subroutine lib_math_set_theory_make_unique_8_byte(set)
        implicit none
        ! dummy
        integer(kind = 8), dimension(:), allocatable, intent(inout) :: set

        ! auxiliary
        type(hash_list_type_8_byte), dimension(:), allocatable :: hash_list
        integer(kind = 8) :: max_value
        integer(kind = 8) :: hash
        integer :: counter
        integer :: i
        integer :: ii

        i = size(set, 1)

        if (i .gt. 0) then
            if (int(i * MARGIN) .lt. 64) then
                allocate(hash_list(64))
            else
                allocate(hash_list(int(i * MARGIN)))
            end if

            hash_list(:)%value = -1
            hash_list(:)%hash_runs = 0

            max_value = size(hash_list)

            counter = 0

            ! hash the set
            do i = lbound(set, 1), ubound(set, 1)
                ! hash runs
                hash = hash_fnv1a(set(i), max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hash_list(hash)%hash_runs .eq. 0) then
                        ! new entry at hash_list
                        hash_list(hash)%value = set(i)
                        hash_list(hash)%hash_runs = ii
                        counter = counter + 1
                        exit
                    else if (hash_list(hash)%value .ne. set(i)) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, ii)
                        hash = hash_fnv1a(hash, max_value)
                    else if (hash_list(hash)%value .eq. set(i)) then
                        ! uindex duplicate found
                        exit
                    else
                        print *, "lib_math_set_theory_get_union: ERROR"
                        print *, "  hash: set"
                        exit
                    end if
                end do

                if (ii .gt. MAXIMUM_NUMBER_OF_HASH_RUNS) then
                    print *, "lib_math_set_theory_make_unique: ERROR"
                    print *, "  Maximum number of hash runs exceeded."
                    print *, "  i: ", i
                    print *, "  set: ", set(i)
                end if
            end do

            if (counter .gt. 0) then
                deallocate(set)
                allocate(set(counter))

                counter = 0
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        set(counter) = hash_list(i)%value
                    end if
                end do
            end if

        end if
    end subroutine lib_math_set_theory_make_unique_8_byte

!    function lib_tree_hierarchy_make_uindex_list_unique(uindex_list) result(rv)
!            implicit none
!            ! dummy
!            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list
!            type(lib_tree_universal_index), dimension(:), allocatable :: rv
!
!            ! auxiliary
!            type(hash_list_type), dimension(size(uindex_list) * LIB_TREE_HIERARCHY_MARGIN) :: hash_list
!            type(lib_tree_universal_index) :: buffer_uindex
!            integer(kind=4) :: max_value
!            integer(kind=4) :: hash
!            integer(kind=UINDEX_BYTES) :: counter
!            integer(kind=UINDEX_BYTES) :: i
!            integer(kind=2) :: ii
!
!            if (size(uindex_list) .gt. 0) then
!
!                hash_list(:)%value = -1
!                hash_list(:)%hash_runs = 0
!
!                max_value = size(hash_list)
!
!                counter = 0
!                do i=1, size(uindex_list)
!                    ! hash runs
!                    do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
!                        hash = hash_fnv1a(uindex_list(i)%n, max_value)
!                        if (hash_list(hash)%hash_runs .eq. 0) then
!                            hash_list(hash)%value = uindex_list(i)%n
!                            hash_list(hash)%hash_runs = ii
!                            counter = counter + 1
!                            exit
!                        else if (hash_list(hash)%value .ne. uindex_list(i)%n) then
!                            hash = IEOR(hash, int(ii, 4))
!                            hash = hash_fnv1a(hash, max_value)
!                        end if
!                    end do
!                end do
!
!                allocate (rv(counter))
!
!                counter = 0
!                buffer_uindex%l = uindex_list(1)%l
!                do i=1, size(hash_list)
!                    if (hash_list(i)%hash_runs .ne. 0) then
!                        counter = counter + 1
!                        buffer_uindex%n = hash_list(i)%value
!                        rv(counter) = buffer_uindex
!                    end if
!                end do
!
!            else
!                rv = uindex_list
!            end if
!        end function lib_tree_hierarchy_make_uindex_list_unique
    function lib_math_set_theory_test_functions() result(rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        if (.not. test_lib_math_set_theory_get_intersection()) rv = rv + 1


    contains
        function test_lib_math_set_theory_get_intersection() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            integer, dimension(:), allocatable :: a
            integer, dimension(:), allocatable :: b
            integer, dimension(:), allocatable :: res
            integer, dimension(:), allocatable :: ground_truth

            allocate(a, source = (/1, 2, 4, 4, 6/))
            allocate(b, source = (/4, 5, 6, 6/))
            allocate(ground_truth, source = (/6, 4/))

            res = get_intersection(a, b)

            rv = .true.
            if (size(ground_truth, 1) .eq. size(res, 1)) then
                print *, "test_lib_math_set_theory_get_intersection"
                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. res(i)) then
                        print *, "  ", i, ": OK"
                    else
                        print *, "  ", i, ": FAILED"
                        print *, "    res: ", res(i)
                        print *, "     gt: ", ground_truth(i)
                    end if
                end do
            else
                print *, "test_lib_math_set_theory_get_intersection: ERROR"
                print *, "  res and ground_truth don't have the same size"
                rv = .false.
            end if

        end function test_lib_math_set_theory_get_intersection

    end function lib_math_set_theory_test_functions

end module lib_math_set_theory
