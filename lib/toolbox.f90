
module toolbox
    implicit none
    public :: reallocate

    interface reallocate
        module procedure reallocate_1d_integer_1
        module procedure reallocate_1d_integer_2
        module procedure reallocate_1d_integer_4
        module procedure reallocate_1d_integer_8
        !module procedure reallocate_1d_integer_16
    end interface

    interface concatenate
        module procedure lib_tree_hf_concatenate_1d_integer_1_array
        module procedure lib_tree_hf_concatenate_1d_integer_2_array
        module procedure lib_tree_hf_concatenate_1d_integer_4_array
        module procedure lib_tree_hf_concatenate_1d_integer_8_array
        !module procedure lib_tree_hf_concatenate_1d_integer_16_array
        module procedure lib_tree_hf_concatenate_1d_integer_1_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_2_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_4_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_8_array_single
        !module procedure lib_tree_hf_concatenate_1d_integer_16_array_single
    end interface

contains
    ! reallocate a 1-dimensional array
    !
    ! Arguments
    ! ----
    !   a: 1-dimensional integer array
    !       original array, will be replaced with the resized array
    !   n: positiv integer
    !       number of additional array elements
    !
    ! change log
    ! ----
    !    - ni_new changed to relative length (n)
    !
    ! source: https://gist.github.com/ponderomotion/3527522
    !! Daniel Fletcher 2012
    !! module for increasing array sizes dynamically
    !! currently new indices must be larger than old
    SUBROUTINE reallocate_1d_integer(a,n)
        implicit none

        INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: n
        INTEGER :: ni_old

        ni_old = SIZE(a)

        ALLOCATE(temp(ni_old+n))

        temp(1:ni_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_1d_integer

    subroutine reallocate_1d_integer_1(a,n)
        implicit none

        integer(kind=1),dimension(:),allocatable,intent(inout) :: a
        integer(kind=1),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate( a(n) )
        end if

    end subroutine

    subroutine reallocate_1d_integer_2(a,n)
        implicit none

        integer(kind=2),dimension(:),allocatable,intent(inout) :: a
        integer(kind=2),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate( a(n) )
        end if

    end subroutine

    subroutine reallocate_1d_integer_4(a,n)
        implicit none

        integer(kind=4),dimension(:),allocatable,intent(inout) :: a
        integer(kind=4),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate( a(n) )
        end if
    end subroutine

    subroutine reallocate_1d_integer_8(a,n)
        implicit none

        integer(kind=8),dimension(:),allocatable,intent(inout) :: a
        integer(kind=8),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if (allocated(a)) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate( a(n) )
        end if

    end subroutine

    !subroutine reallocate_1d_integer_16(a,n)
    !    implicit none
    !
    !    integer(kind=16),dimension(:),allocatable,intent(inout) :: a
    !    integer(kind=16),dimension(:),allocatable :: temp
    !    integer,intent(in) :: n
    !    integer :: ni_old
    !
    !    if ( allocated(a) ) then
    !        ni_old = size(a)
    !
    !        allocate(temp(ni_old+n))
    !
    !        temp(1:ni_old) = a
    !
    !        call move_alloc(temp,a)
    !    else
    !        allocate( a(n) )
    !    end if
    !
    !end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_1_array(a, b)
        implicit none
        integer(kind=1), dimension(:), allocatable, intent(inout) :: a
        integer(kind=1), dimension(:), allocatable, intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, size(b))

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_2_array(a, b)
        implicit none
        integer(kind=2), dimension(:), allocatable, intent(inout) :: a
        integer(kind=2), dimension(:), allocatable, intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, size(b))

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_4_array(a, b)
        implicit none
        integer(kind=4), dimension(:), allocatable, intent(inout) :: a
        integer(kind=4), dimension(:), allocatable, intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, size(b))

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_8_array(a, b)
        implicit none
        integer(kind=8), dimension(:), allocatable, intent(inout) :: a
        integer(kind=8), dimension(:), allocatable, intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, size(b))

        a(size_a_org+1:) = b
    end subroutine

    !subroutine lib_tree_hf_concatenate_1d_integer_16_array(a, b)
    !    implicit none
    !    integer(kind=16), dimension(:), allocatable, intent(inout) :: a
    !    integer(kind=16), dimension(:), allocatable, intent(in) :: b
    !
    !    ! auxiliaray
    !    integer :: size_a_org
    !
    !    if (allocated(a)) then
    !        size_a_org = size(a)
    !    else
    !        size_a_org = 0
    !    end if
    !
    !    call reallocate(a, size(b))
    !
    !    a(size_a_org+1:) = b
    !end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_1_array_single(a, b)
        implicit none
        integer(kind=1), dimension(:), allocatable, intent(inout) :: a
        integer(kind=1), intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, 1)

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_2_array_single(a, b)
        implicit none
        integer(kind=2), dimension(:), allocatable, intent(inout) :: a
        integer(kind=2), intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, 1)

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_4_array_single(a, b)
        implicit none
        integer(kind=4), dimension(:), allocatable, intent(inout) :: a
        integer(kind=4), intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, 1)

        a(size_a_org+1:) = b
    end subroutine

    subroutine lib_tree_hf_concatenate_1d_integer_8_array_single(a, b)
        implicit none
        integer(kind=8), dimension(:), allocatable, intent(inout) :: a
        integer(kind=8), intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if (allocated(a)) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if

        call reallocate(a, 1)

        a(size_a_org+1:) = b
    end subroutine

    !subroutine lib_tree_hf_concatenate_1d_integer_16_array_single(a, b)
    !    implicit none
    !    ! dummy
    !    integer(kind=16), dimension(:), allocatable, intent(inout) :: a
    !    integer(kind=16), intent(in) :: b
    !
    !    ! auxiliaray
    !    integer :: size_a_org
    !
    !    if (allocated(a)) then
    !        size_a_org = size(a)
    !    else
    !        size_a_org = 0
    !    end if
    !
    !    call reallocate(a, 1)
    !
    !    a(size_a_org+1:) = b
    !end subroutine
end module toolbox
