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


module toolbox
    implicit none
    public :: reallocate

    interface reallocate
        module procedure reallocate_1d_integer_1
        module procedure reallocate_1d_integer_2
        module procedure reallocate_1d_integer_4
        module procedure reallocate_1d_integer_8
#ifdef __GFORTRAN__
        module procedure reallocate_1d_integer_16
#endif
    end interface

    interface concatenate
        module procedure lib_tree_hf_concatenate_1d_integer_1_array
        module procedure lib_tree_hf_concatenate_1d_integer_2_array
        module procedure lib_tree_hf_concatenate_1d_integer_4_array
        module procedure lib_tree_hf_concatenate_1d_integer_8_array
#ifdef __GFORTRAN__
        module procedure lib_tree_hf_concatenate_1d_integer_16_array
#endif
        module procedure lib_tree_hf_concatenate_1d_integer_1_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_2_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_4_array_single
        module procedure lib_tree_hf_concatenate_1d_integer_8_array_single
#ifdef __GFORTRAN__
        module procedure lib_tree_hf_concatenate_1d_integer_16_array_single
#endif
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

#ifdef __GFORTRAN__
    subroutine reallocate_1d_integer_16(a,n)
        implicit none

        integer(kind=16),dimension(:),allocatable,intent(inout) :: a
        integer(kind=16),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if ( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate( a(n) )
        end if

    end subroutine
#endif

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

#ifdef __GFORTRAN__
    subroutine lib_tree_hf_concatenate_1d_integer_16_array(a, b)
        implicit none
        integer(kind=16), dimension(:), allocatable, intent(inout) :: a
        integer(kind=16), dimension(:), allocatable, intent(in) :: b

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
#endif

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

#ifdef __GFORTRAN__
    subroutine lib_tree_hf_concatenate_1d_integer_16_array_single(a, b)
        implicit none
        ! dummy
        integer(kind=16), dimension(:), allocatable, intent(inout) :: a
        integer(kind=16), intent(in) :: b

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
#endif

    !               interpolation
    !                  |
    !  y ^             v
    !    |   o   o   o Y o   o <-- data point
    !    |
    !    |   o   o   o Y o   o
    !    |
    !    |   o   o   o Y o   o
    !    |
    !    |---o---o---o-X-o---o-> x
    !                  ^
    !                  x-value
    !
    ! Argument
    ! ----
    !   data: double precision, dimension(:,:)
    !       2-dimensional data set
    !   x_column: integer
    !       x column selector
    !       IMPORTANT: "data" is sorted in ascending order in relation to this column
    !   x_value: double precsison
    !       x value
    !
    ! Returns
    ! ----
    !   rv: double precision, dimension(:)
    !       interpolated values
    subroutine data_interpolation(data, x_column, x_value, rv)
        implicit none
        ! dummy
        double precision, dimension(:,:), intent(in) :: data
        integer, intent(in) :: x_column
        double precision, intent(in) :: x_value

        double precision, dimension(size(data,2)) :: rv

        ! auxiliary
        integer :: row
        integer :: column
        double precision, dimension(size(data,2)) :: value_1
        double precision, dimension(size(data,2)) :: value_2
        double precision, dimension(2) :: point_1
        double precision, dimension(2) :: point_2

        if (x_column .ge. lbound(data, 2) &
            .and. x_column .le. ubound(data, 2)) then

            do row = lbound(data, 1)+1, ubound(data, 1)-1
                if (x_value .lt. data(row, x_column)) then
!                    rv = data(row-1, :)
!                    exit

                    value_1 = data(row-1, :)
                    value_2 = data(row, :)

                    point_1(1) = value_1(x_column)
                    point_2(1) = value_2(x_column)
                    do column=1, size(data,2)
                        point_1(2) = value_1(column)
                        point_2(2) = value_2(column)

                        rv(column) = linear_interpolation(point_1, point_2, x_value)
                    end do
                    exit
                end if
            end do
        end if

    end subroutine

    ! formula f(x) = y = mx + b
    !
    ! Argument
    ! ----
    !   point_1: double precision, dimension(2)
    !       x and y value of point_1
    !   point_2: double precision, dimension(2)
    !       x and y value of point_2
    !   x_value: double precision
    !       x value of the evaluation point
    !
    ! Returns
    ! ----
    !   y: double precision
    !       calculated y value of the evaluation point
    !
    function linear_interpolation(point_1, point_2, x_value) result (y)
        implicit none
        ! dummy
        double precision, dimension(2), intent(in) :: point_1
        double precision, dimension(2), intent(in) :: point_2
        double precision, intent(in) :: x_value

        double precision :: y

        ! auxiliaray
        double precision :: m
        double precision :: b

        m = ( point_2(2) - point_1(2) ) / (point_2(1) - point_1(1))
        b = point_2(2) - m * point_2(1)

        y = m * x_value +b

    end function linear_interpolation
end module toolbox
