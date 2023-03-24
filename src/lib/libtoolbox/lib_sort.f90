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

#define _UINDEX_BYTES_ 8

module lib_sort
    implicit none

    private

    public :: lib_sort_hpsort
    public :: lib_sort_hpsort_integer

    public :: lib_sort_test_functions

    integer(kind=1), parameter, public :: LIB_HPSORT_INTEGER_KIND = _UINDEX_BYTES_
    integer(kind=1), parameter, public :: LIB_HPSORT_OLD_POS_INTEGER_KIND = 4

    contains

    subroutine lib_sort_hpsort_integer(n, ra, ra_old_position)
        INTEGER(kind=LIB_HPSORT_OLD_POS_INTEGER_KIND) n
        integer(kind=LIB_HPSORT_INTEGER_KIND), intent(inout) :: ra(n)
        integer(kind=LIB_HPSORT_OLD_POS_INTEGER_KIND), intent(inout) :: ra_old_position(n)
        !Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
        !input; ra is replaced on output by its sorted rearrangement.
        INTEGER(kind=LIB_HPSORT_INTEGER_KIND) i,ir,j,l
        integer(kind=LIB_HPSORT_INTEGER_KIND) rra

        integer(kind=LIB_HPSORT_OLD_POS_INTEGER_KIND) :: buffer_old_position

        do i=1, n
            ra_old_position(i) = int(i, kind=LIB_HPSORT_OLD_POS_INTEGER_KIND)
        end do

        if (n.lt.2) return
        !The index l will be decremented from its initial value down to 1 during the “hiring” (heap
        !creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
        !down to 1 during the “retirement-and-promotion” (heap selection) phase.
        l=n/2+1
        ir=n
        10 continue
            if(l.gt.1)then      !Still in hiring phase.
                l=l-1
                rra=ra(l)
                buffer_old_position = ra_old_position(l)
            else                !In retirement-and-promotion phase.
                rra=ra(ir)      !Clear a space at end of array.
                buffer_old_position = ra_old_position(ir)
                ra(ir)=ra(1)    !Retire the top of the heap into it.
                ra_old_position(ir) = ra_old_position(1)
                ir=ir-1         !Decrease the size of the corporation.
                if(ir.eq.1)then !Done with the last promotion.
                    ra(1)=rra   !The least competent worker of all!
                    ra_old_position(1) = buffer_old_position
                    return
                endif
            endif
            i=l                 !Whether in the hiring phase or promotion phase, we here
            j=l+l               !set up to sift down element rra to its proper level.
        20 if(j.le.ir)then     !“Do while j.le.ir:”
                if(j.lt.ir)then
                    if(ra(j).lt.ra(j+1))j=j+1   !Compare to the better underling.
                endif
                if(rra.lt.ra(j))then            !Demote rra.
                    ra(i)=ra(j)
                    ra_old_position(i) = ra_old_position(j)
                    i=j
                    j=j+j
                else                            !This is rra’s level. Set j to terminate the sift-down.
                    j=ir+1
                endif
                goto 20
            endif
            ra(i)=rra                           !Put rra into its slot.
            ra_old_position(i) = buffer_old_position
        goto 10
    end subroutine lib_sort_hpsort_integer

    subroutine lib_sort_hpsort(n, ra, ra_old_position)
        INTEGER n
        REAL ra(n)
        integer ra_old_position(n)
        !Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
        !input; ra is replaced on output by its sorted rearrangement.
        INTEGER i,ir,j,l
        REAL rra

        integer :: buffer_old_position

        do i=1, n
            ra_old_position(i) = i
        end do

        if (n.lt.2) return
        !The index l will be decremented from its initial value down to 1 during the “hiring” (heap
        !creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
        !down to 1 during the “retirement-and-promotion” (heap selection) phase.
        l=n/2+1
        ir=n
        10 continue
            if(l.gt.1)then      !Still in hiring phase.
                l=l-1
                rra=ra(l)
                buffer_old_position = ra_old_position(l)
            else                !In retirement-and-promotion phase.
                rra=ra(ir)      !Clear a space at end of array.
                buffer_old_position = ra_old_position(ir)
                ra(ir)=ra(1)    !Retire the top of the heap into it.
                ra_old_position(ir) = ra_old_position(1)
                ir=ir-1         !Decrease the size of the corporation.
                if(ir.eq.1)then !Done with the last promotion.
                    ra(1)=rra   !The least competent worker of all!
                    ra_old_position(1) = buffer_old_position
                    return
                endif
            endif
            i=l                 !Whether in the hiring phase or promotion phase, we here
            j=l+l               !set up to sift down element rra to its proper level.
        20 if(j.le.ir)then     !“Do while j.le.ir:”
                if(j.lt.ir)then
                    if(ra(j).lt.ra(j+1))j=j+1   !Compare to the better underling.
                endif
                if(rra.lt.ra(j))then            !Demote rra.
                    ra(i)=ra(j)
                    ra_old_position(i) = ra_old_position(j)
                    i=j
                    j=j+j
                else                            !This is rra’s level. Set j to terminate the sift-down.
                    j=ir+1
                endif
                goto 20
            endif
            ra(i)=rra                           !Put rra into its slot.
            ra_old_position(i) = buffer_old_position
        goto 10
    end subroutine lib_sort_hpsort

    ! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm.
    !
    ! Argument
    ! ----
    !   n: integer
    !       number of elements
    !   ra: integer, dimension(n)
    !       array to sort
    !
    ! Reference: Numerical Recipes in Fortran 77
    SUBROUTINE hpsort(n,ra)
        INTEGER n
        REAL ra(n)
        !Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
        !input; ra is replaced on output by its sorted rearrangement.
        INTEGER i,ir,j,l
        REAL rra
        if (n.lt.2) return
        !The index l will be decremented from its initial value down to 1 during the “hiring” (heap
        !creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
        !down to 1 during the “retirement-and-promotion” (heap selection) phase.
        l=n/2+1
        ir=n
        10 continue
            if(l.gt.1)then  !Still in hiring phase.
                l=l-1
                rra=ra(l)
            else                !In retirement-and-promotion phase.
                rra=ra(ir)      !Clear a space at end of array.
                ra(ir)=ra(1)    !Retire the top of the heap into it.
                ir=ir-1         !Decrease the size of the corporation.
                if(ir.eq.1)then !Done with the last promotion.
                    ra(1)=rra   !The least competent worker of all!
                    return
                endif
            endif
            i=l                 !Whether in the hiring phase or promotion phase, we here
            j=l+l               !set up to sift down element rra to its proper level.
        20 if(j.le.ir)then     !“Do while j.le.ir:”
                if(j.lt.ir)then
                    if(ra(j).lt.ra(j+1))j=j+1   !Compare to the better underling.
                endif
                if(rra.lt.ra(j))then            !Demote rra.
                    ra(i)=ra(j)
                    i=j
                    j=j+j
                else                            !This is rra’s level. Set j to terminate the sift-down.
                    j=ir+1
                endif
                goto 20
            endif
            ra(i)=rra                           !Put rra into its slot.
        goto 10
    END SUBROUTINE

    ! ----- test functions -----
    function lib_sort_test_functions() result(error_counter)
        implicit none
        ! dummy
        integer :: error_counter

        error_counter = 0

        if (.NOT. test_hpsort()) then
            error_counter = error_counter + 1
        end if
        if (.NOT. test_lib_sort_hpsort()) then
            error_counter = error_counter + 1
        end if
        if (.NOT. test_lib_sort_hpsort_integer()) then
            error_counter = error_counter + 1
        end if

        contains

        function test_lib_sort_hpsort_integer() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=LIB_HPSORT_OLD_POS_INTEGER_KIND), parameter :: n = 5
            integer(kind=LIB_HPSORT_INTEGER_KIND), dimension(n) :: array
            integer(kind=LIB_HPSORT_OLD_POS_INTEGER_KIND), dimension(n) :: array_old_position
            integer(kind=LIB_HPSORT_INTEGER_KIND), dimension(n) :: ground_truth_array_sorted

            integer :: wrong_elements
            integer :: i

            array = (/ 1, 5, 2, 7, 3 /)
            wrong_elements = 0
            ground_truth_array_sorted = (/ 1, 2, 3, 5, 7 /)

            call lib_sort_hpsort_integer(n, array, array_old_position)

            do i=1, n
                if (array(i) .ne. ground_truth_array_sorted(i)) then
                    wrong_elements = wrong_elements + 1
                end if
            end do

            if (wrong_elements .eq. 0) then
                rv = .true.
                print *, "test_lib_sort_hpsort_integer: ", "OK"
            else
                rv = .false.
                print *, "test_lib_sort_hpsort_integer: ", "FAILED"
            end if

        end function test_lib_sort_hpsort_integer

        function test_lib_sort_hpsort() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer, parameter :: n = 5
            real, dimension(n) :: array
            integer, dimension(n) :: array_old_position
            real, dimension(n) :: ground_truth_array_sorted

            integer :: wrong_elements
            integer :: i

            array = (/ 1, 5, 2, 7, 3 /)
            wrong_elements = 0
            ground_truth_array_sorted = (/ 1, 2, 3, 5, 7 /)

            call lib_sort_hpsort(n, array, array_old_position)

            do i=1, n
                if (array(i) .ne. ground_truth_array_sorted(i)) then
                    wrong_elements = wrong_elements + 1
                end if
            end do

            if (wrong_elements .eq. 0) then
                rv = .true.
                print *, "test_lib_sort_hpsort: ", "OK"
            else
                rv = .false.
                print *, "test_lib_sort_hpsort: ", "FAILED"
            end if

        end function test_lib_sort_hpsort

        function test_hpsort() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer, parameter :: n = 5
            real, dimension(n) :: array
            real, dimension(n) :: ground_truth_array_sorted

            integer :: wrong_elements
            integer :: i

            array = (/ 1, 5, 2, 7, 3 /)
            wrong_elements = 0
            ground_truth_array_sorted = (/ 1, 2, 3, 5, 7 /)

            call hpsort(n, array)

            do i=1, n
                if (array(i) .ne. ground_truth_array_sorted(i)) then
                    wrong_elements = wrong_elements + 1
                end if
            end do

            if (wrong_elements .eq. 0) then
                rv = .true.
                print *, "test(_lib_sort)_test_hpsort: ", "OK"
            else
                rv = .false.
                print *, "test(_lib_sort)_test_hpsort: ", "FAILED"
            end if

        end function test_hpsort

    end function lib_sort_test_functions

end module lib_sort
