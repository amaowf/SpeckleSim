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

module lib_math_auxiliaries
    implicit none

    private

    public :: isinf
    public :: isinf_pos
    public :: isinf_neg
    public :: get_inf_pos
    public :: get_inf_neg

    public :: test_lib_math_auxiliaries

    interface isinf
        module procedure lib_math_real_4_isinf
        module procedure lib_math_real_8_isinf
    end interface

    interface isinf_pos
        module procedure lib_math_real_4_is_inf_pos
        module procedure lib_math_real_8_is_inf_pos
    end interface

    interface isinf_neg
        module procedure lib_math_real_4_is_inf_neg
        module procedure lib_math_real_8_is_inf_neg
    end interface

    interface get_inf_pos
!        module procedure lib_math_real_4_get_inf_pos
        module procedure lib_math_real_8_get_inf_pos
    end interface

    interface get_inf_neg
!        module procedure lib_math_real_4_get_inf_neg
        module procedure lib_math_real_8_get_inf_neg
    end interface

    contains

    function lib_math_real_4_isinf(value) result (rv)
        implicit none
        ! dummy
        real(kind=4) :: value

        logical :: rv

        ! auxiliary
        real(kind=4) :: inf

        inf = huge(inf)

        if (value .gt. inf &
            .or. value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_4_is_inf_pos(value) result (rv)
        implicit none
        ! dummy
        real(kind=4) :: value

        logical :: rv

        ! auxiliary
        real(kind=4) :: inf

        inf = huge(inf)

        if (value .gt. inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_4_is_inf_neg(value) result (rv)
        implicit none
        ! dummy
        real(kind=4) :: value

        logical :: rv

        ! auxiliary
        real(kind=4) :: inf

        inf = huge(inf)

        if (value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

!    function lib_math_real_4_get_inf_pos() result (rv)
!        implicit none
!        ! dummy
!        real(kind=4) :: rv
!
!        rv = huge(rv) + 1
!    end function
!
!    function lib_math_real_4_get_inf_neg() result (rv)
!        implicit none
!        ! dummy
!        real(kind=4) :: rv
!
!        rv = -huge(rv) - 1
!    end function

    function lib_math_real_8_isinf(value) result (rv)
        implicit none
        ! dummy
        real(kind=8) :: value

        logical :: rv

        ! auxiliary
        real(kind=8) :: inf

        inf = huge(inf)

        if (value .gt. inf &
            .or. value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_8_is_inf_pos(value) result (rv)
        implicit none
        ! dummy
        real(kind=8) :: value

        logical :: rv

        ! auxiliary
        real(kind=8) :: inf

        inf = huge(inf)

        if (value .gt. inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_8_is_inf_neg(value) result (rv)
        implicit none
        ! dummy
        real(kind=8) :: value

        logical :: rv

        ! auxiliary
        real(kind=8) :: inf

        inf = huge(inf)

        if (value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_8_get_inf_pos() result (rv)
        use ieee_arithmetic
        implicit none
        ! dummy
        real(kind=8) :: rv

        rv = ieee_value(rv,  ieee_positive_inf)
    end function

    function lib_math_real_8_get_inf_neg() result (rv)
        use ieee_arithmetic
        implicit none
        ! dummy
        real(kind=8) :: rv

        rv = ieee_value(rv,  ieee_negative_inf)
    end function

    function test_lib_math_auxiliaries() result (rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        if (.not. test_lib_math_real_4_isinf()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_real_8_isinf()) then
            rv = rv + 1
        end if


        print *, "-------------test_lib_math_auxiliaries----------------"
        if (rv == 0) then
            print *, "test_lib_math_auxiliaries tests: OK"
        else
            print *, rv,"test_lib_math_auxiliaries test(s) FAILED"
        end if
        print *, "------------------------------------------------------"

        contains

            function test_lib_math_real_4_isinf() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 3

                ! auxiliary
                integer :: i
                real(kind=4), dimension(d) :: value
                logical, dimension(d) :: res
                logical, dimension(d) :: ground_truth_res

                value(1) = huge(value)
                value(1) = value(1) + value(1)
                value(2) = -value(1)
                value(3) = huge(value)

                ground_truth_res = (/ .true., .true., .false. /)

                do i=1, d
                    res(i) = lib_math_real_4_isinf(value(i))
                end do

                rv = .true.
                print *, "test_lib_math_real_4_isinf:"
                do i=1, d
                    if (res(i) .neqv. ground_truth_res(i)) then
                        rv = .false.
                        print *, "  ", i, ": FAILED"
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_real_4_isinf

            function test_lib_math_real_8_isinf() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 3

                ! auxiliary
                integer :: i
                real(kind=8), dimension(d) :: value
                logical, dimension(d) :: res
                logical, dimension(d) :: ground_truth_res

                value(1) = huge(value)
                value(1) = value(1) + value(1)
                value(2) = -value(1)
                value(3) = huge(value)

                ground_truth_res = (/ .true., .true., .false. /)

                do i=1, d
                    res(i) = lib_math_real_8_isinf(value(i))
                end do

                rv = .true.
                print *, "test_lib_math_real_8_isinf:"
                do i=1, d
                    if (res(i) .neqv. ground_truth_res(i)) then
                        rv = .false.
                        print *, "  ", i, ": FAILED"
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_real_8_isinf

    end function test_lib_math_auxiliaries

end module lib_math_auxiliaries
