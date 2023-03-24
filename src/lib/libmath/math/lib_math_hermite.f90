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

module lib_math_hermite
    use fsu_hermite_polynomial
    use fsu_test_values
    implicit none

    private

    public :: lib_math_hermite_polynomial

    public :: lib_math_hermite_test_functions


    interface lib_math_hermite_polynomial
        module procedure lib_math_hermite_polynomial_one_x
        module procedure lib_math_hermite_polynomial_array_x
    end interface

    contains

        ! Argument
        ! ----
        !   x: double precision
        !       evaluation point
        !   fnu: integer
        !       lowest order of the polynomial to compute, fnu >= 0
        !   n: integer
        !       number of polynomial, n >= 1
        !
        ! Returns
        ! ----
        !   h: double precision(fnu:fnu+x-1)
        !       values of the Hermite polynomial at the point x
        !
        subroutine lib_math_hermite_polynomial_one_x(x, fnu, n, h)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            double precision, dimension(:), allocatable :: h

            ! auxiliary
            integer ( kind = 4 ) m_n
            real ( kind = 8 ), dimension(1,0:fnu+n-1) :: m_p

            m_n = fnu + n - 1

            call h_polynomial_value ( 1, m_n, (/ x /), m_p )

            if (allocated(h)) deallocate(h)

            allocate(h(fnu:fnu+n-1))

            h = m_p(1, fnu:fnu+n-1)

        end subroutine lib_math_hermite_polynomial_one_x

        ! Argument
        ! ----
        !   x: double precision, dimension(:)
        !       evaluation points
        !   fnu: integer
        !       lowest order of the polynomial to compute, fnu >= 0
        !   n: integer
        !       number of polynomial, n >= 1
        !
        ! Returns
        ! ----
        !   h: double precision(:,fnu:fnu+x-1)
        !       values of the Hermite polynomial at the point x
        !
        subroutine lib_math_hermite_polynomial_array_x(x, fnu, n, h)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: x
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            double precision, dimension(:,:), allocatable :: h

            ! auxiliary
            integer ( kind = 4 ) m_m
            integer ( kind = 4 ) m_n
            real ( kind = 8 ), dimension(:,:), allocatable :: m_p

            allocate(m_p(lbound(x, 1):ubound(x, 1),0:fnu+n-1))

            m_m = size(x)

            if (allocated(h)) deallocate(h)
            allocate(h(lbound(x,1):ubound(x,1), fnu:fnu+n-1))

            m_n = fnu + n - 1
            if (m_n .eq. 0) then
                h = 1
            else
                call h_polynomial_value ( m_m, m_n, x, m_p )

                h = m_p(lbound(x, 1):ubound(x, 1), fnu:fnu+n-1)
            end if

        end subroutine lib_math_hermite_polynomial_array_x

        function lib_math_hermite_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 1d-15

            rv = 0

            if (.not. test_lib_math_hermite_polynomial_one_x()) rv = rv + 1


            print *, "-------------lib_math_hermite_test_functions-------------"
            if (rv == 0) then
                print *, "lib_math_hermite_test_functions tests: OK"
            else
                print *, rv,"lib_math_hermite_test_functions test(s) FAILED"
            end if
            print *, "----------------------------------------------------------"

            contains

                function test_lib_math_hermite_polynomial_one_x() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: n_data
                    integer :: n
                    double precision :: x
                    double precision, dimension(:), allocatable :: fx
                    double precision :: ground_truth_fx

                    double precision :: buffer


                    rv = .true.
                    print *, "test_lib_math_hermite_polynomial_one_x"

                    n_data = 0
                    do
                        call hermite_poly_phys_values ( n_data, n, x, ground_truth_fx )

                        if (n_data .eq. 0) then
                            exit
                        end if

                        call lib_math_hermite_polynomial_one_x(x, n, 1, fx)

                        buffer = ground_truth_fx - fx(lbound(fx,1))

                        if (abs(buffer) .le. ground_truth_e) then
                            print *, "  ", n_data, " :  OK"
                        else
                            print *, "  ", n_data, " :  FAILED"
                            print *, "    fx = ", fx
                            print *, "    GT = ", ground_truth_fx
                            rv = .false.
                        end if

                    end do

                end function test_lib_math_hermite_polynomial_one_x

        end function

end module lib_math_hermite
