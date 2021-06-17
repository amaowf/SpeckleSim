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

module lib_math_solver
    use solver_conjugate_gradient_method
    implicit none

    private

    ! public functions
    public :: lib_math_solver_conjugate_gradient_method

    public :: lib_math_solver_test_functions

    contains

        ! formula: Ax = b
        !
        ! Arguments
        ! ----
        !   x: double complex, dimension(:)
        !       initial guess
        !   matrix_a
        !
        function lib_math_solver_conjugate_gradient_method(matrix_a, x_initial, b, result_precision, max_iterations) result(x)
            implicit none
            ! dummy
            double complex, dimension(:, :), allocatable, intent(inout) :: matrix_a
            double complex, dimension(:), allocatable, intent(in) :: x_initial
            double complex, dimension(:), allocatable, intent(inout) :: b

            double precision, intent(in), optional :: result_precision
            integer, intent(in), optional :: max_iterations

            double complex, dimension(:), allocatable :: x

            ! auxiliary
            double complex :: res
            double complex, dimension(:), allocatable :: m_x
            double precision :: m_result_precision
            integer :: m_max_iterations

            if (present(result_precision)) then
                m_result_precision = result_precision
            else
                m_result_precision = 10D-6
            end if

            if (present(max_iterations)) then
                m_max_iterations = max_iterations
            else
                m_max_iterations = 100
            end if

            allocate( m_x(lbound(x_initial, 1):ubound(x_initial, 1)) )
            m_x = x_initial

            res = dm_bcg(m_x, m_result_precision, m_max_iterations,matrix_a,b)

            if (aimag(res) .gt. m_result_precision) then
                print *, "lib_math_solver_conjugate_gradient_method: WARNING"
                print *, "  res: ", aimag(res)
                print *, "  iterations: ", real(res)
                print *, "  max. iterations: ", m_max_iterations

            end if
            call move_alloc(m_x, x)

        end function


        function lib_math_solver_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! parameter
            double precision :: ground_truth_error = 10D-12

            rv = 0

            if (.not. test_lib_math_solver_conjugate_gradient_method()) then
                rv = rv + 1
            end if


            print *, "-------------lib_math_solver_test_functions----------------"
            if (rv == 0) then
                print *, "lib_math_solver_test_functions tests: OK"
            else
                print *, rv,"lib_math_solver_test_functions test(s) FAILED"
            end if
            print *, "-----------------------------------------------------------"

            contains

            function test_lib_math_solver_conjugate_gradient_method() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                double precision :: buffer
                double complex, dimension(:, :), allocatable :: matrix_a
                double complex, dimension(:), allocatable :: x_initial
                double complex, dimension(:), allocatable :: x
                double complex, dimension(:), allocatable :: b
                double complex, dimension(:), allocatable :: b_test

                allocate(matrix_a(2,2))
                allocate(x_initial(2))
                allocate(b(2))
                allocate(b_test(2))

                matrix_a = reshape((/dcmplx(4, 2), dcmplx(1, 0.2), &
                                     dcmplx(1, 4), dcmplx(3, 1)/), shape(matrix_a))

                b = (/ dcmplx(1, 0), dcmplx(2, 0) /)
                x_initial  = (/ dcmplx(0.1,0), dcmplx(0.5,0) /)

                x = lib_math_solver_conjugate_gradient_method(matrix_a, x_initial, b)

                b_test = matmul(matrix_a, x)

                rv = .true.

                print *, "test_lib_math_solver_conjugate_gradient_method"
                do i=1, size(b)
                    buffer = real(b(i)) - real(b_test(i))
                    if (buffer .gt. ground_truth_error) then
                        print *, "  ", i, "  real: FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, "  real: OK"
                    end if

                    buffer = aimag(b(i)) - aimag(b_test(i))
                    if (buffer .gt. ground_truth_error) then
                        print *, "  ", i, "  cmplx: FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, "  cmplx: OK"
                    end if
                end do

            end function test_lib_math_solver_conjugate_gradient_method

        end function lib_math_solver_test_functions

end module lib_math_solver

