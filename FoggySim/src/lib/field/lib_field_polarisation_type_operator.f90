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

module lib_field_polarisation_type_operator
    use libmath
    use lib_field_polarisation_type
    implicit none

    private

    public :: operator(*)

    public :: lib_field_polarisation_operator_test_functions


    interface operator(*)
        module procedure scalar_field_dot_jones_vector
        module procedure jones_vector_dot_scalar_field

        module procedure field_dot_jones_vector
        module procedure jones_vector_dot_field

        module procedure jones_matrix_dot_vector
        module procedure jones_matrix_dot_matrix
    end interface

    contains

        function scalar_field_dot_jones_vector(lhs, rhs) result(rv)
            implicit none
            ! dummy
            double complex, intent(in) :: lhs
            type(jones_vector_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = dcmplx(0, 0)

        end function

        function jones_vector_dot_scalar_field(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(jones_vector_type), intent(in) :: lhs
            double complex, intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv = scalar_field_dot_jones_vector(rhs, lhs)

        end function

        function field_dot_jones_vector(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type(jones_vector_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs%x
            rv%y = lhs%z * rhs%y
            rv%z = dcmplx(0, 0)

        end function

        function jones_vector_dot_field(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(jones_vector_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv = field_dot_jones_vector(rhs, lhs)

        end function

        function jones_matrix_dot_vector(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(jones_matrix_type), intent(in) :: lhs
            type(jones_vector_type), intent(in) :: rhs

            type(jones_vector_type) :: rv

            rv%x = lhs%m_11 * rhs%x + lhs%m_12 * rhs%y
            rv%y = lhs%m_21 * rhs%x + lhs%m_22 * rhs%y
        end function

        function jones_matrix_dot_matrix(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type(jones_matrix_type), intent(in) :: lhs
            type(jones_matrix_type), intent(in) :: rhs

            type(jones_matrix_type) :: rv

            rv%m_11 = lhs%m_11 * rhs%m_11 + lhs%m_12 * rhs%m_21
            rv%m_12 = lhs%m_11 * rhs%m_12 + lhs%m_12 * rhs%m_22

            rv%m_21 = lhs%m_21 * rhs%m_11 + lhs%m_22 * rhs%m_21
            rv%m_22 = lhs%m_21 * rhs%m_12 + lhs%m_22 * rhs%m_22
        end function

        function lib_field_polarisation_operator_test_functions() result(rv)
            use libmath
            implicit none
            ! dummy
            integer :: rv

            double precision, parameter :: ground_truth_e = 1D-14

            rv = 0

            if (.not. test_jones_matrix_dot_vector()) rv = rv + 1

            if (rv == 0) then
                print *, "lib_field_polarisation_operator_test_functions tests: OK"
            else
                print *, rv,"lib_field_polarisation_operator_test_functions test(s) FAILED"
            end if

            contains

                function test_jones_matrix_dot_vector() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    type(jones_matrix_type) :: m
                    type(jones_vector_type) :: j
                    type(jones_vector_type) :: j_res
                    type(jones_vector_type) :: ground_truth_j_res

                    double complex :: buffer

                    j = lib_field_polarisation_jones_vector_get_linear_rot(PI / 4d0)
                    m = lib_field_polarisation_jones_matrix_get_x_polariser()

                    j_res = m * j

                    ground_truth_j_res%x = dcmplx(1d0 / sqrt(2d0), 0)
                    ground_truth_j_res%y = dcmplx(0, 0)


                    rv = .true.
                    print *, "test_jones_matrix_dot_vector:"

                    buffer = ground_truth_j_res%x - j_res%x

                    if(real(buffer) .lt. ground_truth_e) then
                        print *, "  Re[j_x]: OK"
                    else
                        print *, "  Re[j_x]: FAILED"
                        print *, "        Re[j] = ", real(j_res%x)
                        print *, "       Re[GT] = ", real(ground_truth_j_res%x)
                        rv = .false.
                    end if

                    if(aimag(buffer) .lt. ground_truth_e) then
                        print *, "  Im[j_x]: OK"
                    else
                        print *, "  Im[j_x]: FAILED"
                        print *, "        Re[j] = ", aimag(j_res%x)
                        print *, "       Re[GT] = ", aimag(ground_truth_j_res%x)
                        rv = .false.
                    end if

                    buffer = ground_truth_j_res%y - j_res%y

                    if(real(buffer) .lt. ground_truth_e) then
                        print *, "  Re[j_y]: OK"
                    else
                        print *, "  Re[j_y]: FAILED"
                        print *, "        Re[j] = ", real(j_res%y)
                        print *, "       Re[GT] = ", real(ground_truth_j_res%y)
                        rv = .false.
                    end if

                    if(aimag(buffer) .lt. ground_truth_e) then
                        print *, "  Im[j_y]: OK"
                    else
                        print *, "  Im[j_y]: FAILED"
                        print *, "        Re[j] = ", aimag(j_res%y)
                        print *, "       Re[GT] = ", aimag(ground_truth_j_res%y)
                        rv = .false.
                    end if

                end function test_jones_matrix_dot_vector

        end function
end module lib_field_polarisation_type_operator
