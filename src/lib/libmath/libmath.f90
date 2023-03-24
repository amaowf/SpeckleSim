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

module libmath
    use lib_math_public
    implicit none
!    public :: libmath_hello
!
!    interface libmath_hello
!        module procedure hello
!    end interface

    contains

        function test_lib_math() result (error_counter)
            implicit none
            ! dummy
            integer :: error_counter

            error_counter = 0

            error_counter = error_counter + test_lib_math_auxiliaries()
            error_counter = error_counter + lib_math_factorial_test_functions()
            error_counter = error_counter + lib_math_type_operator_test_functions()
            error_counter = error_counter + lib_math_bessel_test_functions()
            error_counter = error_counter + lib_math_legendre_test_functions()
            error_counter = error_counter + lib_math_hermite_test_functions()
            error_counter = error_counter + lib_math_solver_test_functions()
#ifdef __GFORTRAN__
!            error_counter = error_counter + lib_math_wigner_test_functions()
!            error_counter = error_counter + lib_math_vector_spherical_harmonics_test_functions()
#endif
            error_counter = error_counter + lib_math_gamma_distribution_test_functions()
            error_counter = error_counter + lib_math_set_theory_test_functions()
            error_counter = error_counter + lib_math_triangle_test_functions()

            print *, "-------------test_lib_math----------------"
            if (error_counter == 0) then
                print *, "test_lib_math tests: OK"
            else
                print *, error_counter,"test_lib_math test(s) FAILED"
            end if
            print *, "------------------------------------------"
        end function

end module
