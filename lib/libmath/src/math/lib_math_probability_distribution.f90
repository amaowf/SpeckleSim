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
! Created on Thu Mai 12 10:00:51 2020
!
! @author: Max Daiber-Huppert
!

module lib_math_probability_distribution
    use lib_math_type
    implicit none

    private

    public :: lib_math_gamma_distribution
    public :: lib_math_gamma_distribution_deirmedjian

    public :: lib_math_gamma_distribution_test_functions

    interface lib_math_gamma_distribution
        module procedure lib_math_gamma_distribution_deirmedjian
    end interface

    contains

        ! Deirmendjian modified gamma distribution
        !
        ! Argument
        ! ----
        !   parameter: lib_math_gamma_distribution_deirmedjian_type
        !       parameter of the probability distribution
        !   x: real(kind=lib_math_type_kind)
        !       Argument of the gamma distribution function
        !
        ! Returns
        ! ----
        !   rv: real(kind=lib_math_type_kind)
        !       Probability value at x
        !
        ! Reference: Comparison of laser beam propagation at 785 nm and 1550 nm in fog and haze for optical wireless communications,
        !            Kim, Isaac I. and {McArthur}, Bruce and Korevaar, Eric J.,
        !            eq. (8)
        function lib_math_gamma_distribution_deirmedjian(parameter, x) result(rv)
            implicit none
            ! dummy
            type(lib_math_gamma_distribution_deirmedjian_type), intent(in) :: parameter
            real(kind=lib_math_type_kind), intent(in) :: x

            real(kind=lib_math_type_kind) :: rv

            ! auxiliary

            rv = parameter%a * x**parameter%alpha * exp(-parameter%b * x**parameter%gamma)
        end function

        function lib_math_gamma_distribution_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliary
            double precision, parameter :: ground_truth_error = 1d-7

            rv = 0

            if (.not. test_lib_math_gamma_distribution_deirmedjian()) rv = rv + 1

            contains

            function test_lib_math_gamma_distribution_deirmedjian() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                double precision :: x
                double precision :: res
                double precision :: res_ground_truth
                type(lib_math_gamma_distribution_deirmedjian_type) :: para
                double precision :: buffer

                para%a = 0.027d0
                para%alpha = 3d0
                para%gamma = 1d0
                para%b = 0.3d0

                x = 10.89d0

                res_ground_truth = 1.32925202638267

                res = lib_math_gamma_distribution_deirmedjian(para, x)

                rv = .true.
                buffer = res_ground_truth - res
                if (abs(buffer) .lt. ground_truth_error) then
                    print *, "test_lib_math_gamma_distribution_deirmedjian: OK"
                else
                    print *, "test_lib_math_gamma_distribution_deirmedjian: ERROR"
                    print *, "  res = ", res
                    print *, "   gt = ", res_ground_truth
                end if

            end function test_lib_math_gamma_distribution_deirmedjian

        end function lib_math_gamma_distribution_test_functions
end module lib_math_probability_distribution
