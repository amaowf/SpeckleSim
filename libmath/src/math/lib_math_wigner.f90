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

module lib_math_wigner
    use fwigxjpf
    implicit none

    private

    ! fwigxjpf lib functions
    public :: fwig_table_init
    public :: fwig_temp_init
    public :: fwig_thread_temp_init

    public :: fwig_temp_free
    public :: fwig_table_free

    ! functions
    public :: lib_math_wigner_3j

    public :: lib_math_wigner_test_functions

    contains

        ! Argument
        ! ----
        !
        function lib_math_wigner_3j(j1, j2, j3, m1, m2, m3, double_integer) result(rv)
            implicit none
            ! dummy
            integer, intent(in) :: j1
            integer, intent(in) :: j2
            integer, intent(in) :: j3
            integer, intent(in) :: m1
            integer, intent(in) :: m2
            integer, intent(in) :: m3
            logical, intent(in), optional :: double_integer

            real(kind=8) :: rv


            if (present(double_integer)) then
                if (double_integer) then
                    rv = fwig3jj(j1, j2, j3, &
                                 m1, m2, m3)
                else
                    rv = fwig3jj(2_4*j1, 2_4*j2, 2_4*j3, &
                                 2_4*m1, 2_4*m2, 2_4*m3)
                end if
            else
                rv = fwig3jj(2_4*j1, 2_4*j2, 2_4*j3, &
                             2_4*m1, 2_4*m2, 2_4*m3)
            end if

        end function

        function lib_math_wigner_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! parameter
            real(kind=8), parameter :: ground_truth_e = 1D-14

            rv = 0

            if (.not. test_lib_math_wigner_3j()) then
                rv = rv + 1
            end if

            print *, "-------------lib_math_wigner_test_functions----------------"
            if (rv == 0) then
                print *, "lib_math_wigner_test_functions tests: OK"
            else
                print *, rv,"lib_math_wigner_test_functions test(s) FAILED"
            end if
            print *, "-----------------------------------------------------------"

            contains

                function test_lib_math_wigner_3j() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: n = 5

                    ! auxiliary
                    integer(kind=4) :: j1
                    integer(kind=4) :: j2
                    integer(kind=4), dimension(n) :: j3
                    integer(kind=4) :: m1
                    integer(kind=4) :: m2
                    integer(kind=4) :: m3
                    logical :: double_integer

                    real(kind=8), dimension(n) :: res
                    real(kind=8), dimension(n) :: ground_truth
                    real(kind=8) :: buffer

                    integer :: i
                    integer(kind=4) :: j_max


                    double_integer = .false.

                    ! Reference: Efficient Evaluation of Vector Translation Coefficients in Multiparticle Light-Scattering Theories
                    !            Table 5
                    j1 = 98
                    j2 = 115
                    j3 = (/ 213, 212, 211, 210, 209/)

                    m1 = -69
                    m2 = -100
                    m3 = 169

                    ! less precise (case: ground_truth_e = 10**-14)
                    ! -> error =~ -1.2 * 10**-14
!                    ground_truth(1) = 0.247807608975D-02
!                    ground_truth(2) = 0.693884894548D-02
!                    ground_truth(3) = 0.119894874087D-01
!                    ground_truth(4) = 0.137159251054D-01
!                    ground_truth(5) = 0.878978201272D-02

                    ! Values were generated with mathematica
                    !
                    ! source code:
                    !   >>> j1 = 98;
                    !   >>> j2 = 115;
                    !   >>> j3 = {213, 212, 211, 210, 209};
                    !   >>>
                    !   >>> m1 = -69;
                    !   >>> m2 = -100;
                    !   >>> m3 = 169;
                    !   >>>
                    !   >>> For[i = 1, i <= Length[j3], i ++,
                    !   >>>  erg = ThreeJSymbol[{j1, m1}, {j2, m2}, {Part[j3, i], m3}];
                    !   >>>  (*Print[Part[j3, i], "  : ", N[erg, 16]];*)
                    !   >>>
                    !   >>>  Print[Part[j3, i], "  : ", ScientificForm[N[erg], 16]];
                    !   >>> ]
                    ground_truth(1) = 2.478076089751141D-3
                    ground_truth(2) = 6.938848945479547D-3
                    ground_truth(3) = 1.198948740868428D-2
                    ground_truth(4) = 1.371592510538986D-2
                    ground_truth(5) = 8.78978201272337D-3


                    j_max = max(j1, j2, maxval(j3))
                    call fwig_table_init(2 * j_max, 3_4)
                    call fwig_temp_init(2 * j_max)

                    do i=1, n
                        res(i) = lib_math_wigner_3j(j1, j2, j3(i), m1, m2, m3, double_integer)
                    end do

                    rv = .true.
                    print *, "test_lib_math_wigner_3j:"
                    do i=1, n
                        buffer = res(i) - ground_truth(i)
                        if (abs(buffer) .gt. ground_truth_e) then
                            print *, "  ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  ", i, ": OK"
                        end if
                    end do

                end function test_lib_math_wigner_3j

        end function lib_math_wigner_test_functions

end module lib_math_wigner
