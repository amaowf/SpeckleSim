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

module lib_mie_ms_solver_interface_helper_functions
    use libmath
    implicit none

    private

    public :: lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array
    public :: lib_mie_ms_solver_hf_add_list_list_cmplx_at_array
    public :: lib_mie_ms_solver_hf_get_list_list_cmplx_from_array

    public :: lib_mie_ms_solver_interface_hf_helper_functions

    contains

        ! Argument
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   b_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        ! Returns
        ! ----
        !   array: double complex, dimension(:)
        !
        subroutine lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array(a_nm, b_nm, element_no, n_range, array)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            double complex, dimension(:), intent(inout) :: array

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            double complex, dimension(:), allocatable :: buffer_array

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            call make_array(a_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            last = first + counter - 1
            array(first:last) = buffer_array

            call make_array(b_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            first = last + 1
            last = first + counter - 1
            array(first:last) = buffer_array

        end subroutine lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array

        ! Argument
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   b_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        ! Returns
        ! ----
        !   array: double complex, dimension(:)
        !
        subroutine lib_mie_ms_solver_hf_add_list_list_cmplx_at_array(a_nm, b_nm, element_no, n_range, array)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            double complex, dimension(:), intent(inout) :: array

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            double complex, dimension(:), allocatable :: buffer_array

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            call make_array(a_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            last = first + counter - 1
            array(first:last) =array(first:last) + buffer_array

            call make_array(b_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            first = last + 1
            last = first + counter - 1
            array(first:last) = array(first:last) + buffer_array

        end subroutine lib_mie_ms_solver_hf_add_list_list_cmplx_at_array

        ! Argument
        ! ----
        !   array: double complex, dimension(:)
        !       array of all a_nm and b_nm coefficients
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        subroutine lib_mie_ms_solver_hf_get_list_list_cmplx_from_array(array, element_no, n_range, a_nm, b_nm)
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: array
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: a_nm
            type(list_list_cmplx), intent(inout) :: b_nm

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            last = first + counter - 1
            call make_list(array(first:last), &
                           n_range(1), n_range(2)-n_range(1)+1, &
                           a_nm)
            call remove_zeros(a_nm)

            first = last + 1
            last = first + counter - 1
            call make_list(array(first:last), &
                           n_range(1), n_range(2)-n_range(1)+1, &
                           b_nm)
            call remove_zeros(b_nm)

        end subroutine lib_mie_ms_solver_hf_get_list_list_cmplx_from_array

        function lib_mie_ms_solver_interface_hf_helper_functions() result(error)
            implicit none
            ! dummy
            integer :: error

            error = 0

            if (.not. test_mie_ms_solver_hf()) error = error + 1


            contains

            function test_mie_ms_solver_hf() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: n
                integer :: m

                integer :: no
                integer :: counter

                type(list_list_cmplx), dimension(:), allocatable :: a_nm
                type(list_list_cmplx), dimension(:), allocatable :: b_nm

                double complex, dimension(:), allocatable :: array

                integer, dimension(2) :: n_range

                n_range = (/ 1, 3 /)

                counter = ( 1 + n_range(2))**2 - n_range(1)**2
                no = 3

                allocate(array(2*counter * no))
                allocate(a_nm(no))
                allocate(b_nm(no))

                do i=1, no
                    call init_list(a_nm(i), n_range(1), n_range(2) - n_range(1) + 1, dcmplx(i,0))
                    call init_list(b_nm(i), n_range(1), n_range(2) - n_range(1) + 1, dcmplx(i,0))

                    call lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array(a_nm(i), b_nm(i), i, n_range, array)
                    call lib_mie_ms_solver_hf_add_list_list_cmplx_at_array(a_nm(i), b_nm(i), i, n_range, array)

                    call lib_mie_ms_solver_hf_get_list_list_cmplx_from_array(array, i, n_range, a_nm(i), b_nm(i))
                end do


                rv = .true.
                print *, "test_mie_ms_solver_hf"
                do i = 1, no
                    print *, "  i = ", i
                    print *, "  a_nm"
                    do n = n_range(1) , n_range(2)
                        print *, "    n = ", n
                        do m = -n, n
                            if (a_nm(i)%item(n)%item(m) .eq. dcmplx(2 * i, 0)) then
                                print *, "    m = ", m, ": OK"
                            else
                                print *, "    m = ", m, ": FAILED"
                                rv = .false.
                            end if
                        end do
                    end do

                    print *, "  b_nm"
                    do n = n_range(1) , n_range(2)
                        print *, "    n = ", n
                        do m = -n, n
                            if (b_nm(i)%item(n)%item(m) .eq. dcmplx(2 * i, 0)) then
                                print *, "    m = ", m, ": OK"
                            else
                                print *, "    m = ", m, ": FAILED"
                                rv = .false.
                            end if
                        end do
                    end do
                end do

            end function test_mie_ms_solver_hf

        end function
end module lib_mie_ms_solver_interface_helper_functions
