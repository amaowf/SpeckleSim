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

module lib_field
    use libmath
    use file_io
    implicit none

    private

    public :: lib_field_export
    public :: lib_field_test_functions

    contains



        function lib_field_export(e_field, h_field, path) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:,:), intent(in) :: e_field
            type(cartesian_coordinate_cmplx_type), dimension(lbound(e_field, 1):ubound(e_field, 1), &
                                                             lbound(e_field, 2):ubound(e_field, 2)),&
                                                   intent(in) :: h_field
            character(len=*), intent(in) :: path

            logical :: rv

            ! auxiliary
            integer :: i
            integer :: ii

            integer :: u

            type(cartesian_coordinate_cmplx_type), dimension(:,:), allocatable :: h_field_conjg

            double precision, dimension(:, :), allocatable :: e_field_real_x
            double precision, dimension(:, :), allocatable :: e_field_real_y
            double precision, dimension(:, :), allocatable :: e_field_real_z

            double precision, dimension(:, :), allocatable :: h_field_real_x
            double precision, dimension(:, :), allocatable :: h_field_real_y
            double precision, dimension(:, :), allocatable :: h_field_real_z

            type(cartesian_coordinate_real_type), dimension(: , :), allocatable :: poynting
            type(cartesian_coordinate_cmplx_type) :: buffer_cartesian_cmplx
            double precision, dimension(:, :), allocatable :: poynting_abs

            allocate(e_field_real_x(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(poynting(lbound(e_field, 1):ubound(e_field, 1), &
                              lbound(e_field, 2):ubound(e_field, 2)))

            allocate(e_field_real_y(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(e_field_real_z(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(h_field_real_x(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(h_field_real_y(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(h_field_real_z(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))

            allocate(poynting_abs(lbound(e_field, 1):ubound(e_field, 1), &
                                  lbound(e_field, 2):ubound(e_field, 2)))

            allocate(h_field_conjg(lbound(e_field, 1):ubound(e_field, 1), &
                                   lbound(e_field, 2):ubound(e_field, 2)))

            do i = lbound(e_field, 1), ubound(e_field, 1)
                do ii = lbound(e_field, 2), ubound(e_field, 2)
                    e_field_real_x(i, ii) = real(e_field(i, ii)%x)
                    e_field_real_y(i, ii) = real(e_field(i, ii)%y)
                    e_field_real_z(i, ii) = real(e_field(i, ii)%z)

                    h_field_real_x(i, ii) = real(h_field(i, ii)%x)
                    h_field_real_y(i, ii) = real(h_field(i, ii)%y)
                    h_field_real_z(i, ii) = real(h_field(i, ii)%z)

                    ! calculate the Poynting vector: S = E x H*
                    ! eq. 43
                    h_field_conjg(i, ii)%x = conjg(h_field(i, ii)%x)
                    h_field_conjg(i, ii)%y = conjg(h_field(i, ii)%y)
                    h_field_conjg(i, ii)%z = conjg(h_field(i, ii)%z)

                    buffer_cartesian_cmplx = cross_product(e_field(i, ii), h_field_conjg(i, ii))

                    poynting(i, ii)%x = real(buffer_cartesian_cmplx%x) / 2D0
                    poynting(i, ii)%y = real(buffer_cartesian_cmplx%y) / 2D0
                    poynting(i, ii)%z = real(buffer_cartesian_cmplx%z) / 2D0

                    poynting_abs(i, ii) = abs(poynting(i, ii))
                end do
            end do

            ! --- wirte to PPM ---
            ! e field
            u = 99
            open(unit=u, file= trim(path) // "e_field_x.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_x)
            close(u)

            open(unit=u, file=trim(path) // "e_field_y.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "e_field_z.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "e_field_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "e_field_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "e_field_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_z, logarithmic = .true.)
            close(u)

            ! h field
            u = 99
            open(unit=u, file= trim(path) // "h_field_x.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_x)
            close(u)

            open(unit=u, file= trim(path) // "h_field_y.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "h_field_z.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "h_field_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "h_field_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "h_field_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_z, logarithmic = .true.)
            close(u)

            ! Poynting
            u = 99
            open(unit=u, file= trim(path) // "poynting_x.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%x)
            close(u)

            open(unit=u, file= trim(path) // "poynting_y.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%y)
            close(u)

            open(unit=u, file= trim(path) // "poynting_z.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "poynting_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "poynting_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "poynting_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%z, logarithmic = .true.)
            close(u)

            ! Poynting abs
            open(unit=u, file= trim(path) // "poynting_abs.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting_abs)
            close(u)
            open(unit=u, file= trim(path) // "poynting_abs_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting_abs, logarithmic = .true.)
            close(u)

            ! --- wirte to CSV ---
            ! e field
            u = 99
            open(unit=u, file= trim(path) // "e_field_x.csv", status='unknown')
            rv = write_csv(u, e_field_real_x)
            close(u)

            open(unit=u, file=trim(path) // "e_field_y.csv", status='unknown')
            rv = write_csv(u, e_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "e_field_z.csv", status='unknown')
            rv = write_csv(u, e_field_real_z)
            close(u)

            ! h field
            u = 99
            open(unit=u, file= trim(path) // "h_field_x.csv", status='unknown')
            rv = write_csv(u, h_field_real_x)
            close(u)

            open(unit=u, file= trim(path) // "h_field_y.csv", status='unknown')
            rv = write_csv(u, h_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "h_field_z.csv", status='unknown')
            rv = write_csv(u, h_field_real_z)
            close(u)

            ! Poynting
            u = 99
            open(unit=u, file= trim(path) // "poynting_x.csv", status='unknown')
            rv = write_csv(u, poynting(:,:)%x)
            close(u)

            open(unit=u, file= trim(path) // "poynting_y.csv", status='unknown')
            rv = write_csv(u, poynting(:,:)%y)
            close(u)

            open(unit=u, file= trim(path) // "poynting_z.csv", status='unknown')
            rv = write_csv(u, poynting(:,:)%z)
            close(u)

            ! Poynting abs
            open(unit=u, file= trim(path) // "poynting_abs.csv", status='unknown')
            rv = write_csv(u, poynting_abs)
            close(u)

        end function lib_field_export

        function lib_field_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if ( .not. test_lib_field_export() ) rv = rv + 1

            if (rv == 0) then
                print *, "lib_field_test_functions tests: OK"
            else
                print *, rv,"lib_field_test_functions test(s) FAILED"
            end if

            contains

            function test_lib_field_export() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii
                type(spherical_coordinate_cmplx_type), dimension(:,:), allocatable :: e_field_spherical
                type(spherical_coordinate_cmplx_type), dimension(:,:), allocatable :: h_field_spherical

                type(cartesian_coordinate_cmplx_type), dimension(:,:), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:,:), allocatable :: h_field
                character(len=20) :: path

                type(cartesian_coordinate_real_type) :: d

                path = "temp/test_field"

                allocate(e_field_spherical(5,5))
                allocate(h_field_spherical(5,5))

                allocate(e_field(5,5))
                allocate(h_field(5,5))

                e_field_spherical(:,:)%rho = dcmplx(0, 0)
                e_field_spherical(:,:)%theta = dcmplx(0, 0)
                e_field_spherical(:,:)%phi = dcmplx(0, 0)

                h_field_spherical(:,:)%rho = dcmplx(0, 0)
                h_field_spherical(:,:)%theta = dcmplx(0, 0)
                h_field_spherical(:,:)%phi = dcmplx(0, 0)

                e_field_spherical(2,3)%rho = dcmplx(1, 0)
                e_field_spherical(2,3)%theta = dcmplx(0, 0)
                e_field_spherical(2,3)%phi = dcmplx(0, 0)

                d%x = 0
                d%y = 0
                d%z = -3 * unit_mu

                do i = lbound(e_field_spherical, 1), ubound(e_field_spherical, 1)
                    do ii = lbound(e_field_spherical, 2), ubound(e_field_spherical, 2)
                        e_field(i, ii) = make_cartesian(e_field_spherical(i,ii), d)
                        h_field(i, ii) = make_cartesian(h_field_spherical(i,ii), d)
                    end do
                end do

                rv = lib_field_export(e_field, h_field, "temp/test_field/")


            end function
        end function
end module lib_field
