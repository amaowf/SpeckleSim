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

module file_io

    private

    ! public functions
    public :: operator( .f. )
    public :: file_exists
    public :: write_csv
    public :: read_csv
    public :: write_ppm_p2
    public :: write_ppm_p3

    public :: write_csv_cmplx_array_1d

    public :: test_file_io

    ! interface
    interface operator( .f. )
      module procedure file_exists
    end interface

    interface write_csv
        module procedure write_csv_real
        module procedure write_csv_real_array
        module procedure write_csv_cmplx
        module procedure write_csv_cmplx_array_1d
        module procedure write_csv_cmplx_array_2d
        module procedure write_csv_real_cmplx
    end interface

    interface read_csv
        module procedure file_io_csv_run_get_double
    end interface

    interface write_ppm_p3
        module procedure write_ppm_p3_single
        module procedure write_ppm_p3_RGB_single

        module procedure write_ppm_p3_double
        module procedure write_ppm_p3_RGB_double

        module procedure write_ppm_p3_image
    end interface

    interface color_map_temperature
        module procedure color_map_temperature_single
        module procedure color_map_temperature_double
    end interface

    ! parameter
    integer, parameter :: MAX_LINES = 65535

    contains

    function file_exists(filename) result(res)
      implicit none
      character(len=*),intent(in) :: filename
      logical                     :: res

      ! Check if the file exists
      inquire( file=trim(filename), exist=res )
    end function

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !
    !
    function write_ppm_p2(u, img) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double precision, allocatable, dimension(:, :), intent(inout) :: img

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ' '
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        integer :: i
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(img, 1), size(img, 2)) :: img_discretised

        img_max_value = maxval(img)
        img_min_value = minval(img)

        img_discretised = int( (img - img_min_value) * real(max_value) / (img_max_value - img_min_value), kind=2 )


        write(u, *) "P2"
        write(u, *) size(img, 1), delimiter, size(img, 2)
        write(u, *) max_value
        write(u, *) "# The part above is the header"
        write(u, *) "# P2   means this is a graymap image in ASCII"
        write(u, *) "# ", size(img, 1), delimiter, size(img, 2), "  is the width and height of the image in pixels"
        write(u, *) "# pixel value range: 0 - ", max_value
        write(u, *) "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, *) "# The part below is image data"


        do i=1, size(img, 2)
                write(u, *) img_discretised(:, i)
        end do

        rv = .true.

    end function write_ppm_p2

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   imgR: double precision, 2D array
    !       red channel
    !   imgG: double precision, 2D array
    !       green channel
    !   imgB: double precision, 2D array
    !       blue channel
    !   rescale: logical, opional (std = .false.)
    !       rescales the pixel values to 0-255
    !   hdr: logical, optional (std = .false.)
    !       use high dynamic range
    !
    function write_ppm_p3_RGB_single(u, imgR, imgG, imgB, rescale, hdr) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        real, dimension(:, :), intent(inout) :: imgR
        real, dimension(lbound(imgR,1):ubound(imgR,1), lbound(imgR,2):ubound(imgR,2)), intent(inout) :: imgG
        real, dimension(lbound(imgR,1):ubound(imgR,1), lbound(imgR,2):ubound(imgR,2)), intent(inout) :: imgB
        logical, intent(in), optional :: rescale
        logical, intent(in), optional :: hdr

        logical :: rv

        ! parameter
        character(len=*), parameter :: formatstr = '(10(I3, x, I3, x, I3, x))'
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        logical :: m_rescale
        logical :: m_hdr

        integer :: i
        integer :: ii
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(imgR, 1)*3, size(imgR, 2)) :: img_discretised_color

        integer(kind=2), dimension(3) :: buffer

        m_rescale = .false.
        if (present(rescale)) m_rescale = rescale

        m_hdr = .false.
        if (present(hdr)) m_hdr = hdr

        img_max_value = max(maxval(imgR), maxval(imgG), maxval(imgB))
        img_min_value = min(minval(imgR), minval(imgG), minval(imgB))

        write(u, '(A)') "P3"
        write(u, '(I6.1, X, I6.1)') size(imgR, 1), size(imgR, 2)
        write(u, '(I3)') max_value
        write(u, '(A)') "# The part above is the header"
        write(u, '(A)') "# P3   means this is a RGB color image in ASCII"
        write(u, '(A, I6.1, A, I6.1, A)') "# ", size(imgR, 1), " x ", size(imgR, 2), &
                                          "  is the width and height of the image in pixels"
        write(u, '(A, I3)') "# pixel value range for each color: 0 - ", max_value
        write(u, '(A, EN17.2, A, EN17.2)') "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, '(A)') "# The part below is image data: RGB triplets"

        rv = .true.

        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i=1, size(imgR, 1)
            do ii=1, size(imgR, 2)
                if (m_hdr) then
                    ! reference: Karsten Frenner: F2
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int(log(imgR(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int(log(imgG(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int(log(imgB(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                else if (m_rescale) then
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int((imgR(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int((imgG(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int((imgB(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                else
                    if (imgR(i,ii) .gt. 255) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(max_value, 2)
                    else if (imgR(i,ii) .lt. 0) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(0, 2)
                    else
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(imgR(i, ii), 2)
                    end if

                    if (imgG(i,ii) .gt. 255) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(max_value, 2)
                    else if (imgG(i,ii) .lt. 0) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(0, 2)
                    else
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(imgG(i, ii), 2)
                    end if

                    if (imgB(i,ii) .gt. 255) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(max_value, 2)
                    else if (imgB(i,ii) .lt. 0) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(0, 2)
                    else
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(imgB(i, ii), 2)
                    end if
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        if (rv) then
            do i=1, size(img_discretised_color, 2)
                write(u, formatstr) img_discretised_color(:, i)
            end do
        end if

    end function write_ppm_p3_RGB_single

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   imgR: double precision, 2D array
    !       red channel
    !   imgG: double precision, 2D array
    !       green channel
    !   imgB: double precision, 2D array
    !       blue channel
    !   rescale: logical, opional (std = .false.)
    !       rescales the pixel values to 0-255
    !   hdr: logical, optional (std = .false.)
    !       use high dynamic range
    !
    function write_ppm_p3_RGB_double(u, imgR, imgG, imgB, rescale, hdr) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double precision, dimension(:, :), intent(inout) :: imgR
        double precision, dimension(lbound(imgR,1):ubound(imgR,1), lbound(imgR,2):ubound(imgR,2)), intent(inout) :: imgG
        double precision, dimension(lbound(imgR,1):ubound(imgR,1), lbound(imgR,2):ubound(imgR,2)), intent(inout) :: imgB
        logical, intent(in), optional :: rescale
        logical, intent(in), optional :: hdr

        logical :: rv

        ! parameter
        character(len=*), parameter :: formatstr = '(10(I3, x, I3, x, I3, x))'
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        logical :: m_rescale
        logical :: m_hdr

        integer :: i
        integer :: ii
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(imgR, 1)*3, size(imgR, 2)) :: img_discretised_color

        integer(kind=2), dimension(3) :: buffer

        m_rescale = .true.
        if (present(rescale)) m_rescale = rescale

        m_hdr = .false.
        if (present(hdr)) m_hdr = hdr

        img_max_value = max(maxval(imgR), maxval(imgG), maxval(imgB))
        img_min_value = min(minval(imgR), minval(imgG), minval(imgB))

        write(u, '(A)') "P3"
        write(u, '(I6.1, X, I6.1)') size(imgR, 1), size(imgR, 2)
        write(u, '(I3)') max_value
        write(u, '(A)') "# The part above is the header"
        write(u, '(A)') "# P3   means this is a RGB color image in ASCII"
        write(u, '(A, I6.1, A, I6.1, A)') "# ", size(imgR, 1), " x ", size(imgR, 2), &
                                          "  is the width and height of the image in pixels"
        write(u, '(A, I3)') "# pixel value range for each color: 0 - ", max_value
        write(u, '(A, EN17.2, A, EN17.2)') "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, '(A)') "# The part below is image data: RGB triplets"

        rv = .true.

        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i=1, size(imgR, 1)
            do ii=1, size(imgR, 2)
                if (m_hdr) then
                    ! reference: Karsten Frenner: F2
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int(log(imgR(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int(log(imgG(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int(log(imgB(i, ii) + 1) / log(img_max_value) * max_value, kind=2)
                else if (m_rescale) then
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int((imgR(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int((imgG(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int((imgB(i, ii)-img_min_value) / img_max_value * max_value, kind=2)
                else
                    if (imgR(i,ii) .gt. 255) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(max_value, 2)
                    else if (imgR(i,ii) .lt. 0) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(0, 2)
                    else
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(imgR(i, ii), 2)
                    end if

                    if (imgG(i,ii) .gt. 255) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(max_value, 2)
                    else if (imgG(i,ii) .lt. 0) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(0, 2)
                    else
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(imgG(i, ii), 2)
                    end if

                    if (imgB(i,ii) .gt. 255) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(max_value, 2)
                    else if (imgB(i,ii) .lt. 0) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(0, 2)
                    else
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(imgB(i, ii), 2)
                    end if
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        if (rv) then
            do i=1, size(img_discretised_color, 2)
                write(u, formatstr) img_discretised_color(:, i)
            end do
        end if

    end function write_ppm_p3_RGB_double

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   imgR: double precision, 2D array
    !       red channel
    !   imgG: double precision, 2D array
    !       green channel
    !   imgB: double precision, 2D array
    !       blue channel
    !   rescale: logical, opional (std = .false.)
    !       rescales the pixel values to 0-255
    !   hdr: logical, optional (std = .false.)
    !       use high dynamic range
    !
    function write_ppm_p3_image(u, image, rescale, hdr) result (rv)
        use lib_image_public
        implicit none
        ! dummy
        integer, intent(in) :: u
        type(image_type), intent(in) :: image
        logical, intent(in), optional :: rescale
        logical, intent(in), optional :: hdr

        logical :: rv

        ! parameter
        character(len=*), parameter :: formatstr = '(10(I3, x, I3, x, I3, x))'
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        logical :: m_rescale
        logical :: m_hdr

        integer :: i
        integer :: ii
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(image%pixel, 1)*3, size(image%pixel, 2)) :: img_discretised_color

        integer(kind=2), dimension(3) :: buffer

        m_rescale = .true.
        if (present(rescale)) m_rescale = rescale

        m_hdr = .false.
        if (present(hdr)) m_hdr = hdr

        img_max_value = max(maxval(image%pixel(:,:)%RGB(1)), maxval(image%pixel(:,:)%RGB(2)), maxval(image%pixel(:,:)%RGB(3)))
        img_min_value = min(minval(image%pixel(:,:)%RGB(1)), minval(image%pixel(:,:)%RGB(2)), minval(image%pixel(:,:)%RGB(3)))

        write(u, '(A)') "P3"
        write(u, '(I6.1, X, I6.1)') size(image%pixel, 1), size(image%pixel, 2)
        write(u, '(I3)') max_value
        write(u, '(A)') "# The part above is the header"
        write(u, '(A)') "# P3   means this is a RGB color image in ASCII"
        write(u, '(A, I6.1, A, I6.1, A)') "# ", size(image%pixel, 1), " x ", size(image%pixel, 2), &
                                          "  is the width and height of the image in pixels"
        write(u, '(A, I3)') "# pixel value range for each color: 0 - ", max_value
        write(u, '(A, EN17.2, A, EN17.2)') "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, '(A)') "# The part below is image data: RGB triplets"

        rv = .true.

        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i=1, size(image%pixel, 1)
            do ii=1, size(image%pixel, 2)
                if (m_hdr) then
                    ! reference: Karsten Frenner: F2
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int(log(image%pixel(i, ii)%RGB(1) + 1) &
                                                               / log(img_max_value) * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int(log(image%pixel(i, ii)%RGB(2) + 1) &
                                                               / log(img_max_value) * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int(log(image%pixel(i, ii)%RGB(3) + 1) &
                                                               / log(img_max_value) * max_value, kind=2)
                else if (m_rescale) then
                    ! red
                    img_discretised_color(3*(i-1)+1, ii) = int((image%pixel(i, ii)%RGB(1)-img_min_value) &
                                                                / img_max_value * max_value, kind=2)
                    ! green
                    img_discretised_color(3*(i-1)+2, ii) = int((image%pixel(i, ii)%RGB(2)-img_min_value) &
                                                                / img_max_value * max_value, kind=2)
                    ! blue
                    img_discretised_color(3*(i-1)+3, ii) = int((image%pixel(i, ii)%RGB(3)-img_min_value) &
                                                                / img_max_value * max_value, kind=2)
                else
                    if (image%pixel(i, ii)%RGB(1) .gt. 255) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(max_value, 2)
                    else if (image%pixel(i, ii)%RGB(1) .lt. 0) then
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(0, 2)
                    else
                        ! red
                        img_discretised_color(3*(i-1)+1, ii) = int(image%pixel(i, ii)%RGB(1), 2)
                    end if

                    if (image%pixel(i, ii)%RGB(2) .gt. 255) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(max_value, 2)
                    else if (image%pixel(i, ii)%RGB(2) .lt. 0) then
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(0, 2)
                    else
                        ! green
                        img_discretised_color(3*(i-1)+2, ii) = int(image%pixel(i, ii)%RGB(2), 2)
                    end if

                    if (image%pixel(i, ii)%RGB(3) .gt. 255) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(max_value, 2)
                    else if (image%pixel(i, ii)%RGB(3) .lt. 0) then
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(0, 2)
                    else
                        ! blue
                        img_discretised_color(3*(i-1)+3, ii) = int(image%pixel(i, ii)%RGB(3), 2)
                    end if
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        if (rv) then
            do i=1, size(img_discretised_color, 2)
                write(u, formatstr) img_discretised_color(:, i)
            end do
        end if

    end function write_ppm_p3_image

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   img: double precision, 2D array
    !   color_map: integer, opional (std = 1)
    !       - 1: temperatur
    !
    !
    function write_ppm_p3_single(u, img, color_map, logarithmic) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        real, dimension(:, :), intent(inout) :: img
        integer, intent(in), optional :: color_map
        logical, intent(in), optional :: logarithmic

        logical :: rv

        ! parameter
        character(len=*), parameter :: formatstr = '(10(I3, x, I3, x, I3, x))'
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        integer :: i
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(img, 1)*3, size(img, 2)) :: img_discretised_color

        img_max_value = maxval(img)
        img_min_value = minval(img)

        write(u, '(A)') "P3"
        write(u, '(I6.1, X, I6.1)') size(img, 1), size(img, 2)
        write(u, '(I3)') max_value
        write(u, '(A)') "# The part above is the header"
        write(u, '(A)') "# P3   means this is a RGB color image in ASCII"
        write(u, '(A, I6.1, A, I6.1, A)') "# ", size(img, 1), " x ", size(img, 2), &
                                          "  is the width and height of the image in pixels"
        write(u, '(A, I3)') "# pixel value range for each color: 0 - ", max_value
        write(u, '(A, EN17.2, A, EN17.2)') "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, '(A)') "# The part below is image data: RGB triplets"

        rv = .true.

        if (present(color_map)) then
            if (color_map .eq. 1) then
                 call color_map_temperature(img, max_value, img_discretised_color, logarithmic=logarithmic)
            else
                rv = .false.
            end if
        else
            call color_map_temperature(img, max_value, img_discretised_color, logarithmic=logarithmic)
        end if

        if (rv) then
            do i=1, size(img_discretised_color, 2)
                write(u, formatstr) img_discretised_color(:, i)
            end do
        end if

    end function write_ppm_p3_single

    ! writes two columns into a ppm file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   img: double precision, 2D array
    !   color_map: integer, opional (std = 1)
    !       - 1: temperatur
    !
    !
    function write_ppm_p3_double(u, img, color_map, logarithmic) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double precision, dimension(:, :), intent(inout) :: img
        integer, intent(in), optional :: color_map
        logical, intent(in), optional :: logarithmic

        logical :: rv

        ! parameter
        character(len=*), parameter :: formatstr = '(10(I3, x, I3, x, I3, x))'
        integer(kind=2), parameter :: max_value = 255

        ! auxiliary
        integer :: i
        double precision :: img_max_value
        double precision :: img_min_value

        integer(kind=2), dimension(size(img, 1)*3, size(img, 2)) :: img_discretised_color

        img_max_value = maxval(img)
        img_min_value = minval(img)

        write(u, '(A)') "P3"
        write(u, '(I6.1, X, I6.1)') size(img, 1), size(img, 2)
        write(u, '(I3)') max_value
        write(u, '(A)') "# The part above is the header"
        write(u, '(A)') "# P3   means this is a RGB color image in ASCII"
        write(u, '(A, I6.1, A, I6.1, A)') "# ", size(img, 1), " x ", size(img, 2), &
                                          "  is the width and height of the image in pixels"
        write(u, '(A, I3)') "# pixel value range for each color: 0 - ", max_value
        write(u, '(A, EN17.2, A, EN17.2)') "# data point value range: ", img_min_value, " - ", img_max_value
        write(u, '(A)') "# The part below is image data: RGB triplets"

        rv = .true.

        if (present(color_map)) then
            if (color_map .eq. 1) then
                 call color_map_temperature(img, max_value, img_discretised_color, logarithmic=logarithmic)
            else
                rv = .false.
            end if
        else
            call color_map_temperature(img, max_value, img_discretised_color, logarithmic=logarithmic)
        end if

        if (rv) then
            do i=1, size(img_discretised_color, 2)
                write(u, formatstr) img_discretised_color(:, i)
            end do
        end if

    end function write_ppm_p3_double

    subroutine color_map_temperature_single(img, max_value, rv, logarithmic)
        implicit none
        ! dummy
        real, dimension(:, :), intent(in) :: img
        integer(kind=2), intent(in) :: max_value
        logical, intent(in), optional :: logarithmic

        integer(kind=2), dimension(size(img, 1)*3, size(img, 2)), intent(inout) :: rv

        ! auxiliary
        integer :: i
        integer :: ii
        double precision :: img_max_abs_value
        integer(kind=2) :: buffer

        img_max_abs_value = max(abs(maxval(img)), abs(minval(img)))

        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i=1, size(img, 1)
            do ii=1, size(img, 2)
                if (img(i, ii) .gt. 0) then
                    if (present(logarithmic)) then
                        if (logarithmic) then
                            buffer = int( log(img(i,ii)+1) * real(max_value) / log(img_max_abs_value+1), kind=2)
                        else
                            buffer = int( img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                        end if
                    else
                        buffer = int( img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                    end if

                    ! red
                    rv(3*(i-1)+1, ii) = max_value
                    ! green
                    rv(3*(i-1)+2, ii) = max_value - buffer
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value - buffer
                else if (img(i, ii) .eq. 0) then
                    ! red
                    rv(3*(i-1)+1, ii) = max_value
                    ! green
                    rv(3*(i-1)+2, ii) = max_value
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value
                else
                    if (present(logarithmic)) then
                        if (logarithmic) then
                            buffer = int( log(-img(i,ii)+1) * real(max_value) / log(img_max_abs_value+1), kind=2)
                        else
                            buffer = int( -img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                        end if
                    else
                        buffer = int( -img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                    end if

!                    buffer = int( img(i,ii) * real(max_value) / img_min_value, kind=2)
                    ! red
                    rv(3*(i-1)+1, ii) = max_value - buffer
                    ! green
                    rv(3*(i-1)+2, ii) = max_value - buffer
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value
                end if
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine

    subroutine color_map_temperature_double(img, max_value, rv, logarithmic)
        implicit none
        ! dummy
        double precision, dimension(:, :), intent(in) :: img
        integer(kind=2), intent(in) :: max_value
        logical, intent(in), optional :: logarithmic

        integer(kind=2), dimension(size(img, 1)*3, size(img, 2)), intent(inout) :: rv

        ! auxiliary
        integer :: i
        integer :: ii
        double precision :: img_max_abs_value
        integer(kind=2) :: buffer

        img_max_abs_value = max(abs(maxval(img)), abs(minval(img)))

        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i=1, size(img, 1)
            do ii=1, size(img, 2)
                if (img(i, ii) .gt. 0) then
                    if (present(logarithmic)) then
                        if (logarithmic) then
                            buffer = int( log(img(i,ii)+1) * real(max_value) / log(img_max_abs_value+1), kind=2)
                        else
                            buffer = int( img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                        end if
                    else
                        buffer = int( img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                    end if

                    ! red
                    rv(3*(i-1)+1, ii) = max_value
                    ! green
                    rv(3*(i-1)+2, ii) = max_value - buffer
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value - buffer
                else if (img(i, ii) .eq. 0) then
                    ! red
                    rv(3*(i-1)+1, ii) = max_value
                    ! green
                    rv(3*(i-1)+2, ii) = max_value
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value
                else
                    if (present(logarithmic)) then
                        if (logarithmic) then
                            buffer = int( log(-img(i,ii)+1) * real(max_value) / log(img_max_abs_value+1), kind=2)
                        else
                            buffer = int( -img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                        end if
                    else
                        buffer = int( -img(i,ii) * real(max_value) / img_max_abs_value, kind=2)
                    end if

!                    buffer = int( img(i,ii) * real(max_value) / img_min_value, kind=2)
                    ! red
                    rv(3*(i-1)+1, ii) = max_value - buffer
                    ! green
                    rv(3*(i-1)+2, ii) = max_value - buffer
                    ! blue
                    rv(3*(i-1)+3, ii) = max_value
                end if
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_real(u, header, c1, c2, c3, c4) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        double precision, dimension(:) :: c1
        double precision, dimension(size(c1)), optional :: c2
        double precision, dimension(size(c1)), optional :: c3
        double precision, dimension(size(c1)), optional :: c4

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i))
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (present(c2) .and. .not. present(c3) .and. .not. present(c4)) then
!                write(u, '(F2.15,A,ES2.15)') c1(i), delimiter, c2(i)
                write(u, *) c1(i), delimiter, c2(i)
            else if (present(c2) .and. present(c3) .and. .not. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i)
                write(u, *) c1(i), delimiter, c2(i), delimiter, c3(i)
            else if (present(c2) .and. present(c3) .and. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i), delimiter, c4(i)
                write(u, *) c1(i), delimiter, c2(i), delimiter, c3(i), delimiter, c4(i)
            else
                write(u, '(ES2.15)') c1(i)
            end if
        end do

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_real_array(u, c, header) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double precision, dimension(:,:) :: c
        character(len=*), dimension(:), optional :: header

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        integer :: ii
        character(len=511) :: dummy_str

        if ( present(header) ) then
            dummy_str = trim(header(1))
            do i=2, size(header)
                dummy_str = trim(dummy_str) // delimiter // trim(header(i))
            end do
            write(u, *) trim(dummy_str)
        end if

        do i = lbound(c, 1), ubound(c, 1)
            ii = lbound(c, 2)
            write(u, '(ES25.16)', advance='no') c(i, ii)
            do ii = lbound(c, 2)+1, ubound(c, 2)
                write(u, '(A, ES25.16)', advance='no') delimiter, c(i, ii)
            end do
            write(u, *) '' ! line break
        end do

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_cmplx(u, header, c1, c2, c3) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        complex(kind=8), dimension(:) :: c1
        complex(kind=8), dimension(size(c1)), optional :: c2
        complex(kind=8), dimension(size(c1)), optional :: c3

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i))
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (present(c2) .and. .not. present(c3)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c2) .and. present(c3)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i))
            else
!                write(u, '(ES2.15,A,ES2.15)') real(c1(i)), delimiter, aimag(c2(i))
                write(u, *) real(c1(i)), delimiter, aimag(c2(i))
            end if
        end do

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_cmplx_array_1d(u, c, header) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double complex, dimension(:), intent(in) :: c
        character(len=*), dimension(:), intent(in), optional :: header

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str
        character(len=25), dimension(2) :: m_header

        m_header(1) = "index"
        m_header(2) = "value"

        if ( present(header) ) then
            dummy_str = trim(header(1))
            do i=2, size(header)
                dummy_str = trim(dummy_str) // delimiter // trim(header(i)) // "(real)" // delimiter // trim(header(i)) // "(cmplx)"
            end do
            write(u, *) trim(dummy_str)

            do i = lbound(c, 1), ubound(c, 1)
                write(u, '(I6, A, ES25.16, A, ES25.16)') i, delimiter, real(c(i)), delimiter, aimag(c(i))
            end do
        else
            dummy_str = trim(m_header(1))
            do i=2, size(m_header)
                dummy_str = trim(dummy_str) // delimiter // trim(m_header(i)) // "(real)" // delimiter &
                                                         // trim(m_header(i)) // "(cmplx)"
            end do
            write(u, *) trim(dummy_str)

            do i = lbound(c, 1), ubound(c, 1)
                write(u, '(I6, A, ES25.16, A, ES25.16)') i, delimiter, real(c(i)), delimiter, aimag(c(i))
            end do
        end if

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_cmplx_array_2d(u, c, header) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        double complex, dimension(:,:) :: c
        character(len=*), dimension(:), optional :: header

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        integer :: ii
        character(len=10000) :: dummy_str
        character(len=25), dimension(2) :: m_header
        character(len=6) :: int_str



        if ( present(header) ) then
            dummy_str = trim(header(1))
            do i=2, size(header)
                dummy_str = trim(dummy_str) // delimiter // trim(header(i)) // "(real)" // delimiter // trim(header(i)) // "(cmplx)"
            end do
            write(u, *) trim(dummy_str)

            do i = lbound(c, 1), ubound(c, 1)
                ii = lbound(c, 2)
                write(u, '(ES25.16, A, ES25.16)', advance='no') real(c(i, ii)), delimiter, aimag(c(i, ii))
                do ii = lbound(c, 2)+1, ubound(c, 2)
                    write(u, '(A, ES25.16, A, ES25.16)', advance='no') delimiter, real(c(i, ii)), delimiter, aimag(c(i, ii))
                end do
                write(u, *) '' ! line break
            end do
        else
            m_header(1) = "row \ column"
            m_header(2) = " "
            dummy_str = trim(m_header(1))
            do i=lbound(c, 2), ubound(c, 2)
                write(int_str, '(I6)') i
                dummy_str = trim(dummy_str) // delimiter // int_str // trim(m_header(2)) // "(real)" // delimiter &
                                                         // int_str // trim(m_header(2)) // "(cmplx)"
            end do
            write(u, *) trim(dummy_str)

            do i = lbound(c, 1), ubound(c, 1)
                ii = lbound(c, 2)
                write(u, '(I6, A, ES25.16, A, ES25.16)', advance='no') i, delimiter, real(c(i, ii)), delimiter, aimag(c(i, ii))
                do ii = lbound(c, 2)+1, ubound(c, 2)
                    write(u, '(A, ES25.16, A, ES25.16)', advance='no') delimiter, real(c(i, ii)), delimiter, aimag(c(i, ii))
                end do
                write(u, *) '' ! line break
            end do
        end if

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_real_cmplx(u, header, c1, c2, c3, c4) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        real(kind=8), dimension(:) :: c1
        complex(kind=8), dimension(size(c1)) :: c2
        complex(kind=8), dimension(size(c1)), optional :: c3
        complex(kind=8), dimension(size(c1)), optional :: c4

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i)) // "(real)" // delimiter // trim(header(i)) // "(cmplx)"
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (.not. present(c3) .and. .not. present(c4)) then
                write(u, *) c1(i), delimiter, real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c3) .and. .not. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) c1(i), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i))
            else if (present(c3) .and. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) c1(i), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i)), delimiter, &
                            real(c4(i)), delimiter, aimag(c4(i))
            end if
        end do

        rv = .true.

    end function

    ! reads a csv file with double entries
    !
    ! Arguments
    ! ----
    !   file_name: character(len=*)
    !       name of the csv file
    !   header_lines: integer, optional (std: 1)
    !       number of header lines
    subroutine file_io_csv_run_get_double(file_name, columns, data, header_lines, data_lines)
        implicit none
        ! dummy
        character(len=*) :: file_name
        integer :: columns
        integer, optional :: header_lines
        integer, optional :: data_lines

        double precision, dimension(:,:), allocatable, intent(inout) :: data

        ! auxiliary
        integer :: m_nlines
        integer :: m_header_lines
        integer :: m_data_lines
        integer, dimension(2) :: m_line_range

        ! --- init ---
        m_header_lines = 1
        if (present(header_lines)) then
            m_header_lines = header_lines
        end if

        m_data_lines = -1
        if (present(data_lines)) then
            m_data_lines = data_lines
        end if


        ! get number of lines
        m_nlines = file_io_csv_get_number_of_lines(file_name)

        if (m_nlines .lt. m_header_lines) then
            print *, "file_io_csv_run_get_double: ERROR"
            print *, "  no. of lines < header_lines"
            print *, "  no. of lines: ", m_nlines
            print *, "  header lines: ", m_header_lines
            return
        end if

        if (m_data_lines .gt. 0 .and.&
            m_data_lines .le. m_nlines-m_header_lines) then
            m_line_range = (/ m_header_lines + 1, m_header_lines + m_data_lines /)
        else
            m_line_range = (/ m_header_lines + 1, m_nlines /)
        end if

        ! read data
        m_nlines = file_io_csv_read_double(file_name, m_line_range, columns, data)

    end subroutine file_io_csv_run_get_double

    ! Counts the number of lines in a text file
    !
    ! Argument
    ! ----
    !   file_name: character(len=*)
    !       file name
    !
    ! Result
    ! ----
    !   lines: integer
    !       number of lines of the file
    !       specific values:
    !           -1: open statement is NOT successful
    !           -2: the maximum number of rows has been reached (nax = 65535)
    !
    function file_io_csv_get_number_of_lines(file_name) result(lines)
        implicit none
        ! dummy
        character(len=*) :: file_name

        integer :: lines

        ! auxiliary
        integer :: io
        integer :: unit

        lines = 0

        open (newunit = unit, file = file_name)

        do
            if (lines .gt. MAX_LINES) then
                lines = -2
                exit
            end if

            read(unit,*,iostat=io)


            if (io > 0) then
                print *, "file_io_csv_get_number_of_lines:"
                print *, "  Check input.  Something was wrong"
                print *, "  line: ", lines+1
                lines = -1
                exit
            else if (io < 0) then
                ! end of file
                exit
            else
                lines = lines + 1
            end if
        end do

        close(unit)

    end function

    ! Counts the number of lines in a text file
    !
    ! Argument
    ! ----
    !   unit: integer
    !       unit number, range 9-99
    !
    ! Result
    ! ----
    !   lines: integer
    !       number of lines of the file
    !       specific values:
    !           -1: open statement is NOT successful
    !           -2: the maximum number of rows has been reached (nax = 65535)
    !
    function file_io_csv_read_double(file_name, line_range, number_of_columns, data) result(lines)
        implicit none
        ! dummy
        character(len=*) :: file_name
        integer, dimension(2) :: line_range
        integer :: number_of_columns

        double precision, dimension(:,:), allocatable, intent(inout) :: data
        integer :: lines

        ! auxiliary
        integer :: unit
        integer :: io

        if (allocated(data)) then
            deallocate(data)
        end if
        allocate(data(line_range(1):line_range(2), number_of_columns))

        lines = 0

        open (newunit = unit, file = file_name)

        do
            if (lines .gt. MAX_LINES) then
                lines = -2
                exit
            end if

            if (lines .lt. line_range(1)-1) then
                ! last line < line_range(1) - 1
                read(unit,*,iostat=io)
            else if (lines .ge. line_range(2)) then
                ! last line >= line_range(2)
                exit
            else
                read(unit,*,iostat=io) data(lines+1, :)
            end if


            if (io > 0) then
                print *, "file_io_csv_get_number_of_lines:"
                print *, "  Check input.  Something was wrong"
                print *, "  line: ", lines+1
                lines = -1
                exit
            else if (io < 0) then
                ! end of file
                exit
            else
                lines = lines + 1
            end if
        end do

        close(unit)

    end function file_io_csv_read_double

    subroutine test_file_io()
        implicit none

        call test_write_ppm_p2
        call test_write_ppm_p3
        call test_write_ppm_p3_log

        contains

            subroutine test_write_ppm_p2()
                implicit none

                ! auxiliary
                integer :: u
                double precision, allocatable, dimension(:, :) :: img

                logical :: rv

                allocate (img(10,3))

                img = reshape((/ 0,0,0,1,2,3,4,5,6,7, &
                                 0,0,0,0,0,0,0,0,0,0, &
                                 0,0,0,0,0,0,0,0,0,0 /), shape(img))

                u = 99
                open(unit=u, file="test.ppm", status='unknown')
                rv = write_ppm_p2(u, img)
                close(u)

            end subroutine test_write_ppm_p2

            subroutine test_write_ppm_p3()
                implicit none

                ! auxiliary
                integer :: u
                double precision, allocatable, dimension(:, :) :: img

                logical :: rv

                allocate (img(10,3))

                img = reshape((/ 0,0,0,1,2,3,4,5,6,7, &
                                 0,0,0,0,0,0,0,0,0,0, &
                                 -7,-6,-5,-4,-3,-2,-1,0,0,0 /), shape(img))

                u = 99
                open(unit=u, file="test.ppm", status='unknown')
                rv = write_ppm_p3(u, img)
                close(u)

            end subroutine test_write_ppm_p3

            subroutine test_write_ppm_p3_log()
                implicit none

                ! auxiliary
                integer :: u
                double precision, allocatable, dimension(:, :) :: img

                logical :: rv

                allocate (img(10,3))

!                img = reshape((/ 0,0,0,1,2,3,4,5,6,7, &
!                                 0,0,0,0,0,0,0,0,0,0, &
!                                 -7,-6,-5,-4,-3,-2,-1,0,0,0 /), shape(img))

                img(:,1) = (/ 0,0,0,1,2,3,4,5,6,7 /)
                img(:,2) = (/ 0,0,0,0,0,0,0,0,0,0 /)
                img(:,3) = (/ -7,-6,-5,-4,-3,-2,-1,0,0,0 /)

                u = 99
                open(unit=u, file="test.ppm", status='unknown')
                rv = write_ppm_p3(u, img, logarithmic=.true.)
                close(u)

            end subroutine test_write_ppm_p3_log

    end subroutine

end module
