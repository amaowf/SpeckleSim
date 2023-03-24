module lib_image
    use lib_colour
    use lib_image_type
    implicit none

    private

    public :: sum
    public :: lib_image_generate_mask_binary_random
    public :: lib_image_generate_mask_binary_random_pattern
    public :: lib_image_set_pixel

    interface sum
        module procedure lib_image_get_sum
    end interface

contains

    ! Argument
    ! ----
    !   rows: integer
    !       rows of the mask image
    !   columns: integer
    !       columns of the mask image
    !   threshold: real, optional (std = 0.5)
    !       threshold level for selecting colour 1 or 2
    !   colour_1: colour_type
    !       colour of the pixel if the random number is lower than *thrashold*
    !   colour_2: colour_type
    !       colour of the pixel if the random number is greater than or equal to *thrashold*
    !
    ! Returns
    ! ----
    !   mask: image_type
    !       binary random mask
    function lib_image_generate_mask_binary_random(rows, columns, threshold, colour_1, colour_2) result(mask)
        !$  use omp_lib
        implicit none
        ! dummy
        integer, intent(in) :: rows
        integer, intent(in) :: columns
        real, intent(in), optional :: threshold
        type(colour_type), intent(in), optional :: colour_1
        type(colour_type), intent(in), optional :: colour_2

        type(image_type) :: mask

        ! auxiliary
        integer :: m_r
        integer :: m_c
        real :: m_ran
        real :: m_threshold
        type(colour_type) :: m_colour_1
        type(colour_type) :: m_colour_2

        if (present(threshold)) then
            m_threshold = threshold
        else
            m_threshold = 0.5
        end if

        if (present(colour_1)) then
            m_colour_1 = colour_1
        else
            m_colour_1 = make_colour_RGB(1., 1., 1.)
        end if

        if (present(colour_2)) then
            m_colour_2 = colour_2
        else
            m_colour_2 = make_colour_RGB(0., 0., 0.)
        end if

        allocate(mask%pixel(rows, columns))

        call random_seed()

        !$OMP PARALLEL DO PRIVATE(m_r, m_c, m_ran)
        do m_r = 1, rows
            do m_c = 1, columns
                call random_number(m_ran)
                if (m_ran .ge. m_threshold) then
                    mask%pixel(m_r, m_c) = m_colour_1
                else
                    mask%pixel(m_r, m_c) = m_colour_2
                end if
            end do
        end do
        !$OMP END PARALLEL DO

    end function lib_image_generate_mask_binary_random

    ! Argument
    ! ----
    !   pattern: image_type
    !       pixels of the same color are treated equally in the randomly generated mask
    !   threshold: real, optional (std = 0.5)
    !       threshold level for selecting colour 1 or 2
    !   colour_1: colour_type, optional (std: RGB = (1, 1, 1))
    !       colour of the pixel if the random number is lower than *thrashold*
    !   colour_2: colour_type, optional (std: RGB = (0, 0, 0))
    !       colour of the pixel if the random number is greater than or equal to *thrashold*
    !
    ! Returns
    ! ----
    !   mask: image_type
    !       binary random mask
    function lib_image_generate_mask_binary_random_pattern(pattern, threshold, colour_1, colour_2) result(mask)
        !$  use omp_lib
        implicit none
        ! dummy
        type(image_type), intent(in) :: pattern
        real, intent(in), optional :: threshold
        type(colour_type), intent(in), optional :: colour_1
        type(colour_type), intent(in), optional :: colour_2

        type(image_type) :: mask

        ! auxiliary
        integer :: m_r
        integer :: m_c
        real, dimension(:), allocatable :: ran
        real :: m_threshold
        type(colour_type) :: m_colour_1
        type(colour_type) :: m_colour_2
        double precision, dimension(2) :: m_r_range
        double precision, dimension(2) :: m_g_range
        double precision, dimension(2) :: m_b_range
        integer, dimension(2) :: m_range
        integer :: buffer

        if (present(threshold)) then
            m_threshold = threshold
        else
            m_threshold = 0.5
        end if

        if (present(colour_1)) then
            m_colour_1 = colour_1
        else
            m_colour_1 = make_colour_RGB(1., 1., 1.)
        end if

        if (present(colour_2)) then
            m_colour_2 = colour_2
        else
            m_colour_2 = make_colour_RGB(0., 0., 0.)
        end if

        allocate(mask%pixel(lbound(pattern%pixel, 1):ubound(pattern%pixel, 1), &
                            lbound(pattern%pixel, 2):ubound(pattern%pixel, 2)))

        call random_seed()

        m_r_range(1) = minval(pattern%pixel(:,:)%RGB(1))
        m_r_range(2) = maxval(pattern%pixel(:,:)%RGB(1))

        m_g_range(1) = minval(pattern%pixel(:,:)%RGB(2))
        m_g_range(2) = maxval(pattern%pixel(:,:)%RGB(2))

        m_b_range(1) = minval(pattern%pixel(:,:)%RGB(3))
        m_b_range(2) = maxval(pattern%pixel(:,:)%RGB(3))

        m_range(1) = int(min(m_r_range(1), m_g_range(1), m_b_range(1)))
        m_range(2) = int(max(m_r_range(2), m_g_range(2), m_b_range(2)))

        allocate(ran(m_range(1):m_range(2)))
        call random_number(ran)

        !$OMP PARALLEL DO PRIVATE(m_r, m_c, buffer)
        do m_r = lbound(pattern%pixel, 1), ubound(pattern%pixel, 1)
            do m_c = lbound(pattern%pixel, 2), ubound(pattern%pixel, 2)
                buffer = int(pattern%pixel(m_r, m_c)%RGB(1))
                if (ran(buffer) .lt. m_threshold) then
                    mask%pixel(m_r, m_c)%RGB(1) = m_colour_1%RGB(1)
                else
                    mask%pixel(m_r, m_c)%RGB(1) = m_colour_2%RGB(1)
                end if

                buffer = int(pattern%pixel(m_r, m_c)%RGB(2))
                if (ran(buffer) .lt. m_threshold) then
                    mask%pixel(m_r, m_c)%RGB(2) = m_colour_1%RGB(2)
                else
                    mask%pixel(m_r, m_c)%RGB(2) = m_colour_2%RGB(2)
                end if

                buffer = int(pattern%pixel(m_r, m_c)%RGB(3))
                if (ran(buffer) .lt. m_threshold) then
                    mask%pixel(m_r, m_c)%RGB(3) = m_colour_1%RGB(3)
                else
                    mask%pixel(m_r, m_c)%RGB(3) = m_colour_2%RGB(3)
                end if
            end do
        end do
        !$OMP END PARALLEL DO

    end function

    ! Argument
    ! ----
    !   image: image_type
    !       image
    !   channel: integer
    !       RGB (1, 2, 3) colour channel of the image
    !
    ! Returns
    ! ----
    !   res: double precision
    !       sum of all pixel values of the selectet colour channel
    !
    function lib_image_get_sum(image, channel) result (res)
        implicit none
        ! dummy
        type(image_type), intent(in) :: image
        integer, intent(in) :: channel

        double precision :: res

        if (channel .ge. 1 .or. channel .le. 3) then
            res = sum(image%pixel(:,:)%RGB(channel))
        else
            print *, "lib_image_get_sum: ERROR"
            print *, "  channel not defined: ", channel
        end if

    end function lib_image_get_sum

    subroutine lib_image_set_pixel(image, colour, range_r, range_c)
        implicit none
        ! dummy
        type(image_type), intent(inout) :: image
        type(colour_type), intent(in) :: colour
        integer, dimension(2), intent(in), optional :: range_r
        integer, dimension(2), intent(in), optional :: range_c

        ! auxiliary
        integer :: m_r
        integer :: m_c
        integer, dimension(2) :: m_range_r
        integer, dimension(2) :: m_range_c

        if (present(range_r)) then
            m_range_r = range_r
        else
            m_range_r(1) = lbound(image%pixel, 1)
            m_range_r(2) = ubound(image%pixel, 1)
        end if

        if (present(range_c)) then
            m_range_c = range_c
        else
            m_range_c(1) = lbound(image%pixel, 2)
            m_range_c(2) = ubound(image%pixel, 2)
        end if

        !$OMP PARALLEL DO PRIVATE(m_r, m_c)
        do m_r = m_range_r(1), m_range_r(2)
            do m_c = m_range_c(1), m_range_c(2)
                image%pixel(m_r, m_c) = colour
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine lib_image_set_pixel

end module lib_image
