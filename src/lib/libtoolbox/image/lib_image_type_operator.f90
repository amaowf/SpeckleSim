module lib_image_type_operator
    use lib_colour
    use lib_image_type
    implicit none

    private

    public :: operator(*)
    public :: make_image_RGB

    interface operator(*)
        module procedure lib_image_type_operator_img_mul_img
    end interface

    interface make_image_RGB
        module procedure lib_image_make_image_RGB_real
        module procedure lib_image_make_image_RGB_dble
    end interface

contains

    ! Argument
    ! ----
    !   lhs: image_type
    !   rhs: image_type
    !       same size as lhs
    !
    ! Returns
    ! ----
    !   res: image_type
    !
    function lib_image_type_operator_img_mul_img(lhs, rhs) result(res)
        implicit none
        ! dummy
        type(image_type), intent(in) :: lhs
        type(image_type), intent(in) :: rhs

        type(image_type) :: res

        ! auxiliary
        integer :: m_r
        integer :: m_c

        if (lbound(lhs%pixel, 1) .eq. lbound(rhs%pixel, 1) &
            .and. lbound(lhs%pixel, 1) .eq. lbound(rhs%pixel, 1) &
            .and. lbound(lhs%pixel, 2) .eq. lbound(rhs%pixel, 2) &
            .and. lbound(lhs%pixel, 2) .eq. lbound(rhs%pixel, 2) &
            ) then

            allocate(res%pixel(lbound(lhs%pixel, 1):ubound(lhs%pixel, 1), &
                               lbound(lhs%pixel, 2):ubound(lhs%pixel, 2)))

            !$OMP PARALLEL DO PRIVATE(m_r, m_c)
            do m_r = lbound(lhs%pixel, 1), ubound(lhs%pixel, 1)
                do m_c = lbound(lhs%pixel, 2), ubound(lhs%pixel, 2)
                    res%pixel(m_r, m_c) = lhs%pixel(m_r, m_c) * rhs%pixel(m_r, m_c)
                end do
            end do
            !$OMP END PARALLEL DO

        else
            print *, "lib_image_type_operator_img_mul_img: ERROR"
            print *, "  no match of the array dimensions"
        end if

    end function lib_image_type_operator_img_mul_img

    ! Argument
    ! ----
    !   r: double precision, dimension(:,:), allocatable
    !       red channel of the image
    !   g: double precision, dimension(:,:), allocatable
    !       green channel of the image
    !   b: double precision, dimension(:,:), allocatable
    !       blue channel of the image
    !
    ! Returns
    ! ----
    !   img: image_type
    !       image with pixels of the the colour_type
    function lib_image_make_image_RGB_real(r, g, b) result(img)
        implicit none
        ! dummy
        real, dimension(:,:), allocatable, intent(in) :: r
        real, dimension(:,:), allocatable, intent(in) :: g
        real, dimension(:,:), allocatable, intent(in) :: b

        type(image_type) :: img

        ! auxiliary
        integer :: m_r
        integer :: m_c

        if (size(r) .eq. size(g) .and. size(r) .eq. size(b)) then
            allocate(img%pixel(lbound(r, 1):ubound(r, 1), &
                               lbound(r, 2):ubound(r, 2)))
            !$OMP PARALLEL DO PRIVATE(m_r, m_c)
            do m_r = lbound(r, 1), ubound(r, 1)
                do m_c = lbound(r, 2), ubound(r, 2)
                    img%pixel(m_r, m_c) = make_colour_RGB(r(m_r, m_c), g(m_r, m_c), b(m_r, m_c))
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end function lib_image_make_image_RGB_real

    ! Argument
    ! ----
    !   r: double precision, dimension(:,:), allocatable
    !       red channel of the image
    !   g: double precision, dimension(:,:), allocatable
    !       green channel of the image
    !   b: double precision, dimension(:,:), allocatable
    !       blue channel of the image
    !
    ! Returns
    ! ----
    !   img: image_type
    !       image with pixels of the the colour_type
    function lib_image_make_image_RGB_dble(r, g, b) result(img)
        implicit none
        ! dummy
        double precision, dimension(:,:), allocatable, intent(in) :: r
        double precision, dimension(:,:), allocatable, intent(in) :: g
        double precision, dimension(:,:), allocatable, intent(in) :: b

        type(image_type) :: img

        ! auxiliary
        integer :: m_r
        integer :: m_c

        if (size(r) .eq. size(g) .and. size(r) .eq. size(b)) then
            allocate(img%pixel(lbound(r, 1):ubound(r, 1), &
                               lbound(r, 2):ubound(r, 2)))
            !$OMP PARALLEL DO PRIVATE(m_r, m_c)
            do m_r = lbound(r, 1), ubound(r, 1)
                do m_c = lbound(r, 2), ubound(r, 2)
                    img%pixel(m_r, m_c) = make_colour_RGB(r(m_r, m_c), g(m_r, m_c), b(m_r, m_c))
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end function lib_image_make_image_RGB_dble

    ! Argument
    ! ----
    !   lhs: image_type
    !   rhs: image_type
    !       same size as lhs
    !
    ! Returns
    ! ----
    !   res: image_type
    !
!    function lib_image_type_operator_img_mul_RGB(img, R, G, B) result(res)
!        implicit none
!        ! dummy
!        type(image_type), intent(in) :: lhs
!        type(image_type), intent(in) :: rhs
!
!        type(image_type) :: res
!
!        ! auxiliary
!        integer :: m_r
!        integer :: m_c
!
!        if (lbound(img%pixel, 1) .eq. lbound(%pixel, 1) &
!            .and. lbound(lhs%pixel, 1) .eq. lbound(rhs%pixel, 1) &
!            ) then
!
!            allocate(res%pixel(lbound(lhs%pixel, 1):ubound(lhs%pixel, 1), &
!                               lbound(lhs%pixel, 2):ubound(lhs%pixel, 2)))
!
!            do m_r = lbound(lhs%pixel, 1), ubound(lhs%pixel, 1)
!                do m_c = lbound(lhs%pixel, 2), ubound(lhs%pixel, 2)
!                    res%pixel(m_r, m_c) = lhs%pixel(m_r, m_c) * rhs%pixel(m_r, m_c)
!                end do
!            end do
!
!        else
!            print *, "lib_image_type_operator_img_mul_img: ERROR"
!            print *, "  no match of the array dimensions"
!        end if
!
!    end function lib_image_type_operator_img_mul_img

end module lib_image_type_operator
