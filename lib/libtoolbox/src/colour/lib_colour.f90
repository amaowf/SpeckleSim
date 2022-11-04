module lib_colour
    use lib_colour_type
    implicit none

    private

    public :: make_colour_RGB
    public :: set_colour_RGB
    public :: get_colour_RGB
    public :: operator (*)

    interface set_colour_RGB
        module procedure lib_colour_set_RGB_real
        module procedure lib_colour_set_RGB_real_components
        module procedure lib_colour_set_RGB_dble
        module procedure lib_colour_set_RGB_dble_components
    end interface

    interface make_colour_RGB
        module procedure lib_colour_make_RGB_real
        module procedure lib_colour_make_RGB_real_components
        module procedure lib_colour_make_RGB_dble
        module procedure lib_colour_make_RGB_dble_components
    end interface

    interface get_colour_RGB
        module procedure lib_colour_get_RGB_real
        module procedure lib_colour_get_RGB_dble
    end interface

    interface operator (*)
        module procedure lib_colour_colour_mul_skalar_real
        module procedure lib_colour_colour_mul_skalar_dble

        module procedure lib_colour_skalar_mul_colour_real
        module procedure lib_colour_skalar_mul_colour_dble

        module procedure lib_colour_colour_mul_colour_RGB
    end interface

    contains

        ! Argument
        ! ----
        !   rgb: real, dimension(3)
        !       vector with the RGB colour values [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        subroutine lib_colour_set_RGB_real(colour, rgb)
            implicit none
            ! dummy
            type(colour_type), intent(inout) :: colour
            real, dimension(3), intent(in) :: rgb

            colour%rgb = dble(rgb)
        end subroutine

        ! Argument
        ! ----
        !   red: real
        !       red component [0, 255]
        !   green: real
        !       green component [0, 255]
        !   blue: real
        !       blue component [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        subroutine lib_colour_set_RGB_real_components(colour, red, green, blue)
            implicit none
            ! dummy
            type(colour_type), intent(inout) :: colour
            real, intent(in) :: red
            real, intent(in) :: green
            real, intent(in) :: blue

            colour%rgb = dble((/red, green, blue/))
        end subroutine

        ! Argument
        ! ----
        !   rgb: double precision, dimension(3)
        !       vector with the RGB colour values [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        subroutine lib_colour_set_RGB_dble(colour, rgb)
            implicit none
            ! dummy
            type(colour_type), intent(inout) :: colour
            double precision, dimension(3), intent(in) :: rgb

            colour%rgb = dble(rgb)
        end subroutine

        ! Argument
        ! ----
        !   red: double precision
        !       red component [0, 255]
        !   green: double precision
        !       green component [0, 255]
        !   blue: double precision
        !       blue component [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        subroutine lib_colour_set_RGB_dble_components(colour, red, green, blue)
            implicit none
            ! dummy
            type(colour_type), intent(inout) :: colour
            double precision, intent(in) :: red
            double precision, intent(in) :: green
            double precision, intent(in) :: blue

            colour%rgb = dble((/red, green, blue/))
        end subroutine

        ! Argument
        ! ----
        !   colour: colour_type
        !       Color object from which the RGB colors are extracted.
        !
        ! Retunrs
        ! ----
        !   rv: real, dimension(3)
        !       vector with the RBG colour values [0, 255]
        subroutine lib_colour_get_RGB_real(colour, rv)
            implicit none
            ! dummy
            type(colour_type), intent(in) :: colour

            real, dimension(3), intent(out) :: rv

            rv = real(colour%rgb)
        end subroutine

        ! Argument
        ! ----
        !   colour: colour_type
        !       Color object from which the RGB colors are extracted.
        !
        ! Retunrs
        ! ----
        !   rv: real, dimension(3)
        !       vector with the RBG colour values [0, 255]
        subroutine lib_colour_get_RGB_dble(colour, rv)
            implicit none
            ! dummy
            type(colour_type), intent(in) :: colour

            double precision, dimension(3), intent(out) :: rv

            rv = dble(colour%rgb)
        end subroutine

        ! Argument
        ! ----
        !   rgb: real, dimension(3)
        !       vector with the RGB colour values [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        function lib_colour_make_RGB_real(rgb) result(colour)
            implicit none
            ! dummy
            real, dimension(3), intent(in) :: rgb

            type(colour_type) :: colour

            call set_colour_RGB(colour, rgb)
        end function

        ! Argument
        ! ----
        !   red: real
        !       red component [0, 255]
        !   green: real
        !       green component [0, 255]
        !   blue: real
        !       blue component [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        function lib_colour_make_RGB_real_components(red, green, blue) result(colour)
            implicit none
            ! dummy
            real, intent(in) :: red
            real, intent(in) :: green
            real, intent(in) :: blue

            type(colour_type) :: colour

            call set_colour_RGB(colour, red, green, blue)
        end function

        ! Argument
        ! ----
        !   rgb: double precision, dimension(3)
        !       vector with the RGB colour values [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        function lib_colour_make_RGB_dble(rgb) result(colour)
            implicit none
            ! dummy
            double precision, dimension(3), intent(in) :: rgb

            type(colour_type) :: colour

            call set_colour_RGB(colour, rgb)
        end function

        ! Argument
        ! ----
        !   red: double precision
        !       red component [0, 255]
        !   green: double precision
        !       green component [0, 255]
        !   blue: double precision
        !       blue component [0, 255]
        !
        ! Returns
        ! ----
        !   colour: colour_type
        !       Color object where the RGB colours are set.
        function lib_colour_make_RGB_dble_components(red, green, blue) result(colour)
            implicit none
            ! dummy
            double precision, intent(in) :: red
            double precision, intent(in) :: green
            double precision, intent(in) :: blue

            type(colour_type) :: colour

            call set_colour_RGB(colour, red, green, blue)
        end function

        function lib_colour_colour_mul_skalar_real(c, s) result(res)
            implicit none
            ! dummy
            type(colour_type), intent(in) :: c
            real, intent(in) :: s

            type(colour_type) :: res

            ! auxiliary
            double precision, dimension(3) :: rgb

            call get_colour_RGB(c, rgb)

            rgb = rgb * dble(s)

            call lib_colour_set_RGB_dble(res, rgb)

        end function lib_colour_colour_mul_skalar_real

        function lib_colour_colour_mul_skalar_dble(c, s) result(res)
            implicit none
            ! dummy
            type(colour_type), intent(in) :: c
            double precision, intent(in) :: s

            type(colour_type) :: res

            ! auxiliary
            double precision, dimension(3) :: rgb

            call get_colour_RGB(c, rgb)

            rgb = rgb * dble(s)

            call lib_colour_set_RGB_dble(res, rgb)

        end function lib_colour_colour_mul_skalar_dble

        function lib_colour_skalar_mul_colour_real(s, c) result(res)
            implicit none
            ! dummy
            real, intent(in) :: s
            type(colour_type), intent(in) :: c

            type(colour_type) :: res

            res = c * s

        end function lib_colour_skalar_mul_colour_real

        function lib_colour_skalar_mul_colour_dble(s, c) result(res)
            implicit none
            ! dummy
            double precision, intent(in) :: s
            type(colour_type), intent(in) :: c

            type(colour_type) :: res

            res = c * s

        end function lib_colour_skalar_mul_colour_dble

        function lib_colour_colour_mul_colour_RGB(lhs, rhs) result(res)
            implicit none
            ! dummy
            type(colour_type), intent(in) :: lhs
            type(colour_type), intent(in) :: rhs

            type(colour_type) :: res

            res%RGB = lhs%RGB * rhs%RGB

        end function lib_colour_colour_mul_colour_RGB

end module lib_colour

