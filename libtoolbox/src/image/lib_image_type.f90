module lib_image_type
    use lib_colour_type
    implicit none

    type image_type
        type(colour_type), dimension(:,:), allocatable :: pixel
    end type

end module lib_image_type
