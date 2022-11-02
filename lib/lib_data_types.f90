module lib_data_types
    implicit none

    integer, parameter :: element_kind = 4

    ! list element
    type lib_tree_element
        integer(kind=element_kind)   :: start
        integer(kind=element_kind)   :: end
        integer(kind=element_kind)   :: parent_element
    end type lib_tree_element

!    type tree_parent
!        type(tree_element), dimension(8) :: children
!        type(ectree_element), dimension() :: neighbor
!    end type tree_level



end module lib_data_types
