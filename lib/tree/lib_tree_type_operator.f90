! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

module lib_tree_type_operator
    use lib_tree_type
    implicit none

    interface operator (+)
        module procedure lib_tree_type_spatial_point_operator_add
    end interface

    interface operator (-)
        module procedure lib_tree_type_spatial_point_operator_sub
    end interface

    interface operator (*)
        module procedure lib_tree_type_spatial_point_operator_mul_scalar
        module procedure lib_tree_type_spatial_point_operator_scalar_product
    end interface

    interface operator (.eq.)
        module procedure lib_tree_type_compare_uindex_eq
    end interface

    interface operator (.ne.)
        module procedure lib_tree_type_compare_uindex_ne
    end interface

    interface abs
        module procedure lib_tree_type_spatial_point_operator_abs
    end interface

    contains

    function lib_tree_type_spatial_point_operator_add(lhs, rhs) result(rv)
        implicit none
        ! dummy
        type (lib_tree_spatial_point), intent(in) :: lhs
        type (lib_tree_spatial_point), intent(in) :: rhs
        type (lib_tree_spatial_point) :: rv

        ! auxilary
        integer :: i

        do i=1, TREE_DIMENSIONS
            rv%x(i) = lhs%x(i) + rhs%x(i)
        end do
    end function

    function lib_tree_type_spatial_point_operator_sub(lhs, rhs) result(rv)
        implicit none
        ! dummy
        type (lib_tree_spatial_point), intent(in) :: lhs
        type (lib_tree_spatial_point), intent(in) :: rhs
        type (lib_tree_spatial_point) :: rv

        ! auxilary
        integer :: i

        do i=1, TREE_DIMENSIONS
            rv%x(i) = lhs%x(i) - rhs%x(i)
        end do
    end function

    function lib_tree_type_spatial_point_operator_mul_scalar(lhs, rhs) result(rv)
        use lib_tree_type
        implicit none
        ! dummy
        type (lib_tree_spatial_point), intent(in) :: lhs
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, intent(in):: rhs
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, intent(in):: rhs
#endif
        type (lib_tree_spatial_point) :: rv

        ! auxilary
        integer :: i

        do i=1, TREE_DIMENSIONS
            rv%x(i) = lhs%x(i) * rhs
        end do
    end function

    function lib_tree_type_spatial_point_operator_scalar_mul(lhs, rhs) result(rv)
        implicit none
        ! dummy
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, intent(in):: lhs
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, intent(in):: lhs
#endif
        type (lib_tree_spatial_point), intent(in) :: rhs
        type (lib_tree_spatial_point) :: rv

        ! auxilary
        integer :: i

        do i=1, TREE_DIMENSIONS
            rv%x(i) = lhs * rhs%x(i)
        end do
    end function

    function lib_tree_type_spatial_point_operator_scalar_product(lhs, rhs) result(rv)
        implicit none
        ! dummy
        type (lib_tree_spatial_point), intent(in) :: lhs
        type (lib_tree_spatial_point), intent(in) :: rhs
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision :: rv
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real :: rv
#endif
        ! auxilary
        integer :: i

        rv = 0
        do i=1, TREE_DIMENSIONS
            rv = rv + rhs%x(i) * lhs%x(i)
        end do
    end function

    function lib_tree_type_spatial_point_operator_abs(rhs) result(rv)
        implicit none
        ! dummy
        type (lib_tree_spatial_point), intent(in) :: rhs
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision :: rv
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real :: rv
#endif
        ! auxilary
        integer :: i

        rv = 0
        do i=1, TREE_DIMENSIONS
            rv = rv + rhs%x(i)
        end do

        rv = sqrt(rv)

    end function

    function lib_tree_type_compare_uindex_eq(lhs, rhs) result(rv)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: lhs
        type(lib_tree_universal_index), intent(in) :: rhs
        logical :: rv

        rv = .true.

        if (lhs%n .ne. rhs%n) then
            rv = .false.
        end if

        if (lhs%l .ne. rhs%l) then
            rv = .false.
        end if

    end function

    function lib_tree_type_compare_uindex_ne(lhs, rhs) result(rv)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: lhs
        type(lib_tree_universal_index), intent(in) :: rhs
        logical :: rv

        rv = .true.

        if ((lhs%n .eq. rhs%n) .and. &
            (lhs%l .eq. rhs%l)) then
            rv = .false.
        end if

    end function
end module lib_tree_type_operator
