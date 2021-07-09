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
! Created on Thu May 28 10:09:51 2020
!
! @author: Max Daiber-Huppert
!

!#define DEBUG
!#define TRACE

module lib_tree_traversal
    use libmath
    use lib_tree
    use lib_tree_type
    implicit none

    private

    public :: lib_tree_traversal_get_traversed_boxes
    public :: lib_tree_traversal_get_callback_functions

    public :: lib_tree_traversal_test_functions

    public :: lib_tree_traversal_callback_type

    integer, parameter :: LIST_ALLOCATION_SETP = 10

    interface
        subroutine lib_tree_traversal_octree_get_box_size(uindex, box_size, box_min, box_max)
            use libmath
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(cartesian_coordinate_real_type), intent(out) :: box_size
            type(cartesian_coordinate_real_type), intent(out) :: box_min
            type(cartesian_coordinate_real_type), intent(out) :: box_max
        end subroutine

        function lib_tree_traversal_octree_is_terminal(uindex, data) result (rv)
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_tree_traversal_data), intent(inout) :: data
            logical :: rv
        end function

        subroutine lib_tree_traversal_travesed_box(uindex, data)
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_tree_traversal_data), intent(inout) :: data
        end subroutine

    end interface

    type lib_tree_traversal_callback_type
        procedure(lib_tree_traversal_octree_get_box_size), pointer, nopass :: octree_get_box_size => null()
        procedure(lib_tree_traversal_octree_is_terminal), pointer, nopass :: is_box_terminal => null()
        procedure(lib_tree_traversal_travesed_box), pointer, nopass :: traversed_box => null()
    end type

    integer(kind=1), parameter :: END = 8

contains

    ! Argument
    ! ----
    !   ray_origin: cartesian_coordinate_real_type
    !       origin of the ray
    !   ray_direction: cartesian_coordinate_real_type
    !       unit length direction vector
    !   callback: lib_tree_traversal_callback_type
    !       If *lib_tree* is used to build the tree hierarchy, only
    !       *is_box_terminal* has to be defined.
    !       This implies a box the an edge length of 1 in the first quadrant.
    !       With a vertex in the origin. The remaining functions can be set
    !       by *lib_tree_traversal_get_callback_functions*.
    !
    ! Returns
    ! ----
    !   data: lib_tree_traversal_data
    !       contains a list of all fields that the ray passes through
    subroutine lib_tree_traversal_get_traversed_boxes(callback, data)
        implicit none
        ! dummy
        type(lib_tree_traversal_callback_type), intent(in) :: callback
        type(lib_tree_traversal_data), intent(inout) :: data

        ! auxiliary
        type(lib_tree_spatial_point) :: buffer
        type(cartesian_coordinate_real_type) :: buffer_2
        type(lib_tree_universal_index), dimension(:), allocatable :: traversed_box_list


        buffer%x = data%ray_origin
        buffer = lib_tree_get_scaled_point(buffer)
        buffer_2 = buffer%x

        call ray_parameter(buffer_2, data%ray_direction, callback, data)

        if (allocated(data%traversed_box_list)) then
            allocate(traversed_box_list(1:data%index))

            traversed_box_list = data%traversed_box_list(1:data%index)

            call move_alloc(traversed_box_list, data%traversed_box_list)

        end if

    end subroutine lib_tree_traversal_get_traversed_boxes

    ! Argument
    ! ----
    !   ray_origin: cartesian_coordinate_real_type
    !       origin of the ray
    !   ray_direction: cartesian_coordinate_real_type
    !       unit length direction vector
    !   callback: lib_tree_traversal_callback_type
    !       If *lib_tree* is used to build the tree hierarchy, only
    !       *is_box_terminal* has to be defined.
    !       This implies a box the an edge length of 1 in the first quadrant.
    !       With a vertex in the origin. The remaining functions can be set
    !       by *lib_tree_traversal_get_callback_functions*.
    !
    ! Returns
    ! ----
    !   data: lib_tree_traversal_data
    !       contains a list of all fields that the ray passes through
    subroutine lib_tree_traversal_traverse_tree(callback, data)
        implicit none
        ! dummy
        type(lib_tree_traversal_callback_type), intent(in) :: callback
        type(lib_tree_traversal_data), intent(inout) :: data

        ! auxiliary
        type(lib_tree_traversal_data) :: m_data
        type(lib_tree_spatial_point) :: buffer

        m_data = data

        buffer%x = m_data%ray_origin
        buffer = lib_tree_get_scaled_point(buffer)
        m_data%ray_origin = buffer%x

        call ray_parameter(m_data%ray_origin, m_data%ray_direction, callback, m_data)

        if (allocated(m_data%traversed_box_list)) then
            if (allocated(data%traversed_box_list)) deallocate(data%traversed_box_list)
            allocate(data%traversed_box_list(m_data%index))
            data%traversed_box_list = m_data%traversed_box_list(1:m_data%index)
        else
            if (allocated(data%traversed_box_list)) deallocate(data%traversed_box_list)
        end if
        call move_alloc(m_data%traversed_box_list, data%traversed_box_list)

    end subroutine lib_tree_traversal_traverse_tree

    ! Argument
    ! ----
    !   ray_origin: cartesian_coordinate_real_type
    !       origin of the ray
    !   ray_direction: cartesian_coordinate_real_type
    !       unit length direction vector
    !   callback: lib_tree_traversal_callback_type
    !       If *lib_tree* is used to build the tree hierarchy, only
    !       *is_box_terminal* has to be defined.
    !       This implies a box the an edge length of 1 in the first quadrant.
    !       With a vertex in the origin. The remaining functions can be set
    !       by *lib_tree_traversal_get_callback_functions*.
    !
    ! Returns
    ! ----
    !   data: lib_tree_traversal_data
    !       contains a list of all fields that the ray passes through
    !
    ! Reference: An Efﬁcient Parametric Algorithm for Octree Traversal, Jorge Revelles
    !            pseudocode
    subroutine ray_parameter(ray_origin, ray_direction, callback, data)
        implicit none
        ! dummy
        type(cartesian_coordinate_real_type), intent(in) :: ray_origin
        type(cartesian_coordinate_real_type), intent(in) :: ray_direction
        type(lib_tree_traversal_callback_type), intent(in) :: callback
        type(lib_tree_traversal_data), intent(inout) :: data

        ! auxiliary
        integer(kind=1) :: a
        type(lib_tree_universal_index) :: m_uindex
        type(cartesian_coordinate_real_type) :: m_ray_origin
        type(cartesian_coordinate_real_type) :: m_ray_direction
        type(cartesian_coordinate_real_type) :: m_box_size
        type(cartesian_coordinate_real_type) :: m_box_min
        type(cartesian_coordinate_real_type) :: m_box_max

        type(cartesian_coordinate_real_type) :: m_t0
        type(cartesian_coordinate_real_type) :: m_t1

        logical :: buffer
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_list

        m_ray_origin = ray_origin
        m_ray_direction = ray_direction
        m_uindex%l = 0
        m_uindex%n = 0
        if (associated(callback%octree_get_box_size)) then
            call callback%octree_get_box_size(m_uindex, m_box_size, m_box_min, m_box_max)
        else
            print *, "lib_tree_traversal / ray_parameter: ERROR"
            print *, "  callback%octree_get_box_size isn't associated"
            return
        end if

        if (.not. associated(callback%is_box_terminal)) then
            print *, "lib_tree_traversal / ray_parameter: ERROR"
            print *, "  callback%is_box_terminal isn't associated"
            return
        end if

        if (.not. associated(callback%traversed_box)) then
            print *, "lib_tree_traversal / ray_parameter: ERROR"
            print *, "  callback%traversed_box isn't associated"
            return
        end if

        a = 0_1

        if (m_ray_direction%x .lt. 0d0) then
            m_ray_origin%x = m_box_size%x - m_ray_origin%x
            m_ray_direction%x = - m_ray_direction%x
            a = xor(a, 4_1)
        end if

        if (m_ray_direction%y .lt. 0d0) then
            m_ray_origin%y = m_box_size%y - m_ray_origin%y
            m_ray_direction%y = - m_ray_direction%y
            a = xor(a, 2_1)
        end if

        if (m_ray_direction%z .lt. 0d0) then
            m_ray_origin%z = m_box_size%z - m_ray_origin%z
            m_ray_direction%z = - m_ray_direction%z
            a = xor(a, 1_1)
        end if

        m_t0%x = (m_box_min%x - m_ray_origin%x) / m_ray_direction%x
        m_t1%x = (m_box_max%x - m_ray_origin%x) / m_ray_direction%x

        m_t0%y = (m_box_min%y - m_ray_origin%y) / m_ray_direction%y
        m_t1%y = (m_box_max%y - m_ray_origin%y) / m_ray_direction%y

        m_t0%z = (m_box_min%z - m_ray_origin%z) / m_ray_direction%z
        m_t1%z = (m_box_max%z - m_ray_origin%z) / m_ray_direction%z

        if (max(m_t0%x, m_t0%x, m_t0%z) .ge. min(m_t1%x, m_t1%y, m_t1%z)) then
            return
        end if

        ! --- set ray origin outside the octree area
        if (m_ray_origin%x .ge. m_box_min%x &
            .and. m_ray_origin%y .ge. m_box_min%y &
            .and. m_ray_origin%z .ge. m_box_min%z) then

            m_ray_origin = m_ray_origin - m_ray_direction * abs(m_box_size)
        end if
        ! ~~~

        m_t0%x = (m_box_min%x - m_ray_origin%x) / m_ray_direction%x
        m_t1%x = (m_box_max%x - m_ray_origin%x) / m_ray_direction%x

        m_t0%y = (m_box_min%y - m_ray_origin%y) / m_ray_direction%y
        m_t1%y = (m_box_max%y - m_ray_origin%y) / m_ray_direction%y

        m_t0%z = (m_box_min%z - m_ray_origin%z) / m_ray_direction%z
        m_t1%z = (m_box_max%z - m_ray_origin%z) / m_ray_direction%z

        ! exclude rays along a box face if m_t*%* = 0 / 0
        if (isnan(m_t0%x) .or. isnan(m_t1%x) &
            .or. isnan(m_t0%y) .or. isnan(m_t1%y) &
            .or. isnan(m_t0%z) .or. isnan(m_t1%z)) then
            return
        end if

        if (max(m_t0%x, m_t0%x, m_t0%z) .lt. min(m_t1%x, m_t1%y, m_t1%z)) then
            call callback%traversed_box(m_uindex, data)
            buffer = proc_subtree(m_t0, m_t1, m_uindex, callback, a, m_ray_origin, m_ray_direction, data)

            if (allocated(data%traversed_box_list)) then
                allocate(buffer_list(data%index))
                buffer_list = data%traversed_box_list(1:data%index)
                call move_alloc(buffer_list, data%traversed_box_list)
            end if
        end if

    end subroutine ray_parameter

    ! Reference: An Efﬁcient Parametric Algorithm for Octree Traversal, Jorge Revelles
    !            pseudocode
    recursive function proc_subtree(t0, t1, uindex, callback, a, ray_origin, ray_direction, data) result(rv)
        implicit none
        ! dummy
        type(cartesian_coordinate_real_type), intent(in) :: t0
        type(cartesian_coordinate_real_type), intent(in) :: t1
!        type(cartesian_coordinate_real_type), intent(in) :: tm
        type(lib_tree_universal_index), intent(in) :: uindex
        type(lib_tree_traversal_callback_type), intent(in) :: callback
        integer(kind=1), intent(in) :: a
        type(cartesian_coordinate_real_type), intent(in) :: ray_origin
        type(cartesian_coordinate_real_type), intent(in) :: ray_direction
        type(lib_tree_traversal_data), intent(inout) :: data

        logical :: rv

        ! auxiliary
        integer(kind=1) :: n
        integer(kind=1) :: m_current_box
        type(lib_tree_universal_index) :: m_uindex
        type(cartesian_coordinate_real_type) :: m_t0
        type(cartesian_coordinate_real_type) :: m_t1
        type(cartesian_coordinate_real_type) :: m_tm
        type(cartesian_coordinate_real_type) :: m_box_size
        type(cartesian_coordinate_real_type) :: m_box_min
        type(cartesian_coordinate_real_type) :: m_box_max

        m_uindex = uindex

        m_t0 = t0
        m_t1 = t1

        rv = .true.

        if (m_t1%x .lt. 0d0 &
            .or. m_t1%y .lt. 0d0 &
            .or. m_t1%z .lt. 0d0) then
#ifdef TRACE
            print *, "TRACE proc_subtree: t1 .lt. 0"
#endif
            rv = .false.
            return
        end if

        if (.not. lib_tree_traversal_is_on_ray_path(uindex, data%ray_origin, data%ray_direction)) then
#ifdef TRACE
            print *, "TRACE proc_subtree: lib_tree_traversal_is_on_ray_path .eqv. .false."
            print *, "  uindex: ", uindex%l, uindex%n
#endif
            rv = .false.
            return
        end if


        if (callback%is_box_terminal(uindex, data)) then
#ifdef TRACE
            print *, "TRACE proc_subtree: is terminal"
            print *, "  uindex: ", uindex%l, uindex%n
#endif
            ! box has no child boxes
            call callback%traversed_box(m_uindex, data)
            rv = .false.
            return
        end if

        call callback%octree_get_box_size(m_uindex, m_box_size, m_box_min, m_box_max)

        ! eq. (12)
        if (isinf_neg(m_t0%x) .and. isinf_pos(m_t1%x)) then
            if (ray_origin%x .lt. (m_box_min%x + m_box_max%x) / 2) then
                m_tm%x = get_inf_pos()
            else
                m_tm%x = get_inf_neg()
            end if
        else
            m_tm%x = 0.5d0 * (m_t0%x + m_t1%x)
        end if

        if (isinf_neg(m_t0%y) .and. isinf_pos(m_t1%y)) then
            if (ray_origin%y .lt. (m_box_min%y + m_box_max%y) / 2) then
                m_tm%y = get_inf_pos()
            else
                m_tm%y = get_inf_neg()
            end if
        else
            m_tm%y = 0.5d0 * (m_t0%y + m_t1%y)
        end if

        if (isinf_neg(m_t0%z) .and. isinf_pos(m_t1%z)) then
            if (ray_origin%z .lt. (m_box_min%z + m_box_max%z) / 2) then
                m_tm%z = get_inf_pos()
            else
                m_tm%z = get_inf_neg()
            end if
        else
            m_tm%z = 0.5d0 * (m_t0%z + m_t1%z)
        end if

        m_current_box = first_box(m_t0, m_tm)
#ifdef TRACE
        print *, "m_current_box = ", m_current_box
#endif

        do
            select case(m_current_box)
                case(0)
                    m_current_box = new_box(m_tm%x, 4_1, m_tm%y, 2_1, m_tm%z, 1_1)

                    m_t0%x = t0%x
                    m_t0%y = t0%y
                    m_t0%z = t0%z

                    m_t1%x = m_tm%x
                    m_t1%y = m_tm%y
                    m_t1%z = m_tm%z

                    n = a

                case(1)
                    m_current_box = new_box(m_tm%x, 5_1, m_tm%y, 3_1, t1%z, END)

                    m_t0%x = t0%x
                    m_t0%y = t0%y
                    m_t0%z = m_tm%z

                    m_t1%x = m_tm%x
                    m_t1%y = m_tm%y
                    m_t1%z = t1%z

                    n = xor(1_1, a)
                case(2)
                    m_current_box = new_box(m_tm%x, 6_1, t1%y, END, m_tm%z, 3_1)

                    m_t0%x = t0%x
                    m_t0%y = m_tm%y
                    m_t0%z = t0%z

                    m_t1%x = m_tm%x
                    m_t1%y = t1%y
                    m_t1%z = m_tm%z

                    n = xor(2_1, a)
                case(3)
                    m_current_box = new_box(m_tm%x, 7_1, t1%y, END, t1%z, END)

                    m_t0%x = t0%x
                    m_t0%y = m_tm%y
                    m_t0%z = m_tm%z

                    m_t1%x = m_tm%x
                    m_t1%y = t1%y
                    m_t1%z = t1%z

                    n = xor(3_1, a)
                case(4)
                    m_current_box = new_box(t1%x, END, m_tm%y, 6_1, m_tm%z, 5_1)

                    m_t0%x = m_tm%x
                    m_t0%y = t0%y
                    m_t0%z = t0%z

                    m_t1%x = t1%x
                    m_t1%y = m_tm%y
                    m_t1%z = m_tm%z

                    n = xor(4_1, a)
                case(5)
                    m_current_box = new_box(t1%x, END, m_tm%y, 7_1, t1%z, END)

                    m_t0%x = m_tm%x
                    m_t0%y = t0%y
                    m_t0%z = m_tm%z

                    m_t1%x = t1%x
                    m_t1%y = m_tm%y
                    m_t1%z = t1%z

                    n = xor(5_1, a)
                case(6)
                    m_current_box = new_box(t1%x, END, t1%y, END, m_tm%z, 7_1)

                    m_t0%x = m_tm%x
                    m_t0%y = m_tm%y
                    m_t0%z = t0%z

                    m_t1%x = t1%x
                    m_t1%y = t1%y
                    m_t1%z = m_tm%z

                    n = xor(6_1, a)
                case(7)
                    m_current_box = END

                    m_t0%x = m_tm%x
                    m_t0%y = m_tm%y
                    m_t0%z = m_tm%z

                    m_t1%x = t1%x
                    m_t1%y = t1%y
                    m_t1%z = t1%z

                    n = xor(7_1, a)
                case default
                    print *, "lib_tree_traversal / proc_subtree: ERROR"
                    print *, "  unexpected current_box: ", m_current_box
            end select
#ifdef TRACE
            print *, "new box: ", m_current_box
            print *, "n: ", n
            print *, "l: ", uindex%l + 1_1
#endif
            m_uindex%n = uindex%n * int(2, kind=UINDEX_BYTES)**(TREE_DIMENSIONS) + n
            m_uindex%l = uindex%l + 1_1

            if (proc_subtree(m_t0, m_t1, m_uindex, callback, a, ray_origin, ray_direction, data)) then
!                m_uindex%n = uindex%n * int(2, kind=UINDEX_BYTES)**(TREE_DIMENSIONS) + xor(n, a)
                call callback%traversed_box(m_uindex, data)
            end if

            if (m_current_box .ge. END) then
                exit
            end if
        end do

    end function proc_subtree

    ! Reference: An Efﬁcient Parametric Algorithm for Octree Traversal, Jorge Revelles
    !            Table 1 & 2
    function first_box(t0, tm) result(rv)
        implicit none
        ! dummy
        type(cartesian_coordinate_real_type), intent(in) :: t0
        type(cartesian_coordinate_real_type), intent(in) :: tm
        integer(kind=1) :: rv

        integer, parameter :: XY = 0
        integer, parameter :: XZ = 1
        integer, parameter :: YZ = 2

        ! auxiliary
        integer :: entry_plane

        rv = 0

        ! Table 2
        entry_plane = XY
        if ( t0%x .gt. t0%z) then
            entry_plane = YZ
            if (t0%y .gt. t0%x) then
                entry_plane = XZ
            end if
        else if ( t0%y .gt. t0%z) then
            entry_plane = XZ
        end if

        ! Table 1
        select case(entry_plane)
            case(XY)
#ifdef TRACE
                print *, "Entry Plane: XY"
#endif
                if (tm%x .lt. t0%z) then
                    rv = or(rv, 4_1) !or(rv, 1_1)
                end if

                if (tm%y .lt. t0%z) then
                    rv =  or(rv, 2_1) !or(rv, 2_1)
                end if
            case(XZ)
#ifdef TRACE
                print *, "Entry Plane: XZ"
#endif
                if (tm%x .lt. t0%y) then
                    rv = or(rv, 4_1) !PAPER-ORG: or(rv, 1_1)
                end if

                if (tm%z .lt. t0%y) then
                    rv = or(rv, 1_1) !PAPER-ORG: or(rv, 4_1)
                end if
            case(YZ)
#ifdef TRACE
                print *, "Entry Plane: YZ"
#endif
                if (tm%y .lt. t0%x) then
                    rv = or(rv, 2_1) !PAPER-ORG: or(rv, 2_1)
                end if

                if (tm%z .lt. t0%x) then
                    rv = or(rv, 1_1) !PAPER-ORG: or(rv, 4_1)
                end if
            case default
                print *, "lib_tree_traversal / first_box: ERROR"
                print *, "  unexpected value for entry_plane: ", entry_plane
        end select
    end function first_box

    ! Reference: An Efﬁcient Parametric Algorithm for Octree Traversal, Jorge Revelles
    !            Table 3 and §3.2 Obtaining the Next box
    function new_box(a1, i1, a2, i2, a3, i3) result(rv)
        implicit none
        ! dummy
        double precision, intent(in) :: a1
        integer(kind=1), intent(in) :: i1
        double precision, intent(in) :: a2
        integer(kind=1), intent(in) :: i2
        double precision, intent(in) :: a3
        integer(kind=1), intent(in) :: i3

        integer(kind=1) :: rv

        !auxiliary
        double precision :: buffer

        rv = i3
        buffer = a3
        if (a1 .lt. buffer) then
            rv = i1
            buffer = a1
        end if

        if (a2 .lt. buffer) then
            rv = i2
        end if

    end function new_box

    !
    function lib_tree_traversal_is_on_ray_path(uindex, ray_origin, ray_direction) result(rv)
        use lib_tree
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        type(cartesian_coordinate_real_type), intent(in) :: ray_origin
        type(cartesian_coordinate_real_type), intent(in) :: ray_direction

        logical :: rv

        ! auxiliary
        type(cartesian_coordinate_real_type) :: point
        type(lib_tree_spatial_point) :: buffer_point
        double precision :: buffer


        buffer_point%x = ray_origin
        if (lib_tree_is_in_box(uindex, buffer_point)) then
            rv = .true.
        else
            buffer_point = lib_tree_get_centre_of_box(uindex)
            point = buffer_point%x
            point = point - ray_origin

            buffer = dot_product(point, ray_direction)

            if (buffer .gt. 0d0) then
                rv = .true.
            else
                rv = .false.
            end if
        end if

    end function lib_tree_traversal_is_on_ray_path

    ! Argument
    ! ----
    !   uindex: lib_tree_universal_index
    !       universal index of the box
    !
    ! Retruns
    ! ----
    !   box_size: cartesian_coordinate_real_type
    !       expansion of the box in the x, y and z direction
    !   box_min: cartesian_coordinate_real_type
    !       coordinate of the box vertex with the shortest distance to the origin
    !   box_max: cartesian_coordinate_real_type
    !       coordinate of the box vertex with the longest distance to the origin
    subroutine lib_tree_traversal_get_box_size_std(uindex, box_size, box_min, box_max)
            use libmath
            use lib_tree_type
            use lib_tree
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(cartesian_coordinate_real_type), intent(out) :: box_size
            type(cartesian_coordinate_real_type), intent(out) :: box_min
            type(cartesian_coordinate_real_type), intent(out) :: box_max

            ! auxiliary
            type(lib_tree_spatial_point) :: centre
            double precision :: edge_length


            centre = lib_tree_get_centre_of_box(uindex)
            edge_length = lib_tree_get_box_edge_length(uindex%l)

            box_size = make_cartesian(edge_length, edge_length, edge_length)

            box_min = make_cartesian(centre%x(1) - edge_length / 2d0, &
                                      centre%x(2) - edge_length / 2d0, &
                                      centre%x(3) - edge_length / 2d0)

            box_max = make_cartesian(centre%x(1) + edge_length / 2d0, &
                                      centre%x(2) + edge_length / 2d0, &
                                      centre%x(3) + edge_length / 2d0)

        end subroutine lib_tree_traversal_get_box_size_std

        ! This function adds the *uindex* to the list contained in *data*.
        !
        ! Argument
        ! ----
        !   uindex: lib_tree_universal_index
        !       universal index of the box through which the ray passes
        !   data: lib_tree_traversal_data
        !       contains the dynamically allocated list of the boxes passed through
        subroutine lib_tree_traversal_travesed_box_std(uindex, data)
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_tree_traversal_data), intent(inout) :: data

            ! auxiliary
            type(lib_tree_universal_index), dimension(:), allocatable :: buffer_list

#ifdef DEBUG
            print *, "---"
            print *, "traversed box (l, n):", uindex%l, uindex%n
#endif

            if (allocated(data%traversed_box_list)) then
                if (ubound(data%traversed_box_list, 1) .le. data%index) then
                    allocate(buffer_list(ubound(data%traversed_box_list, 1) + LIST_ALLOCATION_SETP))
                    buffer_list(1:ubound(data%traversed_box_list, 1)) = data%traversed_box_list
                    call move_alloc(buffer_list, data%traversed_box_list)
                end if
            else
                allocate(data%traversed_box_list(LIST_ALLOCATION_SETP))
                data%index = 0
            end if

            data%traversed_box_list(data%index + 1) = uindex

            data%index = data%index + 1

        end subroutine lib_tree_traversal_travesed_box_std

        ! Argument
        ! ----
        !   uindex: lib_tree_universal_index
        !       universal index of the box
        !
        ! Returns
        ! ----
        !   rv: logical
        !       true: box is terminal -> no further child boxes
        !       false: box is NOG terminal -> there are further child boxes
        function lib_tree_traversal_is_terminal_test(uindex, data) result (rv)
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_tree_traversal_data), intent(inout) :: data
            logical :: rv

!            print *, "lib_tree_traversal_is_terminal_test:"
!            print *, "  Please define an own function to check if a box is terminal."
!            print *, "  This is only a test function."

            if (uindex%l .ge. 2) then
                rv = .true.
            else
                rv = .false.
            end if

        end function

        ! Resturns
        !   rv: lib_tree_traversal_callback_type
        !       set of standard callback functions
        !       - octree_get_box_size: lib_tree_traversal_get_box_size
        !       - traversed_box: lib_tree_traversal_travesed_box
        !       HINT:
        !           By calling *lib_tree_traversal* functions outside this module,
        !           *is_box_terminal* must point to a function specified for the
        !           use case.
        function lib_tree_traversal_get_callback_functions() result(rv)
            implicit none
            ! dummy
            type(lib_tree_traversal_callback_type) :: rv

            rv%octree_get_box_size => lib_tree_traversal_get_box_size_std
            rv%traversed_box => lib_tree_traversal_travesed_box_std
            rv%is_box_terminal => lib_tree_traversal_is_terminal_test

        end function

    function lib_tree_traversal_test_functions() result(rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        if (.not. test_lib_tree_traversal_get_traversed_boxes()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_2()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_3()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_4()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_5()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_6()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_7()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_8()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_9()) rv = rv + 1
        if (.not. test_lib_tree_traversal_get_traversed_boxes_10()) rv = rv + 1

    contains

        function test_lib_tree_traversal_get_traversed_boxes() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(6) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 33

            ground_truth(3)%l = 1
            ground_truth(3)%n = 4

            ground_truth(4)%l = 2
            ground_truth(4)%n = 40

            ground_truth(5)%l = 2
            ground_truth(5)%n = 41

            ground_truth(6)%l = 1
            ground_truth(6)%n = 5

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     | o   |     |  o: ray
            !           |-------------------------
            !           |--> x
            !
            data%ray_origin = make_cartesian(0.55d0, 0.1d0, 0.4d0)
            data%ray_direction = make_spherical(1d0, PI / 32d0, PI / 32d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes

        function test_lib_tree_traversal_get_traversed_boxes_2() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(7) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 52

            ground_truth(3)%l = 2
            ground_truth(3)%n = 48

            ground_truth(4)%l = 1
            ground_truth(4)%n = 6

            ground_truth(5)%l = 2
            ground_truth(5)%n = 20

            ground_truth(6)%l = 2
            ground_truth(6)%n = 16

            ground_truth(7)%l = 1
            ground_truth(7)%n = 2

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |<-- ray
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x
            !
            data%ray_origin = make_cartesian(1.1d0, 0.6d0, 0.15d0)
            data%ray_direction = make_spherical(1d0, PI / 2.01d0, 0.999*PI) ! make_cartesian(-1d0, 0d0, 0d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_2:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_2: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_2

        function test_lib_tree_traversal_get_traversed_boxes_3() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(7) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 50

            ground_truth(3)%l = 2
            ground_truth(3)%n = 48

            ground_truth(4)%l = 1
            ground_truth(4)%n = 6

            ground_truth(5)%l = 2
            ground_truth(5)%n = 34

            ground_truth(6)%l = 2
            ground_truth(6)%n = 32

            ground_truth(7)%l = 1
            ground_truth(7)%n = 4

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x         | ray
            !                          v
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x
            !
            data%ray_origin = make_cartesian(0.6d0, 1.1d0, 0.15d0)
            data%ray_direction = make_cartesian(0d0, -1d0, 0d0) ! make_spherical(1d0, PI / 2.01d0, 3.01d0/2d0*PI)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_3:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_3: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_3

        function test_lib_tree_traversal_get_traversed_boxes_4() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(7) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 32

            ground_truth(3)%l = 2
            ground_truth(3)%n = 34

            ground_truth(4)%l = 1
            ground_truth(4)%n = 4

            ground_truth(5)%l = 2
            ground_truth(5)%n = 48

            ground_truth(6)%l = 2
            ground_truth(6)%n = 50

            ground_truth(7)%l = 1
            ground_truth(7)%n = 6

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x         ^
            !                          | ray
            !
            data%ray_origin = make_cartesian(0.6d0, -0.1d0, 0.15d0)
            data%ray_direction = make_spherical(1d0, PI / 2d0, PI/1.99d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_4:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_4: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_4

        function test_lib_tree_traversal_get_traversed_boxes_5() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(5) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 4

            ground_truth(3)%l = 2
            ground_truth(3)%n = 0

            ground_truth(4)%l = 2
            ground_truth(4)%n = 2

            ground_truth(5)%l = 1
            ground_truth(5)%n = 0

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x      ^
            !                        \ ray
            !
            data%ray_origin = make_cartesian(0.55d0, -0.1d0, 0.15d0)
            data%ray_direction = make_spherical(1d0, PI / 2d0, 3d0 / 4d0 * PI)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_5:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_5: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_5

        function test_lib_tree_traversal_get_traversed_boxes_6() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(5) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 4

            ground_truth(3)%l = 2
            ground_truth(3)%n = 0

            ground_truth(4)%l = 2
            ground_truth(4)%n = 2

            ground_truth(5)%l = 1
            ground_truth(5)%n = 0

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     |
            !            ------0-----------4------
            !            |  0 r|  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x      ^
            !                    r:  \ ray
            !
            data%ray_origin = make_cartesian(0.3d0, 0.15d0, 0.15d0)
            data%ray_direction = make_spherical(1d0, PI / 2d0, 3d0 / 4d0 * PI)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_6:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_6: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_6

        function test_lib_tree_traversal_get_traversed_boxes_7() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(9) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 34

            ground_truth(3)%l = 1
            ground_truth(3)%n = 4

            ground_truth(4)%l = 2
            ground_truth(4)%n = 6

            ground_truth(5)%l = 1
            ground_truth(5)%n = 0

            ground_truth(6)%l = 2
            ground_truth(6)%n = 20

            ground_truth(7)%l = 2
            ground_truth(7)%n = 16

            ground_truth(8)%l = 2
            ground_truth(8)%n = 18

            ground_truth(9)%l = 1
            ground_truth(9)%n = 2

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |r    |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x      ^
            !                    r:  \ ray
            !
            data%ray_origin = make_cartesian(0.6d0, 0.3d0, 0.15d0)
            data%ray_direction = make_spherical(1d0, PI / 2d0, 3d0 / 4d0 * PI)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_7:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_7: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_7

        function test_lib_tree_traversal_get_traversed_boxes_8() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            |     |     |     |     r
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x      ^
            !                    r: | ray
            !
            data%ray_origin = make_cartesian(1d0, 0.3d0, 0.15d0)
            data%ray_direction = make_cartesian(0d0, 1d0, 0d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (.not. allocated(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_8: OK"
            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_8: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_8

        function test_lib_tree_traversal_get_traversed_boxes_9() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            r     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x      ^
            !                    r: | ray
            !
            data%ray_origin = make_cartesian(0d0, 0.3d0, 0.15d0)
            data%ray_direction = make_cartesian(0d0, 1d0, 0d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (.not. allocated(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_9: OK"

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_9: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_9

        function test_lib_tree_traversal_get_traversed_boxes_10() result(rv)
            use lib_tree
            use lib_tree_type_operator
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            type(lib_tree_traversal_callback_type) :: callback

            type(lib_tree_traversal_data) :: data

            type(lib_tree_data_element), dimension(:), allocatable :: element_list

            type(lib_tree_universal_index), dimension(7) :: ground_truth

            ground_truth(1)%l = 0
            ground_truth(1)%n = 0

            ground_truth(2)%l = 2
            ground_truth(2)%n = 2

            ground_truth(3)%l = 2
            ground_truth(3)%n = 6

            ground_truth(4)%l = 1
            ground_truth(4)%n = 0

            ground_truth(5)%l = 2
            ground_truth(5)%n = 34

            ground_truth(6)%l = 2
            ground_truth(6)%n = 38

            ground_truth(7)%l = 1
            ground_truth(7)%n = 4

            allocate(element_list, source = lib_tree_get_diagonal_test_dataset(10_8, 1_1, 1_1))

            call lib_tree_constructor(element_list)

            callback = lib_tree_traversal_get_callback_functions()
            !            -------------
            !            |  2  |  6  |
            !            |     |     |
            !            ------0------
            !            |  0  |  4  |
            !         y ^|     |     |
            !           |-------------
            !           |--> x
            !
            !            -------------------------
            !            |  18 |  22 |  50 |  54 |
            !            |     |     |     |     |
            !            ------2-----------6------
            !            |  16 |  20 |  48 |  52 |
            !            |     |     |     |     |
            !            -------------------------
            !            |  2  |  6  |  34 |  38 |
            !            r     |     |     |     |
            !            ------0-----------4------
            !            |  0  |  4  |  32 |  36 |
            !         y ^|     |     |     |     |
            !           |-------------------------
            !           |--> x
            !                    r:  --> ray
            !
            data%ray_origin = make_cartesian(0d0, 0.3d0, 0.15d0)
            data%ray_direction = make_cartesian(1d0, 0d0, 0d0)

            call lib_tree_traversal_get_traversed_boxes(callback, data)

            rv = .true.

            if (size(ground_truth) .eq. size(data%traversed_box_list)) then
                print *, "test_lib_tree_traversal_get_traversed_boxes_10:"

                do i = lbound(ground_truth, 1), ubound(ground_truth, 1)
                    if (ground_truth(i) .eq. data%traversed_box_list(i)) then
                        print *, "  ", i, " : OK"
                    else
                        rv = .false.
                        print *, "  ", i, " : FAILED"
                        print *, "  box (l, n): ",  data%traversed_box_list(i)%l, data%traversed_box_list(i)%n
                        print *, "   gt (l, n): ",  ground_truth(i)%l, ground_truth(i)%n
                    end if
                end do

            else
                rv = .false.
                print *, "test_lib_tree_traversal_get_traversed_boxes_10: FAILED"
            end if

            call lib_tree_destructor()

        end function test_lib_tree_traversal_get_traversed_boxes_10

    end function lib_tree_traversal_test_functions

end module lib_tree_traversal
