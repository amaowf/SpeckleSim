module lib_math_triangle
    use lib_math_type
    use lib_math_type_operator
    use lib_math_constants
    implicit none

    private

    public :: set_normal
    public :: get_normal

    public :: lib_math_triangle_circumscribed_circle
    public :: lib_math_triangle_intersection_ray

    public :: lib_math_triangle_test_functions

    ! --- interfaces ---
    interface set_normal
        module procedure lib_math_triangle_set_normal
    end interface

    interface get_normal
        module procedure lib_math_plane_get_normal
    end interface

contains

    subroutine lib_math_triangle_set_normal(triangle)
        implicit none
        ! dummy
        type(triangle_type), intent(inout) :: triangle

        triangle%normal = lib_math_plane_get_normal(triangle%vertex)
    end subroutine lib_math_triangle_set_normal

    ! Argument
    ! ----
    !   vertex: cartesian_coordinate_real_type, dimension(3)
    !       three points on a plane
    ! Results
    ! ----
    !   rv: cartesian_coordinate_real_type
    !       plane normal, right-hand-rule
    function lib_math_plane_get_normal(vertex) result(rv)
        implicit none
        ! dummy
        type(cartesian_coordinate_real_type), dimension(3), intent(in) :: vertex
        type(cartesian_coordinate_real_type) :: rv


        rv = cross_product(vertex(2) - vertex(1),&
                           vertex(3) - vertex(1))
        rv = unit_vector(rv)
    end function lib_math_plane_get_normal

    ! Argument
    ! ----
    !   triangle: triangle_type
    !       triangle in tree-dimensional space
    !       HINT: the normal of the triangle must be given
    !
    ! Returns
    ! ----
    !   circumcentre: cartesian_coordinate_real_type
    !       centre of the circumscribed circle of the triangle
    !   circumradius: real
    !       radius of the circumscribed circle of the triangle
    subroutine lib_math_triangle_circumscribed_circle(triangle, circumcentre, circumradius)
        implicit none
        ! dummy
        type(triangle_type), intent(in) :: triangle

        type(cartesian_coordinate_real_type), intent(out) :: circumcentre
        real(kind = lib_math_type_kind), intent(out) :: circumradius

        ! auxiliary
        integer :: i
        type(spherical_coordinate_real_type) :: n

        type(cartesian_coordinate_real_type), dimension(3) :: vertex

        real(kind = lib_math_type_kind), dimension(3, 2) :: triangle_2D
        real(kind = lib_math_type_kind), dimension(2) :: circumcentre_2D

        type(cartresian_coordinate_rot_matrix_type) :: r

        ! coordiante transformation to the triangle plane
        n = triangle%normal
!        r = lib_math_get_matrix_rot(n%phi, n%theta, n%phi)
        r = lib_math_get_matrix_rot(-(n%phi - PI/2), -n%theta, 0d0)

        do i = lbound(triangle%vertex, 1), ubound(triangle%vertex, 1)
            vertex(i) = r * triangle%vertex(i)
        end do

        triangle_2D(1, :) = (/ vertex(1)%x, vertex(1)%y /)
        triangle_2D(2, :) = (/ vertex(2)%x, vertex(2)%y /)
        triangle_2D(3, :) = (/ vertex(3)%x, vertex(3)%y /)

        ! calculation of *circumcentre* and *circumradius*
        call lib_math_triangle_circumscribed_circle_2D(triangle_2D, circumcentre_2D, circumradius)

        circumcentre = make_cartesian(circumcentre_2D(1), circumcentre_2D(2), vertex(1)%z)

        r = lib_math_get_matrix_rot((n%phi - PI/2), n%theta, 0d0)
        circumcentre = r * circumcentre

    end subroutine lib_math_triangle_circumscribed_circle

    ! Argument
    ! ----
    !   triangle: real, dimension(3, 2)
    !       triangle in two-dimensional space
    !       HINT:
    !                  x1 y1
    !           array: x2 y2
    !                  x3 y3
    !
    ! Returns
    ! ----
    !   circumcentre: real(kind = lib_math_type_kind), dimension(2)
    !       centre of the circumscribed circle of the triangle
    !       HINT:
    !           array: x_0 y_0
    !   circumradius: real
    !       radius of the circumscribed circle of the triangle
    !
    ! Reference
    ! ----
    !   https://mathworld.wolfram.com/Circumcircle.html
    !   Equations: (4), (6), (7), (8), (11), (12), (13)
    subroutine lib_math_triangle_circumscribed_circle_2D(triangle, circumcentre, circumradius)
        implicit none
        ! dummy
        real(kind = lib_math_type_kind), dimension(3, 2), intent(in) :: triangle

        real(kind = lib_math_type_kind), dimension(2), intent(out) :: circumcentre
        real(kind = lib_math_type_kind), intent(out) :: circumradius

        ! auxiliary
        real(kind = lib_math_type_kind) :: a
        real(kind = lib_math_type_kind) :: a2
        real(kind = lib_math_type_kind) :: b_x
        real(kind = lib_math_type_kind) :: b_y
        real(kind = lib_math_type_kind) :: c

        real(kind = lib_math_type_kind) :: x_1_pw_2_y_1_pw_2
        real(kind = lib_math_type_kind) :: x_2_pw_2_y_2_pw_2
        real(kind = lib_math_type_kind) :: x_3_pw_2_y_3_pw_2

        real(kind = lib_math_type_kind) :: x_1_y_2
        real(kind = lib_math_type_kind) :: x_2_y_3
        real(kind = lib_math_type_kind) :: x_3_y_1

        real(kind = lib_math_type_kind) :: y_1_x_2
        real(kind = lib_math_type_kind) :: y_2_x_3
        real(kind = lib_math_type_kind) :: y_3_x_1

        x_1_y_2 = triangle(1, 1) * triangle(2, 2)
        x_2_y_3 = triangle(2, 1) * triangle(3, 2)
        x_3_y_1 = triangle(3, 1) * triangle(1, 2)

        y_1_x_2 = triangle(1, 2) * triangle(2, 1)
        y_2_x_3 = triangle(2, 2) * triangle(3, 1)
        y_3_x_1 = triangle(3, 2) * triangle(1, 1)

        x_1_pw_2_y_1_pw_2 = triangle(1, 1) * triangle(1, 1) + triangle(1, 2) * triangle(1, 2)
        x_2_pw_2_y_2_pw_2 = triangle(2, 1) * triangle(2, 1) + triangle(2, 2) * triangle(2, 2)
        x_3_pw_2_y_3_pw_2 = triangle(3, 1) * triangle(3, 1) + triangle(3, 2) * triangle(3, 2)

        ! eq. (4)
        a = x_1_y_2 + x_3_y_1 + x_2_y_3 &
            - y_1_x_2 - y_2_x_3 - y_3_x_1

        ! eq. (6)
        b_x = - (x_1_pw_2_y_1_pw_2 * triangle(2, 2) + x_2_pw_2_y_2_pw_2 * triangle(3, 2) + x_3_pw_2_y_3_pw_2 * triangle(1, 2) &
                 - x_3_pw_2_y_3_pw_2 * triangle(2, 2) - x_2_pw_2_y_2_pw_2 * triangle(1, 2) - x_1_pw_2_y_1_pw_2 * triangle(3, 2))

        ! eq. (7)
        b_y = x_1_pw_2_y_1_pw_2 * triangle(2, 1) + x_2_pw_2_y_2_pw_2 * triangle(3, 1) + x_3_pw_2_y_3_pw_2 * triangle(1, 1) &
              - x_3_pw_2_y_3_pw_2 * triangle(2, 1) - x_2_pw_2_y_2_pw_2 * triangle(1, 1) - x_1_pw_2_y_1_pw_2 * triangle(3, 1)

        ! eq. (8)
        c = - (x_1_pw_2_y_1_pw_2 * x_2_y_3 + x_2_pw_2_y_2_pw_2 * x_3_y_1 + x_3_pw_2_y_3_pw_2 * x_1_y_2 &
               - x_3_pw_2_y_3_pw_2 * y_1_x_2 - x_2_pw_2_y_2_pw_2 * y_3_x_1 - x_1_pw_2_y_1_pw_2 * y_2_x_3)

        ! eq. (11), (12)
        a2 = 2 * a
        circumcentre(1) = - b_x / a2
        circumcentre(2) = - b_y / a2

        ! eq. (13)
        circumradius = sqrt(b_x * b_x + b_y * b_y - 4 * a * c) / abs(a2)

    end subroutine lib_math_triangle_circumscribed_circle_2D

    ! Argument
    ! ----
    !   ray_origin: cartesian_coordinate_real_type
    !       origin of a ray
    !   ray_direction: cartesian_coordinate_real_type
    !       direction of a ray
    !   triangle: triangle_type
    !       triangle, which must be checked whether the ray and the triangle have an intersection.
    !
    ! Resturns
    ! ----
    !   intersection_point: cartesian_coordinate_real_type
    !       If the return value *rv* is true, the intersection of the ray and the triangle was calculated
    !   rv: logical
    !       true: there is an intersetion between the triangle and the ray
    !       false: there is no point of intersection
    !
    ! Reference: Fast, Minimum Storage RayTriangle Intersection, Moller and Trumbore
    !            and http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
    function lib_math_triangle_intersection_ray(ray_origin, &
        ray_direction, &
        triangle, &
        intersection_point) result(rv)
        implicit none
        ! dummy
        type(cartesian_coordinate_real_type), intent(in) :: ray_origin
        type(cartesian_coordinate_real_type), intent(in) :: ray_direction
        type(triangle_type), intent(in) :: triangle

        type(cartesian_coordinate_real_type), intent(out) :: intersection_point
        logical :: rv

        ! parameter
        real(lib_math_type_kind), parameter :: EPSILON = 0.0000001;

        ! auxiliary
        type(cartesian_coordinate_real_type) :: edge1, edge2, h, s, q
        real(kind=lib_math_type_kind) ::  det, inf_det, u, v, t

        edge1 = triangle%vertex(2) - triangle%vertex(1)
        edge2 = triangle%vertex(3) - triangle%vertex(1)
        h = cross_product(ray_direction, edge2)
        det = dot_product(edge1, h)
        if (det > -EPSILON .and. det < EPSILON) then
            rv =  .false.    ! This ray is parallel to this triangle.
            return
        end if
        inf_det = 1.0d0/det
        s = ray_origin - triangle%vertex(1)
        u = inf_det * dot_product(s, h)
        if (u .lt. 0.0d0 .or. u .gt. 1.0d0) then
            rv = .false.
            return
        end if
        q = cross_product(s, edge1)
        v = inf_det * dot_product(ray_direction, q)
        if (v .lt. 0.0d0 .or. (u + v) .gt. 1.0d0) then
            rv = .false.
            return
        end if
        ! At this stage we can compute t to find out where the intersection point is on the line.
        t = inf_det * dot_product(edge2, q)
        if (t .gt. EPSILON) then! ray intersection
            intersection_point = ray_origin + ray_direction * t
            rv = .true.
            return
        else ! This means that there is a line intersection but not a ray intersection.
            rv = .false.
            return
        end if
    end function lib_math_triangle_intersection_ray

    function lib_math_triangle_test_functions() result(rv)
        implicit none
        ! dummy
        integer :: rv

        real(kind = lib_math_type_kind), parameter :: ground_truth_error = 1d-14

        rv = 0
        if (.not. test_lib_math_triangle_circumscribed_circle_2D()) rv = rv + 1
        if (.not. test_lib_math_triangle_circumscribed_circle()) rv = rv + 1
        if (.not. test_lib_math_triangle_intersection_ray()) rv = rv + 1

    contains

        function test_lib_math_triangle_circumscribed_circle_2D() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer :: i
            real(kind = lib_math_type_kind), dimension(3, 2) :: triangle
            real(kind = lib_math_type_kind), dimension(2) :: circumcentre
            real(kind = lib_math_type_kind) :: circumradius

            real(kind = lib_math_type_kind), dimension(2) :: gt_circumcentre
            real(kind = lib_math_type_kind) :: gt_circumradius

            triangle(1, :) = (/ -1, 0 /) + 1
            triangle(2, :) = (/ 0, 1 /) + 1
            triangle(3, :) = (/ 1, 0 /) + 1

            gt_circumcentre(:) = (/ 1, 1 /)
            gt_circumradius = 1

            call lib_math_triangle_circumscribed_circle_2D(triangle, circumcentre, circumradius)

            rv = .true.

            print *, "test_lib_math_triangle_circumscribed_circle_2D"
            if (circumradius .eq. gt_circumradius) then
                print *, "  radius: OK"
            else
                print *, "  radius: FAILED"
                print *, "  radius: ", circumradius
                print *, "      gt: ", gt_circumradius
                rv = .false.
            end if

            do i = lbound(gt_circumcentre, 1), ubound(gt_circumcentre, 1)
                if (circumcentre(i) .eq. gt_circumcentre(i)) then
                    print *, "  coordinate component ", i, " : OK"
                else
                    print *, "  coordinate component ", i, " : FAILED"
                    print *, "  value: ", circumcentre(i)
                    print *, "     gt: ", gt_circumcentre(i)
                end if
            end do

        end function test_lib_math_triangle_circumscribed_circle_2D

        function test_lib_math_triangle_circumscribed_circle() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(triangle_type) :: triangle
            type(cartesian_coordinate_real_type) :: circumcentre
            real(kind = lib_math_type_kind) :: circumradius
            type(cartesian_coordinate_real_type) :: offset

            type(cartesian_coordinate_real_type) :: ground_truth_circumcentre
            real(kind = lib_math_type_kind) :: ground_truth_circumradius

            real(kind = lib_math_type_kind) :: buffer

            offset = 10d0 * make_cartesian(1d0, 1d0, 1d0)

            triangle%vertex(1) = make_cartesian(-1d0, 0d0, 0d0) + offset
            triangle%vertex(2) = make_cartesian( 0d0, 0d0, 1d0) + offset
            triangle%vertex(3) = make_cartesian( 1d0, 0d0, 0d0) + offset

            call set_normal(triangle)

            ground_truth_circumcentre = offset !make_cartesian(1d0, 1d0, 1d0)
            ground_truth_circumradius = 1d0

            call lib_math_triangle_circumscribed_circle(triangle, circumcentre, circumradius)

            rv = .true.

            print *, "test_lib_math_triangle_circumscribed_circle"

            buffer = ground_truth_circumradius - circumradius
            if (abs(buffer) .lt. ground_truth_error) then
                print *, "  radius: OK"
            else
                print *, "  radius: FAILED"
                print *, "  radius: ", circumradius
                print *, "      gt: ", ground_truth_circumradius
                rv = .false.
            end if

            buffer = ground_truth_circumcentre%x - circumcentre%x
            if (abs(buffer) .lt. ground_truth_error) then
                print *, "  coordinate component x : OK"
            else
                print *, "  coordinate component x : FAILED"
                print *, "  value: ", circumcentre%x
                print *, "     gt: ", ground_truth_circumcentre%x
                rv = .false.
            end if

            buffer = ground_truth_circumcentre%y - circumcentre%y
            if (abs(buffer) .lt. ground_truth_error) then
                print *, "  coordinate component y : OK"
            else
                print *, "  coordinate component y : FAILED"
                print *, "  value: ", circumcentre%y
                print *, "     gt: ", ground_truth_circumcentre%y
                rv = .false.
            end if

            buffer = ground_truth_circumcentre%z - circumcentre%z
            if (abs(buffer) .lt. ground_truth_error) then
                print *, "  coordinate component z : OK"
            else
                print *, "  coordinate component z : FAILED"
                print *, "  value: ", circumcentre%z
                print *, "     gt: ", ground_truth_circumcentre%z
                rv = .false.
            end if

        end function test_lib_math_triangle_circumscribed_circle

        function test_lib_math_triangle_intersection_ray() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(cartesian_coordinate_real_type):: ray_origin
            type(cartesian_coordinate_real_type):: ray_direction
            type(triangle_type) :: triangle

            type(cartesian_coordinate_real_type) :: intersection_point
            type(cartesian_coordinate_real_type) :: intersection_point_ground_truth

            logical :: res
            logical :: res_ground_truth

            ray_origin = make_cartesian(-1d0, 0d0, 0d0)
            ray_direction = make_spherical(1d0, PI / 2d0, 0d0)

            triangle%vertex(1) = make_cartesian(1d0, 0d0, 1d0)
            triangle%vertex(2) = make_cartesian(1d0, 1d0, -1d0)
            triangle%vertex(3) = make_cartesian(1d0, -1d0, -1d0)

            res_ground_truth = .true.
            intersection_point_ground_truth = make_cartesian(1d0, 0d0, 0d0)

            res = lib_math_triangle_intersection_ray(ray_origin, ray_direction, triangle, intersection_point)


            if (res .eqv. res_ground_truth) then
                rv = .true.
                print *, "test_lib_math_triangle_intersection_ray: OK"
            else
                rv = .false.
                print *, "test_lib_math_triangle_intersection_ray: ERROR"
                print *, "  intersection point: ", intersection_point%x, intersection_point%y, intersection_point%z
                print *, "        ground truth: ", intersection_point_ground_truth%x, intersection_point_ground_truth%y, &
                    intersection_point_ground_truth%z
            end if

        end function test_lib_math_triangle_intersection_ray

    end function lib_math_triangle_test_functions

end module lib_math_triangle
