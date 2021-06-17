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

module lib_scene_generator
    use libmath
    use lib_scene_type
    implicit none

    private

    ! public functions
    public :: lib_scene_gerator_hcp_lattice
    public :: lib_scene_generator_hcp_lattice_fill_cuboid
    public :: lib_scene_generator_hcp_lattice_fill_sphere

    public :: lib_scene_generator_test_functions


    contains

        ! Argument
        ! ----
        !   i, j, k: integer
        !       sphere index starting at 0 for the x-, y- and z-coordinates
        !   sphere_radius: doubel precision
        !       radius of one sphere [m]
        !
        ! Returns
        ! ----
        !   rv: type(cartesian_coordinate_real_type)
        !       coordinate point of the lattice for a hexagonal close packing
        !
        ! Reference: http://mathworld.wolfram.com/HexagonalClosePacking.html
        !            https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres <-- todo: replace
        !
        function lib_scene_gerator_hcp_lattice(i, j, k, sphere_radius) result(rv)
            implicit none
            ! dummy
            integer, intent(in) :: i
            integer, intent(in) :: j
            integer, intent(in) :: k
            double precision, intent(in) :: sphere_radius

            type(cartesian_coordinate_real_type) :: rv

            rv%x = dble(2 * i + mod(dble(j + k), 2d0))
            rv%y = sqrt(3d0) * (dble(j) + mod(dble(k), 2d0) / 3d0)
            rv%z = 2d0 * sqrt(6d0) * dble(k) / 3d0

            rv = sphere_radius * rv

        end function lib_scene_gerator_hcp_lattice

        ! Argument
        ! ----
        !   size_x, size_y, size_z: doubel precision
        !       edge length of the cuboid with a hexagonal close packing scheme
        !   sphere_radius: double precision
        !       radius of the spheres at the lattice coordinate points
        !
        ! Returns
        ! ----
        !   rv: type(lib_scene_object_hcp_cuboid)
        !       an arangement of sphers with a hexagonal cloase packing (hcp) in a cuboid
        !
        !   <---      size_x     --->
        !    . . . . . . . . . . . .  ^
        !   . . . . . . . . . . . . . |
        !    . . . . . z . . . . . .  |
        !   . . . . . . ^ . . . . . . |
        !    . . . .K_o |. . . . . .  |
        !   . . . . . . --> x . . . . | size_z
        !    . . . . . . . . . . . .  |
        !   . . . . . . . . . . . . . |
        !    . . . . . . . . . . . .  |
        !   . . . . . . . . . . . . . |
        !    . . . . . . . . . . . .  v
        !
        !
        !   K_o: object coordinate system
        !   .: spheres at the hcp_lattice_coordiantes
        !
        function lib_scene_generator_hcp_lattice_fill_cuboid(d_o_j, size_x, size_y, size_z, sphere_radius) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type) :: d_o_j
            double precision, intent(in) :: size_x
            double precision, intent(in) :: size_y
            double precision, intent(in) :: size_z
            double precision, intent(in) :: sphere_radius

            type(lib_scene_object_hcp_cuboid_type) :: rv

            ! auxiliary
            integer :: x
            integer :: i
            integer :: j
            integer :: k

            integer :: max_i
            integer :: max_j
            integer :: max_k

            double precision :: buffer

            type(cartesian_coordinate_real_type) :: offset

            rv%d_o_j = d_o_j

            buffer = (size_z - 2d0 * sphere_radius) / sphere_radius
            buffer = buffer * 3d0 / ( 2d0 * sqrt(6d0) )
            max_k = max(int(floor(  buffer  )), 0)

            buffer = (size_y - 2d0 * sphere_radius) / sphere_radius
            buffer = buffer / sqrt(3d0)
            buffer = buffer - mod(dble(max_k), 2d0) / 3d0
            max_j = max(int(floor(  buffer  )), 0)

            buffer = (size_x - 2d0 * sphere_radius) / sphere_radius
            buffer = buffer - mod(dble(max_j + max_k), 2d0)
            buffer = buffer / 2d0
            max_i = max(int(floor(  buffer  )), 0)

            offset = make_cartesian(0d0, 0d0, 0d0)
            offset = offset - lib_scene_gerator_hcp_lattice(max_i, max_j, max_k, sphere_radius) / 2d0

            allocate(rv%hcp_lattice_coordiantes((max_i+1) * (max_j+1) * (max_k+1)))

            x = 0
            do i = 0, max_i
                do j = 0, max_j
                    do k = 0, max_k
                        x = x + 1
                        rv%hcp_lattice_coordiantes(x) = lib_scene_gerator_hcp_lattice(i, j, k, sphere_radius) &
                                                        + offset
                    end do
                end do
            end do

        end function lib_scene_generator_hcp_lattice_fill_cuboid


        ! Argument
        ! ----
        !   size_x, size_y, size_z: doubel precision
        !       edge length of the cuboid with a hexagonal close packing scheme
        !   sphere_radius: double precision
        !       radius of the sphere containing the lattice spheres (long distance order)
        !   lattice_sphere_radius: double precision
        !       radius of the spheres at the lattice coordinate points
        !   spherical_cap_point: type(cartesian_coordinate_real_type)
        !       point on the plane intersecting the sphere
        !   spherical_cap_normal: type(cartesian_coordinate_real_type)
        !       plane normal vector, in this direction the lattice spheres are removed
        !
        !
        ! Returns
        ! ----
        !   rv: type(lib_scene_object_hcp_cuboid)
        !       an arangement of sphers with a hexagonal cloase packing (hcp) in a cuboid
        !
        !   <---      size_x     --->
        !    x x x x . . . . x x x x  ^
        !    x x  . . . . . . . x x x |
        !    x . . . . z . . . . . x  |
        !   x . . . . . ^ . . . . . x |
        !    . . . .K_o |. . . . . .  |
        !   . . . . . . --> x . . . . | size_z
        !    . . . . . . . . . . . .  |
        !   x . . . . . . . . . . . x |
        !    x . . . . . . . . . . x  |
        !   x x x . . . . . . . x x x |
        !    x x x x . . . . x x x x  v
        !
        !   K_o: object coordinate system
        !   .: spheres at the hcp_lattice_coordiantes inside the sphere with the radius "sphere_radius"
        !   x: spheres at the hcp_lattice_coordiantes outside the sphere with the radius "sphere_radius"
        !
        !
        ! with spherical cap
        !
        !   <---      size_x     --->
        !    x x x x . . \ x x x x x  ^
        !    x x  . . . . \ x x x x x |
        !    x . . . . z . \ x x x x  |
        !   x . . . . . ^ . \ x x x x |
        !    . . . .K_o |. . \ x x x  |
        !   . . . . . . --> x \ x x x | size_z
        !    . . . . . . . . . \ x x  |
        !   x . . . . . . . . . \ x x |
        !    x . . . . . . . . . \ x  |
        !   x x x . . . . . . . x \ x |
        !    x x x x . . . . x x x \  v
        !
        !
        function lib_scene_generator_hcp_lattice_fill_sphere(d_o_j, sphere_radius, lattice_sphere_radius, &
                                                             cap_plane_point, cap_plane_normal) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type) :: d_o_j
            double precision, intent(in) :: sphere_radius
            double precision, intent(in) :: lattice_sphere_radius

            type(cartesian_coordinate_real_type), intent(in), optional :: cap_plane_point
            type(cartesian_coordinate_real_type), intent(in), optional :: cap_plane_normal

            type(lib_scene_object_hcp_sphere_type) :: rv

            ! auxiliary
            integer :: i

            integer, dimension(2):: i_range

            type(spherical_coordinate_real_type) :: spherical_coord


            logical :: use_spherical_cap
            type(cartesian_coordinate_real_type) :: m_cap_plane_point
            type(cartesian_coordinate_real_type) :: m_cap_plane_normal
            double precision :: distance

            rv%d_o_j = d_o_j

            if (present(cap_plane_point)) m_cap_plane_point = cap_plane_point
            if (present(cap_plane_normal)) m_cap_plane_normal = cap_plane_normal / abs(cap_plane_normal)

            use_spherical_cap = .false.
            if (present(cap_plane_point) .and. present(cap_plane_normal)) use_spherical_cap = .true.

            rv%hcp_cuboid = lib_scene_generator_hcp_lattice_fill_cuboid(d_o_j, &
                                                                        2d0 * sphere_radius, &
                                                                        2d0 * sphere_radius, &
                                                                        2d0 * sphere_radius, &
                                                                        lattice_sphere_radius)

            i_range(1) = lbound(rv%hcp_cuboid%hcp_lattice_coordiantes, 1)
            i_range(2) = ubound(rv%hcp_cuboid%hcp_lattice_coordiantes, 1)


            allocate (rv%inside_sphere(i_range(1):i_range(2)))

            !$OMP PARALLEL DO PRIVATE(i, spherical_coord, distance)
            do i = i_range(1), i_range(2)
                spherical_coord = rv%hcp_cuboid%hcp_lattice_coordiantes(i)
                if ( spherical_coord%rho .lt. sphere_radius - lattice_sphere_radius / 2d0) then
                    if (use_spherical_cap) then
                        distance = dot_product(rv%hcp_cuboid%hcp_lattice_coordiantes(i) - cap_plane_point, &
                                               cap_plane_normal)
                        if (distance .lt. 0d0) then
                            rv%inside_sphere(i) = .true.
                        else
                            rv%inside_sphere(i) = .false.
                        end if
                    else
                        rv%inside_sphere(i) = .true.
                    end if
                else
                    rv%inside_sphere(i) = .false.
                end if
            end do
            !$OMP END PARALLEL DO

        end function lib_scene_generator_hcp_lattice_fill_sphere

        function lib_scene_generator_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_scene_generator_hcp_lattice_fill_cuboid()) rv = rv + 1
            if (.not. test_lib_scene_generator_hcp_lattice_fill_sphere()) rv = rv + 1

            print *, "-------------test_lib_scene_generator----------------"
            if (rv == 0) then
                print *, "test_lib_scene_generator tests: OK"
            else
                print *, rv,"test_lib_scene_generator test(s) FAILED"
            end if
            print *, "------------------------------------------------------"

            contains

                function test_lib_scene_generator_hcp_lattice_fill_cuboid() result(rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    integer, dimension(2) :: range_i

                    double precision, dimension(:,:), allocatable :: data_list
                    type(cartesian_coordinate_real_type), dimension(:), allocatable :: coordinates

                    type(cartesian_coordinate_real_type) :: d_o_j

                    double precision :: size_x
                    double precision :: size_y
                    double precision :: size_z

                    double precision :: sphere_radius

                    type(lib_scene_object_hcp_cuboid_type) :: hcp_cuboid

                    integer :: u
                    character(len=50) :: filename

                    size_x = 8
                    size_y = 8
                    size_z = 8

                    sphere_radius = 0.5
                    d_o_j = make_cartesian(0d0, 0d0, 0d0)
                    hcp_cuboid = lib_scene_generator_hcp_lattice_fill_cuboid(d_o_j, size_x, size_y, size_z, sphere_radius)

                    range_i(1) = lbound(hcp_cuboid%hcp_lattice_coordiantes, 1)
                    range_i(2) = ubound(hcp_cuboid%hcp_lattice_coordiantes, 1)

                    allocate(data_list(range_i(1):range_i(2), 3))

                    do i = range_i(1), range_i(2)
                        data_list(i, 1) = hcp_cuboid%hcp_lattice_coordiantes(i)%x
                        data_list(i, 2) = hcp_cuboid%hcp_lattice_coordiantes(i)%y
                        data_list(i, 3) = hcp_cuboid%hcp_lattice_coordiantes(i)%z
                    end do

                    u = 99
                    filename = "temp/hcp_cuboid.csv"
                    open(unit=u, file=trim(filename), status='unknown')
                    rv = write_csv(u, data_list)
                    close(u)

                end function test_lib_scene_generator_hcp_lattice_fill_cuboid

                function test_lib_scene_generator_hcp_lattice_fill_sphere() result(rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    integer, dimension(2) :: range_i

                    double precision, dimension(:,:), allocatable :: data_list
                    type(cartesian_coordinate_real_type), dimension(:), allocatable :: coordinates

                    type(cartesian_coordinate_real_type) :: d_o_j

                    double precision :: sphere_radius
                    double precision :: lattice_sphere_radius

                    type(cartesian_coordinate_real_type) :: cap_plane_point
                    type(cartesian_coordinate_real_type) :: cap_plane_normal

                    type(lib_scene_object_hcp_sphere_type) :: hcp_sphere

                    integer :: u
                    character(len=50) :: filename

                    sphere_radius = 0.5 * unit_mu
                    lattice_sphere_radius = 15.5 * unit_nm
                    d_o_j = make_cartesian(0d0, 0d0, 0d0)

                    cap_plane_point = make_cartesian(0d0, 0d0, 0d0)
                    cap_plane_normal = make_cartesian(0d0, 0d0, -1d0)

                    hcp_sphere = lib_scene_generator_hcp_lattice_fill_sphere(d_o_j, sphere_radius, lattice_sphere_radius, &
                                                                             cap_plane_point, cap_plane_normal)

                    coordinates = list_filter(hcp_sphere%hcp_cuboid%hcp_lattice_coordiantes, &
                                              hcp_sphere%inside_sphere)

                    range_i(1) = lbound(coordinates, 1)
                    range_i(2) = ubound(coordinates, 1)

                    allocate(data_list(range_i(1):range_i(2), 3))

                    do i = range_i(1), range_i(2)
                        data_list(i, 1) = coordinates(i)%x
                        data_list(i, 2) = coordinates(i)%y
                        data_list(i, 3) = coordinates(i)%z
                    end do

                    u = 99
                    filename = "temp/hcp_sphere.csv"
                    open(unit=u, file=trim(filename), status='unknown')
                    rv = write_csv(u, data_list)
                    close(u)

                end function test_lib_scene_generator_hcp_lattice_fill_sphere

        end function lib_scene_generator_test_functions
end module lib_scene_generator
