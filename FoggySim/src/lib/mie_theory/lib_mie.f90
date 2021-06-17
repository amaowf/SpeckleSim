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

module lib_mie
    use libmath
    use lib_constants
    implicit none

    private

    public :: lib_mie_get_field

    contains

        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle [radian]
        !   phi: double precision
        !       azimuthal angle [radian]
        !   r: double precision
        !       distance
        !   e_field_0: double precision
        !       amplitude of the incident wave
        !   lambda: double precision
        !       wave length
        !   n_medium: double precision
        !       refractive index of the medium
        !   r_particle: double precision
        !       radius of the sphere
        !   n_particle: double precision
        !       refractive index of the sphere
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   alpha: double precision, optional(std: 0)
        !       incident angle with respect to the z-axis
        !       codomain: 0..Pi
        !   beta: double precision, optional(std: 0)
        !       angle between the x axis and the projection of the wave vector on the x-y plane
        !       codomain: 0..2Pi
        !
        ! Results
        ! ----
        !   field_s: type(spherical_coordinate_cmplx_type), dimension(2)
        !       1: scattered e field
        !       2: scattered h field
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 17
        !
        function lib_mie_get_field(theta, phi, r, &
                                    e_field_0, wave_length_0, &
                                    n_medium, &
                                    p0_nm, q0_nm, z_selector) &
                                result (field_0)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            type(list_list_cmplx),intent(in) :: p0_nm
            type(list_list_cmplx),intent(in) :: q0_nm
            integer(kind=1), intent(in), optional :: z_selector

            type(spherical_coordinate_cmplx_type), dimension(2) :: field_0

            ! auxiliary
            integer(kind=4) :: i
            type(spherical_coordinate_cmplx_type) :: e_field_incident_0
            type(spherical_coordinate_cmplx_type) :: h_field_incident_0

            double precision :: k0 ! wave number = 2 pi / lambda

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            integer(kind=4) :: m
            integer(kind=4) :: n

            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: m_z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(size(p0_nm%item, 1)) :: e_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_e_field_n_incident_0
            type(spherical_coordinate_cmplx_type), dimension(size(p0_nm%item, 1)) :: h_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_h_field_n_incident_0

            type(list_list_logical) :: calc_order_m

            ! --- init ---
            ! standard values
            m_z_selector = 1
            if (present(z_selector)) m_z_selector = int(z_selector, 1)

            k0 = 2.0_8 * PI / wave_length_0

            n_range = (/ lbound(p0_nm%item, 1), ubound(p0_nm%item, 1) /)

            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1, .true.)

            call lib_math_vector_spherical_harmonics_components(theta, phi, r, &
                                                                k0 * n_medium, &
                                                                n_range, m_z_selector, &
                                                                M_nm, N_nm)
            ! eq. (5)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real)
            do n=n_range(1), n_range(2)
                do m=-n, n
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_real = abs(e_field_0) * real((2*n+1), kind=8) &
                                      * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

                        if (isnan(buffer_real)) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
                            print *, "  buffer_real is NaN"
                            print * , "  n = ", n
                            print * , "  m = ", m
                        end if
                        e_field_nm%item(n)%item(m) = buffer_real * cmplx(0,1, kind=8)**n

                        if(buffer_real .eq. 0.0) then
                            calc_order_m%item(n)%item(m) = .false.
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO

            ! eq. (17)
            !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx, buffer_e_field_n_incident_0, buffer_h_field_n_incident_0)
            do n= n_range(1), n_range(2)
                i = n - n_range(1) + 1

                e_field_n_incident_0(i)%rho = cmplx(0,0)
                e_field_n_incident_0(i)%theta = cmplx(0,0)
                e_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_e_field_n_incident_0%rho = cmplx(0,0)
                buffer_e_field_n_incident_0%theta = cmplx(0,0)
                buffer_e_field_n_incident_0%phi = cmplx(0,0)

                h_field_n_incident_0(i)%rho = cmplx(0,0)
                h_field_n_incident_0(i)%theta = cmplx(0,0)
                h_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_h_field_n_incident_0%rho = cmplx(0,0)
                buffer_h_field_n_incident_0%theta = cmplx(0,0)
                buffer_h_field_n_incident_0%phi = cmplx(0,0)

                do m=-n, n
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                        buffer_e_field_n_incident_0 = buffer_cmplx &
                                                      * (N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                         +M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                        buffer_h_field_n_incident_0 = e_field_nm%item(n)%item(m) &
                                             * (M_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                +N_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                        if (isnan(real(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%theta)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%theta)) ) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
                            print *, "  buffer_e_field_n_incident_0 is NaN"
                            print * , "  n = ", n
                            print * , "  m = ", m
                        end if

                        e_field_n_incident_0(i) = e_field_n_incident_0(i) + buffer_e_field_n_incident_0
                        h_field_n_incident_0(i) = h_field_n_incident_0(i) + buffer_h_field_n_incident_0
                    end if

                end do
            end do
            !$OMP END PARALLEL DO

            e_field_incident_0%theta = cmplx(0,0,kind=8)
            e_field_incident_0%phi = cmplx(0,0,kind=8)
            e_field_incident_0%rho = cmplx(0,0,kind=8)

            h_field_incident_0%theta = cmplx(0,0,kind=8)
            h_field_incident_0%phi = cmplx(0,0,kind=8)
            h_field_incident_0%rho = cmplx(0,0,kind=8)

            do i=n_range(2)-n_range(1)+1, 1, -1
                e_field_incident_0 = e_field_incident_0 + e_field_n_incident_0(i)
                h_field_incident_0 = h_field_incident_0 + h_field_n_incident_0(i)
            end do

            field_0(1) = (-1D0) * e_field_incident_0

            ! omega = k * v
            ! v = c0 / n_medium
            ! mu = 1 <-- definition
            !
            !   k / (omega * mu) = k / ( k * v * 1 )
            ! = k / ( k * c0 / n_medium )
            ! = n_medium / c0
            field_0(2) = (-1D0) * h_field_incident_0 * n_medium / real(const_c0, kind=8)

        end function lib_mie_get_field
end module lib_mie
