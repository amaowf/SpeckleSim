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

module lib_mie_type_functions
    use libmath
    use lib_mie_type
    use lib_mie_ss_helper_functions
    implicit none

    contains

        ! Argument
        ! ----
        !   lambda_0: double precision
        !       vacuum wave length
        !   n_medium: double precision
        !       refractive index of the medium
        !   r_particle: double precision
        !       radius of the particel
        !   n_particle: double complex
        !       refractive index of the particle
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /)
        !
        ! Returns
        ! ----
        !   sphere_parameter: type(lib_mie_sshere_parameter_type)
        function lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                              r_particle, n_particle,&
                                              n_range) &
                                            result (sphere_parameter)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: r_particle
            double complex, intent(in) :: n_particle
            integer, dimension(2) :: n_range

            type(lib_mie_sphere_parameter_type) :: sphere_parameter

            ! auxiliary
            double precision :: size_parameter

            size_parameter = 2 * PI * n_medium * r_particle / lambda_0

            sphere_parameter%n_range = n_range
            sphere_parameter%size_parameter = size_parameter
            sphere_parameter%radius = r_particle
            sphere_parameter%refractive_index = n_particle

            allocate(sphere_parameter%a_n%item(n_range(1):n_range(2)))
            allocate(sphere_parameter%b_n%item(n_range(1):n_range(2)))

            if (aimag(n_particle) .eq. 0d0) then
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, real(n_particle)/n_medium, n_range, &
                                                       sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            else
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, n_particle/n_medium, n_range, &
                                                           sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            end if
        end function lib_mie_type_func_get_sphere_parameter

        ! Skatch
        ! ----
        !             _________
        !             ___k^____
        !             ____|____
        !             _________ plane wave
        !                 z
        !                 ^
        !             K_i |
        !                 --> x
        !                ^
        !               /
        !           z  /d_0_i
        !           ^ /
        !       K_0 |/
        !           --> x
        !
        ! K_0: world coordinate system
        ! K_i: illumination coordinate system
        !
        !
        ! Argument
        ! ----
        !   type: integer
        !       illumination type
        !       - 1: plane wave
        !   lambda_0: double precision
        !       vacuum wave length
        !   e_field_0: double precision
        !       magnitude of the oscilating electrical field (peak value)
        !   k: type(cartesian_coordinate_real_type)
        !       sets only the direction of the wave vector, the length of this vector is calculated internally as follows:
        !       |k| = 2 Pi / lambda
        !   d_0_i: type(cartesian_coordinate_real_type)
        !       position of the illumination coordinate system respect to the world coordinate system
        !
        ! Returns
        ! ----
        !   illumination: lib_mie_illumination_parameter
        !
        function lib_mie_type_func_get_plane_wave_illumination(lambda_0, refractive_index_medium, e_field_0, &
                                                               g, k, d_0_i) result(illumination)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda_0
            double precision, intent(in) :: refractive_index_medium
            double precision, intent(in) :: e_field_0
            double precision, dimension(:), intent(in) :: g
            type(cartesian_coordinate_real_type), dimension(lbound(g, 1):ubound(g, 1)), intent(in) :: k
            type(cartesian_coordinate_real_type), dimension(lbound(g, 1):ubound(g, 1)), intent(in) :: d_0_i


            type(lib_mie_illumination_parameter) :: illumination

            ! auxiliary
            integer :: i
            type(spherical_coordinate_real_type) :: k_spherical

            illumination%lambda_0 = lambda_0
            illumination%e_field_0 = e_field_0

            allocate(illumination%plane_wave(lbound(g, 1):ubound(g, 1)))

            do i = lbound(g, 1), ubound(g, 1)
                k_spherical = k(i)

                if (k_spherical%rho .gt. 0) then
                    illumination%plane_wave(i)%g = g(i)
                    illumination%plane_wave(i)%d_0_i = d_0_i(i)
                    illumination%plane_wave(i)%beam_parameter%wave_length_0 = lambda_0
                    illumination%plane_wave(i)%beam_parameter%theta = k_spherical%theta
                    illumination%plane_wave(i)%beam_parameter%phi = k_spherical%phi

                    illumination%plane_wave(i)%beam_parameter%e_field_0 = e_field_0
                    illumination%plane_wave(i)%beam_parameter%polarisation%x = 1
                    illumination%plane_wave(i)%beam_parameter%polarisation%y = 0
                    illumination%plane_wave(i)%beam_parameter%convention = 2

                    illumination%plane_wave(i)%beam_parameter%refractive_index_medium = refractive_index_medium
                else
                    print *, "lib_mie_type_func_get_illumination: ERROR"
                    print *, "  abs(k(i=", i, ")) = 0"
                end if
            end do

        end function lib_mie_type_func_get_plane_wave_illumination

end module lib_mie_type_functions
