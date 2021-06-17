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

module lib_mie_ms_data_container
    use lib_mie_type
    implicit none

    type(lib_mie_simulation_parameter_type) :: simulation_data

    contains

        subroutine lib_mie_ms_data_container_destructor()
            implicit none

            if (allocated(simulation_data%sphere_list)) deallocate(simulation_data%sphere_list)
            if (allocated(simulation_data%sphere_parameter_list)) deallocate(simulation_data%sphere_parameter_list)
            if (allocated(simulation_data%illumination%plane_wave)) deallocate(simulation_data%illumination%plane_wave)
        end subroutine

end module lib_mie_ms_data_container
