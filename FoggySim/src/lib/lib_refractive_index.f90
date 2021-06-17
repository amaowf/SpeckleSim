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

module lib_refractive_index
    use toolbox
    use file_io
    implicit none

!    integer, dimension()
!
!    contains
!
!        function lib_refractive_index(medium, wave_length) result (rv)
!            implicit none
!            ! dummy
!            integer, intent(in) :: medium
!            double precision, intent(in) :: wave_length
!
!            double complex :: rv
!
!            ! auxiliary
!
!            call data_interpolation(data_refractive_index_particle, 1, lambda / unit_mu, &
!                                                    data_refractive_index_particle_interpolation)
!            rv = dcmplx(data_refractive_index_particle_interpolation(2), &
!                                data_refractive_index_particle_interpolation(3))
!
!
!        end function
end module lib_refractive_index
