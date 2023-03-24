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

module lib_constants
    implicit none

    ! --- parameter ---

    ! speed of light in vacuum
    ! https://physics.nist.gov/cgi-bin/cuu/Value?c
    integer, parameter :: const_c0 =  299792458_4 ! m/s

    ! vacuum electric permittivity
    ! https://physics.nist.gov/cgi-bin/cuu/Value?ep0|search_for=permitivity
    double precision, parameter :: const_epsilon_0 = 8.8541878128d-12! [F m-1]

    ! characteristic impedance of vacuum
    ! https://physics.nist.gov/cgi-bin/cuu/Value?z0|search_for=impedance
    double precision, parameter :: const_z_0 = 376.730313668d0 ! [Ohm]

end module lib_constants
