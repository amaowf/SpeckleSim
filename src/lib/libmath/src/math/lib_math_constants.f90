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

module lib_math_constants
    implicit none

    double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

    double precision, parameter :: unit_m =  1.d0
    double precision, parameter :: unit_mm = 1.d-3
    double precision, parameter :: unit_mu = 1.d-6
    double precision, parameter :: unit_nm = 1.d-9

end module lib_math_constants
