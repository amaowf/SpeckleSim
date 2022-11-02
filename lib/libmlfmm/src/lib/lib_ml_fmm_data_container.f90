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

module lib_ml_fmm_data_container
    use  ml_fmm_type
    implicit none

    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the X-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the Y-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_v
end module lib_ml_fmm_data_container
