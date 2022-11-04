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

module lib_data_types
    implicit none

    integer, parameter :: element_kind = 4

    ! list element
    type lib_tree_element
        integer(kind=element_kind)   :: start
        integer(kind=element_kind)   :: end
        integer(kind=element_kind)   :: parent_element
    end type lib_tree_element

!    type tree_parent
!        type(tree_element), dimension(8) :: children
!        type(ectree_element), dimension() :: neighbor
!    end type tree_level



end module lib_data_types
