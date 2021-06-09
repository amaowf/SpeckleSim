!    Copyright (C) 2021  Liwei Fu <liwei.fu@ito.uni-stuttgart.de>
!
!    This file is part of SpeckleSim.
!
!    SpeckleSim is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SpeckleSim is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.!
! 
! @author: Liwei Fu
	
	module lib_sie_constants
   
	implicit none 
	
	!dp : double precision
	!It works only for double precision so that it is compatible with different libraries
	integer, parameter, public :: dp = 8 
	
	real(dp), parameter, public :: c0 = 2.99792458e+8
	real(dp), parameter, public :: eps_0 = 8.854187817e-12
	real(dp), parameter, public :: my_0 = 1.2566370614e-6
	complex(dp), parameter :: im = (0.0, 1.0)
		
end module lib_sie_constants
