module lib_sie_constants
   
	implicit none 
	
	! dp : precision
	!It works only for double precision, in order to coordinate with different libraries
	integer, parameter, public :: dp = 8 
	
	real(dp), parameter, public :: c0 = 2.99792458e+8
	real(dp), parameter, public :: eps_0 = 8.854187817e-12
	real(dp), parameter, public :: my_0 = 1.2566370614e-6
	complex(dp), parameter :: im = (0.0, 1.0)
		
end module lib_sie_constants
