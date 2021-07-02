module input_tri
  
    
    use libmath
	
	implicit none  
	integer, parameter, public :: dp = 8
    !real(kind = dp), parameter, public :: PI = 3.141592653589793238462	
    real(kind = dp), parameter, public :: c0 = 2.99792458e+8
    real(kind = dp), parameter, public :: eps_0 = 8.854187817e-12
    real(kind = dp), parameter, public :: my_0 = 1.2566370614e-6
    complex(kind = dp), parameter :: im = (0.0, 1.0)	
!**************************************************************************
    ! Input parameters    
	real(kind = dp), parameter, public :: Lambda = 600e-9 !meter, 600nm, 500THz     	
	real(kind = dp), parameter :: theta_in = 0*PI/180 !incident angle for Gaussian beam
	real(kind = dp), parameter :: phi_in = 0.0 !incident angle for Gaussian beam
	complex(kind = dp), dimension(3), public :: E_in_hat = (/cos(theta_in)*cos(phi_in), cos(theta_in)*sin(phi_in), -sin(theta_in)/) !Incident E-field vector  (p-polarized)   
	real(kind = dp), dimension(3), public :: k_in_hat = (/sin(theta_in)*cos(phi_in), sin(theta_in)*sin(phi_in), cos(theta_in)/)		
	
	complex(kind = dp), parameter, public :: k0 = 2*PI/Lambda*(1.0 + im*0.0)
	real(kind = dp), parameter, public :: Omega = k0*c0		
	complex(kind = dp), parameter, public :: eps_r1 = 1.0 + im*0!*Epsilon_r of the space outside the object	
    complex(kind = dp), parameter, public :: eps_r2 =  (-17.2, -0.498)! (-8.0, - 1.66)!Au !!Ag *Epsilon_r of the object	
	complex(kind = dp), parameter, public :: eps_1 = eps_0*eps_r1
	complex(kind = dp), parameter, public :: eps_2 = eps_0*eps_r2
	complex(kind = dp), parameter, public :: k1 = k0*sqrt(eps_r1)
	complex(kind = dp), parameter, public :: k2 = k0*sqrt(eps_r2)
	real(kind = dp), parameter, public ::	my_1 = my_0
	real(kind = dp), parameter, public ::	my_2 = my_0

	!---------For Gaussian beam----------------------------------
	real(kind = dp), parameter, public :: w0 = 300e-6	
	character(len = 50), public :: illumination
	character(len = 50), public :: object_type
	
	!For flat-triangle meshed element:
	integer, public :: Nc, Ne
	integer, public :: n_disc
	real(kind = dp), public :: D_factor
	
	!---------------Input and Output ---------------------------
	
	character(len = 80), public :: Result_SEH
	character(len = 50), public :: Name_NodeVector 
    character(len = 50), public :: Name_Element 
	!real(kind = dp), allocatable, dimension(:, :), public :: Node_Vector
	
	contains
	
	
end module input_tri
