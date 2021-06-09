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
	
module lib_sie_solver_gmres_interface
	use libmath		
	use lib_sie_tri_mlfmm_interface
	use lib_sie_tri_data_container	
		
	use time_stamp
	implicit none

   private
		
	public :: run_ml_fmm_GMRES_tri
   logical :: m_use_ml_fmm
		
   type(solver_gmres_parameter_type) :: m_gmres_parameter
   type(solver_gmres_callback_type) :: m_gmres_callback
	
	integer(kind = 1), parameter :: unit_error	 = 115
   contains
	
	subroutine run_ml_fmm_GMRES_tri()
		implicit none
		
		call set_parameters_tri()		
		call lib_sie_ml_fmm_precalculation_tri()
		call lib_sie_solver_gmres_run()
		!call lib_sie_tri_solver_normal_II()
	end subroutine

      ! Argument
      ! ----
      !   use_ml_fmm: logical, optional (std: .false.)
      !       enables the ML-FMM
	subroutine lib_sie_solver_constructor(use_ml_fmm)
		implicit none
      ! dummy
      logical, intent(in), optional :: use_ml_fmm

      !auxiliary
		! if any
      m_use_ml_fmm = .false.
      if (present(use_ml_fmm)) m_use_ml_fmm = use_ml_fmm

      ! init
      m_gmres_parameter = lib_sie_solver_gmres_get_parameter()!
      if (m_use_ml_fmm) then
            m_gmres_callback%calc_vector_b => lib_sie_ml_fmm_calculate_vector_b_tri				
				m_gmres_callback%preconditioner_left => lib_sie_ml_fmm_preconditioner_left
      end if
				
      m_gmres_callback%get_vector_b => lib_sie_ml_fmm_get_vector_b_tri
      m_gmres_callback%get_vector_x_init => lib_sie_ml_fmm_get_init_vector_x_tri
      m_gmres_callback%save_vector_x => lib_sie_ml_fmm_save_vector_x_tri
		m_gmres_callback%internal%backward_error => zgmres_backward_error	
   end subroutine lib_sie_solver_constructor
		
	subroutine zgmres_backward_error(iteration, backward_error_arnoldi, backward_error_true)
		use lib_sie_tri_data_container
		implicit none
			
      ! dummy
      integer, intent(in) :: iteration
      double precision, intent(in) :: backward_error_arnoldi
      double precision, intent(in) :: backward_error_true

		write (unit_error, '(1(i10, tr5), 1(es19.12, tr5))') iteration, backward_error_arnoldi
		  
	end subroutine zgmres_backward_error		
		
	!
	!GMRES_ORTHOGONALIZATION_SCHEME_MGS = 0
	!GMRES_ORTHOGONALIZATION_SCHEME_IMGS = 1
	!GMRES_ORTHOGONALIZATION_SCHEME_CGS = 2
	!GMRES_ORTHOGONALIZATION_SCHEME_ICGS = 3
	!GMRES_PRECONDITIONING_LEFT = 1
	!GMRES_PRECONDITIONING_NO = 0
	function lib_sie_solver_gmres_get_parameter()result(rv)
		implicit none
			
		type(solver_gmres_parameter_type) :: rv
			
		rv%use_initial_guess = .true.				
		rv%max_iterations = max_iterations!12500
		rv%restart = restart! 12500
		rv%convergence_tolerance = convergence_tolerance! 1D-3
		rv%orthogonalization_scheme = GMRES_ORTHOGONALIZATION_SCHEME_IMGS
		rv%use_recurence_formula_at_restart = .false.
		rv%residual_calc_explicitly = .true.
		rv%no_of_elements_vector_x = m_pairs*2
		if (precondition_left) then
			rv%preconditioning = GMRES_PRECONDITIONING_LEFT
		else 
			rv%preconditioning = GMRES_PRECONDITIONING_NO		
		end if
		
	end function lib_sie_solver_gmres_get_parameter

   subroutine lib_sie_solver_gmres_run()
		use lib_math_solver_GMRES			
		use lib_sie_tri_data_container
      implicit none
			
		logical :: use_ml_fmm
		integer :: i
		!character(len = 100) :: file_out	
		real :: test_start_sub, test_finish_sub
		real(dp) :: GMRES_time
      ! WALL-time
		INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub
		   
		call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)
		call timestamp()
		use_ml_fmm = .true.			
		call lib_sie_solver_constructor(use_ml_fmm)
		
		print*, 'm_gmres_parameter%max_iterations', m_gmres_parameter%max_iterations
		if (m_gmres_parameter%preconditioning .eq. 1) then
			file_out_error = trim(file_out_error)//'_LeftP.txt'	
		else 
			file_out_error = trim(file_out_error)//'_NoP.txt'	
		end if		
		
		open (unit = unit_error, file = file_out_error, action="write", status = 'replace')					
			
		call lib_math_solver_gmres_run(m_gmres_parameter, m_gmres_callback, .true.)

		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)
			
		GMRES_time = (test_count_finish_sub-test_count_start_sub) &
                                                      / real(test_count_rate_sub)/60
		  		
		print '("GMRES-Time = ",f15.3," mins.")', GMRES_time
		
		deallocate(preCalculated_TL%TL, preCalculated_TL%my_LT)
		
		
		open (unit=206, file = file_name_output_I, action="write",status = 'replace')		
			write (206, *) 'n_el,', n_el
			write (206, *) 'm_pairs', m_pairs
			do i = 1, m_pairs*2
				write (206, '(201(es19.12, tr5))') simulation_data%vector_x(i)
			end do
		close(206)		
		close(unit_error)	
		
				
		if (pre_types%object .eq. 'sphere') then			
			open (unit = 206,file = file_out_parameters, action="write", status = 'replace')!!
				write (206, *) 'object_type =                  ', pre_types%object
				write (206, *) 'formulation_type =             ', pre_types%formulation
				write (206, *) 'illumination_type =            ', pre_types%illumination
				write (206, *) 'beam_waist(nm) (Gaussian beam) =   ', beam_waist*1e+9
				write (206, *) 'Wavelength (nm) =          ', illumination_p%lambda*1.0e+9
				write (206, *) 'sphere number =            ', number_objects
				write (206, *) 'Sphere diameter (nm) =     ', tri_sp_parameters(1:number_objects)%D*1.0e+9
				write (206, *) 'Ndisc =                    ', tri_sp_parameters(1:number_objects)%n_disc
				write (206, *) 'N_el =                    ', n_el
				write (206, *) 'Sphere centroid (m) =      ', tri_sp_parameters(1:number_objects)%centroid
				write (206, *) 'eps_r2 =                   ', eps_r2
				write (206, *) 'evaluation_type =			  ', pre_types%evaluation
				write (206, *) 'evaluation_parameter%dim_c(nm) =   ', evaluation_parameter%dim_c
				write (206, *) 'evaluation_parameter%dim_a(nm) =   ', evaluation_parameter%dim_a*1.0e+9
				write (206, *) 'evaluation_parameter%dim_b(nm) =   ', evaluation_parameter%dim_b*1.0e+9
				write (206, *) 'illumination_p%theta_in =          ',	illumination_p%theta_in*180/PI
				write (206, *) 'illumination_p%phi_in =          ',	illumination_p%phi_in*180/PI
				write (206, *) 'Sphere n_disc =                ', tri_sp_parameters%n_disc
				write (206, *) 'Sphere sampling point =        ', truncation_number
				write (206, *) 'GMRES_time (min.) =               ', GMRES_time
				write (206, *) 'Iteration procedure =          ', m_gmres_parameter%orthogonalization_scheme				
				write (206, *) 'convergence_tolerance =        ', m_gmres_parameter%convergence_tolerance
				write (206, *) 'l_max =                        ', box_at_l%lmax
				write (206, *) 'dc at l_max =                  ', box_at_l%dc
				write (206, *) 'tree_s_opt  =                  ', tree_s_opt				
			close(206)
	else if (pre_types%object .eq. 'surface') then
		open (unit = 206,file = file_out_parameters, action="write", status = 'replace')!!!
			write (206, *) 'object_type =                  ', pre_types%object
			write (206, *) 'formulation_type =             ', pre_types%formulation
			write (206, *) 'illumination_type =            ', pre_types%illumination
			write (206, *) 'beam_waist(nm) (Gaussian beam) =   ', beam_waist*1e+9
			write (206, *) 'Wavelength (nm) =            ', illumination_p%lambda*1.0e+9
			write (206, *) 'object number =             ', number_objects
			write (206, *) 'object dimension I (nm) =   ', tri_sf_parameters(1:number_objects)%D(1)*1.0e+9
			write (206, *) 'object dimension II (nm) =  ', tri_sf_parameters(1:number_objects)%D(2)*1.0e+9
			write (206, *) 'object centroid (m) =       ', tri_sf_parameters(1:number_objects)%centroid
			write (206, *) 'N_el =                    ', n_el
			write (206, *) 'eps_r2 =                    ', eps_r2
			write (206, *) 'evaluation_type =			  ', pre_types%evaluation
			write (206, *) 'evaluation_parameter%dim_c(nm) =   ', evaluation_parameter%dim_c
			write (206, *) 'evaluation_parameter%dim_a(nm) =   ', evaluation_parameter%dim_a*1.0e+9
			write (206, *) 'evaluation_parameter%dim_b(nm) =   ', evaluation_parameter%dim_b*1.0e+9
			write (206, *) 'illumination_p%theta_in =          ',	illumination_p%theta_in*180/PI
			write (206, *) 'illumination_p%phi_in =          ',	illumination_p%phi_in*180/PI
			write (206, *) 'Sphere sampling point =        ', truncation_number_default
			write (206, *) 'GMRES_time (min.) =               ', GMRES_time
			write (206, *) 'Iteration procedure =          ', m_gmres_parameter%orthogonalization_scheme				
			write (206, *) 'convergence_tolerance =        ', m_gmres_parameter%convergence_tolerance
			write (206, *) 'l_max =                        ', box_at_l%lmax
			write (206, *) 'dc at l_max =                  ', box_at_l%dc
			write (206, *) 'tree_s_opt  =                  ', tree_s_opt
			write (206, *) 'm_pairs  =                  ', m_pairs
		close(206)		
	end if
			
			
   end subroutine

end module lib_sie_solver_gmres_interface