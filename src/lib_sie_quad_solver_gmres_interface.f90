module lib_sie_quad_solver_gmres_interface
	use libmath
	
	use lib_sie_quad_mlfmm_interface	
	use lib_sie_quad_data_container		
	use time_stamp
	implicit none
	
	private
	public :: run_ml_fmm_GMRES_quad
	
	logical :: m_use_ml_fmm
	type(solver_gmres_parameter_type) :: m_gmres_parameter
   type(solver_gmres_callback_type) :: m_gmres_callback	
	integer(kind = 1), parameter :: unit_error	 = 115
	
   contains	
	
	subroutine run_ml_fmm_GMRES_quad()
		implicit none
		
		call set_parameters_quad()		
		call lib_sie_ml_fmm_precalculation_quad()
		call lib_sie_solver_gmres_run_quad()
		
	end subroutine
	
	subroutine lib_sie_solver_constructor_quad(use_ml_fmm)
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
            m_gmres_callback%calc_vector_b => lib_sie_ml_fmm_calculate_vector_b_quad				
      end if
				
      m_gmres_callback%get_vector_b => lib_sie_ml_fmm_get_vector_b_quad
		m_gmres_callback%get_vector_x_init => lib_sie_ml_fmm_get_init_vector_x_quad
      m_gmres_callback%save_vector_x => lib_sie_ml_fmm_save_vector_x_quad
		m_gmres_callback%internal%backward_error => zgmres_backward_error	
			
   end subroutine lib_sie_solver_constructor_quad
	
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
			rv%max_iterations = 100
			rv%restart = 50
			rv%convergence_tolerance = 1D-3
			rv%orthogonalization_scheme = GMRES_ORTHOGONALIZATION_SCHEME_IMGS
			rv%use_recurence_formula_at_restart = .false.
			rv%residual_calc_explicitly = .true.
			if (pre_types%calculation .eq. 'Gmres_quad')then
				rv%no_of_elements_vector_x = m_edge*2
			else
				rv%no_of_elements_vector_x = m_pairs*2
			end if
			
			rv%preconditioning = GMRES_PRECONDITIONING_NO
			
	end function lib_sie_solver_gmres_get_parameter	
	
	subroutine zgmres_backward_error(iteration, backward_error_arnoldi, backward_error_true)
		use lib_sie_tri_data_container
		implicit none
			
      ! dummy
      integer, intent(in) :: iteration
      double precision, intent(in) :: backward_error_arnoldi
      double precision, intent(in) :: backward_error_true

		write (unit_error, '(1(i10, tr5), 1(es19.12, tr5))') iteration, backward_error_arnoldi
		  
	end subroutine zgmres_backward_error	
		
	subroutine lib_sie_solver_gmres_run_quad()
			use lib_math_solver_GMRES
         implicit none
			
			logical :: use_ml_fmm
			integer :: i
			character (len=40) :: file_out			
			real :: test_start_sub, test_finish_sub
			real(dp) :: GMRES_time
        ! WALL-time
			INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub
		   
			call system_clock(test_count_start_sub, test_count_rate_sub)
			call cpu_time(test_start_sub)
			call timestamp()
			use_ml_fmm = .true.			
			
			call lib_sie_solver_constructor_quad(use_ml_fmm)			
			
			print*, 'm_gmres_parameter%max_iterations', m_gmres_parameter%max_iterations
			
			open (unit = unit_error, file = 'iteration_error.txt', action="write", status = 'replace')					
			
			call lib_math_solver_gmres_run(m_gmres_parameter, m_gmres_callback, .true.)
  
			call cpu_time(test_finish_sub)
         call system_clock(test_count_finish_sub, test_count_rate_sub)
			
			GMRES_time = (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)
		  		
			print '("   GMRES-Time = ",f10.3," seconds.")', GMRES_time
																		 
			file_out = 'Result_SEH_quad_'//trim(adjustl(pre_types%formulation))//'.txt'	
			open (unit=206, file = file_out, action="write",status = 'replace')		
				write (206, *) 'n_el,', n_el
				write (206, *) 'm_edge', m_edge
				do i = 1, m_edge*2
					write (206, '(201(es19.12, tr5))') simulation_data%vector_x(i)
				end do
			close(206)
			
			close(unit_error)
		
      end subroutine

end module lib_sie_quad_solver_gmres_interface