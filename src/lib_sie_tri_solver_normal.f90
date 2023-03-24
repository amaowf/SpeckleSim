  !Copyright (C) 2021  Liwei Fu <liwei.fu@ito.uni-stuttgart.de>
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
	
module lib_sie_tri_solver_normal
	!Using this version, objects with different material parameters 
	!can be calculated. 

	use omp_lib
	use lib_sie_math	
	use lib_sie_constants
	use lib_sie_data_container
	use lib_sie_type_function 
	use time_stamp
	
	use lib_sie_tri_data_container
	use lib_sie_tri_calculation_mod
	
	implicit none	
	
	private	
	
	public :: run_normal_tri_calculation	
	public :: lib_sie_tri_solver_normal_II
	public :: test_conical_illumination
	public :: test_Gaussian_illumination
	
	complex(dp), dimension(:), allocatable :: b_vec!
	character(len = 10) :: calculation_step 
		
	contains 
		
	subroutine run_normal_tri_calculation()	
		implicit none
		if (calc_p(6) .eq. 1) then
			calculation_step  = 'solver_I' 
			call lib_sie_tri_solver_normal_I()		
		else if (calc_p(6) .eq. 2) then
			calculation_step  = 'solver_II' 			
			call lib_sie_tri_solver_normal_II()
		else if (calc_p(6) .eq. 3) then
			calculation_step  = 'solver_III' 
			call lib_sie_tri_solver_normal_I()
			call lib_sie_tri_solver_normal_II()
		else 
		
			print*, 'Not a proper calculation step type!'
			call exit
		end if
	end subroutine run_normal_tri_calculation
  
	subroutine lib_sie_tri_solver_normal_I()
				
		implicit none
		integer :: m, n, counter, m_edge, p, q
		character (len=40) :: file_out
		
		real(dp) :: r_ref
		complex(dp), dimension(:,:), allocatable :: D_mat
		complex(dp) :: D_mat_tmp(2, 2)
		
		real :: test_start_sub, test_finish_sub, impedance_time, LU_time
        ! WALL-time
		INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub
		
		call set_parameters_tri()
		call data_import_tri()
		
		n_el = size(struc_tri%elements)
		m_pairs = size(struc_tri%neighbours)			
		m_edge = 2*m_pairs !1536
		print*,'m_pairs =', m_pairs
		print*,'n_el =', n_el
		
		call precalculation_fn(struc_tri)			
		allocate(D_mat(m_edge, m_edge)) !
		allocate(b_vec(m_edge))		
		call get_vector_b_tri(b_vec)		
	 
		print*,'Begin to calculate step-I'
		call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)		
		
		counter = 0	 
		!$omp parallel private (m, n, p, q)  &
		!$omp private (r_ref, D_mat_tmp) &
		!$omp shared (m_pairs, struc_tri, b_vec, D_mat, counter)				
		!$omp do
		!PRINT*, 'Number of threads =', omp_get_num_threads()
		do m = 1, m_pairs ! 			
			counter = counter+1
			if (m_pairs >10000) then		
				if (modulo(counter,500).eq.0) then
					print *, 'pair nr', counter, 'of total', m_pairs
				end if		
			else
				if (modulo(counter,100).eq.0) then
					print *, 'pair nr', counter, 'of total', m_pairs
				end if	
			end if		
			p = struc_tri%neighbours(m)%element(1)
			q = struc_tri%neighbours(m)%element(2)	
			r_ref = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)
			do n = 1, m_pairs				
				D_mat_tmp = lib_sie_tri_get_impedance_edges(m, n, r_ref)
				D_mat(m, n) = D_mat_tmp(1, 1)
				D_mat(m, n + m_pairs) = D_mat_tmp(1, 2)
				D_mat(m + m_pairs, n) = D_mat_tmp(2, 1)
				D_mat(m + m_pairs, n + m_pairs) = D_mat_tmp(2, 2)							
			end do !loop n
		end do !loop m
		!$omp end do
		!$omp end parallel		
		
		print*, 'Impedance elements are done'
		call timestamp()
		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)
		impedance_time = (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)/60		
		print '("Time for impedance calculation = ",f13.5," minutes")', impedance_time 
				  
      call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)
		
		call solver_LU(m_edge, D_mat, b_vec)
				
		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)
  
		LU_time = (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)/60
		
		print '("Time for LU Factorization = ",f13.5," minutes.")', LU_time		
		
		call set_evaluation_parameters_tri()		
			
		open (unit=206, file = file_name_output_I, action="write",status = 'replace')!!!!       	
			write (206, *) 'n_el               ', n_el
			write (206, *) 'm_pairs            ', m_pairs
			do m = 1, m_edge
				write (206, '(201(es19.12, tr5))') real(b_vec(m)), aimag(b_vec(m)) 
			end do 	
		close(206)		
		
		print*, 'Finished step-I ', file_name_output_I
		
		!file_out = 'Calculation_parameters_normal.txt'	
		if (pre_types%object .eq. 'sphere') then
			
		   write(file_out, '(I8.8)') tri_sp_parameters(1)%n_disc		
			file_out_parameters = 'ndisc_'//trim(file_out)//'_'//trim(file_out_parameters)
			
			open (unit = 206,file = file_out_parameters, action="write", status = 'replace')!!!
			
				write (206, *) 'object_type =                  ', pre_types%object
				write (206, *) 'formulation_type =             ', pre_types%formulation
				write (206, *) 'illumination_type =            ', pre_types%illumination
				write (206, *) 'beam_waist(um) (Gaussian beam) =   ', beam_waist*1e+6
			
				write (206, *) 'Wavelength (nm) =          ', illumination_p%lambda*1.0e+9
				write (206, *) 'sphere number =            ', number_objects
				write (206, *) 'Sphere diameter (um) =     ', tri_sp_parameters(1:number_objects)%D*1.0e+6
				write (206, *) 'Ndisc =                    ', tri_sp_parameters(1:number_objects)%n_disc
				write (206, *) 'N_el =                    ', n_el
				write (206, *) 'Sphere centroid (m) =      ', tri_sp_parameters(1:number_objects)%centroid				
				write (206, *) 'eps_r2 =                   ', eps_r2
				write (206, *) 'evaluation_type =			  ', pre_types%evaluation
				write (206, *) 'evaluation_parameter%dim_c(um) =   ', evaluation_parameter%dim_c
				write (206, *) 'evaluation_parameter%dim_a(um) =   ', (evaluation_parameter%dim_a(2)-evaluation_parameter%dim_a(1))*1.0e+6
				write (206, *) 'evaluation_parameter%dim_b(um) =   ', (evaluation_parameter%dim_b(2)- evaluation_parameter%dim_b(1))*1.0e+6
				write (206, *) 'illumination_p%theta_in =          ',	illumination_p%theta_in*180/PI
				write (206, *) 'illumination_p%phi_in =          ',	illumination_p%phi_in*180/PI
				write (206, *) 'Impedance time (min.)         ', impedance_time 
				write (206, *) 'LU time (min.)                ', LU_time
				write (206, *) 'm_pairs                 ', m_pairs
			close(206)
		else if (pre_types%object .eq. 'surface') then
		
			open (unit = 206,file = file_out_parameters, action="write", status = 'replace')!!!			
			
				write (206, *) 'object_type =                  ', pre_types%object
				write (206, *) 'formulation_type =             ', pre_types%formulation
				write (206, *) 'illumination_type =            ', pre_types%illumination
				write (206, *) 'beam_waist(nm) (Gaussian beam) =   ', beam_waist*1e+9
				write (206, *) 'Wavelength (nm) =           ', illumination_p%lambda*1.0e+9
				write (206, *) 'object number =             ', number_objects
				write (206, *) 'object dimension I (um) =   ', tri_sf_parameters(1:number_objects)%D(1)*1.0e+6
				write (206, *) 'object dimension II (um) =  ', tri_sf_parameters(1:number_objects)%D(2)*1.0e+6
				write (206, *) 'object centroid (m) =       ', tri_sf_parameters(1:number_objects)%centroid
				write (206, *) 'N_el =                      ', n_el
				write (206, *) 'eps_r2 =                    ', eps_r2
				write (206, *) 'evaluation_type =			  ', pre_types%evaluation
				write (206, *) 'evaluation_parameter%dim_c(um) =   ', evaluation_parameter%dim_c
				write (206, *) 'evaluation_parameter%dim_a(um) =   ',  (evaluation_parameter%dim_a(2)-evaluation_parameter%dim_a(1))*1.0e+6
				write (206, *) 'evaluation_parameter%dim_b(um) =   ',  (evaluation_parameter%dim_b(2)-evaluation_parameter%dim_b(1))*1.0e+6
				write (206, *) 'illumination_p%theta_in =          ',	illumination_p%theta_in*180/PI
				write (206, *) 'illumination_p%phi_in =          ',	illumination_p%phi_in*180/PI
				write (206, *) 'Impedance time(min.)           ', impedance_time 
				write (206, *) 'LU time (min.)                 ', LU_time
				write (206, *) 'm_pairs                 ', m_pairs
			close(206)		
		end if
		
	end subroutine lib_sie_tri_solver_normal_I
	
	subroutine lib_sie_tri_solver_normal_II()	
		implicit none		
		type(vector_c) , dimension(:), allocatable :: E_out, H_out!		
		
		!dummy
		real(dp), dimension(2) ::  bb_vec
		integer :: in_unit
		integer :: m, m_edge, n_el !			
		character (len=40) :: dummy !file_SEH,
		
		real :: test_start_sub, test_finish_sub
        ! WALL-time
		INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub
    	
		! Open the calculated coefficient from BEM
		!-----------------------------------------
		if (calc_p(6) .eq. 2 ) then			
			
			call set_parameters_tri()
			call data_import_tri()
			call precalculation_fn(struc_tri)
	
			m_edge = 2*m_pairs			
			allocate(b_vec(2*m_pairs))
			in_unit = 115
			open(unit = in_unit, file = trim(file_name_output_I))
				read(in_unit, *) dummy, n_el
				read(in_unit, *) dummy, m_pairs
				
				do m = 1, m_pairs*2
					read(in_unit, *) bb_vec
					b_vec(m) = bb_vec(1) + im*bb_vec(2)
				end do
			close(in_unit)
		else if ((pre_types%calculation .eq. 'Gmres_tri') .and. (calc_p(6) .eq. 3 )) then
			b_vec = simulation_data%vector_x		
		else if ((pre_types%calculation .eq. 'Normal_tri') .and. (calculation_step .ne. 'solver_II')) then	
			continue 
		else
			print*, 'No matching calculation type!'
			call exit
		end if
				
		call set_evaluation_parameters_tri()
		call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)
		
		print*, 'Begin to calculate step-II'		
		call lib_sie_tri_observation_field_calculation(b_vec, E_out, H_out)	
		
		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)

		print '("Field calculation time = ",f10.3," min.")', &
					(test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)	/60
		print*, 'Finished step-II'
	end subroutine lib_sie_tri_solver_normal_II

	subroutine lib_sie_tri_observation_field_calculation(S_EH, E_out, H_out)	!S_EH
		implicit none
		
		type(vector_c), dimension(:), allocatable, intent(out) :: E_out, H_out
		complex(dp), dimension(:), allocatable, intent(in) :: S_EH 
		!real(dp) :: r1, r2
		
		!dummy
		integer :: mm, Nop, n, ngp, counter		
		type(point) :: r_local
		type(vector_c), dimension(2) :: sum_aa
		type(vector_c) :: Hfield_in, Efield_in
		type(lib_sie_evaluation_point_type), dimension(:), allocatable :: input_field		
		real(dp), dimension(:), allocatable :: Energy_scat, theta_scatter
		real(dp) :: dx_a, center(3)
				
		counter=0 
		if (allocated(E_out)) then
			deallocate(E_out)
		end if
		if (allocated(H_out)) then
			deallocate(H_out)
		end if
		if (allocated(Energy_scat)) then
			deallocate(Energy_scat)
		end if
		
		if (allocated(Energy_scat)) then
			deallocate(Energy_scat)
		end if
		
		call get_evaluation_points_domain_new(evaluation_parameter, eps_r1, eps_r2, calc_p, geometry_p, surface_center, r_media)		
		
		call lib_sie_tri_input_field_calculation(input_field)			
		
		Nop = size(input_field)				
		allocate(E_out(Nop), H_out(Nop), Energy_scat(Nop))
				
		if ((pre_types%evaluation .eq. 'BRDF_n' ) .or. (pre_types%evaluation .eq. 'BRDF_p') )then
			allocate(theta_scatter(evaluation_parameter%N_dim(1)))			
			dx_a = (evaluation_parameter%dim_a(2) - evaluation_parameter%dim_a(1))/(evaluation_parameter%N_dim(1)-1)
			do mm = 1, evaluation_parameter%N_dim(1)
				theta_scatter(mm) = evaluation_parameter%dim_a(1) + dx_a*(mm-1)				
				!Intensity_GB_surface(m) =  !should be normalized by the incident energy.
			end do			
		end if			
		
		if ((pre_types%evaluation .eq. 'rcs_p') .or. (pre_types%evaluation .eq. 'rcs_n')) then		
			!$omp parallel private (mm, sum_aa, Efield_in, Hfield_in, ngp) &  
			!$omp shared (m_pairs, Nop, counter, E_out, H_out, S_EH, struc_tri, Energy_scat, r_media) 			
			!$omp do
			do mm = 1, Nop ! for all r
				E_out(mm)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if				
				ngp = 3				
				call Integration_scattered_field_tri(r_media(mm)%point, struc_tri, ngp, S_EH, Sum_aa) !
				E_out(mm)%vector =  sum_aa(1)%vector ! 
				H_out(mm)%vector =  sum_aa(2)%vector ! 
				Energy_scat(mm) = dot_product(E_out(mm)%vector, E_out(mm)%vector)*4*pi*vec_len(r_media(mm)%point)**2! &
				counter = counter + 1
			end do
			!$omp end do
			!$omp end parallel  
		else if ((pre_types%evaluation .eq. 'BRDF_n' ) .or. (pre_types%evaluation .eq. 'BRDF_p')) then		
			!$omp parallel private (mm, n, sum_aa, Efield_in, Hfield_in, ngp) &  
			!$omp shared (m_pairs, Nop, counter, E_out, H_out, S_EH, struc_tri, Energy_scat, r_media) 		
			!$omp do
			do mm = 1, Nop ! for all r
				E_out(mm)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if				
				ngp = 3
				call Integration_scattered_field_tri(r_media(mm)%point, struc_tri, ngp, S_EH, Sum_aa)
				E_out(mm)%vector =  sum_aa(1)%vector ! 
				H_out(mm)%vector =  sum_aa(2)%vector ! 							
				Energy_scat(mm) = dot_product(E_out(mm)%vector, E_out(mm)%vector)*vec_len(r_media(mm)%point)**2! &	
				Energy_scat(mm)= Energy_scat(mm)/(cos(illumination_p%theta_in)*cos(theta_scatter(mm)))
				counter = counter + 1
			end do
			!$omp end do
			!$omp end parallel  
			
		else if (pre_types%evaluation .eq. 'cly_field') then
			!$omp parallel private (mm, n, sum_aa, Efield_in, Hfield_in, ngp) &  
			!$omp shared (m_pairs, Nop, counter, E_out, H_out, r_media, S_EH, struc_tri)  &
			!$omp shared (input_field) !, nearfield_D
			!$omp do
			do mm = 1, Nop ! for all r
				E_out(mm)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if
				ngp = 3
				call Integration_scattered_field_tri(r_media(mm)%point, struc_tri, ngp, S_EH, sum_aa)					
				E_out(mm)%vector = input_field(mm)%e_field%vector*total_field + sum_aa(1)%vector
				H_out(mm)%vector = input_field(mm)%h_field%vector*total_field + sum_aa(2)%vector										
				counter = counter + 1
			end do	
			!$omp end do
			!$omp end parallel  		
		else
			!$omp parallel private (mm, n, sum_aa, Efield_in, Hfield_in,  ngp) &  
			!$omp shared (m_pairs, Nop, counter, E_out, H_out, r_media, S_EH, struc_tri)  &
			!$omp shared (input_field) !
			!$omp do
			
			do mm = 1, Nop ! for all r	
				E_out(mm)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)						
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if				
			!	
				if (r_media(mm)%media .eq. 2)then				
					ngp = r_media(mm)%ngp
				
					call Integration_scattered_field_tri_within_medium(r_media(mm), struc_tri, ngp, S_EH, sum_aa)
					if (abs(sum_aa(1)%vector(1)) .lt. 1.0e-20) then
						sum_aa(1)%vector(1) = (0.0, 0.0)
						sum_aa(2)%vector(1) = (0.0, 0.0)
					end if
					if (abs(sum_aa(1)%vector(2)) .lt. 1.0e-20) then
						sum_aa(1)%vector(2) = (0.0, 0.0)
						sum_aa(2)%vector(2) = (0.0, 0.0)
					end if
					
					if (abs(sum_aa(1)%vector(3)) .lt. 1.0e-20) then
						sum_aa(1)%vector(3) = (0.0, 0.0)
						sum_aa(2)%vector(3) = (0.0, 0.0)
					end if
					
					E_out(mm)%vector = sum_aa(1)%vector
					H_out(mm)%vector = sum_aa(2)%vector						
				else		
				   ngp = r_media(mm)%ngp
					!if ((abs(abs(r_media(mm)%point(1)) - 1.0e-6) .lt. 1e-15) .and. (abs(abs(r_media(mm)%point(3)) - 1.0e-6) .lt. 1e-15)) then
					
					call Integration_scattered_field_tri(r_media(mm)%point, struc_tri, ngp, S_EH, sum_aa)					
					
					E_out(mm)%vector = input_field(mm)%e_field%vector*total_field + sum_aa(1)%vector
					H_out(mm)%vector = input_field(mm)%h_field%vector*total_field + sum_aa(2)%vector	
					
					!print*, 'r_local=', r_media(mm)%point
					!!print*, 'E_out(mm)%vector=', E_out(mm)%vector
					!end if
					
				end if				
				counter = counter + 1				
			end do	
			!$omp end do
			!$omp end parallel  		
		
		end if		
				
		call timestamp()
		
		if ((pre_types%evaluation .eq. 'rcs_n') .or. (pre_types%evaluation .eq. 'rcs_p') .or. &		
			(pre_types%evaluation .eq. 'BRDF_p') .or. (pre_types%evaluation .eq. 'BRDF_n'))then
			
			open (unit = 203, file = file_name_output_Efield, action = "write",status = 'replace')
				do mm = 1, Nop        
					write (203, '(1001(E19.12, tr3))') Energy_scat(mm)        
				end do
			close(203)
		else
			
			open (unit = 203, file = file_name_output_Efield, action = "write",status = 'replace')
			do mm = 1, Nop        
					write (203, '(1001(E19.12, tr3))')  (real(E_out(mm)%vector(n)), n= 1, 3), (imag(E_out(mm)%vector(n)), n= 1, 3), (real(H_out(mm)%vector(n)), n= 1, 3), (imag(H_out(mm)%vector(n)), n= 1, 3)
					end do
			close(203)
					
		end if
		
		if (calc_p(1) .ne. 4)then
			if (allocated(struc_tri%points))then
				deallocate(struc_tri%points)
			end if
		
			if (allocated(struc_tri%x_c))then
				deallocate(struc_tri%x_c)
			end if
		
			if (allocated(struc_tri%midpoint))then
				deallocate(struc_tri%midpoint)
			end if	
		
			if (allocated(struc_tri%neighbours))then
				deallocate(struc_tri%neighbours)
			end if
				
			if (allocated(struc_tri%midpoint_edge))then
				deallocate(struc_tri%midpoint_edge)
			end if
				
			if (allocated(struc_tri%elements))then
				deallocate(struc_tri%elements)
			end if				
		
			if (allocated(b_vec))then
				deallocate(b_vec)
			end if				
		end if
		!open (unit = 203, file = 'test_domain.txt', action = "write",status = 'replace')
		!	do mm = 1, nop        
		!		write (203, '(1001(e19.12, tr3))') 	real(r_media(mm)%eps_r)
		!	end do
		!close(203)
		
		!open (unit = 203, file = 'test_domain.txt', action = "write",status = 'replace')
		!	do mm = 1, nop        
		!		write (203, '(1001(e19.12, tr3))') 	(real(input_field(mm)%e_field%vector(n)), n=1,3)
		!	end do
		!close(203)
		
		!
		!!Only input field 
		!********************************		!
		!open (unit = 203, file = file_name_output_Efield, action = "write",status = 'replace')
		!	do m = 1, Nop        
		!		write (203, '(1001(E19.12, tr3))') 	(real(input_field(m)%e_field%vector(n)), n= 1, 3), &
		!		(imag(input_field(m)%e_field%vector(n)), n= 1, 3)	
		!	end do
		!close(203)
		
		return
	end subroutine
	
	!Test for the conical_illumination(illumination, p_obj, r_local, M_obj, E_surf)	
	subroutine test_conical_illumination(p_obj)
		type(lib_sie_parameter_objective), intent(in) :: p_obj
		type(vector_c) :: E_surf
		real(dp) :: r_local(3), M_obj 
		integer :: Nop, m, n, in_unit, counter, Nx, Nz
		type(point), dimension(:), allocatable :: r_local_tt
		type(vector_c), dimension(:), allocatable :: E_out
		
		!dummy
		character (len=40) :: dummy !		
		real(dp), dimension(2) ::  bb_vec
		real(dp) :: z0, x0, y0, dx, dz
		!
		Nx = 201
		Nz = 201
		
		Nop = Nx*Nz
		dx = 4.0e-6/(Nx-1)
		dz = 4.0e-6/(Nz-1)
		
		allocate(r_local_tt(Nop), E_out(Nop))
		counter = 0!
		
		do m = 1, Nz
			do n = 1, Nx
				counter = (m-1)*Nz + n
				x0 = -2.0e-6 + dx*(n-1)
				!y0 = -1.0e-6 + 0.02e-6*m
				y0 = 0.0
				z0 = -2.0e-6 + dz*(m-1)
				r_local_tt(counter)%point(1:3) = (/x0, y0, z0/)
			end do
		end do
					
		do m = 1, Nop ! for all r
			r_local(1:3) =  r_local_tt(m)%point(1:3)		
			call conical_illumination(p_obj, r_local, illumination_p,  E_out(m))
		end do
		
		
		open (unit = 203, file = 'test_focus_field.txt', action = "write",status = 'replace')
			do m = 1, Nop        
				write (203, '(1001(E19.12, tr3))') (real(E_out(m)%vector(n)), n= 1, 3), (imag(E_out(m)%vector(n)), n= 1, 3)
			end do
		close(203)
		
	end subroutine	
	
	subroutine test_Gaussian_illumination()
		!type(lib_sie_parameter_objective), intent(in) :: p_obj
		type(vector_c) :: E_surf
		real(dp) :: r_local(3), M_obj 
		integer :: Nop, m, n, in_unit, counter, Nx, Nz
		type(point), dimension(:), allocatable :: r_local_tt
		type(vector_c), dimension(:), allocatable :: E_out
		
		!dummy
		character (len=40) :: dummy !file_SEH,		
		real(dp), dimension(2) ::  bb_vec
		real(dp) :: z0, x0, y0, dx, dz
		complex(dp) :: kd
		
		kd = k1
		!
		Nx = 201
		Nz = 201
		
		Nop = Nx*Nz
		dx = 4.0e-6/(Nx-1)
		dz = 4.0e-6/(Nz-1)
		
		allocate(r_local_tt(Nop), E_out(Nop))
		counter = 0!
		
		do m = 1, Nz
			do n = 1, Nx
				counter = (m-1)*Nz + n
				x0 = -2.0e-6 + dx*(n-1)
				!y0 = -1.0e-6 + 0.02e-6*m
				y0 = 0.0
				z0 = -2.0e-6 + dz*(m-1)
				r_local_tt(counter)%point(1:3) = (/x0, y0, z0/)
			end do
		end do
					
		do m = 1, Nop ! for all r
			r_local(1:3) =  r_local_tt(m)%point(1:3)		
			call Gaussian_beam(beam_waist, r_local, illumination_p, E_out(m), kd)!(p_obj, r_local, illumination_p,  E_out(m))				
		end do
				
		open (unit = 203, file = 'test_Gaussian_field.txt', action = "write",status = 'replace')
			do m = 1, Nop        
				write (203, '(1001(E19.12, tr3))') (real(E_out(m)%vector(n)), n= 1, 3), (imag(E_out(m)%vector(n)), n= 1, 3)
			end do
		close(203)
		
	end subroutine	
	
	
	!			
	!	
	
	
end module lib_sie_tri_solver_normal