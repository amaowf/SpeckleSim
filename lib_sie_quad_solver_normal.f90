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
	
!***********************************************************************
!  PURPOSE:  Solving surface tangential fields using MOM + Galerkin method and
!            8-node 10-edge quadrilateral elements!
!****************************************************************************
module lib_sie_quad_solver_normal
	use omp_lib
	use time_stamp
	use lib_sie_math
	use lib_sie_constants
	use lib_sie_data_container
	use lib_sie_type_function
 !
	use lib_sie_quad_singularity
	use lib_sie_quad_calculation_mod
	use lib_sie_quad_data_container	
		 
	implicit none
	
	private
	
	public :: run_normal_quad_calculation	
	character(len = 10) :: calculation_step 
	complex(dp), dimension(:), allocatable :: V_EH!
		
	contains
	
	subroutine run_normal_quad_calculation()	
		implicit none

		if (calc_p(6) .eq. 1) then
			calculation_step  = 'solver_I' 			
			call lib_sie_quad_solver_normal_I()
		else if (calc_p(6) .eq. 2) then			
			calculation_step  = 'solver_II' 
			call lib_sie_quad_solver_normal_II()			
		else if (calc_p(6) .eq. 3) then
			calculation_step  = 'solver_III'
			call lib_sie_quad_solver_normal_I()
			call lib_sie_quad_solver_normal_II()			
		else 
			print*, 'Not a proper calculation step type!'
			call exit
		end if	
	end subroutine run_normal_quad_calculation
    
	subroutine lib_sie_quad_solver_normal_I()
		implicit none    
		integer :: m, s, n, t, p, N_LU
		integer :: ii, jj, counter
		integer :: factor_II, pI_a, pI_q
		integer :: ngp, loc_corner(3), loc_edge(3)
    
		!Parameters used during the matrix calculation		
		complex(dp), dimension(:, :), allocatable :: Sum_KA, Sum_KDE, Sum_KDH		
		complex(dp) :: eta_a
		character(len = 20) :: type_calculation
		complex(dp), dimension(3, 10) :: Sum_aa
		character (len=80) :: file_out
	
		integer  :: N_gesv 
		external    ZGESV ! 
		external    CGESV
				
		complex(dp), dimension(:, :), allocatable :: M_mm, M_em, M_me	
		complex(dp), dimension(:, :), allocatable :: MM		
		 
		real :: test_start_sub, test_finish_sub, impedance_time, LU_time
		
		! WALL-time
		INTEGER :: test_count_start_sub, test_count_finish_sub, test_count_rate_sub
	
		!Initialize variables	
		print*, 'Initialize the structure!'
		call set_parameters_quad()
		call data_import_quad()    !Input the meshing data of the object
		n_el = size(struc_quad%elements)
		print*, 'n_el =', n_el
		
		call Neighbours_EdgeCorner_Full(struc_quad, n_el) !Find out the overlapped edges and corners
		call Sorting_edge(struc_quad, M_edge, Edge_I, Edge_II) 
		!
		N_gesv = 10*n_el
	 
		allocate (Sum_KA(M_edge, N_gesv))
		allocate (Sum_KDE(M_edge, N_gesv))
		allocate (Sum_KDH(M_edge, N_gesv))
	 
		allocate (M_me(M_edge, M_edge))   
		allocate (M_em(M_edge, M_edge))   
		allocate (M_mm(M_edge, M_edge))   
		allocate (MM(2*M_edge, 2*M_edge))		
		Sum_KA(1:m_edge,1:n_gesv) = (0.0, 0.0) 
		Sum_KDE(1:m_edge,1:n_gesv) = (0.0, 0.0)             
		Sum_KDH(1:m_edge,1:n_gesv) = (0.0, 0.0) 
		
		call get_vector_b_quad(edge_I, V_EH)
		
		Sum_aa(1:3, 1:10) = 0.0d0
  
		print*, 'Cores this program has access to:', omp_get_num_procs( ) 		
		print*, 'M_edge=', M_edge
		print*, 'Begin to calculate step-I: surface current'
		call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)
		counter = 0
		
		!$omp parallel private(n, p, s, t) & !
		!$omp private(ngp, Sum_aa) &
		!$omp private(type_calculation, loc_corner, loc_edge) &
		!$omp shared(struc_quad, M_edge, Sum_KA, Sum_KDE, Sum_KDH, counter) 	
		!$omp do
		 do s = 1, M_edge				
			counter = counter+1
			if (M_edge >1000) then		
				if (modulo(counter,500).eq.0) then
					print *, 'pair nr', counter, 'of total', M_edge
				end if		
			else
				if (modulo(counter,100).eq.0) then
					print *, 'pair nr', counter, 'of total', M_edge
				end if	
			end if	
			do t = 1, n_el !Element q  
				call impedance_edge_versus_element(s, t, sum_aa)
				do n = 1, 10
					p = (t-1)*10 + n
					Sum_KA(s, p) =  Sum_aa(1, n)
					Sum_KDE(s, p) = Sum_aa(2, n) 
					Sum_KDH(s, p) = Sum_aa(3, n) !combined KB and KC
				end do ! n-loop			
		  end do ! t-loop   
		end do ! s-loop 
		!$omp end do
		!$omp end parallel
		! 
		print*, 'Impedance elements are done'
		
		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)
		impedance_time = (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)/60		
		print '("  Time for impedance calculation = ",f15.5," mins.")', impedance_time 
				  
      call system_clock(test_count_start_sub, test_count_rate_sub)
		call cpu_time(test_start_sub)
	 
		! Matrix assembly          
		!-----------------------------------------------------------------
		do s = 1, M_edge			
			do t = 1, M_edge
				ii = edge_I(t)%edge_p(1)
				pI_a = edge_I(t)%edge_p(2) 	
				jj = edge_I(t)%edge_p(3)
				pI_q = edge_I(t)%edge_p(4)
								
				factor_II = edge_I(t)%edge_p(5)
				p = (ii-1)*10 + pI_a
				n = (jj-1)*10 + pI_q
				
				!jj could be 0 since 5th and 10th edges are not overlaped
				if (factor_II .NE. 0) then   					
					M_mm(s, t) = Sum_KA(s, p) + Sum_KA(s, n)*factor_II
					M_me(s, t) = Sum_KDH(s, p) + Sum_KDH(s, n)*factor_II
					M_em(s, t) = Sum_KDE(s, p) + Sum_KDE(s, n)*factor_II
				else   
					M_mm(s, t) = Sum_KA(s, p)  !
					M_me(s, t) = Sum_KDH(s, p) 
					M_em(s, t) = Sum_KDE(s, p) 
				end if  
			end do
		end do
		deallocate (Sum_KA)
		deallocate (Sum_KDE)
		deallocate (Sum_KDH)  
			
		if (pre_types%formulation .eq. 'PMCHWT')then!
			MM(1 : M_edge, 1 : M_edge) = M_mm ! M_mm = M_ee
			MM(1 : M_edge, M_edge+1 : 2*M_edge) = M_em
			MM(M_edge + 1 : 2*M_edge, 1:M_edge) = M_me
			MM(M_edge + 1 : 2*M_edge, M_edge+1 : 2*M_edge) = M_mm			
		else if (pre_types%formulation .eq. 'ICTF')then!
			eta_a = (eta_1 + eta_2)/2
			MM(1 : M_edge, 1 : M_edge) = M_mm/eta_a
			MM(1 : M_edge, M_edge+1 : 2*M_edge) = M_em/eta_a
			MM(M_edge + 1 : 2*M_edge, 1:M_edge) = M_me*eta_a
			MM(M_edge + 1 : 2*M_edge, M_edge+1 : 2*M_edge) = M_mm*eta_a
		else if (pre_types%formulation .eq. 'MCTF')then!	
			MM(1 : M_edge, 1 : M_edge) = M_mm
			MM(1 : M_edge, M_edge+1 : 2*M_edge) = M_em
			MM(M_edge + 1 : 2*M_edge, 1:M_edge) = M_me*eta_1*eta_2
			MM(M_edge + 1 : 2*M_edge, M_edge+1 : 2*M_edge) = M_mm*eta_1*eta_2
		else 
			print*, "Not a proper formulation!"
		end if
		
		deallocate (M_mm)
		deallocate (M_me)
		deallocate (M_em)
  
		N_LU = 2*M_edge      
		call solver_LU(N_LU, MM, V_EH) 
		call cpu_time(test_finish_sub)
      call system_clock(test_count_finish_sub, test_count_rate_sub)		
		LU_time = (test_count_finish_sub-test_count_start_sub) &
                                                       / real(test_count_rate_sub)/60
		print*, 'Linear-system is solved'
		print '(" Time for LU factorization = ",f15.5," mins.")', LU_time
  
		!Save the surface current coefficients
		!-------------------------------------------------------------		
		open (unit = 206,file = file_name_output_I, action="write",status = 'replace')!!!!   
		do m = 1, 2*M_edge
			write (206, '(201(es19.12, tr5))') real(V_EH(m)), aimag(V_EH(m)) 
		end do    
		close(206)
		
		if (pre_types%object .eq. 'sphere') then	
			
			open (unit = 206,file = file_out_parameters, action="write",status = 'replace')!		
			write (206,'(A, f9.4)') 'Wavelength (nm) =                            ', illumination_p%lambda*1.0e+9
			write (206,'(A, f9.4)') 'Sphere diameter (nm) =                       ', quad_sp_parameters(1:number_objects)%D*1.0e+9
			write (206,'(A, I6)')   'Number of ements: Nt =                       ', quad_sp_parameters(1:number_objects)%Nt
			write (206,'(A, I6)')   'Number of nodes (corners): Np =              ', quad_sp_parameters(1:number_objects)%Np
			write (206,'(A, 3f9.4)')'Centroid of the spheres                      ', quad_sp_parameters(1:number_objects)%centroid
			write (206, '(A, f9.4)')'Impedance time(mins.)                            ', impedance_time 
			write (206, *) 'eps_r2 =                    ', eps_r2
			write (206, "(A, es19.12)") 'LU time (mins.)                    ', LU_time	
			write (206, *) 'm_edge  =                  ', m_edge
			close(206)
				
		else if (pre_types%object .eq. 'surface') then			
			
			open (unit = 206,file = file_out_parameters, action="write",status = 'replace')!!	
			
			write (206, '(A, f9.4)') 'Wavelength (nm) =              ', illumination_p%lambda*1.0e+9
			write (206, "(A, f9.4)") 'size of the surface (nm) =     ', quad_sf_parameters(1:number_objects)%D(1)*1.0e+9
			write (206, "(A, I6)")   'Number of ements: Nt =         ', quad_sf_parameters(1:number_objects)%Nt
			write (206, "(A, I6)")   'Number of nodes (corners): Np =', quad_sf_parameters(1:number_objects)%Np
			write (206, "(A, 3f9.4)")'Centroid of the sphere         ', quad_sf_parameters(1:number_objects)%centroid
			write (206, "(A, f10.4)")'Impedance time(mins.)              ', impedance_time 
			write (206, *) 'eps_r2 =                    ', eps_r2
			write (206, "(A, es19.12)") 'LU time (mins.)                          ', LU_time
			write (206, *) 'm_edge  =                  ', m_edge
			close(206)
							
		else
			print*, 'Not a defined objects'
		end if
		
	end subroutine lib_sie_quad_solver_normal_I
	
	subroutine lib_sie_quad_solver_normal_II()
		implicit none  

		type(vector_c) , dimension(:), allocatable :: E_out, H_out!
		complex(dp), allocatable, dimension(:, :) :: S_EH		
		integer :: i, j, n, m 
		integer:: factor_II, pI_a, pII_a
		integer :: m_edge, in_unit
		
		character (len=50) :: file_name		
		real(dp), dimension(2) ::  bb_vec !
	
		if (calculation_step .eq. 'solver_II') then	

			call set_parameters_quad()	
			call data_import_quad()			
			n_el = size(struc_quad%elements)
			call Neighbours_EdgeCorner(struc_quad, n_el)
			call Sorting_edge(struc_quad, M_edge, Edge_I, Edge_II)
			print*, 'M_edge = ', M_edge
		
			allocate (V_EH(2*M_edge))
			allocate (S_EH(n_el*2, 10))
			
			!Import the results from Program_I
			!----------------------------------------------------
			in_unit = 30
			
			open(unit = in_unit, file=file_name_output_I,status='old',action='read')
				do n = 1, M_Edge*2			
					read(in_unit,*) bb_vec
					V_EH(n) = bb_vec(1) + im*bb_vec(2)
				end do
			close(in_unit)
		else
			allocate (S_EH(n_el*2, 10))
			call Sorting_edge(struc_quad, M_edge, Edge_I, Edge_II)
		end if
		
		do m = 1, M_edge        
			i = edge_I(m)%edge_p(1) !element number
			pI_a = edge_I(m)%edge_p(2) !edge number			
			j = edge_I(m)%edge_p(3) !element index which has overlapped edge with edge-s
			pII_a = edge_I(m)%edge_p(4) !edge index overlapped with the edge-s 
			factor_II = edge_I(m)%edge_p(5) !edge direction 
			S_EH(i, pI_a) = V_EH(m)
			S_EH(i + n_el, pI_a) = V_EH(m + m_edge)
			if (j .NE. 0) then
				S_EH(j, pII_a) = S_EH(i, pI_a)*factor_II
				S_EH(j + n_el, pII_a) = S_EH(i + n_el, pI_a)*factor_II
			end if    
		end do
	
		!Input parameters for calculating fields in space 
		!--------------------------------------------------------------------
		call set_evaluation_parameters_quad()		
		call timestamp()
		call lib_sie_quad_observation_field_calculation(S_EH, E_out, H_out)!
		call timestamp()
		print*, 'Step-II is done'
		
	end subroutine 
		
	subroutine lib_sie_quad_observation_field_calculation(S_EH, E_out, H_out)!
		implicit none
	
		type(vector_c), dimension(:), allocatable, intent(out) :: E_out, H_out
		complex(dp), dimension(:, :), intent(in) :: S_EH 
		
		!dummy
		real(dp), dimension(8, 3) :: Element_q ! 		
		integer :: m, Nop, n, ngp, counter, s, object_index
		integer(kind = 1) :: tot_field
		type(point) :: r_loc
		type(point), dimension(:), allocatable :: r_local
		type(vector_c), dimension(2) :: sum_a	
		type(lib_sie_evaluation_point_type), dimension(:), allocatable :: input_field		
		real(dp), dimension(:), allocatable :: Energy_scat,theta_scatter
		character(len = 50) :: file_name 
		real(dp) :: dx_a
		
		counter=0 
		
		call get_evaluation_points(evaluation_parameter, pre_types%evaluation, r_local)
		
		input_field = lib_sie_quad_input_field_calculation(r_local)	
		tot_field = evaluation_parameter%tot_field
		Nop = size(input_field)
		allocate(E_out(Nop), H_out(Nop), Energy_scat(Nop))
		print*, 'Nop=', Nop
		
		if (pre_types%evaluation .eq. 'BRDF')then
			allocate(theta_scatter(evaluation_parameter%N_dim(1)))			
			dx_a = (evaluation_parameter%dim_a(2) - evaluation_parameter%dim_a(1))/(evaluation_parameter%N_dim(1)-1)
			do m = 1, evaluation_parameter%N_dim(1)
				theta_scatter(m) = evaluation_parameter%dim_a(1) + dx_a*(m-1)
			end do
		end if
		
		if ((pre_types%evaluation .eq. 'rcs_p') .or. (pre_types%evaluation .eq. 'rcs_n')) then
		
			!$omp parallel private( m, n, s, Element_q) & !
			!$omp shared (counter, Nop, E_out, H_out, struc_quad, S_EH, tot_field) & !, nearfield_distance
			!$omp private (sum_a, r_loc, ngp, object_index)
			!$omp do
    
			do m = 1, Nop
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if
				E_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				H_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)			
				r_loc =  input_field(m)%coordinate
				ngp = 5
				do s = 1, n_el !
					do n = 1, 8      
						Element_q(n, 1:3) = struc_quad%node_vector(struc_quad%elements(s)%vertices(n))%point(1:3)
					end do 
					object_index = 1
					call Integration_scattered_field_quad(r_loc%point, Element_q, ngp, S_EH(s, :), S_EH(s + n_el, :), sum_a)
					E_out(m)%vector = E_out(m)%vector + sum_a(1)%vector
					H_out(m)%vector = H_out(m)%vector + sum_a(2)%vector
				end do !n-loop	
				Energy_scat(m) = dot_product(E_out(m)%vector, E_out(m)%vector)*4*pi*vec_len(r_loc%point)**2 !&
						!/dot_product(input_field(m)%e_field%vector, input_field(m)%e_field%vector)! this might not work for Gaussian beam
				counter = counter + 1
			end do !    
			!$omp end do  
			!$omp end parallel	
		elseif (pre_types%evaluation .eq. 'BRDF') then
		
			!$omp parallel private( m, n, s, Element_q) & !
			!$omp shared (counter, Nop, E_out, H_out, struc_quad, S_EH, tot_field) & !, nearfield_distance
			!$omp private (sum_a, r_loc, ngp, object_index)
			!$omp do
    
			do m = 1, Nop
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if
				E_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				H_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)			
				r_loc =  input_field(m)%coordinate
				ngp = 5
				do s = 1, n_el !
					do n = 1, 8      
						Element_q(n, 1:3) = struc_quad%node_vector(struc_quad%elements(s)%vertices(n))%point(1:3)
					end do 
					object_index = 1
					call Integration_scattered_field_quad(r_loc%point, Element_q, ngp, S_EH(s, :), S_EH(s + n_el, :), sum_a)
					E_out(m)%vector = E_out(m)%vector + sum_a(1)%vector
					H_out(m)%vector = H_out(m)%vector + sum_a(2)%vector
				end do !n-loop	
				Energy_scat(m) = dot_product(E_out(m)%vector, E_out(m)%vector)*4*pi*vec_len(r_loc%point)**2 !&
						!/dot_product(input_field(m)%e_field%vector, input_field(m)%e_field%vector)!	
				Energy_scat(m)= Energy_scat(m)/(cos(illumination_p%theta_in)*cos(theta_scatter(m)))
				
				counter = counter + 1
			end do !    
			!$omp end do  
			!$omp end parallel
		else
			!$omp parallel private( m, n, s, Element_q) & !
			!$omp shared (counter, Nop, E_out, H_out, struc_quad, S_EH, tot_field) & !, nearfield_distance
			!$omp private (sum_a, r_loc, ngp, object_index)
			!$omp do
    
			do m = 1, Nop
				if (mod(counter, Nop/10).eq. 0 .and. counter .ne. 0)then
					write (*, "(i5,a1)")int(ceiling(counter*100.0/Nop)),'%'
				end if
				E_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				H_out(m)%vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)			
				r_loc =  input_field(m)%coordinate			
				ngp = input_field(m)%ngp
			
				if (ngp .eq. 50) then					
					E_out(m)%vector(1:3) = (0.0d0, 0.0d0)
					H_out(m)%vector(1:3) = (0.0d0, 0.0d0)	
				else
					do s = 1, n_el !
						do n = 1, 8      
							Element_q(n, 1:3) = struc_quad%node_vector(struc_quad%elements(s)%vertices(n))%point(1:3)
						end do 						
						call Integration_scattered_field_quad(r_loc%point, Element_q, ngp, S_EH(s, :), S_EH(s + n_el, :), sum_a)
						E_out(m)%vector = E_out(m)%vector + sum_a(1)%vector
						H_out(m)%vector = H_out(m)%vector + sum_a(2)%vector
					end do !n-loop							
					E_out(m)%vector = input_field(m)%e_field%vector*tot_field - E_out(m)%vector
					H_out(m)%vector = input_field(m)%h_field%vector*tot_field - H_out(m)%vector
				end if	
				counter = counter + 1
			end do !    
			!$omp end do  
			!$omp end parallel
			end if
			
		if ((pre_types%evaluation .eq. 'rcs_n') .or. (pre_types%evaluation .eq. 'rcs_p'))then
		   
			open (unit = 203, file = file_name_output_II, action = "write",status = 'replace')
				do m = 1, Nop        
					write (203, '(1001(E19.12, tr3))') Energy_scat(m)        
				end do
			close(203)
		else
			open (unit = 203, file = file_name_output_II, action = "write",status = 'replace')
				do m = 1, Nop        
					write (203, '(1001(E19.12, tr3))') (real(E_out(m)%vector(n)), n= 1, 3), (imag(E_out(m)%vector(n)), n= 1, 3)         
				end do
			close(203)
		end if 			
		return
	end subroutine
	
end module lib_sie_quad_solver_normal