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
! several code lines are from AAS, RUNE ØISTEIN's master thesis
! “Electromagnetic Scattering: ” Norwegian University of Science and Technology, 2012
!  https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/246792
	
module lib_sie_tri_calculation_mod
	
	!use libmath
	use lib_sie_math	
	use lib_sie_constants
	use lib_sie_data_container
	use lib_sie_tri_data_container
	use lib_sie_integration_points
	
	implicit none
	
	private
	public :: cross_c
	public :: initialize_structure_tri
	public :: find_neighbours	
	public :: fn_parameter
	public :: fn_edge_center
	public :: Integration_scattered_field_tri
	public :: Integration_scattered_field_tri_within_medium
	public :: normal_integration_new	
	public :: quadrature_tri
	public :: vec_len
	public :: R_simplex
	public :: R_simplex_new
	public :: get_vector_b_tri
	public :: lib_sie_tri_get_impedance_edges
	public :: singular_integration_new
	public :: recursive_triangulation	
	public :: lib_sie_tri_input_field_calculation
	public :: precalculation_fn
	public :: data_import_tri
	public :: formulation_coefficient
	!test functions
	!public :: test_Iq_L	
	!public :: test_K3q
	!public :: test_K4q
	
	contains	
		
	subroutine data_import_tri()    !
      implicit none
      integer :: j,  p, m, q, s	
		character(len = 30) :: file_name
		
		allocate(ne_arr(number_objects))
		
		if (allocated(m_pairs_arr)) then
			deallocate(m_pairs_arr)
		end if
		
		allocate(m_pairs_arr(number_objects))
		allocate(np_arr(number_objects))
		
      select case (pre_types%object)			
		case('sphere')	
			do j = 1, number_objects
				call recursive_triangulation(struc_tri_spheres(j), tri_sp_parameters(j)%n_disc, tri_sp_parameters(j)%D*0.5)
				
				call initialize_structure_tri(struc_tri_spheres(j))
				
				m_pairs_arr(j) = size(struc_tri_spheres(j)%neighbours)
				ne_arr(j) = size(struc_tri_spheres(j)%elements)				
				np_arr(j) = size(struc_tri_spheres(j)%points)
			end do
			
			if (allocated(struc_tri%points))then
				deallocate(struc_tri%points)
			end if
			allocate(struc_tri%points(sum(np_arr)))	
				
			if (allocated(struc_tri%neighbours))then
				deallocate(struc_tri%neighbours)
			end if				
			allocate(struc_tri%neighbours(sum(m_pairs_arr)))	
				
			if (allocated(struc_tri%midpoint_edge))then
				deallocate(struc_tri%midpoint_edge)
			end if				
			allocate(struc_tri%midpoint_edge(sum(m_pairs_arr)))	
				
			if (allocated(struc_tri%elements))then
				deallocate(struc_tri%elements)
			end if				
			allocate(struc_tri%elements(sum(ne_arr)))	
				
			if (allocated(struc_tri%midpoint))then
				deallocate(struc_tri%midpoint)
			end if
			allocate(struc_tri%midpoint(sum(ne_arr)))	
				
			p = 0
			q = 0
			s = 0
			do j = 1, number_objects
				if (j .eq. 1) then
					do m = 1, m_pairs_arr(j)							
						struc_tri%neighbours(m)%corner = struc_tri_spheres(j)%neighbours(m)%corner														
						struc_tri%neighbours(m)%element = struc_tri_spheres(j)%neighbours(m)%element					
						struc_tri%midpoint_edge(m)%point = struc_tri_spheres(j)%midpoint_edge(m)%point + tri_sp_parameters(j)%centroid%point
					end do
				else
					p = p + m_pairs_arr(j-1)
					q = q + np_arr(j-1)
					s = s + ne_arr(j-1)
					do m = 1, m_pairs_arr(j)
						struc_tri%neighbours(m + p)%corner = struc_tri_spheres(j)%neighbours(m)%corner + q
						struc_tri%neighbours(m + p)%element = struc_tri_spheres(j)%neighbours(m)%element + s	
						struc_tri%midpoint_edge(m + p)%point = struc_tri_spheres(j)%midpoint_edge(m)%point + tri_sp_parameters(j)%centroid%point
					end do
				end if					
			end do
			
			p = 0	
			q = 0
			do j = 1, number_objects
				if (j .eq. 1) then
					do m = 1, ne_arr(j)
						allocate(struc_tri%elements(m)%corners(3))
						struc_tri%elements(m)%corners = struc_tri_spheres(j)%elements(m)%corners
						struc_tri%midpoint(m)%point = struc_tri_spheres(j)%midpoint(m)%point + tri_sp_parameters(j)%centroid%point	
					end do
				else
					p = p + ne_arr(j-1)						
					q = q + np_arr(j-1)
					do m = 1, ne_arr(j)
						allocate(struc_tri%elements(m + p)%corners(3))
						struc_tri%elements(m + p)%corners = struc_tri_spheres(j)%elements(m)%corners + q
						struc_tri%midpoint(m + p)%point = struc_tri_spheres(j)%midpoint(m)%point + tri_sp_parameters(j)%centroid%point	
					end do
				end if	
			end do
			!--------------
			p = 0	
			do j = 1, number_objects
				if (j .eq. 1) then
					do m = 1, tri_sp_parameters(j)%n_disc
						struc_tri%points(m)%point = struc_tri_spheres(j)%points(m)%point + tri_sp_parameters(j)%centroid%point
					end do
				else
					p = p + tri_sp_parameters(j-1)%n_disc
					do m = 1, tri_sp_parameters(j)%n_disc
						struc_tri%points(m + p)%point = struc_tri_spheres(j)%points(m)%point + tri_sp_parameters(j)%centroid%point
					end do
				end if	
			end do
			
			! save the meshed data from a sphere
			!file_name = 'tt_sphere_n258.txt'
			!!!!print*, 'file_name =', file_name
   !!
			!open (unit=206, file = file_name, action="write",status = 'replace')
			!	do j = 1, sum(ne_arr)
			!		write (206, '(201(i10, tr5))') struc_tri%elements(j)%corners
			!	end do
			!close(206)
   !
			!file_name = 'pp_sphere_n258.txt'
			!
			!open (unit=206, file = file_name, action="write",status = 'replace')			
			!	do j = 1, sum(np_arr)
			!		write (206, '(201(es19.12, tr5))') struc_tri%points(j)%point
			!	end do
			!close(206)	
			!call exit
			
		case ('surface')			
			do j = 1, number_objects
			
				call initialize_structure_tri(struc_tri_surfaces(j))	
				
				m_pairs_arr(j) = size(struc_tri_surfaces(j)%neighbours)
				ne_arr(j) = size(struc_tri_surfaces(j)%elements)
				np_arr(j) = size(struc_tri_surfaces(j)%points)
			end do			
			
			if (allocated(struc_tri%points))then
				deallocate(struc_tri%points)
			end if
			allocate(struc_tri%points(sum(np_arr)))	
				
			if (allocated(struc_tri%neighbours))then
				deallocate(struc_tri%neighbours)
			end if				
			allocate(struc_tri%neighbours(sum(m_pairs_arr)))	
				
			if (allocated(struc_tri%midpoint_edge))then
				deallocate(struc_tri%midpoint_edge)
			end if				
			allocate(struc_tri%midpoint_edge(sum(m_pairs_arr)))	
				
			if (allocated(struc_tri%elements))then
				deallocate(struc_tri%elements)
			end if				
			allocate(struc_tri%elements(sum(ne_arr)))	
				
			if (allocated(struc_tri%midpoint))then
				deallocate(struc_tri%midpoint)
			end if
			allocate(struc_tri%midpoint(sum(ne_arr)))	
				
			p = 0
			q = 0
			s = 0
			do j = 1, number_objects
				if (j .eq. 1) then
					!$omp parallel private (m) &  
					!$omp shared (struc_tri, struc_tri_surfaces, j) !
					!$omp do
					do m = 1, m_pairs_arr(j)							
						struc_tri%neighbours(m)%corner = struc_tri_surfaces(j)%neighbours(m)%corner														
						struc_tri%neighbours(m)%element = struc_tri_surfaces(j)%neighbours(m)%element					
						struc_tri%midpoint_edge(m)%point = struc_tri_surfaces(j)%midpoint_edge(m)%point !+ &
							!tri_sf_parameters(j)%centroid%point
					end do
					!$omp end do
					!$omp end parallel  
				else
					p = p + m_pairs_arr(j-1)
					q = q + np_arr(j-1)
					s = s + ne_arr(j-1)
					!$omp parallel private (m) &  
					!$omp shared (struc_tri, struc_tri_surfaces, tri_sf_parameters, j, p, q, s) !
					!$omp do
					do m = 1, m_pairs_arr(j)
						struc_tri%neighbours(m + p)%corner = struc_tri_surfaces(j)%neighbours(m)%corner + q
						struc_tri%neighbours(m + p)%element = struc_tri_surfaces(j)%neighbours(m)%element + s	
						struc_tri%midpoint_edge(m + p)%point = struc_tri_surfaces(j)%midpoint_edge(m)%point + &
							tri_sf_parameters(j)%centroid%point
					end do
					!$omp end do
					!$omp end parallel  					
				end if					
			end do
			
			p = 0	
			q = 0
			do j = 1, number_objects
				if (j .eq. 1) then
					!$omp parallel private (m) &  
					!$omp shared (struc_tri, ne_arr, struc_tri_surfaces, j) !
					!$omp do				
					do m = 1, ne_arr(j)
						allocate(struc_tri%elements(m)%corners(3))
						struc_tri%elements(m)%corners = struc_tri_surfaces(j)%elements(m)%corners
						struc_tri%midpoint(m)%point = struc_tri_surfaces(j)%midpoint(m)%point !+ &
							!tri_sf_parameters(j)%centroid%point					
					end do
					!$omp end do
					!$omp end parallel  							
				else
					p = p + ne_arr(j-1)						
					q = q + np_arr(j-1)
					!$omp parallel private (m) &  
					!$omp shared (struc_tri, struc_tri_surfaces, tri_sf_parameters, j, p, q) !
					!$omp do
					
					do m = 1, ne_arr(j)
						allocate(struc_tri%elements(m + p)%corners(3))
						struc_tri%elements(m + p)%corners = struc_tri_surfaces(j)%elements(m)%corners + q
						struc_tri%midpoint(m + p)%point = struc_tri_surfaces(j)%midpoint(m)%point + &
							tri_sf_parameters(j)%centroid%point
					end do
					!$omp end do
					!$omp end parallel  								
				end if	
			end do
			!--------------
			p = 0	
			do j = 1, number_objects			
				if (j .eq. 1) then
					!$omp parallel private (m) &  
					!$omp shared (struc_tri, struc_tri_surfaces, np_arr) !
					!$omp do								
					do m = 1, np_arr(j)
						struc_tri%points(m)%point = struc_tri_surfaces(j)%points(m)%point! + tri_sf_parameters(j)%centroid%point
					end do
					!$omp end do
					!$omp end parallel 					
				else
					p = p + ne_arr(j-1)
					!$omp parallel private (m) &  
					!$omp shared (j, p, struc_tri, tri_sf_parameters, struc_tri_surfaces, np_arr) !
					!$omp do
					do m = 1, np_arr(j)
						struc_tri%points(m + p)%point = struc_tri_surfaces(j)%points(m)%point + tri_sf_parameters(j)%centroid%point
					end do
					!$omp end do
					!$omp end parallel					
				end if	
			end do  				
				
		end select	
		deallocate(ne_arr, np_arr)
		return
	end subroutine data_import_tri	

	! when new structure is added, do not forget to change the range of capc_p(5)
	! since surface and sphere has different functions to calculate fn
	subroutine precalculation_fn(struct)
		implicit none
		integer :: m, ot, j, p		
		type(structure_tri) :: struct
		type(fn_rwg) :: fn
				
		m_pairs = size(struct%neighbours)		
		if (allocated(structure_edges)) then
			deallocate(structure_edges)
		end if
		
		allocate(structure_edges(m_pairs))
		if (calc_p(5) .eq. 1) then		
			p = 0
			do j = 1, number_objects
				if (j .eq. 1)then
					do m = 1, m_pairs_arr(j)
						do ot = 1, 2							
							!call fn_parameter_type_sphere(struct, m, ot, tri_sp_parameters(j)%centroid%point, fn)
							call fn_parameter_type(struct, m, ot, tri_sp_parameters(j)%centroid%point, fn)
							structure_edges(m)%fn_ot(ot)%fn = fn	
						end do			
					end do
				else 
					p = p + m_pairs_arr(j-1)
					do m = 1, m_pairs_arr(j)						
						do ot = 1, 2
							!call fn_parameter_type_sphere(struct, m+p, ot, tri_sp_parameters(j)%centroid%point, fn)					
							call fn_parameter_type(struct, m, ot, tri_sp_parameters(j)%centroid%point, fn)
							structure_edges(m+p)%fn_ot(ot)%fn = fn	
						end do			
					end do
				end if
			end do				
		else !if (calc_p(5) .eq. 2)then
			do m = 1, m_pairs
				do ot = 1, 2
					!tri_sf_parameters(1) is hard coded for test, has to be changed. 
					if (number_objects .eq. 1) then 
						!call fn_parameter_type_sphere(struct, m, ot, tri_sf_parameters(1)%centroid%point, fn) 
						call fn_parameter_type(struct, m, ot, tri_sf_parameters(1)%centroid%point, fn) 
					else 
						print*, 'Not implemented for multiple surfaces yet'
					end if
					structure_edges(m)%fn_ot(ot)%fn = fn	
				end do			
			end do				
		end if 		
		
	end subroutine 
	
	subroutine	formulation_coefficient(a_lc, d_lc)
		complex(dp), intent(out) :: a_lc(2), d_lc(2) ! for 
		real(dp) :: x, y
		x = 1.0
		y = 0.50
		
		a_lc(1) = 2*x*eta_1/(eta_1 + y*eta_2)
		a_lc(2) = 2*x*eta_2/(eta_1 + y*eta_2)
		d_lc(1) = (eta_1 + y*eta_2)/(2*x*eta_1)
		d_lc(2) = (eta_1 + y*eta_2)/(2*x*eta_2)
		
		!Based on MCTF
		
		!a_lc(1) = x*eta_1
		!a_lc(2) = y*eta_2
		!d_lc(1) = y*eta_2
		!d_lc(2) = x*eta_1
		
		return
	end subroutine formulation_coefficient
	
!!*********** 07.09.2022 for test. Combine the equation tegether.		
	subroutine get_vector_b_tri(vector_b)
		implicit none
      ! dummy		
		double complex, dimension(:), allocatable, intent(inout) :: vector_b
		complex(dp), dimension(4) :: sum_EH!
		complex(dp) :: a_lc(2), d_lc(2) ! for 
		integer :: m, ngp
	
		ngp = ng
		
		call formulation_coefficient(a_lc, d_lc)
		
		select case (pre_types%formulation)
			 case ('PMCHWT')
				do m = 1, m_pairs	
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = Sum_EH(1)
					vector_b(m + m_pairs) = Sum_EH(2)
				end do
			 case ('PMCHWT-S')
				do m = 1, m_pairs	
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = Sum_EH(1)
					vector_b(m + m_pairs) = Sum_EH(2)*eta_1
				end do
			 case ('MCTF')												
				do m = 1, m_pairs	
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = Sum_EH(1)
					vector_b(m + m_pairs) = Sum_EH(2)*eta_1*eta_2
				end do
			case ('JMCFIE')
				do m = 1, m_pairs	
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = 1/eta_1*Sum_EH(1) + Sum_EH(4)
					vector_b(m + m_pairs) = -Sum_EH(3) + eta_1*Sum_EH(2)
				end do		
			case ('MCTF2')! CTF was removed, instead MCTF2 was built
				do m = 1, m_pairs							
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = a_lc(1)/eta_1*Sum_EH(1)
					vector_b(m + m_pairs) = d_lc(1)*eta_1*Sum_EH(2)
				end do
			case ('ICTF')				
				do m = 1, m_pairs						
					call incident_field_tri_general(m, ngp, Sum_EH) !
					!eta_a = (eta_1 + eta_2)/2
					vector_b(m) = Sum_EH(1)/eta_a  
					vector_b(m + m_pairs) = Sum_EH(2)*eta_a					
				end do
			case ('CNF')
				do m = 1, m_pairs						
					call incident_field_tri_general(m, ngp, Sum_EH) !
					vector_b(m) = Sum_EH(4)
					vector_b(m + m_pairs) = -Sum_EH(3)
				end do		
			case('MNMF')			
				do m = 1, m_pairs						
					call incident_field_tri_general(m, ngp, Sum_EH) !					
					vector_b(m) = Sum_EH(4)*my_1/(my_1 + my_2)
					vector_b(m + m_pairs) = -Sum_EH(3)*eps_1/(eps_1 + eps_2)
				end do
		end select
!		
		!open (unit=206, file = 'vector_b_MCTF2.txt', action="write",status = 'replace')				
		!		do m = 1, m_pairs*2
		!			write (206, '(201(es19.12, tr5))') vector_b(m)
		!		end do
		!close(206)
	end subroutine get_vector_b_tri
		
		
	function lib_sie_tri_get_impedance_edges(m, n, r_ref)result(D_mat_tmp)
		use lib_sie_tri_data_container
		implicit none
		
		integer, intent(in) :: m, n
		real(dp), intent(in) :: r_ref
		complex(dp), dimension(2,2) :: D_mat_tmp		

		
		!dummy
		real(dp) :: r_ob 
		integer :: t, s, ngp, p, q
		complex(dp), dimension(:), allocatable :: zz
		character(len = 20) :: singular_type
				
		ngp = ng
		D_mat_tmp(1:2, 1:2) = (0.0, 0.0)
		select case(pre_types%formulation)
			case ('PMCHWT')
			do t = 1, 2
				do s = 1, 2!
					p = struc_tri%neighbours(m)%element(t)					
					q = struc_tri%neighbours(n)%element(s)	
					r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)
  
					if (r_ob .gt. 1.5*r_ref) then						
						singular_type = 'normal'							
						call tri_impedance_building_PMCHWT(m, n, t, s, ngp, zz, singular_type)
					else							
						singular_type = 'singular'
						call tri_impedance_building_PMCHWT(  m, n, t, s, ngp, zz, singular_type)
					end if						
					!--------------------------------------------------------
					D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
					D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
					D_mat_tmp(2, 1) = -D_mat_tmp(1, 2) 
					D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(3)
				end do
			end do
			case ('ICTF')
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)					
						q = struc_tri%neighbours(n)%element(s)	
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then									
							singular_type = 'normal'							
							call tri_impedance_building_ICTF(  m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'								
							call tri_impedance_building_ICTF(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = D_mat_tmp(2, 1) + zz(3)
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(4)
					end do
				end do
		!--------------------------
			case ('MCTF')
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)					
						q = struc_tri%neighbours(n)%element(s)	
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then									
							singular_type = 'normal'							
							call tri_impedance_building_PMCHWT(m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'								
							call tri_impedance_building_PMCHWT(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = -D_mat_tmp(1, 2)*eta_1*eta_2
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(3)*eta_1*eta_2
					end do
				end do
				
			case ('PMCHWT-S') !not complete
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)
						q = struc_tri%neighbours(n)%element(s)
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then									
							singular_type = 'normal'							
							call tri_impedance_building_PMCHWT(  m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'								
							call tri_impedance_building_PMCHWT(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = -D_mat_tmp(2, 1) + zz(2)
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(3)
					end do
				end do	
				D_mat_tmp(1, 2) = D_mat_tmp(1, 2)*eta_1
				D_mat_tmp(2, 1) = D_mat_tmp(2, 1)*eta_1
				D_mat_tmp(2, 2) = D_mat_tmp(2, 2)*eta_1**2
		!---------------------
			case ('JMCFIE')
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)					
						q = struc_tri%neighbours(n)%element(s)	
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then	
							singular_type = 'normal'							
							call tri_impedance_building_JMCFIE(  m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'
							call tri_impedance_building_JMCFIE(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = D_mat_tmp(2, 1) + zz(3)
						D_mat_tmp(2, 2) = D_mat_tmp(1, 1)						
					end do
				end do
					
			case ('MCTF2') ! CTF was removed, instead MCTF2 was built
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)
						q = struc_tri%neighbours(n)%element(s)
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then	
							singular_type = 'normal'							
							call tri_impedance_building_MCTF2(  m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'
							call tri_impedance_building_MCTF2(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = D_mat_tmp(2, 1) + zz(3)
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(4)						
					end do
				end do	
					
			case ('CNF')
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)					
						q = struc_tri%neighbours(n)%element(s)	
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then	
							singular_type = 'normal'							
							call tri_impedance_building_CNF(m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'
							call tri_impedance_building_CNF(m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = D_mat_tmp(2, 1) + zz(3)
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(4)
					end do
				end do
			case ('MNMF')
				do t = 1, 2
					do s = 1, 2!
						p = struc_tri%neighbours(m)%element(t)					
						q = struc_tri%neighbours(n)%element(s)	
						r_ob = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)  
						if (r_ob .gt. 1.5*r_ref) then
							singular_type = 'normal'							
							call tri_impedance_building_MNMF(  m, n, t, s, ngp, zz, singular_type)							
						else							
							singular_type = 'singular'
							call tri_impedance_building_MNMF(  m, n, t, s, ngp, zz, singular_type)	
						end if
						D_mat_tmp(1, 1) = D_mat_tmp(1, 1) + zz(1)
						D_mat_tmp(1, 2) = D_mat_tmp(1, 2) + zz(2)
						D_mat_tmp(2, 1) = D_mat_tmp(2, 1) + zz(3)
						D_mat_tmp(2, 2) = D_mat_tmp(2, 2) + zz(4)
					end do
				end do
			end select
		end function
	
	subroutine recursive_triangulation(struct, n_discr, R_factor)
     ! Creating symmetric triangular discretization recursively
     ! n=6,18,58,258,1026, ...
      type(structure_tri), intent(inout) :: struct
      integer :: n_discr, i, ne, np_curr, ne_curr, lim!, split_num
      real(dp) :: R_factor

      print *, 'Creating discretization recursively with #points = ', n_discr
		if (allocated(struct%points))then
			deallocate(struct%points)
		end if
      allocate(struct%points(n_discr))
		
      struct%points(1)%point=(/ 1.0, 0.0, 0.0 /)*R_factor
      struct%points(2)%point=(/-1.0, 0.0, 0.0 /)*R_factor
      struct%points(3)%point=(/ 0.0, 1.0, 0.0 /)*R_factor
      struct%points(4)%point=(/ 0.0,-1.0, 0.0 /)*R_factor
      struct%points(5)%point=(/ 0.0, 0.0, 1.0 /)*R_factor
      struct%points(6)%point=(/ 0.0, 0.0,-1.0 /)*R_factor

      ne=2*(n_discr-2) !int(8*4**(split_num-1))
		
		if (allocated(struct%elements))then
			deallocate(struct%elements)
		end if 
      allocate(struct%elements(ne))
      do i=1,ne		
			allocate(struct%elements(i)%corners(3))
      end do
      struct%elements(1)%corners=(/ 1, 5, 3/)
      struct%elements(2)%corners=(/ 3, 5, 2/)
      struct%elements(3)%corners=(/ 2, 5, 4/)
      struct%elements(4)%corners=(/ 4, 5, 1/)
      struct%elements(5)%corners=(/ 1, 3, 6/)
      struct%elements(6)%corners=(/ 3, 2, 6/)
      struct%elements(7)%corners=(/ 2, 4, 6/)
      struct%elements(8)%corners=(/ 4, 1, 6/)

      np_curr=6
      ne_curr=8

      lim=int(log10(ne/6.0)/log10(4.0))
      do i=1, lim
			call recursive_split(struct, np_curr, ne_curr, R_factor)
			ne_curr=ne_curr*4
			np_curr=(ne_curr/2)+2
      end do
		
		do i = 1, size(struct%points)
			struct%points(i)%point(1) = struct%points(i)%point(1)
			struct%points(i)%point(2) = struct%points(i)%point(2)
			struct%points(i)%point(3) = struct%points(i)%point(3)
		end do
		
		return
	end subroutine
	
	subroutine initialize_structure_tri(struct)
		type(structure_tri), intent(inout) :: struct
				
		call find_neighbours(struct)
      call find_midpoint_element(struct)    
		call find_midpoint_edge(struct)
		
	end subroutine

	subroutine recursive_split(struct, n, ne, Rf)
		integer :: n, ne, i, counter, &
            ind_a, ind_b, ind_c
      type(structure_tri) :: struct
      type(Element_tri), dimension(:), allocatable :: el_list
      real(dp), intent(in) :: Rf !Scaling factor
      real(dp) :: a(3),b(3),c(3)
      logical :: lo(3)
      				
		allocate(el_list(4*ne))
      counter=n+1
		
      do i=1,ne
         lo=(/ .false., .false., .false. /)
         a=normalize(0.5*(struct%points(struct%elements(i)%corners(1))&
         %point + struct%points(struct%elements(i)%corners(3))%point ))*Rf
         b=normalize(0.5*(struct%points(struct%elements(i)%corners(1))&
         %point+struct%points(struct%elements(i)%corners(2))%point ))*Rf
         c=normalize(0.5*(struct%points(struct%elements(i)%corners(2))&
         %point+struct%points(struct%elements(i)%corners(3))%point ))*Rf
         ind_a=point_in_list(a,struct%points)
         lo(1)=ind_a.eq. -1
         if (lo(1))then ! a is not present in list
               struct%points(counter)%point=a
               ind_a=counter
               counter=counter+1
         end if

         ind_b=point_in_list(b,struct%points)
         lo(1)=ind_b.eq. -1
         if (lo(1))then ! a is not present in list
               struct%points(counter)%point=b
               ind_b=counter
               counter=counter+1
         end if

         ind_c=point_in_list(c,struct%points)
         lo(1)=ind_c.eq. -1
         if (lo(1))then ! a is not present in list
               struct%points(counter)%point=c
               ind_c=counter
               counter=counter+1
         end if
         struct%elements(ne+(i-1)*3+1)%corners=&
         (/ struct%elements(i)%corners(1), ind_b, ind_a/)
         struct%elements(ne+(i-1)*3+2)%corners=&
         (/ ind_b, struct%elements(i)%corners(2), ind_c/)
         struct%elements(ne+(i-1)*3+3)%corners=&
         (/ ind_a, ind_c, struct%elements(i)%corners(3)/)
         struct%elements(i)%corners=(/ind_a, ind_b, ind_c/)
      end do
    end subroutine

    function point_in_list(p, list)
     ! Function returning the index of the point p in list
     ! and returns -1 if it is not in the list
        type(Point), dimension(:), allocatable, intent(in) :: list
        real(dp) :: temp(3)
        real(dp), intent(in) :: p(3)
        integer :: point_in_list
        integer :: i
		  
        point_in_list=-1
        do i=1,size(list)
            temp=list(i)%point-p
            temp(1)=sqrt(temp(1)**2+temp(2)**2+temp(3)**2)
            if (temp(1)<1e-30)then
                point_in_list=i
            return
            end if
        end do
        return
	end function
  
	subroutine find_midpoint_element(struct)
		 ! Subroutine finding all centroids of the elements in struct
		 ! and saving the resulting points in struct%midpoints.
		 ! Compatible with only triangular elements.
		 integer :: ne, i
		 type(structure_tri), intent(inout) :: struct
    
		 real(dp) :: v1(3),v2(3),v3(3), mid(3)

		 ne=size(struct%elements)
		 if (allocated(struct%midpoint))then
			deallocate(struct%midpoint)
		 endif
		 
		 allocate(struct%midpoint(ne))
		 do i=1,ne
		 v1=struct%points(struct%elements(i)%corners(1))%point
		 v2=struct%points(struct%elements(i)%corners(2))%point
		 v3=struct%points(struct%elements(i)%corners(3))%point
		 mid=(v1+v2+v3)/3.0
		 struct%midpoint(i)%point=mid
		 end do
    end subroutine
	 
	 subroutine fn_edge_center(struct, m, ot, p_edgec)
		type(structure_tri), intent(in) :: struct
		integer, intent(in) :: m, ot
		real(dp), intent(out) :: p_edgec(3)
		real(dp), dimension(3) :: p1, p2, pc
		
		integer :: j, nv, p
		type(Pair) :: pair_p
		
      pair_p = struct%neighbours(m)
      do j = 1, 3                
         if (struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(1) &
               .and. struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(2))then
               pc(1:3)=struct%points(struct%elements(pair_p%element(ot))&
               %corners(j))%point ! vertex of the triangle
			nv = j
         end if
      end do
		p = struct%neighbours(m)%element(ot) !element index
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
			p_edgec = (p1 + p2)/2 !midpoint of the common edge
      return		
	end subroutine 
	
	subroutine find_midpoint_edge(struct)
		 ! Subroutine finding all centroids of the elements in struct
		 ! and saving the resulting points in struct%midpoints.
		 ! Compatible with only triangular elements.
		 integer :: n_pair, m, ot
		 type(structure_tri), intent(inout) :: struct
		 
		 n_pair=size(struct%neighbours)
		 if (allocated(struct%midpoint_edge))then
			deallocate(struct%midpoint_edge)
		 endif
		 ot = 1
		 allocate(struct%midpoint_edge(n_pair))
		 do m = 1, n_pair
			call fn_edge_center(struct, m, ot, struct%midpoint_edge(m)%point)
		 end do
    end subroutine

	 ! The proper name should be find pairs
	 ! The total number of neighbours is the number of pairs
    subroutine find_neighbours(struct)
        ! Subroutine finding all neighbours in struct and
        ! saving the result in struct%neighbours.
        ! Compatible with arbitrary number of corners.
        integer :: ne, i, j, d1,d2, nn, &
        out_arr(2)
        type(Pair), dimension(:), allocatable :: temp_neighbours
        type(structure_tri), intent(inout) :: struct
        type(structure_tri) :: temp_struct

        nn=0
        ne=size(struct%elements)
        if (allocated(struct%neighbours))then
            deallocate(struct%neighbours)
        endif

        if (allocated(temp_struct%neighbours))then
            deallocate(temp_struct%neighbours)
        endif

        allocate(temp_neighbours(3*ne/2))
        allocate(temp_struct%neighbours(3*ne/2))

        do i = 1, ne
            do j = i+1,ne 
                d1 = size(struct%elements(i)%corners)
                d2 = size(struct%elements(j)%corners)
                out_arr = common_2(struct%elements(i)%corners, &
                        struct%elements(j)%corners, d1, d2 )
                if ( out_arr(2) /= -1 )then
                    nn = nn + 1
                    temp_struct%neighbours(nn)%element=(/ i, j /)
                    temp_struct%neighbours(nn)%corner=(/out_arr(1:2)/)
                end if
            end do
        end do
        allocate(struct%neighbours(nn))
        struct%neighbours(1:nn) = temp_struct%neighbours(1:nn)
        deallocate(temp_struct%neighbours)
        return
    end subroutine

    function common_2(arr1, arr2, dim1, dim2)
    ! function checking whether arr1 and arr2 have 2 common elements.
    ! If true, the function returns the two common elements
    ! If false, the function returns (-1, -1)
    ! Compatible with arbitrary number of dimensions
    integer, intent(in) :: dim1, dim2
    integer :: i, j, n_equal, common_2(2), n_different
    integer, intent(in) :: arr1(dim1), arr2(dim2)
    common_2=(/-1,-1/)

    n_equal=0
    do i=1,dim1
        n_different=0
        do j=1,dim2
            if (arr1(i) == arr2 (j))then
                n_equal=n_equal+1
                common_2(n_equal)=arr1(i)
            endif
        end do
    end do
    return
    end function common_2

	function normalize(v)
		real(dp), intent(in) :: v(3)
		real(dp) :: normalize(3), l		
		l = sqrt(v(1)**2+v(2)**2+v(3)**2)
		normalize(1:3)=v(1:3)/l		
		return
	end function normalize
			
	!Including V_nH, tangential and normal components
	subroutine incident_field_tri_general(m, ngp, Sum_EH)
		complex(dp), intent(out) ::  Sum_EH(4)
		integer, intent(in) :: ngp, m 
		
		type(fn_rwg) :: fm
		integer :: i, ot
		real(dp) :: R_a(3), r_fm(3)!
		real(dp), dimension(100) :: a, b, w
		complex(dp) :: k_hat(3)
		type(vector_c) :: Ec_in, Hc_in
		
		
		k_hat = illumination_p%k_in*(1 + im*0)	
		
		call Quadrature_tri(ngp, a, b, w)

		Sum_EH(1:4) = (/(0.0, 0.0),  (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		if (pre_types%illumination == 'Gaussian')  then
			do ot = 1, 2	
				fm = structure_edges(m)%fn_ot(ot)%fn
				do i = 1, ngp				
					R_a = r_simplex_new(a(i), b(i), fm%corners)
					r_fm = (R_a - fm%corners(3)%point)*fm%df*0.5
					call Gaussian_beam(beam_waist, R_a, illumination_p, Ec_in, k1)
					Hc_in%vector = cross_c(k_hat, Ec_in%vector)/(my_1*c0)
					Sum_EH(1) = Sum_EH(1) + dot_product(r_fm, Ec_in%vector)*w(i)*fm%area!
					Sum_EH(2) = Sum_EH(2) + dot_product(r_fm, Hc_in%vector)*w(i)*fm%area!
					Sum_EH(3) = Sum_EH(3) + dot_product(r_fm, cross_rc(fm%nq, Ec_in%vector))*w(i)*fm%area!
					Sum_EH(4) = Sum_EH(4) + dot_product(r_fm, cross_rc(fm%nq, Hc_in%vector))*w(i)*fm%area!
				end do !loop i		
			end do
		else if ((pre_types%illumination == 'Plane') .or. (pre_types%illumination == 'Conical_plane')) then
			do ot = 1, 2	
				fm = structure_edges(m)%fn_ot(ot)%fn
				!print*, 'illumination_p%k_in=', illumination_p%k_in
				do i = 1, ngp				
					R_a = r_simplex_new(a(i), b(i), fm%corners)
					r_fm = (R_a - fm%corners(3)%point)*fm%df*0.5
					Ec_in%vector = illumination_p%E_in*exp(-im*dot_product(k1*illumination_p%k_in, R_a))	
					Hc_in%vector = cross_c(k_hat, Ec_in%vector)/(my_1*c0)
					Sum_EH(1) = Sum_EH(1) + dot_product(r_fm, Ec_in%vector)*w(i)*fm%area!
					Sum_EH(2) = Sum_EH(2) + dot_product(r_fm, Hc_in%vector)*w(i)*fm%area!
					Sum_EH(3) = Sum_EH(3) + dot_product(r_fm, cross_rc(fm%nq, Ec_in%vector))*w(i)*fm%area!
					Sum_EH(4) = Sum_EH(4) + dot_product(r_fm, cross_rc(fm%nq, Hc_in%vector))*w(i)*fm%area!
				end do !loop i		
			end do
			
		else if (pre_types%illumination == 'Conical') then
			do ot = 1, 2	
				fm = structure_edges(m)%fn_ot(ot)%fn
				do i = 1, ngp				
					R_a = r_simplex_new(a(i), b(i), fm%corners)
					r_fm = (R_a - fm%corners(3)%point)*fm%df*0.5				
					call conical_illumination(p_obj, R_a, illumination_p, Ec_in)				
					Hc_in%vector = cross_c(k_hat, Ec_in%vector)/(my_1*c0)
					Sum_EH(1) = Sum_EH(1) + dot_product(r_fm, Ec_in%vector)*w(i)*fm%area!
					Sum_EH(2) = Sum_EH(2) + dot_product(r_fm, Hc_in%vector)*w(i)*fm%area!
					Sum_EH(3) = Sum_EH(3) + dot_product(r_fm, cross_rc(fm%nq, Ec_in%vector))*w(i)*fm%area!
					Sum_EH(4) = Sum_EH(4) + dot_product(r_fm, cross_rc(fm%nq, Hc_in%vector))*w(i)*fm%area!
				end do !loop i		
			end do
		end if
		return
	end subroutine incident_field_tri_general
		
	
	! Normal vector calculation, the original one	
	! The different lies in the propagation direciotn of k0 in calc_p(7)
	function norm_vec_r_new(corner_fn, mid, center)!result(norm_vec_r)
	 ! Function returning the normal vector of the surface 
	 ! (zero imaginary part) of a discretization element
	 ! pointing outwards from the origin when p1,p2,p3 are
	 ! corners of the element.
		type(point), intent(in) :: corner_fn(3)
		real(dp), intent(in) :: mid(3), center(3)		
		real(dp) :: norm_vec_r_new(3)
		real(dp) :: v1(3), v2(3), vc(3)		
		integer :: i
		if (calc_p(7) .eq. 1)then
			v1 = corner_fn(1)%point-corner_fn(2)%point
			v2 = corner_fn(3)%point-corner_fn(2)%point
		else 
			v1 = corner_fn(3)%point-corner_fn(1)%point
			v2 = corner_fn(2)%point-corner_fn(1)%point	
		end if
		norm_vec_r_new = cross_r(v1,v2)/vec_len(cross_r(v1,v2))		
		return
	end function
	
	!function norm_vec_r(corner_fn, mid, center)!result(norm_vec_r)
	! ! Function returning the normal vector
	! ! (zero imaginary part) of a discretization element
	! ! pointing outwards from the origin when p1,p2,p3 are
	! ! corners of the element.
	!	type(point), intent(in) :: corner_fn(3)
	!	real(dp), intent(in) :: mid(3), center(3)
	!	
	!	real(dp) :: norm_vec_r(3)
	!	real(dp) :: v1(3), v2(3), vc(3)
	!	
	!	integer :: i
	!	
	!	v1 = corner_fn(3)%point-corner_fn(1)%point
	!	v2 = corner_fn(2)%point-corner_fn(1)%point		
	!	norm_vec_r = cross_r(v1,v2)/vec_len(cross_r(v1,v2))
	!	
	!	return
	!end function
	 
	subroutine m_vector(n_vec, corner_fn, m_v)
        real(dp), dimension(3), intent(in) :: n_vec!p1, p2, p3, 
		  type(point), intent(in) :: corner_fn(3)
        real(dp), intent(out) :: m_v(3, 3)		  
        m_v(1, :) = cross_r(n_vec, corner_fn(3)%point-corner_fn(2)%point) &
			/vec_len(cross_r(n_vec, corner_fn(3)%point-corner_fn(2)%point))
        m_v(2, :) = cross_r(n_vec, corner_fn(1)%point-corner_fn(3)%point) &
		   /vec_len(cross_r(n_vec, corner_fn(1)%point-corner_fn(3)%point))
        m_v(3, :) = cross_r(n_vec, corner_fn(2)%point-corner_fn(1)%point) &
		   /vec_len(cross_r(n_vec, corner_fn(2)%point-corner_fn(1)%point))        
        return
    end subroutine 
	
	subroutine test_Iq_L()
		implicit none
		real(dp) :: p1(3), p2(3), Ra(3), ll
		integer :: q 
		q = -1
		p1 = (/1.0, 2.0, 3.0/)
		p2 = (/-1.0, -2.0, -3.0/)
		Ra = (/0.1, 0.2, -0.3/) 		
		ll = Iq_L(q, Ra, p1, p2)
		print*, 'll =', ll
	end subroutine
	
	!Equation (33) in I. Hänninen's paper, PIER 63, 243-278 (2006)
	function Iq_L(q, Ra, p1, p2)
      real(dp), intent(in) :: Ra(3), p1(3), p2(3)
		integer, intent(in) :: q   
      real(dp) :: IL_n1, IL_p1, IL_p3, IL_p5, Iq_L, &
			Sum_p, Sum_n, Sp, Sn, Rp, Rn, S_vec(3) !negative and positive for n and p
      real(dp) :: Dif_p, Dif_n, R0      
 
		s_vec = (p2-p1)/vec_len(p2-p1)		
		Sp = dot_product((p2 - Ra), s_vec)
		Sn = dot_product((p1 - Ra), s_vec)		
		Rn = vec_len(Ra - p1)
		Rp = vec_len(Ra - p2)
		
		Sum_p = Rp + Sp !eq. (41-43) in Hänninen
		Sum_n = Rn + Sn
		Dif_p = Rp - Sp
		Dif_n = Rn - Sn
		R0 = sqrt(Rp**2 -Sp**2)
                        
		if (R0 .eq. 0.0) then
				IL_n1 = log(Sp/Sn)
		else  
			if (abs(Dif_p) > abs(Sum_n)) then
					IL_n1 = log(Dif_n/Dif_p)
			else
					IL_n1 = log(Sum_p/Sum_n)
			end if
		end if             
        
		IL_p1 = 0.5*R0**2*IL_n1 + 0.5*(Sp*Rp-Sn*Rn) !q = 1
		IL_p3 = 0.75*R0**2*IL_p1 + 0.25*(Sp*Rp**3-Sn*Rn**3) !q = 3
		IL_p5 = 5.0/6.0*R0**2*IL_p3 + 1.0/6.0*(Sp*Rp**5-Sn*Rn**5) !q = 5
		
		if (q == -1) then
			Iq_L = IL_n1
		else if (q == 1)then
			Iq_L = IL_p1
		else if (q == 3)then
			Iq_L = IL_p3
		else if (q == 5)then
			Iq_L = IL_p5
		else 
		print*,	'q is not in the given range!'
		end if
		return   
	end function Iq_L
	
	subroutine test_K3q(struct)
		implicit none
		type(Structure_tri), intent(in) :: struct		
		type(fn_rwg) :: fm, fn
		integer :: m, n, ngp1, ngp2, tt, ll		
		real(dp), dimension(20) :: w1, a1, b1
		real(dp) :: w2(20), a2(20), b2(20), TRF
		real(dp) :: Rq(3), Ra(3), r_fn(3)
		real(dp) :: R(3), I_K3q(3)
		integer :: i, j,  q
		real(dp) :: sum_a(3)
		logical :: Iqs_L
		
		Iqs_L = .false.		
		sum_a = (/0.0, 0.0, 0.0/)
		
		m = 1
		tt = 2
		n=4
		ll = 1

		ngp1 = 6
		call Quadrature_tri(ngp1, a1, b1, w1)
		
		ngp2 = 12
		call Quadrature_tri(ngp2, a2, b2, w2) 
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
				
		do i = 3, 3
			Ra = r_simplex_new(a1(i), b1(i), fm%corners)			
			do j = 1, ngp2 
				Rq = r_simplex_new(a2(j), b2(j), fn%corners)
				r_fn = 0.5*(Rq-fn%corners(3)%point)*fn%df
				R = Ra - Rq
				TRF = w2(j)*fn%area	
				sum_a = sum_a + 3*TRF !q = 1
			end do !loop j
		end do  !loop i
		print*, 'sum_a', sum_a
		q = 1
		I_K3q = K3q(q, Ra, fn, Iqs_L)
		print*, 'I_k3q=', I_K3q		
	end subroutine
	
	subroutine test_K2q(struct)
		implicit none
		type(Structure_tri), intent(in) :: struct		
		type(fn_rwg) :: fm, fn
		integer :: m, n, ngp1, ngp2, tt, ll		
		real(dp), dimension(20) :: w1, a1, b1
		real(dp) :: w2(20), a2(20), b2(20), TRF
		real(dp) :: Rq(3), Ra(3), r_fn(3)
		real(dp) :: R(3), I_K2q(3)
		integer :: i, j,  q
		real(dp) :: sum_a(3)
		logical :: Iqs_L
		
		Iqs_L = .false.		
		sum_a = (/0.0, 0.0, 0.0/)
		
		m = 1
		tt = 2
		n=4
		ll = 1

		ngp1 = 6
		call Quadrature_tri(ngp1, a1, b1, w1)
		
		ngp2 = 12
		call Quadrature_tri(ngp2, a2, b2, w2) 		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
				
		do i = 3, 3
			Ra = r_simplex_new(a1(i), b1(i), fm%corners)			
			do j = 1, ngp2 
				Rq = r_simplex_new(a2(j), b2(j), fn%corners)
				r_fn = 0.5*(Rq-fn%corners(3)%point)*fn%df
				R = Ra - Rq
				TRF = w2(j)*fn%area	
				sum_a = sum_a + vec_len(R)*r_fn*TRF !q = 1
			end do !loop j
		end do  !loop i
		print*, 'sum_a', sum_a
		q = 1
		I_K2q = K2q(q, Ra, fn, Iqs_L)
		print*, 'I_k2q=', I_K2q		
	end subroutine
	
	subroutine test_K4q(struct)
		implicit none
		type(Structure_tri), intent(in) :: struct		
		type(fn_rwg) :: fm, fn
		integer :: m, n, ngp1, ngp2, tt, ll		
		real(dp), dimension(20) :: w1, a1, b1
		real(dp) :: w2(20), a2(20), b2(20), TRF
		real(dp) :: Rq(3), Ra(3), r_fn(3)
		real(dp) :: Rn(3), R(3), I_K4q(3)
		integer :: i, j,  q
		real(dp) :: sum_a(3)
		logical :: Iqs_L
		
		Iqs_L = .false.		
		sum_a = (/0.0, 0.0, 0.0/)
		
		m = 1
		tt = 2
		n = 12
		ll = 1

		ngp1 = 6
		call Quadrature_tri(ngp1, a1, b1, w1)		
		ngp2 = 12
		call Quadrature_tri(ngp2, a2, b2, w2) 	
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
				
		do i = 3, 3
			Ra = r_simplex_new(a1(i), b1(i), fm%corners)			
			do j = 1, ngp2 
				Rq = r_simplex_new(a2(j), b2(j), fn%corners)
				r_fn = 0.5*(Rq-fn%corners(3)%point)*fn%df				
				R = Ra - Rq				
				Rn = R/vec_len(R)
				TRF = w2(j)*fn%area				
				sum_a = sum_a - cross_r(Rn, r_fn)*TRF !q = 1, very good
			end do !loop j
		end do  !loop i
		print*, 'sum_a', sum_a
		q = 1
		I_K4q = K4q(q, Ra, fn, Iqs_L)
		print*, 'I_k4q=', I_K4q		
	end subroutine
	
	function Iq_s(q, R_a, fn, Iqs_l)
		real(dp), dimension(3), intent(in) :: R_a
		type(fn_rwg), intent(in) :: fn
		integer, intent(in) :: q		
		real(dp) :: Iq_s, m_r(3, 3), r(3), nul
		real(dp) :: Omeg, x, y, sum_IL, Is_n3, Is_n1, Is_p1, h  !rp is the r' in Eq.(45). Hänninen's paper
		real(dp), dimension(3, 3) :: arr
		real(dp), dimension(3) :: t_v, IL !n_q,
		integer :: i, n
		logical :: Iqs_l
		
		type(edge_tri), dimension(3) :: edge_fn
		
		nul = 0.0
		do i = 1, 3
			arr(i, :) = (fn%corners(i)%point - R_a)/vec_len(fn%corners(i)%point - R_a)
      end do
		
		Is_n3 = nul
		x = 1 + dot_product(arr(1,:), arr(2,:)) + dot_product(arr(1,:), arr(3,:)) + dot_product(arr(2,:), arr(3,:))
		r = cross_r(arr(2,:), arr(3,:))
		y = abs(dot_product(arr(1,:), r))
		
		Omeg = 2*atan2(y, x)				
		h = dot_product(fn%nq, (R_a - fn%mid))
		
		if (h == nul) then
			Is_n3 = nul
		else 
			if (Iqs_l .eqv. .true.) then
				Is_n3 = nul
			else if ((Omeg <= PI) .and. (Omeg > -PI)) then
				Is_n3 = -1/h*Omeg
			else
				print*, 'Omega is outside the range'			
				print*, 'Ra =', R_a
				print*, 'p1 =', fn%corners(1)%point
				print*, 'p2 =', fn%corners(2)%point		
				print*, 'p3 =', fn%corners(3)%point
			end if	
		end if
		
		call m_vector(fn%nq, fn%corners, m_r)
		call edge_vectors_assigned(fn%corners, edge_fn)
		n = -1
		sum_IL = 0.0
		do i = 1, 3
			t_v(i) = dot_product(m_r(i, :), R_a-(edge_fn(i)%p_t + edge_fn(i)%p_s)/2)  !
			IL(i) = Iq_L(n, R_a, edge_fn(i)%p_t, edge_fn(i)%p_s)
			sum_IL = sum_IL + t_v(i)*IL(i)
		end do		
		Is_n1 = 1/(real(n)+2)*(n*h**2*Is_n3 - sum_IL) 
		
		n = 1
		sum_IL = 0.0
		do i = 1, 3
			IL(i) = Iq_L(n, R_a, edge_fn(i)%p_t, edge_fn(i)%p_s)		        
			sum_IL = sum_IL + t_v(i)*IL(i)
		end do 
		Is_p1 = 1.0/(real(n)+2.0)*(n*h**2*Is_n1 - sum_IL )
  
		if ((Iqs_l .eqv. .true.) .and. (q > -1)) then
			Iq_s = 0.0
		else 
			if (q ==-3) then 
				Iq_s = Is_n3
			else if (q==-1) then
				Iq_s = Is_n1
			else if (q == 1) then
				Iq_s = Is_p1
			else 
				print*, 'q is beyond the given range!'
			end if
		end if
        return    
    end function Iq_s
	
	function K2q(q, R_a, fn, Iqs_L)
      real(dp), intent(in) :: R_a(3)
		type(fn_rwg), intent(in) :: fn
      integer, intent(in) :: q		
      logical, intent(in) :: Iqs_L	
		real(dp) :: K2q(3)
		
		!dummy
      integer :: i   				
      real(dp), dimension(3) :: rho, IL_arr!, n_q
      real(dp) :: Iqs_q, sum_IL(3), h, m_v(3, 3)
		type(edge_tri), dimension(3) :: edge_fn
		
		K2q(1:3) = 0.0
		
		h = dot_product(fn%nq, (R_a - fn%mid))
		call m_vector(fn%nq, fn%corners, m_v)
		
		call edge_vectors_assigned(fn%corners, edge_fn) 		
		do i = 1, 3
			IL_arr(i) = Iq_L(q+2, R_a, edge_fn(i)%p_s, edge_fn(i)%p_t)
		end do
		
		Iqs_q = Iq_s(q, R_a, fn, Iqs_L)		
		sum_IL = 0.0
		do i = 1, 3
			sum_IL = sum_IL + IL_arr(i)*m_v(i, :)
		end do 
		rho = R_a - h*fn%nq!
		K2q = (1/(real(q)+2)*sum_IL + (rho - fn%corners(3)%point)*Iqs_q)*fn%df*0.5
		return
    end function K2q
    !
    function K4q(q, R_a, fn, Iqs_L)
      real(dp), intent(in) :: R_a(3)
		type(fn_rwg), intent(in) :: fn
      logical, intent(in) :: Iqs_L
		
      real(dp) :: I_arr(3), K4q(3), K3q(3)!n_q(3),
		
		!dummy
		type(edge_tri), dimension(3) :: edge_fn  
		real(dp) :: Is, h, m_v(3, 3), sum_tmp(3)
      integer :: i, q
		
		sum_tmp = (/0.0, 0.0, 0.0/)
		
		if (Iqs_L) then
			K4q = (/0.0, 0.0, 0.0/)
		else 
			call m_vector(fn%nq, fn%corners, m_v)
			call edge_vectors_assigned(fn%corners, edge_fn) 		
			do i = 1, 3
				I_arr(i) = Iq_L(q, R_a, edge_fn(i)%p_s, edge_fn(i)%p_t)				
				sum_tmp(:) = sum_tmp(:) +  m_v(i, :)*I_arr(i) 
			end do
			h = dot_product(fn%nq, (R_a - fn%mid))
			Is = Iq_s(q-2, R_a, fn, Iqs_L)
			K3q = sum_tmp - h*real(q)*fn%nq*Is!
			K4q = -cross_r((R_a-fn%corners(3)%point), K3q)*fn%df*0.5!
		end if
		return		
	end function K4q
	
	function K3q(q, R_a, fn, Iqs_L)
      real(dp), intent(in) :: R_a(3)
		type(fn_rwg), intent(in) :: fn
      logical, intent(in) :: Iqs_L		
      real(dp) :: I_arr(3), K3q(3)!n_q(3), 
		
		!dummy
		type(edge_tri), dimension(3) :: edge_fn  
		real(dp) :: Is, h, m_v(3, 3), sum_tmp(3)
      integer :: i, q
		
		sum_tmp = (/0.0, 0.0, 0.0/)
		call m_vector(fn%nq, fn%corners, m_v)
		call edge_vectors_assigned(fn%corners, edge_fn) 		
		do i = 1, 3
			I_arr(i) = Iq_L(q, R_a, edge_fn(i)%p_s, edge_fn(i)%p_t)				
			sum_tmp(:) = sum_tmp(:) +  m_v(i, :)*I_arr(i) 
		end do
		h = dot_product(fn%nq, (R_a - fn%mid))
		Is = Iq_s(q-2, R_a, fn, Iqs_L)
		K3q = sum_tmp - h*real(q)*fn%nq*Is!		
      return	
	end function K3q
    !
   subroutine Diver_Green(km, R_norm, Greenf, DGf)
      implicit none
      complex(dp), intent(in) :: km
      real(dp), intent(in) :: R_norm
      complex(dp), intent(out) :: DGf, Greenf
      complex(dp) :: Exp_f
        
      Exp_f = exp(-im*km*R_norm) ! It should be norm(R). m: medium
      DGf = -(1.0 + im*km*R_norm)*Exp_f/(4*PI*R_norm**2)!r' by Mitschang, without the unit vector
      Greenf = Exp_f/(4*PI*R_norm) !
      return 
   end subroutine Diver_Green
    
   !Smoothed part of Green's function
	subroutine Diver_Green_s(km, R_vl, Gs, DGs)
		implicit none
		real(dp), intent(in) :: R_vl
      complex(dp), intent(in) :: km 
      complex(dp), intent(out) :: DGs, Gs
      complex(dp) :: Exp_f
      
      Exp_f = exp(-im*km*R_vl) !
		if (R_vl == 0.0) then
			Gs = -im*km/(4*PI)
			DGs = km**2/(8*PI) 
		else 
			Gs = (Exp_f-1)/(4*PI*R_vl) + km**2*R_vl/(8*PI) !
			DGs = (1.0 -(1.0 + im*km*R_vl)*Exp_f)/(4*PI*R_vl**2) + km**2/(8*PI) 
		end if
		return 
    end subroutine Diver_Green_s 
   
    subroutine fn_parameter(struct, m, ot, p1, p2, pc, lng, area)
      implicit none
		type(structure_tri), intent(in) :: struct
		integer, intent(in) :: m, ot
      real(dp), dimension(3), intent(out) :: p1, p2, pc
      real(dp), intent(out) :: lng, area
      integer :: j, nv, p
		type(Pair) :: pair_p
		
      pair_p = struct%neighbours(m)
      do j = 1, 3                
         if (struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(1) &
               .and. struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(2))then
               pc(1:3)=struct%points(struct%elements(pair_p%element(ot))&
               %corners(j))%point ! vertex of the triangle
			nv = j
         end if
      end do
		p = struct%neighbours(m)%element(ot) !element index
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
			lng = vec_len(p1-p2)
			area = 0.5*vec_len(cross_r((p1-pc), (p2-pc))) 
			return
    end subroutine fn_parameter
	 
		
	!subroutine fn_parameter_type_sphere(struct, m, ot, center, fn)
 !     implicit none
	!	type(structure_tri), intent(in) :: struct
	!	integer, intent(in) :: m, ot
	!	real(dp), intent(in) :: center(3) !centroid of the sphere or any ojbect
	!	
 !     real(dp), dimension(3) :: p1, p2, pc !, intent(out)
 !     real(dp) :: lng, area
 !     integer :: j, nv, p
	!	type(Pair) :: pair_p
	!	type(fn_rwg), intent(out) :: fn
	!	
 !     pair_p = struct%neighbours(m)
 !     do j = 1, 3
 !        if (struct%elements(pair_p%element(ot))%corners(j) &
 !              /= pair_p%corner(1) &
 !              .and. struct%elements(pair_p%element(ot))%corners(j) &
 !              /= pair_p%corner(2))then
 !              pc(1:3)=struct%points(struct%elements(pair_p%element(ot))&
 !              %corners(j))%point ! vertex of the triangle
	!		nv = j
 !        end if
 !     end do
	!	
	!	p = struct%neighbours(m)%element(ot) !element index
	!	if (nv == 1) then		
	!		p1 = struct%points(struct%elements(p)%corners(2))%point
	!		p2 = struct%points(struct%elements(p)%corners(3))%point
	!	else if (nv == 2) then
	!		p1 = struct%points(struct%elements(p)%corners(3))%point
	!		p2 = struct%points(struct%elements(p)%corners(1))%point
	!	else 
	!		p1 = struct%points(struct%elements(p)%corners(1))%point
	!		p2 = struct%points(struct%elements(p)%corners(2))%point
	!	end if		
	!		lng = vec_len(p1-p2)
	!		area = 0.5*vec_len(cross_r((p1-pc), (p2-pc))) 
	!		fn%corners(1)%point = p1
	!		fn%corners(2)%point = p2
	!		fn%corners(3)%point = pc
	!		fn%area = area
	!		fn%length = lng
	!		fn%sign = -1*(((ot-1)*ot)-1)
	!		fn%df = fn%sign*lng/area!
	!		fn%mid = struct%midpoint(p)%point					
	!		fn%nq = norm_vec_r_sphere(fn%corners, fn%mid, center)
	!		return    
	! end subroutine fn_parameter_type_sphere
	 
	subroutine fn_parameter_type(struct, m, ot, center, fn)
      implicit none
		type(structure_tri), intent(in) :: struct
		integer, intent(in) :: m, ot
		real(dp), intent(in) :: center(3) !centroid of the sphere or any ojbect
		
      real(dp), dimension(3) :: p1, p2, pc !, intent(out)
      real(dp) :: lng, area
      integer :: j, nv, p
		type(Pair) :: pair_p
		type(fn_rwg), intent(out) :: fn
		
      pair_p = struct%neighbours(m)
      do j = 1, 3
         if (struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(1) &
               .and. struct%elements(pair_p%element(ot))%corners(j) &
               /= pair_p%corner(2))then
               pc(1:3)=struct%points(struct%elements(pair_p%element(ot))&
               %corners(j))%point ! vertex of the triangle
			nv = j
         end if
      end do
		
		p = struct%neighbours(m)%element(ot) !element index
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
			lng = vec_len(p1-p2)
			area = 0.5*vec_len(cross_r((p1-pc), (p2-pc))) 
			fn%corners(1)%point = p1
			fn%corners(2)%point = p2
			fn%corners(3)%point = pc
			fn%area = area
			fn%length = lng
			fn%sign = -1*(((ot-1)*ot)-1)
			fn%df = fn%sign*lng/area!
			fn%mid = struct%midpoint(p)%point					
			fn%nq = norm_vec_r_new(fn%corners, fn%mid, center)
			return    
	 end subroutine fn_parameter_type
		 
	 subroutine test_address_csr_element()
		integer, dimension(:), allocatable :: b_csr, ia, ja
		integer :: i, j, b
		
		allocate(b_csr(13), ja(1), ia(6))
		
		b_csr =  (/1, -1, -3, -2, 5, 4, 6, 4, -1, 2, 7, 8, -5/)
		ja  = (/1, 2, 4, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5/)
		ia = (/1, 4, 6, 9, 12, 14/)		
		i = 5 
		j = 5		
		call address_elements_in_ncsr_format(i, j, ia, ja, b_csr, b)
		print*, 'b=', b
	end subroutine 
	
	subroutine address_elements_in_ncsr_format(i, j, ia, ja, a_csr, x)
		integer, intent(in) :: i, j
		integer, dimension(:), intent(in) :: a_csr
		
		integer, intent(out) :: x
		integer, dimension(:), intent(in) :: ia, ja
		integer :: is, it, m		

		x = 0
		is = ia(i)
		it = ia(i+1)-1
		do m = is, it
			if (ja(m) .eq. j) then
				x = a_csr(m)
			else				
			 cycle	
			end if		
		end do
	end subroutine
	
	subroutine address_elements_in_zcsr_format(i, j, ia, ja, a_csr, x)
		integer, intent(in) :: i, j
		complex(kind=dp), dimension(:), intent(in) :: a_csr	
		complex(dp), intent(out) :: x	
		integer, dimension(:), intent(in) :: ia, ja
		integer :: is, it, m
		
		x = (0.0, 0.0)
		is = ia(i)
		it = ia(i+1)-1

		do m = is, it
			if (ja(m) .eq. j) then
				x = a_csr(m)
			else				
			 cycle	
			end if		
		end do
	end subroutine
	
	subroutine normal_integration_new(m, n, tt, ll, ngp, sum_a, sum_b, sum_c)
		implicit none		
		integer, intent(in) :: m, n, ngp, tt, ll		
		complex(dp), intent(out) :: sum_a, sum_b, sum_c
		
		real(dp) :: TRF, c_tmp
		real(dp) :: w(20), a(20), b(20)
		real(dp) :: Rq(3), Ra(3), r_fm(3), r_fn(3), dfm, dfn, Rn(3), R(3) 
		complex(dp) :: f1, f2, DG_1, DG_2, G_1, G_2, c1, c2
		
		type(fn_rwg) :: fm, fn
		integer :: i, j
		
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		call Quadrature_tri(ngp, a, b, w)
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		dfm = fm%sign*fm%length!      
		dfn = fn%sign*fn%length
		
		do i = 1, ngp
			Ra = r_simplex_new(a(i), b(i), fm%corners)
			r_fm = 0.5*(Ra-fm%corners(3)%point)*dfm
			do j = 1, ngp 
				Rq = r_simplex_new(a(j), b(j), fn%corners)
				r_fn = 0.5*(Rq-fn%corners(3)%point)*dfn
				R = Ra - Rq				
				Rn = R/vec_len(R)
				call Diver_Green(k1, vec_len(R), G_1, DG_1)
				call Diver_Green(k2, vec_len(R), G_2, DG_2)						
				TRF = w(i)*w(j)
				sum_a = sum_a + dot_product(r_fm, cross_r(r_fn, Rn)*(DG_1 + DG_2))*TRF 
				
				c_tmp = dot_product(r_fm, r_fn)
				c1 = (c_tmp - 1/(k1*k1)*dfm*dfn)*G_1
				c2 = (c_tmp - 1/(k2*k2)*dfm*dfn)*G_2
				f1 = im*Omega*my_1*c1
				f2 = im*Omega*my_2*c2
			
				sum_b = sum_b + (f1 + f2)*TRF !D_11  					
				f1 = im*Omega*eps_1*c1
				f2 = im*Omega*eps_2*c2
				sum_c = sum_c + (f1 + f2)*TRF !D_22 
			end do !loop j
		end do  !loop i
	end subroutine normal_integration_new	
		
	subroutine init_result_inner_JMCFIE(rv)
		implicit none
		type(result_inner_JMCFIE), intent(inout) :: rv
		
		rv%IK = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		rv%IT_a = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		rv%IT_bn = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		rv%IT_bt = (0.0, 0.0)
	end subroutine
	
	subroutine init_element_JMCFIE(rv)
		implicit none
		type(result_element_JMCFIE), intent(inout) :: rv		
		rv%Kn = (0.0, 0.0)
		rv%Kt = (0.0, 0.0)
		rv%Tt = (0.0, 0.0)
		rv%Tn = (0.0, 0.0)	
	end subroutine
	
	subroutine tri_impedance_building_PMCHWT(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		integer, intent(in) :: ngp
		character(len = 20), intent(in) :: singular_type
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		
		!dummy		
		complex(dp) :: tmp_kt2, tmp_kt1
		integer :: md, pp, qq
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm
		logical :: Iqs_L, edge_ad
		allocate(zz(3))
		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		if (singular_type .eq. 'singular') then	
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)
			if ((pp .eq. qq)) then
	
				Iqs_L = .true.
				md = 1	
				rv_1 = tri_outer_integration_singular(fm, fn, ngp, md, Iqs_L)
				md = 2
				rv_2 = tri_outer_integration_singular(fm, fn, ngp, md, Iqs_L)			
			else 
			   Iqs_L = .false.
				edge_ad = .false.		
				call edge_adjacent_tri(m, n, tt, ll, edge_ad)
				md = 1
				rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
				md = 2
				rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)		
				
		      if (edge_ad) then
					md = 1
					tmp_kt1 = tri_outer_integration_edge_adjacent(  fm, fn, md)
					md = 2
					tmp_kt2 = tri_outer_integration_edge_adjacent(  fm, fn, md)
					rv_1%Kt = tmp_kt1
					rv_2%Kt = tmp_kt2
				end if				
			end if
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if
		zz(1) = eta_1*rv_1%Tt + eta_2*rv_2%Tt !z_11	
		zz(2) = -(rv_1%Kt + rv_2%Kt)  !z_12, just as the equation in Mitschang 94,  
	   zz(3) = 1/eta_1*rv_1%Tt + 1/eta_2*rv_2%Tt !z_22
		return
	end subroutine tri_impedance_building_PMCHWT
	
	subroutine tri_impedance_building_JMCFIE(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		integer, intent(in) :: ngp
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		character(len=20), intent(in) :: singular_type
		
		! dummy
		integer :: md, pp, qq
		complex(dp) :: tmp_kt2, tmp_kt1
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm		
		logical :: Iqs_L,  edge_ad
		
		allocate(zz(3))
		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
				
		if (singular_type .eq. 'singular') then			
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)					
			if (pp .eq. qq) then
				Iqs_L = .true.
				md = 1		
				rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
				md = 2
				rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)			
			else 
			   Iqs_L = .false.
				edge_ad = .false.		
				call edge_adjacent_tri(m, n, tt, ll, edge_ad)
				md = 1
				rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
				md = 2
				rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
				
			   if (edge_ad) then
					md = 1
					tmp_kt1 = tri_outer_integration_edge_adjacent(  fm, fn, md)
					md = 2
					tmp_kt2 = tri_outer_integration_edge_adjacent(  fm, fn, md)					
					rv_1%Kt = tmp_kt1
					rv_2%Kt = tmp_kt2
				end if				
			end if
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if		
		zz(1) = rv_1%Kn - rv_2%Kn + rv_1%Tt + rv_2%Tt !z_11
		zz(2) = -(1/eta_1*(rv_1%Kt - rv_1%Tn) + 1/eta_2*(rv_2%Kt + rv_2%Tn)) !z_12
		zz(3) = -eta_1*(rv_1%Tn - rv_1%Kt) + eta_2*(rv_2%Kt + rv_2%Tn) !z_21
		return
	end subroutine tri_impedance_building_JMCFIE
	
		
	!The same as ICTF, just for test so that a compact expression can be used.
	subroutine tri_impedance_building_MCTF2(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		integer, intent(in) :: ngp
		character(len=20), intent(in) :: singular_type
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		!dummy
		integer :: md, pp, qq		
		complex(dp) :: a_lc(2), d_lc(2)
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm
		logical :: Iqs_L
		
		allocate(zz(4))
		
		call formulation_coefficient(a_lc, d_lc)
		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		if (singular_type .eq. 'singular') then
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)					
			if ((pp == qq)) then
				Iqs_L = .true.
			else
				Iqs_L = .false.
			end if		
			md = 1		
			rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
			md = 2
			rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)			
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if
		zz(1) = rv_1%Tt*a_lc(1) + rv_2%Tt*a_lc(2) !z_11
		zz(2) = -(rv_1%Kt*a_lc(1)/eta_1 + rv_2%Kt*a_lc(2)/eta_2)!z_12
		zz(3) = eta_1*d_lc(1)*rv_1%Kt + rv_2%Kt*eta_2*d_lc(2) !z_21
		zz(4)	= d_lc(1)*rv_1%Tt + d_lc(2)*rv_2%Tt
		
	end subroutine tri_impedance_building_MCTF2	
	
	subroutine tri_impedance_building_ICTF(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		integer, intent(in) :: ngp
		character(len=20), intent(in) :: singular_type
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		!dummy
		integer :: md, pp, qq
		!complex(dp) ::  eta_a !
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm
		logical :: Iqs_L
		
		allocate(zz(4))
		
		!eta_a = (eta_1 + eta_2)/2
		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		if (singular_type .eq. 'singular') then
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)					
			if ((pp == qq)) then
				Iqs_L = .true.
			else
				Iqs_L = .false.
			end if		
			md = 1		
			rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
			md = 2
			rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)			
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if
		zz(1) = (eta_1*rv_1%Tt + eta_2*rv_2%Tt)/eta_a !z_11
		zz(2) = -(rv_1%Kt + rv_2%Kt)/eta_a !z_12
		zz(3) = (rv_1%Kt + rv_2%Kt)*eta_a !z_21
		zz(4) = (rv_1%Tt/eta_1 + rv_2%Tt/eta_2)*eta_a !z_22
		
	end subroutine tri_impedance_building_ICTF
	
	subroutine tri_impedance_building_CNF(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		integer, intent(in) :: ngp
		character(len=20), intent(in) :: singular_type
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		
		!dummy	
		integer :: md, pp, qq
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm
		logical :: Iqs_L
		
		allocate(zz(3))
	
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		if (singular_type .eq. 'singular') then
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)
			if ((pp == qq)) then
				Iqs_L = .true.
			else
				Iqs_L = .false.
			end if		
			md = 1		
			rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
			md = 2
			rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if
		
		zz(1) = rv_1%Kn - rv_2%Kn   !z_11
		zz(2) = rv_1%Tn/eta_1 - rv_2%Tn/eta_2  !z_12
		zz(3) = eta_2*rv_2%Tn - eta_1*rv_1%Tn  !z_21
	end subroutine tri_impedance_building_CNF
		
	subroutine tri_impedance_building_MNMF(m, n, tt, ll, ngp, zz, singular_type)
		implicit none
		integer, intent(in) :: m, n, tt, ll
		complex(dp), dimension(:), allocatable, intent(out) :: zz
		integer, intent(in) :: ngp
		character(len=20), intent(in) :: singular_type

		!dummy
		integer :: md, pp, qq		
		complex(dp) :: b1, b2, c1, c2
		type(result_element_JMCFIE) :: rv_1, rv_2
		type(fn_rwg) :: fn, fm
		logical :: Iqs_L
		
		allocate(zz(4))
		b1 = eps_1/(eps_1 + eps_2)
		b2 = eps_2/(eps_1 + eps_2)
		c1 = my_1/(my_1 + my_2)
		c2 = my_2/(my_1 + my_2)
		
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		if (singular_type .eq. 'singular') then	
			pp = struc_tri%neighbours(m)%element(tt)
			qq = struc_tri%neighbours(n)%element(ll)
			
			if ((pp == qq)) then
				Iqs_L = .true.
			else
				Iqs_L = .false.
			end if		
			md = 1		
			rv_1 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
			md = 2
			rv_2 = tri_outer_integration_singular(  fm, fn, ngp, md, Iqs_L)
		else
			md = 1		
			rv_1 = tri_outer_integration(  fm, fn, ngp, md)
			md = 2
			rv_2 = tri_outer_integration(  fm, fn, ngp, md)
		end if
		zz(1) = b1*rv_1%Kn - b2*rv_2%Kn !z_11
		zz(2) = b1/eta_1*rv_1%Tn - b2/eta_2*rv_2%Tn  !z_12
		zz(3) = eta_2*c2*rv_2%Tn - eta_1*c1*rv_1%Tn   !z_21
		zz(4) = c1*rv_1%Kn - c2*rv_2%Kn
	end subroutine tri_impedance_building_MNMF
	
	function tri_inner_integration_edge_adjacent(r_a, fn, ngp, kd)result(Ik)
		implicit none
		integer, intent(in) :: ngp
		real(dp), intent(in) :: r_a(3)
		complex(dp), intent(in) :: kd		
		type(fn_rwg) :: fn
		complex(dp) :: Ik(3)
		
		!dummy
		integer :: j 
		real(dp) :: R(3), Rn(3), r_q(3), r_fn(3), TRF
		real(dp), dimension(20) :: a, b, w		
		complex(dp) :: G_d, DG_d
		
		Ik = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		call Quadrature_tri(ngp, a, b, w)
		do j = 1, ngp 
			r_q = r_simplex_new(a(j), b(j), fn%corners)
			r_fn = 0.5*(r_q-fn%corners(3)%point)*fn%df
			R = r_a - r_q
			Rn = R/vec_len(R)
			call Diver_Green(kd, vec_len(R), G_d, DG_d)
			TRF = w(j)*fn%area			
			Ik = Ik + cross_r(r_fn, Rn)*DG_d*TRF !
		end do !loop j
	end function tri_inner_integration_edge_adjacent
	
	!partly according to Ergül 2009
	!IEEE Trans. Antenna and Propagation
	function tri_inner_integration(r_a, fn, ngp, kd)result(rv)
		implicit none
		integer, intent(in) :: ngp
		real(dp), intent(in) :: r_a(3)
		complex(dp), intent(in) :: kd		
		type(fn_rwg) :: fn
		type(result_inner_JMCFIE) :: rv
		
		!dummy
		integer :: j 
		real(dp) :: R(3), Rn(3), r_q(3), r_fn(3), TRF
		real(dp), dimension(20) :: a, b, w		
		complex(dp) :: G_d, DG_d
		
		call Quadrature_tri(ngp, a, b, w)  
		call init_result_inner_JMCFIE(rv)
		
		do j = 1, ngp 
			r_q = r_simplex_new(a(j), b(j), fn%corners)
			r_fn = 0.5*(r_q-fn%corners(3)%point)*fn%df
			R = r_a - r_q
			Rn = R/vec_len(R)
			call Diver_Green(kd, vec_len(R), G_d, DG_d)
			TRF = w(j)*fn%area			
			rv%IK = rv%IK + cross_r(r_fn, Rn)*DG_d*TRF !
			rv%IT_a = rv%IT_a + r_fn*G_d*TRF
			rv%IT_bn = rv%IT_bn + fn%df*Rn*DG_d*TRF
			rv%IT_bt = rv%IT_bt + fn%df*G_d*TRF
		end do !loop j
	end function tri_inner_integration
	
	function tri_inner_singular_integration_smooth(r_a, fn, ngp, kd, Iqs_L)result(rv)
		implicit none
		integer, intent(in) :: ngp
		real(dp), intent(in) :: r_a(3)
		complex(dp), intent(in) :: kd		
		logical, intent(in) :: Iqs_L
		type(fn_rwg), intent(in) :: fn
		type(result_inner_JMCFIE) :: rv
		
		!dummy
		integer :: j 
		real(dp) :: R(3), Rn(3), r_q(3), r_fn(3), TRF
		real(dp), dimension(20) :: a, b, w
		
		complex(dp) :: Gs_d, DGs_d
		
		call Quadrature_tri(ngp, a, b, w)  
		call init_result_inner_JMCFIE(rv)
		
		do j = 1, ngp 
			r_q = r_simplex_new(a(j), b(j), fn%corners)
			r_fn = 0.5*(r_q-fn%corners(3)%point)*fn%df
			R = r_a - r_q
			Rn = R/vec_len(R)
			call Diver_Green_s(kd, vec_len(R), Gs_d, DGs_d)
			TRF = w(j)*fn%area	
			if (Iqs_L .eqv. .true.) then
				rv%IK = (0.0, 0.0)
				rv%IT_bn = (0.0, 0.0)
			else !
				rv%IK = rv%IK + cross_r(r_fn, Rn)*DGs_d*TRF !
				rv%IT_bn = rv%IT_bn + fn%df*Rn*DGs_d*TRF
			end if
			rv%IT_a = rv%IT_a + r_fn*Gs_d*TRF
			rv%IT_bt = rv%IT_bt + fn%df*Gs_d*TRF
		end do !loop j
	end function tri_inner_singular_integration_smooth
	
	function tri_inner_singular_integration_singular(r_a, fn, kd, Iqs_L)result(rv)
		implicit none
		
		real(dp), intent(in) :: r_a(3)
		complex(dp), intent(in) :: kd		
		logical, intent(in) :: Iqs_L
		type(fn_rwg), intent(in) :: fn
		type(result_inner_JMCFIE) :: rv
		
		!dummy
		integer :: q!, ngp
		real(dp) :: k2_n1(3), k2_p1(3), k4_n1(3), k4_p1(3)
		real(dp) :: Iqs_n1, Iqs_p1, k3_n1(3), k3_p1(3)
		
		call init_result_inner_JMCFIE(rv)	
		q = -1
		Iqs_n1 = Iq_s(q, r_a, fn, Iqs_L)
		K2_n1 = K2q(q, r_a, fn, Iqs_L)!		
		K3_n1 = K3q(q, r_a, fn, Iqs_L)!			
		K4_n1 = K4q(q, r_a, fn, Iqs_L)!
			
		q = 1
		Iqs_p1 = Iq_s(q, r_a, fn, Iqs_L)
		K2_p1 = K2q(q, r_a, fn, Iqs_L)!
		K3_p1 = K3q(q, r_a, fn, Iqs_L)!
		K4_p1 = K4q(q, r_a, fn, Iqs_L)!
				
		if (Iqs_L) then
			rv%IK = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		else !
			rv%IK = 1/(4*PI)*K4_n1 - kd**2/(8*PI)*K4_p1
		end if
		rv%IT_a =  1/(4*PI)*K2_n1 - kd**2/(8*PI)*K2_p1
		rv%IT_bn = -(fn%df/(4*PI)*K3_n1 - kd**2*fn%df/(8*PI))*K3_p1
		rv%IT_bt = fn%df/(4*PI)*Iqs_n1-kd**2*fn%df/(8*PI)*Iqs_p1		
		return
	end function tri_inner_singular_integration_singular
	
	function quadrature_integral(r_a, fn, ngp)result(k3)
		type(fn_rwg) :: fn
		real(dp) :: r_a(3), k2(3), r_q(3), r_fn(3), k3(3)
		integer :: ngp, j
		
		!dummy
		real(dp), dimension(20) :: a, b, w
		k2 = (/0.0, 0.0, 0.0/)
		k3 = (/0.0, 0.0, 0.0/)
		call Quadrature_tri(ngp, a, b, w) 
		do j = 1, ngp
			r_q = r_simplex_new(a(j), b(j), fn%corners)
			r_fn = 0.5*(r_q-fn%corners(3)%point)*fn%df
			k2 = k2 + (vec_len(r_a - r_q)*r_fn)*w(j)*fn%area
		end do
		return
	end function
	
	function tri_outer_integration_edge_adjacent(fm, fn, md)result(Kt)
		implicit none
		integer, intent(in) :: md !md is the medium index
		type(fn_rwg), intent(in) :: fn, fm
				
		!dummy
		integer :: j, ngp 
		real(dp) :: r_a(3), r_fn(3), TRF, r_fm(3)
		real(dp), dimension(20) :: a, b, w
		real(dp) :: x_nt
		complex(dp) :: kd, Kt
		complex(dp) :: Ik(3)
		
		ngp = 12
		call Quadrature_tri(ngp, a, b, w) 
		if (md .eq. 1)then
			kd = k1
		else 
			kd = k2
		end if
		Kt = (0.0, 0.0)
		x_nt = 0.0
		do j = 1, ngp 
			r_a = r_simplex_new(a(j), b(j), fm%corners)
			r_fm = 0.5*(r_a-fm%corners(3)%point)*fm%df
			r_fn = 0.5*(r_a-fn%corners(3)%point)*fn%df
			TRF = w(j)*fm%area
			Ik = tri_inner_integration_edge_adjacent(r_a, fn, ngp, kd)			
			Kt = Kt + dot_product(r_fm, Ik)*TRF
			x_nt = x_nt - dot_product(r_fn, cross_r(fm%nq, r_fm))*TRF
		end do !loop j
		Kt = Kt + 0.5*(-1)**md*cmplx(x_nt, 0.0)		
	end function tri_outer_integration_edge_adjacent
	
	function tri_outer_integration_singular(fm, fn, ngp, md, Iqs_L)result(rv_out)
		implicit none
		integer, intent(in) :: ngp, md !md is the index of the medium
		logical, intent(in) :: Iqs_L
		type(fn_rwg), intent(in) :: fn, fm
		type(result_element_JMCFIE) :: rv_out
		
		!dummy
		integer :: j 
		real(dp) :: r_a(3), r_fn(3), TRF, r_fm(3)
		real(dp), dimension(20) :: a, b, w
		real(dp) :: x_nn, x_nt
		complex(dp) :: kd
		complex(dp) :: tmp_s, tmp_v(3)
		type(result_inner_JMCFIE) :: rv_in, rv_sm, rv
		
		call Quadrature_tri(ngp, a, b, w)  
		call init_element_JMCFIE(rv_out)
		
		if (md .eq. 1)then
			kd = k1
		else 
			kd = k2
		end if
		
		x_nn = 0.0
		x_nt = 0.0
		do j = 1, ngp 
			r_a = r_simplex_new(a(j), b(j), fm%corners)
			r_fm = 0.5*(r_a-fm%corners(3)%point)*fm%df
			r_fn = 0.5*(r_a-fn%corners(3)%point)*fn%df
			TRF = w(j)*fm%area
			
			rv_in = tri_inner_singular_integration_singular(r_a, fn, kd, Iqs_L)
			rv_sm = tri_inner_singular_integration_smooth(r_a, fn, ngp, kd, Iqs_L)
			
			rv%It_bn = rv_in%It_bn + rv_sm%It_bn
			rv%It_bt = rv_in%It_bt + rv_sm%It_bt
			rv%IK = rv_in%IK + rv_sm%IK
			rv%It_a = rv_in%It_a + rv_sm%It_a
			
			rv_out%Kn = rv_out%Kn + dot_product(r_fm, cross_rc(fm%nq, rv%IK))*TRF !
			rv_out%Kt = rv_out%Kt + dot_product(r_fm, rv%IK)*TRF
			
			tmp_v = cross_rc(fm%nq, kd*rv%It_a) + cross_rc(fm%nq, 1.0/kd*rv%It_bn)
			rv_out%Tn = rv_out%Tn + im*TRF*dot_product(r_fm, tmp_v)
						
			tmp_s = 1.0/kd*fm%df*rv%IT_bt
			rv_out%Tt = rv_out%Tt + im*TRF*(kd*dot_product(r_fm, rv%It_a) - tmp_s)
			x_nt = x_nt - dot_product(cross_r(fm%nq, r_fm), r_fn)*TRF
			x_nn = x_nn - dot_product(r_fm, r_fn)*TRF			
		end do !loop j
		rv_out%Kn = rv_out%Kn + 0.5*(-1)**md*cmplx(x_nn, 0.0)
		rv_out%Kt = rv_out%Kt + 0.5*(-1)**md*cmplx(x_nt, 0.0)		
	end function tri_outer_integration_singular
	
	function tri_outer_integration(fm, fn, ngp, md)result(rv)
		implicit none
		
		integer, intent(in) :: ngp, md !md is the index of the medium		
		type(fn_rwg), intent(in) :: fn, fm
		type(result_element_JMCFIE) :: rv		
		
		!dummy
		integer :: j 
		real(dp) :: r_a(3), r_fn(3), TRF, r_fm(3)
		real(dp), dimension(20) :: a, b, w
		real(dp) :: x_nn, x_nt
		complex(dp) :: kd
		complex(dp) :: tmp_s, tmp_v(3)
		type(result_inner_JMCFIE) :: rv_in
		
		call Quadrature_tri(ngp, a, b, w)  
		call init_element_JMCFIE(rv)
		
		if (md .eq. 1)then
			kd = k1
		else 
			kd = k2
		end if
		
		x_nn = 0.0
		x_nt = 0.0
		do j = 1, ngp 
			r_a = r_simplex_new(a(j), b(j), fm%corners)
			r_fm = 0.5*(r_a-fm%corners(3)%point)*fm%df
			r_fn = 0.5*(r_a-fn%corners(3)%point)*fn%df
			TRF = w(j)*fm%area
			rv_in = tri_inner_integration(r_a, fn, ngp, kd)
			rv%Kn = rv%Kn + dot_product(r_fm, cross_rc(fm%nq, rv_in%IK))*TRF !
			rv%Kt = rv%Kt + dot_product(r_fm, rv_in%IK)*TRF
			
			tmp_v = cross_rc(fm%nq, (kd*rv_in%It_a + 1.0/kd*rv_in%It_bn))
			rv%Tn = rv%Tn + im*TRF*dot_product(r_fm, tmp_v)
			
			tmp_s = 1.0/kd*fm%df*rv_in%IT_bt
			rv%Tt = rv%Tt + im*TRF*(kd*dot_product(r_fm, rv_in%It_a) - tmp_s)
			x_nt = x_nt - dot_product(cross_r(fm%nq, r_fm), r_fn)*TRF
			x_nn = x_nn - dot_product(r_fm, r_fn)*TRF			
		end do !loop j
		rv%Kn = rv%Kn + 0.5*(-1)**md*cmplx(x_nn, 0.0)
		rv%Kt = rv%Kt + 0.5*(-1)**md*cmplx(x_nt, 0.0)
			
	end function tri_outer_integration	
	
	subroutine Singular_Integration(m, n, tt, ll, ngp, suma, sumb, sumc)
		implicit none
		integer, intent(in) :: m, n, ngp, tt, ll		
		complex(dp), intent(out) :: suma, sumb, sumc
		
		integer :: i, j, q, pp, qq
		real(dp) :: TRF, TRF_m
		real(dp), dimension(100) :: w(100), a(100), b(100)
		real(dp) :: Rq(3), Rn(3), Ra(3), R(3), r_fm(3), r_fn(3)!
		real(dp) :: k2_n1(3), k2_p1(3), k4_n1(3), k4_p1(3)
		real(dp) :: pc, Iqs_n1, Iqs_p1
		type(fn_rwg) :: fm, fn
		
		complex(dp) :: suma_s, sumb_s, sumc_s
		complex(dp) :: Gs_1, Gs_2, DGs_1, DGs_2 
		complex(dp) :: cc(4), f1, f2, c3
		logical :: Iqs_L
		
      call Quadrature_tri(ngp, a, b, w)
      suma = (0.0, 0.0)
      sumb = (0.0, 0.0)
      sumc = (0.0, 0.0)
      suma_s = (0.0, 0.0)
      sumb_s = (0.0, 0.0)
      sumc_s = (0.0, 0.0)
		
		k2_p1  = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		
		pp = struc_tri%neighbours(m)%element(tt)
		qq = struc_tri%neighbours(n)%element(ll)
					
		if ((pp == qq)) then
			Iqs_L = .true.
		else
			Iqs_L = .false.
		end if
		
		do i = 1, ngp
			Ra = r_simplex_new(a(i), b(i), fm%corners)
			r_fm = (Ra - fm%corners(3)%point)*fm%df*0.5
			do j = 1, ngp 
				Rq = r_simplex_new(a(j), b(j), fn%corners)
				r_fn = (Rq - fn%corners(3)%point)*fn%df*0.5
				R = Ra - Rq					
				Rn = R/vec_len(R)
				call Diver_Green_s(k1, vec_len(R), Gs_1, DGs_1)
				call Diver_Green_s(k2, vec_len(R), Gs_2, DGs_2)				
				TRF = w(i)*w(j)*fn%area*fm%area
				if (Iqs_L .eqv. .true.) then
					suma = (0.0, 0.0)
				else !K_12s singular integration of K_12	
					suma = suma + dot_product(r_fm, cross_r(r_fn, Rn)*(DGs_1 + DGs_2))*TRF								
				end if				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				f1 = im*(Omega*my_1*dot_product(r_fm, r_fn) - 1/(Omega*eps_1)*fm%df*fn%df)*Gs_1 
				f2 = im*(Omega*my_2*dot_product(r_fm, r_fn) - 1/(Omega*eps_2)*fm%df*fn%df)*Gs_2 
				sumb = sumb + (f1 + f2)*TRF !D_11 
                        
				f1 = im*(Omega*eps_1*dot_product(r_fm, r_fn) - 1/(Omega*my_1)*fm%df*fn%df)*Gs_1 
				f2 = im*(Omega*eps_2*dot_product(r_fm, r_fn) - 1/(Omega*my_2)*fm%df*fn%df)*Gs_2 
				sumc = sumc + (f1 + f2)*TRF !
			end do ! loop j	
         
			TRF_m = w(i)*fm%area
			pc = fm%df*fn%df/Omega
			!print*, 'pc= ', pc
			q = -1
			k2_n1 = K2q(q, Ra, fn, Iqs_L)
			k4_n1 = K4q(q, Ra, fn, Iqs_L)			
			Iqs_n1 = Iq_s(q, Ra, fn, Iqs_L) !
					
			q = 1 
			k2_p1 = K2q(q, Ra, fn, Iqs_L)
			k4_p1 = K4q(q, Ra, fn, Iqs_L) !due to the divergence to r
			Iqs_p1 = Iq_s(q, Ra, fn, Iqs_L) !
			
			c3 = -(k1**2 + k2**2)/(8*PI)		
			if (Iqs_L .eqv. .true.) then
				suma_s = (0.0, 0.0)
			else
				suma_s = suma_s + (1/(2*PI)*dot_product(r_fm, k4_n1) + c3*dot_product(r_fm, k4_p1))*TRF_m !K_21s
			end if			
			call Coefficient_D11s(  pc, cc)	!D_11s			
			sumb_s = sumb_s + (dot_product(r_fm, (cc(1)*k2_n1 + cc(3)*k2_p1)) + (cc(2)*Iqs_n1 + cc(4)*Iqs_p1))*TRF_m !
			
			call Coefficient_D22s(  pc, cc)	!D_22s					
			sumc_s = sumc_s + (dot_product(r_fm, (cc(1)*k2_n1 + cc(3)*k2_p1)) + (cc(2)*Iqs_n1 + cc(4)*Iqs_p1))*TRF_m ! 			
		end do !loop i

		suma = (suma + suma_s)
		sumb = (sumb + sumb_s)
		sumc = (sumc + sumc_s)		
	end subroutine Singular_Integration  
	
	! singular integration = smooth_part + singular_part
	subroutine singular_integration_part_smooth(fm, fn, ngp, suma, sumb, sumc, Iqs_L)	
		implicit none
		integer, intent(in) :: ngp
		complex(dp), intent(out) :: suma, sumb, sumc		
		integer :: i, j
		real(dp) :: TRF
		real(dp), dimension(100) :: w(100), a(100), b(100)
		real(dp) :: Rq(3), Rn(3), Ra(3), R(3)
		real(dp) :: r_fm(3), r_fn(3) 		
				
		complex(dp) :: Gs_1, Gs_2, DGs_1, DGs_2 
		complex(dp) :: f1, f2 
		logical, intent(in) :: Iqs_L
		type(fn_rwg), intent(in) :: fm, fn
				
      call Quadrature_tri(ngp, a, b, w)
      suma = (0.0, 0.0)
      sumb = (0.0, 0.0)
      sumc = (0.0, 0.0)
      
		do i = 1, ngp
			Ra = r_simplex_new(a(i), b(i), fm%corners)
			r_fm = (Ra-fm%corners(3)%point)*fm%df*0.5 
			do j = 1, ngp 
				Rq = r_simplex_new(a(j), b(j), fn%corners)
				r_fn = (Rq - fn%corners(3)%point)*fn%df*0.5
				R = Ra - Rq			
				Rn = R/vec_len(R)
				call Diver_Green_s(k1, vec_len(R), Gs_1, DGs_1)
				call Diver_Green_s(k2, vec_len(R), Gs_2, DGs_2)				
				TRF = w(j)*fn%area*w(i)*fm%area
				if (Iqs_L .eqv. .true.) then
					suma = (0.0, 0.0)
				else !K_12s singular integration of K_12	
					suma = suma + dot_product(r_fm, cross_r(r_fn, Rn)*(DGs_1 + DGs_2))*TRF								
				end if		
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				f1 = im*(Omega*my_1*dot_product(r_fm, r_fn) - 1/(Omega*eps_1)*fm%df*fn%df)*Gs_1 
				f2 = im*(Omega*my_2*dot_product(r_fm, r_fn) - 1/(Omega*eps_2)*fm%df*fn%df)*Gs_2 
				sumb = sumb + (f1 + f2)*TRF !D_11 
                        
				f1 = im*(Omega*eps_1*dot_product(r_fm, r_fn) - 1/(Omega*my_1)*fm%df*fn%df)*Gs_1 
				f2 = im*(Omega*eps_2*dot_product(r_fm, r_fn) - 1/(Omega*my_2)*fm%df*fn%df)*Gs_2 
				sumc = sumc + (f1 + f2)*TRF !
			end do ! loop j	
		end do
	end subroutine singular_integration_part_smooth
	
	! Singular integration = smooth_part + singular_part
	! This is singular part without edge adjacent
	subroutine singular_integration_part_singular(fm, fn, ngp, Iqs_L, suma_s, sumb_s, sumc_s)
		implicit none
		integer, intent(in) ::ngp
		type(fn_rwg), intent(in) :: fm, fn
		logical, intent(in) :: Iqs_L		
		complex(dp), intent(out) :: suma_s, sumb_s, sumc_s
		
		!dummy
		integer :: i, q
		real(dp) :: TRF_m
		real(dp), dimension(100) :: w(100), a(100), b(100)
		real(dp) :: Ra(3), r_fm(3)!, r_fp(3)
		real(dp) :: k2_n1(3), k2_p1(3), k4_n1(3), k4_p1(3) 
		real(dp) :: pc, Iqs_n1, Iqs_p1
		
		complex(dp) :: cc(4), c3
		
      call Quadrature_tri(ngp, a, b, w)     
      suma_s = (0.0, 0.0)
      sumb_s = (0.0, 0.0)
      sumc_s = (0.0, 0.0)
		!
		do i = 1, ngp
			Ra = r_simplex_new(a(i), b(i), fm%corners)
			r_fm = (Ra-fm%corners(3)%point)*fm%df*0.5    
			TRF_m = w(i)*fm%area
			pc = fm%df*fn%df/Omega
			q = -1
			k2_n1 = K2q(q, Ra, fn, Iqs_L)
			k4_n1 = K4q(q, Ra, fn, Iqs_L)
			Iqs_n1 = Iq_s(q, Ra, fn, Iqs_L) !
					
			q = 1 
			k2_p1 = K2q(q, Ra, fn, Iqs_L)
			k4_p1 = K4q(q, Ra, fn, Iqs_L)!			
			Iqs_p1 = Iq_s(q, Ra, fn, Iqs_L) !			
			c3 = -(k1**2 + k2**2)/(8*PI)
			if (Iqs_L .eqv. .true.) then
				suma_s = (0.0, 0.0)
			else
				! KA + KB
				suma_s = suma_s + (1/(2*PI)*dot_product(r_fm, k4_n1) + c3*dot_product(r_fm, k4_p1))*TRF_m !K_21s
			end if
			
			call Coefficient_D11s(  pc, cc)	!D_11s			
			sumb_s = sumb_s + (dot_product(r_fm, (cc(1)*k2_n1 + cc(3)*k2_p1)) + (cc(2)*Iqs_n1 + cc(4)*Iqs_p1))*TRF_m !
			
			call Coefficient_D22s(  pc, cc)	!D_22s
			sumc_s = sumc_s + (dot_product(r_fm, (cc(1)*k2_n1 + cc(3)*k2_p1)) + (cc(2)*Iqs_n1 + cc(4)*Iqs_p1))*TRF_m !
		end do !loop i  
	end subroutine singular_integration_part_singular
	
	subroutine singular_integration_new(m, n, tt, ll, ngp, suma, sumb, sumc)
		integer, intent(in) :: m, n, ngp, tt, ll
		complex(dp) :: suma_s, sumb_s, sumc_s, suma_tmp	
		complex(dp) :: suma_sin, sumb_sin, sumc_sin	
		complex(dp), intent(out) :: suma, sumb, sumc
		integer :: qq, pp
		logical :: Iqs_L		
		logical :: edge_ad
		type(fn_rwg) :: fm, fn
				
		fm = structure_edges(m)%fn_ot(tt)%fn
		fn = structure_edges(n)%fn_ot(ll)%fn
		pp = struc_tri%neighbours(m)%element(tt)
		qq = struc_tri%neighbours(n)%element(ll)
		if ((pp == qq)) then
			Iqs_L = .true.
		else
			Iqs_L = .false.
		end if
		
		edge_ad = .false.
		if (.not. Iqs_L) then 
			call edge_adjacent_tri(m, n, tt, ll, edge_ad)
		end if
				
		call singular_integration_part_smooth(  fm, fn, ngp, suma_s, sumb_s, sumc_s, Iqs_L)		
		if (edge_ad) then			
			call singular_integration_part_singular(  fm, fn, ngp, Iqs_L, suma_sin, sumb_sin, sumc_sin)
			call singular_integration_edge_adjacent_ngp(  fm, fn, suma_tmp)
			suma = suma_tmp
			
			sumb = sumb_s + sumb_sin
			sumc = sumc_s + sumc_sin
		else 
			call singular_integration_part_singular(  fm, fn, ngp, Iqs_L, suma_sin, sumb_sin, sumc_sin)
			suma = suma_s + suma_sin
			sumb = sumb_s + sumb_sin
			sumc = sumc_s + sumc_sin
		end if				
	end subroutine singular_integration_new	
	
	subroutine edge_adjacent_tri(m, n, tt, ll, edge_ad)
		integer, intent(in) :: m, n, tt, ll
		logical, intent(out) :: edge_ad	
		integer :: counter, i, j
		
		type(Pair):: pair_m, pair_n
		pair_m = struc_tri%neighbours(m)
		pair_n = struc_tri%neighbours(n)
		counter = 0
      do j = 1, 3
			do i = 1, 3
         if (struc_tri%elements(pair_m%element(tt))%corners(j) &
               .eq. struc_tri%elements(pair_n%element(ll))%corners(i))then					
				counter = counter + 1         
			else 			
				continue
			end if
			end do
      end do
		if (counter .eq. 2) then
			edge_ad = .true.
		else if (counter == 3) then
			edge_ad = .false.
		end if	
		!print*, 'counter =', counter
	end subroutine edge_adjacent_tri
	
	subroutine lib_sie_tri_input_field_calculation(field_input)!result()	
		implicit none
		
		integer :: n1, n2, m, n, ngp!
		
		real(dp) :: r_local_tmp(3), I_00 
		type(vector_c) :: E_00
		
		type(lib_sie_evaluation_point_type), dimension(:), allocatable, intent(out) :: field_input		
	
		n1 = evaluation_parameter%N_dim(1)
		n2 = evaluation_parameter%N_dim(2)
		
		if (allocated(field_input)) then
			deallocate(field_input)
		end if
		
		allocate(field_input(n1*n2))
	
		if (pre_types%illumination == 'Gaussian') then
			do m = 1, n1*n2	
				call Gaussian_beam(beam_waist, r_media(m)%point, illumination_p, field_input(m)%e_field, k1)				
				field_input(m)%h_field%vector = cross_c(k1*illumination_p%k_in, field_input(m)%e_field%vector)/(my_1*omega)
				ngp = ngp_near_field_distance_tri(r_media(m)%point, nearfield_distance)
				r_media(m)%ngp = ngp
			end do
		else if (pre_types%illumination == 'Conical') then			
			do m = 1, n1*n2
				call conical_illumination(p_obj, r_media(m)%point, illumination_p, field_input(m)%e_field)	
				ngp = ngp_near_field_distance_tri(r_media(m)%point, nearfield_distance)
				r_media(m)%ngp = ngp
				field_input(m)%h_field%vector = cross_c(k1*illumination_p%k_in, Field_input(m)%e_field%vector)/(my_1*omega)
			end do			
		else if ((pre_types%illumination == 'Plane') .or. (pre_types%illumination == 'Conical_plane')) then
			do m = 1, n1*n2
				field_input(m)%e_field%vector = illumination_p%E_in*exp(-im*dot_product(k1*illumination_p%k_in, r_media(m)%point))
				ngp = ngp_near_field_distance_tri(r_media(m)%point, nearfield_distance)
				Field_input(m)%h_field%vector = cross_c(k1*illumination_p%k_in, Field_input(m)%e_field%vector)/(my_1*omega)
				r_media(m)%ngp = ngp
			end do	
		else 
			print*, 'Not a proper illumination method'
			call exit
		end if
		
		!open (unit = 203, file = 'input_field.txt', action = "write",status = 'replace')
		!	do m = 1, n1*n2        
		!		write (203, '(1001(E19.12, tr3))') 	(real(field_input(m)%e_field%vector(n)), n= 1, 3), &
		!		(imag(field_input(m)%e_field%vector(n)), n= 1, 3),	(real(field_input(m)%h_field%vector(n)), n= 1, 3), &
		!		(imag(field_input(m)%h_field%vector(n)), n= 1, 3)	
		!	end do
		!close(203)
		!deallocate(r_local)
		return
	end subroutine
	
	function ngp_near_field_distance_tri(r_loc, distance_nf)result(ngp)
		integer :: ngp
		integer :: p, ngp_arr(number_objects)		
		real(dp) :: r_loc(3)
		real(dp), dimension(2), intent(in) :: distance_nf
		real(dp) :: tmp_r		
		
		!calc_p(5) = 4, 1
		!if (pre_types%object .eq. 'sphere')then
		if (calc_p(5) .eq. 1)then !sphere from SpeckleSim
			do p = 1, number_objects
				tmp_r = (vec_len(r_loc - tri_sp_parameters(p)%centroid%point) - tri_sp_parameters(p)%D*0.5)		
				if ((tmp_r .le. distance_nf(1)) .and. (tmp_r .gt. 0.0) ) then
					ngp_arr(p) = 12
				else if ((tmp_r .le. distance_nf(2)) .and. (tmp_r .gt. distance_nf(1))) then
					ngp_arr(p) = 6
				else if (tmp_r .gt. distance_nf(2)) then				
					ngp_arr(p) = 3
				else 
					ngp_arr(p) = 6
				end if
				ngp = maxval(ngp_arr)
			end do
		else if (calc_p(5) .eq. 4)then !sphere from COMSOL			
				p=1
				tmp_r = (vec_len(r_loc - tri_sf_parameters(p)%centroid%point) - tri_sf_parameters(p)%D(1)*0.5)
				if ((tmp_r .le. distance_nf(1)) .and. (tmp_r .gt. 0.0) ) then
					ngp_arr(p) = 12
				else if ((tmp_r .le. distance_nf(2)) .and. (tmp_r .gt. distance_nf(1))) then
					ngp_arr(p) = 6
				else if (tmp_r .gt. distance_nf(2)) then				
					ngp_arr(p) = 3
				else 
					ngp_arr(p) = 6
					!print*, 'ngp=6', ngp_arr(p)
				end if
				ngp = maxval(ngp_arr)			
		!else if (pre_types%object .eq. 'surface')then		
		else if ((calc_p(5) .eq. 2) .or. (calc_p(5) .eq. 3))then			
			ngp = 3
		else
			ngp = 3
			!print*, 'not a properly defined object type'			
		end if		
		return
	end function
	
	!Just using more integration points	
	subroutine singular_integration_edge_adjacent_ngp(fm, fn, sum_a)
		implicit none
		complex(dp), intent(out) :: sum_a
		type(fn_rwg), intent(in) :: fm, fn
		
		!dummy
		real(dp) :: TRF
		real(dp), dimension(100) :: w(100), a(100), b(100)
		real(dp) :: Rq(3), Ra(3), r_fm(3), r_fn(3), Rn(3), R(3) 
		complex(dp) :: DG_1, DG_2, G_1, G_2
		integer :: i, j, ngp
		
		ngp = 12		
		sum_a = (0.0, 0.0)	
		call Quadrature_tri(ngp, a, b, w) 
		do i = 1, ngp
			Ra = r_simplex_new(a(i), b(i), fm%corners)
			r_fm = (Ra-fm%corners(3)%point)*fm%df*0.5
			do j = 1, ngp 
				Rq = r_simplex_new(a(j), b(j), fn%corners)
				r_fn = (Rq-fn%corners(3)%point)*fn%df*0.5
				R = Ra - Rq				
				Rn = R/vec_len(R)
				call Diver_Green(k1, vec_len(R), G_1, DG_1)
				call Diver_Green(k2, vec_len(R), G_2, DG_2)						
				TRF = w(i)*w(j)*fn%area*fm%area
				sum_a = sum_a + dot_product(r_fm, cross_r(r_fn, Rn)*(DG_1 + DG_2))*TRF !		
			end do !loop j
		end do  !loop i
	end subroutine singular_integration_edge_adjacent_ngp
	
	subroutine edge_vectors_assigned(corner_fm, edge_fm)
		type(edge_tri), dimension(3), intent(out) :: edge_fm
		type(point), intent(in) :: corner_fm(3)
		edge_fm(1)%p_s = corner_fm(3)%point
		edge_fm(1)%p_t = corner_fm(2)%point
		
		edge_fm(2)%p_s = corner_fm(1)%point
		edge_fm(2)%p_t = corner_fm(3)%point
		
		edge_fm(3)%p_s = corner_fm(2)%point
		edge_fm(3)%p_t = corner_fm(1)%point
	end subroutine
	
	subroutine line_integration_nodes(edge_fm, ngp, xx, ww)
		implicit none
		type(edge_tri), intent(in) :: edge_fm
		integer, intent(in) :: ngp
		real(dp), dimension(:), allocatable, intent(out) ::xx, ww
		
		!dummy
		real(dp) :: sp, sn, s_vec(3)
		real(dp) :: p1(3), p2(3)		
		real(dp), dimension(:), allocatable :: xx_t, ww_t 
		
		p1 = edge_fm%p_s
		p2 = edge_fm%p_t		
		s_vec = (p2-p1)/vec_len(p2-p1)		
		sp = dot_product(p2, s_vec)
		sn = dot_product(p1, s_vec)	
		call legendre_handle(ngp, sn, sp, xx_t, ww_t)
		xx = xx_t
		ww = ww_t
	end subroutine

	subroutine Coefficient_D11s(P, c_arr)
		real(dp), intent(in) :: P
		complex(dp), dimension(4), intent(out) :: c_arr
		c_arr(1) = Im*Omega/(4*PI)*(my_1 + my_2)
		c_arr(2) = -Im*P/(4*PI)*(1/eps_1 + 1/eps_2)
		c_arr(3) = -Im*Omega/(8*PI)*(my_1*k1**2 + my_2*k2**2)
		c_arr(4) = Im*P/(8*PI)*(k1**2/eps_1 + k2**2/eps_2)
		return
	end subroutine Coefficient_D11s
		
	subroutine Coefficient_D22s(Pc, c_arr)
		real(dp), intent(in) :: Pc		
		complex(dp), dimension(4), intent(out) :: c_arr
		c_arr(1) = Im*Omega/(4*PI)*(eps_1 + eps_2)
		c_arr(2) = -Im*Pc/(4*PI)*(1/my_1 + 1/my_2)
		c_arr(3) = -Im*Omega/(8*PI)*(eps_1*k1**2 + eps_2*k2**2)
		c_arr(4) = Im*Pc/(8*PI)*(k1**2/my_1 + k2**2/my_2)
		return
	end subroutine Coefficient_D22s			
   
	function r_simplex(alpha, beta, r1, r2, r3)
      implicit none
      real(dp), intent(in) :: r1(3), r2(3), r3(3), alpha, beta       
      real(dp) :: r_simplex(3)
      r_simplex(1) = (1-alpha - beta)*r1(1) + alpha*r2(1) + beta*r3(1)
      r_simplex(2) = (1-alpha - beta)*r1(2) + alpha*r2(2) + beta*r3(2)
      r_simplex(3) = (1-alpha - beta)*r1(3) + alpha*r2(3) + beta*r3(3)
      return    
	end function 
	 
	function r_simplex_new(alpha, beta, corner_fn)result(rv)
      implicit none
      real(dp), intent(in) :: alpha, beta       
		type(point), intent(in) :: corner_fn(3)
      real(dp) :: rv(3)
      rv(1) = (1-alpha - beta)*corner_fn(1)%point(1) + alpha*corner_fn(2)%point(1) + beta*corner_fn(3)%point(1)
      rv(2) = (1-alpha - beta)*corner_fn(1)%point(2) + alpha*corner_fn(2)%point(2) + beta*corner_fn(3)%point(2)
      rv(3) = (1-alpha - beta)*corner_fn(1)%point(3) + alpha*corner_fn(2)%point(3) + beta*corner_fn(3)%point(3)
      return    
	end function 
	
	!Calculation of the scattered fields!
	subroutine Integration_scattered_field_tri(r_a, struct, ngp, SEH, Sum_a)
      complex(dp), dimension(:), intent(in) ::  SEH        
      real(dp), intent(in) :: r_a(3)
		type(structure_tri), intent(in) :: struct
		integer, intent(in) :: ngp
		
      type(Vector_c), dimension(2), intent(out) :: Sum_a
		
      integer i, t, pos_neg, m
      real(dp) :: R_norm, Len_m, area_m, dfm 
		complex(dp) :: GradG(3), K_tmp(3), Greenf, DG_SS, f_mat(3)
		real(dp) :: r_q(3), p1(3), p2(3), vm_t(3)
		real(dp), dimension(100) :: a, b, w	
		integer :: n_object, m_tmp
		
      call Quadrature_tri(ngp, a, b, w)
      do t = 1, 2
         Sum_a(t)%Vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		end do  
		
		do m = 1, m_pairs
			do t = 1, 2				
				call fn_parameter(struct, m, t, p1, p2, vm_t, Len_m, area_m)
				pos_neg = -1*(((t-1)*t)-1)            
				dfm = pos_neg*(Len_m)	!/area_m
				do i = 1, ngp
					r_q = r_simplex(a(i), b(i), p1, p2, vm_t)
					f_mat=0.5*dfm*(r_q-vm_t) !*area_m
					R_norm = vec_len(r_a - r_q)
					call Diver_Green(k1, R_norm, Greenf, DG_SS)
					GradG = (r_a - r_q)/R_norm*DG_SS					
      
					K_tmp = im*Omega*my_1*SEH(m)*f_mat*Greenf - &
						cross_c(SEH(m+m_pairs)*f_mat, GradG) - 1/(im*Omega*eps_1)*SEH(m)*dfm*GradG*0
					Sum_a(1)%Vector = Sum_a(1)%Vector - K_tmp*w(i)!TRF
					K_tmp = im*Omega*eps_1*SEH(m+m_pairs)*f_mat*Greenf &
						- 1/(im*Omega*my_1)*SEH(m+m_pairs)*dfm*GradG + SEH(m)*cross_c(f_mat, GradG)
					Sum_a(2)%Vector = Sum_a(2)%Vector - K_tmp*w(i)!TRF
				end do !loop i
			end do !loop t						
		end do
		!print*, 'Efield, Sum_a(1) =', Sum_a(1)
      return
	end subroutine Integration_scattered_field_tri	
	
subroutine Integration_scattered_field_tri_within_medium(r_a, struct, ngp, SEH, Sum_a)
      complex(dp), dimension(:), intent(in) ::  SEH 
      !type(point), intent(in) :: r_a
		type(structure_tri), intent(in) :: struct
		integer, intent(in) :: ngp
		type(evaluation_r_media), intent(in):: r_a
      type(Vector_c), dimension(2), intent(out) :: Sum_a
		
      integer i, t, pos_neg, mm
      real(dp) :: R_norm, Len_m, area_m, dfm 
		complex(dp) :: GradG(3), Ke_tmp(3), Kh_tmp(3), Greenf, DG_SS, f_mat(3)
		complex(dp) :: SE, SH, km, eps_r
		real(dp) :: r_q(3), p1(3), p2(3), vm_t(3)
		real(dp), dimension(100) :: a, b, w	
		integer :: n_object, m_tmp	
				
      call Quadrature_tri(ngp, a, b, w)
		
      do t = 1, 2
         Sum_a(t)%Vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      end do  
		eps_r =	eps_0*r_a%eps_r
		km = k0*sqrt(r_a%eps_r)
		do mm = 1, m_pairs		
			SE = -1.0*SEH(mm)
			SH = -1.0*SEH(mm + m_pairs)
			do t = 1, 2
				call fn_parameter(struct, mm, t, p1, p2, vm_t, Len_m, area_m)				
				pos_neg = -1*(((t-1)*t)-1) 
				dfm = pos_neg*(Len_m)	!/area_m
				do i = 1, ngp
					r_q = r_simplex(a(i), b(i), p1, p2, vm_t)
					f_mat=0.5*dfm*(r_q-vm_t) !*area_m
					R_norm = vec_len(r_a%point - r_q)
					call Diver_Green(km, R_norm, Greenf, DG_SS)
					GradG = (r_a%point - r_q)/R_norm*DG_SS					
					
					Ke_tmp = im*Omega*my_0*SE*f_mat*Greenf - &
						cross_c(SH*f_mat, GradG) - 1/(im*Omega*eps_r)*SE*dfm*GradG	
					Sum_a(1)%Vector = Sum_a(1)%Vector - Ke_tmp*w(i)!TRF				
					
					Kh_tmp = im*Omega*eps_r*SH*f_mat*Greenf &
						- 1/(im*Omega*my_0)*SH*dfm*GradG + SE*cross_c(f_mat, GradG)
					Sum_a(2)%Vector = Sum_a(2)%Vector - Kh_tmp*w(i)!TRF
				end do            
			end do !loop i 
		end do
      return
	end subroutine Integration_scattered_field_tri_within_medium	
	 
 end module lib_sie_tri_calculation_mod