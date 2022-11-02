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
	
module lib_sie_quad_data_container	
	use lib_sie_type_function
	use lib_sie_data_container
		
	implicit none
	
	!surface related variables
	integer :: Np, Nt
	real(dp), dimension(:, :), allocatable :: pp
	integer, dimension(:, :), allocatable :: tt
	real(dp) :: surface_length_quad(2)
	type(quad_surface_parameter), dimension(:), allocatable :: quad_sf_parameters !
	type(structure_quad), dimension(:), allocatable :: struc_quad_surfaces
	
		
	type(quad_sphere_parameter), dimension(:), allocatable :: quad_sp_parameters !
	type(structure_quad), dimension(:), allocatable :: struc_quad_spheres
	
	!common variables
	character(len = 20) :: Result_SEH	
	type(structure_quad), public :: struc_quad			
	type(edges_parameter), dimension(:), allocatable :: edge_I
	type(element_edges), dimension(:), allocatable :: edge_II
	
	character(len = 50) :: File_NodeVector
	character(len = 20) :: File_ElementIndicies
	
	contains
	
	subroutine set_parameters_quad()				
		call set_global_parameters_quad()		
		call set_local_material_parameters_quad()
	end subroutine
	
	subroutine set_global_parameters_quad()
		k0 = 2*PI/illumination_p%lambda
		omega = 2*PI*c0/illumination_p%lambda 		
		my_2 = my_0
		my_1 = my_0
		number_objects = 1			
		if (pre_types%object .eq. 'sphere') then			
			allocate(quad_sp_parameters(number_objects))
			allocate(struc_quad_spheres(number_objects))					
			call set_sphere_parameters_quad()
			eps_r2 = quad_sp_parameters(1)%eps_r2
			eps_r1 = quad_sp_parameters(1)%eps_r1
			call set_local_material_parameters_quad()			
			
		else if (pre_types%object .eq. 'surface') then	
			allocate(quad_sf_parameters(number_objects))	
			allocate(struc_quad_surfaces(number_objects))
			call set_surface_parameters_quad()	
			eps_r2 = quad_sf_parameters(1)%eps_r2
			eps_r1 = quad_sf_parameters(1)%eps_r1
			call set_local_material_parameters_quad()
		else 
			print*, 'Not a proper object type.'
		end if
	end subroutine
	
	subroutine set_surface_parameters_quad()
		implicit none
		integer :: m
		type(point) :: shift
				
		call sorting_quad()
		
		shift%point = (/0.0, 0.0, 1300.0e-9/)		
		do m = 1, number_objects			
			quad_sf_parameters(m)%D(1:2) = geometry_p(1) !
			quad_sf_parameters(m)%D(3) = 0.0
			quad_sf_parameters(m)%centroid%point = (m-1)*shift%point
			quad_sf_parameters(m)%np = Np 
			quad_sf_parameters(m)%nt =	Nt 
			quad_sf_parameters(m)%eps_r2 =  eps_r2_main !(-16.075, -0.442342)! (2.465, 0.0)! (-17.2, -0.498)!
			quad_sf_parameters(m)%eps_r1 = (1.0, 0.0)!(-17.2, -0.498)!
		end do	
	end subroutine
		
	subroutine set_sphere_parameters_quad()
		implicit none
		integer :: m
		type(point) :: shift
		character(len = 20) :: str1
		
		shift%point = (/1500.0e-9, 0.0, 0.0/)	
		
		do m = 1, number_objects
			quad_sp_parameters(m)%D = 400.0e-9
			quad_sp_parameters(m)%centroid%point = shift%point*(m-1)
			quad_sp_parameters(m)%nt =32!   168! 512!  2688 !  672  ! 128!
			quad_sp_parameters(m)%np =98!   506! 1538! 8066 !  2018 ! 386!
			quad_sp_parameters(m)%eps_r2 = eps_r2_main !(-16.075, -0.442342)! 
			quad_sp_parameters(m)%eps_r1 = (1.0, 0.0)
			
			write(str1,'(I5)') quad_sp_parameters(m)%Nt
			File_NodeVector = 'E'//trim(adjustl(str1))//'_ppQua.txt'			
			File_ElementIndicies = 'E'//trim(adjustl(str1))//'_ttQua.txt'				
		end do
	end subroutine	
	
	subroutine set_local_material_parameters_quad()
		implicit none		
		k0 = 2*PI/illumination_p%lambda
		omega = 2*PI*c0/illumination_p%lambda 		
		my_2 = my_0
		my_1 = my_0		
		k1 = k0*sqrt(eps_r1)
		k2 = k0*sqrt(eps_r2)
		eps_1 = eps_0*eps_r1
		eps_2 = eps_0*eps_r2
		eta_1 = sqrt(my_1/eps_1)
		eta_2 = sqrt(my_2/eps_2)
	end subroutine
	
	subroutine set_evaluation_parameters_quad()
		implicit none
		real(dp) :: phi, r_Far
		
		evaluation_parameter%tot_field = total_field
		
		if ((pre_types%evaluation .eq. 'rcs_p') .or.  (pre_types%evaluation .eq. 'rcs_n') .or. &
			(pre_types%evaluation .eq. 'BRDF_p') .or. (pre_types%evaluation .eq. 'BRDF_n')) then			
			select case (pre_types%evaluation)
				case('rcs_p')
					phi = 0 !
				case('rcs_n')
					phi = PI/2 !
				case('BRDF_p')	
					phi = 0 !
				case('BRDF_n')	
					phi = PI/2 !					
			end select
		
			r_far = 5.0e+9
			!theta
			evaluation_parameter%dim_a(1) = theta_start
			evaluation_parameter%dim_a(2) = theta_end
			evaluation_parameter%dim_b(1) = phi! phi_start
			evaluation_parameter%dim_b(2) = phi! phi_end			
			evaluation_parameter%N_dim(1) = 301
			evaluation_parameter%N_dim(2) = 1
			if (pre_types%object .eq. 'sphere') then
				evaluation_parameter%dim_c = sum(quad_sp_parameters(1:number_objects)%D)*r_Far !
			else if (pre_types%object .eq. 'surface') then				
				evaluation_parameter%dim_c = (sum(quad_sf_parameters(1:number_objects)%D(1)) &
				+ sum(quad_sf_parameters(1:number_objects)%D(2)))*r_Far !
			else 
				print*, 'Wrong calculation type'
				call exit				
			end if
		else
			evaluation_parameter%N_dim(1) = 201
			evaluation_parameter%N_dim(2) = 201		
			if (pre_types%object .eq. 'sphere') then
				evaluation_parameter%dim_a(1) = -sum(quad_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_a(2) =  sum(quad_sp_parameters(1:number_objects)%D)*2.5		
				evaluation_parameter%dim_b(1) = -sum(quad_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_b(2) =  sum(quad_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_c = 0.0				
			else if (pre_types%object .eq. 'surface') then
				evaluation_parameter%dim_a(1) = -sum(quad_sf_parameters(1:number_objects)%D(1))*2.5
				evaluation_parameter%dim_a(2) =  sum(quad_sf_parameters(1:number_objects)%D(1))*2.5		
				evaluation_parameter%dim_b(1) = -sum(quad_sf_parameters(1:number_objects)%D(2))*2.5
				evaluation_parameter%dim_b(2) =  sum(quad_sf_parameters(1:number_objects)%D(2))*2.5
				evaluation_parameter%dim_c = 0.0 !position of the third dimension					
			else 
				print*, 'Not a well defined object_type!'
				call exit 
			end if	
		end if
	end subroutine
	
	subroutine data_import_quad()    
		implicit none    		
		integer, dimension(:, :), allocatable :: Elements_indicies
		real(dp), dimension(:, :), allocatable :: Node_Vector
		integer :: j, m, q, p
		
		select case (pre_types%object)		
		case ('sphere')
			allocate(struc_quad%elements(sum(quad_sp_parameters(1:number_objects)%Nt)))
			allocate(struc_quad%node_vector(sum(quad_sp_parameters(1:number_objects)%Np)))
			print*, 'node numbers:', sum(quad_sp_parameters(1:number_objects)%Np)
			p = 0 
			q = 0
			do m = 1, number_objects				
				allocate(node_vector(1:4, 1:quad_sp_parameters(m)%Np)) 
				allocate(elements_indicies(1:10, 1:quad_sp_parameters(m)%Nt))			
				
				open(unit = 80, file = File_NodeVector, status = 'old', action='read') !Small_
				read(80, *) node_vector
				close(80)
				
				!Array of node index for each element
				open(unit = 81, file = File_ElementIndicies, status = 'old', action='read') !Small_
				read(81, *) elements_indicies
				close(81) 
				
				if (m .gt. 1)then
					p = p + quad_sp_parameters(m-1)%Np				
					do j = 1, quad_sp_parameters(m)%Np !pp						
						struc_quad%node_vector(j + p)%point = node_vector(2:4, j)*quad_sp_parameters(m)%D + &
								quad_sp_parameters(m)%centroid%point !for sphere, here is 2:4
					end do	
				else
					do j = 1, quad_sp_parameters(m)%Np !pp
						struc_quad%node_vector(j)%point = node_vector(2:4, j)*quad_sp_parameters(m)%D + &
								quad_sp_parameters(m)%centroid%point !for sphere, here is 2:4
					end do
				end if
								
				if (m .gt. 1) then
					q = q + quad_sp_parameters(m-1)%Nt
					do j = 1, quad_sp_parameters(m)%Nt !tt
					   allocate(struc_quad%elements(j + q)%vertices(8))
						struc_quad%elements(j + q)%vertices = elements_indicies(3:10, j) + p
					end do
				else 	
					do j = 1, quad_sp_parameters(m)%Nt !tt				
						 allocate(struc_quad%elements(j)%vertices(8))
						struc_quad%elements(j)%vertices = elements_indicies(3:10, j) 					   	
					end do				
				end if	
				
				deallocate(node_vector)
				deallocate(elements_indicies)				
			end do
			!file_name = 'three_quad_sphere_tt.txt'
			!print*, 'file_name =', file_name
   !
			!open (unit=206, file = file_name, action="write",status = 'replace')
			!	do j = 1, sum(quad_sp_parameters(1:number_objects)%Nt)
			!		write (206, '(10(i10, tr5))') j, 2, struc_quad%elements(j)%vertices
			!	end do
			!close(206)
   !
			!file_name ='three_quad_sphere_pp.txt'
			!
			!open (unit=206, file = file_name, action="write",status = 'replace')			
			!	do j = 1, sum(quad_sp_parameters(1:number_objects)%Np)
			!		write (206, '(3(es19.12, tr5))') struc_quad%node_vector(j)
			!	end do
			!close(206)	
		case ('surface')
			allocate(struc_quad%elements(sum(quad_sf_parameters(1:number_objects)%Nt)))
			allocate(struc_quad%node_vector(sum(quad_sf_parameters(1:number_objects)%Np)))
			print*, 'node numbers:', sum(quad_sf_parameters(1:number_objects)%Np)
			p = 0 
			q = 0
			do m = 1, number_objects
		
				allocate(node_vector(1:4, 1:quad_sf_parameters(m)%Np)) 
				allocate(elements_indicies(1:10, 1:quad_sf_parameters(m)%Nt))			
				
				do p = 1, quad_sf_parameters(m)%Np
					node_vector(1, p) = p					
				end do
				node_vector(2:4, :) = pp(1:3, :)				
				elements_indicies = tt
				
				!open(unit = 80, file = File_NodeVector, status = 'old', action='read') !Small_
				!read(80, *) node_vector
				!close(80)
    !   
				!!Array of node index for each element
				!open(unit = 81, file = File_ElementIndicies, status = 'old', action='read') !Small_
				!read(81, *) elements_indicies
				!close(81) 				
				
				if (m .gt. 1)then
					p = p + quad_sf_parameters(m-1)%Np				
					do j = 1, quad_sf_parameters(m)%Np !pp						
						struc_quad%node_vector(j + p)%point = node_vector(2:4, j) + &
								quad_sf_parameters(m)%centroid%point !for sphere, here is 2:4
					end do	
				else
					do j = 1, quad_sf_parameters(m)%Np !pp
						struc_quad%node_vector(j)%point = node_vector(2:4, j) + &
								quad_sf_parameters(m)%centroid%point !for sphere, here is 2:4
					end do
				end if
								
				if (m .gt. 1) then
					q = q + quad_sf_parameters(m-1)%Nt
					do j = 1, quad_sf_parameters(m)%Nt !
					   allocate(struc_quad%elements(j + q)%vertices(8))
						struc_quad%elements(j + q)%vertices = elements_indicies(3:10, j) + p
					end do
				else 	
					do j = 1, quad_sf_parameters(m)%Nt !tt		
						allocate(struc_quad%elements(j)%vertices(8))
						struc_quad%elements(j)%vertices = elements_indicies(3:10, j) 
					end do
				end if					
			
				deallocate(node_vector)
				deallocate(elements_indicies)				
			end do		
		end select
		
		return		
   end subroutine data_import_quad
	
	subroutine sorting_quad()
		implicit none
		
		integer, dimension(:, :), allocatable :: tt_temp, ss 
		integer :: i, j
		integer :: Nl(2), N_element(2)
		real(dp), dimension(:, :), allocatable :: zz 
		real :: dL(2), n_tmp(2)
		character(len = 60) :: file_c
	
		open(unit = 92, file = file_name_surface, status = 'old', action='read')
			read(92, *) surface_length_quad(1:2)
			read(92, *) n_tmp(1:2)
		close(92)
		
		n_element = int(n_tmp)
		
		Nl = N_element*2 + 1 
				
		allocate(zz(1: Nl(1), 1: Nl(2)))
	
		open(unit = 92, file = file_name_surface, status = 'old', action='read')
			read(92, *) surface_length_quad(1:2)
			read(92, *) n_tmp(1:2)			
			read(92, *) ((zz(i, j), i = 1, Nl(1)),j = 1,Nl(2))
		close(92)
		 
		dL(1) = surface_length_quad(1)/(Nl(1)-1)
		dL(2) = surface_length_quad(2)/(Nl(2)-1)
		
		Np = 0
		Nt = 0
		Nl = N_element*2 + 1

		allocate (ss(Nl(1), Nl(2)))
		allocate (pp(3, Nl(1)*Nl(2)))
 
		do i = 1, Nl(1)
			do j = 1, Nl(2) 
				Np = Np + 1 ! all the points
				ss(i, j) = Np      
				pp(1, Np) = -surface_length_quad(1)/2+dL(1)*(i-1)
				pp(2, Np) = -surface_length_quad(2)/2+dL(2)*(j-1)
				pp(3, Np) = zz(i, j) 
			end do
		end do   
		
		allocate (tt_temp(1:10, 1:Np))
	
		do j = 1, Nl(1)-2, 2
			do i = 1, Nl(2)-2, 2      
				Nt = Nt + 1 ! for tt     
				tt_temp(1, Nt) = Nt
				tt_temp(2, Nt) = 2
				tt_temp(3, Nt) = ss(j, i) !the first node 
				tt_temp(4, Nt) = ss(j+1, i) !the second
				tt_temp(5, Nt) = ss(j+2, i) 
				tt_temp(6, Nt) = ss(j+2, i+1)
				tt_temp(7, Nt) = ss(j+2, i+2)
				tt_temp(8, Nt) = ss(j+1, i+2)
				tt_temp(9, Nt) = ss(j, i+2)
				tt_temp(10, Nt) = ss(j, i+1) 
			end do 
		end do
		
		allocate(tt(1:10, 1:Nt))
 
		tt(1:10, 1:Nt) = tt_temp(1:10, 1:Nt)
		
		file_c = trim('tt_quad_'//adjustl(file_name_surface))		
		open (unit = 203, file = file_c, action = "write",status = 'replace')
			write (203, '(10(i8, tr3))') (tt(:, j), j = 1, Nt)   
		close(203)
  !
		file_c = trim('pp_quad_'//adjustl(file_name_surface))
		open (unit=200,file =	file_c, action="write",status = 'replace')!   
			write (200, '(i5, tr5, e19.12, tr5, e19.12, tr5, e19.12, tr5)') (j, pp(:, j), j = 1, Np)   
		close(200)		
		
		deallocate(zz, tt_temp, ss)		
		return
		
	end subroutine Sorting_quad
	
	subroutine sorting_quad_new()
		implicit none
		
		integer, dimension(:, :), allocatable :: tt_temp, ss 
		integer :: i, j
		integer :: Nl(2), N_element(2), Nop
		real(dp), dimension(:, :), allocatable :: zz 
		real :: n_tmp(2)
		character(len = 60) :: file_c
	
		open(unit = 92, file = file_name_surface, status = 'old', action='read')
			read(92, *) surface_length_quad(1:2)
			read(92, *) n_tmp(1:2)
		close(92)
		
		Nl = int(n_tmp)		
		Nop = Nl(2)*Nl(1)		
		allocate(zz(3, Nop))
	
		open(unit = 92, file = file_name_surface, status = 'old', action='read')
			read(92, *) surface_length_quad(1:2)
			read(92, *) n_tmp(1:2)			
			read(92, *) ((zz(i, j), i = 1, 3),j = 1, Nop)
		close(92)
				
		Np = 0
		Nt = 0
		
		allocate (ss(Nl(1), Nl(2)))
		allocate (pp(3, Nl(1)*Nl(2)))
 
		do i = 1, Nl(2) !y
			do j = 1, Nl(1) !x 
				Np = Np + 1 ! all the points
				ss(j, i) = Np      				
				pp(1, Np) = zz(1, Np) 
				pp(2, Np) = zz(2, Np) 
				pp(3, Np) = zz(3, Np) 
			end do
		end do   
		
		allocate (tt_temp(1:10, 1:Np))
	
		do j = 1, Nl(1)-2, 2
			do i = 1, Nl(2)-2, 2      
				Nt = Nt + 1 ! for tt     
				tt_temp(1, Nt) = Nt
				tt_temp(2, Nt) = 2
				tt_temp(3, Nt) = ss(j, i) !the first node 
				tt_temp(4, Nt) = ss(j+1, i) !the second
				tt_temp(5, Nt) = ss(j+2, i) 
				tt_temp(6, Nt) = ss(j+2, i+1)
				tt_temp(7, Nt) = ss(j+2, i+2)
				tt_temp(8, Nt) = ss(j+1, i+2)
				tt_temp(9, Nt) = ss(j, i+2)
				tt_temp(10, Nt) = ss(j, i+1) 
			end do 
		end do
		
		allocate(tt(1:10, 1:Nt))
 
		tt(1:10, 1:Nt) = tt_temp(1:10, 1:Nt)
		
		file_c = trim('tt_quad_'//adjustl(file_name_surface))		
		open (unit = 203, file = file_c, action = "write",status = 'replace')
			write (203, '(10(i8, tr3))') (tt(:, j), j = 1, Nt)   
		close(203)
  !
		file_c = trim('pp_quad_'//adjustl(file_name_surface))
		open (unit=200,file =	file_c, action="write",status = 'replace')!   
			write (200, '(i5, tr5, e19.12, tr5, e19.12, tr5, e19.12, tr5)') (j, pp(:, j), j = 1, Np)   
		close(200)		
		
		deallocate(zz, tt_temp, ss)		
		return
		
	end subroutine Sorting_quad_new
	
end module lib_sie_quad_data_container