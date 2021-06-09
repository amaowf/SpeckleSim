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

module lib_sie_tri_data_container
	use ml_fmm_type
	
	use lib_sie_math
	use lib_sie_constants
	use lib_sie_data_container
	use lib_sie_type_function
	
	implicit none	
	
	integer, dimension(:), allocatable :: m_pairs_arr
		
	type(tri_sphere_parameter), dimension(:), allocatable, public :: tri_sp_parameters !	
	type(tri_surface_parameter),dimension(:), allocatable, public :: tri_sf_parameters !
	type(structure_tri), dimension(:), allocatable, public :: struc_tri_spheres
	type(structure_tri), dimension(:), allocatable, public :: struc_tri_surfaces	
	
	type(structure_tri) :: struc_tri	
	real(dp) :: Rfx, Rfy
	
	
	real(dp), dimension(3) :: centroid	
	type(edges_in_structure), dimension(:), allocatable :: structure_edges	
	integer, dimension(:), allocatable :: np_arr, ne_arr!
	
	contains
	
	subroutine set_parameters_tri()		
		call set_global_parameters_tri()		
		call set_local_material_parameters_tri()	
	end subroutine
	
	subroutine set_surface_parameters_tri()
		implicit none
		integer :: m
		type(point) :: shift
		
		shift%point = (/0.0, 0.0, 800.0e-9/)
		do m = 1, number_objects			
			call read_pp_tt_file_tri(struc_tri_surfaces(m))
			tri_sf_parameters(m)%D(1) = surface_length
			tri_sf_parameters(m)%D(2) = surface_length
			tri_sf_parameters(m)%centroid%point = (m-1)*shift%point	
			tri_sf_parameters(m)%eps_r2 = eps_r2_main !(2.25, 0.0)! (-16.075, -0.442342)!(-17.2, -0.498)! 
			tri_sf_parameters(m)%eps_r1 = (1.0, 0.0)
		end do	
		!eps_r2 = (-17.2, -0.498) !(-8.0, - 1.66) !(2.50, 0.0)!
	end subroutine
	
	subroutine set_sphere_parameters_tri()
		implicit none
		integer :: m
		type(point) :: shift
		
		shift%point = (/1000.0e-9, 0.0, 0.0/)
		do m = 1, number_objects		
			tri_sp_parameters(m)%D = 500e-9
			tri_sp_parameters(m)%centroid%point = shift%point*(m-1)
			tri_sp_parameters(m)%eps_r2 = eps_r2_main! (-16.075, -0.442342)!(2.465, 0.0)! 
			tri_sp_parameters(m)%eps_r1 = (1.0, 0.0)
			!n_disc = 4x(4x6-6)-6 
			!n_disc = 4x(4x(4x6-6)-6)-6....
			!n_disc = 6, 18, 66, 1026, 4098, 16386, 65538, 80500 
			!n_disc = 66(nel=128), 258(nel = 512), 1026(nel = 2048)
			!n_disc = 4098(nel = 8192), 16386(nel=32768)
			tri_sp_parameters(m)%n_disc =  258! 4098! 66! 65538 !,  16386!  16386 ! 1026			
		end do	
		
		!when the spheres have different size and disretization
		!tri_sp_parameters(1)%D = 400e-9
		!tri_sp_parameters(2)%D = 600e-9
		!tri_sp_parameters(3)%D = 800e-9
		!tri_sp_parameters(4)%D = 600e-9
		!tri_sp_parameters(5)%D = 400e-9		
		!tri_sp_parameters(1)%centroid%point = (/1400.0e-9, 0.0, 0.0/) 
		!tri_sp_parameters(2)%centroid%point = (/800.0e-9, 0.0, 0.0/) 
		!tri_sp_parameters(3)%centroid%point = (/0.0, 0.0, 0.0/) 
		!tri_sp_parameters(4)%centroid%point = (/-800.0e-9, 0.0, 0.0/) 
		!tri_sp_parameters(5)%centroid%point = (/-1400.0e-9, 0.0, 0.0/) 		
		!
		!tri_sp_parameters(1)%n_disc = 258 !
		!tri_sp_parameters(2)%n_disc = 258 !
		!tri_sp_parameters(3)%n_disc = 1026 !
		!tri_sp_parameters(4)%n_disc = 258 !
		!tri_sp_parameters(5)%n_disc = 258 !
		
		
	end subroutine
	
	subroutine set_evaluation_parameters_tri()	
		implicit none		
		real(dp) :: phi
		
		evaluation_parameter%tot_field = total_field
		
		if ((pre_types%evaluation .eq. 'rcs_p') .or.  (pre_types%evaluation .eq. 'rcs_n') .or. (pre_types%evaluation .eq. 'BRDF')) then			
			select case (pre_types%evaluation)
				case('rcs_p')
					phi = 0 !
				case('rcs_n')
					phi = PI/2 !
				case('BRDF')	
					phi = 0 !
			end select
			
			!theta
			evaluation_parameter%dim_a(1) = theta_start!
			evaluation_parameter%dim_a(2) = theta_end!
			evaluation_parameter%dim_b(1) = phi_start!
			evaluation_parameter%dim_b(2) = phi_end!			
			evaluation_parameter%N_dim(1) = 301 !for theta
			evaluation_parameter%N_dim(2) = 1!for phi
			if (pre_types%object .eq. 'sphere') then
				evaluation_parameter%dim_c = -sum(tri_sp_parameters(1:number_objects)%D)*5e+9 !
			else if (pre_types%object .eq. 'surface') then				
				evaluation_parameter%dim_c = -(sum(tri_sf_parameters(1:number_objects)%D(1)) &
				+ sum(tri_sf_parameters(1:number_objects)%D(2)))*5e+9 !
			else 
				print*, 'Wrong calculation type'
				call exit				
			end if	
		
		else
			evaluation_parameter%N_dim(1) = 201
			evaluation_parameter%N_dim(2) = 201
			if (pre_types%object .eq. 'sphere') then
				evaluation_parameter%dim_c = 0.0 !the position of the third dimension
				evaluation_parameter%dim_a(1) = -sum(tri_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_a(2) = sum(tri_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_b(1) = -sum(tri_sp_parameters(1:number_objects)%D)*2.5
				evaluation_parameter%dim_b(2) = sum(tri_sp_parameters(1:number_objects)%D)*2.5
			else if (pre_types%object .eq. 'surface') then					
				evaluation_parameter%dim_c = 0.0e-6 !the position of the third dimension
				evaluation_parameter%dim_a(1) = -sum(tri_sf_parameters(1:number_objects)%D(1))*2.5
				evaluation_parameter%dim_a(2) = sum(tri_sf_parameters(1:number_objects)%D(1))*2.5
				evaluation_parameter%dim_b(1) = -sum(tri_sf_parameters(1:number_objects)%D(2))*2.5
				evaluation_parameter%dim_b(2) = sum(tri_sf_parameters(1:number_objects)%D(2))*2.5
			else
				print*, 'Wrong calculation type!'
				call exit
			end if
		end if		
	end subroutine
	
	subroutine set_global_parameters_tri()
		implicit none
		
		k0 = 2*PI/illumination_p%lambda
		omega = 2*PI*c0/illumination_p%lambda 		
		my_2 = my_0
		my_1 = my_0
		number_objects = 1			
		if (pre_types%object .eq. 'sphere') then			
			allocate(tri_sp_parameters(number_objects))
			allocate(struc_tri_spheres(number_objects))
			call set_sphere_parameters_tri()	
					
			eps_r2 = tri_sp_parameters(1)%eps_r2
			eps_r1 = tri_sp_parameters(1)%eps_r1
			call set_local_material_parameters_tri()
			
		else if (pre_types%object .eq. 'surface') then
			if (allocated(tri_sf_parameters)) then
				deallocate(tri_sf_parameters)
			end if
			if (allocated(struc_tri_surfaces)) then
				deallocate(struc_tri_surfaces)
			end if
			allocate(tri_sf_parameters(number_objects))
			allocate(struc_tri_surfaces(number_objects))
			call set_surface_parameters_tri()
			
			eps_r2 = tri_sf_parameters(1)%eps_r2
			eps_r1 = tri_sf_parameters(1)%eps_r1
			call set_local_material_parameters_tri()
		else 
			print*, 'Not a proper object type.'
		end if
		
		if (data_structure_type .eq. 'top_down') then
			tree_l_max = 4	
			tree_bounding_box(1)%x = (/-1.0, -1.0, -1.0/)*(illumination_p%lambda)*2**(tree_l_max-1)
			tree_bounding_box(2)%x = (/1.0, 1.0, 1.0/)*illumination_p%lambda*2**(tree_l_max-1)!
						
		end if 
	end subroutine
	
	subroutine set_local_material_parameters_tri()		
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
	
	subroutine read_pp_tt_file_tri(struct)	
		implicit none
		type(structure_tri) :: struct		
		integer, dimension(:, :), allocatable :: Element_NodeIndices	
      real(dp), dimension(:, :), allocatable :: Node_coordinates				
		integer :: Nt, Np, j
		
		call read_file_array_size_returned_tri(File_NodeCoordinates, Np)
				
		allocate(Node_coordinates(3, Np))
		allocate(struct%points(Np))

		open(unit = 81, file = File_NodeCoordinates, status = 'old', action='read') !			
			read(81, *) Node_coordinates
		close(81) 
				
		do j = 1, Np 
			struct%points(j)%point = Node_coordinates(:, j) !
		end do
		
		call read_file_array_size_returned_tri(File_Element_NodeIndices, Nt)
		allocate(Element_NodeIndices(3, Nt))  
		allocate(struct%elements(Nt))
   
		!pp: array of cartesian node coordinates x, y and z
		open(unit = 80, file = File_Element_NodeIndices, status = 'old', action='read') !
			read(80, *) Element_NodeIndices
		close(80)
		
		do j = 1, Nt !tt, node indices
			allocate(struct%elements(j)%corners(3)) ! 
			struct%elements(j)%corners = Element_NodeIndices(:, j)
		end do
		deallocate(Element_NodeIndices)
		deallocate(Node_coordinates)
	end subroutine 
	
	subroutine read_file_array_size_returned_tri(file_name, Np)
		implicit none
		character(len = 100), intent(in) :: file_name
		
		real(dp) :: wert
		integer, intent(out) :: Np		
		integer :: status_open           ! Rueckgabewert aus iostat bei open
      integer :: status_read
		
		Np = 0
		open(unit=25, file = file_name, status='old', action='read', iostat=status_open)
			oeffnen_1: if ( status_open == 0 ) then
				einlese_schleife_1: do
				read (25, *, iostat=status_read) wert
					if ( status_read /= 0 ) exit
						Np = Np + 1
					end do einlese_schleife_1!
			readif_1: if ( status_read > 0 ) then
						write(*,*) 'Beim Lesen von Zeile ', Np + 1, &
							' ist ein Fehler aufgetreten'
						else ! status_read < 0
							write(*,*)
							write(*,*) ' Hinweis: das Dateiende wurde erreicht'
							write(*,*) ' => N = ', Np, 'nodes or elements'
						end if readif_1
					else oeffnen_1
						write(*,*) 'Beim OEffnen der Datei trat &
										&  Systemfehler Nr. ', status_open,' auf'
					end if oeffnen_1
		close( unit=25 )  !Datei schliessen
	end subroutine read_file_array_size_returned_tri
	
	!The mesh file from COMSOL has to be treated 
	!Especially, the index of the corner in tt file begins at zero,
	!which is not programmed in this simulator
	subroutine read_file_meshed_comsol()
		implicit none
		
		character(len = 100) :: file_name_out, file_name_in
		integer :: in_unit,  mm, out_unit 
		integer, allocatable, dimension(:, :) :: tt
		real(kind = 8), allocatable, dimension(:, :) :: pp
		
		intrinsic :: max
		integer :: Nt, Np
	
		in_unit = 20
		out_unit = 30	
		
		write(*,*) 'Give the file name of the surface'
		read(*,'(A)') file_name_surface
		
		file_name_in = trim(adjustl(file_name_surface))//'_pp.mphtxt'	
	
		file_name_out = 'pp_'//trim(adjustl(file_name_surface))//'.txt'
	
		call read_file_array_size_returned_tri(file_name_in, Np)
	
		print*, 'Np=', Np
	
		allocate(pp(3, Np))
	
		open(unit = in_unit, file = trim(file_name_in))	
		do mm = 1, Np
			read(in_unit, *) pp(1:3, mm)		
		end do
		close(in_unit)	
		
		open(unit = out_unit, file = trim(file_name_out))
		do mm = 1, Np
		 !  tmp = pp(1, mm)
			!pp(1, mm) = pp(3, mm)
			!pp(3, mm) = pp(2, mm)
			!pp(2, mm) = tmp		
			pp(3, mm) = pp(3, mm)	!reverse the direction of the surface
			write(out_unit, *) pp(:, mm)		
		end do	
		close(out_unit)	
	
	file_name_in = trim(adjustl(file_name_surface))//'_tt.mphtxt'
	
	file_name_out = 'tt_'//trim(adjustl(file_name_surface))//'.txt'
	
	call read_file_array_size_returned_tri(file_name_in, Nt)
	
	print*, 'Nt=', Nt	
	allocate(tt(3, Nt))	
	open(unit = in_unit, file = trim(file_name_in))	
	do mm = 1, Nt
		read(in_unit, *) tt(1:3, mm)		
	end do
	close(in_unit)	
	!	
	open(unit = in_unit, file = trim(file_name_out))	
	do mm = 1, Nt
		write(in_unit, *) tt(1:3, mm)+1		
	end do
	close(in_unit)	
	
	end subroutine
	
end module lib_sie_tri_data_container