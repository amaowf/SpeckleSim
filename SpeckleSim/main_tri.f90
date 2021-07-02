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
	
program main_tri
	use lib_sie_solver_gmres_interface
	use lib_sie_tri_solver_normal	
	use lib_sie_quad_solver_normal	
	use lib_sie_data_container
	use lib_sie_tri_data_container
	
	implicit none
	integer :: file_n, file_start_n
	
	character(len=100) :: str1, str2	
	
	!The mesh data from COMSOL has to be modified slightly
	!Especially for the start number of indexing (
	!call read_file_meshed_comsol()
		
	call initalize_calculation()
	call start_sie_calculation()	
	
	contains
			
	subroutine set_calculation_parameters()
		implicit none
		
		!calc_p(1) : illumination type
		!1- plane wave
		!2- Gaussian beam
		!3- Focused Guassian beam
		calc_p(1) = 2
		
		!When Gaussian beam
		if (calc_p(1) .eq. 2) then
			beam_waist = 1.5e-6 !waist radius
		end if
		
		!calc_p(2) : calculation type
		!1- tri, triangular element, normal calculation (LU-decomposition) 
		!2- quad, quadrilateral element, normal calculation (LU-decomposition)
		!3- Gmres_tri, MLFMM using GMRES iterative process
		calc_p(2) = 1
		
		!calc_p(3) : formulation type
		!1- 'PMCHWT' 
		!2- 'MCTF'
		!3- 'ICTF'		
		calc_p(3) = 1
		
		!calc_p(4) : field evaluation plane/method
		!1- 'xz' 
		!2- 'xy'
		!3- 'yz'
		!4- 'rcs_n',radar cross section n-polarization 
		!5- 'rcs_p', radar cross section p-polarization 
		!6- 'BRDF', bidirectional refelction distribution function
		calc_p(4) = 5
		
		!calc_p(5) : object type
		!1- 'sphere'
		!2- 'surface'
		!3- 'incident' !only incident field is calculated		
		calc_p(5) = 2
		
		if (calc_p(5) .eq. 2)then
			surface_length = 4.0e-6	
		end if
		
		!calc_p(6) : field evaluation
		!1- surface current or field
		!2- observation field
		!3- both are calculated subsequently	
		calc_p(6) = 2
		
		!Relative permittivity of the surface (object)
		!environment is air
		eps_r2_main =  (-16.075, -0.442342)!(-17.2, -0.498)! (2.25, 0.0)!
		total_field = 1 ! '1' total field including incident field is calculted
		              ! '0' only scattered field is calculated

		!only for the case of rcs or BRDF
		if (calc_p(4) .gt. 3)then
			theta_start = -75.0/180.0*PI 
			theta_end = 75.0/180.0*PI 
			phi_start = 0.0 
			phi_end = 0.0
		end if
	end subroutine
	
	!Initalize GMRES calculations
	subroutine initalize_gmres()	
	
		x_initial = 0.0				
		tree_s_opt = 7
		truncation_number_default = 15
		max_iterations = 12000
		
		!when restart is not wished
		!set restart large
		restart = 12000
		convergence_tolerance = 1.0e-3
		precondition_left =  .true. !.false. !
		
		!There is a bug for top_down case
		data_structure_type =  'bottom_up'!'top_down'! 
			
		if (data_structure_type .eq. 'top_down') then
			tree_l_max_predefined = 4			
		end if
		
	end subroutine
		
	!Convention of exp(iwt) for the EM wave
	!Illuminaiton parameters
	subroutine set_illumination_parameters()
		implicit none
		
		real(dp) :: theta, phi		
		illumination_p%lambda = 600e-9		
		illumination_p%theta_in = 30.0*PI/180
		illumination_p%phi_in = 0.0*PI/180
		theta = illumination_p%theta_in
		phi = illumination_p%phi_in		
		illumination_p%k_in = &
			(/sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)/)	
		illumination_p%E_in = &
			(/cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)/)
		
	end subroutine
	
	subroutine initalize_calculation()		
		
		call set_calculation_parameters()
		call set_illumination_parameters()
		call set_calculation_types()
	
		if (calc_p(5) .eq. 1) then	
			pre_types%object = 'sphere'
			call name_outputfile_single()
			
		else if (calc_p(5) .eq. 2) then	
			pre_types%object = 'surface'	
			write(*,*) 'Give the number of files to be calculated'
			read(*,'(I3.3)') file_n			
			file_n = file_n -1 
			
			if (file_n .eq. 0)then
				file_start_n = 1
			else
				write(*,*) 'Give the start-number of files'
				read(*,'(I3.3)') file_start_n			
			end if
			
			print*, '                           '			
			print*, 'Calculation for a surface'
			print*, '-------------------------'
			!The input name should not include pp_ por tt_
			!If only one file to be calculated, the input name should include .txt
			!if multiple files to be calculated, 
			!the input name should not include '_001.txt'!
			write(*,*) 'Give the file name of the surface'
			read(*,'(A)') file_name_surface
		else
			print*, 'Not a proper object type'
		end if	
			
	end subroutine	
	
	!No setting is necessary
	subroutine start_sie_calculation()
		implicit none
		integer :: m
		
		if (calc_p(2) .eq. 1) then
			print*, '                                       '	
			print*, 'Normal calculation using triangular elements'
			print*, '--------------------------------------------'			
			if (file_n .gt. 1) then
				do m = file_start_n, file_n+file_start_n
					print*, 'The index of the file is', m
					call name_outputfile_multiple(m)
					call run_normal_tri_calculation()
				end do
			else 
				call name_outputfile_single()
				call run_normal_tri_calculation()
			end if 	
		else if (calc_p(2) .eq. 2) then		
			print*, '                                       '	
			print*, 'Normal calculation using quadrilateral elements'
			print*, '-------------------------------------------------'					
			if (file_n .gt. 1) then
				do m = file_start_n, file_n + file_start_n
					call name_outputfile_multiple(m)
					call run_normal_quad_calculation()
				end do
			else 
				call name_outputfile_single()
				call run_normal_quad_calculation()
			end if 
		else if (calc_p(2) .eq. 3) then			
			print*, '                                       '	
			print*, 'GMRES solver using triangular elements'
			print*, '------------------------------------------------'		
			call initalize_gmres()
			if (file_n .gt. 1) then
				do m = 1, file_n
					call name_outputfile_multiple(m)
					if (calc_p(6) .eq. 2) then
						call lib_sie_tri_solver_normal_II()
					else if (calc_p(6) .eq. 3) then
						call run_ml_fmm_GMRES_tri()
						call lib_sie_tri_solver_normal_II()
					else			
						call run_ml_fmm_GMRES_tri()
					end if
				end do
			else 
				call name_outputfile_single()
				call initalize_gmres()
				if (calc_p(6) .eq. 2) then
					call lib_sie_tri_solver_normal_II()
				else if (calc_p(6) .eq. 3) then
					call run_ml_fmm_GMRES_tri()
					call lib_sie_tri_solver_normal_II()
				else			
					call run_ml_fmm_GMRES_tri()
				end if			
			end if 
		end if		
	end subroutine
	
	!No setting is necessary
	subroutine set_calculation_types()
		!
		if (calc_p(1) .eq. 1) then	
			pre_types%illumination = 'plane'! 
			print*, '                              '			
			print*, 'Plane wave illumination'
			print*, '------------------------'
		else if (calc_p(1) .eq. 2) then	
			pre_types%illumination = 'Gaussian'!
			print*, '                              '
			print*, 'Gaussian beam illumination'
			print*, '------------------------'
		else if (calc_p(1) .eq. 3) then		
			pre_types%illumination = 'Focused' !  
			print*, '                              '			
			print*, 'Focused beam illumination'
			print*, '------------------------'			
		else
			print*, 'Not a proper illumination type'
		end if
		
		if (calc_p(2) .eq. 1) then	
			pre_types%calculation = 'Normal_tri'
		else if (calc_p(2) .eq. 2) then	
			pre_types%calculation = 'Normal_quad'
		else if (calc_p(2) .eq. 3) then		
			pre_types%calculation = 'Gmres_tri'			
		else
			print*, 'Not a proper calculation type'
		end if
		
		if (calc_p(3) .eq. 1) then	
			pre_types%formulation = 'PMCHWT'
		else if (calc_p(3) .eq. 2) then	
			pre_types%formulation = 'MCTF'
		else if (calc_p(3) .eq. 3) then		
			pre_types%formulation = 'ICTF'
		else if (calc_p(3) .eq. 4) then	
		!Not exact, maybe there is a bug
			pre_types%formulation = 'CTF'
		else
			print*, 'Not a proper formulation type'
		end if

		if (calc_p(4) .eq. 1) then	
			pre_types%evaluation = 'xz'
		else if (calc_p(4) .eq. 2) then	
			pre_types%evaluation = 'xy'
		else if (calc_p(4) .eq. 3) then		
			pre_types%evaluation = 'yz'
		else if (calc_p(4) .eq. 4) then	
			pre_types%evaluation = 'rcs_n'
		else if (calc_p(4) .eq. 5) then		
			pre_types%evaluation = 'rcs_p'
		else if (calc_p(4) .eq. 6) then		
			pre_types%evaluation = 'BRDF'			
		else
			print*, 'Not a proper evaluation type'
		end if
		
	end subroutine
	
	!Only multiple surfaces are calculated in one run
	subroutine name_outputfile_multiple(nn)
		integer, intent(in) :: nn

		character(len=8) :: file_nr 
		write(file_nr, '(I3.3)') nn			
		if (calc_p(2) .eq. 2)then !quad
			str1 = trim(adjustl(pre_types%formulation))//'_'//trim(file_nr)//'.txt'	
			str2 = trim(adjustl(pre_types%evaluation))//'_'//trim(file_nr)//'.txt'
			file_name_surface = trim(adjustl(file_name_surface))//'_'//trim(file_nr)//'.txt'
			file_name_output_I = 'Result_SEH_quad_'//trim(str1)	
			file_name_output_II = 'Result_Efield_quad_'//trim(str2)
		else !GMRES and tri
			str1 = trim(adjustl(pre_types%formulation))//'_'//trim(file_nr)//'.txt'	
			str2 = trim(adjustl(pre_types%evaluation))//'_'//trim(file_nr)//'.txt'				
			File_NodeCoordinates = 'pp_'//trim(adjustl(file_name_surface))//'_'//trim(file_nr)//'.txt'	
			File_Element_NodeIndices = 'tt_'//trim(adjustl(file_name_surface))//'_'//trim(file_nr)//'.txt'			
			file_name_output_I = 'Result_SEH_tri_'//trim(str1)
			file_name_output_II = 'Result_Efield_tri_'//trim(str2)					
			if (calc_p(5) .eq. 3)then
				str1 = trim(adjustl(pre_types%formulation))//'_'//trim(file_nr)
				file_out_error = 'iteration_error_'//trim(str1)								
				file_out_parameters = 'Calculation_parameters_GMRES_'//trim(str1)
			else !calc_p(5) = 1
				file_out_parameters = 'Calculation_parameters_normal_'//trim(str1)
			end if	
		end if	!	
	end subroutine
	
	subroutine name_outputfile_single()	
		character(len = 10) :: str3
		str1 = trim(adjustl(pre_types%formulation))//'.txt'	
		str2 = trim(adjustl(pre_types%evaluation))//'.txt'	
		
		if (calc_p(2) .eq. 2)then !quad
			str3 = 'quad'
			File_NodeCoordinates = 'pp_'//trim(adjustl(file_name_surface))//'.txt'	
			File_Element_NodeIndices = 'tt_'//trim(adjustl(file_name_surface))//'.txt'	
			file_out_parameters = 'Calculation_parameters_normal_'//trim(str3)//'.txt'
		else
			str3 = 'tri'!//trim(adjustl(pre_types%formulation))
			File_NodeCoordinates = 'pp_'//trim(adjustl(file_name_surface))
			File_Element_NodeIndices = 'tt_'//trim(adjustl(file_name_surface))
			if (calc_p(2) .eq. 1)then !tri normal
				file_out_parameters = 'Calculation_parameters_normal_'//trim(str3)//'.txt'
			else !(calc_p(2) .eq. 3), GMRES 												
				file_out_parameters = 'Calculation_parameters_GMRES_'//trim(str3)//'.txt'
				file_out_error = 'iteration_error_'//trim(str3)
			end if!			
		end if		
		file_name_output_I = 'Result_SEH_'//trim(str3)//'_'//trim(str1)
		file_name_output_II = 'Result_Efield_'//trim(str3)//'_'//trim(str2)
		print*, 'file_name_output_II =', file_name_output_II
	end subroutine
	
end program