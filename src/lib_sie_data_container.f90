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
	
module lib_sie_data_container
	use lib_sie_constants
	use lib_sie_type_function
	use ml_fmm_type	
	
	implicit none	
	
	!public variables
	real(dp) :: k0 
	real(dp) :: Omega	
	complex(dp) :: eps_1, eps_2
	complex(dp) :: eps_r1, eps_r2, eps_r2_main
	complex(dp) :: k1, k2
	real(dp) ::	my_1, my_2 
	complex(dp) ::	eta_1, eta_2
	integer(kind = 1) :: total_field
	
	integer :: calc_p(8)
	character(100) :: file_name_surface	
	character(len = 100) :: File_NodeCoordinates
	character(len = 100) :: File_Element_NodeIndices	
	character(100) :: file_name_output_I
	character(100) :: file_name_output_Efield
	character(100) :: file_name_output_Hfield	
	character (len=100) :: file_out_error
	character (len=100) :: file_out_parameters
	character (len=20) :: field_output
		
	real(dp) :: theta_start, theta_end, phi_start, phi_end
	real(dp) :: geometry_p(3)
	real(dp) :: surface_center(3)
	real(dp), dimension(:, :), allocatable :: surface_center_arr
	real(dp) :: position_c, dim_a_min, dim_a_max, dim_b_min, dim_b_max
	integer :: sampling_na, sampling_nb, n_disc
	
	type(evaluation_r_media), dimension(:), allocatable :: r_media
	!--------------------------------------------------------------------------	
	real(dp), public :: beam_waist 
	type(lib_sie_illumination_parameter), public :: illumination_p
	type(preset_type), public :: pre_types	
	type(lib_sie_parameter_objective), public :: p_obj
	
	type(lib_sie_evaluation_parameter_type), public :: evaluation_parameter
	real(dp), dimension(2), parameter :: nearfield_distance = (/200.0e-9, 500.0e-9/)
	integer, public :: number_objects 
	integer, parameter :: ng = 3
	integer, parameter :: ng4 = 4
	
	integer :: n_el, m_pairs, m_edge !m_edge for quadrilateral elements		
	
	!for x in GMRES
	real(dp), public :: x_initial
	integer(kind=1) :: m_tree_l_max
	type(box_size_at_l) :: box_at_l	
	integer(kind=1), public :: tree_l_max_predefined
	integer :: truncation_number_default
	integer(kind = 4) :: tree_s_opt	
	character(len = 20), public :: data_structure_type
	logical :: precondition_left
	integer :: max_iterations, restart
	real(dp) :: convergence_tolerance
	
	
	!Parameters for MLFMM	
	type(lib_sie_simulation_parameter_type) :: simulation_data
	type(truncation_hl), dimension(:), allocatable :: truncation_number	
	type(lib_ml_fmm_coefficient), dimension(:), allocatable :: preCalculated_RT
	type(lib_ml_fmm_coefficient), dimension(:, :), allocatable :: preCalculated_RT_quad
	type(lib_ml_fmm_coefficient) :: fxy	
	
	integer(kind=1) :: tree_l_max	
	double precision, dimension(:,:), allocatable :: k_hat_p
	double precision, dimension(:), allocatable :: TRF_p	
	double precision, dimension(:,:), allocatable :: transfer_theta
	double precision, dimension(:,:), allocatable :: transfer_phi
	complex(dp), dimension(:), allocatable :: diagonal_element
	type(lib_tree_spatial_point), dimension(2) :: tree_bounding_box 
	type(lib_sie_mlfmm_TL_coefficient) :: preCalculated_TL			
	complex(dp), public :: eta_a
	
end module