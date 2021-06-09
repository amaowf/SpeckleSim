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
	
module lib_sie_quad_calculation_mod

	use lib_sie_math
	use lib_sie_constants
	use lib_sie_quad_data_container
	use lib_sie_integration_points
	
	implicit none
	
	private
	
	public :: lib_sie_quad_input_field_calculation
	public :: ngp_near_field_distance_quad
	public :: phi_range_4Feld
	public :: r_shape	
	public :: f_edge_new
	public :: n_differential
	public :: rho_2_4feld
	public :: data_import_quad
	public :: paras_st_a
	public :: paras_st_b
	public :: diver_green
	public :: neighbours_edgecorner_full
	public :: find_edge_index	
	public :: incident_int
	public :: incident_gauss_quad
	public :: integration_scattered_field_quad	
	public :: sorting_edge	
	public :: neighbours_edgecorner
	public :: derived_element_parameters_TypeB
	public :: derived_element_parameters_TypeA
	public :: r_local_parametric
	public :: get_vector_b_quad	
	public :: find_midpoint_edge_quad
	
	!type(point), dimension(:), allocatable :: r_local
	
	type derived_element_parameters_TypeB
		real(dp), dimension(10) :: n_rot_N
		real(dp) :: JJ
		type(point) :: nq		
		type(point), dimension(10) :: Nqi
	end type
	
	type derived_element_parameters_TypeA
		real(dp) :: JJ
		type(point), dimension(10) :: Naj
	end type
	
   contains
	
	function lib_sie_quad_input_field_calculation(r_local)result(field_input)	

		implicit none
		type(point), dimension(:), allocatable, intent(in) :: r_local
		integer :: n1, n2, m, ngp
		
		complex(dp) :: k1
		
		type(lib_sie_evaluation_point_type), dimension(:), allocatable :: field_input		
		
		n1 = evaluation_parameter%N_dim(1)
		n2 = evaluation_parameter%N_dim(2)
		
		allocate(field_input(n1*n2))
		k1 = k0		
		
		do m = 1, n1*n2
			if (pre_types%illumination == 'Gaussian')then
				call Gaussian_beam(beam_waist, r_local(m)%point, illumination_p, field_input(m)%e_field, k1)
			else if (pre_types%illumination == 'Focused')then
				call Gaussian_beam_focused(beam_waist, r_local(m)%point, field_input(m)%e_field, k1)
			else if ((pre_types%illumination == 'Plane') .or. (pre_types%illumination == 'plane'))then
				field_input(m)%e_field%vector = illumination_p%E_in*exp(-im*dot_product(k1*illumination_p%k_in, r_local(m)%point))				
			else 
				print*, 'Not a proper illumination method'
				call exit
			end if
			ngp = ngp_near_field_distance_quad(r_local(m), nearfield_distance)
			Field_input(m)%h_field%vector = cross_c(k1*illumination_p%k_in, Field_input(m)%e_field%vector)/(my_1*c0)
			field_input(m)%coordinate = r_local(m)
			field_input(m)%ngp = ngp
		end do
		return
	end function
		
	subroutine get_vector_b_quad(Edge_I, vector_b)
		complex(kind = dp), dimension(:), allocatable, intent(out) ::  vector_b		
		type(edges_parameter), dimension(:), intent(in) :: edge_I
		!integer, dimension(:, :), intent(in) :: Edge_I		
		integer :: pI_a, ii, s, n
		real(dp), dimension(8, 3) :: Element_a		
		complex(kind = dp) :: k_hat(3)
		complex(dp) :: eta_a!
		complex(kind = dp), dimension(:, :), allocatable :: Sum_REH
		integer :: ngp 		
		
		ngp = 6
		allocate(vector_b(2*m_edge), Sum_REH(m_edge, 2))
		do s= 1, m_edge			
			ii = Edge_I(s)%edge_p(1)! element index of the s-th edge !
			pI_a = Edge_I(s)%edge_p(2) ! edge index in element-ii
			do n = 1, 8 
				Element_a(n, 1:3) = struc_quad%node_vector(struc_quad%elements(ii)%Vertices(n))%point(1:3)
			end do	
			k_hat = illumination_p%k_in*(1 + im*0)
			if (pre_types%illumination == 'Gaussian') then
				call incident_Gauss_quad(Element_a, pI_a, ngp, Sum_REH(s, :))
			else if ((pre_types%illumination == 'Plane') .or. (pre_types%illumination == 'plane')) then
				call incident_int(Element_a, pI_a, ngp, Sum_REH(s, :))
			else 
				print*, 'Not a proper illumination type'
			end if
		end do
		
		select case (pre_types%formulation)
			case ('PMCHWT')
				vector_b(1 : m_edge) = sum_REH(:, 1)
				vector_b(m_edge + 1 : 2*m_edge) = sum_REH(:, 2)
			case ('MCTF')
				vector_b(1 : m_edge) = sum_REH(:, 1)
				vector_b(m_edge + 1 : 2*m_edge) = sum_REH(:, 2)*eta_1*eta_2				
			case ('ICTF')
				print*, 'ICIF'
				eta_a = (eta_1 + eta_2)/2
				vector_b(1 : m_edge) = sum_REH(:, 1)/eta_a
				vector_b(m_edge + 1 : 2*m_edge) = sum_REH(:, 2)*eta_a
		end select
		deallocate(sum_REH)
		return
	end subroutine get_vector_b_quad
	
	function ngp_near_field_distance_quad(r_loc, distance_nf)result(ngp)
		integer :: ngp
		integer :: p, ngp_arr(number_objects)
		type(point) :: r_loc
		real(dp), dimension(2), intent(in) :: distance_nf
		real(dp) :: tmp_r		
		
		if (pre_types%object .eq. 'sphere')then
			do p = 1, number_objects
				tmp_r = (vec_len(r_loc%point - quad_sp_parameters(p)%centroid%point) - quad_sp_parameters(p)%D*0.5)			
				if ((tmp_r .le. distance_nf(1)) .and. (tmp_r .gt. 0.0) ) then
					ngp_arr(p) = 10
				else if ((tmp_r .le. distance_nf(2)) .and. (tmp_r .gt. distance_nf(1))) then
					ngp_arr(p) = 8
				else if (tmp_r .gt. distance_nf(2)) then
					ngp_arr(p) = 5
				else 
					ngp_arr(p) = 50
				end if
				ngp = maxval(ngp_arr)
			end do
		else if (pre_types%object .eq. 'surface')then
			ngp = 5 ! firstly consider only a rough surface
			!do p = 1, number_objects
			!	tmp_r = (vec_len(r_loc%point - quad_msf_parameters(p)%centroid%point) - quad_msf_parameters(p)%D*0.5)			
			!	if ((tmp_r .le. distance_nf(1)) .and. (tmp_r .gt. 0.0) ) then
			!		ngp_arr(p) = 10
			!	else if ((tmp_r .le. distance_nf(2)) .and. (tmp_r .gt. distance_nf(1))) then
			!		ngp_arr(p) = 8
			!	else if (tmp_r .gt. distance_nf(2)) then
			!		ngp_arr(p) = 5
			!	else 
			!		ngp_arr(p) = 0
			!	end if
			!	ngp = maxval(ngp_arr)
			!end do
		else
			print*, 'not a properly defined object type'
			ngp = 0
		end if
		return
	end function	
	
	subroutine find_midpoint_element_quad(struct)
		implicit none
		integer :: ne, t, n
		type(structure_quad), intent(inout) :: struct    
		real(kind = dp) :: mid(3), ksi, eta
		
		real(kind = dp), dimension(8, 3) :: Element_q
		ne=size(struct%elements)
		if (allocated(struct%midpoint))then
		deallocate(struct%midpoint)
		endif
		ksi = 0.0
		eta = 0.0
		allocate(struct%midpoint(ne))
		do t = 1, ne
			do n = 1, 8
				Element_q(n, 1:3) = struct%node_vector(struct%elements(t)%Vertices(n))%point(1:3)
			end do
			mid = r_local_parametric(ksi, eta, Element_q)
			struct%midpoint(t)%point=mid
		end do
   end subroutine
	
	subroutine find_midpoint_edge_quad(struct, m_edge)
		implicit none
		
		integer, intent(in) :: m_edge		
		type(structure_quad), intent(inout) :: struct    
		real(kind = dp), dimension(3) :: point_I, point_II
		real(kind = dp), dimension(8, 3) :: Element_q
		
		integer :: edge(10, 2), pI_a, node_I, node_II, ii, t, n
		real(kind = dp) :: mid(3), ksi, eta
		
		edge(1, 1:2) = (/1, 2/) ! 
      edge(2, 1:2) = (/2, 3/) ! 
      edge(3, 1:2) = (/7, 6/) ! 
      edge(4, 1:2) = (/6, 5/)
      edge(5, 1:2) = (/8, 4/)
      edge(6, 1:2) = (/1, 8/)
      edge(7, 1:2) = (/8, 7/)
      edge(8, 1:2) = (/3, 4/)
      edge(9, 1:2) = (/4, 5/)
      edge(10, 1:2) = (/2, 6/)
		
		if (allocated(struct%midpoint))then
			deallocate(struct%midpoint)
		endif
		
		allocate(struct%midpoint(m_edge))
		do t = 1, m_edge
			ii = edge_I(t)%edge_p(1)
			pI_a = edge_I(t)%edge_p(2)	
			
			node_I = edge(pI_a, 1)
			node_II = edge(pI_a, 2)
			
			point_I = struct%node_vector(struct%elements(ii)%Vertices(node_I))%point(1:3)
			point_II = struct%node_vector(struct%elements(ii)%Vertices(node_II))%point(1:3)
			struct%midpoint(t)%point = (point_I + point_II)/2
			
			!do n = 1, 8
			!	Element_q(n, 1:3) = struct%node_vector(struct%elements(ii)%Vertices(n))%point(1:3)
			!end do
			!mid = r_local_parametric(ksi, eta, Element_q)
			!struct%midpoint(t)%point=mid
			
		end do
   end subroutine find_midpoint_edge_quad
	
	
	!Finding the global coordinates at (ksi, eta)	
	!Eq. (3.7) and (3.8)
   function R_shape(ksi, eta)	
		integer, parameter :: NP = 8 ! Points number      		
		real(kind = dp), dimension(8) :: R_shape !The coordinate of one piont at ksi and eta
		real(kind = dp), intent(in) :: ksi, eta
		R_shape(1) = 0.25 * (1 - ksi) * (1 - eta) * (-1 - eta - ksi)
		R_shape(2) = 0.5 * (1 - ksi**2) * (1 - eta)
		R_shape(3) = 0.25 * (1 + ksi) * (1 - eta) * (-1 - eta + ksi)
		R_shape(4) = 0.5 * (1 - eta**2) * (1 + ksi)
		R_shape(5) = 0.25 * (1 + ksi) * (1 + eta) * (-1 + eta + ksi)
		R_shape(6) = 0.5 * (1 - ksi**2) * (1 + eta)
		R_shape(7) = 0.25 * (1 - ksi) * (1 + eta) * (-1 + eta - ksi)
		R_shape(8) = 0.5 * (1 - eta**2) * (1 - ksi) 
		return
	end function R_shape 
	
	!Differential of the shape function	R_shape
	function N_differential(ksi, eta)	
		real(kind = dp), dimension(8, 2) :: N_differential
		real(kind = dp), intent(in) :: ksi, eta
		real(kind = dp), dimension(8) :: Nd_ksi, Nd_eta!, Ndd_ksi, Ndd_eta
		Nd_ksi(1) =  0.25*(1 - eta)*(2*ksi+eta)
		Nd_ksi(2) = -ksi*(1 - eta)
		Nd_ksi(3) = 0.25*(1 - eta)*(2*ksi-eta)
		Nd_ksi(4) = 0.50*(1 - eta**2) 
		Nd_ksi(5) = 0.25*(1 + eta)*(2*ksi+eta)
		Nd_ksi(6) = - ksi*(1 + eta)
		Nd_ksi(7) = 0.25*(1 + eta)*(2*ksi-eta)
		Nd_ksi(8) = -0.50*(1- eta**2) 
   
		Nd_eta(1) = 0.25*(1 - ksi)*(2*eta + ksi)
		Nd_eta(2) = -0.5*(1 - ksi**2) 
		Nd_eta(3) = 0.25*(1 + ksi)*(2*eta - ksi)
		Nd_eta(4) = -eta*(1 + ksi)
		Nd_eta(5) = 0.25*(1 + ksi)*(2*eta + ksi)
		Nd_eta(6) = 0.5*(1-ksi**2) 
		Nd_eta(7) = 0.25*(1-ksi)*(2*eta - ksi)
		Nd_eta(8) = -eta*(1 - ksi)     
		N_differential(:, 1) = Nd_ksi(:)
		N_differential(:, 2) = Nd_eta(:)
		return
	end function N_differential
	
	!The second order differential of the shape function R_shape
	function Nd_differential(ksi, eta) !	
		real(kind = dp), dimension(8, 4) :: Nd_differential
		real(kind = dp), intent(in) :: ksi, eta
		real(kind = dp), dimension(8) :: Nd_ksi_ksi, Nd_eta_ksi, Nd_ksi_eta, Nd_eta_eta

		Nd_ksi_ksi(1) = 0.5*(1 - eta)
		Nd_ksi_ksi(2) = -(1 - eta)
		Nd_ksi_ksi(3) = 0.5*(1 - eta)
		Nd_ksi_ksi(4) = 0.0
		Nd_ksi_ksi(5) = 0.5*(1 + eta)
		Nd_ksi_ksi(6) = -(1 + eta)
		Nd_ksi_ksi(7) = 0.5*(1 + eta)
		Nd_ksi_ksi(8) = 0.0
	
		Nd_ksi_eta(1) = 0.25*(-2*ksi -2*eta + 1)
		Nd_ksi_eta(2) = ksi
		Nd_ksi_eta(3) = 0.25*(2*eta - 2*ksi - 1)
		Nd_ksi_eta(4) = -eta
		Nd_ksi_eta(5) = 0.25*(2*ksi + 2*eta + 1)
		Nd_ksi_eta(6) = -ksi
		Nd_ksi_eta(7) = 0.25*(2*ksi - 2*eta - 1)
		Nd_ksi_eta(8) = eta 
     
		Nd_eta_ksi(1) = 0.25*(-2*ksi - 2*eta + 1)
		Nd_eta_ksi(2) = ksi
		Nd_eta_ksi(3) = 0.25*(2*eta - 2*ksi - 1)
		Nd_eta_ksi(4) = -eta
		Nd_eta_ksi(5) = 0.25*(2*ksi + 2*eta + 1)
		Nd_eta_ksi(6) = -ksi
		Nd_eta_ksi(7) = -0.25*(2*eta - 2*ksi + 1)
		Nd_eta_ksi(8) = eta

		Nd_eta_eta(1) = 0.5*(1 - ksi)
		Nd_eta_eta(2) = 0.0
		Nd_eta_eta(3) = 0.5*(1 + ksi)
		Nd_eta_eta(4) = -(1 + ksi)
		Nd_eta_eta(5) = 0.5*(1 + ksi)
		Nd_eta_eta(6) = 0.0
		Nd_eta_eta(7) = 0.5*(1- ksi)
		Nd_eta_eta(8) = -(1 - ksi) 

		Nd_differential(:, 1) = Nd_ksi_ksi(:)
		Nd_differential(:, 2) = Nd_ksi_eta(:)
		Nd_differential(:, 3) = Nd_eta_ksi(:)
		Nd_differential(:, 4) = Nd_eta_eta(:)
	end function Nd_differential
	
	!Coefficients of the 10 edge-functions and its diffential
	!Eqs. (3.51a-3.51j) in Huber's disseration
   function f_edge_new(ksi,eta)    
      real(kind = dp), dimension(10, 3) :: f_edge_new
      real(kind = dp), dimension(10) :: df_ksi, df_eta
      real(kind = dp), intent(in) :: ksi, eta
		
      f_edge_new(1, 1) = -0.25*(1-eta)*(2*ksi+eta)
      f_edge_new(2, 1) = 0.25*(1-eta)*(2*ksi-eta)
      f_edge_new(3, 1) = -0.25*(1 + eta)*(2*ksi-eta)
      f_edge_new(4, 1) = 0.25*(1 + eta)*(2*ksi + eta)
      f_edge_new(5, 1) = 0.5*(1 - eta**2)
      f_edge_new(6, 1) = -0.25*(1 - ksi)*(2*eta + ksi)
      f_edge_new(7, 1) = 0.25*(1 - ksi)*(2*eta - ksi)
      f_edge_new(8, 1) = -0.25*(1 + ksi)*(2*eta - ksi)
      f_edge_new(9, 1) = 0.25*(1 + ksi)*(2*eta + ksi)
      f_edge_new(10, 1) = 0.5*(1 - ksi**2)
		
     
      df_ksi(1) =  -0.5*(1 - eta)!f_dksi
      df_ksi(2) =  0.5*(1 - eta)!f_dksi
      df_ksi(3) =  -0.5*(1 + eta)!f_dksi
      df_ksi(4) =  0.5*(1 + eta)!f_dksi
      df_ksi(5) = 0
      
      df_ksi(6) = 0.25*(2*eta + 2*ksi - 1) !f_dksi
      df_ksi(7) = 0.25*(2*ksi - 2*eta - 1) !f_dksi
      df_ksi(8) = 0.25*(2*ksi - 2*eta + 1) !f_dksi
      df_ksi(9) = 0.25*(2*ksi + 2*eta + 1) !f_dksi
      df_ksi(10) = -ksi !d
      
      df_eta(1) = 0.25*(2*ksi + 2*eta - 1)
      df_eta(2) = 0.25*(2*eta - 2*ksi - 1)
      df_eta(3) = 0.25*(2*eta - 2*ksi + 1)
      df_eta(4) = 0.25*(2*ksi + 2*eta + 1)
      df_eta(5) = -eta !df_deta
      
      df_eta(6) = -0.5*(1 - ksi)
      df_eta(7) = 0.5*(1 - ksi)
      df_eta(8) = -0.25*(1 + ksi)
      df_eta(9) = 0.25*(1 + ksi)
      df_eta(10) = 0 !df_deta
      
      f_edge_new(:, 2) = df_ksi(:)
      f_edge_new(:, 3) = df_eta(:)
      return    
   end function f_edge_new
	
	!Namely, Ni = fi*vi, i= 1, 10, fi = f_edge_new 
	!Ev: dr/deta, Eu: dr/dksi, !EEu: Eq.(A.6-A.7)		
   function Edge_vector_new(ksi, eta, Element_q) 
		implicit none
      real(kind = dp) :: EEu(3), EEv(3) 
      real(kind = dp) :: Normal_n(3)
      real(kind = dp) :: f_temp(10,3), V_temp(3), V2(3), V1(3)  
      real(kind = dp) :: Edge_vector_new(10, 3) 
      real(kind = dp), dimension (8, 2) :: Nod_d  
      integer :: i, j
      real(kind = dp), intent(in) :: Element_q(8,3)
      real(kind = dp), intent(in) ::  ksi, eta 
        
		!Differential of the shape function
      Nod_d = N_differential(ksi, eta) 		
      do j = 1, 3
         EEu(j) = dot_product(Nod_d(1:8, 1),  Element_q(1:8, j))
         EEv(j) = dot_product(Nod_d(1:8, 2),  Element_q(1:8, j)) 
      end do
      Normal_n = cross_r(EEu, EEv)/vec_len(cross_r(EEu, EEv))
        
      f_temp = f_edge_new(ksi,eta) 
		!Eq. 3.49a in Huber's diss
      do i = 1, 5  
         V_temp = cross_r(EEv, Normal_n)
         V1 = V_temp/dot_product(EEu, V_temp)      
         Edge_vector_new(i, :) = f_temp(i,1)*V1(:)
      end do
		!Eq. 3.49b
      do i = 6, 10  
         V_temp = cross_r(Normal_n, EEu)
         V2 =  V_temp/dot_product(EEu, cross_r(EEv, Normal_n))
         Edge_vector_new(i, :) = f_temp(i,1)*V2(:)   
      end do
      return 
   end function Edge_vector_new  !

	!Eq. (3.57) in Huber
	function n_rot_N(ksi, eta, JJ) 
		real(kind = dp), intent(in) :: ksi, eta, JJ
		real(kind = dp), dimension(10) :: n_rot_N
		integer :: i
		real(kind = dp), dimension(10,3):: f_temp  
		f_temp = f_edge_new(ksi,eta)
		!!JJ = EEu x EEv
		do i = 1, 5     
			n_rot_N(i) = -f_temp(i,3)/JJ !df/deta  
		end do  
      
		do i = 6, 10     
			n_rot_N(i) = f_temp(i,2)/JJ !df/dksi 
		end do  
	end function n_rot_N
	 
   function Orts_coordinate(rn, ra, boundary, error)
   ! Find out the original ort vector when ksi and eta in the parametric system are given. 
        implicit none
        real(kind = dp), intent(in) :: rn(8, 3), ra(3)          
        real(kind = dp) :: dksi, ksi, eta, deta, boundary, error, r(3), Orts_coordinate(3), ksi_a, eta_a
        integer :: i, j, n     
        dksi = 2*boundary/100 
        deta = dksi
        Orts_coordinate = (/0.0, 0.0, 0.0/)
        
outer: do i = 1, 101
            ksi = 0-boundary+dksi*(i-1)
            ksi_a= ksi
            do j = 1, 101
                eta = -boundary+deta*(j-1)
                eta_a = eta
            do n = 1, 3
                r(n) = dot_product(R_shape(ksi_a, eta_a),  rn(1:8, n)) 
            end do
            if ((abs(r(1)-ra(1)) < error) .and. (abs(r(2)-ra(2)) < error)) then
            exit outer 
            end if
            end do          
        end do outer
        Orts_coordinate(1) = ksi 
        Orts_coordinate(2) = eta 
        return
   end function Orts_coordinate   
	
	!R_p(ksi, eta): local vector in the parametric coordinate system
	!Eq.(3.8) in Huber
	!Normal vector is a function of (ksi, eta)
	subroutine paras_st_a(mu_s, nu_s, Element_s, delement_pa)
		implicit none
		real(kind = dp), intent(in) :: mu_s, nu_s, Element_s(8, 3)
		integer :: n
		real(kind = dp) :: f_diff(8, 2), v_temp(3), f_temp(10, 3), Eu(3), Ev(3), Normal_n(3)
		type(derived_element_parameters_TypeA), intent(out) :: delement_pa
		
		f_diff = N_differential(mu_s, nu_s)
		do n = 1, 3
			Eu(n) = dot_product(f_diff(1:8, 1),  Element_s(1:8, n)) !dr/dksi
			Ev(n) = dot_product(f_diff(1:8, 2),  Element_s(1:8, n)) !dr/deta
		end do        
		delement_pa%JJ = vec_len(cross_r(Eu, Ev)) 
		Normal_n = cross_r(Eu, Ev)/delement_pa%JJ
		f_temp = f_edge_new(mu_s, nu_s)
		do n = 1, 10                
			if (n<=5) then
				V_temp = cross_r(Ev, Normal_n)/delement_pa%JJ
				delement_pa%Naj(n)%point = f_temp(n, 1)*V_temp
			else                
				V_temp = cross_r(Normal_n, Eu)/delement_pa%JJ
				delement_pa%Naj(n)%point= f_temp(n, 1)*V_temp
			end if
		end do
		return            
	end subroutine paras_st_a  
	
	function R_local_parametric(mu_s, nu_s, Element_s)result(rv)
		implicit none
		real(kind = dp), intent(in) :: mu_s, nu_s, Element_s(8, 3)
		real(kind = dp) :: ff(8), rv(3)
		integer :: n
		
		ff = R_shape(mu_s, nu_s)
		do n = 1, 3
			rv(n) = dot_product(ff(1:8),  Element_s(1:8, n))
		end do
   end function R_local_parametric
	 
	subroutine paras_st_b(mu_s, nu_s, Element_s, delement_pq)   
      implicit none
		type(derived_element_parameters_TypeB), intent(out) :: delement_pq
      real(kind = dp), intent(in) :: mu_s, nu_s, Element_s(8, 3)
      integer :: n
      real(kind = dp) :: f_diff(8, 2), v_temp(3), f_temp(10, 3), r_tmp(3), Eu(3), Ev(3)
      
      f_diff = N_differential(mu_s, nu_s)       
      do n = 1, 3            
         Eu(n) = dot_product(f_diff(1:8, 1),  Element_s(1:8, n)) !dr/dksi
         Ev(n) = dot_product(f_diff(1:8, 2),  Element_s(1:8, n)) !dr/deta
      end do
		r_tmp = cross_r(Eu, Ev)
      delement_pq%JJ = vec_len(r_tmp)
      delement_pq%nq%point = r_tmp/delement_pq%JJ
      f_temp = f_edge_new(mu_s, nu_s)
      do n = 1, 10                
			if (n<=5) then
				V_temp = cross_r(Ev, delement_pq%nq%point)/delement_pq%JJ
				delement_pq%Nqi(n)%point = f_temp(n, 1)*V_temp
            delement_pq%n_rot_N(n) = -f_temp(n, 3)/delement_pq%JJ!*vec_len(Ev)
         else                
				V_temp = cross_r(delement_pq%nq%point, Eu)/delement_pq%JJ					
            delement_pq%Nqi(n)%point = f_temp(n, 1)*V_temp
				delement_pq%n_rot_N(n) = f_temp(n, 2)/delement_pq%JJ!*vec_len(Eu)
         end if
      end do
      return            
   end subroutine paras_st_b
  
	subroutine Diver_Green(km, R_norm, Greenf, GG_SS)
		implicit none
		complex(kind = dp), intent(in) :: km
		real(kind = dp), intent(in) :: R_norm
		complex(kind = dp), intent(out) :: GG_SS, Greenf
		complex(kind = dp) :: Exp_f

		Exp_f = exp(-im*km*R_norm) ! It should be vec_len(R). m: medium
		GG_SS = -(1 + im*km*R_norm)*Exp_f/(4*PI*R_norm**2)!r' in Mitschang's paper is here Rq 
		Greenf = Exp_f/(4*PI*R_norm) !
		return 
	end subroutine Diver_Green
	
	!For the two edge adjecent elements. 
	!Only the first two nodes of the edge are necessary
	subroutine find_edge_index(S_edge, S_edge_index)
      implicit none        
      
      integer, dimension(6), intent(in) :: S_edge
      integer, dimension(2), intent(out) :: S_edge_index
      integer :: edge(10, 2), m, n, j, kn, pI_a, pII_a, A_edge(6)    
        
      edge(1, 1:2) = (/1, 2/) ! 
      edge(2, 1:2) = (/2, 3/) ! 
      edge(3, 1:2) = (/7, 6/) ! 
      edge(4, 1:2) = (/6, 5/)
      edge(5, 1:2) = (/8, 4/)
      edge(6, 1:2) = (/1, 8/)
      edge(7, 1:2) = (/8, 7/)
      edge(8, 1:2) = (/3, 4/)
      edge(9, 1:2) = (/4, 5/)
      edge(10, 1:2) = (/2, 6/)
		
		A_edge = S_edge
      if ((A_edge(2) - A_edge(1)) .NE. 1 ) then !the third and second node on the edge
         m = A_edge(2) 
         A_edge(2) = A_edge(3) 
         A_edge(3) = m  
         m = A_edge(5) 
         A_edge(5) = A_edge(6) 
         A_edge(6) = m 
      end if 
		do j = 1, 10
			if ((A_edge(1) .eq. Edge(j, 1)) .and. (A_edge(2) .eq. Edge(j, 2))) then
				pI_a = j !Edge index of element-s
				do kn = 1, 10
					if ((A_edge(4) .eq. Edge(kn, 1)) .and. (A_edge(5) .eq. Edge(kn, 2))) then
						pII_a = kn !
						exit
					else if ((A_edge(5) == Edge(kn, 1)) .and. (A_edge(4) == Edge(kn, 2))) then
						pII_a = kn
						exit
					end if
				end do !Loop kn  				
			else if ((A_edge(2) == Edge(j, 1)) .and. (A_edge(1) == Edge(j, 2))) then
				pI_a = j 
				do kn = 1, 10
					if ((A_edge(4) .eq.  Edge(kn, 1)) .and. (A_edge(5) .eq.  Edge(kn, 2))) then
						pII_a = kn						
						exit						
					else if ((A_edge(5) .eq. Edge(kn, 1)) .and. (A_edge(4) .eq. Edge(kn, 2))) then
						pII_a = kn
						exit
					end if
				end do !Loop kn					
			end if  			
      end do
		s_edge_index = (/pI_a, pII_a/)
      return
	end subroutine find_edge_index
	
	subroutine Sorting_edge(struc_temp, N_edge, edge_vs_element, element_vs_edge)!Edge_I, Edge_II,
		implicit none
		type(structure_quad), intent(in) :: struc_TEMP		  
        
		integer :: A_edge(7), edge_p(4)!
		integer, intent(out) :: N_edge 
		integer, dimension (N_el, 10) :: Edge_Order
		integer :: i, j, s, p, q, n_CommonEdge, m		  
		type(element_edges), dimension(:), allocatable, intent(out) :: element_vs_edge
		  type(edges_parameter), dimension(:), allocatable, intent(out) :: edge_vs_element
		  type(edges_parameter), dimension(:), allocatable :: nr_edges_tmp, edge_vs_element_tmp, nr_element_tmp
		  
		  allocate(element_vs_edge(N_el), edge_vs_element_tmp(n_el*10), nr_edges_tmp(10), nr_element_tmp(4))
		
		!allocate(element_vs_edge(N_el), edge_vs_element_tmp(n_el*10), nr_edges_tmp(10))
		  
		do i = 1, 10
		nr_edges_tmp(i)%edge_p(1:5) = 0
		end do
		  
		do i = 1, N_el*10             
			edge_vs_element_tmp(i)%edge_p(1:5) = (/0, 0, 0, 0, 0/) 
		end do   
		  
		do s = 1, N_el !Initialize the edge matrix 
			Edge_Order(s, 1:10) = 1 !
		end do  
      
      !Sort the independent edges
      p = 0
		do s = 1, N_el!
			n_commonedge = size(struc_quad%elements(s)%edge_overlap)
			q = 0
			do j = 1, 10 !Edges in element-s
				if (Edge_order(s, j) .NE. 0) then
					p = p + 1
					q = q + 1
					edge_vs_element_tmp(p)%edge_p = (/s, j, 0, 0, 0/)					
					do i = 1, n_CommonEdge ! 		
						A_edge(1:7) = struc_temp%elements(s)%Edge_overlap(i)%CommonEdge_Elements(2:8)
						if ((A_edge(3) - A_edge(2)) .NE. 1 ) then !
							m = A_edge(3) 
							A_edge(3) = A_edge(4)  
							A_edge(4) = m
							m = A_edge(6) 
							A_edge(6) = A_edge(7) 
							A_edge(7) = m 
						end if
						edge_p = (/j, 0, 0, 0/)					
						call find_edges_in_element(a_edge, edge_p)						
						if (edge_p(2) .NE. 0) then   !                 
							Edge_order(edge_p(2), edge_p(3)) = 0 ! The edge has been used	
							nr_edges_tmp(q)%edge_p(1:5) = (/p, edge_p(1:4)/)
							edge_vs_element_tmp(p)%edge_p = (/s, edge_p(1:4)/)
						else !for the case of j = 5, and 10					
							nr_edges_tmp(q)%edge_p(1:5) = (/p, j, 0, 0, 0/)
						end if
					end do !Loop i-n_CommonEdge
				end if
			end do
			
			!edge_II
			allocate(element_vs_edge(s)%nr_edges(q))
			do i = 1, q
				element_vs_edge(s)%nr_edges(i)%edge_p = nr_edges_tmp(i)%edge_p		
			end do			
		end do !end of s-loop
      
		!edge_I
      N_edge = p
		allocate(edge_vs_element(n_edge))  		
		edge_vs_element(1:n_edge) = edge_vs_element_tmp(1:n_edge)
		deallocate(nr_edges_tmp, edge_vs_element_tmp)
      return		
	end subroutine Sorting_edge
	! 
	subroutine find_edges_in_element(a_edge, edge_p)
		integer, intent(inout) ::  a_edge(7)
		integer, intent(inout) :: edge_p(4) 
		
		integer :: edge_tmp(2), edge(10, 2) 
		integer :: pI_a, pII_a, j, jj
		integer :: factor_II, factor_I
		
		edge(1, 1:2) = (/1, 2/) ! 
      edge(2, 1:2) = (/2, 3/) ! 
      edge(3, 1:2) = (/7, 6/) ! 
      edge(4, 1:2) = (/6, 5/)
      edge(5, 1:2) = (/8, 4/)
      edge(6, 1:2) = (/1, 8/)
      edge(7, 1:2) = (/8, 7/)
      edge(8, 1:2) = (/3, 4/)
      edge(9, 1:2) = (/4, 5/)
      edge(10, 1:2) = (/2, 6/)
		
		jj = 0					
		!This is to change the order of the nodes (1, 7, 8) into (1, 8, 7)
		j = edge_p(1)		
      if ((A_edge(2) .eq. Edge(j, 1)) .and. (A_edge(3) .eq. Edge(j, 2))) then
         pI_a = j !Edge index of element-s
         factor_I = 1
			edge_tmp = (/A_edge(5), A_edge(6)/)
			call sorting_edge_sub(edge, edge_tmp)
			jj = a_edge(1)
			pII_a = edge_tmp(1)
			factor_II = edge_tmp(2)	
      else if ((A_edge(3) .eq. Edge(j, 1)) .and. (A_edge(2) .eq. Edge(j, 2))) then 
         pI_a = j                        
         factor_I = -1
			edge_tmp = (/A_edge(5), A_edge(6)/)
			call sorting_edge_sub(edge, edge_tmp)
			jj = a_edge(1)
			pII_a = edge_tmp(1)
			factor_II = edge_tmp(2)
      else if ((A_edge(3) .eq. Edge(j, 1)) .and. (A_edge(4) .eq. Edge(j, 2))) then 
         pI_a = j !
         factor_I = 1 !
			edge_tmp = (/A_edge(6), A_edge(7)/)
			call sorting_edge_sub(edge, edge_tmp)
			jj = a_edge(1)
			pII_a = edge_tmp(1)
			factor_II = edge_tmp(2)
      else if ((A_edge(4) .eq. Edge(j, 1)) .and. (A_edge(3) .eq. Edge(j, 2))) then
         pI_a = j
         factor_I = -1
			edge_tmp = (/A_edge(6), A_edge(7)/)
			call sorting_edge_sub(edge, edge_tmp)
			jj = a_edge(1)
			pII_a = edge_tmp(1)
			factor_II = edge_tmp(2)				
      end if
		if ( (jj .NE. 0)) then   ! 
			edge_p = (/pI_a, jj, pII_a, factor_I*factor_II/)		
		end if		
		return
	end subroutine
	 
	subroutine sorting_edge_sub(edge, edge_tmp)
		implicit none
		integer, intent(in) :: Edge(10, 2)		
		integer, intent(inout) :: edge_tmp(2)
		
		integer :: kn
		
		do kn = 1, 10
			if ((edge_tmp(1) .eq. edge(kn, 1)) .and. (edge_tmp(2) .eq. edge(kn, 2))) then				
				edge_tmp(1) = kn
				edge_tmp(2) = 1
				exit
			else if ((edge_tmp(2) .eq. edge(kn, 1)) .and. (edge_tmp(1) .eq. edge(kn, 2))) then
				edge_tmp(1) = kn
				edge_tmp(2) = -1
				exit
			end if
		end do !Loop kn
		return
		
	 end subroutine
    
    subroutine Neighbours_EdgeCorner(struct, N_element)
    !Find the elements with common edges  
        integer, intent(in) :: N_element !Ne for sphere
        integer :: nn_e, nn_c, i, j, out_arr(3), N_gb1(3), N_gb2(3) 		  
		  
        type(structure_quad) :: struct
        type(Edge_Adjacent), dimension(:), allocatable :: temp_adjacent_edge    
        type(Vertex_Adjacent), dimension(:), allocatable :: temp_adjacent_corner
   
        allocate(temp_adjacent_edge(N_element))      
        allocate(temp_adjacent_corner(4*N_element))
        do i = 1, N_element
            nn_e = 0  
            nn_c = 0  
            do j = i + 1, N_element   
                out_arr = Common_EdgeCorner(struct%elements(i)%Vertices, &
                struct%elements(j)%Vertices, N_gb1, N_gb2)                 
                if ( out_arr(3) .ne. -1 ) then !Edge overlap
                    nn_e = nn_e + 1           
                    temp_adjacent_edge(nn_e)%CommonEdge_Elements = (/i, j, N_gb1(1), N_gb1(2), N_gb1(3), N_gb2(1), N_gb2(2), N_gb2(3)/) 
                    temp_adjacent_edge(nn_e)%CommonEdge_Nodes = out_arr(1:3) !
                else if ((out_arr(1) .ne. -1) .and. (out_arr(3) .eq. -1)) then !corner overlap
                    nn_c = nn_c +1                     
                    temp_adjacent_corner(nn_c)%CommonCorner_Elements = (/i, j, N_gb1(1), N_gb1(2), N_gb2(1), N_gb2(2)/)
                    temp_adjacent_corner(nn_c)%CommonCorner = out_arr(1) 
                end if      
            end do         
            allocate(struct%elements(i)%edge_overlap(nn_e))
            allocate(struct%elements(i)%corner_overlap(nn_c)) 
            do j = 1, nn_e
                struct%elements(i)%edge_overlap(j) = temp_adjacent_edge(j)
            end do
            do j = 1, nn_c
                struct%elements(i)%corner_overlap(j) = temp_adjacent_corner(j)
            end do
        end do !Loop i
        return
    end subroutine Neighbours_EdgeCorner   
    !
    subroutine Neighbours_EdgeCorner_Full(struct, N_element)
    !Find the elements with common edges  
        integer, intent(in) :: N_element !
        integer :: nn_e, nn_c, i, j, out_arr(3), N_gb1(3), N_gb2(3) !
        type(structure_quad) :: struct
        type(Edge_Adjacent), dimension(:), allocatable :: temp_adjacent_edge    
        type( Vertex_Adjacent), dimension(:), allocatable :: temp_adjacent_corner
   
        allocate(temp_adjacent_edge(N_element))      
        allocate(temp_adjacent_corner(4*N_element))
        do i = 1, N_element
            nn_e = 0  
            nn_c = 0  
            do j = 1, N_element 					
                if (j .NE. i) then				
                    out_arr = Common_EdgeCorner(struct%elements(i)%Vertices, &
									struct%elements(j)%Vertices, N_gb1, N_gb2)
						  
                    if ( out_arr(3) /= -1 ) then !Edge overlap
                        nn_e = nn_e + 1           
                        temp_adjacent_edge(nn_e)%CommonEdge_Elements = (/i, j, N_gb1(1), N_gb1(2), N_gb1(3), N_gb2(1), N_gb2(2), N_gb2(3)/) !
                        temp_adjacent_edge(nn_e)%CommonEdge_Nodes = out_arr(1:3) !It seems it is useless
                    else if ((out_arr(1) /= -1) .and. (out_arr(3) == -1)) then !corner overlap
                        nn_c = nn_c +1                     
                        temp_adjacent_corner(nn_c)%CommonCorner_Elements = (/i, j, N_gb1(1), N_gb1(2), N_gb2(1), N_gb2(2)/)
                        temp_adjacent_corner(nn_c)%CommonCorner = out_arr(1) 
                    end if  
                else
                end if                
            end do         
            allocate(struct%elements(i)%edge_overlap(nn_e))
            allocate(struct%elements(i)%corner_overlap(nn_c)) 
   
            do j = 1, nn_e
                struct%elements(i)%edge_overlap(j) = temp_adjacent_edge(j)
            end do
            do j = 1, nn_c
                struct%elements(i)%corner_overlap(j) = temp_adjacent_corner(j)
            end do
        end do !Loop i  
		  deallocate(temp_adjacent_edge, temp_adjacent_corner)
        return
    end subroutine Neighbours_EdgeCorner_Full 
         
    function Common_EdgeCorner(arr_a, arr_q, Nr_a, Nr_q) ! Check whether arr_a and arr_q have two common edges or corners
        integer, intent(out) :: Nr_a(3), Nr_q(3)
        integer :: i, j, n_equal, Common_EdgeCorner(3)
        integer, intent(in) :: arr_a(8), arr_q(8) !node numbers of elements, ddi: element number, arr has 8 nodes.
        Common_EdgeCorner = (/-1, -1, -1/) 
        n_equal = 0  
        Nr_a = (/0, 0, 0/)
        Nr_q = (/0, 0, 0/)      
        do i = 1, 8
            do j = 1, 8
                if (arr_a(i) == arr_q(j)) then  
                    n_equal = n_equal + 1            
                    Common_EdgeCorner(n_equal) = arr_q(j) 
                    if (n_equal == 1) then 
                       Nr_a(1) = i !node-index of the overlapped corner or edge in the first element                               
                       Nr_q(1) = j !node index in the second element
                    else if (n_equal == 2) then 
                       Nr_a(2) = i                   
                       Nr_q(2) = j
                    else if (n_equal == 3) then 
                       Nr_a(3) = i 
                       Nr_q(3) = j                  
                    end if
                end if         
            end do      
        end do
        return
    end function  Common_EdgeCorner    
      
    function rho_2(mu, nu, phi, m)
        real(kind = dp), intent(in) :: mu, nu, phi
        integer, intent(in) :: m
        real(kind = dp) :: rho_2a(8), rho_2        
            rho_2a(1) = abs((1-mu)/cos(phi))
            rho_2a(2) = ((1-nu)/sin(PI-phi))        
            rho_2a(3) = abs(rho_2a(2))
            rho_2a(4) = abs((1+mu)/cos(PI-phi))
            rho_2a(5) = rho_2a(4)         
            rho_2a(6) = abs((1+nu)/cos(1.5*PI-phi))
            rho_2a(7) = rho_2a(6)                 
            rho_2a(8) = rho_2a(1)          
            rho_2 = rho_2a(m)
        return 
    end function rho_2
   
    function rho_2_4Feld(mu, nu, phi, m)
        real(kind = dp), intent(in) :: mu, nu, phi
        integer, intent(in) :: m
        real(kind = dp) :: rho_2a(4), rho_2_4Feld
            rho_2a(1) = ((1-nu)/sin(PI-phi))
            rho_2a(2) = abs((1+mu)/cos(PI-phi))
            rho_2a(3) = abs((1+nu)/cos(1.5*PI-phi))
            rho_2a(4) = abs((1-mu)/cos(phi))
            rho_2_4Feld = rho_2a(m)
        return 
    end function rho_2_4Feld
    
    function rho_2_3Feld(mu, phi, m)
        real(kind = dp), intent(in) :: mu, phi
        integer, intent(in) :: m
        real(kind = dp) :: rho_2a(3), rho_2_3Feld
            rho_2a(1) = (1-mu)/cos(phi)
            rho_2a(2) = 2/sin(phi)
            rho_2a(3) = (1+mu)/cos(pi- phi)
            rho_2_3Feld = rho_2a(m)
        return 
    end function rho_2_3Feld

    function phi_range(mu, nu)
      real(kind = dp), intent(in) :: mu, nu   
      real(kind = dp) :: phi_range(9)
         phi_range(1) = 0.0*PI  
         phi_range(3) = PI/2
         phi_range(5) = PI
         phi_range(7) = 3*PI/2
         phi_range(9) = 2*PI 
         phi_range(2) = atan((1-nu)/(1-mu))   
         phi_range(4) = PI - atan((1-nu)/(1+mu))                    
         phi_range(6) = PI + atan((1+nu)/(1+mu))
         phi_range(8) = 2*PI - atan((1+nu)/(1-mu)) 
      return     
   end function  phi_range

   function phi_range_4Feld(mu, nu)
      real(kind = dp), intent(in) :: mu, nu   
      real(kind = dp) :: phi_range_4Feld(5)
         phi_range_4Feld(1) = atan((1-nu)/(1-mu))   
         phi_range_4Feld(2) = PI - atan((1-nu)/(1+mu))                    
         phi_range_4Feld(3) = PI + atan((1+nu)/(1+mu))
         phi_range_4Feld(4) = 2*PI - atan((1+nu)/(1-mu)) 
         phi_range_4Feld(5) = 2*PI + atan((1-nu)/(1-mu))   
      return     
   end function  phi_range_4Feld
   
   function phi_range_3Feld(mu)
      real(kind = dp), intent(in) :: mu
      real(kind = dp) :: phi_range_3Feld(4)        
         phi_range_3Feld(1) = 0
         phi_range_3Feld(2) = atan(2/(1-mu))   
         phi_range_3Feld(3) = PI - atan(2/(1+mu)) 
         phi_range_3Feld(4) = PI 
      return     
   end function  phi_range_3Feld
   
   
	SUBROUTINE M22INV (A, AINV, OK_FLAG)
      implicit none
      real(kind = dp), dimension(2,2), intent(in)  :: A
      real(kind = dp), dimension(2,2), intent(out) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      real(kind = dp), parameter :: EPS = 1.0D-10
      real(kind = dp) :: DET
      real(kind = dp), dimension(2,2) :: COFACTOR

      DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
      if (ABS(DET) .LE. EPS) then
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         return
      end if

      COFACTOR(1,1) = +A(2,2)
      COFACTOR(1,2) = -A(2,1)
      COFACTOR(2,1) = -A(1,2)
      COFACTOR(2,2) = +A(1,1)

      AINV = transpose(COFACTOR) / DET
      OK_FLAG = .TRUE.
		return
	end subroutine M22INV	
	
	subroutine incident_int(Element_aa, pI_aa, ngp_a, Sum_EH)
        complex(kind = dp), intent(out) ::  Sum_EH(2)
        complex(kind = dp) :: Ec_in(3), Hc_in(3)        
        integer, intent(in) :: pI_aa, ngp_a
        integer i, j, n
        real(kind = dp), dimension(10, 3) :: Nedge_a
        real(kind = dp) :: R_a(3), Eu_a(3), Ev_a(3), Nod_d(8,2), TRF
        real(kind = dp) :: eta_a, ksi_a, Naj_v(3)
        real(kind = dp), intent(in) :: Element_aa(8,3)

        real(kind = dp), dimension(100) :: x, w
			
        call xw_GPoint_1D(x, w, ngp_a)
        Sum_EH(1:2) = (/(0.0, 0.0),  (0.0, 0.0)/)
        do j = 1, ngp_a     !Loop for integration of Element_a  
            eta_a = x(j)
            do i = 1, ngp_a  
                ksi_a = x(i) 
                Nod_d = N_differential(ksi_a, eta_a)
                do n = 1, 3
                    R_a(n)  = dot_product(R_shape(ksi_a, eta_a),  Element_aa(1:8, n)) 
                    Eu_a(n) = dot_product(Nod_d(1:8, 1),  Element_aa(1:8, n))
                    Ev_a(n) = dot_product(Nod_d(1:8, 2),  Element_aa(1:8, n)) 
                end do 
                Nedge_a = Edge_vector_new(ksi_a, eta_a, Element_aa)
                TRF = w(i)*w(j)*vec_len(cross_r(Eu_a, Ev_a))
                Ec_in = illumination_p%E_in*exp(-im*dot_product(k1*illumination_p%k_in, R_a)) !
                Hc_in = cross_c(illumination_p%k_in*(1+im*0.0), Ec_in)/(my_1*c0)
                Naj_v = Nedge_a(pI_aa, :)
                Sum_EH(1) = Sum_EH(1) + dot_product(Naj_v, Ec_in)*TRF
                Sum_EH(2) = Sum_EH(2) + dot_product(Naj_v, Hc_in)*TRF      
            end do
        end do   
    end subroutine incident_int    
    
	subroutine incident_Gauss_quad(Element_aa, pI_aa, ngp_a, Sum_EH)
		complex(kind = dp), intent(out) ::  Sum_EH(2)
		
      type(vector_c) :: Ec_in, Hc_in
      integer, intent(in) :: pI_aa, ngp_a
      integer i, j, n
		real(kind = dp), dimension(10, 3) :: Nedge_a
		real(kind = dp) :: R_a(3), Eu_a(3), Ev_a(3), Nod_d(8,2), TRF
		real(kind = dp) :: eta_a, ksi_a, Naj_v(3)
		real(kind = dp), intent(in) :: Element_aa(8,3)

		real(kind = dp), dimension(100) :: x, w
		complex(kind = dp) :: k_hat(3)
				
      k_hat = illumination_p%k_in !unit vector
      call xw_GPoint_1D(x, w, ngp_a)
      Sum_EH(1:2) = (/(0.0, 0.0),  (0.0, 0.0)/)
      do j = 1, ngp_a     !Loop for integration of Element_a  
         eta_a = x(j)
         do i = 1, ngp_a  
            ksi_a = x(i) 
            Nod_d = N_differential(ksi_a, eta_a)
            do n = 1, 3
               R_a(n)  = dot_product(R_shape(ksi_a, eta_a),  Element_aa(1:8, n)) 
               Eu_a(n) = dot_product(Nod_d(1:8, 1),  Element_aa(1:8, n))
               Ev_a(n) = dot_product(Nod_d(1:8, 2),  Element_aa(1:8, n)) 
            end do 
				!print*, 'beam_waist', beam_waist
				call Gaussian_beam(beam_waist, R_a, illumination_p, Ec_in, k1) ! 
            Nedge_a = Edge_vector_new(ksi_a, eta_a, Element_aa)
            TRF = w(i)*w(j)*vec_len(cross_r(Eu_a, Ev_a))
            Hc_in%vector = cross_c(k_hat, Ec_in%vector)/(my_1*c0)
            Naj_v = Nedge_a(pI_aa, :)
            Sum_EH(1) = Sum_EH(1) + dot_product(Naj_v, Ec_in%vector)*TRF
            Sum_EH(2) = Sum_EH(2) + dot_product(Naj_v, Hc_in%vector)*TRF      
         end do
      end do   
	end subroutine incident_Gauss_quad    
	
   !integration of the field in free space when SE and SH have been calculated
	subroutine Integration_scattered_field_quad(Ra, Element_qq,  ngp_q, SE, SH, Sum_aa)        
		complex(kind = dp), intent(in) ::  SE(10), SH(10)        
		real(kind = dp), intent(in) :: Element_qq(8,3), Ra(3)
		integer, intent(in) :: ngp_q
		type(Vector_c), dimension(2), intent(out) :: Sum_aa
			
		!dummy
		type(point), dimension(:), allocatable :: r_local
      integer i, j, n
      real(kind = dp):: R_norm !
      complex(kind = dp) :: GradG(3), Greenf, GG_SS
      complex(kind = dp) :: Vect_temp(3), KA_temp_1(3), KA(3), KB(3), KC(3) 
       
      real(kind = dp) ::  ff,  Rq(3)! 
      real(kind = dp), dimension(100) :: w, x
		  
		type(derived_element_parameters_TypeB) :: delement_pq
      call xw_GPoint_1D(x, w, ngp_q)       
      do n = 1, 2
         Sum_aa(n)%Vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      end do
		  
      do i = 1,  ngp_q      
         do j = 1, ngp_q      
            call paras_st_B(x(i), x(j), Element_qq, delement_pq)!
				Rq = R_local_parametric(x(i), x(j), Element_qq) 		
            R_norm = vec_len(Ra-Rq)
            call Diver_Green(k1, R_norm, Greenf, GG_SS)
            GradG = (Ra-Rq)/R_norm*GG_SS 
            do n = 1, 10 ! Edge loop
               KA_temp_1 = cross_r(delement_pq%nq%point, delement_pq%Nqi(n)%point)
               Vect_temp = KA_temp_1
					ff = delement_pq%JJ*w(i)*w(j)!
               KA = cross_c(Vect_temp, GradG)*ff
               KB = delement_pq%n_rot_n(n)*GradG*ff
               KC = KA_temp_1*Greenf*ff                    
               Sum_aa(1)%vector = Sum_aa(1)%vector + KA*SE(n) + &
							1.0_dp/(im*Omega*eps_0)*KB*SH(n) - (-im*Omega*my_0)*KC*SH(n)   !sum for edges! !KA_E
               Sum_aa(2)%vector = Sum_aa(2)%vector + KA*SH(n) + &
							1.0_dp/(-im*Omega*my_0)*KB*SE(n) - im*Omega*eps_0*KC*SE(n)  !KB_E
            end do                 
         end do !loop j
      end do !loop i       
      return
   end subroutine Integration_scattered_field_quad
   
end module lib_sie_quad_calculation_mod