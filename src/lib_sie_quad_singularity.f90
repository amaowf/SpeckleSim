module lib_sie_quad_singularity
	
   use lib_sie_math
	use lib_sie_constants
	use lib_sie_type_function
   use lib_sie_integration_points
   use lib_sie_quad_calculation_mod
	use lib_sie_quad_data_container    
    
	implicit none
	 
	public :: impedance_edge_versus_element
	 
	type edge_sum_area
		complex(kind=dp), dimension(3, 10) :: edge_area 
	end type 
	
	private
	
	contains
	
	subroutine impedance_edge_versus_element(s_edge, t_element, sum_aa)
		implicit none
		integer, intent(in) :: s_edge, t_element
		complex(dp), intent(out) :: sum_aa(3, 10) 
		complex(dp) :: sum_a(3, 10)
		integer :: n
		integer :: ii, jj
		integer :: n_CommonEdge, n_CommonCorner
		integer :: pI_a, pI_aa, pI_q, Edge_Rotated_a(2, 10), Edge_Rotated_q(2, 10)    
		integer :: ngp, loc_corner(3),  t_edge_index(2), t_edge(6)
		integer, allocatable, dimension(:) :: t_corner
    
		!Parameters used during the matrix calculation
		real(kind = dp), dimension(8, 3) :: Element_a, Element_q, Element_aa, Element_qq!
		character(len = 20) :: type_calculation
		
		allocate(t_corner(5))
		
		ii = Edge_I(s_edge)%edge_p(1)! element index of the s-th edge !
		pI_a = Edge_I(s_edge)%edge_p(2) ! edge index in element-ii
		pI_q = Edge_I(s_edge)%edge_p(4) ! 	
		n_commonedge = size(struc_quad%elements(ii)%edge_overlap)
		do n = 1, n_CommonEdge
			t_edge = struc_quad%elements(ii)%edge_overlap(n)%CommonEdge_Elements(3:8)			
			jj = struc_quad%elements(ii)%edge_overlap(n)%CommonEdge_Elements(2)
			if (jj .eq. t_element) then
				call find_edge_index(t_edge, t_edge_index)
				exit
			endif  
		end do
					
		!The number of elements having overlapped corners with element-ii	
		n_CommonCorner = size(struc_quad%elements(ii)%corner_overlap)
		do n = 1, n_CommonCorner
			t_corner(1:5) = &
			struc_quad%elements(ii)%corner_overlap(n)%CommonCorner_Elements(2:6)!
			if (t_corner(1) .eq. t_element) then
				exit
			endif
		end do
		
		do n = 1, 8 
			Element_a(n, 1:3) = struc_quad%node_vector(struc_quad%elements(ii)%Vertices(n))%point(1:3)
		end do
											
		do n = 1, 8
			Element_q(n, 1:3) = struc_quad%node_vector(struc_quad%elements(t_element)%Vertices(n))%point(1:3)
		end do
		if (t_element .eq. ii) then 
		
			type_calculation = 'SelfTerm'               
			ngp = 15						
			call Element_integration(element_a, element_q, pI_a, type_calculation, ngp, Sum_aa) 			
		 else if (t_element .eq. jj) then
		 
			type_calculation = 'EdgeA'
			ngp = 5			
			!The edge indices of the two elements have to be renumbered for integration
			call Element_Rotation_edge(t_edge_index(1), Element_a, Element_aa, Edge_Rotated_a)
			call Element_Rotation_edge(t_edge_index(2), Element_q, Element_qq, Edge_Rotated_q)                
			pI_aa = Edge_Rotated_a(1, pI_a) !The edge index of element_a in the roated element
			call Element_integration(element_aa, element_qq, pI_aa, type_calculation, ngp, Sum_a)
			do n = 1, 10
				Sum_aa(:, n) = Sum_a(:, Edge_Rotated_q(1, n))*Edge_Rotated_q(2, n)*Edge_Rotated_a(2, pI_a)
			end do
		else if (t_element .eq. t_corner(1)) then
		
			type_calculation = 'Normal' ! 'VertexA'!
			ngp = 6
			call Element_integration(element_a, element_q, pI_a, type_calculation, ngp, Sum_aa)  
		else !  
			type_calculation = 'Normal'					
			ngp = 4
			call Element_integration(element_a, element_q, pI_a, type_calculation, ngp, Sum_aa)
		end if
		return
	end subroutine
	 
	subroutine Element_integration(element_aa, element_qq, pI_aa, type_cal, ngp_a, Sum_area)
		implicit none   
		integer, intent(in) :: ngp_a, pI_aa
		  
		integer :: i, j, m, n, s, p, N_Field
		integer :: ll, n_max
		integer :: ngp_1, ngp_2, ngp_3, ngp_4         
		character(len = 10), intent(in) :: type_cal        
		real(kind = dp), intent(in) :: Element_aa(8, 3), Element_qq(8, 3)                
		real(kind = dp), dimension(100) :: xx, xw
		
		type(derived_element_parameters_TypeA) :: delement_pa
		real(kind = dp) :: TRF, Naj_v(3), R_a(3)
		real(kind = dp) :: eta_a, ksi_a
        
		complex(kind=dp), dimension(3, 10),  intent(out) :: Sum_area 
		complex(kind=dp), dimension(3, 10) :: Temp 
		  
		do j = 1, 3			
			Sum_area(j, 1:10) = (0.0, 0.0)				
			Temp(j, 1:10) = (0.0, 0.0)
		end do
		
		call xw_GPoint_1D(xx, xw, ngp_a)				
		select case (type_cal) 
		case ("SelfTerm") 
		!print*, 'self term'	
		!Loop for integration of Element_a  
			N_field = 4
			do j = 1, ngp_a 
				eta_a = xx(j)
				do i = 1, ngp_a
					ksi_a = xx(i)
					call paras_st_A(ksi_a, eta_a, Element_aa, delement_pa)   !
					Naj_v = delement_pa%Naj(pI_aa)%point
					TRF = xw(i)*xw(j)*delement_pa%JJ
					Sum_area = Sum_area + Integration_KABCD_Taylor_ST_4FEdge(pI_aa, ksi_a, eta_a, Element_aa, N_Field, ngp_a)*TRF
				end do !Loop i
			end do ! Loop j
		case ("EdgeA")
			ngp_4 = ngp_a
			ngp_3 = ngp_a!
			ngp_2 = ngp_a!
			ngp_1 = ngp_a!
			do s = 1, 2  !
				m = s - 1 
				do ll = 1, 3
					p = ll-1 !(p = l) 
					if (((p == 0) .and. (m == 0)) .or. ((p == 2) .and. (m == 0))) then
						n_max = 4
					else
						n_max = 2
					end if
					do n = 1, n_max
						Sum_area = Sum_area + &
							GL_Int_theta_EA_8Node(p, m, n, Element_aa, Element_qq, pI_aa, ngp_4, ngp_3, ngp_2, ngp_1)
					end do
				end do !Loop s    
			end do !Loop-i    
		case ('VertexA') 
			ngp_4 = ngp_a!
			ngp_3 = ngp_a!
			ngp_2 = ngp_a!
			ngp_1 = ngp_a!
			do m = 1, 2 !
				do  n = 1, 2
					do p = 1, 2                
					Sum_area = Sum_area + &
						GL_Int_ThetaP_VA_8Node(m, n, p, Element_aa, Element_qq, pI_aa, ngp_4, ngp_3, ngp_2, ngp_1) 
					end do
				end do
			end do
		case default !normal
			!print*, 'normal integration'                
			do j = 1, ngp_a !Loop for integration of Element_a  
				eta_a = xx(j)
				do i = 1, ngp_a
					ksi_a = xx(i)
					call paras_st_A(ksi_a, eta_a, Element_aa, delement_pa)                        
					Naj_v = delement_pa%Naj(pI_aa)%point !a loop can be built to save time. 
					TRF = xw(i)*xw(j)*delement_pa%JJ
					R_a = r_local_parametric(ksi_a, eta_a, Element_aa)
					Sum_area = Sum_area + &
						Integration_KABCD_Edge_new(Naj_v, Element_qq, R_a, ngp_a)*TRF
				end do !Loop-i
			end do !Loop-j    
		end select
		return 
	end subroutine Element_integration 
	
	!2021.01.13 to accelerate the calculation
	subroutine Element_integration_new(element_aa, element_qq, edge_arr, type_cal, ngp_a, Sum_area)
		implicit none   
		integer, intent(in) :: ngp_a 
		integer, dimension(:), intent(in) :: edge_arr
		
		!dummy
		integer :: i, j, m, n, s, p, N_Field
		integer :: ll, n_max, pI_aa, n_edge
		integer :: ngp_1, ngp_2, ngp_3, ngp_4         
		character(len = 10), intent(in) :: type_cal        
		real(kind = dp), intent(in) :: Element_aa(8, 3), Element_qq(8, 3)                
		real(kind = dp), dimension(100) :: xx, xw
		
		type(derived_element_parameters_TypeA) :: delement_pa
		real(kind = dp) :: TRF, Naj_v(3), R_a(3)
		real(kind = dp) :: eta_a, ksi_a
      
		type(edge_sum_area), dimension(:), allocatable, intent(out) :: sum_area

		complex(kind=dp), dimension(3, 10) :: Temp 
		
		allocate(sum_area(size(edge_arr)))
		do j = 1, 10
			Sum_area(j)%edge_area(1:3, 1:10) = (0.0, 0.0)							
		end do
		Temp(1:3, 1:10) = (0.0, 0.0)
		
		call xw_GPoint_1D(xx, xw, ngp_a)				
		select case (type_cal) 
		case ("SelfTerm") 
		!print*, 'self term'	
		!Loop for integration of Element_a  
			N_field = 4
			do j = 1, ngp_a 
				eta_a = xx(j)
				do i = 1, ngp_a
					ksi_a = xx(i)
					call paras_st_A(ksi_a, eta_a, Element_aa, delement_pa)   !
					do n_edge = 1, size(edge_arr)
						pI_aa = edge_arr(n_edge)
						Naj_v = delement_pa%Naj(pI_aa)%point
						TRF = xw(i)*xw(j)*delement_pa%JJ
						Sum_area(pI_aa)%edge_area = Sum_area(pI_aa)%edge_area + &
						 Integration_KABCD_Taylor_ST_4FEdge(pI_aa, ksi_a, eta_a, Element_aa, N_Field, ngp_a)*TRF
					end do	
				end do !Loop i
			end do ! Loop j
		case ("EdgeA")
			ngp_4 = ngp_a
			ngp_3 = ngp_a!
			ngp_2 = ngp_a!
			ngp_1 = ngp_a!
			do n_edge = 1, size(edge_arr)
				pI_aa = edge_arr(n_edge)			
				do s = 1, 2  !
					m = s - 1 
					do ll = 1, 3
						p = ll-1 !(p = l) 
						if (((p == 0) .and. (m == 0)) .or. ((p == 2) .and. (m == 0))) then
							n_max = 4
						else
							n_max = 2
						end if
						do n = 1, n_max
							Sum_area(pI_aa)%edge_area = Sum_area(pI_aa)%edge_area + &
								GL_Int_theta_EA_8Node(p, m, n, Element_aa, Element_qq, pI_aa, ngp_4, ngp_3, ngp_2, ngp_1)
						end do
					end do !Loop s    
				end do !Loop-i  
			end do !loop pI_aa
		case ('VertexA') 
			ngp_4 = ngp_a!
			ngp_3 = ngp_a!
			ngp_2 = ngp_a!
			ngp_1 = ngp_a!
			do n_edge = 1, size(edge_arr)
				pI_aa = edge_arr(n_edge)
				do m = 1, 2 !
					do  n = 1, 2
						do p = 1, 2                
						Sum_area(pI_aa)%edge_area = Sum_area(pI_aa)%edge_area + &
							GL_Int_ThetaP_VA_8Node(m, n, p, Element_aa, Element_qq, pI_aa, ngp_4, ngp_3, ngp_2, ngp_1) 
						end do
					end do
				end do
			end do
		case default !normal
			!print*, 'normal integration'                
			do j = 1, ngp_a !Loop for integration of Element_a  
				eta_a = xx(j)
				do i = 1, ngp_a
					ksi_a = xx(i)
					call paras_st_A(ksi_a, eta_a, Element_aa, delement_pa)  
					do n_edge = 1, size(edge_arr)
						pI_aa = edge_arr(n_edge)
						Naj_v = delement_pa%Naj(pI_aa)%point !a loop can be built to save time. 
						TRF = xw(i)*xw(j)*delement_pa%JJ
						R_a = r_local_parametric(ksi_a, eta_a, Element_aa)
						Sum_area(pI_aa)%edge_area = Sum_area(pI_aa)%edge_area + &
							Integration_KABCD_Edge_new(Naj_v, Element_qq, R_a, ngp_a)*TRF
					end do
				end do !Loop-i
			end do !Loop-j    
		end select
		return 
	end subroutine Element_integration_new    
   
	!self term
	function Integration_KABCD_Taylor_ST_4FEdge(pI_a, mu_a, nu_a, Element_q, N_Field, ngp)result(rv)
		real(kind = dp), intent(in) :: mu_a, nu_a
		real(kind = dp), dimension(8, 3), intent(in) :: Element_q  
		integer, intent(in) :: ngp, pI_a, N_Field
      
		!dummy
		integer :: i, j, m, ngp_1, ngp_2        
		real(kind = dp) :: A_vec(3),  a, c, A_amp, Naj_v(3) !JJ_q, JJ_a, n_q(3), n_a(3),
		real(kind = dp), dimension(100) :: x1, w1, x2, w2
		complex(kind = dp) :: ff(5, 10), rv(3, 10)
		complex(kind = dp), dimension(10) :: KA_a, KA_b, KBE_a, KBH_a, KBE_b, KBH_b, KCH, KCE  
		  
		real(kind = dp) :: R_q(3), rho2, rho, phi_b(5), phi, f_AB(2, 10), R_a(3) 
		real(kind = dp) ::  R(3), Eu_a(3), Ev_a(3)!
		real(kind = dp) ::  mu_q, nu_q, TRF, TRF_0  		  
		real(kind = dp) :: f_diff(8, 2)
		  
		type(derived_element_parameters_TypeB):: derived_pq, derived_pa
		  
		f_diff = N_differential(mu_a, nu_a)
       
      do i = 1, 10 !Initialize the area array        
         KA_a(i) = (0.0, 0.0)
         KA_b(i) = (0.0, 0.0)
         KBE_a(i) = (0.0, 0.0)
         KBE_b(i) = (0.0, 0.0)       
         KBH_a(i) = (0.0, 0.0)
         KBH_b(i) = (0.0, 0.0)       
         KCE(i) = (0.0, 0.0)
         KCH(i) = (0.0, 0.0)      
      end do
      !'a': observation point, 'q': source point.
      phi_b = phi_range_4Feld(mu_a, nu_a)
      ngp_1 = ngp
      call xw_GPoint_1D(x1, w1, ngp_1)        
      ngp_2 = 5       
      call xw_GPoint_1D(x2, w2, ngp_2)        
		call paras_st_B(mu_a, nu_a, Element_q, derived_pa)
		  
		R_a = r_local_parametric(mu_a, nu_a, Element_q)
			
		do m = 1, 3           
			Eu_a(m) = dot_product(f_diff(1:8, 1),  Element_q(1:8, m)) !dr/dksi
			Ev_a(m) = dot_product(f_diff(1:8, 2),  Element_q(1:8, m)) !dr/deta
		end do
			
		Naj_v = derived_pa%Nqi(pI_a)%point
      do m = 1, N_Field!
         a = (phi_b(m+1) - phi_b(m))/2 
         c = (phi_b(m+1) + phi_b(m))/2 			
         do i = 1, ngp_1
               phi = c + a*x1(i)  !Eq.(4.12)
               A_vec = - Eu_a*cos(phi) - Ev_a*sin(phi)!
               A_amp = vec_len(A_vec)                
               rho2 = rho_2_4Feld(mu_a, nu_a, phi, m)              
					f_AB = ff_phi_AB_Edge(Naj_v, A_vec, derived_pa)!
               do j = 1, ngp_2 !integration over q-element                      
                  rho = rho2/2*(1+x2(j))
                  mu_q = mu_a + rho * cos(phi)  
                  nu_q = nu_a + rho * sin(phi)
                  call paras_st_B(mu_q, nu_q, Element_q, derived_pq)	!Eu_q, Ev_q,
						R_q = r_local_parametric(mu_q, nu_q, Element_q)
                  R = R_a - R_q  !
                  TRF = 0.5*rho2*w1(i)*w2(j)*a                                       
						ff = ff_ABC_II_Edge(Naj_v, R, derived_pq)
                  KA_a = KA_a + (ff(1, :)*rho - 2*f_AB(1, :)/rho)*TRF!
                  KBE_a = KBE_a + (ff(2, :)*rho - (1 + 1/eps_r2)*f_AB(2, :)/rho)*TRF!!!!!!!
                  KBH_a = KBH_a + (ff(3, :)*rho - 2*f_AB(2, :)/rho)*TRF!2 !!!!!!                    
                  KCE = KCE + ff(4, :)*rho*TRF
                  KCH = KCH + ff(5, :)*rho*TRF
               end do ! Loop j 
               TRF_0 = log(abs(rho2)*A_amp)*w1(i)*a
               KA_b = KA_b + 2*f_AB(1, :)*TRF_0!  
               KBE_b = KBE_b + (1+1/eps_r2)*f_AB(2, :)*TRF_0 !!
               KBH_b = KBH_b + 2*f_AB(2, :)*TRF_0 
         end do ! Loop i
      end do ! Loop m 
      rv(1, 1:10) = (KA_a + KA_b)!
      rv(2, 1:10) = (KBE_a + KBE_b)/(im*Omega*eps_0)-(-im*Omega*my_0*KCE)!  
      rv(3, 1:10) = (KBH_a + KBH_b)/(-im*Omega*my_0)-im*Omega*eps_0*KCH! 				  
      return
	end function Integration_KABCD_Taylor_ST_4FEdge   
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Element rotation for calculating vertex and edge overlapped cases
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
	subroutine Element_rotation_edge(qI, Element_in, Element_out, Edge_Rotation)
      real(kind = dp), intent(in):: Element_in(8, 3) 
      real(kind = dp), intent(out) :: Element_out(8, 3)
      integer, intent(out) :: Edge_Rotation(2, 10)
      integer, intent(in) :: qI 
      integer :: j, jj
      !print*, 'qI in the subroutine = ', qI
      if ((qI == 3) .or. (qI == 4)) then
         do j = 1, 8                
               if (j<= 4) then
                  jj = 4 + j !jj is the output index
                  Element_out(jj, :) = Element_in(j, :)
               else 
                  jj = j - 4
                  Element_out(jj, :) = Element_in(j, :)
               end if                
         end do
         Edge_Rotation(1, :) = (/4, 3, 2, 1, 5, 9, 8, 7, 6, 10/) !edge index of the roated
         Edge_Rotation(2, :) = (/-1, -1, -1, -1, -1, -1, -1, -1, -1, -1/) !change of the edge direction
            
      else if ((qI == 6) .or. (qI == 7)) then
         do j = 1, 8
               if (j < 7) then
                  jj = j + 2
                  Element_out(jj, :) = Element_in(j, :)
               else 
                  jj = j - 6
                  Element_out(jj, :) = Element_in(j, :)
               end if
         end do
         Edge_Rotation(1, :) = (/8, 9, 6, 7, 10, 2, 1, 4, 3, 5/) !The edge index of the roated element
         Edge_Rotation(2, :) = (/1, 1, 1, 1, 1, -1, -1, -1, -1, -1/)
            
      else if ((qI == 8) .or. (qI == 9)) then
         do j = 1, 8
               if (j > 2) then
                  jj = j - 2
                  Element_out(jj, :) = Element_in(j, :)
               else 
                  jj = j + 6
                  Element_out(jj, :) = Element_in(j, :)
               end if
         end do
         Edge_Rotation(1, :) = (/7, 6, 9, 8, 10, 3, 4, 1, 2, 5/)
         Edge_Rotation(2, :) = (/-1, -1, -1, -1, -1, 1, 1, 1, 1, 1/)
            
      else 
         Element_out = Element_in
         Edge_Rotation(1, :) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/) 
         Edge_Rotation(2, :) = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
      end if
      return        
   end subroutine Element_rotation_edge
    
   subroutine Element_Rotation_corner(cI, Element_in, Element_out, Edge_Rotation)    
        real(kind = dp), intent(in):: Element_in(8, 3) 
        real(kind = dp), intent(out) :: Element_out(8, 3)
        integer, intent(out) :: Edge_Rotation(2, 10)
        integer, intent(in) :: cI
        integer :: j, jj
        !print*, 'qI in the subroutine = ', qI         
        if (cI .NE. 1) then
            do j = 1, 8
                if (j <= (8 - cI + 1)) then
                    jj = j + cI - 1
                    Element_out(j, :) = Element_in(jj, :)
                else 
                    jj = j - (8 - cI + 1)
                    Element_out(j, :) = Element_in(jj, :)
                end if
            end do
        else 
            Element_out = Element_in
        end if
        if (cI == 3) then
            Edge_Rotation(1, :) = (/7, 6, 9, 8, 10, 3, 4, 1, 2, 5/)
            Edge_Rotation(2, :) = (/-1, -1, -1, -1, -1, 1, 1, 1, 1, 1/)            
        else if (cI == 5) then            
            Edge_Rotation(1, :) = (/4, 3, 2, 1, 5, 9, 8, 7, 6, 10/)
            Edge_Rotation(2, :) = (/-1, -1, -1, -1, -1, -1, -1, -1, -1, -1/)            
        else if (cI == 7) then
            Edge_Rotation(1, :) = (/8, 9, 6, 7, 10, 2, 1, 4, 3, 5/)
            Edge_Rotation(2, :) = (/1, 1, 1, 1, 1, -1, -1, -1, -1, -1/)
        else
            Edge_Rotation(1, :) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/) 
            Edge_Rotation(2, :) = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
        end if        
        return        
    end subroutine Element_rotation_corner
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 
	 function Integration_KABCD_Edge_new(Naj_v, Element_q, R_a, ngp) result(rs)!
		real(kind = dp), intent(in) :: Naj_v(3), R_a(3)  
      real(kind = dp), dimension(8, 3), intent(in) :: Element_q
      integer, intent(in) :: ngp
      integer :: i, j, ngp_1, ngp_2     
      real(kind = dp), dimension(100) :: x1, w1, x2, w2
      complex(kind = dp) :: ff(5, 10), rs(3, 10)
      complex(kind = dp), dimension(10) :: KA, KBE, KBH, KCH, KCE  
      real(kind = dp) :: R(3), R_q(3) ! Nqi(10, 3),
      real(kind = dp) ::  mu_q, nu_q, TRF
		type(derived_element_parameters_typeB) :: delement_p
		
      do i = 1, 10 !
         KA(i) = (0.0, 0.0)        
         KBE(i) = (0.0, 0.0)    
         KBH(i) = (0.0, 0.0)
         KCE(i) = (0.0, 0.0)
         KCH(i) = (0.0, 0.0)      
      end do
      ngp_1 = ngp
      call xw_GPoint_1D(x1, w1, ngp_1)
      ngp_2 = ngp
      call xw_GPoint_1D(x2, w2, ngp_2)
      do i = 1, ngp_2
         mu_q = x1(i)
         do j = 1, ngp_2 !Integration over q-element
            nu_q = x2(j)
            call paras_st_B(mu_q, nu_q, Element_q, delement_p)	
				R_q = r_local_parametric(mu_q, nu_q, Element_q)
            R = R_a - R_q
            ff = ff_ABC_II_Edge(Naj_v, R, delement_p)
            TRF = w1(i)*w2(j)
            KA = KA + ff(1, :)*TRF
            KBE = KBE + ff(2, :)*TRF
            KBH = KBH + ff(3, :)*TRF
            KCE = KCE + ff(4, :)*TRF
            KCH = KCH + ff(5, :)*TRF
         end do ! Loop j 
      end do ! Loop i       
      rs(1, 1:10) = KA
      rs(2, 1:10) = KBE/(im*Omega*eps_0) - (-im*Omega*my_0*KCE) ! 
      rs(3, 1:10) = KBH/(-im*Omega*my_0)- im*Omega*eps_0*KCH
      return
   end function Integration_KABCD_Edge_new
    
   function Integration_KABCD_VA(m, n, s, theta_p, Element_a, theta_q, Element_q, pI_a, Psi, N1)result(rv) !Equation 74(a)
      real(kind = dp), intent(in) :: theta_p, theta_q, Psi, Element_q(8, 3), Element_a(8, 3)!
      integer, intent(in) :: m, n, s, N1, pI_a! 
      real(kind = dp) :: a, b, am, bm, x
      real(kind = dp), dimension(3) :: Ra, Rq!
      real(kind = dp), dimension(100) :: x1, w1
      complex(kind = dp) :: ff(5, 10), rv(3, 10)
      complex(kind = dp), dimension(10) :: KA, KBE, KBH, KCH, KCE 
        
      real(kind = dp) :: R(3)!
      real(kind = dp) ::  mu_q, nu_q, mu_a, nu_a, TRF
        
     type(derived_element_parameters_typeB) :: delement_pq 
	  type(derived_element_parameters_typeA) :: delement_pa 
		integer :: j   
        
		KA(1:10) = (0.0, 0.0)        
		KBE(1:10) = (0.0, 0.0)    
      KBH(1:10) = (0.0, 0.0)
      KCE(1:10) = (0.0, 0.0)
		KCH(1:10) = (0.0, 0.0)      
        
      a = 0.0
      if (s==1) then! case a
         b = Lp_Psi_VA(m, theta_p, Psi) !Lp 
      else 
         b = Lq_Psi_VA(n, theta_q, Psi)
      end if        
      am = (b-a)/2
      bm = (b+a)/2
      call xw_GPoint_1D(x1, w1, N1)        
		do j = 1, N1
			x = am*x1(j) + bm
			mu_a = -1 + x*cos(Psi)*cos(theta_p)
			nu_a = -1 + x*cos(Psi)*sin(theta_p)            
			mu_q = -1 + x*sin(Psi)*cos(theta_q)
			nu_q = -1 + x*sin(Psi)*sin(theta_q)  
			
			call paras_st_A(mu_a, nu_a, Element_a, delement_pa)	
			call paras_st_B(mu_q, nu_q, Element_q, delement_pq)	
			Rq = r_local_parametric(mu_q, nu_q, Element_q)
			Ra = r_local_parametric(mu_a, nu_a, Element_a)
			R = Ra - Rq 			
			ff = ff_ABC_II_Edge(delement_pa%Naj(pI_a)%point, R, delement_pq)
			TRF = am*w1(j)*cos(Psi)*sin(Psi)*delement_pa%JJ*x**3
			KA = KA + ff(1, :)*TRF!            
			KBE = KBE + ff(2, :)*TRF!!!!!!!
			KBH = KBH + ff(3, :)*TRF!2 !!!!!!
			KCE = KCE + ff(4, :)*TRF
			KCH = KCH + ff(5, :)*TRF       
		end do ! Loop j 
		rv(1, 1:10) = KA
		rv(2, 1:10) = KBE/(im*Omega*eps_0) - (-im*Omega*my_0*KCE )! 
		rv(3, 1:10) = KBH/(-im*Omega*my_0) - im*Omega*eps_0*KCH !      
		return
   end function Integration_KABCD_VA
	
	function ff_ABC_II_Edge(Naj_v, Raq, delement_p)result(rs)
		type(derived_element_parameters_typeB), intent(in) :: delement_p
      real(kind = dp), intent(in) :: Naj_v(3), Raq(3)
      complex(kind = dp) :: rs(5, 10), GG_SS_m, GG_SS_a, ff_1, ff_2, ff_3, Greenf_m, Greenf_a
      complex(kind = dp), dimension(10) :: FF_AS, FF_BES ! 
      complex(kind = dp), dimension(10) :: FF_BHS, FF_CE, FF_CH ! Strong and weak parts
      real(kind = dp) :: aa(3), rr_n(3), JJ, nR
      integer :: n
		  
      nR = vec_len(Raq)
      rr_n = Raq/vec_len(Raq) 
		JJ  = delement_p%JJ
		call Diver_Green(k1, nR, Greenf_a, GG_SS_a)         
      call Diver_Green(k2, nR, Greenf_m, GG_SS_m)	
      do n= 1, 10
			aa = cross_r(delement_p%nq%point, delement_p%Nqi(n)%point) 
			ff_1 = dot_product(Naj_v, cross_r(aa, rr_n))*JJ ! kaji, Eq.(3.63)
			ff_2 = dot_product(Naj_v, rr_n)*delement_p%n_rot_n(n)*JJ ! Eq.(3.64)
			ff_3 = dot_product(Naj_v, aa)*JJ ! kcji, Eq.(3.65) 
			
			! Strongly singular part for ka !ff_1!* + GG_SS_m
			FF_AS(n) = ff_1*(GG_SS_a + GG_SS_m) 
			
			! Strongly singular part for kb
			FF_BES(n) = ff_2*(1/eps_r1*GG_SS_a +1/eps_r2*GG_SS_m) 			
			FF_BHS(n) = ff_2*(GG_SS_a + GG_SS_m) !
			
			! Weakly singular part for kc !
			FF_CE(n) = ff_3*(Greenf_a + Greenf_m) 
			FF_CH(n) = ff_3*(eps_r1*Greenf_a + eps_r2*Greenf_m)!
      end do        
      rs(1, :) = FF_AS
      rs(2, :) = FF_BES
      rs(3, :) = FF_BHS
      rs(4, :) = FF_CE
      rs(5, :) = FF_CH
      return
   end function ff_ABC_II_Edge !
    
	!handle the singularity for self term	
	function ff_phi_AB_Edge(Naj_v, A_vec, derived_pa)result(rv)
		real(kind = dp), dimension(3), intent(in) :: A_vec!
		type(derived_element_parameters_typeB), intent(in) :: derived_pa
		  
		real(kind = dp), intent(in) :: Naj_v(3)!
		real(kind = dp) :: rv(2, 10)
		real(kind = dp) :: aa(3), A_amp, bb(3)
		integer :: n        
		A_amp = vec_len(A_vec)
		bb = -A_vec/A_amp**3/(4*PI) !
		do n = 1, 10
			aa = cross_r(derived_pa%nq%point, derived_pa%Nqi(n)%point)
			rv(1, n) = dot_product(Naj_v, cross_r(aa, bb))*derived_pa%JJ !! 
			rv(2, n) = dot_product(Naj_v, bb)*derived_pa%JJ*derived_pa%n_rot_n(n)!  
		end do
      return   
	end function ff_phi_AB_Edge
   
   !Integration between edge adjacent elements      
    function uv_EA(s, theta, Psi, u, Lambd) 
        real(kind = dp), dimension(2) :: uv_EA
        real(kind = dp) :: uu
        real(kind = dp), intent(in) :: Psi, theta, Lambd, u
        integer, intent(in) :: s
        uu = u
        if (s==1) then !p
            uv_EA(1) = uu !u
            uv_EA(2) = Lambd*sin(Psi)-1 ! v
        else if (s==2) then !q
            uv_EA(1) = Lambd*cos(Psi)*cos(theta)-uu !u
            uv_EA(2) = Lambd*cos(Psi)*sin(theta)-1 ! v
        end if
        return
    end function uv_EA
    
    function Lm_Psi_EA(ll, mm, nn, theta)
        !implicit none
        real (kind = dp), intent(in) :: theta
        integer, intent(in) :: ll, mm, nn
        
        real(kind= dp), dimension(2) :: Lm_Psi_EA
        real(kind = dp):: Psi_theta_1, Psi_theta_2, Psi_12
        
        Psi_theta_1 = atan(cos(theta))
        Psi_theta_2 = atan(sin(theta))       
        Psi_12 = atan(sin(theta))
        
        if ((ll == 0) .AND. (mm == 0)) then          
            if (nn == 1) then
                Lm_Psi_EA(1) = 0 
                Lm_Psi_EA(2) = Psi_theta_1
            else if (nn==2) then
                Lm_Psi_EA(1) = Psi_theta_1
                Lm_Psi_EA(2) = PI/2
            else if (nn==3) then
                Lm_Psi_EA(1) = 0
                Lm_Psi_EA(2) = Psi_theta_2
            else if (nn==4) then
                Lm_Psi_EA(1) = Psi_theta_2
                Lm_Psi_EA(2) = PI/2
            else 
                print*, 'n is too large'
            end if            
        else if ((ll == 0) .AND. (mm == 1)) then
            if (nn == 1) then
                Lm_Psi_EA(1) = Psi_theta_1
                Lm_Psi_EA(2) = PI/2
            else if (nn==2) then
                Lm_Psi_EA(1) = Psi_theta_2
                Lm_Psi_EA(2) = PI/2
            else 
                print*, 'nn is too large'
            end if            
        else if ((ll == 1) .AND. (mm == 0)) then
            if (nn == 1) then
                Lm_Psi_EA(1) = 0.0
                Lm_Psi_EA(2) = Psi_12
            else if (nn==2) then
                Lm_Psi_EA(1) = 0.0
                Lm_Psi_EA(2) = Psi_12
            else 
                print*, 'n is too large'
            end if 
            
        else if ((ll == 1) .AND. (mm == 1)) then     
            if (nn == 1) then
                Lm_Psi_EA(1) = Psi_12
                Lm_Psi_EA(2) = PI/2
            else if (nn==2) then
                Lm_Psi_EA(1) = Psi_12
                Lm_Psi_EA(2) = PI/2
            else 
                print*, 'n is too large'
            end if             
            
        else if ((ll == 2) .AND. (mm == 0)) then    
            if (nn == 1) then
                Lm_Psi_EA(1) = 0
                Lm_Psi_EA(2) = Psi_theta_2
            else if (nn==2) then
                Lm_Psi_EA(1) = Psi_theta_2
                Lm_Psi_EA(2) = PI/2
            else if (nn==3) then
                Lm_Psi_EA(1) = 0
                Lm_Psi_EA(2) = -Psi_theta_1
            else if (nn==4) then
                Lm_Psi_EA(1) = -Psi_theta_1
                Lm_Psi_EA(2) = PI/2
            else 
                print*, 'n is too large'
            end if             
        else if ((ll == 2) .AND. (mm == 1)) then      
            if (nn == 1) then
                Lm_Psi_EA(1) = Psi_theta_2
                Lm_Psi_EA(2) = PI/2
            else if (nn==2) then
                Lm_Psi_EA(1) = -Psi_theta_1
                Lm_Psi_EA(2) = PI/2
            else 
                print*, 'n is too large'
            end if 
        else
            print*, 'll is too large'
        end if
        return
    end function Lm_Psi_EA
    
    function Lm_mu_EA(ll, mm, nn, theta, Psi)
        !implicit none
        real (kind = dp), intent(in) :: theta, Psi
        integer, intent(in) :: ll, mm, nn
        real(kind= dp), dimension(2) :: Lm_mu_EA, Mu_ab
        real(kind = dp) :: u1_Psi, u2_Psi, u1_theta, u2_theta
        
        u1_Psi = 2*cos(theta)*tan(PI/2-Psi)-1
        u2_Psi = 2*cos(theta)*tan(PI/2-Psi)+1
        u1_theta = 2*tan(PI/2-theta)-1
        u2_theta = 2*tan(PI/2-theta)+1 
        
        if ((ll == 0) .AND. (mm == 0)) then          
            if (nn == 1) then
                Mu_ab(1) = -1
                Mu_ab(2) = 1
            else if (nn==2) then
                Mu_ab(1) = -1
                Mu_ab(2) = u1_Psi
            else if (nn==3) then
                Mu_ab(1) = -1
                Mu_ab(2) = u1_theta
            else if (nn==4) then
                Mu_ab(1) = -1
                Mu_ab(2) = u1_Psi
            else 
                print*, 'n is too large'
            end if            
        else if ((ll == 0) .AND. (mm == 1)) then
            if (nn == 1) then
                Mu_ab(1) = u1_Psi
                Mu_ab(2) = 1
            else if (nn==2) then
                Mu_ab(1) = u1_Psi
                Mu_ab(2) = u1_theta
            else 
                print*, 'nn is too large'
            end if            
        else if ((ll == 1) .AND. (mm == 0)) then
            if (nn == 1) then
                Mu_ab(1) = u1_theta
                Mu_ab(2) = 1
            else if (nn==2) then
                Mu_ab(1) = -1
                Mu_ab(2) = u2_theta
            else 
                print*, 'n is too large'
            end if 
            
        else if ((ll == 1) .AND. (mm == 1)) then     
            if (nn == 1) then
                Mu_ab(1) = u1_theta
                Mu_ab(2) = 1
            else if (nn==2) then
                Mu_ab(1) = -1
                Mu_ab(2) = u2_theta
            else 
                print*, 'n is too large'
            end if             
            
        else if ((ll == 2) .AND. (mm == 0)) then    
            if (nn == 1) then
                Mu_ab(1) = u2_theta
                Mu_ab(2) = 1
            else if (nn==2) then
                Mu_ab(1) = u2_Psi
                Mu_ab(2) = 1
            else if (nn==3) then
                Mu_ab(1) = -1
                Mu_ab(2) = 1
            else if (nn==4) then
                Mu_ab(1) = u2_Psi
                Mu_ab(2) = 1
            else 
                print*, 'n is too large'
            end if             
        else if ((ll == 2) .AND. (mm == 1)) then      
            if (nn == 1) then
                Mu_ab(1) = u2_theta
                Mu_ab(2) = u2_Psi
            else if (nn==2) then
                Mu_ab(1) = -1
                Mu_ab(2) = u2_Psi
            else 
                print*, 'n is too large'
            end if 
        else
            print*, 'll is too large'
        end if
        Lm_mu_EA = Mu_ab
        return
    end function Lm_mu_EA
    
	function Lm_Lambd_EA(ll, mm, theta, Psi, u)
		!implicit none
		real (kind = dp), intent(in) :: theta, Psi, u
		real (kind = dp) :: Lm_Lambd_EA
		integer, intent(in) :: ll, mm
        
		if ((ll == 0) .AND. (mm == 0)) then
			Lm_Lambd_EA = (1 + u)/(cos(theta)*cos(Psi))              
		else if ((ll == 0) .AND. (mm == 1)) then
			Lm_Lambd_EA = 2/sin(Psi)
		else if ((ll == 1) .AND. (mm == 0)) then
			Lm_Lambd_EA = 2/(sin(theta)*cos(Psi)) !treat it later
		else if ((ll == 1) .AND. (mm == 1)) then
			Lm_Lambd_EA = 2/sin(Psi)
		else if ((ll == 2) .AND. (mm == 0)) then
			Lm_Lambd_EA = (u - 1)/(cos(theta)*cos(Psi))
		else if ((ll == 2) .AND. (mm == 1)) then 
			Lm_Lambd_EA = 2/sin(Psi)
		end if 
		return
	end function Lm_Lambd_EA
    
	function Lm_theta_EA(ll, mm, nn)
		!implicit none
		real (kind = dp) :: Lm_theta_EA
		integer, intent(in) :: ll, mm, nn       
		if ((ll == 0) .AND. (mm == 0)) then
			if (nn == 1) then
					Lm_theta_EA = 0
			else if (nn==2) then
					Lm_theta_EA = 0
			else if (nn==3) then
					Lm_theta_EA = PI/4
			else if (nn==4) then
					Lm_theta_EA = PI/4
			else
					print*, 'nn is too large'
			end if
		else if ((ll == 0) .AND. (mm == 1)) then
			if (nn==1) then
					Lm_theta_EA = 0
			else if (nn==2) then
					Lm_theta_EA = PI/4
			else 
					print*, 'nn is too large'
			end if          
            
		else if ((ll == 1) .AND. (mm == 0)) then
			if (nn == 1) then
					Lm_theta_EA = PI/4
			else if (nn==2) then
					Lm_theta_EA = PI/2
			else
					print*, 'nn is too large'
			end if  
		else if ((ll == 1) .AND. (mm == 1)) then            
			if (nn==1) then
					Lm_theta_EA = PI/4
			else if (nn==2) then
					Lm_theta_EA = PI/2
			else
					print*, 'nn is too large'        
			end if
            
		else if ((ll == 2) .AND. (mm == 0)) then
			if ((nn==1) .or. (nn==2)) then
					Lm_theta_EA = PI/2
			else if ((nn==3).or. (nn==4)) then
					Lm_theta_EA = 3*PI/4
			else
					print*, 'nn is too large'
			end if
		else if ((ll == 2) .AND. (mm == 1)) then 
			if (nn==1) then
					Lm_theta_EA = PI/2
			else if (nn==2) then               
					Lm_theta_EA = 3*PI/4
			else 
					print*, 'n is too large'
			end if
		end if 
		return
	end function Lm_theta_EA

	function GL_Int_Lambd_EA_8Node(ll, mm, Element_a, Element_q, pI_a, theta, Psi, u, N1)result(rv) 
		real(kind = dp), intent(in) :: theta, Psi, u, Element_a(8, 3), Element_q(8, 3)
		integer, intent(in) :: ll, mm, N1, pI_a
		  
		real(kind = dp), dimension(3) :: Ra, Rq, Naj_v
		real(kind = dp), dimension(100) :: x, w1, x1 !
		real(kind = dp) :: a, b, am, bm!
        
		real(kind = dp) :: R(3)!
		real(kind = dp) :: mu_q, nu_q, mu_a, nu_a, TRF        
		type(derived_element_parameters_TypeB) :: delement_pq
		type(derived_element_parameters_TypeA) :: delement_pa
		  
		complex(kind = dp) :: rv(3, 10), Summ(5, 10), ff(5, 10)
        
		integer :: j 

		mu_a = u
      Summ(1:5, 1:10) = 0.0d0
      a = 0.0
      b = Lm_Lambd_EA(ll, mm, theta, Psi, mu_a)
      am = (b-a)/2
      bm = (b+a)/2        
      call xw_GPoint_1D(x1, w1, N1)
       
      x = am*x1 + bm
      do j = 1, N1 
			nu_a = -1 + x(j)*sin(Psi)            
			mu_q = x(j)*cos(Psi)*cos(theta)-mu_a
			nu_q = x(j)*cos(Psi)*sin(theta) - 1
            
			call paras_st_A(mu_a, nu_a, Element_a, delement_pa) 
			call paras_st_B(mu_q, nu_q, Element_q, delement_pq) 				
			Ra = r_local_parametric(mu_a, nu_a, Element_a)
			Rq = r_local_parametric(mu_q, nu_q, Element_q)
			R = Ra - Rq !
			Naj_v = delement_pa%Naj(pI_a)%point
			ff = ff_ABC_II_Edge(Naj_v, R, delement_pq)            
			TRF = am*w1(j)*delement_pa%JJ*x(j)**2*cos(Psi)
			Summ = Summ + ff*TRF            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do  
      rv(1, 1:10) = Summ(1, 1:10)
      rv(2, 1:10) = Summ(2, 1:10)/(im*Omega*eps_0) - (-im*Omega*my_0*Summ(4, 1:10))
      rv(3, 1:10) = Summ(3, 1:10)/(-im*Omega*my_0) - (im*Omega*eps_0*Summ(5, 1:10)) !        
      return
   end function GL_Int_Lambd_EA_8Node   

   function GL_Int_mu_EA_8Node(ll, mm, nn, Element_a, Element_q, pI_a, theta, Psi, N2, N1) !Equation 74(a)
      real(kind = dp), intent(in) :: theta, Psi, Element_a(8, 3), Element_q(8, 3)
      integer, intent(in) :: ll, mm, nn, N1, N2, pI_a! The first integration layer !s for (p, q)
      real(kind = dp) :: am, bm
      real(kind = dp), dimension(2) :: mu_ab ! Integration limit of mu
      real(kind = dp), dimension(100) :: x, w2, x2 !, Lambda
      complex(kind = dp) :: GL_Int_mu_EA_8Node(3, 10), Summ(3, 10)
      integer :: j ! The first integration layer !s for (p, q)        
        
		Summ(1:3, 1:10) = 0.0d0
      mu_ab = Lm_mu_EA(ll, mm, nn, theta, Psi)
      am = (mu_ab(2) - mu_ab(1))/2
      bm = (mu_ab(2) + mu_ab(1))/2
      call xw_GPoint_1D(x2, w2, N2)        
      x = am*x2 + bm
      do j = 1, N2
         Summ = Summ + &
				am*w2(j)*GL_Int_Lambd_EA_8Node(ll, mm, Element_a, Element_q, pI_a, theta, Psi, x(j), N1)
      end do       
      GL_Int_mu_EA_8Node = Summ
      return
   end function GL_Int_mu_EA_8Node
    
	function GL_Int_Psi_EA_8Node(ll, mm, nn, Element_a, Element_q, pI_a, theta, N3, N2, N1) !Equation 74(a)
		real(kind = dp), intent(in) :: theta, Element_a(8, 3), Element_q(8, 3)
		integer, intent(in) :: ll, mm, nn, pI_a
		real(kind = dp) :: am, bm
		real(kind = dp), dimension(2) :: psi_ab ! Integration limit of mu
		real(kind = dp), dimension(100) :: x, w3, x3 !, Lambda
		complex(kind = dp) :: GL_Int_Psi_EA_8Node(3, 10), Summ(3, 10)
		integer :: j, N1, N2, N3 ! The first integration layer !s for (p, q)                
     
		Summ(1:3, 1:10) = 0.0d0     
      psi_ab = Lm_Psi_EA(ll, mm, nn, theta)        
      am = (psi_ab(2) - psi_ab(1))/2
      bm = (psi_ab(2) + psi_ab(1))/2 
        
      call xw_GPoint_1D(x3, w3, N3)
      x = am*x3 + bm !x =Psi
      do j = 1, N3              
         Summ = Summ + &
			 am*w3(j)*GL_Int_mu_EA_8Node(ll, mm, nn, Element_a, Element_q, pI_a, theta, x(j), N2, N1)
      end do
      GL_Int_Psi_EA_8Node = Summ
      return
	end function GL_Int_Psi_EA_8Node
    ! 
	function GL_Int_theta_EA_8Node(ll, mm, nn, Element_a, Element_q, pI_a, N4, N3, N2, N1) !Equation 74(a)
		real(kind = dp), intent(in) :: Element_a(8, 3), Element_q(8, 3) 
		real(kind = dp) :: am, bm, a, b        
		real(kind = dp), dimension(100) :: x, w4, x4 !, Lambda
		complex(kind = dp) :: GL_Int_theta_EA_8Node(3, 10), Summ(3, 10)
		integer, intent(in) :: mm, ll, nn, N4, N1, N2, N3, pI_a! The first integration layer !s for (p, q)        
		integer :: i, j
        
		do j= 1, 10
			do i = 1, 3
			Summ(i, j) = 0.0d0
			end do
		end do          
		a = Lm_theta_EA(ll, mm, nn)
		b = a + PI/4
		am = (b-a)/2
		bm = (b+a)/2
		call xw_GPoint_1D(x4, w4, N4)            
		x = am*x4 + bm
		do j = 1, N4             
			Summ = Summ + &
				am*w4(j)*GL_Int_Psi_EA_8Node(ll, mm, nn, Element_a, Element_q, pI_a, x(j), N3, N2, N1)
		end do       
		GL_Int_theta_EA_8Node = Summ
		return
	end function GL_Int_theta_EA_8Node     
    
	function Lpq_VA(m, theta) ! Equation (70), 
		real(kind = dp), intent(in) :: theta
		integer, intent(in) :: m 
		real(kind = dp) :: Lpq_VA

		if (m ==1) then
				Lpq_VA = 2/cos(theta) ! Lp1, Lq1
		else if (m==2) then
				Lpq_VA = 2/sin(theta) ! Lq2, Lp2, 
		else 
				print*, 'Value of m overflow'
		end if
		return     
	end function Lpq_VA
     
	function Lp_Psi_VA(m, theta_p, Psi) ! Euqation 74, paper Tambova, .. Polimeridis (2017)                             
		real(kind = dp), intent(in) :: Psi, theta_p ! The second in (74)
		integer, intent(in) :: m        
		real(kind = dp) :: Lp_Psi_VA 
		Lp_Psi_VA = Lpq_VA(m, theta_p)/cos(Psi) !L1(m,n) in (74)
		return     
	end function Lp_Psi_VA

	function Lq_Psi_VA(n, theta_q, Psi) ! Euqation 74, paper Tambova, .. Polimeridis (2017)                                
		real(kind = dp), intent(in) :: theta_q, Psi
		integer, intent(in) :: n        
		real(kind = dp) :: Lq_Psi_VA
        
		Lq_Psi_VA = Lpq_VA(n, theta_q)/sin(Psi) !L1(m,n) in (74)
		return     
	end function Lq_Psi_VA   

	function PsiI_VA(m, n, theta_p, theta_q) !The third equation in (74) !
		real(kind = dp), intent(in) :: theta_p, theta_q
		integer, intent(in) :: m, n
		real(kind = dp) :: PsiI_VA, Lp, Lq 
        
		Lp = Lpq_VA(m, theta_p) 
		Lq = Lpq_VA(n, theta_q) 
		PsiI_VA = atan(Lq/Lp)
		return
	end function PsiI_VA

	!Integration between vertex overlapped elements
	function GL_Int_Psi_VA_8Node(m, n, s, theta_p, Element_p, theta_q, Element_q, pI_a, N2, N1) !Equation 74(a),  theta_p, theta_q, 
		real(kind = dp), intent(in) :: theta_p, theta_q, Element_p(8, 3), Element_q(8, 3)
		integer, intent(in) :: m, n, s, N2, N1, pI_a  ! The second integration loop
		real(kind = dp) :: a, b, am, bm, PsiI, TRF
		real(kind = dp), dimension(100) ::  x, w, x1 !x is Lambda         
		complex(kind = dp) :: GL_Int_Psi_VA_8Node(3, 10), Summ(3, 10) 
		integer :: j, i ! The second integration loop
       
		do j= 1, 10
			do i = 1, 3
			Summ(i, j) = 0.0d0
			end do
		end do        
        
		PsiI = PsiI_VA(m, n, theta_p, theta_q)
		if (s==1) then!case p
			a = 0.0d0
			b = PsiI 
		else
			a = PsiI
			b = PI/2
		end if        
		am = (b-a)/2
		bm = (b+a)/2        
		call xw_GPoint_1D(x, w, N2)        
       
		x1 = am*x + bm
		do j = 1, N2
			TRF = am*w(j)
			Summ = Summ + Integration_KABCD_VA(m, n, s, theta_p, Element_p, theta_q, Element_q, pI_a, x1(j), N1)*TRF
		end do        
		GL_Int_Psi_VA_8Node = Summ         
		return
	end function GL_Int_Psi_VA_8Node    
    
	function GL_Int_ThetaQ_VA_8Node(m, n, s, theta_a, Element_a, Element_q, pI_a, N3, N2, N1) !GL: Gaussian laguerre ! The third loop, N4
		real(kind = dp), intent(in) :: theta_a, Element_a(8, 3), Element_q(8, 3) 
		integer, intent(in) :: N3, N2, N1, m, n, s, pI_a! The third integration loop
		real(kind = dp) :: am, bm
		real(kind = dp), dimension(100) :: x, w, x1
		complex(kind = dp) :: GL_Int_ThetaQ_VA_8Node(3, 10), Summ(3, 10)
		real(kind = 8), dimension(2) :: aa
		integer :: i, j ! The third integration loop
		do j= 1, 10
			do i = 1, 3
					Summ(i, j) = 0.0d0
			end do
		end do         
		aa = Lm_Theta(n)        
		am = (aa(2) - aa(1))/2
		bm = (aa(2) + aa(1))/2        
		call xw_GPoint_1D(x, w, N3)
		x1 = am*x + bm        
		!Theta_a is theta_p in the formula, x1 is theta_q
		do j = 1, N3
			Summ = Summ + GL_Int_Psi_VA_8Node(m, n, s, theta_a, Element_a, x1(j), Element_q, pI_a, N2, N1)*am*w(j)
		end do         
		GL_Int_ThetaQ_VA_8Node = Summ        
		return
	end function GL_Int_ThetaQ_VA_8Node    !
    !
	function GL_Int_ThetaP_VA_8Node(m, n, s, Element_p, Element_q, pI_a, N4, N3, N2, N1) !GL: Gaussian laguerre
		integer, intent(in) :: m, n, s, N4, N3, N2, N1, pI_a! The fourth integration loop
		real(kind = dp), intent(in) :: Element_p(8, 3), Element_q(8, 3)
		real(kind = dp) :: am, bm
		real(kind = dp), dimension(100) :: x, w, x1 !x is 
		real(kind = dp), dimension(2) :: aa
		complex(kind = dp) :: GL_Int_ThetaP_VA_8Node(3, 10), Summ(3, 10)
		integer :: i, j ! The fourth integration loop
		do j= 1, 10
			do i = 1, 3
					Summ(i, j) = 0.0d0
			end do
		end do        
		aa = Lm_Theta(m)
		am = (aa(2) - aa(1))/2
		bm = (aa(2) + aa(1))/2        
		call xw_GPoint_1D(x, w, N4)
		x1 = am*x + bm
		do j = 1, N4
			Summ = Summ + am*w(j)*GL_Int_ThetaQ_VA_8Node(m, n, s, x1(j), Element_p, Element_q, pI_a, N3, N2, N1)
		end do      
		GL_Int_ThetaP_VA_8Node = Summ   
		!print*, 'ThetaP_VA', Summ(1, 1)
		return
	end function GL_Int_ThetaP_VA_8Node    
  
	function Lm_Theta(m)
		integer, intent(in) :: m
		real(kind = dp), dimension(2) :: Lm_Theta
        
		if (m ==1) then
			Lm_Theta(1) = 0
			Lm_Theta(2) = PI/4
			else   
			Lm_Theta(1) = PI/4
			Lm_Theta(2) = PI/2
		end if
		return
	end function Lm_Theta  
     
end module lib_sie_quad_singularity
    