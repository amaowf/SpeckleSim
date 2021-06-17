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

module lib_sie_mlfmm
	use lib_sie_math
	use lib_sie_data_container
	use lib_sie_lookup_table
	use lib_sie_tri_calculation_mod
	
	implicit none
	
	type(csr_format) :: field_csr
	
	contains
	
	function LT3_get_TL_at_loc(r_a, r_b, l)result(TL)
		
		real(dp), intent(in) :: r_a(3), r_b(3)
		real(dp) :: r_ab(3), col_value
		integer, intent(in) :: l
		type(lib_ml_fmm_coefficient) :: TL
		type(LT3_loc_t)	:: my_loc		
		
		r_ab = r_a - r_b	
		my_loc = LT3_get_loc(preCalculated_TL%my_lt(l), r_ab(1), r_ab(2), r_ab(3))		
		col_value = LT3_get_col_at_loc(preCalculated_TL%my_lt(l), 1, my_loc)
		TL%a_nm = preCalculated_TL%TL(l)%item(NINT(col_value))		
	end function LT3_get_TL_at_loc
	
	subroutine precalculation_LT3_TL(k_neighbour, l_max)	
		implicit none
		
		integer(kind = 1), intent(in) :: l_max
		real(dp) :: x_min, x_min_3d(3)  !< Minimum x-coordinate
		real(dp) :: x_max, x_max_3d(3)  !< Maximum x-coordinate
		real(dp) :: d_c
		real(dp), dimension(:, :, :, :), allocatable :: spaced_data
		
		integer(kind=UINDEX_BYTES), intent(in) :: k_neighbour!
		integer :: n_cols
		integer :: i, n_points, n_points_3d(3)
		integer(kind = 1) :: l_min,  l
		
		n_cols = 1		
		! For the case k_neighbour = 1				
		if (k_neighbour == 1 ) then
			n_points = 7			
		else if (k_neighbour == 2) then
			n_points = 11			
		else
			print*, 'Not a proper neighbour input in the LT3 creator'
		end if
		
		allocate(preCalculated_TL%TL(l_max))
		
		if (k_neighbour == 1)then	
			l_min = 2
			allocate(preCalculated_TL%TL(1)%item(1))
		else if (k_neighbour == 2)then
			l_min = 2
			allocate(preCalculated_TL%TL(1)%item(1))
			allocate(preCalculated_TL%TL(2)%item(1))
		end if
		
		allocate(preCalculated_TL%my_lt(l_max))
		do l = l_min, l_max				
			d_c =  lib_tree_get_box_edge_length(l)*lib_tree_scaling_D%x(1)			
			x_min = -(n_points-1)/2*d_c
			x_max = -x_min			
			x_min_3d = (/x_min, x_min, x_min/)
			x_max_3d = (/x_max, x_max, x_max/)
			n_points_3d = (/n_points, n_points, n_points/)
			spaced_data = LT3_create_spaced_data(n_points_3d, n_cols)
			preCalculated_TL%my_LT(l) = LT3_create_from_data(x_min_3d, x_max_3d, spaced_data)		
			call PreCalc_transfer_function(preCalculated_TL%my_LT(l), spaced_data, l)
		end do	
		
		
	!#ifdef _TRACE_
	
        print *, "TRACE: precalculation_LT3_TL ..done"       
		  
	!#endif
		
	end subroutine precalculation_LT3_TL
	
	subroutine PreCalc_transfer_function(my_LT, spaced_data, l)
		
		integer(kind = 1), intent(in) :: l

		real(dp), dimension(:, :, :, :), intent(in) :: spaced_data		
		real(dp), dimension(:), allocatable :: x_data, y_data, z_data
		real(dp) :: r_ab(3)
		
		integer :: i, j, k, nx, ny, nz, t, K_m
		type(LT3_t), intent(in) :: my_LT		
		type(lib_ml_fmm_coefficient) :: TL_tmp
		
		nx = size(spaced_data, 1)
		ny = size(spaced_data, 2)
		nz = size(spaced_data, 3)
		
		K_m = truncation_number(l)%n		
		allocate(x_data(nx), y_data(ny), z_data(nz))	
		x_data = LT_get_xdata(my_LT%x_min(1), my_LT%dx(1), my_LT%n_points(1))
		y_data = LT_get_xdata(my_LT%x_min(2), my_LT%dx(2), my_LT%n_points(2))
		z_data = LT_get_xdata(my_LT%x_min(3), my_LT%dx(3), my_LT%n_points(3))		
		
		if (allocated(preCalculated_TL%TL(l)%item)) then
			deallocate(preCalculated_TL%TL(l)%item)
		end if		
		allocate(preCalculated_TL%TL(l)%item(nx*ny*nz))
		do i = 1, nx
			do j = 1, ny
				do k = 1, nz
					t = k + (j-1)*ny + (i-1)*ny*nx
					r_ab = (/x_data(i), y_data(j), z_data(k)/)				
					if (spaced_data(i, j, k, 1) .gt. 0.0_dp) then
						TL_tmp = TL_km_ml_fmm(r_ab, K_m, l)						
						preCalculated_TL%TL(l)%item(t) = TL_tmp%a_nm
					else
						call init_list_sie(TL_tmp, 2, 1)	
						preCalculated_TL%TL(l)%item(t) = TL_tmp%a_nm
					end if
				end do
			end do
		end do
	end subroutine PreCalc_transfer_function
	
	!r_ab: distance of the two boxes
	!k_m: sampling rate along phi
	function TL_km_ml_fmm(r_ab, K_m, level)result(TL_tilde)
      integer, intent(in) :: K_m    !
		integer(kind = 1), intent(in) :: level
      real(dp), intent(in) :: r_ab(3)
      real(dp) :: r_hat(3), x
		real(dp), dimension(:, :), allocatable :: pm_arr
      complex(dp) :: r_h(2), kd(2)
		
		!dummy
      real (dp), dimension(K_m) :: pm, dummy 
      complex(dp), dimension(:, :), allocatable :: TL
		type(lib_ml_fmm_coefficient) :: TL_tilde
      complex(dp), dimension(K_m) :: hl 
      integer :: j, fnu, m, i, Ns
		integer(kind = 1) :: l_max
		m = K_m
      fnu = 0        
      r_hat = r_ab/vec_len(r_ab)
		kd = (/k1, k2/)
      r_h =  kd*vec_len(r_ab) 
		Ns = 2*K_m*K_m
		if (allocated(TL)) then
			deallocate(TL)
		end if
		
		call init_list_sie(TL_tilde, 2, Ns+2)
		l_max = lib_tree_get_level_max(tree_s_opt)
		
		allocate(pm_arr(Ns+2, m), TL(2, Ns+2)) !! different in TL_km_arr
		TL(1:2, 1:Ns+2) = (0.0, 0.0)   
		call Khat_and_TRF_ExtendedK(K_m)!
      do j = 1, Ns+2
         x = dot_product(r_hat(:), k_hat_p(j, :))
         call lib_math_legendre_polynomial(x, fnu, m, pm, dummy)
         pm_arr(j, :) = pm(:)
      end do
		
		!When k2 is complex, only at the level lmax 
		!k2 will be considered.
		if ((abs(imag(k2)) .gt. 0.0) .and. (level .lt. l_max-1))then
			i = 1
			hl =  lib_math_hankel_spherical_2(r_h(1), fnu, m)
			do j = 1, m			
				TL(1, :) = TL(1, :) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd(i)/(4*PI)**2 ! 
			end do			
			TL_tilde%a_nm%item(i)%item(:) = TL(1, :)
		else
			do i = 1, 2
				hl =  lib_math_hankel_spherical_2(r_h(i), fnu, m) 			
				do j = 1, m			
					TL(i, :) = TL(i, :) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd(i)/(4*PI)**2 ! 
				end do
			end do		
			do i = 1, 2
				TL_tilde%a_nm%item(i)%item(:) = TL(i, :)
			end do
		end if
      return 
		deallocate(pm_arr, TL)
	end function TL_km_ml_fmm
	
	subroutine Khat_and_TRF_ExtendedK(K_m)! 	
		implicit none	
		
		real(dp), dimension(:), allocatable :: xx, ww!, x        
		real(dp) :: d_phi, phi, a, b, dzero
        
		integer :: j, k, i, N_sp
		integer, intent(in) :: K_m
		dzero = 0.0		
		N_sp = 2*K_m*K_m
        
        if (allocated(k_hat_p)) then
                deallocate(k_hat_p)
        end if
        
        if (allocated(TRF_p)) then
                deallocate(TRF_p)
        end if
        
        if (allocated(transfer_theta)) then
                deallocate(transfer_theta)
        end if
        
        if (allocated(transfer_phi)) then
                deallocate(transfer_phi)
        end if
        
		allocate(k_hat_p(N_sp+2, 3))
		allocate(TRF_p(N_sp+2))
		allocate(transfer_theta(N_sp+2, 3), transfer_phi(N_sp+2, 3))
		allocate(xx(K_m))
		allocate(ww(K_m))
		a = -1.0
		b = 1.0		
		call legendre_handle(K_m, a, b, xx, ww)
		
		d_phi = 2*PI/(2*K_m)
		
		do j = 1, K_m		   
			do k = 1, K_m*2
				phi = (k-1)*d_phi 
				i = k + (j-1)*K_m*2				
				k_hat_p(i, :) = (/sin(acos(xx(j)))*cos(phi), sin(acos(xx(j)))*sin(phi), xx(j)/) !,
				transfer_theta(i, 1:3) = (/xx(j)*cos(phi), xx(j)*sin(phi), -sin(acos(xx(j)))/)
				transfer_phi(i, 1:3) = (/-sin(phi), cos(phi), dzero/)
				TRF_p(i) = ww(j)*d_phi				
			end do
			
		end do
		k_hat_p(N_sp+1, :)= (/0.0, 0.0, 1.0/)
		k_hat_p(N_sp+2, :)= (/0.0, 0.0, -1.0/)
		transfer_theta(N_sp+1, :) = (/1.0, 0.0, 0.0/) !For saving fx at theta = 0
		transfer_theta(N_sp+2, :) = (/1.0, 0.0, 0.0/) !For saving fx at theta = 0		
		transfer_phi(N_sp+1, :) = (/0.0, 1.0, 0.0/)
		transfer_phi(N_sp+2, :) = (/0.0, 1.0, 0.0/)
		TRF_p(N_sp+1:N_sp+2) = 0.0
		deallocate(xx)
		deallocate(ww)
   end subroutine Khat_and_TRF_ExtendedK
	
	subroutine fx_fy(fn_scs)		
		use libmath
		implicit none
		type(lib_ml_fmm_coefficient), intent(in) :: fn_scs 	
		integer :: i, N
		call init_list_sie(fxy, 2, 2)
		N = size(Fn_scs%a_nm%item(1)%item)
		do i = 1, 2 
			fxy%a_nm%item(i)%item(1) = Fn_scs%a_nm%item(i)%item(N-1) !fx, Ns+2: theta = 0
			fxy%b_nm%item(i)%item(1) = Fn_scs%b_nm%item(i)%item(N-1) !fy, Ns+2: theta = 0
			!new
			fxy%a_nm%item(i)%item(2) = Fn_scs%a_nm%item(i)%item(N) !fx, Ns+1: theta = pi
			fxy%b_nm%item(i)%item(2) = Fn_scs%b_nm%item(i)%item(N) !fy, Ns+1: theta = pi			
		end do
	end subroutine
	
	subroutine Lagrange_Interpolation_New(pp, K_m, K_n, Fn_scs, f_inter)
		implicit none
		
		integer, intent(in) :: pp, K_m, K_n
		type(lib_ml_fmm_coefficient), intent(in) :: Fn_scs
		type(lib_ml_fmm_coefficient), intent(out) :: f_inter
		
		integer :: i, j, k, tt, Ns_m, Ns_n
		integer :: it, jp, m_phi, n_phi, t, s, kp, kt
				
		real(dp) :: dphi, a, b, theta, phi		
		complex(dp), dimension(:, :), allocatable :: ff_inter_a, ff_inter_b!		
		type(lib_ml_fmm_coefficient) :: C_tmp
		
		complex(dp), dimension(:, :, :), allocatable :: ff_in, ff_out_a, ff_out_b!
		real(dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		!	 
		n_phi = 2*K_n
		m_phi = 2*K_m
		dphi = 2*PI/m_phi		
		
		Ns_n = K_n*K_n*2
		Ns_m = K_m*K_m*2
		
		allocate(ff_in(K_m, 2*K_m, 4))!,  ff_out_a(K_m+2, 2*K_m, 4))!
		allocate(ff_inter_a(n_phi*K_n, 4)) !
		allocate(ff_inter_b(n_phi*K_n, 4)) !
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
		call init_list_sie(f_inter, 4, Ns_n+2)
		call init_list_sie(C_tmp, 2, 2)
				
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				! for xe
				ff_in(i, j, 1)= fn_scs%a_nm%item(1)%item(k) ! ff_in(i, j, 1:2) theta, k1 and k2
				ff_in(i, j, 2)= fn_scs%a_nm%item(2)%item(k) !
				ff_in(i, j, 3)= fn_scs%b_nm%item(1)%item(k) ! ff_in(i, j, 3:4) phi, k1 and k2
				ff_in(i, j, 4)= fn_scs%b_nm%item(2)%item(k) !
			end do
		end do
		
		C_tmp%a_nm%item(1:2) = fn_scs%a_nm%item(1:2)
		C_tmp%b_nm%item(1:2) = fn_scs%b_nm%item(1:2)
		call fx_fy(C_tmp)
		call Sorting_theta_phi_LenReduced_Arr_New(pp, K_m, ff_in, ff_out_a)				
		do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j		
				!for xh
				ff_in(i, j, 1)= fn_scs%a_nm%item(3)%item(k) ! ff_in(i, j, 1:2) theta, k1 and k2 
				ff_in(i, j, 2)= fn_scs%a_nm%item(4)%item(k) !
				ff_in(i, j, 3)= fn_scs%b_nm%item(3)%item(k) ! ff_in(i, j, 3:4) phi, k1 and k2
				ff_in(i, j, 4)= fn_scs%b_nm%item(4)%item(k) !
			end do
		end do
		C_tmp%a_nm%item(1:2) = fn_scs%a_nm%item(3:4)
		C_tmp%b_nm%item(1:2) = fn_scs%b_nm%item(3:4)
		call fx_fy(C_tmp)
		call Sorting_theta_phi_LenReduced_Arr_New(pp, K_m, ff_in, ff_out_b)
		deallocate(ff_in)!
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		a = -1.0 
		b = 1.0
		call legendre_handle(K_n, a, b, xt, x_tmp)	
		call legendre_handle(K_m, a, b, xt_old, x_tmp_old)	
		do i = 1, n_phi
			xp(i) = (i-1)*6*PI/n_phi
		end do		
		do i = 1, m_phi
			xp_old(i) = (i-1)*2*PI/m_phi		
		end do		
		!
		!it: index for the interpolated theta		
		do it = 1, K_n
			theta = acos(xt(it))			
			t = floor((it-1.0)/K_n*K_m) + 1
			kt = pp + t
			call Lagrange_Interpolation_Weights(theta, t, pp, xt_old, a_theta, wk_theta)
			do jp = 1, n_phi			   
				tt = (it-1)*n_phi + jp			
				ff_inter_a(tt, 1:4) = (0.0, 0.0)
				ff_inter_b(tt, 1:4) = (0.0, 0.0)
				phi = (jp-1)*2*PI/n_phi
				s = floor((jp-1.0)/K_n*K_m)+1	
				kp = pp + s
				call Lagrange_Interpolation_Weights(phi, s, pp, xp_old, a_phi, wk_phi)
				do i = 1, 2*pp
					do j = 1, 2*pp						
						ff_inter_a(tt, 1:4) = ff_inter_a(tt, 1:4) + wk_phi(j)*wk_theta(i)*ff_out_a(kt+i-pp, kp+j-pp, 1:4)
						ff_inter_b(tt, 1:4) = ff_inter_b(tt, 1:4) + wk_phi(j)*wk_theta(i)*ff_out_b(kt+i-pp, kp+j-pp, 1:4)						
					end do
				end do				
			end do
		end do
		! for x_e
		f_inter%a_nm%item(1)%item(1:Ns_n)= ff_inter_a(1:Ns_n, 1)
		f_inter%a_nm%item(2)%item(1:Ns_n)= ff_inter_a(1:Ns_n, 2) 
		f_inter%b_nm%item(1)%item(1:Ns_n)= ff_inter_a(1:Ns_n, 3) 
		f_inter%b_nm%item(2)%item(1:Ns_n)= ff_inter_a(1:Ns_n, 4) 
		! for x_h	
		f_inter%a_nm%item(3)%item(1:Ns_n)= ff_inter_b(1:Ns_n, 1)
		f_inter%a_nm%item(4)%item(1:Ns_n)= ff_inter_b(1:Ns_n, 2) 
		f_inter%b_nm%item(3)%item(1:Ns_n)= ff_inter_b(1:Ns_n, 3) 
		f_inter%b_nm%item(4)%item(1:Ns_n)= ff_inter_b(1:Ns_n, 4) 
		!
		deallocate(ff_inter_a, ff_inter_b)!
		!save fx and fy at theta=0 and pi
		do i = 1, 2		
			f_inter%a_nm%item(i)%item(Ns_n+1:Ns_n+2) = fn_scs%a_nm%item(i)%item(Ns_m+1:Ns_m+2)! 
			f_inter%b_nm%item(i)%item(Ns_n+1:Ns_n+2) = fn_scs%b_nm%item(i)%item(Ns_m+1:Ns_m+2) ! 
			f_inter%a_nm%item(i+2)%item(Ns_n+1:Ns_n+2) = fn_scs%a_nm%item(i+2)%item(Ns_m+1:Ns_m+2)! 
			f_inter%b_nm%item(i+2)%item(Ns_n+1:Ns_n+2) = fn_scs%b_nm%item(i+2)%item(Ns_m+1:Ns_m+2) 
		end do
	end subroutine Lagrange_Interpolation_New 
	
	subroutine Lagrange_Interpolation_Weights(x, t, p, tp_arr, tp_arr_ext, wk)
		implicit none
	
		integer, intent(in):: t, p!
		real(dp), dimension(:), allocatable :: xx, wk, tp_arr, tp_arr_ext
		real(dp) :: x
		integer :: K_arr, j, kt, m	
		k_arr = size(tp_arr)
		allocate(xx(2*p))		
		kt = K_arr+t
		xx(1:2*p) = tp_arr_ext(kt+1-p:kt+p)	
		
		do j = 1, 2*p
			wk(j) = 1.0
			do m = 1, 2*p
				if (m .ne. j) then
					wk(j) = (x-xx(m))/(xx(j)-xx(m))*wk(j)
				end if			
			end do		
		end do
	end subroutine Lagrange_Interpolation_Weights

	subroutine Sorting_theta_phi_LenReduced_Arr_New(p, K_m, ff_in, ff_out)
		!use lib_sie_data_container		
		implicit none
		complex(dp), dimension(:, :, :), allocatable, intent(out) :: ff_out
		complex(dp), dimension(:, :, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m, p
		integer :: n_phi, m, n, K_n, j, n_theta
		real(dp), dimension(:), allocatable :: a_theta, a_phi
		complex(dp) :: Max_p(2), Max_t(2)
		real(dp) :: xt, xp
		complex(dp) :: fx_fy(2, 4)
		
		K_n = 2*K_m
		n_phi = K_n+2*p+1
		n_theta = K_m+2*p+2	
		
		call theta_phi_extended(K_m, a_phi, a_theta)	
	
		allocate(ff_out(n_theta, n_phi, 4))
		ff_out=(0.0, 0.0)
		
		! ++ fx and fy at theta = 0 and PI		
		do j = 1, 2!k1, k2
			fx_fy(j, 1:2) = (/fxy%a_nm%item(j)%item(1), fxy%b_nm%item(j)%item(1)/) !theta = 0, phi=0 fx, fy
			!New
			fx_fy(j, 3:4) = (/fxy%a_nm%item(j)%item(2), fxy%b_nm%item(j)%item(2)/) !theta = pi, phi=0 fx, fy
		end do
				
		! ++ f_theta and f_phi at theta = 0
		do m = 1, n_phi !phi
			n = m+K_n-p !the position of phi varied from 0 to 2*pi
			Max_t = (/cos(a_phi(n)), sin(a_phi(n))/)  !f_theta
			Max_p = (/-sin(a_phi(n)), cos(a_phi(n))/) !f_phi			
			do j = 1, 2 !k1 and k2
				ff_out(p+1, m, j) = -dot_product(Max_t, fx_fy(j, 3:4))  !!theta = pi, phi
				ff_out(K_m+p+2, m, j) = dot_product(Max_t, fx_fy(j, 1:2))  !!theta = 0, phi
				ff_out(p+1, m, j+2) = dot_product(Max_p, fx_fy(j, 3:4))    !!theta = pi, phi			
				ff_out(K_m+p+2, m, j+2) = dot_product(Max_p, fx_fy(j, 1:2)) !!theta = 0, phi
			end do
		end do
		
		ff_out(p+2:K_m+p+1, p+1:K_n+p, 1:4) = ff_in(1:K_m, 1:K_n, 1:4)
	
		do m = 1, K_m+2*p+2!
			!print*, 'm=', m
			xt = a_theta(m-p + K_m+1)
			do n = 1, K_n+2*p+1!n_phi
				xp = a_phi(n-p+K_n)
				if (xt<0.0)then
					if (xp<PI .and. xp>=0.0) then
						ff_out(m, n, 1:4) = -ff_out(2*(K_m+p+2)-m, n+K_n/2, 1:4)
					elseif (xp<2*PI .and. xp>=PI) then
						ff_out(m, n, 1:4) = -ff_out(2*(K_m+p+2)-m, n-K_n/2, 1:4)
					endif
				elseif (xt>=PI) then
					if (xp<PI .and. xp>=0.0)then					 
						ff_out(m, n, 1:4) = -ff_out(2*(p+1)-m, n+K_n/2, 1:4)
					elseif (xp<2*PI .and. xp>=PI)then
						ff_out(m, n, 1:4) = -ff_out(2*(p+1)-m, n-K_n/2, 1:4)
					end if
				elseif (xt<PI .and. xt>0.0)then
					if (xp<2*PI .and. xp>=0.0)then
						ff_out(m, n, 1:4) = ff_out(m, n, 1:4)
					end if
				endif
			end do
		end do
	
		do m = 1, k_m+2*p+2
			xt = a_theta(m-p+k_m+1)
			do n = 1, k_n+2*p+1
				xp = a_phi(n-p+k_n)
				if (xp>=2*pi)then
				ff_out(m, n, 1:4) = ff_out(m, n-k_n, 1:4)
				elseif (xp<0)then
				ff_out(m, n, 1:4) = ff_out(m, n+k_n, 1:4) !
				end if
			end do
		end do
	end subroutine Sorting_theta_phi_LenReduced_Arr_New
    
	subroutine Sorting_theta_phi_LenReduced_Arr_II(p, K_m, ff_in, ff_out)
		implicit none
		complex(dp), dimension(:, :, :), allocatable, intent(out) :: ff_out
		complex(dp), dimension(:, :, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m, p
		integer :: n_phi, m, n, K_n, j
		real(dp), dimension(:), allocatable :: a_theta, a_phi
		real(dp) :: Max_p(2), Max_t(2)
		real(dp) :: xt, xp
		complex(dp) :: fxy1(4, 2), fxy2(4, 2)
		
		K_n = 2*K_m
		n_phi = K_n+2*p+1
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		allocate(ff_out(K_m+2*p+2, K_n+2*p+1, 4))
		ff_out=(0.0, 0.0)
		
		do j = 1, 4
			fxy1(j, 1:2) =(/ff_in(K_m+1, 1, j), ff_in(K_m+1, 2, j)/) !theta = 0, k1 and k2
			fxy2(j, 1:2) =(/ff_in(K_m+2, 1, j), ff_in(K_m+2, 2, j)/) !theta = pi, k1, and k2
		end do
		
		! ++ f_theta and f_phi at theta = 0 and PI
		do m = 1, n_phi
			n = m-p+K_n
			Max_t = (/cos(a_phi(n)), sin(a_phi(n))/)  !theta = 0
			Max_p = (/-sin(a_phi(n)), cos(a_phi(n))/) !theta = pi
			do j = 1, 2
				ff_out(p+1, m, j) = -dot_product(Max_t, fxy2(j, :))  !
				ff_out(K_m+p+2, m, j) = dot_product(Max_t, fxy1(j, :))  !
				ff_out(p+1, m, j+2) = dot_product(Max_p, fxy2(j+2, :))  !
				ff_out(K_m+p+2, m, j+2) = dot_product(Max_p, fxy1(j+2, :) ) 
			end do
		end do
		
		ff_out(p+2:K_m+p+1, p+1:K_n+p, 1:4) = ff_in(1:K_m, 1:K_n, 1:4)
		
		do m = 1, K_m+2*p+2!		
			xt = a_theta(m-p + K_m+1)
			do n = 1, K_n+2*p+1!n_phi
				xp = a_phi(n-p+K_n)
				if (xt<0.0)then
					if (xp<PI .and. xp>=0.0) then
						ff_out(m, n, 1:4) = -ff_out(2*(K_m+p+2)-m, n+K_n/2, 1:4)
					elseif (xp<2*PI .and. xp>=PI) then
						ff_out(m, n, 1:4) = -ff_out(2*(K_m+p+2)-m, n-K_n/2, 1:4)
					endif
				elseif (xt>=PI) then
					if (xp<PI .and. xp>=0.0)then					 
						ff_out(m, n, 1:4) = -ff_out(2*(p+1)-m, n+K_n/2, 1:4)
					elseif (xp<2*PI .and. xp>=PI)then
						ff_out(m, n, 1:4) = -ff_out(2*(p+1)-m, n-K_n/2, 1:4)
					end if
				elseif (xt<PI .and. xt>0.0)then
					if (xp<2*PI .and. xp>=0.0)then
						ff_out(m, n, 1:4) = ff_out(m, n, 1:4)
					end if
				endif
			end do
		end do
		do m = 1, k_m+2*p+2
			xt = a_theta(m-p+k_m+1)
			do n = 1, k_n+2*p+1
				xp = a_phi(n-p+k_n)
				if (xp>=2*pi)then
				ff_out(m, n, 1:4) = ff_out(m, n-k_n, 1:4)
				elseif (xp<0)then
				ff_out(m, n, 1:4) = ff_out(m, n+k_n, 1:4) !
				end if
			end do
		end do
	end subroutine Sorting_theta_phi_LenReduced_Arr_II
    
	subroutine theta_phi_extended(K_m, a_phi, a_theta)
	! theta and phi for Lagrange interpolation
	! p is the one-side sampling point number closing the source point. 	
		real(dp), dimension(:), allocatable, intent(out) :: a_theta, a_phi
		integer, intent(in) :: K_m
		real(dp) :: d_phi, a, b
		integer :: Np, m_ph
		real(dp), dimension(:), allocatable :: xx, ww
		integer :: j, k		
				
		m_ph = 3*(2*K_m)
		Np = K_m*m_ph
		d_phi= 6*PI/m_ph
		
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(m_ph))
		a = -1.0
		b = 1.0
		call legendre_handle(K_m, a, b, xx, ww)				
		!Theta from 2*PI to -PI
		do j = 1, K_m
			a_theta(j) = 2*PI - acos(xx(K_m -(j-1))) !between 2*PI ~ PI		
			a_theta(j + K_m + 1) = acos(xx(j)) !theta between PI ~ 0
			a_theta(2*K_m + 2 + j) = -acos(xx(K_m-(j-1))) !between 0 ~ -PI
		end do
		a_theta(K_m+1)= PI
		a_theta(2*K_m+2) = 0.0
		!---------------------------------
		!when k=2*Km + 1, phi = 0 
		!when k=2*(2*Km + 1), phi=2*PI		
		!----------------------------------
		k = 0
		do k = 1, m_ph
			a_phi(k) = -2*PI + d_phi*(k-1) 
		end do
	end subroutine theta_phi_extended
	
	subroutine find_elements_in_field_zcsr(m, n, D_mat_tmp)!
		integer, intent(in) :: m, n !m: row		
		complex(dp), dimension(2,2), intent(out) :: D_mat_tmp
		integer :: is, it, i
		
		is = field_csr%ia(m)
		it = field_csr%ia(m+1)-1			
		do i = is, it		
			if (field_csr%ja(i) .eq. n) then
				D_mat_tmp(1, 1) = field_csr%a(i)
				D_mat_tmp(1, 2) = field_csr%b(i)
				D_mat_tmp(2, 1) = field_csr%c(i)
				D_mat_tmp(2, 2) = field_csr%d(i)
			else
				cycle	
			end if		
		end do
	end subroutine
	
	subroutine Truncation_number_calculation()
		implicit none
		
		! auxiliary
		integer(kind = 1) :: k
		integer :: K_m, n, d0!
		logical :: near_field
		real(dp) :: d_c, error, r_max, dd 
				
		error = 1.0e-4
		d0 = 4
		
		k = m_tree_l_max		
		d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
		box_at_l%lmax = m_tree_l_max
		box_at_l%dc = d_c
		
		if (allocated(truncation_number)) then
			deallocate(truncation_number)
		end if
		
		allocate(truncation_number(k))
		truncation_number(:)%n = 0
		truncation_number(:)%near_field = .false.
		k = 2
		d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)	
		dd = abs(k1)*d_c*sqrt(3.0) !diameter of the sphere 
		n = abs(dd) + 1.8*(d0)**2/3*abs(dd)**(1/3)
		!check why on server n=14 instead of 9
		
		!Assuming that only at the lowest level, k2 is considered. 
		!R_min is the smallest effective distance to consider k2			
		!if (imag(k2) .ne. 0.0) then
		!	k = m_tree_l_max-1
		!	call Truncation_number_near(error, d_c, r_max, near_field, K_m)
		!	truncation_number(m_tree_l_max)%n = K_m
		!	truncation_number(m_tree_l_max)%near_field = near_field
		!	do
		!		if (k .lt. 2) exit
		!		d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
		!		if (d_c*2 .lt. r_max) then 
		!			call Truncation_number_near(error, d_c, r_max, near_field, K_m)
		!			truncation_number(k)%n = K_m
		!			truncation_number(k)%near_field = near_field	
		!		else			
		!			call Truncation_number_far(error, d_c, truncation_number(k+1)%n, K_m) !
		!			truncation_number(k)%n = K_m
		!		end if
		!		k = k-1
		!	end do
		!	print*, 'truncation_number =',  truncation_number
		!else
		!	k = m_tree_l_max
		!	print*, 'k=', k
		!	do
		!		if (k .lt. 2) exit
		!		d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
		!		call Truncation_number_far(error, d_c, truncation_number(k)%n, K_m) !
		!		truncation_number(k)%n = K_m				
		!		k = k-1
		!	end do
		!end if
		
		do k = 1, m_tree_l_max
			truncation_number(k)%n = truncation_number_default! truncation_number(k)%n + 3 		!14!
			!print*, 'k=', k
			!print*, 'truncation_number =', truncation_number(k)%n
		end do
	end subroutine Truncation_number_calculation
end module lib_sie_mlfmm