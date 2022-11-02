module lib_sie_tri_mlfmm
	
	use omp_lib	
	use libmath
	use libmlfmm
	use lib_sie_constants
	use lib_sie_data_container
	use lib_sie_mlfmm
	
	use lib_sie_tri_data_container	
	use lib_sie_tri_calculation_mod	
	use Gauss_Quadrature_Fast

	implicit none	
	private 
    
	public :: Khat_and_TRF_ExtendedK 	
	public :: TL_km_arr_II 	
	public :: init_list_sie
	public :: MLFMM_tri_summation_New
	public :: Lagrange_Interpolation_New	
	public :: fx_fy	
	public :: Truncation_number_near
	!public :: Truncation_number_far
	public :: PreCalculate_RT	
	public :: LT3_get_TL_at_loc
	public :: Find_box_index_for_element
	
	!Test
	public :: Test_PreCalculate_RT_02
	public :: Test_PreCalculate_RT_01	
	public :: Test_lagrange_interpolation
		
	!interface init_list_sie
	!	module procedure init_list_sie_c
	!	module procedure init_list_sie_v
	!end interface
	
	contains
	
	subroutine Khat_and_TRF_II(K_m, K_hat_II, TRF_II)
		
		implicit none	
		integer, intent(in) :: K_m
		real(dp), dimension(:), allocatable :: xx, ww 
		real(dp), dimension(:), allocatable, intent(out) :: TRF_II
		real(dp), dimension(:, :), allocatable, intent(out) :: k_hat_II
		real(dp) :: d_phi, phi, a, b
		integer :: j, k, i, N_sp
		
		N_sp = K_m*K_m*2
		allocate(k_hat_II(N_sp, 3))
		allocate(xx(K_m))
		allocate(ww(K_m))		
		allocate(TRF_II(N_sp))
		a = -1.0
		b = 1.0		
		call legendre_handle(K_m, a, b, xx, ww)
		d_phi = 2*PI/(2*K_m)
		do j = 1, K_m
            do k = 1, 2*K_m
                phi = (k-1)*d_phi 
                i = k + (j-1)*2*K_m
                k_hat_II(i, :) = (/sin(acos(xx(j)))*cos(phi), sin(acos(xx(j)))*sin(phi), xx(j)/)
                TRF_II(i) = ww(j)*d_phi
            end do
		end do
		deallocate(xx)
		deallocate(ww)
	end subroutine Khat_and_TRF_II	

	subroutine Find_box_index_for_element(element_uindex)
		use libtree
		!use lib_sie_data_container
		use libmlfmm		
		
		implicit none		
		type(lib_tree_universal_index) :: uindex
		type(lib_tree_universal_index), dimension(:), allocatable, intent(out) :: element_uindex
      ! auxiliary
      integer(kind=UINDEX_BYTES) :: number_of_boxes, m_element_number
		type(lib_tree_data_element), dimension(:), allocatable :: data_element
		type(lib_tree_spatial_point) :: x_c
		
		integer(kind=4), dimension(:), allocatable :: element_number		
      integer(kind=UINDEX_BYTES) :: i, j, m_tree_l_max      
      integer(kind=1) :: hierarchy_type
		
		if (allocated(struc_tri%x_c)) then 
				deallocate(struc_tri%x_c)
		end if
		  
		allocate(struc_tri%x_c(m_pairs))
		allocate(element_uindex(m_pairs))
		
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		number_of_boxes = size(m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index)	
		uindex%l = m_tree_l_max
		if (m_ml_fmm_hierarchy(m_tree_l_max)%is_hashed) then			
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((uindex%n .ge. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
							!**********************************************							
							data_element = lib_tree_get_domain_e1(uindex, element_number)
							if ((allocated (data_element)) &
									.and. (size(data_element) .gt. 0)) then
									do j = 1, size(element_number)
										m_element_number = element_number(j)
										x_c = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex))
										struc_tri%x_c(m_element_number)%point = x_c%x
										element_uindex(m_element_number) = uindex
									end do								
							end if				
                end if
            end do           
        else
            do i=1, number_of_boxes
					uindex%n = i - int(1, 1)                
					data_element = lib_tree_get_domain_e1(uindex, element_number)
					if ((allocated (data_element)) &
							.and. (size(data_element) .gt. 0)) then
						do j = 1, size(element_number)
							m_element_number = element_number(j)
							x_c = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex))
							struc_tri%x_c(m_element_number)%point = x_c%x
							element_uindex(m_element_number) = uindex
						end do								
					end if
            end do            
        end if
		  if (allocated(data_element)) then 
				deallocate(data_element)
		  end if
		  if (allocated(element_number)) then 
				deallocate(element_number)
		  end if
	end subroutine
		
	subroutine Truncation_number_near(error, d_c, r_near, near_field, K_m)
		
		implicit none		
		real(dp), intent(in) :: error, d_c
		complex(dp) :: kd
		logical, intent(out) :: near_field
		real(dp), dimension(8, 3) :: ra_sub, rb_sub!
		
		real(dp) :: d0 !
		real(dp), intent(out) :: r_near

      integer :: n!
		integer, intent(out) :: K_m
		complex(dp) :: dd
		
		d0 = log10(error)
		call two_cubes_containing_subcubes(d_c*2, ra_sub, rb_sub)
		kd = k0
		dd = kd*d_c*sqrt(3.0) !diagonal of the cube
		n = abs(dd) + 1.8*(d0)**2/3*abs(dd)**(1/3)
		print*, 'n in medium k0', n		
		
		!Using the truncation number from k1 for k2
		!error_max = error
		!r_near = 0.0
		!kd = k2
		!counter = 0
		!!The far most point in cube-a
		!ro = ra_sub(1, :)
		!ra = ro
		!roa = ro-ra		
		!do
		!	if (error_max .gt. error) exit
		!	do j = 1, 8
		! 		rb = rb_sub(j, :)
		!		rp = rb_sub(j, :) + (/1.0, 1.0, 1.0/)*d_c/2
		!		rab = ra - rb
		!		K_m = n
		!		error_tmp = sampling_error(K_m, ra, rb, ro, rp, kd)				
		!		r_tmp = abs(vec_len(rab))
		!		
		!		!initialize r_near
		!		if ((error_tmp .ge. error) .and. (counter .eq. 0)) then !
		!			print*, 'error =', error_tmp
		!			error_max = error_tmp
		!			r_near = r_tmp
		!			counter = counter + 1
		!		end if
		!		
		!		!find the smallest rab suitable for taking k2 into account, at level lmax			
		!		!find the largest error using this n
		!		if ((error_tmp .ge. error) .and. (counter .ge. 1)) then
		!			counter = counter + 1				
		!			if (r_tmp .lt. r_near) then
		!				r_near = r_tmp					
		!			end if
		!			if (error_tmp .gt. error_max) then
		!				error_max = error_tmp					
		!			end if
		!		end if				
		!	end do !loop for points in cube-b
		!	if (error_max .gt. error*5) then
		!		print*, 'not a near-field'
		!		near_field = .false.
		!	elseif (error_max .le. error*5) then
		!		near_field = .true.
		!	end if
		!	n = n-1
		!end do
		K_m = n
		return
	end subroutine Truncation_number_near
	
	!subroutine Truncation_number_far(error, d_c, n_lmin, K_m) !
	!	
	!	implicit none		
	!	real(dp), intent(in) :: error, d_c
	!	integer, intent(in) :: n_lmin
	!	complex(dp) :: kd
	!	complex(dp) :: dd
	!	
	!	real(dp) :: dc, d0
	!	real(dp) :: r_near
	!	real(dp), dimension(8, 3) :: ra_sub, rb_sub
	!	
 !     integer :: n!
	!	integer, intent(out) :: K_m
	!	
	!	d0 = log10(error)			
	!	!print*, 'For medium k2, the maximum rab is lambda/8 or lambda/9'
	!	!print*, 'Lambda to edge length of the smallest box =', NINT(lambda*2/d_c)
	!	
	!	call two_cubes_containing_subcubes(d_c*2, ra_sub, rb_sub)		
	!	kd = k2
	!	dc = d_c
	!	dd = kd*dc*sqrt(3.0) !diameter of the sphere 
	!	n = abs(dd) + 1.8*(d0)**2/3*abs(dd)**(1/3)
	!	K_m = n
	!	!consider n_lmin to avoid larger l having more sampling points
	!	!if (n .le. n_lmin) then
	!	!	n = n_lmin
	!	!	do j = 1, 8
	!	!		rb = rb_sub(j, :)
	!	!		rp = rb_sub(j, :) + (/1.0, 1.0, 1.0/)*dc/2
	!	!		rab = ra - rb
	!	!		rpb = rp - rb				
	!	!		K_m = n
	!	!		error_tmp = sampling_error(K_m, ra, rb, ro, rp, kd)				
	!	!		r_tmp = abs(vec_len(rab))
	!	!		if (error_tmp .ge. error) then !	
	!	!			print*, 'error is =', error_tmp
	!	!			print*, 'Recalculate n_lmin'
	!	!		end if
	!	!	end do
	!	!else
	!	!	ro = ra_sub(1, :)
	!	!	ra = ro
	!	!	roa = ro-ra
	!	!	do
	!	!		do j = 1, 8
	!	!			rb = rb_sub(j, :)
	!	!			rp = rb_sub(j, :) + (/1.0, 1.0, 1.0/)*dc/2
	!	!			rab = ra - rb
	!	!			rpb = rp - rb				
	!	!			K_m = n
	!	!			error_tmp = sampling_error(K_m, ra, rb, ro, rp, kd)				
	!	!			r_tmp = abs(vec_len(rab))
	!	!			if (error_tmp .ge. error) then !
	!	!				!print*, 'the smallest n=', n
	!	!				!print*, 'error =', error_tmp				
	!	!			end if
	!	!		end do		
	!	!		if ((error_tmp .gt. error) .or. (n .eq. n_lmin)) exit
	!	!		n = n-1
	!	!	end do
	!		!K_m = n 
	!	!end if			
	!	return
	!end subroutine Truncation_number_far
	
	subroutine two_cubes_containing_subcubes(dc, ra_sub, rb_sub)
		implicit none			
		real(dp), dimension(8, 3) :: r_unit_sub, r_sub
		real(dp), dimension(8, 3), intent(out) :: ra_sub, rb_sub
		real(dp), dimension(3) :: r_unit, ra, rb, r_shift
		real(dp), intent(in) :: dc
		integer :: i
			
		r_unit = (/0.0, 0.0, 0.0/)
		r_unit_sub(1, :) = (/0.25, 0.25, 0.25/) !center of the local box
		r_unit_sub(2, :) = (/0.25, -0.25, 0.25/)
		r_unit_sub(3, :) = (/0.25, 0.25, -0.25/)
		r_unit_sub(4, :) = (/0.25, -0.25, -0.25/)
		r_unit_sub(5, :) = (/-0.25, 0.25, 0.25/)
		r_unit_sub(6, :) = (/-0.25, -0.25, 0.25/)
		r_unit_sub(7, :) = (/-0.25, 0.25, -0.25/)
		r_unit_sub(8, :) = (/-0.25, -0.25, -0.25/)
		r_sub = r_unit_sub*dc
		r_shift = (/0.5, 0.5, 0.5/)*dc
		!r_shift = (/0.5, 0.0, 0.0/)*dc
		ra = (r_unit + r_shift)
		rb = (r_unit - r_shift)
		do i = 1, 8			
			ra_sub(i, :) = r_sub(i, :) + r_shift(:)
			rb_sub(i, :) = r_sub(i, :) - r_shift(:)
		end do
	end subroutine
	
	!function sampling_error(K_m, ra, rb, ro, rp, kd) result(error)
	!
	!	use lib_sie_data_container
	!	use libmath
	!			
	!	implicit none
	!	integer, intent(in) :: K_m
	!	integer :: n, N3, m
	!	complex(dp) :: tmp, sum_a
	!	complex(dp), intent(in) :: kd
	!	complex(dp), dimension(:), allocatable :: int_2, TRF_t, ff
	!	real(dp) :: error
	!	real(dp), intent(in) :: ra(3), rb(3), rp(3), ro(3)
	!	! intrinsic exp
	!				
	!	n = K_m		
	!	N3 = 2*K_m**2
	!	tmp = exp(-im*kd*vec_len(ro-rp))/(4*pi*vec_len(ro-rp))
	!	call Khat_and_TRF_Legendre(K_m)		
	!	!-------Receive and radiation functions -------
	!	int_2 = TL_km_arr(N3, ra-rb, n, kd)		
	!	TRF_t = Tf_rd(N3, rp-rb, kd) 
	!	ff = Rf_rc(N3, ro-ra, kd) !
	!	sum_a = cmplx(0.0d0, 0.0d0)
	!	do m = 1, N3         
	!		sum_a = sum_a + ff(m)*TRF_t(m)*int_2(m)
	!	end do
	!	error = abs((tmp-sum_a)/tmp)					
	!end function
	
	
	subroutine Lagrange_Interpolation(pp, K_m, K_D, Fn_scs, ff_inter)
		use libmath		
		implicit none
		
		integer, intent(in) :: pp, K_m, K_D
		integer :: i, j, k, K_n, tt
		integer :: it, jp, m_phi, n_phi, t, s, kp, kt
		
		!type(Integration_fn_scs) :: Int_fn_inter, Int_fn
		real(dp) :: dphi, a, b, theta, phi
		
		complex(dp), dimension(:, :), allocatable, intent(out) :: ff_inter
		complex(dp), dimension(:, :), allocatable, intent(in) :: Fn_scs
		complex(dp), dimension(:, :, :), allocatable :: ff_in, ff_out		
		real(dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		
	   K_n = K_m + K_D
		n_phi = 2*K_n
		m_phi = 2*K_m
		dphi = 2*PI/m_phi		
		allocate(ff_in(K_m+2, 2*K_m, 4))!, ff_out(K_m+2, 2*K_m, 4)
		allocate(ff_inter(n_phi*K_n, 4)) 
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_in(i, j, 1:4)= Fn_scs(k, 1:4) !theta in k1		Why no in k2?		
			end do
		end do
		do i = 1, 2		
			ff_in(K_m+1:K_m+2, 1, i) =	  Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i) !theta
			ff_in(K_m+1:K_m+2, 1, i+2) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i) !theta		
			ff_in(K_m+1:K_m+2, 2, i) =   Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i+2)	!phi	
			ff_in(K_m+1:K_m+2, 2, i+2) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i+2)	!phi	
		end do
		
		call Sorting_theta_phi_LenReduced_Arr_II(pp, K_m, ff_in, ff_out)				
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
				ff_inter(tt, 1:4) = (0.0, 0.0)
				phi = (jp-1)*2*PI/n_phi
				s = floor((jp-1.0)/K_n*K_m)+1	
				kp = pp + s
				call Lagrange_Interpolation_Weights(phi, s, pp, xp_old, a_phi, wk_phi)
				do i = 1, 2*pp
					do j = 1, 2*pp						
						ff_inter(tt, 1:4) = ff_inter(tt, 1:4) + wk_phi(j)*wk_theta(i)*ff_out(kt+i-pp, kp+j-pp, 1:4)
					end do
				end do				
			end do
        end do
	end subroutine Lagrange_Interpolation
	
	subroutine MLFMM_tri_summation_New(K_m, f_rl, fD, sum_r)
		
      !sum_a = K1+K2
      !sum_b = L1+L2
      !sum_c = e1/u1*L1 + e2/u2*L2
		! Z_LE =zz(1) = sum_b
		! Z_KH = zz(2)= -sum_a
		! Z_KE = zz(3) = sum_a
		! Z_LH = zz(4) = sum_c
      !Equation 28 in Mitchang's paper Vol. 11, P1383, JOSAA, 1994
      !K, L operators for FMM are Eq.9.25-9.33 in Gibson, p336
		integer, intent(in) :: K_m
      type(lib_ml_fmm_coefficient), intent(in) :: f_rl, fD 
		
		! auxiliary 
		type(lib_ml_fmm_coefficient) :: Z_LE, Z_LH, Z_KH, Z_KE, f_rk
		complex(dp), dimension(2), intent(out) :: sum_r
		complex(dp), dimension(2) :: eps_d, my_d, kd	
		complex(dp) :: sum_a, sum_b, sum_c
		complex(dp) :: a_lc(2), d_lc(2)
      integer :: Ns, i
		
		eps_d = (/eps_1, eps_2/)*im*Omega
		
		call formulation_coefficient(a_lc, d_lc)
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		my_d = (/my_1, my_2/)*im*Omega
		kd = (/k1, k2/)*im		
		Ns = K_m*K_m*2		
	
		call init_list_sie(Z_LE, 2, Ns+2)!Z_LH, Z_KH, Z_KE
		call init_list_sie(Z_LH, 2, Ns+2)
		call init_list_sie(Z_KE, 2, Ns+2)
		call init_list_sie(Z_KH, 2, Ns+2)
		call init_list_sie(f_rk, 2, Ns+2)		
		do i = 1, 2
			f_rk%a_nm%item(i)%item =-f_rl%b_nm%item(i)%item
			f_rk%b_nm%item(i)%item = f_rl%a_nm%item(i)%item
			
			Z_LE%a_nm%item(i)%item = f_rl%a_nm%item(i)%item*fD%a_nm%item(i)%item*my_d(i)!
			Z_LE%b_nm%item(i)%item = f_rl%b_nm%item(i)%item*fD%b_nm%item(i)%item*my_d(i)!
			
			!-sum_a
			Z_KE%a_nm%item(i)%item = -f_rk%a_nm%item(i)%item*fD%a_nm%item(i)%item*kd(i)			
			Z_KE%b_nm%item(i)%item = -f_rk%b_nm%item(i)%item*fD%b_nm%item(i)%item*kd(i)!
			
			Z_LH%a_nm%item(i)%item = f_rl%a_nm%item(i)%item*fD%a_nm%item(i+2)%item*eps_d(i)!
			Z_LH%b_nm%item(i)%item = f_rl%b_nm%item(i)%item*fD%b_nm%item(i+2)%item*eps_d(i)!
			
			Z_KH%a_nm%item(i)%item = f_rk%a_nm%item(i)%item*fD%a_nm%item(i+2)%item*kd(i)
			Z_KH%b_nm%item(i)%item = f_rk%b_nm%item(i)%item*fD%b_nm%item(i+2)%item*kd(i)!
		end do
				
		sum_r(1:2) = (0.0, 0.0)
		do i = 1, 2
			sum_r(1) = sum_r(1) + sum(Z_KH%a_nm%item(i)%item + Z_KH%b_nm%item(i)%item) + sum(Z_LE%a_nm%item(i)%item + Z_LE%b_nm%item(i)%item)
			sum_r(2) = sum_r(2) + sum(Z_KE%a_nm%item(i)%item + Z_KE%b_nm%item(i)%item) + sum(Z_LH%a_nm%item(i)%item + Z_LH%b_nm%item(i)%item)			
		end do			
		
		select case(pre_types%formulation)
			case('PMCHWT')
					sum_r(1) = sum_r(1) 
					sum_r(2) = sum_r(2) 				
			case('PMCHWT-S')
				sum_r(2) = sum_r(2)*eta_1
			case('MCTF')
				sum_r(2) = sum_r(2)*eta_1*eta_2
			case('ICTF')
				sum_r(1) = sum_r(1)/eta_a
				sum_r(2) = sum_r(2)*eta_a				
			!case('MCTF2')
			!	sum_r(1) = sum_r(1)/(eta_1*eta_2)
			!	sum_r(2) = sum_r(2)*(eta_1*eta_2)	
			case('MCTF2')
				sum_r(1) = sum_r(1)*a_lc(1)/eta_1
				sum_r(2) = sum_r(2)*d_lc(1)*eta_1
		end select
		
	end subroutine
	!++++++++

	!precalculate the radiation or transmission coefficient
	!as a function of theta and phi
	!and save for use
	subroutine PreCalculate_RT(truncation_n, struct, preCalc_RT)
		use lib_sie_tri_data_container
		use lib_sie_tri_calculation_mod
		
		implicit none
		type(truncation_hl), dimension(:), intent(in) :: truncation_n
		type(lib_ml_fmm_coefficient), dimension(:), allocatable, intent(out) :: preCalc_RT
		type(structure_tri), intent(in) :: struct
		!dummy
		integer :: K_m, lmax
		integer :: m, ll, Ns!, m_pairs
		type(lib_ml_fmm_coefficient) :: fn_scs		
		real(dp) :: Rc(3)
		
		lmax = lib_tree_get_level_max(tree_s_opt)
		K_m = truncation_n(lmax)%n		
		m_pairs = size(struct%neighbours)
		Ns = K_m*K_m*2
		if (allocated(preCalc_RT)) then
			deallocate(preCalc_RT)
		end if
		
		allocate(preCalc_RT(m_pairs))
		call Khat_and_TRF_ExtendedK(K_m)
		do m = 1, m_pairs			
			call init_list_sie(preCalc_RT(m), 2, Ns+2)
			!should be the same for ll = 1 and 2, has to be checked.
			Rc = struct%x_c(m)%point
			do ll = 1, 2
				call RT_function_SCSArr_New(K_m, struct, m, ll, Rc, fn_scs)
				preCalc_RT(m)%a_nm%item(1) = preCalc_RT(m)%a_nm%item(1) + fn_scs%a_nm%item(1)
				preCalc_RT(m)%a_nm%item(2) = preCalc_RT(m)%a_nm%item(2) + fn_scs%a_nm%item(2) !theta
				preCalc_RT(m)%b_nm%item(1) = preCalc_RT(m)%b_nm%item(1) + fn_scs%b_nm%item(1) !phi
				preCalc_RT(m)%b_nm%item(2) = preCalc_RT(m)%b_nm%item(2) + fn_scs%b_nm%item(2)
			end do
		end do		
	end subroutine PreCalculate_RT
	 
	!Receive function Eq.(9.28) and (9.33) without the prefactors 
	!In Gibson's book, 2nd Version
	subroutine RT_function_SCSArr_New(K_m, struct, n, ll, Rc, fn_scs)		
		implicit none		
		real(dp), intent(in) :: Rc(3)		
		integer, intent(in) :: n, K_m, ll  !sampling point on a sphere				
		
		integer :: i, j, pos_neg, Ns
		real(dp), dimension(100) :: w(100), a(100), b(100)
		real(dp) :: pm1(3), pm2(3), vm_t(3), Roc(3), Ro(3), fm(3)
		complex(dp) :: ff1, ff2, tmp
		complex(dp) :: sum_a(3), sum_b(3)
		real(dp) :: area, Len_m, dfm	
		type(structure_tri), intent(in) :: struct	
		type(lib_ml_fmm_coefficient), intent(out) :: fn_scs
		
		call Quadrature_tri(ng, a, b, w)	
		call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
		
		pos_neg = -1*(((ll-1)*ll)-1)!
		dfm = pos_neg*Len_m 			
		Ns = K_m*2*K_m
		call init_list_sie(fn_scs, 2, Ns+2)
		call Khat_and_TRF_ExtendedK(K_m)
		do i = 1, Ns+2
			sum_a = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			sum_b = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			do j = 1, ng
				Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
				fm = 0.5*(Ro-vm_t)*dfm
				Roc = Ro - Rc
				tmp = dot_product(k_hat_p(i, :), Roc(:))
				ff1 = exp(-im*k1*tmp)*w(j)
				ff2 = exp(-im*k2*tmp)*w(j)
				sum_a = sum_a + fm*ff1
				sum_b = sum_b + fm*ff2
			end do			
			fn_scs%a_nm%item(1)%item(i) = dot_product(transfer_theta(i, 1:3), sum_a) !k1, at Ns+1: fx1
			fn_scs%a_nm%item(2)%item(i) = dot_product(transfer_theta(i, 1:3), sum_b) !k2, at Ns+1: fx2
			fn_scs%b_nm%item(1)%item(i) = dot_product(transfer_phi(i, 1:3), sum_a)   !k1, at Ns+1: fy1
			fn_scs%b_nm%item(2)%item(i) = dot_product(transfer_phi(i, 1:3), sum_b)   !k2, at Ns+1: fy2
		end do !
	end subroutine RT_function_SCSArr_New
 
	function TL_km_arr_II(r_ab, K_m)result(TL)
      integer, intent(in) :: K_m    !Lm is the order limit of the polynomials, N is the sampling points on the sphere
      real(dp), intent(in) :: r_ab(3)
      real(dp) :: r_hat(3), x
		real(dp), dimension(:, :), allocatable :: pm_arr
      complex(dp) :: r_h(2), kd(2)
		!real(dp), dimension(:, :), allocatable :: k_hat_II
		!real(dp), dimension(:), allocatable :: TRF_II
      real (dp), dimension(K_m) :: pm, dummy 
      complex(dp), dimension(:, :), allocatable :: TL
      complex(dp), dimension(K_m) :: hl
		
      integer :: j, fnu, m, i, Ns
		m = K_m
      fnu = 0        
      r_hat = r_ab/vec_len(r_ab)
		kd = (/k1, k2/)
      r_h =  kd*vec_len(r_ab) 
		Ns = 2*K_m*K_m
		
		if (allocated(TL)) then
			deallocate(TL)
		end if
		
		allocate(pm_arr(Ns+2, m), TL(Ns+2, 2))
		TL(1:Ns+2, 1:2) = (0.0, 0.0)   
		call Khat_and_TRF_ExtendedK(K_m)!
      do j = 1, Ns+2
         x = dot_product(r_hat(:), k_hat_p(j, :))
         call lib_math_legendre_polynomial(x, fnu, m, pm, dummy)
         pm_arr(j, :) = pm(:)
      end do
		do i = 1, 2
			hl =  lib_math_hankel_spherical_2(r_h(i), fnu, m) 
			do j = 1, m			
				TL(:, i) = TL(:, i) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd(i)/(4*PI)**2 ! 
			end do
		end do
      return 
		deallocate(pm_arr)
	end function TL_km_arr_II
	
	
	
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!Test functions
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 
	subroutine test_Hl2()
		
		!type(list_list_real) :: hl_1, hl_2 
		complex(dp), dimension(:), allocatable :: hl_2
		complex(dp) :: r_h, kd
		integer :: m, fnu, i, N
		real(dp) :: x
		! with fnu = 1, m=50 for instance, hl is the same as with the calculator for x=1.86, v=50
		! with fnu = 0, m=50 is identical with hl of the order v=49. 
		! But when fnu = 0, it does not work for m = 1 
		fnu = 0 
		m = 30
		N = 100
		kd = k1
		allocate(hl_2(m+1))
		open (unit=206, file = 'Hankel_2_m30.txt', action="write",status = 'replace')				
		do i = 2, N
			x = (0.1 + (i-1)*10/N)*illumination_p%lambda
			!print*, 'x=', x
		   r_h =kd*x
			hl_2 = lib_math_hankel_spherical_2(r_h, fnu, m+1)			
			write (206, '(201(es19.12, tr5))') x, abs(r_h), real(hl_2(m+1)), imag(hl_2(m+1)) !real(r_h), imag(r_h), 			
		end do			
		close(206)
	end subroutine
	
	subroutine Test_PreCalculate_RT_02(struct)
		implicit none
		!use lib_sie_data_container
		type(structure_tri), intent(in) :: struct
		integer :: K_m, K_n, n, m, Ns, i, ngp, j
		complex(dp), dimension(:, :), allocatable :: TL
		complex(dp) :: sum_a, sum_b, sum_c, kd(2), my_d(2), eps_d(2)!, tmp
		type(lib_ml_fmm_coefficient) :: fn_scs_r, fn_scs_t
		integer, dimension(:), allocatable:: n_arr
		
		real(dp) :: ra(3), rb(3), rab(3), d_c		
		!type(lib_tree_universal_index):: uindex		
		
		type(lib_ml_fmm_coefficient) :: C_1, B_i_2, B_i_1, D		
		
		K_m = 13
		m = 1
		Ns = K_m*K_m*2
		ngp = 3
		ra = struct%x_c(m)%point
		
		call init_list_sie(C_1, 4, Ns+2)
		call init_list_sie(B_i_1, 4, Ns+2)
		call init_list_sie(B_i_2, 4, Ns+2)
		call init_list_sie(D, 4, Ns+2)
		
		call Khat_and_TRF_ExtendedK(K_m)
		fn_scs_r = preCalculated_RT(m)
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		j = 3
		d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
		print*, 'd_c', d_c
		
		allocate(n_arr(2))
		!n_arr = (/220, 221, 222/)!index%n = 192, l=3
		!n_arr = (/42, 217, 219, 467, 468/) !index%n = 193
		n_arr = (/73, 529, 530, 644, 646/) !uindex%n = 138
		!n_arr = (/74, 314, 315, 531, 532, 706, 709/) !uindex%n = 139
		!n_arr = (/335, 654, 716/)
		n_arr= (/544, 715/) ! %n = 137
		!n_arr = (/316, 317, 318, 712, 73, 529, 530, 644, 646, 74, 314, 315, 531, 532, 706, 709/)
		do j = 1, size(n_arr)
		n = n_arr(j)
		C_1 = preCalculated_RT(n)
		fn_scs_t = preCalculated_RT(n)
		rb = struct%x_c(n)%point
		rab = ra -rb
		!print*, 'r_ab in FMM', vec_len(rab)
		TL = TL_km_arr_II(rab, K_m)
		kd=(/k1, k2/)*im
		my_d = (/my_1, my_2/)*im*Omega
		eps_d = (/eps_1, eps_2/)*im*Omega
		
		do i = 1, 2
			sum_a = sum_a - sum((-fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) &
				+ fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*kd(i)) !
			sum_b = sum_b + sum((fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) &
				+ fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*my_d(i))!
			sum_c = sum_c + sum((fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) &
				+ fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*eps_d(i))!
		end do
		end do	
		
		!print*, 'fn_scs_r%a_nm%item(1)%item(:) =', fn_scs_t%a_nm%item(1)%item(10)
		!print*, 'fn_scs_r%a_nm%item(2)%item(:) =', fn_scs_t%a_nm%item(2)%item(10)
		!print*, 'fn_scs_r%b_nm%item(1)%item(:) =', fn_scs_t%b_nm%item(1)%item(10)
		!print*, 'fn_scs_r%b_nm%item(2)%item(:) =', fn_scs_t%b_nm%item(2)%item(10)
		!
			
		print*, 'sum_a FMM', sum_a		
		print*, 'sum_b FMM', sum_b
		print*, 'sum_c FMM', sum_c		
		!print*, 'TL(:, 1)', TL(10, 1)
		!print*, 'TL(:, 2)', TL(10, 2)
		!++++++++++++++++++++++++++
		
	end subroutine Test_PreCalculate_RT_02
	
	subroutine Test_PreCalculate_RT_01(struct)
		implicit none
		!use lib_sie_data_container
		type(structure_tri), intent(in) :: struct
		integer :: K_m, n, m, Ns, i, ngp, n_arr(5), j
		complex(dp), dimension(:, :), allocatable :: TL
		complex(dp) :: sum_a, sum_b, sum_c, kd(2), my_d(2), eps_d(2), tmp, tmp_k(2)
		type(lib_ml_fmm_coefficient) :: fn_scs_r, fn_scs_t
		
		real(dp) :: ra(3), rb(3), rab(3)
		
		type(lib_tree_spatial_point) :: r_a, r_b, r_ab, r_o, r_p, r_xy, r_pb, r_oa
		type(lib_tree_universal_index):: uindex
		
		type(lib_ml_fmm_coefficient) :: C_1, B_i_2, B_i_1, D
		
		
		K_m = 13
		!n = 43! 470!469!230!229 ! 
		n_arr = (/43, 229, 230, 469, 470/)
		m = 1
		Ns = K_m*K_m*2
		ngp = 3
		ra = struct%x_c(m)%point
		
		call init_list_sie(C_1, 4, Ns+2)
		call init_list_sie(B_i_1, 4, Ns+2)
		call init_list_sie(B_i_2, 4, Ns+2)
		call init_list_sie(D, 4, Ns+2)
		
		call Khat_and_TRF_ExtendedK(K_m)
		fn_scs_r = preCalculated_RT(m)
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		do j = 1, 1!5
		n = 19!n_arr(j)
		C_1 = preCalculated_RT(n)
		fn_scs_t = preCalculated_RT(n)
		rb = struct%x_c(n)%point
		rab = ra -rb
		!print*, 'r_ab in FMM', vec_len(rab)
		TL = TL_km_arr_II(rab, K_m)
		kd=(/k1, k2/)*im
		my_d = (/my_1, my_2/)*im*Omega
		eps_d = (/eps_1, eps_2/)*im*Omega
		
		do i = 1, 2
			sum_a = sum_a - sum((-fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) + fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*kd(i))
			sum_b = sum_b + sum((fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) + fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*my_d(i))
			sum_c = sum_c + sum((fn_scs_r%a_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%a_nm%item(i)%item(:)) + fn_scs_r%b_nm%item(i)%item(:)*TL(:, i)*conjg(fn_scs_t%b_nm%item(i)%item(:)))*TRF_p(:)*eps_d(i))	
		end do
		end do		
		print*, 'sum_a FMM', sum_a
		print*, 'sum_b FMM', sum_b
		print*, 'sum_c FMM', sum_c		
		!++++++++++++++++++++++++++
		!++++++++++++++++++++++++++
		!upward
		uindex%l = 3
		uindex%n = 19
		r_p = lib_tree_get_centre_of_box(uindex)!child position
		
		uindex = lib_tree_get_parent(uindex)
		!print*, 'Test parent index', uindex
		
		uindex%l = 2
		uindex%n = 10
		r_b = lib_tree_get_centre_of_box(uindex) 		
		r_pb = lib_tree_get_unscaled_point(r_p) - lib_tree_get_unscaled_point(r_b) 		
		
		!downward
		uindex%l = 2
		uindex%n = 63
		r_a = lib_tree_get_centre_of_box(uindex)
		uindex%l = 3
		uindex%n = 504
		r_o = lib_tree_get_centre_of_box(uindex)	!child	
		r_oa = lib_tree_get_unscaled_point(r_o) - lib_tree_get_unscaled_point(r_a)
		r_ab = lib_tree_get_unscaled_point(r_a) - lib_tree_get_unscaled_point(r_b)
		
		r_xy = r_oa + r_ab - r_pb
		
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		TL = TL_km_arr_II(r_ab%x, K_m)
		!upward
		do i = 1, Ns+2
			tmp = dot_product(K_hat_p(i, :), r_pb%x)
			tmp_k = (/exp(im*k1*tmp), exp(im*k2*tmp)/)
			do j = 1, 2
				B_i_1%a_nm%item(j)%item(i) = tmp_k(j)*conjg(C_1%a_nm%item(j)%item(i))*TRF_p(i)*TL(i, j)
				B_i_1%b_nm%item(j)%item(i) = tmp_k(j)*conjg(C_1%b_nm%item(j)%item(i))*TRF_p(i)*TL(i, j)
				B_i_1%a_nm%item(j+2)%item(i) = tmp_k(j)*conjg(C_1%a_nm%item(j)%item(i))*TRF_p(i)*TL(i, j)
				B_i_1%b_nm%item(j+2)%item(i) = tmp_k(j)*conjg(C_1%b_nm%item(j)%item(i))*TRF_p(i)*TL(i, j)
			end do
		end do
	
		!downward
		do i = 1, Ns+2
			tmp = dot_product(K_hat_p(i, :), r_oa%x)
			tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
			do j = 1, 2
				D%a_nm%item(j)%item(i) =  -tmp_k(j)*B_i_1%a_nm%item(j)%item(i)*fn_scs_r%b_nm%item(j)%item(i) !
				D%b_nm%item(j)%item(i) = tmp_k(j)*B_i_1%b_nm%item(j)%item(i)*fn_scs_r%a_nm%item(j)%item(i) !
				D%a_nm%item(j+2)%item(i) = - tmp_k(j)*B_i_1%a_nm%item(j+2)%item(i)*fn_scs_r%b_nm%item(j)%item(i)
				D%b_nm%item(j+2)%item(i) = tmp_k(j)*B_i_1%b_nm%item(j+2)%item(i)*fn_scs_r%a_nm%item(j)%item(i)
			end do
		end do
		
		do j = 1, 2
			sum_a = sum_a + sum(D%a_nm%item(j)%item + D%b_nm%item(j)%item)		
			sum_b = sum_b + sum(D%a_nm%item(j+2)%item + D%b_nm%item(j+2)%item)
		end do 	
		!++++++++++++++
	end subroutine Test_PreCalculate_RT_01	
	
	subroutine Test_lagrange_interpolation(struct)
		type(structure_tri), intent(in) :: struct
		integer :: K_m, K_n, ll, pp, n, i, Ns	
		type(lib_ml_fmm_coefficient) :: fn_scs, fn_scs_tmp		
		type(lib_ml_fmm_coefficient) :: f_inter
		
		real(dp) :: Rc(3), a_arr(16)
		ll = 1
		K_m = 14
		K_n = 11
		pp = 3
		n = 2
		Ns = K_m*K_m*2
		Rc = struct%x_c(n)%point
		
		call init_list_sie(fn_scs_tmp, 4, Ns+2)		
		call Khat_and_TRF_ExtendedK(K_m)
		
		!call RT_function_SCSArr_New(K_m, struct, n, ll, Rc, fn_scs)
		!do i = 1, 2
		!	fn_scs_tmp%a_nm%item(i)= fn_scs%a_nm%item(i)
		!	fn_scs_tmp%b_nm%item(i)= fn_scs%b_nm%item(i)
		!	fn_scs_tmp%a_nm%item(i+2)= fn_scs%a_nm%item(i)
		!	fn_scs_tmp%b_nm%item(i+2)= fn_scs%b_nm%item(i)
		!end do 
		!
		!!To test whether the save format makes mistakes
		!open (unit=206, file = 'A_i_2_km14.txt', action="write",status = 'replace')
		!do i = 1, K_m*K_m*2+2
		!	write (206, '(16(es19.12, tr5))') 			&
		!	real(fn_scs_tmp%a_nm%item(1)%item(i)), imag(fn_scs_tmp%a_nm%item(1)%item(i)), &
		!	real(fn_scs_tmp%a_nm%item(2)%item(i)), imag(fn_scs_tmp%a_nm%item(2)%item(i)), &
		!	real(fn_scs_tmp%b_nm%item(1)%item(i)), imag(fn_scs_tmp%b_nm%item(1)%item(i)), &
		!	real(fn_scs_tmp%b_nm%item(2)%item(i)), imag(fn_scs_tmp%b_nm%item(2)%item(i)), &
		!	real(fn_scs_tmp%a_nm%item(3)%item(i)), imag(fn_scs_tmp%a_nm%item(3)%item(i)), &
		!	real(fn_scs_tmp%a_nm%item(4)%item(i)), imag(fn_scs_tmp%a_nm%item(4)%item(i)), &
		!	real(fn_scs_tmp%b_nm%item(3)%item(i)), imag(fn_scs_tmp%b_nm%item(3)%item(i)), &
		!	real(fn_scs_tmp%b_nm%item(4)%item(i)), imag(fn_scs_tmp%b_nm%item(4)%item(i))			
		!end do
		!close (206)
		
		open (unit=206, file = 'A_i_2_Interface_km14.txt') !A_i_2_km14
		do i = 1, K_m*K_m*2+2
			read (206, *) a_arr
			
			fn_scs_tmp%a_nm%item(1)%item(i)	 = cmplx(a_arr(1), a_arr(2))
			fn_scs_tmp%a_nm%item(2)%item(i)	 = cmplx(a_arr(3), a_arr(4))
			fn_scs_tmp%b_nm%item(1)%item(i)	 = cmplx(a_arr(5), a_arr(6))
			fn_scs_tmp%b_nm%item(2)%item(i)	 = cmplx(a_arr(7), a_arr(8))
			
			fn_scs_tmp%a_nm%item(3)%item(i)	 = cmplx(a_arr(1), a_arr(2)) !9, 10
			fn_scs_tmp%a_nm%item(4)%item(i)	 = cmplx(a_arr(3), a_arr(4))!11, 12
			fn_scs_tmp%b_nm%item(3)%item(i)	 = cmplx(a_arr(5), a_arr(6))!13, 14
			fn_scs_tmp%b_nm%item(4)%item(i)	 = cmplx(a_arr(7), a_arr(8))!15, 16					
		end do		
		close(206)		
				
		print*, 'fn_scs_tmp%a_nm%item(1)%item(1:2)', fn_scs_tmp%a_nm%item(1)%item(1:2)	!
		call Lagrange_Interpolation_New(pp, K_m, K_n, fn_scs_tmp, f_inter)		
		open (unit=206, file = 'Interpolation_finter_kn11_B.txt', action="write",status = 'replace')
		do i = 1, K_n*K_n*2+2
			write (206, '(16(es19.12, tr5))') 			&
			real(f_inter%a_nm%item(1)%item(i)), imag(f_inter%a_nm%item(1)%item(i)), &
			real(f_inter%a_nm%item(2)%item(i)), imag(f_inter%a_nm%item(2)%item(i)), &
			real(f_inter%b_nm%item(1)%item(i)), imag(f_inter%b_nm%item(1)%item(i)), &
			real(f_inter%b_nm%item(2)%item(i)), imag(f_inter%b_nm%item(2)%item(i)), &
			real(f_inter%a_nm%item(3)%item(i)), imag(f_inter%a_nm%item(3)%item(i)), &
			real(f_inter%a_nm%item(4)%item(i)), imag(f_inter%a_nm%item(4)%item(i)), &
			real(f_inter%b_nm%item(3)%item(i)), imag(f_inter%b_nm%item(3)%item(i)), &
			real(f_inter%b_nm%item(4)%item(i)), imag(f_inter%b_nm%item(4)%item(i))			
		end do		
		close(206)
	end subroutine test_lagrange_interpolation   
    
	end module lib_sie_tri_mlfmm