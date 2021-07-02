 module lib_sie_mlfmm_tri
	
	use Calculation_mod_tri
	use disc_mod
	use input_tri
	use libmath
	use omp_lib
	use Gauss_Quadrature_Fast
    use lib_sie_data_container

	implicit none
	! save 
	private 
    
    public :: Khat_and_TRF_ExtendedK 
    public :: TL_km_arr_II
    public :: FMM_aggregation_coefficient_element
    public :: Khat_and_TRF_II
    public :: init_list_sie
	public :: FMM_disaggregation_coefficient_RL
    public :: FMM_disaggregation_coefficient_RK
    public :: MLFMM_tri_Dabc_calculation
    public :: Lagrange_Interpolated
    public :: fxy_interpolation
    
    interface init_list_sie
        module procedure init_list_sie_c
        module procedure init_list_sie_v
    end interface
    
	contains
	
    function fxy_interpolation(Fn_scs) result(fxy)
        complex(kind = dp), dimension(:, :), allocatable, intent(in) :: Fn_scs
        complex(kind = dp), dimension(2, 4) :: fxy
        integer :: Nsp
        Nsp = size(Fn_scs(:, 1))
        fxy(1, 1:4) = Fn_scs(Nsp-1, 1:4)
        fxy(2, 1:4) = Fn_scs(Nsp, 1:4)
    end function
    
	subroutine Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter)
		implicit none
		integer, intent(in) :: pp, K_m, K_D
		integer :: n_cube, i, j, k, K_n, tt
		integer :: it, jp, m_phi, n_phi, t, s, kp, kt
		!type(Structure), intent(in) :: struct	
		!type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn_inter, Int_fn
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_inter
		complex(kind = dp), dimension(:, :), allocatable, intent(in) :: Fn_scs
		complex(kind = dp), dimension(:, :, :), allocatable :: ff_in, ff_out		
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		
	    K_n = K_m + K_D
		n_phi = 2*K_n
		m_phi = 2*K_m
		dphi = 2*PI/m_phi		
		allocate(ff_in(K_m+2, 2*K_m, 4), ff_out(K_m+2, 2*K_m, 4))
		allocate(ff_inter(n_phi*K_n, 4)) 
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_in(i, j, 1:4)= Fn_scs(k, 1:4) !theta in k1				
			end do
		end do
		do i = 1, 2		
			ff_in(K_m+1:K_m+2, 1, i) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i) 
			ff_in(K_m+1:K_m+2, 1, i+2) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i)		
			ff_in(K_m+1:K_m+2, 2, i) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i+2)		
			ff_in(K_m+1:K_m+2, 2, i+2) = Fn_scs(2*K_m*K_m+1:2*K_m*K_m+2, i+2)		
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
	end subroutine Lagrange_Interpolated
        
    subroutine Lagrange_Interpolation_Weights(x, t, p, tp_arr, tp_arr_ext, wk)
		implicit none
		!tp_arr: The original theta or phi array		
		integer, intent(in):: t, p!, K_n, K_m
		real(kind = dp), dimension(:), allocatable :: xx, wk, tp_arr, tp_arr_ext
		real(kind = dp) :: x
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
    
    subroutine Sorting_theta_phi_LenReduced_Arr_II(p, K_m, ff_in, ff_out)
		implicit none
		complex(kind = dp), dimension(:, :, :), allocatable, intent(out) :: ff_out
		complex(kind = dp), dimension(:, :, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m, p
		integer :: n_phi, m, n, K_n, j
		real(kind = dp), dimension(:), allocatable :: a_theta, a_phi
		real(kind = dp) :: Max_p(2), Max_t(2)
		real(kind = dp) :: xt, xp
		complex(kind = 8) :: fx1, fy1, fx2, fy2, fxy1(4, 2), fxy2(4, 2)
		intrinsic cos, sin
		
		K_n = 2*K_m
		n_phi = K_n+2*p+1
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		allocate(ff_out(K_m+2*p+2, K_n+2*p+1, 4))
		ff_out=(0.0, 0.0)
		
	  ! ++ fx and fy at theta = 0 and PI		
		
		do j = 1, 4
		fxy1(j, 1:2) =(/ff_in(K_m+1, 1, j), ff_in(K_m+1, 2, j)/)
		fxy2(j, 1:2) =(/ff_in(K_m+2, 1, j), ff_in(K_m+2, 2, j)/)
		end do
		
		! ++ f_theta and f_phi at theta = 0 and PI
		do m = 1, n_phi
		   n = m-p+K_n
			Max_t = (/cos(a_phi(n)), sin(a_phi(n))/)
			Max_p = (/-sin(a_phi(n)), cos(a_phi(n))/)
			do j = 1, 2
				ff_out(p+1, m, j) = -dot_product(Max_t, fxy2(j, :))  !
				ff_out(K_m+p+2, m, j) = dot_product(Max_t, fxy1(j, :))  !
				ff_out(p+1, m, j+2) = dot_product(Max_p, fxy2(j+2, :))  !
				ff_out(K_m+p+2, m, j+2) = dot_product(Max_p, fxy1(j+2, :) ) 
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
    end subroutine	Sorting_theta_phi_LenReduced_Arr_II
    
    subroutine theta_phi_extended(K_m, a_phi, a_theta)
	! theta and phi for Lagrange interpolation
	! p is the one-side sampling point number closing the source point. 
		real(kind = 8), dimension(:), allocatable :: tmp_theta, tmp_phi
		real(kind = 8), dimension(:), allocatable, intent(out) :: a_theta, a_phi
		integer, intent(in) :: K_m
		real(kind = 8) :: phi, d_phi, a, b
		integer :: Np, m_ph
		real(kind = 8), dimension(:), allocatable :: xx, ww
		integer :: j, k, i
		intrinsic :: sqrt, acos
				
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
        
    subroutine MLFMM_tri_Dabc_calculation(f_rk, f_rl, fD, D_abc)
        !sum_a = K1+K2
        !sum_b = L1+L2
        !sum_c = e1/u1*L1 + e2/u2*L2
        !Equation 28 in Mitchang's paper Vol. 11, P1383, JOSAA, 1994
        !K, L operators for FMM are Eq.9.25-9.33 in Gibson's book, p336
        complex(kind = dp), dimension(:, :), intent(in) :: f_rk, f_rl, fD
        complex(kind = dp), dimension(:, :), allocatable, intent(out) :: D_abc
        integer :: t, s, m
        
        m = size(f_rk(:, 1))
        !print*, 'size m=', m
        allocate(D_abc(3, m))
        
        D_abc(1:3, :) = (0.0, 0.0)
        do t = 1, 2
            do s = 1, 2
                D_abc(1, :) = D_abc(1, :) - (f_rk(:, t)*fD(:, s) + f_rk(:, t+2)*fD(:, s+2))
                D_abc(2, :) = D_abc(2, :) + (f_rl(:, t)*fD(:, s) + f_rl(:, t+2)*fD(:, s+2))
		        D_abc(3, :) = D_abc(3, :) + (eps_1/my_1*f_rl(:, t)*fD(:, s) + eps_2/my_2*f_rl(:, t+2)*fD(:, s+2))
            end do
        end do
        !sum_r(1) = sum_r(1) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs)
        !sum_r(2) = sum_r(2) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs)
    end subroutine
    
        
    subroutine init_list_sie_c(C, m, n)
        use libmlfmm
        type(lib_ml_fmm_v), intent(inout) :: C        
        integer, intent(in) :: m, n
        integer :: i
        if (m==1) then
            if (allocated(C%c)) then
                deallocate(C%c)
            end if
            allocate(C%C(1:n))
            C%C(1:n) = (0.0, 0.0)
        else
            if (allocated(C%a_nm%item)) then
                deallocate(C%a_nm%item)
            end if
            allocate(C%a_nm%item(1:m))
            do i = 1, m
                allocate(C%a_nm%item(i)%item(1:n))
                C%a_nm%item(i)%item(1:n) = (0.0, 0.0)
            end do            
        end if
    end subroutine
    
    subroutine init_list_sie_v(C, m, n)
        use libmlfmm
        type(lib_ml_fmm_coefficient), intent(inout) :: C        
        integer, intent(in) :: m, n
        integer :: i
        if (m==1) then
            if (allocated(C%c)) then
                deallocate(C%c)
            end if
            allocate(C%C(1:n))
            C%C(1:n) = (0.0, 0.0)
        else
            if (allocated(C%a_nm%item)) then
                deallocate(C%a_nm%item)
            end if
            allocate(C%a_nm%item(1:m))
            do i = 1, m
                allocate(C%a_nm%item(i)%item(1:n))
                C%a_nm%item(i)%item(1:n) = (0.0, 0.0)
            end do            
        end if
    end subroutine
    
    subroutine Khat_and_TRF_ExtendedK(K_m)!, K_hat_p, transfer_theta, transfer_phi, TRF_p
		implicit none	
		real(kind = dp), dimension(:), allocatable :: xx, ww!, x        
		real(kind = dp) :: d_phi, phi, a, b, sum_a, dzero
        
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
        allocate(TRF_p(N_sp + 2))
		allocate(transfer_theta(N_sp + 2, 3), transfer_phi(N_sp + 2, 3))
		allocate(xx(K_m))
		allocate(ww(K_m))
		a = -1.0
		b = 1.0		
		call legendre_handle(K_m, a, b, xx, ww)			
		
		d_phi = 2*PI/(2*K_m)
		do j = 1, K_m
            do k = 1, 2*K_m
            phi = (k-1)*d_phi 
            i = k + (j-1)*2*K_m
            k_hat_p(i, :) = (/sin(acos(xx(j)))*cos(phi), sin(acos(xx(j)))*sin(phi), xx(j)/)
            transfer_theta(i, 1:3) = (/xx(j)*cos(phi), xx(j)*sin(phi), -sin(acos(xx(j)))/)
		    transfer_phi(i, 1:3) = (/-sin(phi), cos(phi), dzero/)
		    TRF_p(i) = ww(j)*d_phi				
            end do
		end do
		k_hat_p(N_sp+1, :)= (/0.0, 0.0, 1.0/)
		k_hat_p(N_sp+2, :)= (/0.0, 0.0, -1.0/)
		transfer_theta(N_sp+1, 1:3) = (/1.0, 0.0, 0.0/) !For saving fx at theta = 0
		transfer_theta(N_sp+2, 1:3) = (/1.0, 0.0, 0.0/)!For saving fx at theta = pi
		transfer_phi(N_sp+1, 1:3) = (/0.0, 1.0, 0.0/)
		transfer_phi(N_sp+2, 1:3) = (/0.0, 1.0, 0.0/) !For saving fy at pi
		TRF_p(N_sp+1:N_sp+2) = 0.0
		deallocate(xx)
		deallocate(ww)
    end subroutine Khat_and_TRF_ExtendedK

    subroutine FMM_aggregation_coefficient_element(K_m, m, struct, R_c, ff_out) !
		implicit none		
		integer, intent(in) :: m, K_m		
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3) !Center of the box
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_out
        !dummy
        complex(kind = dp), dimension(:, :), allocatable :: Fn_scs
		integer :: tt
		integer :: Ns
		
		Ns = 2*K_m*K_m
		allocate(ff_out(Ns+2, 4))			
		ff_out(Ns+2, 1:4) = (0.0, 0.0)
		do tt = 1, 2 !
			call RT_function_SCSArr(K_m, struct, m, tt, R_c, Fn_scs)
			Fn_scs = conjg(Fn_scs)
			ff_out(1:Ns+2, 1:4) = ff_out(1:Ns+2, 1:4) + Fn_scs(1:Ns+2, 1:4)
        end do	
    end subroutine FMM_aggregation_coefficient_element
    
    subroutine RT_function_SCSArr(K_m, struct, n, ll, Rc, Fn_scs)
		implicit none		
		real(kind = dp), intent(in) :: Rc(3)		
		integer, intent(in) :: n, K_m  !sampling point on a sphere				
		integer :: i, j, pos_neg, ll, m_pairs, Ns
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roc(3), Ro(3), fm(3)
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)!
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: Fn_scs
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		type(Integration_fn) :: Int_fn
        
		m_pairs = size(struct%neighbours)
		call Quadrature_tri(ng, a, b, w)	
		call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
		pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
		dfm = pos_neg*Len_m 			
		Ns = K_m*2*K_m
		allocate(Int_fn%pair_kk(ll)%Rfk_vec(Ns+2, 2, 3))		
		allocate(Fn_scs(Ns+2, 4))
		do i = 1, Ns+2
			Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			do j = 1, ng
				Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
				fm = 0.5*(Ro-vm_t)*dfm
				Roc = Ro - Rc
				tmp = dot_product(k_hat_p(i, :), Roc(:))
				ff1 = exp(-im*k1*tmp)*w(j)
				ff2 = exp(-im*k2*tmp)*w(j)
				Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) + fm*ff1
				Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) + fm*ff2					
			end do
			Fn_scs(i, 1) = dot_product(transfer_theta(i, 1:3), Int_fn%pair_kk(ll)%Rfk_vec(i, 1, 1:3))	!k1, fx1						
			Fn_scs(i, 2) = dot_product(transfer_theta(i, 1:3), Int_fn%pair_kk(ll)%Rfk_vec(i, 2, 1:3)) !k2, fx2
			Fn_scs(i, 3) = dot_product(transfer_phi(i, 1:3), Int_fn%pair_kk(ll)%Rfk_vec(i, 1, 1:3))  !k1, fy1								
			Fn_scs(i, 4) = dot_product(transfer_phi(i, 1:3), Int_fn%pair_kk(ll)%Rfk_vec(i, 2, 1:3))  !k2, fy2
		end do !			
    end subroutine RT_function_SCSArr
	
    
    !Put m in the function
    function TL_km_arr_III(N, r_ab, K_m, kd)
        integer, intent(in) :: K_m, N    !m is the order limit of the polynomials, N is the sampling points on the sphere
        real(kind = dp), intent(in) :: r_ab(3)
        complex(kind =dp), intent(in) :: kd
        
        real(kind = dp) :: r_hat(3), x
        real(kind = dp), dimension(:, :), allocatable :: pm_arr        
        real(kind = dp), dimension(:, :), allocatable :: k_hat_II
        real(kind = dp), dimension(:), allocatable :: TRF_II
        real (kind = dp), dimension(:), allocatable :: pm, dummy 
        
        complex(kind = dp) :: r_h
        complex(kind = dp), dimension(:), allocatable :: TL_km_arr_III
        complex(kind = dp), dimension(:), allocatable :: hl 
        integer :: j, fnu, Lm
		
        fnu = 0        
        r_hat = r_ab/vec_len(r_ab)
        r_h =  vec_len(r_ab)*Kd  
        
        Lm = abs(r_h) + 10*(abs(r_h))**(1/3)
        !if (Lm > int(abs(r_h))) then
        !    Lm = int(abs(r_h))
        !end if
        
        allocate(hl(Lm), pm(Lm), dummy(Lm)) 
        hl =  lib_math_hankel_spherical_2(r_h, fnu, Lm) 
        
        allocate(pm_arr(2*K_m*K_m, Lm), TL_km_arr_III(2*K_m*K_m))
        
        TL_km_arr_III(1:2*K_m*K_m) = (0.0, 0.0)   
        call Khat_and_TRF_II(K_m, K_hat_II, TRF_II)        
        do j = 1, 2*K_m*K_m
            x = dot_product(r_hat(:), k_hat_II(j, :))
            call lib_math_legendre_polynomial(x, fnu, Lm, pm, dummy)
            pm_arr(j, :) = pm(:)
        end do
		
        do j = 1, Lm
            TL_km_arr_III(:) = TL_km_arr_III(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd/(4*PI)**2 ! 
        end do
        return        
    end function TL_km_arr_III
    
    function TL_km_arr_II(r_ab, m, K_m, kd)
      integer, intent(in) :: K_m, m    !Lm is the order limit of the polynomials, N is the sampling points on the sphere
      real(kind = dp), intent(in) :: r_ab(3)
		complex(kind =dp), intent(in) :: kd
      real(kind = dp) :: r_hat(3), x
		real(kind = dp), dimension(:, :), allocatable :: pm_arr
      complex(kind = dp) :: r_h
		real(kind = dp), dimension(:, :), allocatable :: k_hat_II
		real(kind = dp), dimension(:), allocatable :: TRF_II
      real (kind = dp), dimension(m) :: pm, dummy 
      complex(kind = dp), dimension(:), allocatable :: TL_km_arr_II
      complex(kind = dp), dimension(m) :: hl 
      integer :: j, fnu
		
      fnu = 0        
      r_hat = r_ab/vec_len(r_ab)
      r_h =  vec_len(r_ab)*Kd  
      hl =  lib_math_hankel_spherical_2(r_h, fnu, m) 

		allocate(pm_arr(2*K_m*K_m, m), TL_km_arr_II(2*K_m*K_m))
		TL_km_arr_II(1:2*K_m*K_m) = (0.0, 0.0)   
		call Khat_and_TRF_II(K_m, K_hat_II, TRF_II)
      do j = 1, 2*K_m*K_m
         x = dot_product(r_hat(:), k_hat_II(j, :))
         call lib_math_legendre_polynomial(x, fnu, m, pm, dummy)
         pm_arr(j, :) = pm(:)
      end do
		
      do j = 1, m
         TL_km_arr_II(:) = TL_km_arr_II(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd/(4*PI)**2 ! 
		end do
      return        
	end function TL_km_arr_II
    
    subroutine Khat_and_TRF_II(K_m, K_hat_II, TRF_II)
		implicit none	
		integer, intent(in) :: K_m
		real(kind = dp), dimension(:), allocatable :: xx, ww 
		real(kind = dp), dimension(:), allocatable, intent(out) :: TRF_II
		real(kind = dp), dimension(:, :), allocatable, intent(out) :: k_hat_II
		real(kind = dp) :: d_phi, phi, a, b, sum_a
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
    
    subroutine FMM_disaggregation_coefficient_RK(K_m, m, struct, R_c, d_c, ff_out) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m		
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3)
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, Fn_scs, ff_out, ff_tmp
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		
		integer :: K_m, K_D, K_n, pp, tt
		integer :: Nsp		
		
		Nsp = K_m*2*K_m
        
		allocate(ff_tmp(K_m*K_m*2+2, 4))
		if (allocated(ff_out)) then
            deallocate(ff_out)
        end if
        allocate(ff_out(K_m*K_m*2, 4))
		ff_out(K_m*K_m*2, 4) = (0.0, 0.0)
		do tt = 1, 2
			call RT_function_SCSArr(K_m, struct, m, tt, R_c, Fn_scs)                
			ff_tmp(:, 1) = -Fn_scs(:, 3)
			ff_tmp(:, 3) = Fn_scs(:, 1)
			ff_tmp(:, 2) = -Fn_scs(:, 4)
			ff_tmp(:, 4) = Fn_scs(:, 2)                
			ff_out(1:K_m*K_m*2, 4) = ff_out(1:K_m*K_m*2, 4) + ff_tmp(1:K_m*K_m*2, 4)
		end do
        ff_out(:, 1:2) = ff_out(:, 1:2)*im*k1
        ff_out(:, 3:4) = ff_out(:, 3:4)*im*k2
    end subroutine FMM_disaggregation_coefficient_RK
    
    subroutine FMM_disaggregation_coefficient_RL(K_m, m, struct, R_c, d_c, ff_out) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3)
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, Fn_scs, ff_tmp
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_out
        
		integer :: K_m, K_D, K_n, pp, tt
		integer :: Ns
		Ns = K_m*2*K_m
		
		allocate(ff_tmp(Ns+2, 4))
		allocate(ff_out(Ns, 4))
        
		ff_out(Ns, 4) = (0.0, 0.0)
		do tt = 1, 2
			call RT_function_SCSArr(K_m, struct, m, tt, R_c, Fn_scs)				
			ff_out(1:Ns, 4) = ff_out(1:Ns, 4) + Fn_scs(1:Ns, 4)
        end do            
        
        ff_out(:, 1:2) = ff_out(:, 1:2)*im*Omega*my_1
        ff_out(:, 3:4) = ff_out(:, 3:4)*im*Omega*my_2
        !ff_out = ff_out*im*omega*my_0 !
    end subroutine FMM_disaggregation_coefficient_RL
    
 !   subroutine FMM_Near_MVProduct_Tri(m_pairs, struct, struct_cube, d_c, aa_csr, ib, jb, D_mat)
	!    !Near field matrix in a csr format
	!	integer :: m_cube, n_cube, m, n_edge_number, nn_edge, mm_edge, nd
	!	integer :: m_pairs, n_cube_number, m_cube_number, m_edge_number
	!	type(Structure), intent(in) :: struct
	!	type(Structure_cube), intent(in) :: struct_cube
	!	complex(kind = dp) :: sum_a, sum_c, sum_b		
	!	complex(kind = dp) :: ssum_a, ssum_c, ssum_b				
	!	integer, dimension(:), allocatable :: n_arr_tmp
	!	complex(kind = dp), dimension(:), allocatable :: a_tmp, a_csr, c_csr, b_csr
	!	complex(kind = dp), dimension(:), allocatable, intent(out) :: aa_csr
	!	complex(kind = dp), dimension(:, :), allocatable, intent(out) :: D_mat
	!	integer, dimension(:), allocatable :: ia, ja
	!	integer, dimension(:), allocatable, intent(out) :: ib, jb
	!	
	!	real(kind = dp), intent(in) :: d_c ! 
	!	real(kind = dp) :: r_ob, r_ref
	!	integer :: p, q, ngp, t, s, n, j, n_csr, m_csr!, counter 
	!	
	!	!---------Calculate the number of non-zero elements in the array a_st---------------		
	!	!counter = 0
	!	m_cube_number = size(struct_cube%cubes) 
	!	!print*, 'm_cube_number =', m_cube_number
	!	n_csr = 0
	!	do m = 1, m_cube_number
	!		m_edge_number = size(struct_cube%cubes(m)%edges_in_cube)
	!		n_cube_number = size(struct_cube%cubes(m)%cubes_near) !number of cubes in the near region of cube i
	!		do mm_edge = 1, m_edge_number	
	!			do n = 1, n_cube_number
	!				n_cube = struct_cube%cubes(m)%cubes_near(n)
	!				n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube)
	!				n_csr = n_csr + n_edge_number
	!			end do			
	!		end do
	!	end do
	!	allocate(a_csr(n_csr))
	!	allocate(b_csr(n_csr))
	!	allocate(c_csr(n_csr))
	!	allocate(a_tmp(n_csr))
	!	allocate(ja(n_csr))	
	!	allocate(ia(m_pairs + 1))
	!	allocate(D_mat(2*m_pairs, 2*m_pairs))
	!	allocate(ib(2*m_pairs + 1))
	!	allocate(jb(4*n_csr))
	!	allocate(aa_csr(4*n_csr))	
	!	allocate(n_arr_tmp(2*m_pairs))
	!	!
	!	ia(1:m_pairs+1) = 0
	!	ja(1:n_csr) = 0
	!	n_csr = 0
	!	m_csr = 0
	!	ia(1) = 1
	!	D_mat(2*m_pairs, 2*m_pairs) = (0.0, 0.0)
	!	do m = 1, m_pairs	
	!		n_csr = 0
	!		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
	!		n_cube_number = size(struct_cube%cubes(m_cube)%cubes_near) 		
	!		!print*, 'near n_cube_number =', n_cube_number
	!		do n = 1, n_cube_number
	!			n_cube = struct_cube%cubes(m_cube)%cubes_near(n) !cube index in the near 
	!			n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube
	!			do nn_edge = 1, n_edge_number
	!				n_csr = n_csr + 1
	!				n_arr_tmp(n_csr) = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n					
	!			end do
	!		end do
	!		call sort_edge_index(n_arr_tmp(1:n_csr), n_csr, n_arr_tmp(1:n_csr))		
	!		ja(m_csr + 1 : m_csr + n_csr) = n_arr_tmp(1: n_csr)		
	!	
	!		p = struct%neighbours(m)%elements(1)
	!		q = struct%neighbours(m)%elements(2)	
	!		r_ref = vec_len(struct%midpoints(p)%point - struct%midpoints(q)%point)		
	!		do j = 1, n_csr
	!			n = n_arr_tmp(j)
	!			
	!			sum_a = (0.0, 0.0) 
	!			sum_b = (0.0, 0.0) 
	!			sum_c = (0.0, 0.0) 
	!			do t = 1, 2
	!				do s = 1, 2!
	!					p = struct%neighbours(m)%elements(t)					
	!					q = struct%neighbours(n)%elements(s)						
	!					r_ob = vec_len(struct%midpoints(p)%point - struct%midpoints(q)%point)					
	!					if (r_ob > r_ref) then
	!						ngp = 3				
	!						call Normal_Integration_new(m, n, t, s, ngp, struct, ssum_a, ssum_b, ssum_c)								
	!					else
	!						ngp = 3	
	!						call Singular_Integration(m, n, t, s, ngp, struct, ssum_a, ssum_b, ssum_c)
	!					end if
	!					sum_a = sum_a + ssum_a !
	!					sum_b = sum_b + ssum_b
	!					sum_c = sum_c + ssum_c
	!				end do !loop s 
	!			end do  !loop t
	!			a_csr(m_csr + j) = sum_a !one-based indexing
	!			b_csr(m_csr + j) = sum_b
	!			c_csr(m_csr + j) = sum_c
	!			D_mat(m, n) = sum_b
	!			D_mat(m, n + m_pairs) = -sum_a
	!			D_mat(m + m_pairs, n) = sum_a
	!			D_mat(m + m_pairs, n + m_pairs) = sum_c
	!		end do
	!		m_csr = m_csr + n_csr
	!		ia(m+1) = ia(m) + n_csr
	!		n_csr = 0
	!	end do !loop m		
	!	call Expension_A_CSR(a_csr, b_csr, c_csr, ia, ja, aa_csr, ib, jb)
	!end subroutine FMM_Near_MVProduct_Tri
 !   
	end module lib_sie_mlfmm_tri