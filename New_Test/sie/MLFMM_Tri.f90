 module MLFMM_Tri
	
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
	!public :: cube_definition		
	public :: FMM_Near_MVProduct_Tri	
    public :: FMM_Downsampled_RK
    public :: FMM_Downsampled_RL
	public :: FMM_edge_Tri_RT	
    public :: TL_km_arr_II
	public :: FMM_aggregation_coefficient_element
	public :: FMM_Aggregation
    public :: calculate_km
	 
    contains
    
    	subroutine Khat_and_TRF_ExtendedK(K_m)!, K_hat_p, transfer_theta, transfer_phi, TRF_p
		implicit none	
		real(kind = dp), dimension(:), allocatable :: xx, ww!, x        
		real(kind = dp) :: d_phi, phi, a, b, sum_a, dzero
        
		integer :: j, k, i, N_sp
        integer, intent(in) :: K_m
	    dzero = 0.0		
		N_sp = K_m*2*K_m
        
		allocate(k_hat_p(N_sp + 2, 3))
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
		
		!open (unit=206, file = 'Transfer_theta.txt', action="write",status = 'replace')
		!do i = 1, Km*(2*Km+1)
		!	write (206, '(3(es19.12, tr5))')  transfer_theta(i, 1:2)!, x(i)
		!end do
		!close(206)
		!
		!open (unit=206, file = 'Transfer_phi.txt', action="write",status = 'replace')
		!do i = 1, Km*(2*Km+1)
		!	write (206, '(2(es19.12, tr5))')  transfer_phi(i, 1:2)
		!end do		
		!close(206)
		
		deallocate(xx)
		deallocate(ww)
	end subroutine Khat_and_TRF_ExtendedK
		
	subroutine FMM_aggregation_coefficient_element(K_m, m, struct, R_c, ff_out) !
		implicit none		
		integer, intent(in) :: m, K_m		
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3) !Center of the box
		!complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: Fn_scs, ff_out
		!real(kind = dp), intent(in) :: d_c ! sphere diameter of the highest level				
		integer :: tt
		integer :: Ns
		
		Ns = K_m*2*K_m
		allocate(ff_out(Ns, 4))			
		ff_out(Ns, 4) = (0.0, 0.0)
		do tt = 1, 2 !
			call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)
			Fn_scs = conjg(Fn_scs)
			ff_out(1:K_m*K_m*2, 4) = ff_out(1:K_m*K_m*2, 4) + Fn_scs(1:K_m*K_m*2, 4)
		end do		
    end subroutine FMM_aggregation_coefficient_element
    
	!nd : edge index
	subroutine FMM_Aggregation(interpol, m, struct, R_c, d_c, ff_out) !
		implicit none		
		integer, intent(in) :: m
		logical, intent(in) :: interpol
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3) !Center of the box
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, Fn_scs, ff_out
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the highest level				
		integer :: K_m, K_D, K_n, pp, tt
		integer :: Nsp
		K_m = 12
		K_D = + 3
		pp = 4
		K_n = K_m + K_D
		Nsp = K_n*2*K_n
		if (interpol) then
			allocate(ff_out(Nsp, 4))
			ff_out(1:Nsp, 4) = (0.0, 0.0)
			do tt = 1, 2
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)
				Fn_scs = conjg(Fn_scs)
				call Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter_R)				
				ff_out = ff_out + ff_inter_R
			end do
		else 
			allocate(ff_out(K_m*K_m*2, 4))			
			ff_out(K_m*K_m*2, 4) = (0.0, 0.0)
			do tt = 1, 2 !
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)
				Fn_scs = conjg(Fn_scs)
				ff_out(1:K_m*K_m*2, 4) = ff_out(1:K_m*K_m*2, 4) + Fn_scs(1:K_m*K_m*2, 4)
			end do
		end if		
    end subroutine FMM_Aggregation
    
    function calculate_km(d_c) result(k_m)
        real(kind = dp), intent(in) :: d_c
        integer :: k_m
        k_m = abs(k2*d_c)+1.8*(d0)**(2/3)*(abs(k2)*d_c)**(1/3)
    end function
	
	subroutine FMM_sie_translation(k_m, x_start, x_end, d_c, TL_1, TL_2)
		real(kind = dp) :: x_start(3), x_end(3), Rab(3)
		real(kind = dp) :: d_c
		complex(kind = dp), dimension(:), allocatable, intent(out) :: TL_1, TL_2		
		integer :: K_n, K_D, Nsp, Lm
        integer, intent(in) :: k_m

			K_D = 3
			K_n = K_m + K_D
			Nsp = 2*K_n*K_n
			allocate(TL_1(Nsp), TL_2(Nsp))
			Rab = x_end -x_start !target - source
			Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)
			TL_1 = TL_km_arr_II(Rab, Lm, K_n, k1)
			Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
			TL_2 = TL_km_arr_II(Rab, Lm, K_n, k2)		
	end subroutine FMM_sie_translation	
	
	subroutine FMM_Downsampled_RK(interpol, m, struct, R_c, d_c, ff_out) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m
		logical, intent(in) :: interpol
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3)
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, Fn_scs, ff_out, ff_tmp
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		
		integer :: K_m, K_D, K_n, pp, tt
		integer :: Nsp		
		K_m = 12
		K_D = -3
		pp = 4
		K_n = K_m + K_D
		Nsp = K_n*2*K_n
		
		allocate(ff_tmp(K_m*K_m*2+2, 4))
		
		if (interpol) then
			allocate(ff_out(Nsp, 4))
			ff_out(1:Nsp, 4) = (0.0, 0.0)
			do tt = 1, 2
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)	
				ff_tmp(:, 1) = -Fn_scs(:, 3)
				ff_tmp(:, 3) = Fn_scs(:, 1)
				ff_tmp(:, 2) = -Fn_scs(:, 4)
				ff_tmp(:, 4) = Fn_scs(:, 2)
				call Lagrange_Interpolated(pp, K_m, K_D, ff_tmp, ff_inter_R)
				ff_out = ff_out + ff_inter_R
			end do
		else 
			allocate(ff_out(K_m*K_m*2, 4))			
			ff_out(K_m*K_m*2, 4) = (0.0, 0.0)
			do tt = 1, 2
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)                
				ff_tmp(:, 1) = -Fn_scs(:, 3)
				ff_tmp(:, 3) = Fn_scs(:, 1)
				ff_tmp(:, 2) = -Fn_scs(:, 4)
				ff_tmp(:, 4) = Fn_scs(:, 2)                
				ff_out(1:K_m*K_m*2, 4) = ff_out(1:K_m*K_m*2, 4) + ff_tmp(1:K_m*K_m*2, 4)
			end do
        end if
        ff_out(:, 1:2) = ff_out(:, 1:2)*im*k1
        ff_out(:, 3:4) = ff_out(:, 3:4)*im*k2
    end subroutine FMM_Downsampled_RK
    
    subroutine FMM_Downsampled_RL(interpol, m, struct, R_c, d_c, ff_out) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m
		logical, intent(in) :: interpol
		type(Structure), intent(in) :: struct
		real(kind = dp), intent(in) ::  R_c(3)
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, Fn_scs, ff_out, ff_tmp
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		
		integer :: K_m, K_D, K_n, pp, tt
		integer :: Nsp		
		K_m = 12
		K_D = -3
		pp = 4
		K_n = K_m + K_D
		Nsp = K_n*2*K_n
		
		allocate(ff_tmp(K_m*K_m*2+2, 4))
		
		if (interpol) then
			allocate(ff_out(Nsp, 4))
			ff_out(1:Nsp, 4) = (0.0, 0.0)
			do tt = 1, 2
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)	
				call Lagrange_Interpolated(pp, K_m, K_D, fn_scs, ff_inter_R)
				ff_out = ff_out + ff_inter_R
			end do
		else 
			allocate(ff_out(K_m*K_m*2, 4))			
			ff_out(K_m*K_m*2, 4) = (0.0, 0.0)
			do tt = 1, 2
				call RT_function_SCSArr(struct, m, tt, K_m, R_c, Fn_scs)				
				ff_out(1:K_m*K_m*2, 4) = ff_out(1:K_m*K_m*2, 4) + Fn_scs(1:K_m*K_m*2, 4)
            end do            
        end if
        ff_out(:, 1:2) = ff_out(:, 1:2)*im*Omega*my_1
        ff_out(:, 3:4) = ff_out(:, 3:4)*im*Omega*my_2
        !ff_out = ff_out*im*omega*my_0 !
	end subroutine FMM_Downsampled_RL

	! K_D: increased or reduced interpolation points
	subroutine FMM_edge_RT_Lagrange(m, n, struct, struct_cube, d_c, suma, sumb, sumc) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), dimension(:), allocatable :: sum_a, sum_c, sum_b!
		complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2
		complex(kind = dp), dimension(:), allocatable :: ss1, ss2, tmp
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter_R, ff_inter_T, Fn_scs
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3), dphi 
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
		real(kind = dp), dimension(:), allocatable :: TRF_II
		real(kind = dp), dimension(:, :), allocatable :: k_hat_II
		type(Integration_fn) :: Int_fn, Int_fm
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll, K_m, K_D, K_n, pp
		integer :: Nsp
		intrinsic :: sqrt	
		K_m = 11
		K_D = +3
		pp = 4
		K_n = K_m + K_D		
	
		Nsp = K_n*2*K_n		
		allocate(sum_a(Nsp), tmp(Nsp), sum_b(Nsp), sum_c(Nsp))
		allocate(TL_1(Nsp), TL_2(Nsp), ss1(Nsp), ss2(Nsp))
		allocate(TRF_II(Nsp))
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position !target
		Rb = struct_cube%cubes(n_cube)%cube_position !source	
		Rab = Ra - Rb
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)
		TL_1 = TL_km_arr_II(Rab, Lm, K_n, k1)
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr_II(Rab, Lm, K_n, k2)
		
		call Khat_and_TRF_II(K_n, K_hat_II, TRF_II)
		
		sum_a(1:Nsp) = (0.0, 0.0)
		sum_b(1:Nsp) = (0.0, 0.0)
		sum_c(1:Nsp) = (0.0, 0.0)
		do tt = 1, 2 !
			call RT_function_SCSArr(struct, m, tt, K_m, Ra, Fn_scs)
			call Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter_R)			
			do ll = 1, 2
				call RT_function_SCSArr(struct, n, ll, K_m, Rb, Fn_scs)
				call Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter_T)
				ff_inter_T = conjg(ff_inter_T)
				ss1 = -ff_inter_R(:, 3)*ff_inter_T(:, 1) + ff_inter_R(:, 1)*ff_inter_T(:, 3)
				ss2 = -ff_inter_R(:, 4)*ff_inter_T(:, 2) + ff_inter_R(:, 2)*ff_inter_T(:, 4)	
				
				sum_a = sum_a - im*(TL_1*ss1*k1+TL_2*ss2*k2)*TRF_II!
				ss1 = ff_inter_R(:, 1)*ff_inter_T(:, 1) + ff_inter_R(:, 3)*ff_inter_T(:, 3)
				ss2 = ff_inter_R(:, 2)*ff_inter_T(:, 2) + ff_inter_R(:, 4)*ff_inter_T(:, 4)		
				sum_b = sum_b + im*omega*(my_1*TL_1*ss1 + my_2*TL_2*ss2)*TRF_II
				sum_c = sum_c + im*omega*(eps_1*TL_1*ss1 + eps_2*TL_2*ss2)*TRF_II 											
			end do
		end do		
		suma = sum(sum_a)
		sumb = sum(sum_b)
		sumc = sum(sum_c)
	end subroutine FMM_edge_RT_Lagrange
	
	function TL_km_arr(N_sp, r_ab, m, kd)
      integer, intent(in) :: m, N_sp    !Lm is the order limit of the polynomials, N is the sampling points on the sphere
      real(kind = dp), intent(in) :: r_ab(3)
      real(kind = dp) :: r_hat(3), x, pm_arr(N_sp, m)      
      complex(kind = dp) :: r_h
		!real(kind = dp) :: k_hat(N_sp, 3)
      complex(kind = dp), intent(in) :: kd !complex wavevector
      real (kind = dp), dimension(m) :: pm, dummy 
      complex(kind = dp) :: TL_km_arr(N_sp)
      complex(kind = dp), dimension(m) :: hl 
      integer :: j, fnu
      fnu = 0        
      r_hat = r_ab/vec_len(r_ab)
      r_h =  vec_len(r_ab)*kd  
      hl =  lib_math_hankel_spherical_2(r_h, fnu, m) 
      TL_km_arr(1:N_sp) = (0.0, 0.0) 
      
      do j = 1, N_sp
         x = dot_product(r_hat(:), k_hat(j, :))
         call lib_math_legendre_polynomial(x, fnu, m, pm, dummy)
         pm_arr(j, :) = pm(:)
      end do
      do j = 1, m
         TL_km_arr(:) = TL_km_arr(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*kd/(4*PI)**2 ! 
		end do
      return        
	end function TL_km_arr
	
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
	! 
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
		complex(kind = dp), dimension(:, :), allocatable :: Fn_scs
		complex(kind = dp), dimension(:, :, :), allocatable :: ff_in, ff_out		
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		
		!nd = 11		
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
		!allocate(Int_fn_inter%Tri(ll)%Rfk_scs(n_phi*K_n, 2, 2)) 		
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
	
	subroutine Test_Lagrange_Interpolated_LenReduced_Arr_II(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n, ll, tt
		integer :: it, jp, m_phi, n_phi, t, s, pp, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn_inter, Int_fn
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter, Fn_scs
		complex(kind = dp), dimension(:, :, :), allocatable :: ff_in, ff_out		
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		
		nd = 11
		ll = 1 !index of the triangle in the pair
		pp = 4 !points for interpolation !p should be smaller than Km
		K_m = 12
		K_n = K_m - 2
		n_phi = 2*K_n
		
		m_phi = 2*K_m
		dphi = 2*PI/m_phi		
		allocate(ff_in(K_m+2, 2*K_m, 4), ff_out(K_m+2, 2*K_m, 4), ff_inter(n_phi*K_n, 4)) 
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCSArr(struct, nd, ll, K_m, Ra, Fn_scs)	
		
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
		allocate(Int_fn_inter%Tri(ll)%Rfk_scs(n_phi*K_n, 2, 2)) 		
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
		!Int_fn_inter%Tri(ll)%Rfk_scs(:, 1, 1) = ff_inter(:, 1)
		!Int_fn_inter%Tri(ll)%Rfk_scs(:, 2, 1) = ff_inter(:, 2)
		!Int_fn_inter%Tri(ll)%Rfk_scs(:, 1, 2) = ff_inter(:, 3)
		!Int_fn_inter%Tri(ll)%Rfk_scs(:, 2, 2) = ff_inter(:, 4)
		!
		!open (unit=206, file = 'ff_inter_theta1_real_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta2_real_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi1_real_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi2_real_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta1_imag_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta2_imag_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi1_imag_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi2_imag_Km8_Arr.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi1_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi2_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta1_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta2_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi1_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi2_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta1_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta2_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		!end do
		!close(206)	
		!
	end subroutine Test_Lagrange_Interpolated_LenReduced_Arr_II
	
	subroutine RT_function_SCSArr(struct, n, ll, K_m, Rc, Fn_scs)

		implicit none		
		real(kind = dp), intent(in) :: Rc(3)		
		integer, intent(in) :: n, K_m  !sampling point on a sphere				
		integer :: i, j, pos_neg, ll, m_pairs, Nsp
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roc(3), Ro(3), fm(3)
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)!
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: Fn_scs
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		type(Integration_fn) :: Int_fn
        !real(kind = dp), dimension(:, :), allocatable :: K_hat_p
        !real(kind = dp), dimension(:), allocatable :: TRF_p
        !real(kind = dp), dimension(:, :), allocatable :: transfer_theta, transfer_phi
        
        !call Khat_and_TRF_ExtendedK(K_m, K_hat_p, transfer_theta, transfer_phi, TRF_p)
		call Khat_and_TRF_ExtendedK(K_m)
        
		m_pairs = size(struct%neighbours)
		call Quadrature_tri(ng, a, b, w)	
		call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
		pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
		dfm = pos_neg*Len_m 			
		Nsp = K_m*2*K_m
		allocate(Int_fn%pair_kk(ll)%Rfk_vec(Nsp+2, 2, 3))		
		allocate(Fn_scs(Nsp+2, 4))
		do i = 1, Nsp+2
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
	
	subroutine Sorting_theta_phi(Component, K_m, ff_in, ff_out)
		implicit none
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_out
		complex(kind = dp), dimension(:, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m
		integer :: n_phi, m, n, K_n
		real(kind = 8), dimension(:), allocatable :: a_theta, a_phi
		real(kind = dp) :: xt, xp
		complex(kind = 8) :: fx1, fy1, fx2, fy2
		character(len = 1) :: Component
		!character(len = 50), public :: illumination
		intrinsic cos, sin
		
		fx1 = ff_in(K_m+1, 1)
		fy1 = ff_in(K_m+1, 2)
		fx2 = ff_in(K_m+2, 1)
		fy2 = ff_in(K_m+2, 2)
		
		K_n = 2*K_m
		n_phi = 3*K_n
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		allocate(ff_out(3*K_m+2, 3*K_n))
		ff_out=(0.0, 0.0)
		
		if (Component =='P') then		
			do m = 1, n_phi
				ff_out(K_m+1, m) = -sin(a_phi(m))*fx2 +cos(a_phi(m))*fy2  !theta = PI in the new 		
				ff_out(2*K_m+2, m) = -sin(a_phi(m))*fx1 +cos(a_phi(m))*fy1 !theta = 0 in the new array
			end do
		else if (Component =='T') then
			do m = 1, n_phi
				ff_out(K_m+1, m) = -cos(a_phi(m))*fx2 -sin(a_phi(m))*fy2  !theta = PI, here the problem with fx, and fy		
				ff_out(2*K_m+2, m) = cos(a_phi(m))*fx1+sin(a_phi(m))*fy1
			end do		
		else 
			print*, 'Wrong input component'
		end if		
		ff_out(K_m+2:2*K_m+1, K_n+1:2*K_n) = ff_in(1:K_m, 1:2*K_m)
		
		do m = 1, 3*K_m+2
			xt = a_theta(m)
			do n = 1, n_phi
				xp = a_phi(n)
				if (xt<0.0)then
					if (xp<PI .and. xp>=0.0) then
						ff_out(m, n) = -ff_out(4*K_m+4-m, n+K_n/2)					
					elseif (xp<2*PI .and. xp>=PI) then
						ff_out(m, n) = -ff_out(4*K_m+4-m, n-K_n/2)
					endif
				elseif (xt>PI) then
					if (xp<PI .and. xp>=0.0)then
						ff_out(m, n) = -ff_out(2*(K_m+1)-m, n+K_n/2)
					elseif (xp<2*PI .and. xp>=PI)then
						ff_out(m, n) = -ff_out(2*(K_m+1)-m, n-K_n/2)
					end if 
				elseif (xt<PI .and. xt>0.0)then
					if (xp<2*PI .and. xp>=0.0)then
						ff_out(m, n) = ff_out(m, n)
					end if 				
				endif					
			end do			
		end do
		do m = 1, 3*K_m+2
			xt = a_theta(m)
			do n = 1, n_phi
				xp = a_phi(n)
				if (xp>=2*PI)then
				ff_out(m, n) = ff_out(m, n-K_n)
				elseif (xp<0)then
				ff_out(m, n) = ff_out(m, n+K_n) !Eq.(7.10)
				end if
			end do
		end do	
	end subroutine	
	
	subroutine RT_function_SCS(K_m, struct, n, Ra, Int_fn_scs)
	! R_ab : distance between the two points
		implicit none		
		real(kind = dp), intent(in) :: Ra(3)		
		integer, intent(in) :: n, K_m  !sampling point on a sphere				
		integer :: i, j, pos_neg, ll, m_pairs
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roa(3), Ro(3), fm(3)
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		type(Integration_fn) :: Int_fn
		type(Integration_fn_scs), intent(out) :: Int_fn_scs
        
        !real(kind = dp), dimension(:, :), allocatable :: K_hat_p
        !real(kind = dp), dimension(:), allocatable :: TRF_p
        !real(kind = dp), dimension(:, :), allocatable :: transfer_theta, transfer_phi
                
        !call Khat_and_TRF_ExtendedK(K_m, K_hat_p, transfer_theta, transfer_phi, TRF_p)
        call Khat_and_TRF_ExtendedK(K_m)
        
		m_pairs = size(struct%neighbours)		
				
		call Quadrature_tri(ng, a, b, w)
		do ll = 1, 2
			call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
			pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
			dfm = pos_neg*Len_m 			
			allocate(Int_fn%pair_kk(ll)%Rfk_vec(N_sp+2, 2, 3))
			allocate(Int_fn_scs%Tri(ll)%Rfk_scs(N_sp+2, 2, 2))
			do i = 1, N_sp+2
				Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				do j = 1, ng
					Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
					fm = 0.5*(Ro-vm_t)*dfm
					Roa = Ro - Ra
					tmp = dot_product(k_hat_p(i, :), Roa(:))
					ff1 = exp(-im*k1*tmp)*w(j)
					ff2 = exp(-im*k2*tmp)*w(j)
					Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) + fm*ff1
					Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) + fm*ff2					
				end do
				Int_fn_scs%Tri(ll)%Rfk_scs(i, 1, 1) = dot_product(transfer_theta(i, 1:2), Int_fn%pair_kk(ll)%Rfk_vec(i, 1, 1:2)) !k1, fx1
				Int_fn_scs%Tri(ll)%Rfk_scs(i, 2, 1) = dot_product(transfer_theta(i, 1:2), Int_fn%pair_kk(ll)%Rfk_vec(i, 2, 1:2)) !k2, fx2
				Int_fn_scs%Tri(ll)%Rfk_scs(i, 1, 2) = dot_product(transfer_phi(i, 1:2), Int_fn%pair_kk(ll)%Rfk_vec(i, 1, 1:2))   !k1, fy1
				Int_fn_scs%Tri(ll)%Rfk_scs(i, 2, 2) = dot_product(transfer_phi(i, 1:2), Int_fn%pair_kk(ll)%Rfk_vec(i, 2, 1:2))   !k2, fy2
			end do !	
		end do
		!Rfk_scs(i, m, n)
		!i = N_sp+1 : theta = 0
		!i = N_sp+2 : theta = PI
		!m : medium
		!n : 1 = theta and 2 = phi
		!The last two Rfk_scs saves the fx, and fy at 0 and PI
		!Rfk_scs(, , 1) for fx, Rfk_scs(, , 2) for fy
	end subroutine RT_function_SCS
	
	subroutine Sorting_theta_phi_LenReduced(Component, p, K_m, ff_in, ff_out)
		implicit none
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_out
		complex(kind = dp), dimension(:, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m, p
		integer :: n_phi, m, n, K_n
		real(kind = 8), dimension(:), allocatable :: a_theta, a_phi
		real(kind = dp) :: xt, xp
		complex(kind = 8) :: fx1, fy1, fx2, fy2
		character(len = 1) :: Component
		!character(len = 50), public :: illumination
		intrinsic cos, sin
		
		fx1 = ff_in(K_m+1, 1)
		fy1 = ff_in(K_m+1, 2)
		fx2 = ff_in(K_m+2, 1)
		fy2 = ff_in(K_m+2, 2)
		
		K_n = 2*K_m
		n_phi = K_n+2*p+1
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		allocate(ff_out(K_m+2*p+2, K_n+2*p+1))
		ff_out=(0.0, 0.0)
		
		if (Component =='P') then		
			do m = 1, n_phi
			   n = m-p+K_n
				ff_out(p+1, m) = -sin(a_phi(n))*fx2 +cos(a_phi(n))*fy2  !theta = PI in the new 		
				ff_out(K_m+p+2, m) = -sin(a_phi(n))*fx1 +cos(a_phi(n))*fy1 !theta = 0 in the new array
			end do
		else if (Component =='T') then
			do m = 1, n_phi
			    n = m-p+K_n
				ff_out(p+1, m) = -cos(a_phi(n))*fx2 -sin(a_phi(n))*fy2  !theta = PI, here the problem with fx, and fy		
				ff_out(K_m+p+2, m) = cos(a_phi(n))*fx1+sin(a_phi(n))*fy1
			end do		
		else 
			print*, 'Wrong input component type'
		end if		
		ff_out(p+2:K_m+p+1, p+1:K_n+p) = ff_in(1:K_m, 1:K_n)
		
		do m = 1, K_m+2*p+2!3*K_m+2
			print*, 'm=', m
			xt = a_theta(m-p + K_m+1)
			do n = 1, K_n+2*p+1!n_phi
				xp = a_phi(n-p+K_n)
				if (xt<0.0)then
					if (xp<PI .and. xp>=0.0) then
						ff_out(m, n) = -ff_out(2*(K_m+p+2)-m, n+K_n/2)
					elseif (xp<2*PI .and. xp>=PI) then
						ff_out(m, n) = -ff_out(2*(K_m+p+2)-m, n-K_n/2)
					endif
				elseif (xt>=PI) then
					if (xp<PI .and. xp>=0.0)then					 
						ff_out(m, n) = -ff_out(2*(p+1)-m, n+K_n/2)
					elseif (xp<2*PI .and. xp>=PI)then
						ff_out(m, n) = -ff_out(2*(p+1)-m, n-K_n/2)
					end if
				elseif (xt<PI .and. xt>0.0)then
					if (xp<2*PI .and. xp>=0.0)then
						ff_out(m, n) = ff_out(m, n)
					end if
				endif
			end do
		end do
		do m = 1, k_m+2*p+2
			xt = a_theta(m-p+k_m+1)
			do n = 1, k_n+2*p+1
				xp = a_phi(n-p+k_n)
				if (xp>=2*pi)then
				ff_out(m, n) = ff_out(m, n-k_n)
				elseif (xp<0)then
				ff_out(m, n) = ff_out(m, n+k_n) !eq.(7.10)
				end if
			end do
		end do
	end subroutine	Sorting_theta_phi_LenReduced
	
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
	
	subroutine Sorting_theta_phi_LenReduced_Arr(p, K_m, ff_in, ff_out)
		implicit none
		complex(kind = dp), dimension(:, :, :), allocatable, intent(out) :: ff_out
		complex(kind = dp), dimension(:, :, :), allocatable, intent(in) :: ff_in
		integer, intent(in) :: K_m, p
		integer :: n_phi, m, n, K_n
		real(kind = 8), dimension(:), allocatable :: a_theta, a_phi
		real(kind = dp) :: xt, xp
		complex(kind = 8) :: fx1, fy1, fx2, fy2
		intrinsic cos, sin
		
		K_n = 2*K_m
		n_phi = K_n+2*p+1
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		allocate(ff_out(K_m+2*p+2, K_n+2*p+1, 4))
		ff_out=(0.0, 0.0)
		
		fx1 = ff_in(K_m+1, 1, 3)
		fy1 = ff_in(K_m+1, 2, 3)
		fx2 = ff_in(K_m+2, 1, 3)
		fy2 = ff_in(K_m+2, 2, 3)
		!if (Component =='P') then		
		do m = 1, n_phi
		   n = m-p+K_n
			ff_out(p+1, m, 3) = -sin(a_phi(n))*fx2 +cos(a_phi(n))*fy2  !theta = PI in the new 		
			ff_out(K_m+p+2, m, 3) = -sin(a_phi(n))*fx1 +cos(a_phi(n))*fy1 !theta = 0 in the new array		
		end do
		fx1 = ff_in(K_m+1, 1, 4)
		fy1 = ff_in(K_m+1, 2, 4)
		fx2 = ff_in(K_m+2, 1, 4)
		fy2 = ff_in(K_m+2, 2, 4)
		do m = 1, n_phi
		   n = m-p+K_n
			ff_out(p+1, m, 4) = -sin(a_phi(n))*fx2 +cos(a_phi(n))*fy2  !theta = PI in the new 		
			ff_out(K_m+p+2, m, 4) = -sin(a_phi(n))*fx1 +cos(a_phi(n))*fy1 !theta = 0 in the new array
		end do
		
		fx1 = ff_in(K_m+1, 1, 1)
		fy1 = ff_in(K_m+1, 2, 1)
		fx2 = ff_in(K_m+2, 1, 1)
		fy2 = ff_in(K_m+2, 2, 1)
		do m = 1, n_phi
			n = m-p+K_n
			ff_out(p+1, m, 1) = -cos(a_phi(n))*fx2 -sin(a_phi(n))*fy2  !theta = PI, here the problem with fx, and fy		
			ff_out(K_m+p+2, m, 1) = cos(a_phi(n))*fx1+sin(a_phi(n))*fy1
		end do
		fx1 = ff_in(K_m+1, 1, 2)
		fy1 = ff_in(K_m+1, 2, 2)
		fx2 = ff_in(K_m+2, 1, 2)
		fy2 = ff_in(K_m+2, 2, 2)
		do m= 1, n_phi
			n = m-p+K_n
			ff_out(p+1, m, 2) = -cos(a_phi(n))*fx2 -sin(a_phi(n))*fy2  !theta = PI, here the problem with fx, and fy		
			ff_out(K_m+p+2, m, 2) = cos(a_phi(n))*fx1+sin(a_phi(n))*fy1
		end do
		
		ff_out(p+2:K_m+p+1, p+1:K_n+p, 1:4) = ff_in(1:K_m, 1:K_n, 1:4)
		
		do m = 1, K_m+2*p+2!3*K_m+2
			print*, 'm=', m
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
				ff_out(m, n, 1:4) = ff_out(m, n+k_n, 1:4) !eq.(7.10)
				end if
			end do
		end do
	end subroutine	Sorting_theta_phi_LenReduced_Arr
	
	subroutine Test_Lagrange_Interpolated_LenReduced_Arr(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n, ll, tt
		integer :: it, jp, m_phi, n_phi, t, s, pp, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn, Int_fn_inter
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:, :, :), allocatable :: ff_in, ff_out		
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		
		nd = 11
		ll = 1 !index of the triangle in the pair
		pp = 5 !points for interpolation !p should be smaller than Km
		K_m = 12
		K_n = K_m + 2
		n_phi = 2*K_n
		
		m_phi = 2*K_m
		dphi = 2*PI/m_phi		
		allocate(ff_in(K_m+2, 2*K_m, 4), ff_out(K_m+2, 2*K_m, 4), ff_inter(n_phi*K_n, 4)) 
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(K_m, struct, nd, Ra, Int_fn)	
		
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_in(i, j, 1)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 1) !theta in k1
				ff_in(i, j, 2)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 1) !theta in k2
				ff_in(i, j, 3)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 2) !phi in k1
				ff_in(i, j, 4)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 2)	!phi in k2			
			end do
		end do
		
		ff_in(K_m+1:K_m+2, 1:2, 1)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) !fx and fy, at theta = zero at theta = pi, medium 1
		ff_in(K_m+1:K_m+2, 1:2, 2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) !fx and fy, at theta = zero at theta = pi, medium 2
		ff_in(K_m+1:K_m+2, 1:2, 3)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) 
		ff_in(K_m+1:K_m+2, 1:2, 4)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) 
		a = -1.0 
		b = 1.0		
		call Sorting_theta_phi_LenReduced_Arr(pp, K_m, ff_in, ff_out)		
		
		call theta_phi_extended(K_m, a_phi, a_theta)
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
		allocate(Int_fn_inter%Tri(ll)%Rfk_scs(n_phi*K_n, 2, 2)) 		
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
		Int_fn_inter%Tri(ll)%Rfk_scs(:, 1, 1) = ff_inter(:, 1)
		Int_fn_inter%Tri(ll)%Rfk_scs(:, 2, 1) = ff_inter(:, 2)
		Int_fn_inter%Tri(ll)%Rfk_scs(:, 1, 2) = ff_inter(:, 3)
		Int_fn_inter%Tri(ll)%Rfk_scs(:, 2, 2) = ff_inter(:, 4)
		
		open (unit=206, file = 'ff_inter_theta1_real_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_theta2_real_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_phi1_real_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_phi2_real_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_theta1_imag_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_theta2_imag_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_phi1_imag_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_phi2_imag_Km8_Arr.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		end do
		close(206)
		!!
		!open (unit=206, file = 'ff_in_phi1_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi2_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta1_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta2_real_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi1_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi2_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta1_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta2_imag_Km13_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		!end do
		!close(206)	
		!
	end subroutine Test_Lagrange_Interpolated_LenReduced_Arr
	
	subroutine Test_Lagrange_Interpolated_LenReduced(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n, ll, tt
		integer :: it, jp, m_phi, n_phi, t, s, pp, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn, Int_fn_inter
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:, :), allocatable :: ff_in, ff_out		
		complex(kind = dp), dimension(:), allocatable :: ff_in_theta, ff_in_phi
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		complex(kind = dp), dimension(:, :), allocatable :: ff_in_theta_I, ff_in_phi_I, ff_out_theta_I, ff_out_phi_I 
		complex(kind = dp), dimension(:, :), allocatable :: ff_in_theta_II, ff_in_phi_II, ff_out_theta_II, ff_out_phi_II
		
		nd = 11
		ll = 1 !index of the triangle in the pair
		pp = 5 !points for interpolation !p should be smaller than Km
		K_m = 11
		K_n = K_m + 2
		n_phi = 2*K_n
		
		m_phi = 2*K_m
		dphi = 2*PI/m_phi
		allocate(ff_in_theta_I(K_m+2, 2*K_m), ff_in_theta_II(K_m+2, 2*K_m))
		allocate(ff_in_phi_I(K_m+2, 2*K_m), ff_in_phi_II(K_m+2, 2*K_m))	
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*pp), wk_phi(2*pp))
		
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(K_m, struct, nd, Ra, Int_fn)	
		
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_in_theta_I(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 1) !theta in k1
				ff_in_theta_II(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 1) !theta in k2
				ff_in_phi_I(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 2)
				ff_in_phi_II(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 2)				
			end do
		end do
		
		ff_in_theta_I(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) !fx, at theta = zero at theta = pi 
		ff_in_theta_II(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) !fx, at theta = zero 
		ff_in_phi_I(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) !fx, at theta = zero 
		ff_in_phi_II(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) !fx, at theta = zero 
				
		a = -1.0 
		b = 1.0		
		call Sorting_theta_phi_LenReduced('T', pp, K_m, ff_in_theta_I, ff_out_theta_I)
		call Sorting_theta_phi_LenReduced('T', pp, K_m, ff_in_theta_II, ff_out_theta_II)
		call Sorting_theta_phi_LenReduced('P', pp, K_m, ff_in_phi_I, ff_out_phi_I)
		call Sorting_theta_phi_LenReduced('P', pp, K_m, ff_in_phi_II, ff_out_phi_II)!
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		call legendre_handle(K_n, a, b, xt, x_tmp)	
		call legendre_handle(K_m, a, b, xt_old, x_tmp_old)	
		do i = 1, n_phi
			xp(i) = (i-1)*6*PI/n_phi
		end do		
		do i = 1, m_phi
			xp_old(i) = (i-1)*2*PI/m_phi		
		end do		
		!
		!!it: index for the interpolated theta
		allocate(Int_fn_inter%Tri(ll)%Rfk_scs(n_phi*K_n, 2, 2)) 
		do it = 1, K_n
			theta = acos(xt(it))			
			t = floor((it-1.0)/K_n*K_m) + 1
			kt = pp + t
			call Lagrange_Interpolation_Weights(theta, t, pp, xt_old, a_theta, wk_theta)
			do jp = 1, n_phi			   
				tt = (it-1)*n_phi + jp
				Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1:2, 1:2) = (0.0, 0.0)
				phi = (jp-1)*2*PI/n_phi
				s = floor((jp-1.0)/K_n*K_m)+1	
				kp = pp + s
				call Lagrange_Interpolation_Weights(phi, s, pp, xp_old, a_phi, wk_phi)
				do i = 1, 2*pp
					do j = 1, 2*pp
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 1) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 1)+wk_phi(j)*wk_theta(i)*ff_out_theta_I(kt+i-pp, kp+j-pp) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 1) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 1)+wk_phi(j)*wk_theta(i)*ff_out_theta_II(kt+i-pp, kp+j-pp) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 2) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 2)+ wk_phi(j)*wk_theta(i)*ff_out_phi_I(kt+i-pp, kp+j-pp) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 2) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 2)+ wk_phi(j)*wk_theta(i)*ff_out_phi_II(kt+i-pp, kp+j-pp) !					
					end do
				end do				
			end do
		end do
		
		!open (unit=206, file = 'ff_inter_theta1_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta2_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi1_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi2_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta1_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_theta2_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi1_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 2))
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_inter_phi2_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		!do i = 1, K_n*n_phi
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 2))
		!end do
		!close(206)
		!####################################
		open (unit=206, file = 'ff_in_phi1_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_phi2_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_theta1_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_theta2_real_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  real(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_phi1_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_phi2_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_theta1_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		end do
		close(206)
		
		open (unit=206, file = 'ff_in_theta2_imag_Km8_LenReduced.txt', action="write",status = 'replace')		
		do i = 1, K_m*2*K_m
			write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		end do
		close(206)	
		!
	end subroutine Test_Lagrange_Interpolated_LenReduced
	
	subroutine Test_AKK_Function(struct, struct_cube)
		implicit none
		integer :: nd, n_cube
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn
		real(kind = dp) :: Ra(3)
		complex(kind = dp), dimension(:), allocatable :: ff_in_theta, ff_in_phi
		complex(kind = dp), dimension(:), allocatable:: ff_kp 
		real(kind = 8) :: x, phi
		complex(kind = 8) :: suma, ff_kq
		real(kind = 8), dimension(:), allocatable :: y
		real(kind = 8), dimension(:, :), allocatable :: xw_p, xw_q	
		integer :: m, n, Kp, Kq, Nq, Np, i, j, fnu
		type(list_list_cmplx) :: Ylm_q, Ylm_p			
		intrinsic :: conjg		
		Kq = 11
		Kp = int(Kq*PI)  !new sampling rate		
		Nq = Kq*2*Kq
		Np = Kp*2*Kp
		xw_p = theta_phi_w(Kp) !new 
		xw_q = theta_phi_w(Kq) !original  
		suma= (0.0, 0.0)
		!allocate(Akk(1:Nq, 1:Np))
		allocate(y(0:Kq))
		allocate(ff_kp(Np))
		nd=11	
		y = C_lk(Kp, Kq)
		allocate(ff_in_theta(N_sp+2), ff_in_phi(N_sp+2))
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(Kq, struct, nd, Ra, Int_fn)	
		ff_in_theta = Int_fn%Tri(1)%Rfk_scs(:, 1, 1)
		ff_in_phi = Int_fn%Tri(1)%Rfk_scs(:, 1, 2)
		
		do m = 1, Np
			ff_kp(m) = (0.0, 0.0)
			if (modulo(m,50).eq.0) then
				print *, 'm=', m, 'of total', Np
			end if	
			call spherical_harmonics_Darve(Kq+1 , Ylm_p, xw_p(m, 1), xw_p(m, 2))	
			do n = 1, Nq				
				call spherical_harmonics_Darve(Kq+1, Ylm_q, xw_q(n, 1), xw_q(n, 2))				
				do j = 0, Kq
					do i = -j, j
						ff_kp(m) = ff_kp(m) + conjg(Ylm_p%item(j)%item(i))*Ylm_q%item(j)%item(i)*xw_q(n, 3)*ff_in_theta(n)*y(j)
					end do
				end do
			end do		
		end do
		do i = 1, Np
			suma= suma+ff_kp(i)*xw_p(i, 3)
		end do
		print*, 'sum fp=', suma
		suma = (0.0, 0.0)
		do i = 1, Nq
			suma= suma+ff_in_theta(i)*xw_q(i, 3)
		end do
		print*, 'sum fq=', suma
		
		!print*, 'ff_kp =', ff_kp
		!print*, 'ff_kq =', ff_kq
	end subroutine Test_Akk_Function
	
	subroutine Akk_interpolation()
		implicit none
		real(kind = 8) :: x, phi
		complex(kind = 8) :: tmp
		real(kind = 8), dimension(:), allocatable :: y
		real(kind = 8), dimension(:, :), allocatable :: xw_p, xw_q	
		integer :: m, n, Kp, Kq, Nq, Np, i, j, fnu
		type(list_list_cmplx) :: Ylm_q, Ylm_p			
		intrinsic :: conjg		
		Kp = 11
		Kq = Kp!int(Kp*PI)  !new sampling rate		
		Np = Kp*2*Kp
		Nq = Kq*2*Kq
		xw_p = theta_phi_w(Kp) !original 
		xw_q = theta_phi_w(Kq) ! new 
		allocate(Akk(1:Nq, 1:Np))
		allocate(y(0:Kq))
		y = C_lk(Kp, Kq)
		Akk(1:Nq, 1:Np) = (0.0, 0.0)
		do m = 11, 11!Nq		   			
			call spherical_harmonics_Darve(Kq+1 , Ylm_q, xw_q(m, 1), xw_q(m, 2))	
			do n = 1, Np	
			   Akk(m, n) = (0.0, 0.0)
				call spherical_harmonics_Darve(Kq+1, Ylm_p, xw_p(n, 1), xw_p(n, 2))
				do j = 0, Kp
					do i = -j, j
						Akk(m, n) = Akk(m, n) + Ylm_q%item(j)%item(i)*conjg(Ylm_p%item(j)%item(i))*xw_p(n, 3)*y(j)
					end do
				end do
			end do
		end do		
	end subroutine
	
	function theta_phi_w(K_m)
		real(kind = 8), dimension(:, :), allocatable :: theta_phi_w
		integer, intent(in) :: K_m
		real(kind = 8) :: phi, d_phi, a, b
		integer :: Np, m_ph
		real(kind = 8), dimension(:), allocatable :: xx, ww
		integer :: j, k, i
		intrinsic :: sqrt
		
		m_ph = 2*K_m
		Np = K_m*m_ph		
		allocate(theta_phi_w(Np, 3))		
		allocate(xx(K_m))
		allocate(ww(K_m))
		a = -1.0
		b = 1.0
		call legendre_handle(K_m, a, b, xx, ww)		
		d_phi = 2*PI/m_ph
		do j = 1, K_m
         do k = 1, m_ph
            phi = (k-1)*d_phi 
            i = k + (j-1)*m_ph 
				theta_phi_w(i, 1) = xx(j)
				theta_phi_w(i, 2) = phi
				theta_phi_w(i, 3) = ww(j)*d_phi
         end do
		end do
	end function
	
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
	
	subroutine Test_Lagrange_Interpolated_Function_2Step(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n
		integer :: it, jp, m_phi, n_phi, t, s, p, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:), allocatable :: ff_in_theta, ff_in_phi
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		complex(kind = dp), dimension(:, :), allocatable :: ff_tmp_in, ff_out_theta, ff_out_phi, ff_inter_theta, ff_inter_phi, ff_tmp_t, ff_tmp_p
		
		nd = 11
		p = 3
		K_m = 11
		K_n = K_m + 3
		n_phi = 2*K_n
		m_phi = 2*K_m
		dphi = 2*PI/m_phi
		allocate(ff_inter_theta(K_n, n_phi), ff_inter_phi(K_n, n_phi))
		allocate(ff_in_theta(N_sp+2), ff_in_phi(N_sp+2))
		allocate(ff_tmp_in(K_m+2, 2*K_m))
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*p), wk_phi(2*p))
		allocate(ff_tmp_t(K_n, 3*m_phi), ff_tmp_p(K_n, 3*m_phi))
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(K_m, struct, nd, Ra, Int_fn)	
		ff_in_theta = Int_fn%Tri(1)%Rfk_scs(:, 1, 1)
		ff_in_phi = Int_fn%Tri(1)%Rfk_scs(:, 1, 2)
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_tmp_in(i, j)= ff_in_theta(k)
			end do
		end do		
		ff_tmp_in(K_m+1, 1) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1, 1, 1) !fx, at theta = zero 
		ff_tmp_in(K_m+1, 2) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1, 1, 2) !fy
		ff_tmp_in(K_m+2, 1) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+2, 1, 1) !fx, at theta = pi 
		ff_tmp_in(K_m+2, 2) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+2, 1, 2) !fy, phi
		
		a = -1.0 
		b = 1.0		
		call Sorting_theta_phi('T', K_m, ff_tmp_in, ff_out_theta)
		
		do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_tmp_in(i, j)= ff_in_phi(k)
			end do
		end do		
		
		call Sorting_theta_phi('P', K_m, ff_tmp_in, ff_out_phi)		
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		call legendre_handle(K_n, a, b, xt, x_tmp)	
		call legendre_handle(K_m, a, b, xt_old, x_tmp_old)	
		do i = 1, n_phi
			xp(i) = (i-1)*6*PI/n_phi
		end do
		
		do i = 1, m_phi
			xp_old(i) = (i-1)*2*PI/m_phi		
		end do		
		
		!it: index for the interpolated theta
		do it = 1, K_n
			theta = acos(xt(it))
			t = floor((it-1.0)/K_n*K_m)+1 
			!print*, 't=', t
			kt = K_m + t
			call Lagrange_Interpolation_Weights(theta, t, p, xt_old, a_theta, wk_theta)			
			ff_tmp_t(it, :) = (0.0, 0.0)
			ff_tmp_p(it, :) = (0.0, 0.0)
			do i = 1, 2*p
				ff_tmp_t(it, :) = ff_tmp_t(it, :) + wk_theta(i)*ff_out_theta(kt+i-p, :) !
				ff_tmp_p(it, :) = ff_tmp_p(it, :) + wk_theta(i)*ff_out_phi(kt+i-p, :) !
			end do				
		end do	
		
		do jp = 1, n_phi
			ff_inter_theta(:, jp) = (0.0, 0.0)
			ff_inter_phi(:, jp) = (0.0, 0.0)
			phi = (jp-1)*2*PI/n_phi
			s = floor((jp-1.0)/K_n*K_m)+1	
			kp = m_phi + s
			call Lagrange_Interpolation_Weights(phi, s, p, xp_old, a_phi, wk_phi)
			do j = 1, 2*p
				ff_inter_theta(:, jp) = ff_inter_theta(:, jp) + wk_phi(j)*ff_tmp_t(:, kp+j-p) !
				ff_inter_phi(:, jp) = ff_inter_phi(:, jp) + wk_phi(j)*ff_tmp_p(:, kp+j-p)!!
			end do				
		end do
	end subroutine Test_Lagrange_Interpolated_Function_2Step	
	
	subroutine Test_Lagrange_Interpolated_New(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n, ll, tt
		integer :: it, jp, m_phi, n_phi, t, s, p, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn, Int_fn_inter
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:), allocatable :: ff_in_theta, ff_in_phi
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		complex(kind = dp), dimension(:, :), allocatable :: ff_in_theta_I, ff_in_phi_I, ff_out_theta_I, ff_out_phi_I 
		complex(kind = dp), dimension(:, :), allocatable :: ff_in_theta_II, ff_in_phi_II, ff_out_theta_II, ff_out_phi_II
		
		nd = 11
		ll = 1 !index of the triangle in the pair
		p = 5
		K_m = 11
		K_n = K_m + 2
		n_phi = 2*K_n
		
		m_phi = 2*K_m
		dphi = 2*PI/m_phi
		allocate(ff_in_theta_I(K_m+2, 2*K_m), ff_in_theta_II(K_m+2, 2*K_m))
		allocate(ff_in_phi_I(K_m+2, 2*K_m), ff_in_phi_II(K_m+2, 2*K_m))
		
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*p), wk_phi(2*p))
		
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(K_m, struct, nd, Ra, Int_fn)	
		
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_in_theta_I(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 1)
				ff_in_theta_II(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 1)				
				ff_in_phi_I(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 1, 2)
				ff_in_phi_II(i, j)= Int_fn%Tri(ll)%Rfk_scs(k, 2, 2)				
			end do
		end do
		
		ff_in_theta_I(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) !fx, at theta = zero at theta = pi 
		ff_in_theta_II(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) !fx, at theta = zero 
		ff_in_phi_I(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 1, 1:2) !fx, at theta = zero 
		ff_in_phi_II(K_m+1:K_m+2, 1:2)= Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1:K_m*2*K_m+2, 2, 1:2) !fx, at theta = zero 
		
		a = -1.0 
		b = 1.0		
		call Sorting_theta_phi('T', K_m, ff_in_theta_I, ff_out_theta_I)
		call Sorting_theta_phi('T', K_m, ff_in_theta_II, ff_out_theta_II)
		call Sorting_theta_phi('P', K_m, ff_in_phi_I, ff_out_phi_I)
		call Sorting_theta_phi('P', K_m, ff_in_phi_II, ff_out_phi_II)!
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		call legendre_handle(K_n, a, b, xt, x_tmp)	
		call legendre_handle(K_m, a, b, xt_old, x_tmp_old)	
		do i = 1, n_phi
			xp(i) = (i-1)*6*PI/n_phi
		end do		
		do i = 1, m_phi
			xp_old(i) = (i-1)*2*PI/m_phi		
		end do		
		!
		!!it: index for the interpolated theta
		allocate(Int_fn_inter%Tri(ll)%Rfk_scs(n_phi*K_n, 2, 2)) 
		do it = 1, K_n
			theta = acos(xt(it))			
			t = floor((it-1.0)/K_n*K_m)+1 		
			kt = K_m + t
			call Lagrange_Interpolation_Weights(theta, t, p, xt_old, a_theta, wk_theta)
			do jp = 1, n_phi			   
				tt = (it-1)*n_phi + jp
				Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1:2, 1:2) = (0.0, 0.0)
				phi = (jp-1)*2*PI/n_phi
				s = floor((jp-1.0)/K_n*K_m)+1	
				kp = m_phi + s
				call Lagrange_Interpolation_Weights(phi, s, p, xp_old, a_phi, wk_phi)
				do i = 1, 2*p
					do j = 1, 2*p
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 1) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 1)+wk_phi(j)*wk_theta(i)*ff_out_theta_I(kt+i-p, kp+j-p) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 1) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 1)+wk_phi(j)*wk_theta(i)*ff_out_theta_II(kt+i-p, kp+j-p) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 2) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 1, 2)+ wk_phi(j)*wk_theta(i)*ff_out_phi_I(kt+i-p, kp+j-p) !
						Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 2) = Int_fn_inter%Tri(ll)%Rfk_scs(tt, 2, 2)+ wk_phi(j)*wk_theta(i)*ff_out_phi_II(kt+i-p, kp+j-p) !					
					end do
				end do				
			end do
		end do
		
		open (unit=206, file = 'ff_inter_theta1_real_Km8_new.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 1, 1))
		end do
		close(206)
		
		open (unit=206, file = 'ff_inter_theta2_real_Km8_new.txt', action="write",status = 'replace')		
		do i = 1, K_n*n_phi
			write (206, '(200(es19.12, tr5))')  real(Int_fn_inter%Tri(ll)%Rfk_scs(i, 2, 1))
		end do
		close(206)
		
		!open (unit=206, file = 'ff_in_phi1_imag_Km8.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_phi2_imag_Km8.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 2)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta1_imag_Km8.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 1, 1)) ! 
		!end do
		!close(206)
		!
		!open (unit=206, file = 'ff_in_theta2_imag_Km8.txt', action="write",status = 'replace')		
		!do i = 1, K_m*2*K_m
		!	write (206, '(200(es19.12, tr5))')  imag(Int_fn%Tri(ll)%Rfk_scs(i, 2, 1)) ! 
		!end do
		!close(206)		
		
	end subroutine Test_Lagrange_Interpolated_New
	
	subroutine Test_Lagrange_Interpolated_Function(struct, struct_cube)
		implicit none
		integer :: nd, n_cube, i, j, k, K_m, K_n
		integer :: it, jp, m_phi, n_phi, t, s, p, kp, kt
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		type(Integration_fn_scs) :: Int_fn
		real(kind = dp) :: Ra(3), dphi, a, b, theta, phi
		complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:), allocatable :: ff_in_theta, ff_in_phi
		real(kind = dp), dimension(:), allocatable :: x_tmp, xp, xp_old, a_phi, a_theta, wk_phi, xt, wk_theta, xt_old, x_tmp_old
		complex(kind = dp), dimension(:, :), allocatable :: ff_tmp_in, ff_out_theta, ff_out_phi, ff_inter_theta, ff_inter_phi
		
		nd = 11
		p = 3
		K_m = 11
		K_n = K_m - 3
		n_phi = 2*K_n
		m_phi = 2*K_m
		dphi = 2*PI/m_phi
		allocate(ff_inter_theta(K_n, n_phi), ff_inter_phi(K_n, n_phi))
		allocate(ff_in_theta(N_sp+2), ff_in_phi(N_sp+2))
		allocate(ff_tmp_in(K_m+2, 2*K_m))
		allocate(xp(n_phi), xp_old(m_phi))
		allocate(a_theta(3*K_m + 2))
		allocate(a_phi(3*m_phi))
		allocate(wk_theta(2*p), wk_phi(2*p))
		
		n_cube = struct_cube%pair_cubes(nd)%cube_index!
		Ra = struct_cube%cubes(n_cube)%cube_position
		call RT_function_SCS(K_m, struct, nd, Ra, Int_fn)	
		ff_in_theta = Int_fn%Tri(1)%Rfk_scs(:, 2, 1)
		ff_in_phi = Int_fn%Tri(1)%Rfk_scs(:, 2, 2)
	   do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_tmp_in(i, j)= ff_in_theta(k)
			end do
		end do		
		ff_tmp_in(K_m+1, 1) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1, 2, 1) !fx, at theta = zero 
		ff_tmp_in(K_m+1, 2) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+1, 2, 2) !fy
		ff_tmp_in(K_m+2, 1) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+2, 2, 1) !fx, at theta = pi 
		ff_tmp_in(K_m+2, 2) = Int_fn%Tri(1)%Rfk_scs(K_m*2*K_m+2, 2, 2) !fy, phi
		
		a = -1.0 
		b = 1.0		
		call Sorting_theta_phi('T', K_m, ff_tmp_in, ff_out_theta)
		do i = 1, K_m
			do j = 1, 2*K_m
				k = (i-1)*2*K_m + j
				ff_tmp_in(i, j)= ff_in_phi(k)
			end do
		end do		
		!
		call Sorting_theta_phi('P', K_m, ff_tmp_in, ff_out_phi)		
		
		call theta_phi_extended(K_m, a_phi, a_theta)
		call legendre_handle(K_n, a, b, xt, x_tmp)	
		call legendre_handle(K_m, a, b, xt_old, x_tmp_old)	
		do i = 1, n_phi
			xp(i) = (i-1)*6*PI/n_phi
		end do
		
		do i = 1, m_phi
			xp_old(i) = (i-1)*2*PI/m_phi		
		end do		
		
		!it: index for the interpolated theta
		do it = 1, K_n
			theta = acos(xt(it))
			t = floor((it-1.0)/K_n*K_m)+1 
			!print*, 't=', t
			kt = K_m + t
			call Lagrange_Interpolation_Weights(theta, t, p, xt_old, a_theta, wk_theta)
			do jp = 1, n_phi
				ff_inter_theta(it, jp) = (0.0, 0.0)
				ff_inter_phi(it, jp) = (0.0, 0.0)
				phi = (jp-1)*2*PI/n_phi
				s = floor((jp-1.0)/K_n*K_m)+1	
				kp = m_phi + s
				call Lagrange_Interpolation_Weights(phi, s, p, xp_old, a_phi, wk_phi)
				do i = 1, 2*p
					do j = 1, 2*p
						ff_inter_theta(it, jp) = ff_inter_theta(it, jp) + wk_phi(j)*wk_theta(i)*ff_out_theta(kt+i-p, kp+j-p) !
						ff_inter_phi(it, jp) = ff_inter_phi(it, jp) + wk_phi(j)*wk_theta(i)*ff_out_phi(kt+i-p, kp+j-p) !
					end do
				end do				
			end do
		end do
		
		open (unit=206, file = 'ff_inter_phi_real_Km8.txt', action="write",status = 'replace')		
		do i = 1, K_n
			write (206, '(200(es19.12, tr5))')  (real(ff_inter_phi(i, j)), j= 1,n_phi) ! 
		end do
		close(206)
		
		!open (unit=206, file = 'ff_in_phi_real_Km8.txt', action="write",status = 'replace')		
		!do i = 1, K_m
		!	write (206, '(200(es19.12, tr5))')  (real(ff_tmp_in(i, j)), j= 1,m_phi) ! 
		!end do
		!close(206)
	end subroutine 	
	
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
	
	!subroutine Sorting_theta_phi(Component, K_m, ff_in, ff_out)
	!	implicit none
	!	complex(kind = dp), dimension(:, :), allocatable, intent(out) :: ff_out
	!	complex(kind = dp), dimension(:, :), allocatable, intent(in) :: ff_in
	!	integer, intent(in) :: K_m
	!	integer :: n_phi, m, n, K_n
	!	real(kind = 8), dimension(:), allocatable :: a_theta, a_phi
	!	real(kind = dp) :: xt, xp
	!	complex(kind = 8) :: fx1, fy1, fx2, fy2
	!	character(len = 1) :: Component
	!	!character(len = 50), public :: illumination
	!	intrinsic cos, sin
	!	
	!	fx1 = ff_in(K_m+1, 1)
	!	fy1 = ff_in(K_m+1, 2)
	!	fx2 = ff_in(K_m+2, 1)
	!	fy2 = ff_in(K_m+2, 2)
	!	
	!	K_n = 2*K_m
	!	n_phi = 3*K_n
	!	
	!	call theta_phi_extended(K_m, a_phi, a_theta)
	!	allocate(ff_out(3*K_m+2, 3*K_n))
	!	ff_out=(0.0, 0.0)
	!	
	!	if (Component =='P') then		
	!		do m = 1, n_phi
	!			ff_out(K_m+1, m) = -sin(a_phi(m))*fx2 +cos(a_phi(m))*fy2  !theta = PI in the new 		
	!			ff_out(2*K_m+2, m) = -sin(a_phi(m))*fx1 +cos(a_phi(m))*fy1 !theta = 0 in the new array
	!		end do
	!	else if (Component =='T') then
	!		do m = 1, n_phi
	!			ff_out(K_m+1, m) = -cos(a_phi(m))*fx2 -sin(a_phi(m))*fy2  !theta = PI, here the problem with fx, and fy		
	!			ff_out(2*K_m+2, m) = cos(a_phi(m))*fx1+sin(a_phi(m))*fy1
	!		end do		
	!	else 
	!		print*, 'Wrong input component'
	!	end if		
	!	ff_out(K_m+2:2*K_m+1, K_n+1:2*K_n) = ff_in(1:K_m, 1:2*K_m)
	!	
	!	do m = 1, 3*K_m+2
	!		xt = a_theta(m)
	!		do n = 1, n_phi
	!			xp = a_phi(n)
	!			if (xt<0.0)then
	!				if (xp<PI .and. xp>=0.0) then
	!					ff_out(m, n) = -ff_out(4*K_m+4-m, n+K_n/2)					
	!				elseif (xp<2*PI .and. xp>=PI) then
	!					ff_out(m, n) = -ff_out(4*K_m+4-m, n-K_n/2)
	!				endif
	!			elseif (xt>PI) then
	!				if (xp<PI .and. xp>=0.0)then
	!					ff_out(m, n) = -ff_out(2*(K_m+1)-m, n+K_n/2)
	!				elseif (xp<2*PI .and. xp>=PI)then
	!					ff_out(m, n) = -ff_out(2*(K_m+1)-m, n-K_n/2)
	!				end if 
	!			elseif (xt<PI .and. xt>0.0)then
	!				if (xp<2*PI .and. xp>=0.0)then
	!					ff_out(m, n) = ff_out(m, n)
	!				end if 				
	!			endif					
	!		end do			
	!	end do
	!	do m = 1, 3*K_m+2
	!		xt = a_theta(m)
	!		do n = 1, n_phi
	!			xp = a_phi(n)
	!			if (xp>=2*PI)then
	!			ff_out(m, n) = ff_out(m, n-K_n)
	!			elseif (xp<0)then
	!			ff_out(m, n) = ff_out(m, n+K_n) !Eq.(7.10)
	!			end if
	!		end do
	!	end do
	!	!
	!	!open (unit=206, file = 'f_in.txt', action="write",status = 'replace')		
	!	!do m = 1, K_m
	!	!	write (206, '(200(es19.12, tr5))')  (real(ff_in(m, n)), n= 1, 2*K_m)
	!	!end do
	!	!close(206)
	!	!
	!	!open (unit=206, file = 'a_theta.txt', action="write",status = 'replace')		
	!	!do m = 1, 3*K_m+2
	!	!	write (206, '(200(es19.12, tr5))')  a_theta(m)			
	!	!end do
	!	!close(206)
	!	!
	!	!open (unit=206, file = 'a_phi.txt', action="write",status = 'replace')		
	!	!do m = 1, 3*K_n
	!	!	write (206, '(200(es19.12, tr5))')  a_phi(m)			
	!	!end do
	!	!close(206)
	!end subroutine	
	

	!With theta = PI and 0 at N_sp+1 and N_sp+2 for interpolation
	!subroutine FMM_edge_Tri_RT_ExtendedK(m, n, struct, struct_cube, d_c, suma, sumb, sumc) !	
	!!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
	!	implicit none		
	!	integer, intent(in) :: m, n
	!	type(Structure), intent(in) :: struct
	!	type(Structure_cube), intent(in) :: struct_cube
	!	complex(kind = dp), dimension(N_sp+2) :: sum_a, sum_c, sum_b!
	!	complex(kind = dp), dimension(N_sp+2, 2) :: Rf, Tf
	!	complex(kind = dp), dimension(N_sp+2) :: TL_1, TL_2
	!	complex(kind = dp), dimension(N_sp+2, 3) :: ss1, ss2, R1, R2, T1, T2
	!	complex(kind = dp), dimension(N_sp+2, 2, 3) :: Rf_vec, Tf_vec
	!	complex(kind = dp), intent(out) :: suma, sumb, sumc
	!	complex(kind = dp), dimension(N_sp+2) :: f1, g1, f2, g2
	!	real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
	!	real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
	!	type(Integration_fn) :: Int_fn, Int_fm
	!	integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
	!	intrinsic :: sqrt	
	!	
	!	m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
	!	n_cube = struct_cube%pair_cubes(n)%cube_index
	!	Ra = struct_cube%cubes(m_cube)%cube_position
	!	Rb = struct_cube%cubes(n_cube)%cube_position	
	!	Rab = Ra - Rb
	!	
	!	Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
	!	TL_1 = TL_km_arr_ExtendedK(Rab, Lm, k1)
	!	Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
	!	TL_2 = TL_km_arr_ExtendedK(Rab, Lm, k2)
	!	
	!	sum_a(1:N_sp+2) = (0.0, 0.0)
	!	sum_b(1:N_sp+2) = (0.0, 0.0)
	!	sum_c(1:N_sp+2) = (0.0, 0.0)
	!	
	!	call RT_function_ExtendedK(struct, m, Ra, Int_fm)
	!	call RT_function_ExtendedK(struct, n, Rb, Int_fn)
	!	do tt = 1, 2 !Attention, sometimes, t and l might not be smaller than 2			
	!		Rf_vec = Int_fm%pair_kk(tt)%Rfk_vec
	!		Rf = Int_fm%pair_kk(tt)%Rfk
	!		R1(1:N_sp+2, 1:3) = Rf_vec(1:N_sp+2, 1, 1:3)
	!		R2(1:N_sp+2, 1:3) = Rf_vec(1:N_sp+2, 2, 1:3)	
 !
	!		do ll = 1, 2			
	!			Tf_vec = Int_fn%pair_kk(ll)%Rfk_vec
	!			Tf = Int_fn%pair_kk(ll)%Rfk				
	!			Tf_vec = conjg(Tf_vec)
	!			Tf = conjg(Tf)	
	!			T1(1:N_sp+2, 1:3) = Tf_vec(1:N_sp+2, 1, 1:3)
	!			T2(1:N_sp+2, 1:3) = Tf_vec(1:N_sp+2, 2, 1:3)					
	!			
	!			ss1 = cross_rc_arr(N_sp+2, k_hat_p, R1)*im*k1!
	!			ss2 = cross_rc_arr(N_sp+2, k_hat_p, R2)*im*k2!	
	!			
	!			sum_a = sum_a -  (TL_1*dot_product_arr(N_sp+2, ss1, T1) + &
	!										TL_2*dot_product_arr(N_sp+2, ss2, T2))*TRF_p !
	!			f1 = omega*dot_product_arr(N_sp+2, R1, T1)
	!			g1 = Rf(:, 1)*Tf(:, 1)/omega
	!			f2 = omega*dot_product_arr(N_sp+2, R2, T2)
	!			g2 = Rf(:, 2)*Tf(:, 2)/omega	
	!			sum_b = sum_b + (TL_1*(my_1*f1 - g1/eps_1) + TL_2*(my_2*f2 - g2/eps_2))*TRF_p*im
	!			sum_c = sum_c + (TL_1*(eps_1*f1 - g1/my_1) + TL_2*(eps_2*f2 - g2/my_2))*TRF_p*im	
	!		end do
	!	end do		
	!	suma = sum(sum_a)
	!	sumb = sum(sum_b)
	!	sumc = sum(sum_c)
	!end subroutine FMM_edge_Tri_RT_ExtendedK
	!
	
	!
	!subroutine RT_function_ExtendedK(struct, n, Ra, Int_fn)
	!! R_ab : distance between the two points
	!	implicit none		
	!	real(kind = dp), intent(in) :: Ra(3)		
	!	integer, intent(in) :: n  !sampling point on a sphere				
	!	integer :: i, j, pos_neg, ll, m_pairs
	!	real(kind = dp), dimension(100) :: w(100), a(100), b(100)
	!	real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roa(3), Ro(3), fm(3)
	!	complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
	!	real(kind = dp) :: area, Len_m, dfm		
	!	type(Structure), intent(in) :: struct	
	!	type(Integration_fn), intent(out) :: Int_fn
	!	m_pairs = size(struct%neighbours)		
	!			
	!	call Quadrature_tri(ng, a, b, w)
	!	do ll = 1, 2
	!		call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
	!		pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
	!		dfm = pos_neg*Len_m 			
	!		allocate(Int_fn%pair_kk(ll)%Rfk_vec(N_sp+2, 2, 3))
	!		allocate(Int_fn%pair_kk(ll)%Rfk(N_sp+2, 2))
	!		do i = 1, N_sp+2
	!			Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!			Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!			Int_fn%pair_kk(ll)%Rfk(i, 1:2) = (0.0, 0.0)						
	!			do j = 1, ng
	!				Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
	!				fm = 0.5*(Ro-vm_t)*dfm
	!				Roa = Ro - Ra
	!				tmp = dot_product(k_hat_p(i, :), Roa(:))
	!				ff1 = exp(-im*k1*tmp)*w(j)
	!				ff2 = exp(-im*k2*tmp)*w(j)
	!				Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) + fm*ff1
	!				Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) + fm*ff2
	!				Int_fn%pair_kk(ll)%Rfk(i, 1) = Int_fn%pair_kk(ll)%Rfk(i, 1) + dfm*ff1
	!				Int_fn%pair_kk(ll)%Rfk(i, 2) = Int_fn%pair_kk(ll)%Rfk(i, 2) + dfm*ff2
	!			end do
	!		end do !		
	!	end do
	!end subroutine RT_function_ExtendedK
	
	subroutine RT_function(struct, n, Ra, Int_fn)
	! R_ab : distance between the two points
		implicit none		
		real(kind = dp), intent(in) :: Ra(3)		
		integer, intent(in) :: n  !sampling point on a sphere				
		integer :: i, j, pos_neg, ll, m_pairs
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roa(3), Ro(3), fm(3)
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		type(Integration_fn), intent(out) :: Int_fn
		m_pairs = size(struct%neighbours)		
				
		call Quadrature_tri(ng, a, b, w)
		do ll = 1, 2
			call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
			pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
			dfm = pos_neg*Len_m 			
			allocate(Int_fn%pair_kk(ll)%Rfk_vec(N_sp, 2, 3))
			allocate(Int_fn%pair_kk(ll)%Rfk(N_sp, 2))
			do i = 1, N_sp
				Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				Int_fn%pair_kk(ll)%Rfk(i, 1:2) = (0.0, 0.0)						
				do j = 1, ng
					Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
					fm = 0.5*(Ro-vm_t)*dfm
					Roa = Ro - Ra
					tmp = dot_product(k_hat(i, :), Roa(:))
					ff1 = exp(-im*k1*tmp)*w(j)
					ff2 = exp(-im*k2*tmp)*w(j)
					Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) + fm*ff1
					Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) + fm*ff2
					Int_fn%pair_kk(ll)%Rfk(i, 1) = Int_fn%pair_kk(ll)%Rfk(i, 1) + dfm*ff1
					Int_fn%pair_kk(ll)%Rfk(i, 2) = Int_fn%pair_kk(ll)%Rfk(i, 2) + dfm*ff2
				end do
			end do !		
		end do
	end subroutine RT_function
	
	subroutine Khat_and_TRF_legendre(K_m)
		implicit none	
		real(kind = dp), dimension(:), allocatable :: xx, ww
		real(kind = dp) :: d_phi, phi, a, b, sum_a
		integer :: j, k, i!
		integer, intent(in) :: K_m
		N_sp = K_m*m_phi
		allocate(k_hat(N_sp, 3))
		allocate(xx(K_m))
		allocate(ww(K_m))
		allocate(TRF(N_sp))
		a = -1.0
		b = 1.0		
		call legendre_handle(K_m, a, b, xx, ww)		
		d_phi = 2*PI/m_phi		
		do j = 1, K_m
         do k = 1, m_phi
            phi = (k-1)*d_phi 
            i = k + (j-1)*m_phi
            k_hat(i, :) = (/sin(acos(xx(j)))*cos(phi), sin(acos(xx(j)))*sin(phi), xx(j)/)
				TRF(i) = ww(j)*d_phi
         end do
		end do
		deallocate(xx)
		deallocate(ww)
	end subroutine Khat_and_TRF_legendre
	
	subroutine spherical_harmonics_Darve(Lm, Ylm, x, phi)
		implicit none	
		integer(kind = 4), intent(in) :: Lm
		double precision, intent(in) :: phi, x
		type(list_list_real) :: pm, pd
		type(list_list_cmplx), intent(inout) :: Ylm
		integer(kind = 4) :: fnu, m, n, p
		logical :: cs_phase 
		real(kind = dp) :: tmp 
		fnu = 0		
		cs_phase = .false. 		
		call init_list(Ylm, fnu, Lm)	
		call lib_math_associated_legendre_polynomial_range(x, fnu, Lm, pm, pd, cs_phase)	
		pm = 1/sqrt(4*PI)*pm	
		do n = fnu, fnu + Lm-1
			do m = 0, n
				tmp = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n, m)
				Ylm%item(n)%item(m) = (-1)**m*sqrt((2*n+1)*tmp)*pm%item(n)%item(m)*exp(im*m*phi)
				Ylm%item(n)%item(-m) = conjg(Ylm%item(n)%item(m))*(-1)**m
			end do
		end do	
		deallocate(pm%item)
		return
	end subroutine 
	

	
	subroutine test_Ylm()
		implicit none
		type(list_list_cmplx) :: Ylm
		integer :: i, j, Ns, Lm, m, n
		real(kind = 8) :: theta, dtheta, dphi, phi
		complex(kind = 8) :: sum_a, sum_b
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)		
		Ns = 100
		dtheta = PI/100
		dphi = 2*PI/100
		Lm = 5		
		n = 2
		m = 1 
		theta = 0.50
		phi = 0.1
		call spherical_harmonics_Darve(Lm + 1, Ylm, cos(theta), phi)
		sum_a = Ylm%item(n)%item(0)
		!sum_b = 1/sqrt(4*PI)
		!sum_b = sqrt(3/(4*PI))*cos(theta)
		!sum_b = -sqrt(3/(8*PI))*sin(theta)*exp(im*phi) !Y1,+1
		sum_b = sqrt(5/(16*PI))*(3*cos(theta)*cos(theta) -1) !Y2,0
		!sum_b = -sqrt(15/(8*PI))*(sin(theta)*cos(theta)*exp(im*phi))!Y2,+1
		!sum_b = sqrt(35/(64*PI))*(sin(theta)*sin(theta)*sin(theta)*exp(-im*3*phi)) !Y3,-3
		print*, 'sum_a=', sum_a
		print*, 'sum_b=', sum_b
	end subroutine
	
	function C_lk(Lp, Lq)
		integer, intent(in) :: Lp, Lq
		real(kind = dp) :: C_lk(0:Lq)
		integer :: j
		do j = 0, Lq
			if (j <= Lp) then
				C_lk(j) = 1.0				
			else if ((j<= Lq) .and. (j>=Lp)) then
				c_lk(j) = cos((j-Lp)*PI/(2*Lp*(PI-1.0)))*cos((j-Lp)*PI/(2*Lp*(PI-1.0)))
			end if
		end do		
		return
	end function
	
	subroutine test_Plm()
		implicit none	
		integer(kind = 4) :: Lm
		double precision :: phi, x
		type(list_list_real) :: pm, pd
		integer(kind = 4) :: fnu, m, n, p, i, tmp_2
		logical :: cs_phase 
		real(kind = dp), dimension(:), allocatable :: xx, ww
		real(kind = dp) :: a, b, suma, tmp		
		
		real (kind = dp), dimension(:), allocatable :: pl, dummy 		
		Lm = 5
		a = -1.0
		b = 1.0		
		suma = 0.0
		call legendre_handle(Lm, a, b, xx, ww)				
		fnu = 1		
		cs_phase = .true. 		
		x = 0.9		
		call init_list(pm, fnu, Lm + 1)			
		do i = 1, Lm
			call lib_math_associated_legendre_polynomial_range(xx(i), fnu, Lm + 1, pm, pd, cs_phase)
		end do		
				
		tmp = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(3, 2) !this is the (1+m)!/(l-m)!
		tmp = tmp*pm%item(3)%item(2)
		
		allocate(pl(Lm),dummy(Lm))		
		do i = 1, Lm
			call lib_math_legendre_polynomial(xx(i), fnu, Lm, pl, dummy)			
			print*, 'i = ', i
			print*, 'xx(i) =', xx(i)			
			print*, 'pl(i) =', pl(i)			
		end do
		print*, 'Analytical p1=', xx(1)
		print*, 'Analytical p2=', 0.5*(3*xx(2)**2-1)
		print*, 'Analytical p5=', 1.0/8.0*(63*xx(5)**5-70*(xx(5))**3 + 15*xx(5))	
	end subroutine 
	
	subroutine test_Ylm_II()
		implicit none
		type(list_list_cmplx) :: Ylm
		integer :: i, j, Ns, Lm, m, n
		real(kind = 8) :: theta, dtheta, dphi, phi, a, b
		real(kind = dp), dimension(:), allocatable :: xx_p, ww_p
		real(kind = dp), dimension(:), allocatable :: xx_t, ww_t
		complex(kind = 8) :: sum_a, sum_b
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)				
		dtheta = PI/100
		dphi = 2*PI/100
		
		Lm = 11
		n = Lm -2
		m = (n-1)
		
		allocate(xx_t(Lm+1), ww_t(Lm+1))
		allocate(xx_p(Lm+1), ww_p(Lm+1))
		a = -1
		b = 1
		call legendre_handle(Lm+1, a, b, xx_t, ww_t)		
		a = 0
		b = 2*PI		
		dphi = 2*PI/Ns
		
		call legendre_handle(Ns, a, b, xx_p, ww_p)		
		do i = 1, Lm+1
			do j = 1, Ns	
				phi = xx_p(j)
				dphi = ww_p(j)
				call spherical_harmonics_Darve(Lm+1, Ylm, xx_t(i), phi)
				sum_a = sum_a + Ylm%item(n)%item(m)*conjg(Ylm%item(n)%item(m))*ww_t(i)*ww_p(j)
				sum_b = sum_b + ww_t(i)*ww_p(j)
			end do
		end do
		print*, 'sum_a=', sum_a
		print*, 'sum_b=', sum_b
	end subroutine
	
	subroutine test_Ylm_I()
		implicit none
		type(list_list_cmplx) :: Ylm
		integer :: i, j, Ns, Lm
		real(kind = 8) :: theta, dtheta, dphi, phi
		complex(kind = 8) :: sum_a, sum_b
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)		
		Ns = 100
		dtheta = PI/100
		dphi = 2*PI/100
		Lm = 15		
		do i = 1, Ns
			theta = (i-1)*dtheta
			do j = 1, Ns
				phi = (j-1)*dphi
				call spherical_harmonics_Darve(Lm + 1, Ylm, cos(theta), phi)
				sum_a = sum_a + Ylm%item(15)%item(-11)*conjg(Ylm%item(15)%item(-11))*sin(theta)*dtheta*dphi
				sum_b = sum_b + sin(theta)*dtheta*dphi
			end do
		end do
		print*, 'sum_a=', sum_a
		print*, 'sum_b=', sum_b
	end subroutine	

	subroutine Khat_and_TRF_legendre_new(K_m)
		implicit none	
		real(kind = dp), dimension(:, :), allocatable :: xw
		real(kind = dp) :: d_phi, phi, a, b, sum_a
		integer :: j!, m_phi	
        integer, intent(in) :: K_m
		!m_phi = 2*K_m+1
		N_sp = K_m*m_phi
		allocate(k_hat(N_sp, 3))
		allocate(TRF(N_sp))		
		xw = theta_phi_w(K_m)
		do j = 1, N_sp         
         k_hat(j, :) = (/sin(acos(xw(j, 1)))*cos(xw(j, 2)), sin(acos(xw(j, 1)))*sin(xw(j, 2)), xw(j, 1)/)
			TRF(j) = xw(j, 3)
		end do		
		deallocate(xw)
	end subroutine Khat_and_TRF_legendre_new	
	
	subroutine FMM_edge_Tri_RT_kd(m, n, struct, struct_cube, d_c, suma, sumb, sumc) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), dimension(N_sp) :: sum_a, sum_c, sum_b!, intent(out)
		!complex(kind = dp), dimension(N_sp) :: tmp_1, tmp_2
		complex(kind = dp), dimension(N_sp) :: Rf1, Tf1, Rf2, Tf2
		complex(kind = dp), dimension(N_sp) :: TL_1, TL_2
		complex(kind = dp), dimension(N_sp, 3) :: ss1, ss2, R1, R2, T1, T2
		complex(kind = dp), dimension(N_sp, 3) :: Rf1_vec, Tf1_vec, Rf2_vec, Tf2_vec
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		complex(kind = dp), dimension(N_sp) :: f1, g1, f2, g2
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
		type(Integration_fn_II) :: Int_fn1, Int_fn2, Int_fm1, Int_fm2
		!-------------for test use-----------------
		!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
		!       Rq
		!       o
		!       | 
		!	     |
		!	  Ra o-----------------------o Rb
		!                                |
		!                                |
		!                                o Rp
		!------------------------------------------
		
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
		intrinsic :: sqrt	
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		sum_a(1:N_sp) = (0.0, 0.0)
		sum_b(1:N_sp) = (0.0, 0.0)
		sum_c(1:N_sp) = (0.0, 0.0)
		
		call RT_function_kd(m, struct, k1, Ra, Int_fm1)
		call RT_function_kd(m, struct, k2, Ra, Int_fm2)
		call RT_function_kd(n, struct, k1, Rb, Int_fn1)
		call RT_function_kd(n, struct, k2, Rb, Int_fn2)
		
		do tt = 1, 2 !Attention, sometimes, t and l might be smaller than 2			
			Rf1_vec = Int_fm1%pair_ks(tt)%Rfs_vec
			Rf1 = Int_fm1%pair_ks(tt)%Rfs
			!print*, 'Rf1_vec =', Rf1_vec(1, :)
			!print*, 'Rf2_vec =', Rf2_vec(1, :)
			
			Rf2_vec = Int_fm2%pair_ks(tt)%Rfs_vec
			Rf2 = Int_fm2%pair_ks(tt)%Rfs			
			do ll = 1, 2
				Tf1_vec = Int_fn1%pair_ks(ll)%Rfs_vec
				Tf1 = Int_fn1%pair_ks(ll)%Rfs
				Tf1_vec = conjg(Tf1_vec)
				Tf1 = conjg(Tf1)
				
				Tf2_vec = Int_fn2%pair_ks(ll)%Rfs_vec
				Tf2 = Int_fn2%pair_ks(ll)%Rfs
				Tf2_vec = conjg(Tf2_vec)
				Tf2 = conjg(Tf2)
				ss1 = cross_rc_arr(N_sp, k_hat, Rf1_vec)*im*k1!
				ss2 = cross_rc_arr(N_sp, k_hat, Rf2_vec)*im*k2!						
				sum_a = sum_a - (TL_1*dot_product_arr(N_sp, ss1, Tf1_vec) + &
				TL_2*dot_product_arr(N_sp, ss2, Tf2_vec))*TRF
				
				f1 = omega*dot_product_arr(N_sp, Rf1_vec, Tf1_vec)
				g1 = Rf1*Tf1/omega
				f2 = omega*dot_product_arr(N_sp, Rf2_vec, Tf2_vec)
				g2 = Rf2*Tf2/omega	
				sum_b = sum_b + ((TL_1*(my_1*f1 - g1/eps_1) + TL_2*(my_2*f2 - g2/eps_2))*TRF*im)
				sum_c = sum_c + ((TL_1*(eps_1*f1 - g1/my_1) + TL_2*(eps_2*f2 - g2/my_2))*TRF*im)	
			end do
		end do
		suma = sum(sum_a)
		sumb = sum(sum_b)
		sumc = sum(sum_c)
	end subroutine FMM_edge_Tri_RT_kd
	
	subroutine FMM_cube_Tri_arr(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, rr) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n_cube, m_pairs		
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp) :: sum_a, sum_c, sum_b
		character(len = 1), intent(in) :: transa
		
		complex(kind = dp) :: sum_r(2), sum_t(2), Term_I, Term_II
		complex(kind = dp),  intent(out) :: rr(2)
		complex(kind = dp), dimension(:, :), allocatable :: ss1, ss2
		complex(kind = dp), dimension(:), allocatable :: f1, g1, f2, g2 		
		complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2, x1, tmp_1, tmp_2
 
		complex(kind = dp), dimension(:, :), allocatable :: Rfk, Tfk, R1, R2, T1, T2 
		complex(kind = dp), dimension(:, :, :), allocatable :: Rfk_vec, Tfk_vec
		real(kind = dp) ::  Ra(3), Rb(3), Rm(3), Rn(3), R(3), Rab(3)
		
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		integer :: j, m_cube, tt, ll, nn_edge, nd, n_edge_number, Lm
		
		allocate(Rfk_vec(N_sp, 2, 3), Tfk_vec(N_sp, 2, 3)) 
		allocate(TL_1(N_sp), TL_2(N_sp), Rfk(N_sp, 2), Tfk(N_sp, 2))
		allocate(tmp_1(N_sp), tmp_2(N_sp)) 
		allocate(ss1(N_sp, 3), ss2(N_sp, 3), T1(N_sp, 3))
		allocate(T2(N_sp, 3), R1(N_sp, 3), R2(N_sp, 3))
		allocate(f1(N_sp), f2(N_sp), g2(N_sp), g1(N_sp)) 
		
		sum_r(1:2) = (0.0, 0.0)
		sum_t(1:2) = (0.0, 0.0)
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located			
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
	   Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube		
		do nn_edge = 1, n_edge_number		
			nd = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n				
			sum_a = (0.0, 0.0)
			sum_b = (0.0, 0.0)
			sum_c = (0.0, 0.0)
			do tt = 1, 2 !
				call Receivef_PreC(struct, m, tt, Rfk_vec, Rfk)
				Rm = struct%midpoints(struct%neighbours(m)%elements(tt))%point
				do ll = 1, 2
					call Receivef_PreC(struct, nd, ll, Tfk_vec, Tfk)
					Rn = struct%midpoints(struct%neighbours(nd)%elements(ll))%point							
					Tfk_vec = conjg(Tfk_vec)
					Tfk = conjg(Tfk)
					R1(1:N_sp, 1:3) = Rfk_vec(1:N_sp, 1, 1:3)
					R2(1:N_sp, 1:3) = Rfk_vec(1:N_sp, 2, 1:3)
					T1(1:N_sp, 1:3) = Tfk_vec(1:N_sp, 1, 1:3)
					T2(1:N_sp, 1:3) = Tfk_vec(1:N_sp, 2, 1:3)
					
					tmp_1 = exp(im*k1*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab))) !(1.0, 0.0)!
					tmp_2 = exp(im*k2*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab)))!(1.0, 0.0)!
					ss1 = cross_rc_arr(N_sp, k_hat, R1)!*im*k1*tmp_1
					ss2 = cross_rc_arr(N_sp, k_hat, R2)!*im*k2*tmp_2
					do j = 1, N_sp
						ss1(j, :) = ss1(j, :)*tmp_1(j)*im*k1
						ss2(j, :) = ss2(j, :)*tmp_2(j)*im*k2
					end do
					sum_a = sum_a - sum((TL_1*dot_product_arr(N_sp, ss1, T1) + &
					TL_2*dot_product_arr(N_sp, ss2, T2))*TRF)
						
					f1 = omega*dot_product_arr(N_sp, R1, T1)
					g1 = Rfk(:, 1)*Tfk(:, 1)/omega
					f2 = omega*dot_product_arr(N_sp, R2, T2)
					g2 = Rfk(:, 2)*Tfk(:, 2)/omega	
					sum_b = sum_b + sum((TL_1*(my_1*f1 - g1/eps_1)*tmp_1 + TL_2*(my_2*f2 - g2/eps_2)*tmp_2)*TRF*im)
					sum_c = sum_c + sum((TL_1*(eps_1*f1 - g1/my_1)*tmp_1 + TL_2*(eps_2*f2 - g2/my_2)*tmp_2)*TRF*im)					
				end do
			end do
			sum_r(1) = sum_r(1) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs)
			sum_r(2) = sum_r(2) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
			sum_t(1) = sum_t(1) + sum_b*x1(nd) + sum_a*x1(nd + m_pairs)
			sum_t(2) = sum_t(2) - sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
		end do		
		if (transa == 'T') then
			rr = sum_t
		else if (transa == 'N') then
			rr = sum_r
		else 
			print*, 'Not a proper multiplication type!'
		end if		
	end subroutine FMM_cube_Tri_arr	
	
	!subroutine spherical_harmonics_Gibson(l, m, Ylm, x, phi)
	!	implicit none	
	!	integer(kind = 4), intent(in) :: m, l
	!	double precision, intent(in) :: phi, x
	!	double precision, dimension(:, :), allocatable :: pm, pd
	!	double precision :: y
	!	complex(kind = dp), dimension(:), allocatable, intent(out) :: Ylm
	!	integer(kind = 4) :: fnu, n
	!	logical :: cs_phase 
	!	allocate(Ylm(l))
	!	fnu = 1		
	!	cs_phase = .false. 		
	!	n = abs(m)
	!	y = abs(x)
	!	call  lib_math_associated_legendre_polynomial_with_negative_m(x, n, fnu, l, pm, pd, cs_phase)
	!	
	!	if (m>=0) then						
	!		Ylm(:) = 1/sqrt(4*PI)*pm(2, :)*exp(im*m*phi) 		
	!	else			
	!		Ylm(:) = 1/sqrt(4*PI)*pm(1, :)*exp(im*m*phi) 
	!	end if
	!	print*, ''		
	!	return
	!end subroutine 
		
	
	subroutine RT_function_kd(n, struct, kd, Ra, Int_fn)
	! R_ab : distance between the two points
		implicit none		
		real(kind = dp), intent(in) :: Ra(3)	
		complex(kind = dp), intent(in) :: kd !complex wave number 
		integer, intent(in) :: n  !edge number
		integer :: i, j, pos_neg, ll, m_pairs
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roa(3), Ro(3), fm(3)
		complex(kind = dp) :: ff, tmp
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		type(Integration_fn_II), intent(out) :: Int_fn
		
		m_pairs = size(struct%neighbours)		
				
		call Quadrature_tri(ng, a, b, w)
		do ll = 1, 2
			call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
			pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
			dfm = pos_neg*Len_m 			
			allocate(Int_fn%pair_ks(ll)%Rfs_vec(N_sp, 3))
			allocate(Int_fn%pair_ks(ll)%Rfs(N_sp))
			do i = 1, N_sp
				Int_fn%pair_ks(ll)%Rfs_vec(i, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
				Int_fn%pair_ks(ll)%Rfs(i) = (0.0, 0.0)
				do j = 1, ng
					Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
					fm = 0.5*(Ro-vm_t)*dfm
					Roa = Ro - Ra
					tmp = dot_product(k_hat(i, :), Roa(:))
					ff = exp(-im*kd*tmp)*w(j)
					Int_fn%pair_ks(ll)%Rfs_vec(i, :) = Int_fn%pair_ks(ll)%Rfs_vec(i, :) + fm*ff
					Int_fn%pair_ks(ll)%Rfs(i) = Int_fn%pair_ks(ll)%Rfs(i) + dfm*ff
				end do
			end do !		
		end do
	end subroutine RT_function_kd
	
	!subroutine RT_function(struct, n, Ra, Int_fn)
	!! R_ab : distance between the two points
	!	implicit none		
	!	real(kind = dp), intent(in) :: Ra(3)		
	!	integer, intent(in) :: n  !sampling point on a sphere				
	!	integer :: i, j, pos_neg, ll, m_pairs
	!	real(kind = dp), dimension(100) :: w(100), a(100), b(100)
	!	real(kind = dp) :: pm1(3), pm2(3), vm_t(3), Roa(3), Ro(3), fm(3)
	!	complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
	!	real(kind = dp) :: area, Len_m, dfm		
	!	type(Structure), intent(in) :: struct	
	!	type(Integration_fn), intent(out) :: Int_fn
	!	m_pairs = size(struct%neighbours)		
	!			
	!	call Quadrature_tri(ng, a, b, w)
	!	do ll = 1, 2
	!		call fn_parameter(struct, n, ll, pm1, pm2, vm_t, Len_m, area)
	!		pos_neg = -1*(((ll-1)*ll)-1)!put it outside this routine
	!		dfm = pos_neg*Len_m 			
	!		allocate(Int_fn%pair_kk(ll)%Rfk_vec(N_sp, 2, 3))
	!		allocate(Int_fn%pair_kk(ll)%Rfk(N_sp, 2))
	!		do i = 1, N_sp
	!			Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!			Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!			Int_fn%pair_kk(ll)%Rfk(i, 1:2) = (0.0, 0.0)						
	!			do j = 1, ng
	!				Ro = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
	!				fm = 0.5*(Ro-vm_t)*dfm
	!				Roa = Ro - Ra
	!				tmp = dot_product(k_hat(i, :), Roa(:))
	!				ff1 = exp(-im*k1*tmp)*w(j)
	!				ff2 = exp(-im*k2*tmp)*w(j)
	!				Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 1, :) + fm*ff1
	!				Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) = Int_fn%pair_kk(ll)%Rfk_vec(i, 2, :) + fm*ff2
	!				Int_fn%pair_kk(ll)%Rfk(i, 1) = Int_fn%pair_kk(ll)%Rfk(i, 1) + dfm*ff1
	!				Int_fn%pair_kk(ll)%Rfk(i, 2) = Int_fn%pair_kk(ll)%Rfk(i, 2) + dfm*ff2
	!			end do
	!		end do !		
	!	end do
	!end subroutine RT_function
		
	!function TL_km_arr_ExtendedK(r_ab, m, K_m)
 !     integer, intent(in) :: m    !Lm is the order limit of the polynomials, N is the sampling points on the sphere
 !     real(kind = dp), intent(in) :: r_ab(3)
 !     real(kind = dp) :: r_hat(3), x, pm_arr(N_sp+2, m)
 !     complex(kind = dp) :: r_h
	!	!real(kind = dp) :: k_hat(N_sp, 3)
 !     complex(kind = dp), intent(in) :: K_m
 !     real (kind = dp), dimension(m) :: pm, dummy 
 !     complex(kind = dp) :: TL_km_arr_ExtendedK(N_sp+2)
 !     complex(kind = dp), dimension(m) :: hl 
 !     integer :: j, fnu
	!	
 !     fnu = 0        
 !     r_hat = r_ab/vec_len(r_ab)
 !     r_h =  vec_len(r_ab)*K_m  
 !     hl =  lib_math_hankel_spherical_2(r_h, fnu, m) 
 !     TL_km_arr_ExtendedK(1:N_sp+2) = (0.0, 0.0)     
	!	
 !     do j = 1, N_sp+2
 !        x = dot_product(r_hat(:), k_hat_p(j, :))
 !        call lib_math_legendre_polynomial(x, fnu, m, pm, dummy)
 !        pm_arr(j, :) = pm(:)
 !     end do
	!	
 !     do j = 1, m
 !        TL_km_arr_ExtendedK(:) = TL_km_arr_ExtendedK(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*K_m/(4*PI)**2 ! 
	!	end do
 !     return        
	!end function TL_km_arr_ExtendedK
	!
 !
	!
	!function Rf_rc(r_oa, kd) 
	!	implicit none
	!	integer :: j
	!	complex(kind = 8), dimension(N_sp) :: Rf_rc
	!	real(kind = dp) :: r_oa(3)
	!	complex(kind = dp) :: kd
	!	Rf_rc(:) = exp(-im*kd*(k_hat(:, 1)*r_oa(1) + k_hat(:, 2)*r_oa(2) + k_hat(:, 3)*r_oa(3)))*TRF(:)
	!	print*, ''
	!end function	
	!
	!function Tf_rd(r_pb, kd) 
	!	implicit none
	!	integer :: j
	!	complex(kind = 8), dimension(N_sp) :: Tf_rd
	!	real(kind = dp) :: r_pb(3)
	!	complex(kind = dp) :: kd
	!	
	!	Tf_rd(:) = exp(im*kd*(k_hat(:, 1)*r_pb(1)+ k_hat(:, 2)*r_pb(2) + k_hat(:, 3)*r_pb(3)))		
	!end function	
	!
	!!Implemented according to Gibson's book "The Method of Moments in Electromagnetics"	
	subroutine sort_edge_index(vector, v_size, v_ascendente)
		integer, intent(in) :: v_size
		integer, intent(in) :: vector(v_size)
		integer, intent(out) :: v_ascendente(v_size)
		integer :: tmp(v_size)
		integer ::  i
		
		LOGICAL, DIMENSION(v_size) :: mk
		mk  = .TRUE.
		!vector = (/1, 3, 5, 4, 2, 8, 10, 9, 7, 6/)		 
		DO i = 1, v_size
			tmp(i) = MINVAL(vector,mk)
			mk(MINLOC(vector,mk)) = .FALSE.
		END DO
		v_ascendente = tmp 
		return
		!print*, 'v_ascendente =', v_ascendente
	end subroutine
	!
	!subroutine cube_definition(struct, struct_cube_II, d_box) !, cube_list_level_I
	!	!First test for spheres
	!	implicit none
	!	type(Structure_cube) :: struct_cube_I, struct_cube_tmp
	!	type(Structure_cube), intent(out) :: struct_cube_II 
	!	type(Structure), intent(in) :: struct
	!	integer :: n_pairs, m, n_tmp
	!	integer :: cube_list_level_I, ot, N, tmp_far(1200), tmp_near(1200), tmp_arr(1200)
	!	integer :: i, j, k, Nx, Ny, Nz, n_counter, m_counter, p_counter
	!	real(kind = dp) :: p_edgec(3), r_in, r_out, W_out, r(3), W_box, p2(3)
	!	real(kind = dp), intent(out) :: d_box
	!	real(kind = dp) :: x_min, x_max, y_min, y_max, z_min, z_max, dc
	!	intrinsic sqrt
	!	logical :: con1, con2, con3
	!	
	!	!call Data_Import(struct)
	!	n_pairs = size(struct%neighbours)	
	!	
	!	
	!	W_out = D_factor*1.1 !D_factor: the radius of the phere 			
	!	!print*, 'n_pairs =', n_pairs
	!	Nx = 8 !Temporary solution
	!	N = Nx*Nx*Nx!ceiling(W_out/W_box) !only two level		
	!	print*, 'Number of total cubes', N
	!	n_tmp = 0
	!	W_box = W_out/Nx
	!	allocate(Struct_cube_I%cubes(N))
	!	allocate(Struct_cube_tmp%cubes(N))
	!	allocate(Struct_cube_II%pair_cubes(n_pairs))
	!	do i = 1, Nx			
	!		do j = 1, Nx
	!			do k = 1, Nx
	!			    m = ((i-1)*Nx + (j-1))*Nx + k					
	!				Struct_cube_I%cubes(m)%cube_position = &
	!				(/-W_out/2 + (i - 0.5)*W_box, -W_out/2 + (j - 0.5)*W_box, -W_out/2 + (k - 0.5)*W_box/) !					
	!			end do
	!		end do
	!	end do
	!	!		
	!	ot = 1		
	!	d_box = sqrt(3.0)*W_box
	!	print*, 'd_box =', d_box
	!	n_counter = 0
	!	m_counter = 1
	!	
	!	do k = 1, N
	!		p2 = struct_cube_I%cubes(k)%cube_position
	!		x_min = p2(1) - W_box/2
	!		x_max = p2(1) + W_box/2
	!		y_min = p2(2) - W_box/2
	!		y_max = p2(2) + W_box/2
	!		z_min = p2(3) - W_box/2
	!		z_max = p2(3) + W_box/2
	!		do m = 1, n_pairs
	!			call fn_edge_center(struct, m, ot, p_edgec)	!find the edge center of an element pair
	!			con1 = (x_min < p_edgec(1)) .and. (p_edgec(1) <= x_max)
	!			con2 = (y_min < p_edgec(2)) .and. (p_edgec(2) <= y_max)
	!			con3 = (z_min < p_edgec(3)) .and. (p_edgec(3) <= z_max)
	!			
	!			if ( con1 .and. con2 .and. con3) then
	!				n_counter = n_counter + 1 !number of edges in one cube						
	!				tmp_arr(n_counter) = m !edge index saved	
	!				Struct_cube_II%pair_cubes(m)%cube_index = m_counter !edge 				
	!			end if				
	!		end do
	!		if (n_counter .ne. 0) then				
	!			struct_cube_tmp%cubes(m_counter)%cube_position = p2					
	!			allocate(struct_cube_tmp%cubes(m_counter)%edges_in_cube(n_counter))
	!			struct_cube_tmp%cubes(m_counter)%edges_in_cube(1:n_counter) = tmp_arr(1:n_counter)				
	!			m_counter = m_counter + 1 !index of cubes				
	!		end if
	!		n_counter = 0
	!	end do
	!	!Copy the data in tmp to struct_cube_II
	!	m_counter = m_counter - 1
	!	allocate(struct_cube_II%cubes(m_counter))		
	!	do i = 1, m_counter
	!		!print*, 'box number =', i
	!	    m = size(struct_cube_tmp%cubes(i)%edges_in_cube)
	!		 !print*, 'edges in the box =', m
	!		allocate(struct_cube_II%cubes(i)%edges_in_cube(m))
	!		struct_cube_II%cubes(i)%edges_in_cube(1:m) = struct_cube_tmp%cubes(i)%edges_in_cube(1:m)
	!		struct_cube_II%cubes(i)%cube_position = struct_cube_tmp%cubes(i)%cube_position
	!		n_tmp = n_tmp + m
	!	end do
	!	print*, 'n_tmp =', n_tmp
	!	
	!	!-------finde the cubes in the near regions and save the cube indices to cube_near
	!	n_counter = 0
	!	p_counter = 0
	!	do i = 1, m_counter
	!		do j = 1, m_counter
	!			if (vec_len(struct_cube_tmp%cubes(i)%cube_position - struct_cube_tmp%cubes(j)%cube_position) <= 1.5*d_box) then				
	!				n_counter = n_counter + 1
	!				tmp_near(n_counter) = j					
	!			else
	!				p_counter = p_counter + 1
	!				tmp_far(p_counter) = j
	!			end if					
	!		end do
	!		allocate(struct_cube_II%cubes(i)%cubes_near(n_counter))
	!		allocate(struct_cube_II%cubes(i)%cubes_far(p_counter))
	!		struct_cube_II%cubes(i)%cubes_near(1:n_counter) = tmp_near(1:n_counter)
	!		struct_cube_II%cubes(i)%cubes_far(1:p_counter) = tmp_far(1:p_counter)
	!		n_counter = 0
	!		p_counter = 0
	!	end do	
	!end subroutine cube_definition
	
	subroutine FMM_Near_MVProduct_Tri(m_pairs, struct, struct_cube, d_c, aa_csr, ib, jb, D_mat)
	!Near field matrix in a csr format
		integer :: m_cube, n_cube, m, n_edge_number, nn_edge, mm_edge, nd
		integer :: m_pairs, n_cube_number, m_cube_number, m_edge_number
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp) :: sum_a, sum_c, sum_b		
		complex(kind = dp) :: ssum_a, ssum_c, ssum_b				
		integer, dimension(:), allocatable :: n_arr_tmp
		complex(kind = dp), dimension(:), allocatable :: a_tmp, a_csr, c_csr, b_csr
		complex(kind = dp), dimension(:), allocatable, intent(out) :: aa_csr
		complex(kind = dp), dimension(:, :), allocatable, intent(out) :: D_mat
		integer, dimension(:), allocatable :: ia, ja
		integer, dimension(:), allocatable, intent(out) :: ib, jb
		
		real(kind = dp), intent(in) :: d_c ! 
		real(kind = dp) :: r_ob, r_ref
		integer :: p, q, ngp, t, s, n, j, n_csr, m_csr!, counter 
		
		!---------Calculate the number of non-zero elements in the array a_st---------------		
		!counter = 0
		m_cube_number = size(struct_cube%cubes) 
		!print*, 'm_cube_number =', m_cube_number
		n_csr = 0
		do m = 1, m_cube_number
			m_edge_number = size(struct_cube%cubes(m)%edges_in_cube)
			n_cube_number = size(struct_cube%cubes(m)%cubes_near) !number of cubes in the near region of cube i
			do mm_edge = 1, m_edge_number	
				do n = 1, n_cube_number
					n_cube = struct_cube%cubes(m)%cubes_near(n)
					n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube)
					n_csr = n_csr + n_edge_number
				end do			
			end do
		end do
		allocate(a_csr(n_csr))
		allocate(b_csr(n_csr))
		allocate(c_csr(n_csr))
		allocate(a_tmp(n_csr))
		allocate(ja(n_csr))	
		allocate(ia(m_pairs + 1))
		allocate(D_mat(2*m_pairs, 2*m_pairs))
		allocate(ib(2*m_pairs + 1))
		allocate(jb(4*n_csr))
		allocate(aa_csr(4*n_csr))	
		allocate(n_arr_tmp(2*m_pairs))
		!
		ia(1:m_pairs+1) = 0
		ja(1:n_csr) = 0
		n_csr = 0
		m_csr = 0
		ia(1) = 1
		D_mat(2*m_pairs, 2*m_pairs) = (0.0, 0.0)
		do m = 1, m_pairs	
			n_csr = 0
			m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
			n_cube_number = size(struct_cube%cubes(m_cube)%cubes_near) 		
			!print*, 'near n_cube_number =', n_cube_number
			do n = 1, n_cube_number
				n_cube = struct_cube%cubes(m_cube)%cubes_near(n) !cube index in the near 
				n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube
				do nn_edge = 1, n_edge_number
					n_csr = n_csr + 1
					n_arr_tmp(n_csr) = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n					
				end do
			end do
			call sort_edge_index(n_arr_tmp(1:n_csr), n_csr, n_arr_tmp(1:n_csr))		
			ja(m_csr + 1 : m_csr + n_csr) = n_arr_tmp(1: n_csr)		
		
			p = struct%neighbours(m)%elements(1)
			q = struct%neighbours(m)%elements(2)	
			r_ref = vec_len(struct%midpoints(p)%point - struct%midpoints(q)%point)		
			do j = 1, n_csr
				n = n_arr_tmp(j)
				
				sum_a = (0.0, 0.0) 
				sum_b = (0.0, 0.0) 
				sum_c = (0.0, 0.0) 
				do t = 1, 2
					do s = 1, 2!
						p = struct%neighbours(m)%elements(t)					
						q = struct%neighbours(n)%elements(s)						
						r_ob = vec_len(struct%midpoints(p)%point - struct%midpoints(q)%point)					
						if (r_ob > r_ref) then
							ngp = 3				
							call Normal_Integration_new(m, n, t, s, ngp, struct, ssum_a, ssum_b, ssum_c)								
						else
							ngp = 3	
							call Singular_Integration(m, n, t, s, ngp, struct, ssum_a, ssum_b, ssum_c)
						end if
						sum_a = sum_a + ssum_a !
						sum_b = sum_b + ssum_b
						sum_c = sum_c + ssum_c
					end do !loop s 
				end do  !loop t
				a_csr(m_csr + j) = sum_a !one-based indexing
				b_csr(m_csr + j) = sum_b
				c_csr(m_csr + j) = sum_c
				D_mat(m, n) = sum_b
				D_mat(m, n + m_pairs) = -sum_a
				D_mat(m + m_pairs, n) = sum_a
				D_mat(m + m_pairs, n + m_pairs) = sum_c
			end do
			m_csr = m_csr + n_csr
			ia(m+1) = ia(m) + n_csr
			n_csr = 0
		end do !loop m		
		call Expension_A_CSR(a_csr, b_csr, c_csr, ia, ja, aa_csr, ib, jb)
	end subroutine FMM_Near_MVProduct_Tri
	
	subroutine Expension_A_CSR(a_csr, b_csr, c_csr, ia, ja, aa, ib, jb)
		integer :: n, m_csr, m, m_pairs, q
		complex(kind = 8), dimension(:), intent(out) :: aa
		complex(kind = 8), dimension(:), intent(in) :: a_csr, b_csr, c_csr
		integer, dimension(:), intent(in) :: ia, ja
		integer, dimension(:), allocatable, intent(out) :: ib, jb
		
		m_pairs = size(ia)-1
		m_csr = size(a_csr)
		allocate(ib(2*m_pairs + 1))
		allocate(jb(4*m_csr))	
		
		ib(1) = ia(1)
		do n = 1, m_pairs
			ib(n+1) = ib(n) + 2*(ia(n+1) - ia(n))
			m = ia(n+1)-ia(n)
			aa(ib(n): ib(n+1)-1-m) = b_csr(ia(n): ia(n+1)-1)
			jb(ib(n): ib(n+1)-1-m) = ja(ia(n): ia(n+1)-1)
			aa(ib(n+1)-m : ib(n+1)-1) = -a_csr(ia(n) : (ia(n+1)-1))
			jb(ib(n+1)-m : ib(n+1)-1) = ja(ia(n): ia(n+1)-1) + m_pairs
		end do
		ib(1+2*m_pairs) = ia(1) + 2*m_pairs
		do n = 1, m_pairs	  	 	   
			q = n + m_pairs
			ib(q+1) = ib(q) + 2*(ia(n+1) - ia(n))
			m = ia(n+1)-ia(n)
			aa(ib(q): ib(q+1)-1-m) = a_csr(ia(n): ia(n+1)-1)
			jb(ib(q): ib(q+1)-1-m) = ja(ia(n): ia(n+1)-1)
			aa(ib(q+1)-m : ib(q+1)-1) = c_csr(ia(n) : (ia(n+1)-1))
			jb(ib(q+1)-m : ib(q+1)-1) = ja(ia(n): ia(n+1)-1) + m_pairs
		end do
	end subroutine
	
	subroutine FMM_Far_MVProduct_Tri_OMP(transa, m_pairs, struct, struct_cube, d_c, x1, rr)
	!For iterative solution	
		integer, intent(in) :: m_pairs
		character(len = 1), intent(in) :: transa
		integer :: n_cube, m_cube, n, n_edge_number, nn_edge, nd, n_cube_number, m
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube		
		complex(kind = dp):: sum_r(2)
		complex(kind = dp), dimension(:), allocatable, intent(out) :: rr
		complex(kind = dp), dimension(:), allocatable, intent(in) :: x1	
		real(kind = dp), intent(in) :: d_c !
		allocate(rr(2*m_pairs))
		
		call timestamp()
		!$omp parallel private (m, sum_r, m_cube, n_cube, n, n_cube_number) &
		!$omp shared (m_pairs, struct, struct_cube, x1, d_c, rr, transa)
		!$omp do
		!j = omp_get_thread_num ( )
		!write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', j
		do m = 1, m_pairs
			rr(m) = (0.0, 0.0)			
			rr(m + m_pairs) = (0.0, 0.0)
			m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
			n_cube_number = size(struct_cube%cubes(m_cube)%cubes_far)	
			do n = 1, n_cube_number 
				n_cube = struct_cube%cubes(m_cube)%cubes_far(n) !cube index in the far region of the cube m_cube
				!call FMM_cube_Tri_II(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, sum_r)
				call FMM_cube_Tri_arr_PreCal(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, sum_r) 
				rr(m) = rr(m) +  sum_r(1)
				rr(m + m_pairs) = rr(m + m_pairs) +  sum_r(2)
			end do
		end do
		!$omp end do
		!$omp end parallel			
		call timestamp()
		print*, 'The new method'
		return		
	end subroutine FMM_Far_MVProduct_Tri_OMP
		
	subroutine FMM_cube_Tri_arr_PreCal(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, rr) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n_cube, m_pairs		
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp) :: sum_a, sum_c, sum_b
		character(len = 1), intent(in) :: transa
		
		complex(kind = dp) :: sum_r(2), sum_t(2), Term_I, Term_II
		complex(kind = dp),  intent(out) :: rr(2)
		complex(kind = dp), dimension(:, :), allocatable :: ss1, ss2
		complex(kind = dp), dimension(:), allocatable :: f1, g1, f2, g2 		
		complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2, x1, tmp_1, tmp_2
 
		complex(kind = dp), dimension(:, :), allocatable :: Rf, Tf, R1, R2, T1, T2 
		complex(kind = dp), dimension(:, :, :), allocatable :: Rf_vec, Tf_vec
		real(kind = dp) ::  Ra(3), Rb(3), Rm(3), Rn(3), R(3), Rab(3)
		
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level				
		integer :: j, m_cube, tt, ll, nn_edge, nd, n_edge_number, Lm
		
		allocate(Rf_vec(N_sp, 2, 3), Tf_vec(N_sp, 2, 3)) 
		allocate(TL_1(N_sp), TL_2(N_sp), Rf(N_sp, 2), Tf(N_sp, 2))
		allocate(tmp_1(N_sp), tmp_2(N_sp)) 
		allocate(ss1(N_sp, 3), ss2(N_sp, 3), T1(N_sp, 3))
		allocate(T2(N_sp, 3), R1(N_sp, 3), R2(N_sp, 3))
		allocate(f1(N_sp), f2(N_sp), g2(N_sp), g1(N_sp)) 
		
		sum_r(1:2) = (0.0, 0.0)
		sum_t(1:2) = (0.0, 0.0)
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located			
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
	   Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube		
		do nn_edge = 1, n_edge_number		
			nd = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n				
			sum_a = (0.0, 0.0)
			sum_b = (0.0, 0.0)
			sum_c = (0.0, 0.0)
			do tt = 1, 2 !
			   Rf_vec = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec
				Rf = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk
				Rm = struct%midpoints(struct%neighbours(m)%elements(tt))%point
				do ll = 1, 2
					Tf_vec = struct_kk%Integr_fn(nd)%pair_kk(ll)%Rfk_vec
					Tf = struct_kk%Integr_fn(nd)%pair_kk(ll)%Rfk
					Rn = struct%midpoints(struct%neighbours(nd)%elements(ll))%point							
					Tf_vec = conjg(Tf_vec)
					Tf = conjg(Tf)
					R1(1:N_sp, 1:3) = Rf_vec(1:N_sp, 1, 1:3)
					R2(1:N_sp, 1:3) = Rf_vec(1:N_sp, 2, 1:3)
					T1(1:N_sp, 1:3) = Tf_vec(1:N_sp, 1, 1:3)
					T2(1:N_sp, 1:3) = Tf_vec(1:N_sp, 2, 1:3)
					
					!do j = 1, N_sp
						tmp_1 = exp(im*k1*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab))) !(1.0, 0.0)!
						tmp_2 = exp(im*k2*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab)))!(1.0, 0.0)!
						ss1 = cross_rc_arr(N_sp, k_hat, R1)!*im*k1*tmp_1
						ss2 = cross_rc_arr(N_sp, k_hat, R2)!*im*k2*tmp_2
						do j = 1, N_sp
							ss1(j, :) = ss1(j, :)*tmp_1(j)*im*k1
							ss2(j, :) = ss2(j, :)*tmp_2(j)*im*k2
						end do
						sum_a = sum_a - sum((TL_1*dot_product_arr(N_sp, ss1, T1) + &
						TL_2*dot_product_arr(N_sp, ss2, T2))*TRF)
						
						f1 = omega*dot_product_arr(N_sp, R1, T1)
						g1 = Rf(:, 1)*Tf(:, 1)/omega
						f2 = omega*dot_product_arr(N_sp, R2, T2)
						g2 = Rf(:, 2)*Tf(:, 2)/omega	
						sum_b = sum_b + sum((TL_1*(my_1*f1 - g1/eps_1)*tmp_1 + TL_2*(my_2*f2 - g2/eps_2)*tmp_2)*TRF*im)
						sum_c = sum_c + sum((TL_1*(eps_1*f1 - g1/my_1)*tmp_1 + TL_2*(eps_2*f2 - g2/my_2)*tmp_2)*TRF*im)								
					!end do
				end do
			end do		
			
			sum_r(1) = sum_r(1) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs)
			sum_r(2) = sum_r(2) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
			sum_t(1) = sum_t(1) + sum_b*x1(nd) + sum_a*x1(nd + m_pairs)
			sum_t(2) = sum_t(2) - sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
		end do		
		if (transa == 'T') then
			rr = sum_t
		else if (transa == 'N') then
			rr = sum_r
		else 
			print*, 'Not a proper multiplication type!'
		end if		
	end subroutine FMM_cube_Tri_arr_PreCal
	
	subroutine FMM_edge_Tri_New_arr(m, n, struct, struct_cube, d_c, suma, sumb, sumc) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), dimension(N_sp) :: sum_a, sum_c, sum_b!, intent(out)
		complex(kind = dp), dimension(N_sp) :: tmp_1, tmp_2
		complex(kind = dp), dimension(N_sp, 2) :: Rf, Tf
		complex(kind = dp), dimension(N_sp) :: TL_1, TL_2
		complex(kind = dp), dimension(N_sp, 3) :: ss1, ss2, R1, R2, T1, T2
		complex(kind = dp), dimension(N_sp, 2, 3) :: Rf_vec, Tf_vec
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		complex(kind = dp), dimension(N_sp) :: f1, g1, f2, g2
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		

		!-------------for test use-----------------
		!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
		!       Rq
		!       o
		!       | 
		!	     |
		!	  Ra o-----------------------o Rb
		!                                |
		!                                |
		!                                o Rp
		!------------------------------------------
		
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
		intrinsic :: sqrt	
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		sum_a(1:N_sp) = (0.0, 0.0)
		sum_b(1:N_sp) = (0.0, 0.0)
		sum_c(1:N_sp) = (0.0, 0.0)
		do tt = 1, 2 !Attention, sometimes, t and l might not be smaller than 2
			
			Rf_vec = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec
			Rf = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk
			Rm = struct%midpoints(struct%neighbours(m)%elements(tt))%point
			do ll = 1, 2			
				Tf_vec = struct_kk%Integr_fn(n)%pair_kk(ll)%Rfk_vec
				Tf = struct_kk%Integr_fn(n)%pair_kk(ll)%Rfk
				Rn = struct%midpoints(struct%neighbours(n)%elements(ll))%point							
				Tf_vec = conjg(Tf_vec)
				Tf = conjg(Tf)	
				R1(1:N_sp, 1:3) = Rf_vec(1:N_sp, 1, 1:3)
				R2(1:N_sp, 1:3) = Rf_vec(1:N_sp, 2, 1:3)
				T1(1:N_sp, 1:3) = Tf_vec(1:N_sp, 1, 1:3)
				T2(1:N_sp, 1:3) = Tf_vec(1:N_sp, 2, 1:3)
					
				tmp_1 = exp(im*k1*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab))) !(1.0, 0.0)!
				tmp_2 = exp(im*k2*dot_rr_arr(N_sp, k_hat, (Rn-Rm+Rab)))!(1.0, 0.0)!
				ss1 = cross_rc_arr(N_sp, k_hat, R1)*im*k1!
				ss2 = cross_rc_arr(N_sp, k_hat, R2)*im*k2!						
				do j = 1, N_sp
					ss1(j, :) = ss1(j, :)*tmp_1(j)
					ss2(j, :) = ss2(j, :)*tmp_2(j)
				end do
				sum_a = sum_a - (TL_1*dot_product_arr(N_sp, ss1, T1) + &
				TL_2*dot_product_arr(N_sp, ss2, T2))*TRF
				
				f1 = omega*dot_product_arr(N_sp, R1, T1)
				g1 = Rf(:, 1)*Tf(:, 1)/omega
				f2 = omega*dot_product_arr(N_sp, R2, T2)
				g2 = Rf(:, 2)*Tf(:, 2)/omega	
				sum_b = sum_b + (TL_1*(my_1*f1 - g1/eps_1)*tmp_1 + TL_2*(my_2*f2 - g2/eps_2)*tmp_2)*TRF*im
				sum_c = sum_c + (TL_1*(eps_1*f1 - g1/my_1)*tmp_1 + TL_2*(eps_2*f2 - g2/my_2)*tmp_2)*TRF*im	
			end do
		end do
		suma = sum(sum_a)
		sumb = sum(sum_b)
		sumc = sum(sum_c)
	end subroutine FMM_edge_Tri_New_arr	
		 
	subroutine Receivef_PreC(struct, m, tt, Recf_vec, Recf)
	! rc: receive
	! rd: radiation
	! tr: transfer function
	! R_oa :radius around the observation point
	! R_pb: radius around the source point
	! R_ab : distance between the two points
		implicit none
		integer, intent(in) :: m, tt !sampling point on a sphere				
		real(kind = dp) :: R_m(3) 
		
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), R_om(3), R_o(3), fm(3)
		complex(kind = dp), dimension(N_sp, 2, 3), intent(out) :: Recf_vec
		complex(kind = dp), dimension(N_sp, 2), intent(out) :: Recf
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
		real(kind = dp) :: area, Len_m, dfm
		integer :: i, j, n, pos_neg
		type(Structure) :: struct
		
		call Quadrature_tri(ng, a, b, w)
		call fn_parameter_midpoint(struct, m, tt, pm1, pm2, vm_t, R_m, Len_m, area)				
		pos_neg = -1*(((tt-1)*tt)-1)
		dfm = pos_neg*Len_m!
		do i = 1, N_sp				
			Recf_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Recf_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Recf(i, 1:2) = (0.0, 0.0)		
		   do j = 1, ng
				R_o = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
				fm = 0.5*(R_o-vm_t)
				R_om = R_o - R_m			
				tmp = dot_product(k_hat(i, :), R_om(:))
				ff1 = exp(-im*k1*tmp)*dfm*w(j)
				ff2 = exp(-im*k2*tmp)*dfm*w(j)
				Recf_vec(i, 1, :) = Recf_vec(i, 1, :) + fm*ff1
				Recf_vec(i, 2, :) = Recf_vec(i, 2, :) + fm*ff2
				Recf(i, 1) = Recf(i, 1) + ff1
				Recf(i, 2) = Recf(i, 2) + ff2
			end do			
		end do !
		return
	end subroutine Receivef_PreC
	
	subroutine FMM_edge_Tri_PreCal(m, n, struct, struct_cube, d_c, sum_a, sum_b, sum_c) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		!type(Structure_k) :: struct_kk
		
		complex(kind = dp), intent(out) :: sum_a, sum_c, sum_b
		complex(kind = dp) :: tmp_1, tmp_2
		complex(kind = dp), dimension(N_sp, 2) :: Tf, Rf
		complex(kind = dp), dimension(N_sp) :: TL_1, TL_2
		complex(kind = dp), dimension(N_sp, 2, 3) :: Tf_vec, Rf_vec
 
		complex(kind = dp) :: ff1, gg1, ff2, gg2, ss1(3), ss2(3)
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
 
		!-------------for test use-----------------
		!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
		!       Rq
		!       o
		!       | 
		!	     |
		!	  Ra o-----------------------o Rb
		!                                |
		!                                |
		!                                o Rp
		!------------------------------------------
		
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
		intrinsic :: sqrt
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		do tt = 1, 2 !Attention: sometimes, t and l might not be smaller than 2
			!call Receivef_PreC(struct, m, tt, Rfk_vec, Rfk)
			Rm = struct%midpoints(struct%neighbours(m)%elements(tt))%point
			Rf_vec = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec
			Rf = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk
			do ll = 1, 2
				!call Receivef_PreC(struct, n, ll, Tfk_vec, Tfk)
				Rn = struct%midpoints(struct%neighbours(n)%elements(ll))%point	
				Tf_vec = conjg(struct_kk%Integr_fn(n)%pair_kk(ll)%Rfk_vec)
				Tf = conjg(struct_kk%Integr_fn(n)%pair_kk(ll)%Rfk)
				
				do j = 1, N_sp
					tmp_1 = exp(im*k1*dot_product(k_hat(j, :), (Rn-Rm+Rab))) !(1.0, 0.0)!
					tmp_2 = exp(im*k2*dot_product(k_hat(j, :), (Rn-Rm+Rab)))!(1.0, 0.0)!
						
					ss1 = im*k1*(cross_rc(k_hat(j, :), Rf_vec(j, 1, :)))*tmp_1
					ss2 = im*k2*(cross_rc(k_hat(j, :), Rf_vec(j, 2, :)))*tmp_2
						
					sum_a = sum_a - (TL_1(j)*dot_product(conjg(ss1(:)), Tf_vec(j, 1, :)) + &
					TL_2(j)*dot_product(conjg(ss2(:)), Tf_vec(j, 2, :)))*TRF(j)
						
					ff1 = omega*dot_product(conjg(Rf_vec(j, 1, :)), Tf_vec(j, 1, :))
					gg1 = Rf(j, 1)*Tf(j, 1)/omega
					ff2 = omega*dot_product(conjg(Rf_vec(j, 2, :)), Tf_vec(j, 2, :))
					gg2 = Rf(j, 2)*Tf(j, 2)/omega
						
					sum_b = sum_b + (TL_1(j)*(my_1*ff1 - gg1/eps_1)*tmp_1 + TL_2(j)*(my_2*ff2 - gg2/eps_2)*tmp_2)*TRF(j)*im
					sum_c = sum_c + (TL_1(j)*(eps_1*ff1 - gg1/my_1)*tmp_1 + TL_2(j)*(eps_2*ff2 - gg2/my_2)*tmp_2)*TRF(j)*im								
				end do	
			end do
		end do
	end subroutine FMM_edge_Tri_PreCal
	
	subroutine Test_RT_Spherical_Coefficient(struct, struct_cube)
		implicit none
		real(kind = dp) :: Ra(3)		
		integer :: Kn  !sampling point on a sphere				
		integer :: i, j, k, ll, fnu, n_phi, n, m, Np, m_cube, nd
		real(kind = dp), dimension(:, :), allocatable :: xw
		complex(kind = dp) :: r_tmp(3) 
		real(kind = dp) :: dphi	
		type(Structure), intent(in) :: struct	
		type(Structure_cube), intent(in) :: struct_cube	
		
		type(Integration_fn) :: Int_fm
		complex(kind = dp), dimension(N_sp, 2, 3) :: Rf_vec, Tf_vec
		complex(kind = dp), dimension(N_sp, 2) :: Rf
		
		!--------------------------------------------------------------
		type(list_list_cmplx) :: alm_c, alm_d! alm_a, alm_b,
		type(list_list_cmplx) :: Ylm
		intrinsic conjg
		type(list_list_cmplx_vector) :: alm_a, alm_b		
		nd = 2
		Kn = 12	
		n_phi = 2*Kn + 1
		Np = Kn*n_phi
		
		allocate(xw(Np, 3))		
		xw = theta_phi_w(Kn)
		dphi = 2*PI/n_phi
		m_cube = struct_cube%pair_cubes(nd)%cube_index	!
		Ra = struct_cube%cubes(m_cube)%cube_position	
		call RT_Spherical_Coefficient(struct, nd, Ra, alm_a, alm_b, alm_c, alm_d)
		r_tmp = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		do ll = 1, 2! 
			do i = 1, Kn
				do j = 1, n_phi
					k = (i-1)*n_phi + j
					call spherical_harmonics_Darve(Kn+1, Ylm, xw(k, 1), xw(k, 2))
					do n = 0, Kn
						do m = -n, n
							r_tmp = r_tmp + Ylm%item(n)%item(m)*xw(k, 3)*dphi!alm_a%vector(1:3)%item(n)%item(m)*
							!print*, '          '
							!print*, 'r_tmp =', r_tmp
							!alm_b%vector(1:3)%item(n)%item(m) = alm_b%vector(1:3)%item(n)%item(m) + r_tmp(1:3)*conjg(Ylm%item(n)%item(m))*TRF(k)
							!alm_c%item(n)%item(m) = alm_c%item(n)%item(m) + Rf(k, 1)*conjg(Ylm%item(n)%item(m))*TRF(k)
							!alm_d%item(n)%item(m) = alm_d%item(n)%item(m) + Rf(k, 2)*conjg(Ylm%item(n)%item(m))*TRF(k)
						end do
					end do
				end do
			end do
		end do
		print*, 'The up-sampled summation =', r_tmp
	end subroutine Test_RT_Spherical_Coefficient
	
	subroutine RT_Spherical_Coefficient(struct, nd, Ra, alm_a, alm_b, alm_c, alm_d)
	! R_ab : distance between the two points
		implicit none		
		real(kind = dp), intent(in) :: Ra(3)		
		integer, intent(in) :: nd  !sampling point on a sphere				
		integer :: i, j, k, ll, fnu, n_phi, n, m
		real(kind = dp), dimension(:, :), allocatable :: xw
		complex(kind = dp), dimension(N_sp, 1:3) :: R1, R2
		complex(kind = dp) :: r_tmp(3) 
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	

		type(Integration_fn) :: Int_fm
		complex(kind = dp), dimension(N_sp, 2, 3) :: Rf_vec, Tf_vec
		complex(kind = dp), dimension(N_sp, 2) :: Rf
		!--------------------------------------------------------------
		type(list_list_cmplx), intent(out) :: alm_c, alm_d! alm_a, alm_b,
		type(list_list_cmplx) :: Ylm
		intrinsic conjg
		type(list_list_cmplx_vector) :: alm_a, alm_b
		integer :: K_m
        K_m = 12
		fnu = 0		
		n_phi = 2*K_m+1
		call init_list(Ylm, fnu, K_m+1)
		call init_list(alm_a%vector(1), fnu, K_m+1)		
		call init_list(alm_a%vector(2), fnu, K_m+1)
		call init_list(alm_a%vector(3), fnu, K_m+1)
		call init_list(alm_b%vector(1), fnu, K_m+1)		
		call init_list(alm_b%vector(2), fnu, K_m+1)
		call init_list(alm_b%vector(3), fnu, K_m+1)				
		call init_list(alm_c, fnu, K_m+1)		
		call init_list(alm_d, fnu, K_m+1)
		
		call RT_function(struct, nd, Ra, Int_fm)		
		allocate(xw(N_sp, 3))
		xw = theta_phi_w(K_m)		
		do ll = 1, 2! 
			Rf_vec = Int_fm%pair_kk(ll)%Rfk_vec
			Rf = Int_fm%pair_kk(ll)%Rfk
			R1(1:N_sp, 1:3) = Rf_vec(1:N_sp, 1, 1:3)
			R2(1:N_sp, 1:3) = Rf_vec(1:N_sp, 2, 1:3)
			do i = 1, K_m
				do j = 1, n_phi
					k = (i-1)*n_phi + j
					call spherical_harmonics_Darve(K_m+1 , Ylm, xw(k, 1), xw(k, 2))
					do n = 0, K_m
						do m = -n, n
							r_tmp = cross_rc(k_hat(k, :), R1(k, :))
							!r_tmp = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
							alm_a%vector(1:3)%item(n)%item(m) = alm_a%vector(1:3)%item(n)%item(m) + r_tmp(1:3)*conjg(Ylm%item(n)%item(m))*TRF(k)
							r_tmp = cross_rc(k_hat(k, :), R2(k, :))
							alm_b%vector(1:3)%item(n)%item(m) = alm_b%vector(1:3)%item(n)%item(m) + r_tmp(1:3)*conjg(Ylm%item(n)%item(m))*TRF(k)
							alm_c%item(n)%item(m) = alm_c%item(n)%item(m) + Rf(k, 1)*conjg(Ylm%item(n)%item(m))*TRF(k)
							alm_d%item(n)%item(m) = alm_d%item(n)%item(m) + Rf(k, 2)*conjg(Ylm%item(n)%item(m))*TRF(k)
						end do
					end do
				end do
			end do
		end do		
		r_tmp = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)			
	!-----------For the test------------------
		do ll = 1, 2! 		
			do i = 1, K_m
				do j = 1, n_phi
					k = (i-1)*n_phi + j
					r_tmp = r_tmp + TRF(k)!cross_rc(k_hat(k, :), R1(k, :))
				end do				
			end do
		end do	
		print*, 'The original summation =', r_tmp
	end subroutine RT_Spherical_Coefficient
	
	!function TL_km_arr(r_ab, Lm, N_sp, K_m)
 !
 !     integer, intent(in) :: Lm, N_sp     !Lm is the order limit of the polynomials, N is the sampling points on the sphere
 !     real(kind = dp), intent(in) :: r_ab(3)
 !     real(kind = dp) :: r_hat(3), x, pm_arr(N_sp, Lm)
 !     complex(kind = dp) :: r_h
 !     complex(kind = dp), intent(in) :: K_m
 !     real (kind = dp), dimension(Lm) :: pm, dummy 
 !     complex(kind = dp) :: TL_km_arr(N_sp)
 !     complex(kind = dp), dimension(Lm) :: hl 
 !     integer :: j, fnu
	!	
	!	
 !     fnu = 0        
 !     r_hat = r_ab/vec_len(r_ab)
 !     r_h =  vec_len(r_ab)*K_m  
 !     hl =  lib_math_hankel_spherical_2(r_h, fnu, Lm) 
 !     TL_km_arr(1:N_sp) = (0.0, 0.0)         
 !     do j = 1, N_sp
 !        x = dot_product(r_hat, k_hat(j, :))
 !        call lib_math_legendre_polynomial(x, fnu, Lm, pm, dummy)
 !        pm_arr(j, :) = pm(:)
 !     end do
	!	
 !     do j = 1, Lm
 !        TL_km_arr(:) = TL_km_arr(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*K_m/(4*PI)**2 ! 
	!	end do
 !     return        
	!end function TL_km_arr
	
	!subroutine FMM_Far_MVProduct_Tri_OMP(transa, m_pairs, struct, struct_cube, d_c, x1, rr)
	!!For iterative solution	
	!	integer, intent(in) :: m_pairs
	!	character(len = 1), intent(in) :: transa
	!	integer :: n_cube, m_cube, n, n_edge_number, nn_edge, nd, n_cube_number, m
	!	type(Structure), intent(in) :: struct
	!	type(Structure_cube), intent(in) :: struct_cube		
	!	complex(kind = dp):: sum_r(2)
	!	complex(kind = dp), dimension(:), allocatable, intent(out) :: rr
	!	complex(kind = dp), dimension(:), allocatable, intent(in) :: x1	
	!	real(kind = dp), intent(in) :: d_c !
	!	allocate(rr(2*m_pairs))
	!	
	!	call timestamp()
	!	!$omp parallel private (m, sum_r, m_cube, n_cube, n, n_cube_number) &
	!	!$omp shared (m_pairs, struct, struct_cube, x1, d_c, rr, transa)
	!	!$omp do
	!	!j = omp_get_thread_num ( )
	!	!write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', j
	!	do m = 1, m_pairs
	!		rr(m) = (0.0, 0.0)			
	!		rr(m + m_pairs) = (0.0, 0.0)
	!		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
	!		n_cube_number = size(struct_cube%cubes(m_cube)%cubes_far)	
	!		do n = 1, n_cube_number 
	!			n_cube = struct_cube%cubes(m_cube)%cubes_far(n) !cube index in the far region of the cube m_cube
	!			call FMM_cube_Tri_KL_spherical(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, sum_r)
	!			rr(m) = rr(m) +  sum_r(1)
	!			rr(m + m_pairs) = rr(m + m_pairs) +  sum_r(2)
	!		end do
	!	end do
	!	!$omp end do
	!	!$omp end parallel			
	!	call timestamp()
	!	return		
	!end subroutine FMM_Far_MVProduct_Tri_OMP	!	
	
	subroutine Receivef_PreCal(struct)	
	! R_ab : distance between the two points
		implicit none
		integer :: m, tt, m_pairs !sampling point on a sphere				
		integer :: i, j, n, pos_neg
        integer :: K_m
		
		real(kind = dp) :: R_m(3) 		
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: pm1(3), pm2(3), vm_t(3), R_om(3), R_o(3), fm(3)		
		complex(kind = dp) :: ff1, ff2, tmp, tmp_r(3)	
		real(kind = dp) :: area, Len_m, dfm		
		type(Structure), intent(in) :: struct	
		m_pairs = size(struct%neighbours)		
		
		allocate(struct_kk%Integr_fn(m_pairs)) 		!		
		call Quadrature_tri(ng, a, b, w)
        K_m = 12
		print*, 'K_m=', K_m
		do m = 1, m_pairs
		   do tt = 1, 2
				call fn_parameter_midpoint(struct, m, tt, pm1, pm2, vm_t, R_m, Len_m, area)				!Needs to be re-written 
				pos_neg = -1*(((tt-1)*tt)-1)!put it outside this routine
				dfm = pos_neg*Len_m 			
				allocate(struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(N_sp, 2, 3))
				allocate(struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(N_sp, 2))
				do i = 1, N_sp
					struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
					struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
					struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(i, 1:2) = (0.0, 0.0)	
					
					do j = 1, ng
						R_o = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
						fm = 0.5*(R_o-vm_t)
						R_om = R_o - R_m			
						tmp = dot_product(k_hat(i, :), R_om(:))
						ff1 = exp(-im*k1*tmp)*dfm*w(j)
						ff2 = exp(-im*k2*tmp)*dfm*w(j)
						struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 1, :) = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 1, :) + fm*ff1
						struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 2, :) = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk_vec(i, 2, :) + fm*ff2
						struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(i, 1) = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(i, 1) + ff1
						struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(i, 2) = struct_kk%Integr_fn(m)%pair_kk(tt)%Rfk(i, 2) + ff2
					end do						
				end do		
			end do !		
		end do		
	end subroutine Receivef_PreCal
	
	subroutine FMM_edge_Tri_RT_Upsampling(m, n, struct, struct_cube, d_c, sum_a, sum_b, sum_c) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		integer :: Nq, Np, Kn
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), dimension(N_sp) :: sum_a, sum_c, sum_b 
		complex(kind = dp), dimension(:), allocatable ::suma_p, sumb_p, sumc_p
		complex(kind = dp), dimension(N_sp, 2) :: Rf, Tf
		complex(kind = dp), dimension(N_sp) :: TL_1, TL_2
		complex(kind = dp), dimension(N_sp, 3) :: ss1, ss2, R1, R2, T1, T2
		complex(kind = dp), dimension(N_sp, 2, 3) :: Rf_vec, Tf_vec
		!complex(kind = dp), intent(out) :: suma, sumb, sumc
		complex(kind = dp), dimension(N_sp) :: f1, g1, f2, g2
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
		real(kind  = dp), dimension(:), allocatable :: wk_theta, wk_phi
		real(kind  = dp), dimension(:, :), allocatable :: TRF_new
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
		type(Integration_fn) :: Int_fn, Int_fm
		!-------------for test use-----------------
		!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
		!       Rq
		!       o
		!       | 
		!	     |
		!	  Ra o-----------------------o Rb
		!                                |
		!                                |
		!                                o Rp
		!------------------------------------------
		
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
		intrinsic :: sqrt	
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		
		sum_a(1:N_sp) = (0.0, 0.0)
		sum_b(1:N_sp) = (0.0, 0.0)
		sum_c(1:N_sp) = (0.0, 0.0)
		
		call RT_function(struct, m, Ra, Int_fm)
		call RT_function(struct, n, Rb, Int_fn)
		do tt = 1, 2 !Attention, sometimes, t and l might not be smaller than 2			
			Rf_vec = Int_fm%pair_kk(tt)%Rfk_vec
			Rf = Int_fm%pair_kk(tt)%Rfk
			R1(1:N_sp, 1:3) = Rf_vec(1:N_sp, 1, 1:3)
			R2(1:N_sp, 1:3) = Rf_vec(1:N_sp, 2, 1:3)	

			do ll = 1, 2			
				Tf_vec = Int_fn%pair_kk(ll)%Rfk_vec
				Tf = Int_fn%pair_kk(ll)%Rfk				
				Tf_vec = conjg(Tf_vec)
				Tf = conjg(Tf)	
				T1(1:N_sp, 1:3) = Tf_vec(1:N_sp, 1, 1:3)
				T2(1:N_sp, 1:3) = Tf_vec(1:N_sp, 2, 1:3)					
				
				ss1 = cross_rc_arr(N_sp, k_hat, R1)*im*k1!
				ss2 = cross_rc_arr(N_sp, k_hat, R2)*im*k2!		
			
				sum_a = sum_a - (TL_1*dot_product_arr(N_sp, ss1, T1) + &
				TL_2*dot_product_arr(N_sp, ss2, T2))
				
				f1 = omega*dot_product_arr(N_sp, R1, T1)
				g1 = Rf(:, 1)*Tf(:, 1)/omega
				f2 = omega*dot_product_arr(N_sp, R2, T2)
				g2 = Rf(:, 2)*Tf(:, 2)/omega	
				sum_b = sum_b + (TL_1*(my_1*f1 - g1/eps_1) + TL_2*(my_2*f2 - g2/eps_2))*im
				sum_c = sum_c + (TL_1*(eps_1*f1 - g1/my_1) + TL_2*(eps_2*f2 - g2/my_2))*im	
			end do
		end do
		
		!Kn = K_m
		!TRF_new = theta_phi_w(Kn)
		!Np = (2*Kn+1)*Kn	
		!Nq = K_m*(2*K_m+1)
		!allocate(suma_p(Np))
		!suma_p(1:Np)=(0.0, 0.0)
		!print*,  ''
		!print*, 'size(Akk(:, 1))=', size(Akk, 1)
		!print*, 'size(Akk(:, 2))=', size(Akk, 2)		
		!do tt = 1, size(Akk, 1)
		!	suma_p(tt) = suma_p(tt) + dot_product(conjg(Akk(tt, :)), sum_a(:))
		!end do
		!
		!print*, 'sum of sum_a =', sum(sum_a*TRF)
		!print*, 'sum of sum_p =', sum(suma_p(:)*TRF_New(:, 3))

		!open (unit=206, file = 'suma_Km7.txt', action="write",status = 'replace')
		!do tt = 1, K_m*(2*K_m + 1)
		!	write (206, '(2(es19.12, tr5))') real(sum_a(tt)*TRF(tt)), imag(sum_a(tt)*TRF(tt))
		!end do
		!sum_a = (1.0, 0.0)
		
		!++++++++++++++++++++++++++		!
		!open (unit=206, file = 'sump_Km7.txt', action="write",status = 'replace')
		!do tt = 1, Kn*(2*Kn + 1)
		!	write (206, '(2(es19.12, tr5))') real(suma_p(tt)*TRF_new(tt, 3)), imag(suma_p(tt)*TRF_new(tt, 3))
		!end do
		!!+++++++++++++++++++++++++++
		!!
		!print*, 'sum of sum_a=', sum(sum_a*TRF)
		!print*, 'sum of sum_p=', sum(suma_p(:)*TRF_new(:, 3))		
	end subroutine FMM_edge_Tri_RT_Upsampling
	
	subroutine FMM_edge_Tri_RT(m, n, struct, struct_cube, d_c, suma, sumb, sumc) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), dimension(N_sp) :: sum_a, sum_c, sum_b!, intent(out)		
		complex(kind = dp), dimension(N_sp, 2) :: Rf, Tf
		complex(kind = dp), dimension(N_sp) :: TL_1, TL_2
		complex(kind = dp), dimension(N_sp, 3) :: ss1, ss2, R1, R2, T1, T2
		complex(kind = dp), dimension(N_sp, 2, 3) :: Rf_vec, Tf_vec
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		complex(kind = dp), dimension(N_sp) :: f1, g1, f2, g2
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3), Rm(3), Rn(3)
		real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
		type(Integration_fn) :: Int_fn, Int_fm
		!-------------for test use-----------------
		!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
		!       Rq
		!       o
		!       | 
		!	     |
		!	  Ra o-----------------------o Rb
		!                                |
		!                                |
		!                                o Rp
		!------------------------------------------
		
		integer :: j, pos_neg, Lm, m_cube, n_cube, tt, ll
		intrinsic :: sqrt	
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(Rab, Lm, k1)
		!print*, 'TL_1=', TL_1(1:2)
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(Rab, Lm, k2)
		!print*, 'TL_2=', TL_2(1:2)
		
		sum_a(1:N_sp) = (0.0, 0.0)
		sum_b(1:N_sp) = (0.0, 0.0)
		sum_c(1:N_sp) = (0.0, 0.0)
		
		call RT_function(struct, m, Ra, Int_fm)
		call RT_function(struct, n, Rb, Int_fn)
		do tt = 1, 2 !Attention, sometimes, t and l might not be smaller than 2			
			Rf_vec = Int_fm%pair_kk(tt)%Rfk_vec
			Rf = Int_fm%pair_kk(tt)%Rfk
			R1(1:N_sp, 1:3) = Rf_vec(1:N_sp, 1, 1:3)
			R2(1:N_sp, 1:3) = Rf_vec(1:N_sp, 2, 1:3)	

			do ll = 1, 2			
				Tf_vec = Int_fn%pair_kk(ll)%Rfk_vec
				Tf = Int_fn%pair_kk(ll)%Rfk				
				Tf_vec = conjg(Tf_vec)
				Tf = conjg(Tf)	
				T1(1:N_sp, 1:3) = Tf_vec(1:N_sp, 1, 1:3)
				T2(1:N_sp, 1:3) = Tf_vec(1:N_sp, 2, 1:3)					
				
				ss1 = cross_rc_arr(N_sp, k_hat, R1)*im*k1!
				ss2 = cross_rc_arr(N_sp, k_hat, R2)*im*k2!		

				sum_a = sum_a - (TL_1*dot_product_arr(N_sp, ss1, T1) + &
										TL_2*dot_product_arr(N_sp, ss2, T2))*TRF !				
				f1 = omega*dot_product_arr(N_sp, R1, T1)
				g1 = Rf(:, 1)*Tf(:, 1)/omega
				f2 = omega*dot_product_arr(N_sp, R2, T2)
				g2 = Rf(:, 2)*Tf(:, 2)/omega	
				sum_b = sum_b + (TL_1*(my_1*f1 - g1/eps_1) + TL_2*(my_2*f2 - g2/eps_2))*TRF*im
				sum_c = sum_c + (TL_1*(eps_1*f1 - g1/my_1) + TL_2*(eps_2*f2 - g2/my_2))*TRF*im	
			end do
		end do
		suma = sum(sum_a)
		sumb = sum(sum_b)
		sumc = sum(sum_c)
	end subroutine FMM_edge_Tri_RT
	
	end module MLFMM_Tri