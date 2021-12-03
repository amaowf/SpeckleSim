module lib_sie_math
	
	use lib_sie_constants
	use lib_sie_type_function		
	use Gauss_Quadrature_Fast
	
	use libmath
	
	implicit none
	
	real(dp), dimension(:), allocatable, public :: f_boundary_rs
	
	interface init_list_sie
		module procedure init_list_sie_c
		module procedure init_list_sie_v
	end interface
	
	contains
	
	subroutine init_list_sie_c(C, m, n)
      use libmlfmm
      type(lib_ml_fmm_coefficient), intent(inout) :: C        
      integer, intent(in) :: m, n
      integer :: i
      if (m==0) then
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
			!print*, 'size of C in init', size(C%a_nm%item)
         do i = 1, m
               allocate(C%a_nm%item(i)%item(1:n))
               C%a_nm%item(i)%item(1:n) = (0.0, 0.0)
         end do 
 
			if (allocated(C%b_nm%item)) then
               deallocate(C%b_nm%item)
         end if
         allocate(C%b_nm%item(1:m))
         do i = 1, m
               allocate(C%b_nm%item(i)%item(1:n))
               C%b_nm%item(i)%item(1:n) = (0.0, 0.0)
         end do            	
      end if
	end subroutine
	 
	subroutine init_list_sie_v(C, m, n)
      use libmlfmm
      type(lib_ml_fmm_v), intent(inout) :: C        
      integer, intent(in) :: m, n
      integer :: i
      if (m==0) then
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
 
			if (allocated(C%b_nm%item)) then
               deallocate(C%b_nm%item)
         end if
         allocate(C%b_nm%item(1:m))
         do i = 1, m
               allocate(C%b_nm%item(i)%item(1:n))
               C%b_nm%item(i)%item(1:n) = (0.0, 0.0)
         end do 
		end if
	end subroutine
	!
	subroutine Gaussian_beam(w0, r_local, illumination, Es, kd)!, theta, pol, lambda, r0
		real(dp), intent(in) :: r_local(3)
		type(lib_sie_illumination_parameter), intent(in) :: illumination
		type(vector_c), intent(out) :: Es
		real(dp), intent(in) :: w0
		
		!kd: wavevector
		complex(dp), intent(in) :: kd
		
		
		!dummy
		integer :: i
		real(dp) :: zR, ra(3), Rz, wz, phi, A(3, 3), op, expf
				
      zR = PI*w0**2/illumination%lambda
      op = 0.0
		! A: rotation matrix
      A(1, :) = (/cos(illumination%theta_in), op, -sin(illumination%theta_in)/)
      A(2, :) = (/op, cos(illumination%theta_in)/cos(illumination%theta_in), op/)
      A(3, :) = (/sin(illumination%theta_in), op, cos(illumination%theta_in)/)

      do i = 1, 3
         ra(i) = dot_product(A(i, :), r_local(:))
      end do

      if (ra(3) == 0.0) then
         Rz = 0.0
         phi = 0.0
         wz = w0
      else
         Rz = ra(3)*(1 + (zR/ra(3))**2)
         phi = atan(ra(3)/zR)
         wz = w0*sqrt(1+(ra(3)/zR)**2)
      end if

		expf = -(ra(1)**2+ra(2)**2)/wz**2
		if (expf < -100) then
         Es%vector(1:3) = 0.0
		else if (ra(3) == 0.0) then
         Es%vector = illumination%E_in*exp(-(ra(1)**2+ra(2)**2)/wz**2)*exp(im*phi)!
		else
			Es%vector = illumination%E_in*w0/wz*exp(-(ra(1)**2+ra(2)**2)/wz**2)* &
				exp(-im*(kd*ra(3) + kd*(ra(1)**2+ra(2)**2)/(2*Rz)-phi))!
		end if		
		return
	end subroutine Gaussian_beam
	
	!According to T. Pahl 2020 (Optics Express)
	!Using plane wave illumination with a normal incidence 
	subroutine conical_illumination(p_obj, r_local, illumination, E_surf)	
		type(lib_sie_illumination_parameter), intent(in) :: illumination
		real(dp), intent(in) :: r_local(3)
		type(lib_sie_parameter_objective), intent(in) :: p_obj		
		type(vector_c), intent(out) :: E_surf
		
		!dummy
		integer :: nm, m, i, j
		real(dp), dimension(:), allocatable :: jn
		real(dp) :: n_air, k_in(3), kx, ky, kz, ph, E_in(2) !ph: the pupil function, can be extended		
		real(dp) :: theta_min, theta_max, phi_min, phi_max, k_rho, x_j, trf(3),  E_tmp(3)
		
		real(dp) :: k0, dtheta, dphi, theta, phi,filter_ph, factor_m, factor_lp
		integer :: Nx, Ny
		
		complex(dp) :: E_surf_tmp(3), phase_tmp
		nm = 2
		m = 2
		
		allocate(jn(1:m))
		
		n_air = 1.0
		Nx = 21
		Ny = Nx
		
		!Arbitrary factor to obtain the value of maximum E_field around 1.0
		factor_m = 8.0	
		
		theta_min = 0.0
		theta_max = asin(p_obj%NA/n_air);
		
		phi_min = 0
		phi_max = -2*PI
		
		if (illumination%k_in(3) .le. 0.0) then
			factor_lp = -1
		else
			factor_lp = 1
		end if			
		
		dtheta = (theta_max-theta_min)/(Nx-1)
		dphi = (phi_max-phi_min)/(Ny-1)
		
		k0 = PI*2/illumination%lambda		
		E_in = illumination%E_in(1:2)
			
		Ph = 1.0 !for temparory use	
		E_surf%vector(1:3) = (0.0, 0.0)				
		do i = 1, Nx !
			theta = theta_min + dtheta*(i) 
			do j = 1, Ny-1
			   phi = phi_min + dphi*(j-1)
				kx = sin(theta)*cos(phi)				
				ky = sin(theta)*sin(phi)
				kz = cos(theta)*factor_lp
				k_in = (/kx, ky, kz/)*k0								
				k_rho = sqrt(kx**2 + ky**2)*k0
				x_j = k_rho*p_obj%r_ph/p_obj%M				
				jn = lib_math_bessel_spherical_first_kind(x_j, 0, nm)				
							
				call field_rotation_conical_illumination(theta, phi, E_in, E_tmp)
			
				trf = dtheta*dphi*sin(theta)*cos(theta)*E_tmp
				
				!check which order should be in jn()
				filter_ph = 2*PI*p_obj%r_ph/(k_rho*p_obj%M)*jn(1)*Ph
				
				!First using no filter function
				E_surf_tmp = exp(-im*dot_product(k_in, r_local))*k0**2/(4*PI**2)*trf!*filter_ph
				E_surf%vector = E_surf%vector + E_surf_tmp*factor_m				
			end do		
		end do
		deallocate(jn)		
	end subroutine
	
	!According to equation (2) and (5) in Pahl, 2021
	subroutine field_rotation_conical_illumination(theta, phi, Ein_t, E_in)
		real(dp), intent(in) :: theta, phi
		real(dp), intent(in) :: Ein_t(2)
		real(dp) :: Rm(3, 2)
		real(dp), intent(out) :: E_in(3)
		real(dp) :: a, b, c, d
		
		a = sin(theta)		
		b = cos(theta)
		c = sin(phi)
		d = cos(phi)
		
		Rm(1,1) = d**2*b + c**2
		Rm(1,2) = c*d*(b-1)
		Rm(2,1) = c*d*(b-1)
		Rm(2,2) = c**2*b + d**2
		Rm(3,1) = a*d
		Rm(3,2) = a*c
		!
		E_in(1) = Rm(1, 1)*Ein_t(1) + Rm(1, 2)*Ein_t(2)
		E_in(2) = Rm(2, 1)*Ein_t(1) + Rm(2, 2)*Ein_t(2)
		E_in(3) = Rm(3, 1)*Ein_t(1) + Rm(3, 2)*Ein_t(2)
	
	end subroutine
	
	subroutine get_evaluation_points(evaluation_p, r_local, evaluation_type)! 
		type(point), dimension(:), allocatable, intent(out) :: r_local
		character(len = 10), intent(in) :: evaluation_type
		type(lib_sie_evaluation_parameter_type), intent(in) :: evaluation_p

		!dummy
		real(dp ) :: dx_a, dx_b
		real(dp) :: phi, theta !ref_pl, 
		integer :: i, j, n_a, n_b
		
		n_a = evaluation_p%N_dim(1)
		n_b = evaluation_p%N_dim(2)

		allocate(r_local(n_a*n_b))
		
		dx_a = (evaluation_p%dim_a(2) - evaluation_p%dim_a(1))/(n_a-1)
		
		if ((evaluation_type .eq. 'rcs_p') .or. (evaluation_type .eq. 'rcs_n') .or. (evaluation_type .eq. 'BRDF'))then
			dx_b = 0.0
		else
			dx_b = (evaluation_p%dim_b(2) - evaluation_p%dim_b(1))/(n_b-1)
		end if
		
		if (evaluation_type .eq. 'xy')then ! xy-plot		
			do i = 1, n_b
				do j = 1, n_a
					r_local((i-1)*n_a+j)%point= &
					(/evaluation_p%dim_a(2) - dx_a*(j-1), evaluation_p%dim_b(2) - dx_b*(i-1), evaluation_p%dim_c/)
				end do
			end do
		else if (evaluation_type .eq. 'yz')then ! yz-plot      !     
			do i=1, n_b
				do j=1, n_a
					r_local((i-1)*n_a + j)%point=(/evaluation_p%dim_c, evaluation_p%dim_a(2) - dx_a*(j-1), &
					evaluation_p%dim_b(2) - dx_b*(i-1) /)
				end do
			end do
        
		else if (evaluation_type .eq. 'xz')then ! xz-plot
			do i=1, n_b
				do j=1, n_a
					r_local((i-1)*n_a+j)%point=(/evaluation_p%dim_a(1) + dx_a*(j-1), evaluation_p%dim_c, &
					evaluation_p%dim_b(1) + dx_b*(i-1)/)
				end do
			end do
		else if ((evaluation_type .eq. 'rcs_p') .or. (evaluation_type .eq. 'rcs_n'))then ! 
			phi = evaluation_p%dim_b(1)
			i = 1		
			do j = 1, n_a
				theta = evaluation_p%dim_a(1) + dx_a*(j-1)
				r_local((i-1)*n_a+j)%point=(/sin(theta)*cos(phi),  &
						sin(theta)*sin(phi), cos(theta)/)*evaluation_p%dim_c
			end do
		else if ((evaluation_type .eq. 'BRDF_n') .or. (evaluation_type .eq. 'BRDF_p')) then ! 
			phi = evaluation_p%dim_b(1)
			i = 1		
			do j = 1, n_a
				theta = evaluation_p%dim_a(1) + dx_a*(j-1)
				r_local((i-1)*n_a+j)%point=(/sin(theta)*cos(phi),  &
						sin(theta)*sin(phi), cos(theta)/)*evaluation_p%dim_c
			end do
		end if				
		return		
	end subroutine get_evaluation_points
	
	subroutine get_evaluation_points_domain(evaluation_p, eps_r1, eps_r2, calculation_p, geometry_p, r_local)! 		
		type(evaluation_r_media), dimension(:), allocatable, intent(out) :: r_local
		complex(dp), intent(in) :: eps_r1, eps_r2
		integer:: calculation_p(8)
		type(lib_sie_evaluation_parameter_type), intent(in) :: evaluation_p
		real(dp), dimension(:), intent(in) :: geometry_p(2)

		!dummy
		real(dp ) :: dx_a, dx_b, ra(3)
		real(dp) :: phi, theta, wg, hg, dummy
		integer :: i, j, n_a, n_b, n_rs
		real(dp), dimension(:), allocatable :: xx, ff
		
		n_a = evaluation_p%N_dim(1)
		n_b = evaluation_p%N_dim(2)

		allocate(r_local(n_a*n_b))
		r_local%eps_r = eps_r1
		r_local%media = 1
		dx_a = (evaluation_p%dim_a(2) - evaluation_p%dim_a(1))/(n_a-1)
		
		if ((calculation_p(4) .ge. 4) .and. (calculation_p(4) .le. 7))then
			dx_b = 0.0
		else
			dx_b = (evaluation_p%dim_b(2) - evaluation_p%dim_b(1))/(n_b-1)
		end if
		
		if (calculation_p(4) .eq. 2)then ! xy-plot		
			do i = 1, n_b
				do j = 1, n_a
					r_local((i-1)*n_a+j)%point= &
					(/evaluation_p%dim_a(2) - dx_a*(j-1), evaluation_p%dim_b(2) - dx_b*(i-1), evaluation_p%dim_c/)
				end do
			end do
		else if (calculation_p(4) .eq. 3)then ! yz-plot      !     
			do i=1, n_b
				do j=1, n_a
					r_local((i-1)*n_a + j)%point=(/evaluation_p%dim_c, evaluation_p%dim_a(2) - dx_a*(j-1), &
					evaluation_p%dim_b(2) - dx_b*(i-1) /)
				end do
			end do
        
		else if (calculation_p(4) .eq. 1)then ! xz-plot
			if (calculation_p(5) .eq. 3)then !grating
				wg = geometry_p(1)
				hg = geometry_p(2)		
				do i=1, n_b
					do j=1, n_a
						r_local((i-1)*n_a+j)%point(1:3)=(/evaluation_p%dim_a(1) + dx_a*(j-1), evaluation_p%dim_c, &
						evaluation_p%dim_b(1) + dx_b*(i-1)/)
						ra = r_local((i-1)*n_a+j)%point					
						if ((abs(ra(1)) .le. wg/2) .and. (ra(3) .le. 0.0)) then
							r_local((i-1)*n_a+j)%eps_r = eps_r2
							r_local((i-1)*n_a+j)%media = 2							
						elseif ((abs(ra(1)) .ge. wg/2) .and. (ra(3) .le. -hg)) then
							r_local((i-1)*n_a+j)%eps_r = eps_r2		
							r_local((i-1)*n_a+j)%media = 2
						end if
					end do
				end do
			else if ((calculation_p(5) .eq. 1) .or. (calculation_p(5) .eq. 4)) then		!sphere	
				do i=1, n_b
					do j=1, n_a
						r_local((i-1)*n_a+j)%point(1:3)=(/evaluation_p%dim_a(1) + dx_a*(j-1), evaluation_p%dim_c, &
						evaluation_p%dim_b(1) + dx_b*(i-1)/)
						ra = r_local((i-1)*n_a+j)%point! - surface_center !attention: surface center at 0.0 is assumed. 
						if (norm2(ra) .le. geometry_p(1)) then
							r_local((i-1)*n_a+j)%eps_r = eps_r2
							r_local((i-1)*n_a+j)%media = 2
						end if
					end do
				end do			
			else if (calculation_p(5) .eq. 2)then		!rough surface		
				!The file 'f_boundary_rs.dat' is generated simultaneously when the rough surface is generated
				open(unit = 115, file = 'f_boundary_rs.dat', status = 'old', action = 'read')
				allocate(f_boundary_rs(n_a))
				do i = 1, n_a
					read(115, *) dummy, f_boundary_rs(i) 
				end do
				do j=1, n_a !the size of the surface is equal to the surface lengtgh, in x-direction
					do i=1, n_b !z-direction
						r_local((i-1)*n_a+j)%point(1:3)=(/evaluation_p%dim_a(1) + dx_a*(j-1), evaluation_p%dim_c, &
						evaluation_p%dim_b(1) + dx_b*(i-1)/)
						ra = r_local((i-1)*n_a+j)%point					
						if (ra(3) .ge. f_boundary_rs(j)) then						
							r_local((i-1)*n_a+j)%eps_r = eps_r2
							r_local((i-1)*n_a+j)%media = 2							
						end if
					end do
				end do			
			end if
		else if ((calculation_p(4) .ge. 4) .or. (calculation_p(4) .le. 7))then ! 
			phi = evaluation_p%dim_b(1)
			i = 1		
			do j = 1, n_a
				theta = evaluation_p%dim_a(1) + dx_a*(j-1)
				r_local((i-1)*n_a+j)%point=(/sin(theta)*cos(phi),  &				
						sin(theta)*sin(phi), cos(theta)/)*evaluation_p%dim_c
			end do		
		end if			
		return		
	end subroutine get_evaluation_points_domain
	
	function cross_r(v1,v2)!result(rv)
		! Function returning the cross product of
		! two real vectors v1 and v2
		real(dp), intent(in) :: v1(3), v2(3)
		real(dp) ::cross_r(3)
		cross_r(1)=v1(2)*v2(3)-v1(3)*v2(2)
		cross_r(2)=v1(3)*v2(1)-v1(1)*v2(3)
		cross_r(3)=v1(1)*v2(2)-v1(2)*v2(1)
		return
   end function   
	
	function cross_rc(v1,v2)!result(rv)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		real(dp), dimension(3), intent(in) :: v1
		complex(dp), dimension(3), intent(in) :: v2
		complex(dp), dimension(3) :: cross_rc
		
		cross_rc(1)=cmplx(v1(2))*v2(3)-cmplx(v1(3))*v2(2)
		cross_rc(2)=cmplx(v1(3))*v2(1)-cmplx(v1(1))*v2(3)
		cross_rc(3)=cmplx(v1(1))*v2(2)-cmplx(v1(2))*v2(1)
		return
	end function
	
	function vec_len(v) 
   ! Function of a vector norm
      real(dp), intent(in) :: v(3)
      real(dp) :: vec_len
      vec_len = sqrt(v(1)**2+v(2)**2+v(3)**2)     
      return
   end function vec_len
	
	function vec_len_c(v) 
   ! Function of norm for a complex vector
      complex(dp), intent(in) :: v(3)
      real(dp) :: vec_len_c
      vec_len_c = sqrt(abs(v(1))**2+abs(v(2))**2+abs(v(3))**2)     
      return
   end function vec_len_c  
	
	function cross_c(v1,v2)
   ! Function returning the cross product of
   ! two complex vectors v1 and v2
		complex(dp), dimension(3), intent(in) :: v1, v2
		complex(dp), dimension(3) :: cross_c
      cross_c(1)=v1(2)*v2(3)-v1(3)*v2(2)
      cross_c(2)=v1(3)*v2(1)-v1(1)*v2(3)
      cross_c(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
   end function
	
	function cross_rc_arr(n, v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		integer, intent(in) :: n
		real(dp), dimension(n, 3), intent(in) :: v1
		complex(dp), dimension(n, 3), intent(in) :: v2
		complex(dp), dimension(n, 3) :: cross_rc_arr
		
		cross_rc_arr(:, 1)=cmplx(v1(:, 2))*v2(:, 3)-cmplx(v1(:, 3))*v2(:, 2)
		cross_rc_arr(:, 2)=cmplx(v1(:, 3))*v2(:, 1)-cmplx(v1(:, 1))*v2(:, 3)
		cross_rc_arr(:, 3)=cmplx(v1(:, 1))*v2(:, 2)-cmplx(v1(:, 2))*v2(:, 1)
		return
	end function
	
	function dot_rr_arr(n, v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		integer, intent(in) :: n
		real(dp), dimension(n, 3), intent(in) :: v1
		real(dp), dimension(3), intent(in) :: v2
		real(dp), dimension(n) :: dot_rr_arr
		
		dot_rr_arr(:)= v1(:, 1)*v2(1) + v1(:, 2)*v2(2) + v1(:, 3)*v2(3)
		return
	end function
	
	function dot_arr(v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2		
		real(dp), dimension(:, :), intent(in) :: v1
		real(dp), dimension(3), intent(in) :: v2
		real(dp), dimension(:), allocatable :: dot_arr
		integer :: n
		n = size(v1,1)
		allocate(dot_arr(n))
		dot_arr(:)= v1(:, 1)*v2(1) + v1(:, 2)*v2(2) + v1(:, 3)*v2(3)
		return
	end function
	
	function dot_product_arr(n, v1, v2)
		integer, intent(in) :: n
		complex(dp), dimension(n, 3) :: v1, v2
		complex(dp), dimension(n) :: dot_product_arr
		
		dot_product_arr(1:n) = v1(1:n, 1)*v2(1:n, 1) + v1(1:n, 2)*v2(1:n, 2) + v1(1:n, 3)*v2(1:n, 3)
		return
	end function
	
	subroutine solver_LU(Mr, MM_s, VV_EH)
      integer, parameter :: NRHS = 1 
      integer, intent(in) :: Mr !edge number 
      integer :: LDA, LDB
      integer :: info   
      
      complex(dp), dimension(Mr, Mr), intent(in) :: MM_s
      complex(dp), dimension(Mr), intent(in out) :: VV_EH   
      integer ::   IPIV(Mr)
      external   ZGESV ! ZGELS
		external   CGESV ! ZGELS
        LDA = Mr
        LDB = Mr

!     Solve the equations A*X = B.
		if (dp == 4) then
			CALL CGESV(Mr, NRHS, MM_s, LDA, IPIV, VV_EH, LDB, INFO)
		else if (dp == 8) then
			CALL ZGESV(Mr, NRHS, MM_s, LDA, IPIV, VV_EH, LDB, INFO)
		end if
		
		IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagon_al element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
		END IF
    end subroutine solver_LU
	 
	subroutine sort_edge_index(vector, v_size, v_ascendente)
		implicit none
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

end module lib_sie_math