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
	
module lib_sie_math
	
	use lib_sie_constants
	use lib_sie_type_function		
	use Gauss_Quadrature_Fast
	
	use libmath
	
	implicit none
	
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
	
	subroutine get_evaluation_points(evaluation_p, evaluation_type, r_local)! 		
		character(len = 10), intent(in) :: evaluation_type
		type(lib_sie_evaluation_parameter_type), intent(in) :: evaluation_p
		type(point), dimension(:), allocatable, intent(out) :: r_local

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
		else if (evaluation_type .eq. 'BRDF')then ! 
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

	!kd: wavevector
	subroutine Gaussian_beam_focused(w0, r_local, E_00, kd)	
		implicit none		
		
		real(dp), intent(in) :: r_local(3)
		type(vector_c), intent(out) :: E_00
		real(dp), intent(in) :: w0
		complex(dp), intent(in) :: kd
		
		!dummy		
		real (dp) :: f0, f, theta_max, a, x_j
		integer(kind = 4) :: m, nm, i, ngp!
		real(dp), dimension(:), allocatable :: jn, djn, xx, ww
		real(dp) :: c1, c2, rho, phi_r, fw, j0, j1, j2, c3
		real(dp) :: E0, r_max, phi_m, z_max
		real(dp) :: indx_2, n_a
		
		!dummy			
		real(dp) :: r0(3) 
		complex(dp) :: p1, p2, p3, I_00, I_01, I_02, trf, E_tmp(3)    
				
		ngp = 7
		nm = 5
		m = 3
		r_max = 2e-6 !radius of the field to be calculated
		phi_m = 2*PI
		z_max = 5e-6 !half of the calculation range along the z-axis
			
		r0 = (/0.0, 0.0, 0.0/) !center position of the focused point		
		E0 = 1.0
		indx_2 = 1.0
		n_a = 0.8
		a = 0.0
		theta_max = asin(n_a/indx_2)    
		f = 30e-6
		f0 = w0/(f*sin(theta_max))    
		
		allocate(jn(1:m))
		allocate(djn(1:m))    
		allocate(ww(1:ngp))
		allocate(xx(1:ngp))
		call legendre_handle(ngp, a, theta_max, xx, ww)
		!Ex ploarized light focused onto z = 0
		!Eq. (3.49 - 3.66) Novotny, n_ano-Optics
		rho = sqrt((r_local(1) - r0(1))**2 + (r_local(2) - r0(2))**2)
		if (rho .eq. 0.0) then
			phi_r = 0.0
		else 
			phi_r = acos((r_local(1) - r0(1))/rho)		
		end if
		
		I_01 = (0.0, 0.0)
		I_02 = (0.0, 0.0)
		I_00 = (0.0, 0.0)
		do i = 1, ngp        
			c1 = sin(xx(i))
			c2 = cos(xx(i))
			c3 = sin(theta_max)
			fw = exp(-1/f0**2*(c1/C3)**2)             
			x_j = c1*kd*rho
			jn = lib_math_bessel_spherical_first_kind(x_j, 0, nm)
			j0 = jn(1)
			j1 = jn(3)
			j2 = jn(3)              
			trf = fw*sqrt(c2)*exp(im*kd*(r_local(3)-r0(3))*c2)*ww(i)
			I_00 = I_00 + c1*(1 + c2)*trf*j0!
			I_01 = I_01 + c1**2*j1*trf        
			I_02 = I_02 + c1*(1-c2)*j2*trf
		end do        
		p1 = I_00 + I_02*cos(2*phi_r) 
		p2 = I_02*sin(2*phi_r)
		p3 = -2*im*I_01*cos(phi_r)
		E_tmp = (/p1, p2, p3/) 
		E_00%vector = im*kd*f/2*sqrt(1.0/indx_2)*E0*exp(-im*kd*f)*E_tmp
		return
	end subroutine Gaussian_beam_focused
	
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
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
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
		DO i = 1, v_size
			tmp(i) = MINVAL(vector,mk)
			mk(MINLOC(vector,mk)) = .FALSE.
		END DO
		v_ascendente = tmp 
		return		
	end subroutine
	
end module lib_sie_math