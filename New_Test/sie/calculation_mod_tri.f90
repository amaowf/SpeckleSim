module calculation_mod_tri
	
	! Module containing functions useful in the
	! calculation part of the program
	!use omp_lib
	use disc_mod
	use input_tri
	use libmath
	implicit none
	
    contains
	
	function cross_c(v1,v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		complex(kind = dp), dimension(3), intent(in) :: v1, v2
		complex(kind = dp), dimension(3) :: cross_c
		cross_c(1)=v1(2)*v2(3)-v1(3)*v2(2)
		cross_c(2)=v1(3)*v2(1)-v1(1)*v2(3)
		cross_c(3)=v1(1)*v2(2)-v1(2)*v2(1)
		return
	end function
	
	function cross_rc(v1,v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		real(kind = dp), dimension(3), intent(in) :: v1
		complex(kind = dp), dimension(3), intent(in) :: v2
		complex(kind = dp), dimension(3) :: cross_rc
		
		cross_rc(1)=cmplx(v1(2))*v2(3)-cmplx(v1(3))*v2(2)
		cross_rc(2)=cmplx(v1(3))*v2(1)-cmplx(v1(1))*v2(3)
		cross_rc(3)=cmplx(v1(1))*v2(2)-cmplx(v1(2))*v2(1)
		return
	end function
	
	function cross_rc_arr(n, v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		integer, intent(in) :: n
		real(kind = dp), dimension(n, 3), intent(in) :: v1
		complex(kind = dp), dimension(n, 3), intent(in) :: v2
		complex(kind = dp), dimension(n, 3) :: cross_rc_arr
		
		cross_rc_arr(:, 1)=cmplx(v1(:, 2))*v2(:, 3)-cmplx(v1(:, 3))*v2(:, 2)
		cross_rc_arr(:, 2)=cmplx(v1(:, 3))*v2(:, 1)-cmplx(v1(:, 1))*v2(:, 3)
		cross_rc_arr(:, 3)=cmplx(v1(:, 1))*v2(:, 2)-cmplx(v1(:, 2))*v2(:, 1)
		return
	end function
	
	function dot_rr_arr(n, v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2
		integer, intent(in) :: n
		real(kind = dp), dimension(n, 3), intent(in) :: v1
		real(kind = dp), dimension(3), intent(in) :: v2
		real(kind = dp), dimension(n) :: dot_rr_arr
		
		dot_rr_arr(:)= v1(:, 1)*v2(1) + v1(:, 2)*v2(2) + v1(:, 3)*v2(3)
		return
	end function
	
	function dot_arr(v1, v2)
	 ! Function returning the cross product of
	 ! two complex vectors v1 and v2		
		real(kind = dp), dimension(:, :), intent(in) :: v1
		real(kind = dp), dimension(3), intent(in) :: v2
		real(kind = dp), dimension(:), allocatable :: dot_arr
		integer :: n
		n = size(v1,1)
		allocate(dot_arr(n))
		dot_arr(:)= v1(:, 1)*v2(1) + v1(:, 2)*v2(2) + v1(:, 3)*v2(3)
		return
	end function
	
	function dot_product_arr(n, v1, v2)
		integer, intent(in) :: n
		complex(kind = dp), dimension(n, 3) :: v1, v2
		complex(kind = dp), dimension(n) :: dot_product_arr
		
		dot_product_arr(1:n) = v1(1:n, 1)*v2(1:n, 1) + v1(1:n, 2)*v2(1:n, 2) + v1(1:n, 3)*v2(1:n, 3)
		return
	end function
	
	function cross_r(v1,v2)
	 ! Function returning the cross product of
	 ! two real vectors v1 and v2
		real(kind = dp), dimension(3), intent(in) :: v1, v2 
		real(kind = dp), dimension(3) :: cross_r		
		cross_r(1)=v1(2)*v2(3)-v1(3)*v2(2)
		cross_r(2)=v1(3)*v2(1)-v1(1)*v2(3)
		cross_r(3)=v1(1)*v2(2)-v1(2)*v2(1)     
		return
	end function
 
	subroutine incident_plane_tri(struct, m_pair, ngp_a, Sum_EH)
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) ::  Sum_EH(2)
		integer, intent(in) :: m_pair, ngp_a 
		
		complex(kind = dp) :: Ec_in(3), Hc_in(3)
		real(kind = dp) :: lng
		integer :: i, ot, pos_neg
		real(kind = dp) :: R_a(3), p1(3), p2(3), pc(3), fm(3)!, TRF, 
		real(kind = dp) :: area	
		real(kind = dp), dimension(100) :: a, b, w
		complex(kind = dp) :: k_hat(3)
		k_hat = k_in_hat*(1 + im*0)
		
		call Quadrature_tri(ngp_a, a, b, w)
		Sum_EH(1:2) = (/(0.0, 0.0),  (0.0, 0.0)/)
		!area is cancelled out with that in w(i)*area
		do ot = 1, 2
			call fn_parameter(struct, m_pair, ot, p1, p2, pc, lng, area)
			pos_neg = -1*(((ot-1)*ot)-1) 
			do i = 1, ngp_a 
				R_a = r_simplex(a(i), b(i), p1, p2, pc)
				fm = pos_neg*(R_a-pc)*(lng/2)	
				Ec_in = E_in_hat*exp(-im*dot_product(k1*k_in_hat, R_a))	
				Hc_in = cross_c(k_hat, Ec_in)/(my_1*c0)
				Sum_EH(1) = Sum_EH(1) + dot_product(fm, Ec_in)*w(i)!
				Sum_EH(2) = Sum_EH(2) + dot_product(fm, Hc_in)*w(i)!
			end do !loop i		
		end do
		return
	end subroutine incident_plane_tri  
	
	subroutine incident_Gauss_tri(struct, m_pair, ngp_a, Sum_EH)
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) ::  Sum_EH(2)
		integer, intent(in) :: m_pair, ngp_a 
		
		complex(kind = dp) :: Ec_in(3), Hc_in(3)
		real(kind = dp) :: lng
		integer :: i, ot, pos_neg
		real(kind = dp) :: R_a(3), p1(3), p2(3), pc(3), fm(3)!, TRF, 
		real(kind = dp) :: area	
		real(kind = dp), dimension(100) :: a, b, w
		complex(kind = dp) :: k_hat(3)
		k_hat = k_in_hat*(1 + im*0)
		
		call Quadrature_tri(ngp_a, a, b, w)
		Sum_EH(1:2) = (/(0.0, 0.0),  (0.0, 0.0)/)
		!area is cancelled out with that in w(i)*area
		do ot = 1, 2
			call fn_parameter(struct, m_pair, ot, p1, p2, pc, lng, area)		
			pos_neg = -1*(((ot-1)*ot)-1) 
			do i = 1, ngp_a 
				R_a = r_simplex(a(i), b(i), p1, p2, pc)
				fm = pos_neg*(R_a-pc)*(lng/2)	
				call Gaussian_beam_theta(w0, R_a, Ec_in)
				Hc_in = cross_c(k_hat, Ec_in)/(my_1*c0)
				Sum_EH(2) = Sum_EH(2) + dot_product(fm, Hc_in)*w(i)!
			end do !loop i		
		end do
		return
	end subroutine incident_Gauss_tri  
	
	subroutine incident_FocusedGauss_tri(struct, m_pair, ngp_a, Sum_EH)
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) ::  Sum_EH(2)
		integer, intent(in) :: m_pair, ngp_a 
		
		complex(kind = dp) :: Ec_in(3), Hc_in(3)
		real(kind = dp) :: lng
		integer :: i, ot, pos_neg
		real(kind = dp) :: R_a(3), p1(3), p2(3), pc(3), fm(3)!
		real(kind = dp) :: area	
		real(kind = dp), dimension(100) :: a, b, w
		complex(kind = dp) :: k_hat(3)
		k_hat = k_in_hat*(1 + im*0)
		
		call Quadrature_tri(ngp_a, a, b, w)
		Sum_EH(1:2) = (/(0.0, 0.0),  (0.0, 0.0)/)		
		!area is cancelled out with that in w(i)*area
		do ot = 1, 2
			call fn_parameter(struct, m_pair, ot, p1, p2, pc, lng, area)		
			pos_neg = -1*(((ot-1)*ot)-1) 
			do i = 1, ngp_a 
				R_a = r_simplex(a(i), b(i), p1, p2, pc)
				fm = pos_neg*(R_a-pc)*(lng/2)	
				call Focused_Gaussian_beam(w0, R_a, Ec_in)
				Hc_in = cross_c(k_hat, Ec_in)/(my_1*c0)
				Sum_EH(2) = Sum_EH(2) + dot_product(fm, Hc_in)*w(i)!
			end do !loop i		
		end do
		return
	end subroutine incident_FocusedGauss_tri  
 
	function vec_len(v)
		real(kind = dp), intent(in) :: v(3)
		real(kind = dp) :: vec_len
		vec_len = sqrt(v(1)**2+v(2)**2+v(3)**2)
		return
	end function
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Liwei Fu, 22.03.2019
	
	function normal_vec_r(p1, p2, p3)
	 ! Function returning the normal vector
	 ! (zero imaginary part) of a discretization element
	 ! pointing outwards from origo when p1,p2,p3 are
	 ! corners of the element.
		real(kind = dp), intent(in) :: p1(3), p2(3), p3(3)
		real(kind = dp) :: normal_vec_r(3)
		real(kind = dp) :: v1(3), v2(3)
		 
		v1(1:3)=p1(1:3)-p2(1:3) 
		v2(1:3)=p3(1:3)-p2(1:3)				
		normal_vec_r = cross_r(v1,v2)/vec_len(cross_r(v1,v2))
		!mid = (p1 + p2 + p3)/3.0
		!mid_n = mid/vec_len(mid)
		!tmp = dot_product(normal_vec_r, mid)
		!if (tmp < 0)then
		!	print*, 'tmp =', tmp
		!	print*, 'Wrong direction of the normal vector!'
		!	print*, 'p1 =', p1
		!	print*, 'p2 =', p2
		!	print*, 'p3 =', p3
		!end if
		return
	end function
	 
	subroutine m_vector(n_vec, p1, p2, p3, m_v)
        real(kind = dp), dimension(3), intent(in) :: p1, p2, p3, n_vec        
        real(kind = dp), intent(out) :: m_v(3, 3)
        m_v(1, :) = cross_r(n_vec, p3-p2)/vec_len(cross_r(n_vec, p3-p2))
        m_v(2, :) = cross_r(n_vec, p1-p3)/vec_len(cross_r(n_vec, p1-p3))
        m_v(3, :) = cross_r(n_vec, p2-p1)/vec_len(cross_r(n_vec, p2-p1))        
        return
    end subroutine 
	
	function Iq_L(q, Ra, p1, p2)
        real(kind = dp), intent(in) :: Ra(3), p1(3), p2(3)
		integer, intent(in) :: q   
        real(kind = dp) :: IL_n1, IL_p1, IL_p3, IL_p5, Iq_L, &
			Sum_p, Sum_n, Sp, Sn, Rp, Rn, S_vec(3) !negative and positive for n and p
        real(kind = dp) :: Dif_p, Dif_n, R0
        intrinsic :: log
 
        s_vec = (p2-p1)/vec_len(p2-p1)		
        Sp = dot_product((p2 - Ra), s_vec)
        Sn = dot_product((p1 - Ra), s_vec)		
        Rn = vec_len(Ra - p1)
        Rp = vec_len(Ra - p2)
		
        Sum_p = Rp + Sp !eq. (41-43) in Hänninen
        Sum_n = Rn + Sn
        Dif_p = Rp - Sp
        Dif_n = Rn - Sn
        R0 = sqrt(Rp**2 -Sp**2)
                        
        if (R0 .eq. 0.0) then
             IL_n1 = log(Sp/Sn)
        else  
            if (abs(Dif_p) > abs(Sum_n)) then
                IL_n1 = log(Dif_n/Dif_p)
            else
                IL_n1 = log(Sum_p/Sum_n)
            end if
        end if             
        
        IL_p1 = 0.5*R0**2*IL_n1 + 0.5*(Sp*Rp-Sn*Rn) !q = 1
        IL_p3 = 0.75*R0**2*IL_p1 + 0.25*(Sp*Rp**3-Sn*Rn**3) !q = 3
        IL_p5 = 5.0/6.0*R0**2*IL_p3 + 1.0/6.0*(Sp*Rp**5-Sn*Rn**5) !q = 5
		
        if (q == -1) then
            Iq_L = IL_n1
        else if (q == 1)then
            Iq_L = IL_p1
        else if (q == 3)then
            Iq_L = IL_p3
        else if (q == 5)then
            Iq_L = IL_p5
		else 
		print*,	'q is not in the given range!'
        end if
        return   
    end function Iq_L
	
	function Iq_s(q, R_a, p1, p2, p3, Iqs_l)
        real(kind = dp), dimension(3), intent(in) :: R_a, p1, p2, p3
		integer, intent(in) :: q		
        real(kind = dp) :: Iq_s, m_r(3, 3), mid(3), nul
        real(kind = dp) :: Omeg, x, y, sum_IL, Is_n3, Is_n1, Is_p1, h  !rp is the r' in Eq.(45). Hänninen's paper
        real(kind = dp), dimension(3) :: a1, a2, a3, n_q, t_v, IL
        integer :: i, n
        logical :: Iqs_l
		
		nul = 0.0
				
		n_q = normal_vec_r(p1, p2, p3)
        a1 = (p1 - R_a)/vec_len(p1 - R_a)
        a2 = (p2 - R_a)/vec_len(p2 - R_a)
        a3 = (p3 - R_a)/vec_len(p3 - R_a)
		mid = (p1+p2+p3)/3
		
        Is_n3 = nul
        x = 1 + dot_product(a1, a2) + dot_product(a1, a3) + dot_product(a2, a3)
        y = abs(dot_product(a1, cross_r(a2, a3)))
		
		
		
        Omeg = 2*atan2(y, x)				
        h = -dot_product(n_q, (R_a - mid)) !			
		
		if (h == nul) then
			Is_n3 = nul
		else 
			if (Iqs_l .eqv. .true.) then
				Is_n3 = nul
			else if ((Omeg <= PI) .and. (Omeg > -PI)) then
				Is_n3 = 1/h*Omeg
			else
				print*, 'Omega is outside the range'			
				print*, 'Ra =', R_a
				print*, 'p1 =', p1
				print*, 'p2 =', p2		
				print*, 'p3 =', p3		
				!print*, 'Omega=', Omeg
				!print*, 'x =', x
				!print*, 'y =', y 
			end if	
		end if
		
		call m_vector(n_q, p1, p2, p3, m_r)	
		
       	t_v(1) = dot_product(m_r(1, :), R_a-(p2+p3)/2)  !corresponding to the edge index in m_r
		t_v(2) = dot_product(m_r(2, :), R_a-(p3+p1)/2) 
		t_v(3) = dot_product(m_r(3, :), R_a-(p1+p2)/2)         
		
		n = -1
		sum_IL = 0.0
		IL(1) = Iq_L(n, R_a, p2, p3)
		IL(2) = Iq_L(n, R_a, p3, p1) 
		IL(3) = Iq_L(n, R_a, p1, p2)
        do i = 1, 3
			sum_IL = sum_IL + t_v(i)*IL(i)
		end do 
		Is_n1 = 1/(real(n)+2)*(n*h**2*Is_n3 - sum_IL)
 
		n = 1
		sum_IL = 0.0
		IL(1) = Iq_L(n, R_a, p2, p3)
		IL(2) = Iq_L(n, R_a, p3, p1) 
		IL(3) = Iq_L(n, R_a, p1, p2)
		
        do i = 1, 3
			sum_IL = sum_IL + t_v(i)*IL(i)
		end do 		
		Is_p1 = 1/(real(n)+2)*(real(n)*h**2*Is_n1 - sum_IL )
  
		if ((Iqs_l .eqv. .true.) .and. (q > -1)) then
			Iq_s = 0.0
		else 
			if (q ==-3) then 
				Iq_s = Is_n3
			else if (q==-1) then
				Iq_s = Is_n1
			else if (q == 1) then
				Iq_s = Is_p1
			else 
				print*, 'q is beyond the given range!'
			end if
		end if
        return    
    end function Iq_s
	
    function K2q(q, R_a, p1, p2, p3, Iqs_L)
        real(kind = dp), intent(in) :: p1(3), p2(3), p3(3), R_a(3)
        integer, intent(in) :: q		
        logical, intent(in) :: Iqs_L		
        real(kind = dp), dimension(3) :: rho, K2q, IL_arr, n_q
        real(kind = dp) :: Iqs_q, sum_IL(3), mid(3), h, m_v(3, 3)
        integer :: i   
		
        K2q(1:3) = 0.0
		n_q = normal_vec_r(p1, p2, p3)
		mid = (p1+p2+p3)/3
		h = -dot_product(n_q, (R_a - mid)) !
        call m_vector(n_q, p1, p2, p3, m_v)
        
        IL_arr(1) = Iq_L(q+2, R_a, p3, p2)
        IL_arr(2) = Iq_L(q+2, R_a, p1, p3)
        IL_arr(3) = Iq_L(q+2, R_a, p2, p1) 
		
        Iqs_q = Iq_s(q, R_a, p1, p2, p3, Iqs_L)
		
        sum_IL = 0.0
        do i = 1, 3
            sum_IL = sum_IL + IL_arr(i)*m_v(i, :)
        end do 
        rho = R_a + h*n_q !very strange, since the sign in between is oposite to that in Hännien's paper
        K2q = 1/(real(q)+2)*sum_IL + (rho - p3)*Iqs_q
        return
    end function K2q
    !
    function K4q(q, R_a, p1, p2, p3, Iqs_L)
        real(kind = dp), intent(in) :: p1(3), p2(3), p3(3), R_a(3)
        logical, intent(in) :: Iqs_L
		
        real(kind = dp) :: n_q(3), I_arr(3), K4q(3), K3q(3)
        real(kind = dp) :: Is, h, mid(3), m_v(3, 3), sum_tmp(3)
        integer :: i, q
                
        n_q = normal_vec_r(p1, p2, p3)
        call m_vector(n_q, p1, p2, p3, m_v)
		if (Iqs_L .eqv. .true.) then
			K4q = (/0.0, 0.0, 0.0/)
		else 
			!R = R_a - R_q    
			mid = (p1+p2+p3)/3
			h = -dot_product(n_q, (R_a - mid))
			I_arr(1) = Iq_L(q, R_a, p2, p3)!q is not defined yet.
			I_arr(2) = Iq_L(q, R_a, p3, p1)
			I_arr(3) = Iq_L(q, R_a, p1, p2)
			Is = Iq_s(q-2, R_a, p1, p2, p3, Iqs_L)
			sum_tmp = (/0.0, 0.0, 0.0/)
			do i = 1, 3
				sum_tmp(:) = sum_tmp(:) +  m_v(i, :)*I_arr(i) 
			end do			
			K3q = sum_tmp + h*real(q)*n_q*Is !
			if (q==1)then				
				K4q = cross_r((R_a-p3), K3q) ! Without the negative sign in Eq.(76) (Hänninen) yet
			else if (q==-1)then
				K4q = -cross_r((R_a-p3), K3q)!
			else 
				print*, 'q is outside the given range!'
			end if			
		end if
        return		
    end function K4q
    !
    subroutine Diver_Green(km, R_norm, Greenf, DGf)
        implicit none
        complex(kind = dp), intent(in) :: km
        real(kind = dp), intent(in) :: R_norm
        complex(kind = dp), intent(out) :: DGf, Greenf
        complex(kind = dp) :: Exp_f
        intrinsic :: exp
        Exp_f = exp(-im*km*R_norm) ! It should be norm(R). m: medium
        DGf = -(1 + im*km*R_norm)*Exp_f/(4*PI*R_norm**2)!r' by Mitschang, without the direction unit vector
        Greenf = Exp_f/(4*PI*R_norm) !
        return 
    end subroutine Diver_Green
    
    !smoothed part of Green's function
    subroutine Diver_Green_s(km, R_vl, Gs, DGs)
        implicit none
		real(kind = dp), intent(in) :: R_vl
        complex(kind = dp), intent(in) :: km 
        complex(kind = dp), intent(out) :: DGs, Gs
        complex(kind = dp) :: Exp_f
        intrinsic :: exp
        Exp_f = exp(-im*km*R_vl) ! It should be norm(R). m: medium        
		if (R_vl == 0.0) then
        Gs = -im*km/(4*PI)
		DGs = km**2/(8*PI) 
		else 
		Gs = (Exp_f-1)/(4*PI*R_vl) + km**2*R_vl/(8*PI) !
		DGs = (1 -(1 + im*km*R_vl)*Exp_f)/(4*PI*R_vl**2) + km**2/(8*PI) 
		end if
        return 
    end subroutine Diver_Green_s 
   
    subroutine fn_parameter(struct, m, ot, p1, p2, pc, lng, area)
      implicit none
		type(Structure), intent(in) :: struct
		integer, intent(in) :: m, ot
      real(kind = dp), dimension(3), intent(out) :: p1, p2, pc
      real(kind = dp), intent(out) :: lng, area
      integer :: j, nv, p
		type(Pair) :: pair_p
		
        pair_p = struct%neighbours(m)
        do j = 1, 3                
            if (struct%elements(pair_p%elements(ot))%corners(j) &
                /= pair_p%corners(1) &
                .and. struct%elements(pair_p%elements(ot))%corners(j) &
                /= pair_p%corners(2))then
                pc(1:3)=struct%points(struct%elements(pair_p%elements(ot))&
                %corners(j))%point ! vertex of the triangle
				nv = j
            end if
        end do
		p = struct%neighbours(m)%elements(ot) !element index
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
        lng = vec_len(p1-p2)
		area = 0.5*vec_len(cross_r((p1-pc), (p2-pc))) 
        return
    end subroutine fn_parameter
	 
	 subroutine fn_parameter_midpoint(struct, m, ot, p1, p2, pc, rm, lng, area)
      implicit none
		type(Structure), intent(in) :: struct
		integer, intent(in) :: m, ot
      real(kind = dp), dimension(3), intent(out) :: p1, p2, pc, rm
      real(kind = dp), intent(out) :: lng, area
      integer :: j, nv, p
		type(Pair) :: pair_p
		
        pair_p = struct%neighbours(m) !edge index
        do j = 1, 3                
            if (struct%elements(pair_p%elements(ot))%corners(j) &
                /= pair_p%corners(1) &
                .and. struct%elements(pair_p%elements(ot))%corners(j) &
                /= pair_p%corners(2))then
                pc(1:3)=struct%points(struct%elements(pair_p%elements(ot))&
                %corners(j))%point ! vertex of the triangle
				nv = j
            end if
        end do
		p = struct%neighbours(m)%elements(ot) !element index
		rm = struct%midpoints(p)%point
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
        lng = vec_len(p1-p2)
		area = 0.5*vec_len(cross_r((p1-pc), (p2-pc))) 
        return
    end subroutine fn_parameter_midpoint
	 
	 subroutine fn_edge_center(struct, m, ot, p_edgec)
		type(Structure), intent(in) :: struct
		integer, intent(in) :: m, ot
		real(kind = dp), intent(out) :: p_edgec(3)
		real(kind = dp), dimension(3) :: p1, p2, pc
		
		integer :: j, nv, p
		type(Pair) :: pair_p
		
      pair_p = struct%neighbours(m)
      do j = 1, 3                
         if (struct%elements(pair_p%elements(ot))%corners(j) &
               /= pair_p%corners(1) &
               .and. struct%elements(pair_p%elements(ot))%corners(j) &
               /= pair_p%corners(2))then
               pc(1:3)=struct%points(struct%elements(pair_p%elements(ot))&
               %corners(j))%point ! vertex of the triangle
			nv = j
         end if
      end do
		p = struct%neighbours(m)%elements(ot) !element index
		if (nv == 1) then		
			p1 = struct%points(struct%elements(p)%corners(2))%point
			p2 = struct%points(struct%elements(p)%corners(3))%point
		else if (nv == 2) then
			p1 = struct%points(struct%elements(p)%corners(3))%point
			p2 = struct%points(struct%elements(p)%corners(1))%point
		else 
			p1 = struct%points(struct%elements(p)%corners(1))%point
			p2 = struct%points(struct%elements(p)%corners(2))%point
		end if		
			p_edgec = (p1 + p2)/2 !midpoint of the common edge
      return
		
	end subroutine 
  	subroutine Normal_Integration_Compact(m, n, ngp, struct, sum_a, sum_b, sum_c)
		implicit none		
		integer, intent(in) :: m, n, ngp		
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) :: sum_a, sum_b, sum_c
		
		real(kind = dp) :: TRF, Len_m, Len_n, area!, area_n
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: Rq(3), Ra(3), fm(3), fn(3), dfm, dfn, pm1(3), pm2(3), & 
			vm_t(3), vn_l(3), pn1(3), pn2(3), Rn(3), R(3) 
		complex(kind = dp) :: f1, f2, DG_1, DG_2, G_1, G_2
		integer :: i, j,  pos_neg, tt, ll
		intrinsic :: sqrt
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)	    
		
		!In this calculation, all the areas are cancelled out.		
		call Quadrature_tri(ngp, a, b, w)        
		
		do tt = 1, 2
			do ll = 1, 2
				call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
				pos_neg = -1*(((tt-1)*tt)-1)            
				dfm = pos_neg*Len_m!            
				call fn_parameter(struct, n, ll, pn1, pn2, vn_l, Len_n, area)
				pos_neg = -1*(((ll-1)*ll)-1)                
				dfn = pos_neg*Len_n
				do i = 1, ngp
					Ra = r_simplex(a(i), b(i), pm1, pm2, vm_t)
					fm = 0.5*(Ra-vm_t)*dfm
					do j = 1, ngp 
						Rq = r_simplex(a(j), b(j), pn1, pn2, vn_l)
						fn = 0.5*(Rq-vn_l)*dfn
						R = Ra - Rq				
						Rn = R/vec_len(R)
						call Diver_Green(k1, vec_len(R), G_1, DG_1)
						call Diver_Green(k2, vec_len(R), G_2, DG_2)						
						TRF = w(i)*w(j)
						sum_a = sum_a + dot_product(fm, cross_r(fn, Rn)*(DG_1 + DG_2))*TRF !K_21N, very sensitive to Rn
						f1 = im*(Omega*my_1*dot_product(fm, fn) - 1/(Omega*eps_1)*dfm*dfn)*G_1
						f2 = im*(Omega*my_2*dot_product(fm, fn) - 1/(Omega*eps_2)*dfm*dfn)*G_2					
			
						sum_b = sum_b + (f1 + f2)*TRF !D_11  					
						f1 = im*(Omega*eps_1*dot_product(fm, fn) - 1/(Omega*my_1)*dfm*dfn)*G_1 
						f2 = im*(Omega*eps_2*dot_product(fm, fn) - 1/(Omega*my_2)*dfm*dfn)*G_2 
						sum_c = sum_c + (f1 + f2)*TRF !D_22 
					end do !loop j
				end do  !loop i
			end do !ll-loop
		end do !tt-loop
	end subroutine Normal_Integration_Compact
	
  	subroutine Normal_Integration_new(m, n, tt, ll, ngp, struct, sum_a, sum_b, sum_c)
		implicit none		
		integer, intent(in) :: m, n, ngp, tt, ll		
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) :: sum_a, sum_b, sum_c
		
		real(kind = dp) :: TRF, Len_m, Len_n, area!, area_n
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = dp) :: Rq(3), Ra(3), fm(3), fn(3), dfm, dfn, pm1(3), pm2(3), & 
			vm_t(3), vn_l(3), pn1(3), pn2(3), Rn(3), R(3) 
		complex(kind = dp) :: f1, f2, DG_1, DG_2, G_1, G_2
		integer :: i, j,  pos_neg
		intrinsic :: sqrt
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)	    
		
		!In this calculation, all the areas are cancelled out.		
		call Quadrature_tri(ngp, a, b, w)        
		call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
		pos_neg = -1*(((tt-1)*tt)-1)            
		dfm = pos_neg*Len_m!
		!print*, 'dfm=', dfm
            
		call fn_parameter(struct, n, ll, pn1, pn2, vn_l, Len_n, area)
		pos_neg = -1*(((ll-1)*ll)-1)                
		dfn = pos_neg*Len_n
		!print*, 'dfn=', dfn
		
		do i = 1, ngp
			Ra = r_simplex(a(i), b(i), pm1, pm2, vm_t)
			fm = 0.5*(Ra-vm_t)*dfm
			do j = 1, ngp 
				Rq = r_simplex(a(j), b(j), pn1, pn2, vn_l)
				fn = 0.5*(Rq-vn_l)*dfn
				R = Ra - Rq				
				Rn = R/vec_len(R)
				call Diver_Green(k1, vec_len(R), G_1, DG_1)
				call Diver_Green(k2, vec_len(R), G_2, DG_2)						
				TRF = w(i)*w(j)
				sum_a = sum_a + dot_product(fm, cross_r(fn, Rn)*(DG_1 + DG_2))*TRF !K_21N, very sensitive to Rn
				f1 = im*(Omega*my_1*dot_product(fm, fn) - 1/(Omega*eps_1)*dfm*dfn)*G_1
				f2 = im*(Omega*my_2*dot_product(fm, fn) - 1/(Omega*eps_2)*dfm*dfn)*G_2					
			
				sum_b = sum_b + (f1 + f2)*TRF !D_11  					
				f1 = im*(Omega*eps_1*dot_product(fm, fn) - 1/(Omega*my_1)*dfm*dfn)*G_1 
				f2 = im*(Omega*eps_2*dot_product(fm, fn) - 1/(Omega*my_2)*dfm*dfn)*G_2 
				sum_c = sum_c + (f1 + f2)*TRF !D_22 
				
			end do !loop j
		end do  !loop i
	end subroutine Normal_Integration_new
	
	subroutine Singular_Integration_Compact(m, n, ngp, struct, suma, sumb, sumc)
		implicit none
		integer, intent(in) :: m, n, ngp
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		
		integer :: i, j, pos_neg_m, pos_neg_n, q, pp, qq, tt, ll
        real(kind = dp) :: TRF, TRF_m, Len_m, Len_n  
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
        real(kind = dp) :: Rq(3), Rn(3), Ra(3), R(3), fm(3), fn(3), dfm, dfn, pm1(3), pm2(3), vm_t(3), vn_l(3), pn1(3), pn2(3)
		real(kind = dp) :: k2_n1(3), k2_p1(3), k4_n1(3), k4_p1(3), area_n, area_m
		real(kind = dp) :: pc, Iqs_n1, Iqs_p1
		
        complex(kind = dp) :: suma_s, sumb_s, sumc_s
        complex(kind = dp) :: Gs_1, Gs_2, DGs_1, DGs_2 
		complex(kind = dp) :: c1, c2, c3, c4, f1, f2 
		logical :: Iqs_L
		intrinsic :: sqrt
		
      call Quadrature_tri(ngp, a, b, w)
      suma = (0.0, 0.0)
      sumb = (0.0, 0.0)
      sumc = (0.0, 0.0)
      suma_s = (0.0, 0.0)
      sumb_s = (0.0, 0.0)
      sumc_s = (0.0, 0.0)
		
		k2_p1  = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		
		do tt = 1, 2 !
		do ll = 1, 2
			call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area_m)
			pos_neg_m = -1*(((tt-1)*tt)-1)            
			dfm = pos_neg_m*(Len_m/area_m)                 
        
			call fn_parameter(struct, n, ll, pn1, pn2, vn_l, Len_n, area_n)
			pos_neg_n = -1*(((ll-1)*ll)-1)                
			dfn = pos_neg_n*(Len_n/area_n)
		
			pp = struct%neighbours(m)%elements(tt)
			qq = struct%neighbours(n)%elements(ll)
					
			if ((pp == qq)) then
				Iqs_L = .true.
			else
				Iqs_L = .false.
			end if
		
			do i = 1, ngp
				Ra = r_simplex(a(i), b(i), pm1, pm2, vm_t)
				fm = (Ra-vm_t)*dfm/2
		!
				do j = 1, ngp 
					Rq = r_simplex(a(j), b(j), pn1, pn2, vn_l)
					fn = (Rq-vn_l)*dfn/2
					R = Ra - Rq			
					Rn = R/vec_len(R)
					call Diver_Green_s(k1, vec_len(R), Gs_1, DGs_1)
					call Diver_Green_s(k2, vec_len(R), Gs_2, DGs_2)				
					TRF = w(j)*area_n*w(i)*area_m
					if (Iqs_L .eqv. .true.) then
						suma = (0.0, 0.0)
					else !K_12s singular integration of K_12	
						suma = suma + dot_product(fm, cross_r(fn, Rn)*(DGs_1 + DGs_2))*TRF 												
					end if		
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					f1 = im*(Omega*my_1*dot_product(fm, fn) - 1/(Omega*eps_1)*dfm*dfn)*Gs_1 
					f2 = im*(Omega*my_2*dot_product(fm, fn) - 1/(Omega*eps_2)*dfm*dfn)*Gs_2 
					sumb = sumb + (f1 + f2)*TRF !D_11 
                        
					f1 = im*(Omega*eps_1*dot_product(fm, fn) - 1/(Omega*my_1)*dfm*dfn)*Gs_1 
					f2 = im*(Omega*eps_2*dot_product(fm, fn) - 1/(Omega*my_2)*dfm*dfn)*Gs_2 
					sumc = sumc + (f1 + f2)*TRF !
				end do ! loop j	
 
				TRF_m = w(i)*area_m			
				pc = dfm*dfn/Omega
				q = -1
				k2_n1 = K2q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
				k4_n1 = K4q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
				Iqs_n1 = Iq_s(q, Ra, pn1, pn2, vn_l, Iqs_L) !					
				q = 1 
				k2_p1 = K2q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
				k4_p1 = K4q(q, Ra, pn1, pn2, vn_l, Iqs_L)*(-1)*dfn/2!
				Iqs_p1 = Iq_s(q, Ra, pn1, pn2, vn_l, Iqs_L) !
			
				c3 = -(k1**2 + k2**2)/(8*PI)	
				if (Iqs_L .eqv. .true.) then
					suma_s = (0.0, 0.0)
				else
					suma_s = suma_s + (1/(2*PI)*dot_product(fm, k4_n1) + c3*dot_product(fm, k4_p1))*TRF_m !K_21s
				end if			
				call Coefficient_D11s(pc, c1, c2, c3, c4)	!D_11s			
				sumb_s = sumb_s + (dot_product(fm, (c1*k2_n1 + c3*k2_p1)) + (c2*Iqs_n1 + c4*Iqs_p1))*TRF_m !
			
				call Coefficient_D22s(pc, c1, c2, c3, c4)	!D_22s					
				sumc_s = sumc_s + (dot_product(fm, (c1*k2_n1 + c3*k2_p1)) + (c2*Iqs_n1 + c4*Iqs_p1))*TRF_m ! 			
			end do !loop i     
			suma = (suma + suma_s)		
			sumb = (sumb + sumb_s)
			sumc = (sumc + sumc_s)		
		end do 
		end do
    end subroutine Singular_Integration_Compact
	
	subroutine Singular_Integration(m, n, tt, ll, ngp, struct, suma, sumb, sumc)
		implicit none
		integer, intent(in) :: m, n, ngp, tt, ll
		type(Structure), intent(in) :: struct
		complex(kind = dp), intent(out) :: suma, sumb, sumc
		
		integer :: i, j, pos_neg_m, pos_neg_n, q, pp, qq
        real(kind = dp) :: TRF, TRF_m, Len_m, Len_n
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
        real(kind = dp) :: Rq(3), Rn(3), Ra(3), R(3), fm(3), fn(3), dfm, dfn, pm1(3), pm2(3), vm_t(3), vn_l(3), pn1(3), pn2(3)
		real(kind = dp) :: k2_n1(3), k2_p1(3), k4_n1(3), k4_p1(3), area_n, area_m
		real(kind = dp) :: pc, Iqs_n1, Iqs_p1
		
        complex(kind = dp) :: suma_s, sumb_s, sumc_s
        complex(kind = dp) :: Gs_1, Gs_2, DGs_1, DGs_2 
		complex(kind = dp) :: c1, c2, c3, c4, f1, f2 
		logical :: Iqs_L
		intrinsic :: sqrt
		
      call Quadrature_tri(ngp, a, b, w)
      suma = (0.0, 0.0)
      sumb = (0.0, 0.0)
      sumc = (0.0, 0.0)
      suma_s = (0.0, 0.0)
      sumb_s = (0.0, 0.0)
      sumc_s = (0.0, 0.0)
		
		k2_p1  = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
		
      call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area_m)
      pos_neg_m = -1*(((tt-1)*tt)-1)            
      dfm = pos_neg_m*(Len_m/area_m)                 
        
		call fn_parameter(struct, n, ll, pn1, pn2, vn_l, Len_n, area_n)
		pos_neg_n = -1*(((ll-1)*ll)-1)                
		dfn = pos_neg_n*(Len_n/area_n)
		
		pp = struct%neighbours(m)%elements(tt)
		qq = struct%neighbours(n)%elements(ll)
					
		if ((pp == qq)) then
			Iqs_L = .true.
		else
			Iqs_L = .false.
		end if
		
		do i = 1, ngp
			Ra = r_simplex(a(i), b(i), pm1, pm2, vm_t)
			fm = (Ra-vm_t)*dfm/2
   !
			do j = 1, ngp 
				Rq = r_simplex(a(j), b(j), pn1, pn2, vn_l)
				fn = (Rq-vn_l)*dfn/2
				R = Ra - Rq			
				Rn = R/vec_len(R)
				call Diver_Green_s(k1, vec_len(R), Gs_1, DGs_1)
				call Diver_Green_s(k2, vec_len(R), Gs_2, DGs_2)				
				TRF = w(j)*area_n*w(i)*area_m
				if (Iqs_L .eqv. .true.) then
					suma = (0.0, 0.0)
				else !K_12s singular integration of K_12	
					suma = suma + dot_product(fm, cross_r(fn, Rn)*(DGs_1 + DGs_2))*TRF 												
				end if		
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				f1 = im*(Omega*my_1*dot_product(fm, fn) - 1/(Omega*eps_1)*dfm*dfn)*Gs_1 
				f2 = im*(Omega*my_2*dot_product(fm, fn) - 1/(Omega*eps_2)*dfm*dfn)*Gs_2 
				sumb = sumb + (f1 + f2)*TRF !D_11 
                        
				f1 = im*(Omega*eps_1*dot_product(fm, fn) - 1/(Omega*my_1)*dfm*dfn)*Gs_1 
				f2 = im*(Omega*eps_2*dot_product(fm, fn) - 1/(Omega*my_2)*dfm*dfn)*Gs_2 
				sumc = sumc + (f1 + f2)*TRF !
			end do ! loop j	
 
			TRF_m = w(i)*area_m			
			pc = dfm*dfn/Omega
			q = -1

		
			!print*, 'pn1 =', pn1
			!print*, 'pn2 =', pn2			
			!print*, 'vn_l =', vn_l
			
			k2_n1 = K2q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
			k4_n1 = K4q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
			Iqs_n1 = Iq_s(q, Ra, pn1, pn2, vn_l, Iqs_L) !
					
			q = 1 
			k2_p1 = K2q(q, Ra, pn1, pn2, vn_l, Iqs_L)*dfn/2!
			k4_p1 = K4q(q, Ra, pn1, pn2, vn_l, Iqs_L)*(-1)*dfn/2!
			Iqs_p1 = Iq_s(q, Ra, pn1, pn2, vn_l, Iqs_L) !
			
			c3 = -(k1**2 + k2**2)/(8*PI)	
			if (Iqs_L .eqv. .true.) then
				suma_s = (0.0, 0.0)
			else
				suma_s = suma_s + (1/(2*PI)*dot_product(fm, k4_n1) + c3*dot_product(fm, k4_p1))*TRF_m !K_21s
			end if			
			call Coefficient_D11s(pc, c1, c2, c3, c4)	!D_11s			
			sumb_s = sumb_s + (dot_product(fm, (c1*k2_n1 + c3*k2_p1)) + (c2*Iqs_n1 + c4*Iqs_p1))*TRF_m !
			
			call Coefficient_D22s(pc, c1, c2, c3, c4)	!D_22s					
			sumc_s = sumc_s + (dot_product(fm, (c1*k2_n1 + c3*k2_p1)) + (c2*Iqs_n1 + c4*Iqs_p1))*TRF_m ! 			
		end do !loop i            
		suma = (suma + suma_s)		
		sumb = (sumb + sumb_s)
		sumc = (sumc + sumc_s)		
    end subroutine Singular_Integration  
 
	subroutine Integration_FieldFS_tri(r_a, struct, m_pair, ngp_a, SE, SH, Sum_a)
      complex(kind = dp), intent(in) ::  SE, SH        
      real(kind = dp), intent(in) :: r_a(3)
		type(Structure), intent(in) :: struct
		integer, intent(in) :: ngp_a, m_pair  
      type(Vector_c), dimension(2), intent(out) :: Sum_a
		
      integer i, t, pos_neg
      real(kind = dp) :: R_norm, Len_m, area_m, dfm 
		complex(kind = dp) :: GradG(3), K_tmp(3), Greenf, GG_SS, f_mat(3)
		real(kind = dp) :: r_q(3), p1(3), p2(3), vm_t(3)
		real(kind = dp), dimension(100) :: a, b, w	
		
      call Quadrature_tri(ngp_a, a, b, w)
      do t = 1, 2
         Sum_a(t)%Vector = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      end do  
      do t = 1, 2
			call fn_parameter(struct, m_pair, t, p1, p2, vm_t, Len_m, area_m)
			pos_neg = -1*(((t-1)*t)-1)            
			dfm = pos_neg*(Len_m)	!/area_m
			do i = 1, ngp_a
				r_q = r_simplex(a(i), b(i), p1, p2, vm_t)
				f_mat=0.5*dfm*(r_q-vm_t) !*area_m
				R_norm = vec_len(r_a - r_q)
				call Diver_Green(k1, R_norm, Greenf, GG_SS)
				GradG = (r_a - r_q)/R_norm*GG_SS
				!TRF = w(i)*area_m				
				K_tmp = im*Omega*my_1*SE*f_mat*Greenf - cross_c(SH*f_mat, GradG) - 1/(im*Omega*eps_1)*SE*dfm*GradG	
				Sum_a(1)%Vector = Sum_a(1)%Vector - K_tmp*w(i)!TRF				
				K_tmp = im*Omega*eps_1*SH*f_mat*Greenf - 1/(im*Omega*my_1)*SH*dfm*GradG + SE*cross_c(f_mat, GradG)
				Sum_a(2)%Vector = Sum_a(2)%Vector - K_tmp*w(i)!TRF
			end do            
      end do !loop i       
      return
	end subroutine Integration_FieldFS_tri	
	
	subroutine Coefficient_D11s(P, c1, c2, c3, c4)		
		real(kind = dp), intent(in) :: P
		complex(kind = dp), intent(out) :: c1, c2, c3, c4
		intrinsic :: sqrt		
		
		c1 = Im*Omega/(4*PI)*(my_1 + my_2)
		c2 = -Im*P/(4*PI)*(1/eps_1 + 1/eps_2)
		c3 = -Im*Omega/(8*PI)*(my_1*k1**2 + my_2*k2**2)
		c4 = Im*P/(8*PI)*(k1**2/eps_1 + k2**2/eps_2)
		return
	end subroutine Coefficient_D11s
		
	subroutine Coefficient_D22s(Pc, c1, c2, c3, c4)
		real(kind = dp), intent(in) :: Pc		
		complex(kind = dp), intent(out) :: c1, c2, c3, c4		
		intrinsic :: sqrt
 
		c1 = Im*Omega/(4*PI)*(eps_1 + eps_2)
		c2 = -Im*Pc/(4*PI)*(1/my_1 + 1/my_2)
		c3 = -Im*Omega/(8*PI)*(eps_1*k1**2 + eps_2*k2**2)
		c4 = Im*Pc/(8*PI)*(k1**2/my_1 + k2**2/my_2)
		return
	end subroutine Coefficient_D22s		
   
    function r_simplex(alpha, beta, r1, r2, r3)
        implicit none
        real(kind = dp), intent(in) :: r1(3), r2(3), r3(3), alpha, beta       
        real(kind = dp) :: r_simplex(3)
        r_simplex(1) = (1-alpha - beta)*r1(1) + alpha*r2(1) + beta*r3(1)
        r_simplex(2) = (1-alpha - beta)*r1(2) + alpha*r2(2) + beta*r3(2)
        r_simplex(3) = (1-alpha - beta)*r1(3) + alpha*r2(3) + beta*r3(3)
        return    
    end function 
	 
    subroutine Quadrature_tri(ngp, alpha, beta, ww)
        implicit none
        integer, intent(in) :: ngp
        real(kind = dp), dimension(100), intent(out) :: alpha, beta, ww
        integer :: m
        
		  
        do m = 1, 20
            ww(m) = 0.0
            alpha(m) = 0.0
            beta(m) = 0.0
        end do    
        
        if (ngp == 1) then
            alpha(1:ngp) = 0.33333333333333
            beta(1:ngp) = 0.33333333333333
            ww(1:ngp) = 1.000000000
        
        else if (ngp == 3) then
            alpha(1:ngp) = & 
                    (/ 0.16666666666667, &
                       0.16666666666667, &
                       0.66666666666667/)
            beta(1:ngp) = & 
                    (/ 0.16666666666667, &
                       0.66666666666667, &
                       0.16666666666667/)
            ww(1:ngp) = &
                    (/  0.33333333333333, &
                        0.33333333333333, &
                        0.33333333333333/)  
        else if (ngp == 4) then        
            alpha(1:ngp) = & 
                    (/ 0.33333333333333, &
                       0.20000000000000, &
                       0.20000000000000, &
                       0.60000000000000/)
            beta(1:ngp) = & 
                    (/ 0.33333333333333, &
                       0.20000000000000, &
                       0.60000000000000, &
                       0.20000000000000/)
            ww(1:ngp) = &
                    (/-0.56250000000000, &
                       0.52083333333333, &
                       0.52083333333333, &
                       0.52083333333333/)  
           else if (ngp == 6) then        
            alpha(1:ngp) = & 
                    (/ 0.44594849091597, &
                       0.44594849091597, &
                       0.10810301816807, &
                       0.09157621350977, &
                       0.09157621350977, &
                       0.81684757298046/)
            beta(1:ngp) = & 
                    (/ 0.44594849091597, &
                       0.10810301816807, &
                       0.44594849091597, &                       
                       0.09157621350977, &
                       0.81684757298046, &
                       0.09157621350977/)           
            ww(1:ngp) = &
                    (/  0.22338158967801, &
                        0.22338158967801, &
                        0.22338158967801, &
                        0.10995174365532, &
                        0.10995174365532, &            
                        0.10995174365532/) 
            else if (ngp == 7) then
            alpha(1:ngp) = & 
                    (/ 0.33333333333333, &
                       0.47014206410511, &
                       0.47014206410511, &
                       0.05971587178977, &
                       0.10128650732346, &
                       0.10128650732346, &
                       0.79742698535309/)
            beta(1:ngp) = & 
                    (/ 0.33333333333333, &
                       0.47014206410511, &
                       0.05971587178977, &
                       0.47014206410511, &
                       0.10128650732346, &
                       0.79742698535309, &
                       0.10128650732346/)
          
            ww(1:ngp) = &
                    (/  0.22500000000000, &
                        0.13239415278851, &
                        0.13239415278851, &
                        0.13239415278851, &
                        0.12593918054483, &            
                        0.12593918054483, &
                        0.12593918054483 /)
            
            else if (ngp == 12) then
                 alpha(1:ngp) = & 
                (/  0.249286745170910, &
                    0.249286745170910, &
                    0.501426509658180, &
                    0.063089014491500, &
                    0.063089014491500, &
                    0.873821971017000, &
                    0.310352451033780, &
                    0.636502499121400, &
                    0.053145049844820, &
                    0.636502499121400, &
                    0.310352451033780, &
                    0.053145049844820/)
                 beta(1:ngp) = & 
                  (/0.249286745170910, &
                    0.501426509658180, &
                    0.249286745170910, &
                    0.063089014491500, &
                    0.873821971017000, &
                    0.063089014491500, &
                    0.636502499121400, &
                    0.053145049844820, &
                    0.310352451033780, &
                    0.310352451033780, &
                    0.053145049844820, &
                    0.636502499121400/)
                 ww(1:ngp) = &
                 (/ 0.116786275726380, &
                    0.116786275726380, &
                    0.116786275726380, &
                    0.0508449063702100, &
                    0.0508449063702100, &
                    0.0508449063702100, &
                    0.0828510756183700, &
                    0.0828510756183700, &
                    0.0828510756183700, &
                    0.0828510756183700, &
                    0.0828510756183700, &
                    0.0828510756183700 /)
            end if
     end subroutine 
    
 end module calculation_mod_tri