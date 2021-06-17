 module FMM_tri_Spherical	
	
	use Calculation_mod_tri
	use disc_mod
	use input_tri
	use libmath
	use omp_lib
	use Gauss_Quadrature_Fast

	implicit none
	! save 
	private 
	public :: cube_definition	
	public :: test_pl_hf
	public :: FMM_Near_MVProduct_Tri	
	public :: FMM_edge_Tri
	public :: Tf_rd
	public :: Rf_rc
	public :: test_Hl2
	public :: Translation_HF
	public :: TL_km_arr
	public :: FMM_Far_MVProduct_Tri_OMP
	public :: sort_edge_index
	public :: Khat_and_TRF_legendre	
	
	
	integer, parameter :: ng = 3	
	integer, parameter :: Lm = 25
	
	
	real(kind = dp), dimension(:, :), allocatable :: k_hat
	real(kind = dp), dimension(:), allocatable :: TRF
	integer :: N_sp
	
	contains
	
	subroutine test_Hl2()
		
		!type(list_list_real) :: hl_1, hl_2 
		complex(kind = 8), dimension(:), allocatable :: hl_1
		complex(kind = 8) :: r_h, kd
		integer :: m, fnu, i, N
		real(kind = 8) :: x
		fnu = 0
		m = 40
		N = 140
		kd = (0.06 - im*4.14)*2*PI/6.0e-7
		
		open (unit=206, file = 'Hankel_2.txt', action="write",status = 'replace')			
		do i = 1, N
			x = (1.0 + (i-1)*0.1)*6e-7
			r_h = x*k2
			hl_1 =  lib_math_hankel_spherical_2(r_h, fnu, m+1)
			write (206, '(201(es19.12, tr5))') real(r_h), imag(r_h), real(hl_1(m)), imag(hl_1(m)) 			
		end do	
	   print*, 'r_h', r_h	
		print*, 'x', x
		close(206)
	end subroutine
	
	!subroutine test_Plm()
	!	implicit none	
	!	integer(kind = 4) :: Lm
	!	double precision :: phi, x
	!	type(list_list_real) :: pm, pd
	!	integer(kind = 4) :: fnu, m, n, p, i, tmp_2
	!	logical :: cs_phase 
	!	real(kind = dp), dimension(:), allocatable :: xx, ww
	!	real(kind = dp) :: a, b, suma, tmp		
	!	
	!	real (kind = dp), dimension(:), allocatable :: pl, dummy 		
	!	Lm = 5
	!	a = -1.0
	!	b = 1.0		
	!	suma = 0.0
	!	call legendre_handle(Lm, a, b, xx, ww)				
	!	fnu = 1		
	!	cs_phase = .true. 		
	!	x = 0.9		
	!	call init_list(pm, fnu, Lm + 1)			
	!	do i = 1, Lm
	!		call lib_math_associated_legendre_polynomial_range(xx(i), fnu, Lm + 1, pm, pd, cs_phase)
	!	end do		
	!			
	!	tmp = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(3, 2) !this is the (1+m)!/(l-m)!
	!	tmp = tmp*pm%item(3)%item(2)
	!	
	!	allocate(pl(Lm),dummy(Lm))		
	!	do i = 1, Lm
	!		call lib_math_legendre_polynomial(xx(i), fnu, Lm, pl, dummy)			
	!		print*, 'i = ', i
	!		print*, 'xx(i) =', xx(i)			
	!		print*, 'pl(i) =', pl(i)			
	!	end do
	!	print*, 'Analytical p1=', xx(1)
	!	print*, 'Analytical p2=', 0.5*(3*xx(2)**2-1)
	!	print*, 'Analytical p5=', 1.0/8.0*(63*xx(5)**5-70*(xx(5))**3 + 15*xx(5))	
	!end subroutine 
	
	subroutine test_pl_hf()
		real(kind = 8) :: theta
		integer :: j, ll, fnu
		real(kind = 8) :: pl, x, p2, al_j, p4, p3, p5
		real(kind = dp), dimension(:), allocatable :: pl_arr, dummy 	
		
		
		fnu = 0 !can be zero, but the index of pl is shifted. 
		
		
		!theta = 0.5
		x = cos(theta)
		!pl = 0.0
		!p2 = 0.5*(3*x**2-1.0)
		!p4 = 0.125*(35*x*x*x*x - 30.0*x*x + 3.0)
		!p3 = 0.5*(5*x*x*x-3*x)
		!p5 = 0.125*(63*x*x*x*x*x - 70*x*x*x + 15*x)
		ll = 35		
		allocate(pl_arr(ll+1),dummy(ll+1))
		
		!do j = 1, ll
		!	al_j = 4.0/(PI*sqrt((ll+2/PI)**2 - j**2))
		!	pl = pl + al_j*cos(j*theta)	
		!end do
		!pl = pl + 2.0/(PI*ll+2.0)
		
		call lib_math_legendre_polynomial(x, fnu, ll+1, pl_arr, dummy)
		print*, 'pl_0 =', pl_arr(ll)
		
	end subroutine
	!
	function Translation_HF_theta(K_m, theta, r_h) result(TL_HF_theta)
		integer, intent(in) :: K_m
		real(kind = dp), intent(in) :: theta
		real(kind = dp) :: b_j
		integer :: j
		intrinsic :: sqrt
		complex(kind = dp) :: TL_HF, TL_HF_theta
		complex(kind = dp) :: r_h
		
		TL_HF = (0.0, 0.0)
		do j = 1, K_m			
			b_j = 4.0/PI*sqrt((K_m + 1.0)**2 - j**2)
			TL_HF =TL_HF + cmplx(b_j*cos(j*theta), 0.0)
		end do 
		TL_HF = TL_HF + 2/PI*(K_m +1)	!j = 0, delta = 0.5, Eq.(10) in Song2001
		TL_HF = TL_HF*exp(-im*r_h)/(im*r_h)
		TL_HF_theta = TL_HF
	end function
	
	function Translation_HF(N3, r_ab, n, kd)result(TL_HF)
		integer, intent(in) :: N3, n
		complex(kind = dp), intent(in) :: kd
		real(kind = dp), intent(in) :: r_ab(3)
		real(kind = dp) :: r_hat(3), x
		complex(kind = dp), dimension(N3) :: TL_HF
		complex(kind = dp) :: r_h
		integer :: j
		r_hat = r_ab/vec_len(r_ab)
      r_h =  vec_len(r_ab)*kd		
		do j = 1, N3
			x = acos(dot_product(r_hat(:), k_hat(j, :)))
			TL_HF(j) = Translation_HF_theta(n, x, r_h)
		end do		 
	end function
	
	subroutine Khat_and_TRF_legendre()
		implicit none	
		real(kind = dp), dimension(:), allocatable :: xx, ww
		real(kind = dp) :: d_phi, phi, a, b
		integer :: j, k, i, m_phi
	
		m_phi = 2*Lm	
		N_sp = Lm*m_phi     
		allocate(k_hat(N_sp, 3))
		allocate(xx(Lm))
		allocate(ww(Lm))
		allocate(TRF(N_sp))
		a = -1.0
		b = 1.0
		call legendre_handle(Lm, a, b, xx, ww)
		
		d_phi = 2*PI/m_phi
		
		do j = 1, Lm
            do k = 1, m_phi
                phi = (k-1)*d_phi 
                i = k + (j-1)*m_phi
                k_hat(i, :) = (/sin(acos(xx(j)))*cos(phi), sin(acos(xx(j)))*sin(phi), xx(j)/)
				TRF(i) = ww(j)*d_phi
            end do
		end do
		deallocate(xx)
		deallocate(ww)
		print*, 'trf(12)=', TRF(12)
	end subroutine Khat_and_TRF_legendre
	
	function TL_km_arr(Nsp, r_ab, m, kd)
      integer, intent(in) :: m, Nsp    !Lm is the order limit of the polynomials, N is the sampling points on the sphere
      real(kind = 8), intent(in) :: r_ab(3)
      real(kind = 8) :: r_hat(3), x, pm_arr(Nsp, m+1)
      complex(kind = 8) :: r_h
		!real(kind = dp) :: k_hat(Nsp, 3)
      complex(kind = 8), intent(in) :: kd
      real (kind = 8), dimension(m+1) :: pm, dummy 
      complex(kind = 8) :: TL_km_arr(Nsp)
      complex(kind = 8), dimension(:), allocatable :: hl 
      integer(kind = 4) :: j, fnu, i
		
      fnu = 0 
      
      r_hat = r_ab/vec_len(r_ab)
      r_h =  vec_len(r_ab)*kd
      print*, 'r_h =', r_h
		print*, 'm=', m
		allocate(hl(m+1))
      hl =  lib_math_hankel_spherical_2(r_h, fnu, m+1) 
		
		!open (unit=206, file = 'Hankel_2_TL.txt', action="write",status = 'replace')	
		!!	write (206, *) 'fnu=1'
		!do i = 1, m+1
		!	!print*, 'm=', i
		!	!print*, 'hl', hl(i)
		!
		!	write (206, '(201(es19.12, tr5))') real(hl(i)), imag(hl(i)) 
		!end do
		!close(206)
		
      TL_km_arr(1:Nsp) = (0.0, 0.0)     
		
      do j = 1, Nsp
         x = dot_product(r_hat(:), k_hat(j, :))
			!print*, 'j=', j
			!print*, 'x=', abs(x)
         call lib_math_legendre_polynomial(x, fnu, m+1, pm, dummy)
         pm_arr(j, :) = pm(:)
      end do
      
      do j = 1, m+1
			!print*, 'real(hl(j)=', real(hl(j))
         TL_km_arr(:) = TL_km_arr(:) + (real(hl(j))*pm_arr(:, j)*kd/(4*PI)**2)*(-im)**j*(2*(j-1) + 1)! 
    end do
	 
      return
	end function TL_km_arr
	
	function Rf_rc(r_oa, kd) 
		implicit none
		integer :: j
		complex(kind = 8), dimension(N_sp) :: Rf_rc
		real(kind = dp) :: r_oa(3)
		complex(kind = dp) :: kd

		!do j = 1, N_sp
		!	Rf_rc(j) = exp(-im*kd*dot_product(k_hat(j, :), r_oa(:)))*TRF(j)
		!end do
		Rf_rc(:) = TRF(:)*exp(-im*kd*(k_hat(:, 1)*r_oa(1) + k_hat(:, 2)*r_oa(2) + k_hat(:, 3)*r_oa(3)))
		!print*, 'Rf_rc(i)', Rf_rc(15)
	end function	
	
	function Tf_rd(r_pb, kd) 
		implicit none
		integer :: j
		complex(kind = 8), dimension(N_sp) :: Tf_rd
		real(kind = dp) :: r_pb(3)
		complex(kind = dp) :: kd
		
		Tf_rd(:) = exp(im*kd*(k_hat(:, 1)*r_pb(1)+ k_hat(:, 2)*r_pb(2) + k_hat(:, 3)*r_pb(3)))		
	end function	
	
	!Implemented according to Gibson's book "The Method of Moments in Electromagnetics"	
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
	
	subroutine cube_definition(struct, struct_cube_II, d_box) !, cube_list_level_I
		!First test for spheres
		implicit none
		type(Structure_cube) :: struct_cube_I, struct_cube_tmp
		type(Structure_cube), intent(out) :: struct_cube_II 
		type(Structure), intent(in) :: struct
		integer :: n_pairs, m, n_tmp
		integer :: cube_list_level_I, ot, N, tmp_far(1200), tmp_near(1200), tmp_arr(1200)
		integer :: i, j, k, Nx, Ny, Nz, n_counter, m_counter, p_counter
		real(kind = dp) :: p_edgec(3), r_in, r_out, W_out, r(3), W_box, p2(3)
		real(kind = dp), intent(out) :: d_box
		real(kind = dp) :: x_min, x_max, y_min, y_max, z_min, z_max, dc
		intrinsic sqrt
		logical :: con1, con2, con3
		
		!call Data_Import(struct)
		n_pairs = size(struct%neighbours)	
		
		
		W_out = D_factor*1.1 !D_factor: the radius of the phere 			
		!print*, 'n_pairs =', n_pairs
		Nx = 8 !Temporary solution
		N = Nx*Nx*Nx!ceiling(W_out/W_box) !only two level		
		print*, 'Number of total cubes', N
		n_tmp = 0
		W_box = W_out/Nx
		allocate(Struct_cube_I%cubes(N))
		allocate(Struct_cube_tmp%cubes(N))
		allocate(Struct_cube_II%pair_cubes(n_pairs))
		do i = 1, Nx			
			do j = 1, Nx
				do k = 1, Nx
				    m = ((i-1)*Nx + (j-1))*Nx + k					
					Struct_cube_I%cubes(m)%cube_position = &
					(/-W_out/2 + (i - 0.5)*W_box, -W_out/2 + (j - 0.5)*W_box, -W_out/2 + (k - 0.5)*W_box/) !					
				end do
			end do
		end do
		!		
		ot = 1		
		d_box = sqrt(3.0)*W_box
		print*, 'd_box =', d_box
		n_counter = 0
		m_counter = 1
		
		do k = 1, N
			p2 = struct_cube_I%cubes(k)%cube_position
			x_min = p2(1) - W_box/2
			x_max = p2(1) + W_box/2
			y_min = p2(2) - W_box/2
			y_max = p2(2) + W_box/2
			z_min = p2(3) - W_box/2
			z_max = p2(3) + W_box/2
			do m = 1, n_pairs
				call fn_edge_center(struct, m, ot, p_edgec)	!find the edge center of an element pair
				con1 = (x_min < p_edgec(1)) .and. (p_edgec(1) <= x_max)
				con2 = (y_min < p_edgec(2)) .and. (p_edgec(2) <= y_max)
				con3 = (z_min < p_edgec(3)) .and. (p_edgec(3) <= z_max)
				
				if ( con1 .and. con2 .and. con3) then
					n_counter = n_counter + 1 !number of edges in one cube						
					tmp_arr(n_counter) = m !edge index saved	
					Struct_cube_II%pair_cubes(m)%cube_index = m_counter !edge 				
				end if				
			end do
			if (n_counter .ne. 0) then				
				struct_cube_tmp%cubes(m_counter)%cube_position = p2					
				allocate(struct_cube_tmp%cubes(m_counter)%edges_in_cube(n_counter))
				struct_cube_tmp%cubes(m_counter)%edges_in_cube(1:n_counter) = tmp_arr(1:n_counter)				
				m_counter = m_counter + 1 !index of cubes				
			end if
			n_counter = 0
		end do
		!Copy the data in tmp to struct_cube_II
		m_counter = m_counter - 1
		allocate(struct_cube_II%cubes(m_counter))		
		do i = 1, m_counter
			!print*, 'box number =', i
		    m = size(struct_cube_tmp%cubes(i)%edges_in_cube)
			 !print*, 'edges in the box =', m
			allocate(struct_cube_II%cubes(i)%edges_in_cube(m))
			struct_cube_II%cubes(i)%edges_in_cube(1:m) = struct_cube_tmp%cubes(i)%edges_in_cube(1:m)
			struct_cube_II%cubes(i)%cube_position = struct_cube_tmp%cubes(i)%cube_position
			n_tmp = n_tmp + m
		end do
		print*, 'n_tmp =', n_tmp
		
		!-------finde the cubes in the near regions and save the cube indices to cube_near
		n_counter = 0
		p_counter = 0
		do i = 1, m_counter
			do j = 1, m_counter
				if (vec_len(struct_cube_tmp%cubes(i)%cube_position - struct_cube_tmp%cubes(j)%cube_position) <= 1.5*d_box) then				
					n_counter = n_counter + 1
					tmp_near(n_counter) = j					
				else
					p_counter = p_counter + 1
					tmp_far(p_counter) = j
				end if					
			end do
			allocate(struct_cube_II%cubes(i)%cubes_near(n_counter))
			allocate(struct_cube_II%cubes(i)%cubes_far(p_counter))
			struct_cube_II%cubes(i)%cubes_near(1:n_counter) = tmp_near(1:n_counter)
			struct_cube_II%cubes(i)%cubes_far(1:p_counter) = tmp_far(1:p_counter)
			n_counter = 0
			p_counter = 0
		end do	
	end subroutine cube_definition
	
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
		integer :: p, q, ngp, t, s, n, j, n_csr, m_csr, counter 
		
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
				counter = counter +1
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
			!if (m==m_pairs) then
			!print*, 'Edges in the near field counter =', counter
			!end if
			!counter = 0
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
				call FMM_cube_Tri_KL(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, sum_r)
				rr(m) = rr(m) +  sum_r(1)
				rr(m + m_pairs) = rr(m + m_pairs) +  sum_r(2)
			end do
		end do
		!$omp end do
		!$omp end parallel			
		call timestamp()
		return		
	end subroutine FMM_Far_MVProduct_Tri_OMP
	
	subroutine FMM_cube_Tri_KL(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, rr) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n_cube, m_pairs
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp) :: sum_a, sum_c, sum_b
		character(len = 1), intent(in) :: transa
		
		complex(kind = dp) :: sum_r(2), sum_t(2)
		complex(kind = dp),  intent(out) :: rr(2)
		complex(kind = dp) :: ff1, gg1, ff2, gg2, ss1, ss2
		
		complex(kind = dp), dimension(:, :, :), allocatable :: Rfla, Tfla
		complex(kind = dp), dimension(:, :), allocatable :: Rflb, Tflb
		complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2, x1
		complex(kind = dp), dimension(:, :, :), allocatable :: Rfk!, Tfk 
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3) !k_hat(ngp*n_phi, 3), TRF(ngp*n_phi),
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
		
		integer :: j, pos_neg, Lm, m_cube, tt, ll, nn_edge, nd, n_edge_number
		intrinsic :: sqrt
		
		!N_sp = n_theta*n_phi
		
		allocate(Rfla(N_sp, 2, 3), Tfla(N_sp, 2, 3)) 
		allocate(TL_1(N_sp), TL_2(N_sp), Rflb(N_sp, 2), Tflb(N_sp, 2))		
		allocate(Rfk(N_sp, 2, 3)) 
		
		sum_r(1) = (0.0, 0.0)
		sum_r(2) = (0.0, 0.0)
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located			
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(N_sp, Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(N_sp, Rab, Lm, k2)
		
		n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube		
		do nn_edge = 1, n_edge_number		
			nd = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n				
			sum_a = (0.0, 0.0)
			sum_b = (0.0, 0.0)
			sum_c = (0.0, 0.0)
			do tt = 1, 2
				call Receivef_rc_KL(struct, m, tt, Ra, Rfk, Rfla, Rflb)
				do ll = 1, 2
					call Radiationf_rd_KL(struct, nd, ll, Rb, Tfla, Tflb)	!Tfk = Tf1a
					do j = 1, N_sp
					   ss1 = dot_product(conjg(Rfk(j, 1, :)), Tfla(j, 1, :))
						ss2 = dot_product(conjg(Rfk(j, 2, :)), Tfla(j, 2, :))
						
						ff1 = omega*dot_product(conjg(Rfla(j, 1, :)), Tfla(j, 1, :))
						gg1 = Rflb(j, 1)*Tflb(j, 1)/omega
						ff2 = omega*dot_product(conjg(Rfla(j, 2, :)), Tfla(j, 2, :))
						gg2 = Rflb(j, 2)*Tflb(j, 2)/omega
						
						sum_b = sum_b + (TL_1(j)*(my_1*ff1 - gg1/eps_1) + TL_2(j)*(my_2*ff2 - gg2/eps_2))*TRF(j)*im
						sum_c = sum_c + (TL_1(j)*(eps_1*ff1 - gg1/my_1) + TL_2(j)*(eps_2*ff2 - gg2/my_2))*TRF(j)*im		
						sum_a = sum_a + (TL_1(j)*ss1 + TL_2(j)*ss2)*TRF(j)
					end do				
				end do				
			end do
			sum_r(1) = sum_r(1) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs)
			sum_r(2) = sum_r(2) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs) ! for the item with index m_pairs + m
			
			sum_t(1) = sum_t(1) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs)
			sum_t(2) = sum_t(2) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
		end do		
		if (transa == 'T') then
			rr = sum_t
		else if (transa == 'N') then
			rr = sum_r
		else 
			print*, 'Not a proper multiplication type!'
		end if
	end subroutine FMM_cube_Tri_KL
	
	subroutine FMM_edge_Tri(m, n, struct, struct_cube, d_c, sum_a, sum_b, sum_c) !	
	!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
		implicit none		
		integer, intent(in) :: m, n
		type(Structure), intent(in) :: struct
		type(Structure_cube), intent(in) :: struct_cube
		complex(kind = dp), intent(out) :: sum_a, sum_c, sum_b
		complex(kind = dp) :: tmp_1, tmp_2
		complex(kind = dp), dimension(:, :, :), allocatable :: Rfla, Tfla 
		complex(kind = dp), dimension(:, :), allocatable :: Rflb, Tflb	
		complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2
		complex(kind = dp), dimension(:, :, :), allocatable :: Rfk
		!complex(kind = dp), dimension(:), allocatable :: ff1, gg1, ff2, gg2, ss1, ss2
		complex(kind = dp) :: ff1, gg1, ff2, gg2, ss1, ss2
		real(kind = dp) ::  Ra(3), Rb(3), Rab(3) 
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
		
		!N_sp = n_theta*n_phi
		
		allocate(Rfla(N_sp, 2, 3), Tfla(N_sp, 2, 3)) 
		allocate(TL_1(N_sp), TL_2(N_sp), Rflb(N_sp, 2), Tflb(N_sp, 2))
		allocate(Rfk(N_sp, 2, 3))
		
		!allocate(ff1(N_sp), gg1(N_sp), ff2(N_sp), gg2(N_sp), ss1(N_sp), ss2(N_sp))
		
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		
		m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
		n_cube = struct_cube%pair_cubes(n)%cube_index
		Ra = struct_cube%cubes(m_cube)%cube_position
		Rb = struct_cube%cubes(n_cube)%cube_position	
		Rab = Ra - Rb
		
		Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
		TL_1 = TL_km_arr(N_sp, Rab, Lm, k1)
		
		Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
		TL_2 = TL_km_arr(N_sp, Rab, Lm, k2)
		
		sum_a = (0.0, 0.0)
		sum_b = (0.0, 0.0)
		sum_c = (0.0, 0.0)
		do tt = 1, 2
			call Receivef_rc_KL(struct, m, tt, Ra, Rfk, Rfla, Rflb)
			do ll = 1, 2				
				call Radiationf_rd_KL(struct, n, ll, Rb, Tfla, Tflb)	!Tfk = Tf1a				
				!ss1(:)= Rfk(:, 1, 1)*Tfla(:, 1, 1) + Rfk(:, 1, 2)*Tfla(:, 1, 2) + Rfk(:, 1, 3)*Tfla(:, 1, 3) 
				!ss2(:)= Rfk(:, 2, 1)*Tfla(:, 2, 1) + Rfk(:, 2, 2)*Tfla(:, 2, 2) + Rfk(:, 2, 3)*Tfla(:, 2, 3) 
				!ff1(:) = omega*(Rfla(:, 1, 1)*Tfla(:, 1, 1) + Rfla(:, 1, 2)*Tfla(:, 1, 2) + Rfla(:, 1, 3)*Tfla(:, 1, 3))
				!ff2(:) = omega*(Rfla(:, 2, 1)*Tfla(:, 2, 1) + Rfla(:, 2, 2)*Tfla(:, 2, 2) + Rfla(:, 2, 3)*Tfla(:, 2, 3))
				!gg1(:) = Rflb(:, 1)*Tflb(:, 1)/omega
				!gg2(:) = Rflb(:, 2)*Tflb(:, 2)/omega
				!sum_a = sum_a + sum((TL_1(:)*ss1(:) + TL_2(:)*ss2(:))*TRF(:))
				!sum_b = sum_b + sum((TL_1(:)*(my_1*ff1(:) - gg1(:)/eps_1) + TL_2(:)*(my_2*ff2(:) - gg2(:)/eps_2))*TRF(:))*im
				!sum_c = sum_c + sum((TL_1(:)*(eps_1*ff1(:) - gg1(:)/my_1) + TL_2(:)*(eps_2*ff2(:) - gg2(:)/my_2))*TRF(:))*im	
				
				do j = 1, N_sp
					ss1 = dot_product(conjg(Rfk(j, 1, :)), Tfla(j, 1, :))
					ss2 = dot_product(conjg(Rfk(j, 2, :)), Tfla(j, 2, :))
						
					ff1 = omega*dot_product(conjg(Rfla(j, 1, :)), Tfla(j, 1, :))
					gg1 = Rflb(j, 1)*Tflb(j, 1)/omega
					ff2 = omega*dot_product(conjg(Rfla(j, 2, :)), Tfla(j, 2, :))
					gg2 = Rflb(j, 2)*Tflb(j, 2)/omega
						
					sum_b = sum_b + (TL_1(j)*(my_1*ff1 - gg1/eps_1) + TL_2(j)*(my_2*ff2 - gg2/eps_2))*TRF(j)*im
					sum_c = sum_c + (TL_1(j)*(eps_1*ff1 - gg1/my_1) + TL_2(j)*(eps_2*ff2 - gg2/my_2))*TRF(j)*im		
					sum_a = sum_a + (TL_1(j)*ss1 + TL_2(j)*ss2)*TRF(j)
				end do				
			end do				
		end do		
    end subroutine FMM_edge_Tri
	
		
	subroutine Receivef_rc_KL(struct, m, tt, R_a, Rf_rc_k, Rf_rc_La, Rf_rc_Lb)
	! rc: receive
	! rd: radiation
	! tr: transfer function
	! R_oa :radius around the observation point
	! R_pb: radius around the source point
	! R_ab : distance between the two points
		implicit none
		integer, intent(in) :: m, tt !sampling point on a sphere				
		real(kind = dp), intent(in) :: R_a(3) 
		
		real(kind = dp), dimension(100) :: w(100), a(100), b(100)
		real(kind = 8) :: pm1(3), pm2(3), vm_t(3), R_oa(3), R_o(3), fm(3)
		complex(kind = dp), dimension(N_sp, 2, 3), intent(out) :: Rf_rc_La, Rf_rc_k
		complex(kind = dp), dimension(N_sp, 2), intent(out) :: Rf_rc_Lb 
		complex(kind = 8) :: ff1, ff2, tmp, tmp_r(3)
		!complex(kind = 8), dimension(N_sp) :: ff1, ff2
		real(kind = 8) :: area, Len_m, dfm
		integer :: i, j, n, pos_neg
		
		type(Structure) :: struct
		
		call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
		pos_neg = -1*(((tt-1)*tt)-1)            
		dfm = pos_neg*Len_m!				
		call Quadrature_tri(ng, a, b, w)
		do i = 1, N_sp
			Rf_rc_K(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Rf_rc_La(i, 1,:) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Rf_rc_K(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Rf_rc_La(i, 2,:) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			
			Rf_rc_Lb(i, 1:2) = (0.0, 0.0)
			
			do j = 1, ng
				R_o = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
				fm = 0.5*(R_o-vm_t)
				R_oa = R_o - R_a
				tmp = dot_product(k_hat(i, :), R_oa(:))
				ff1 = dfm*exp(-im*k1*tmp)*w(j)!
				ff2 = dfm*exp(-im*k2*tmp)*w(j)!
				
				Rf_rc_La(i, 1, :) = Rf_rc_La(i, 1, :) + fm*ff1
				Rf_rc_Lb(i, 1) = Rf_rc_Lb(i, 1) + ff1
				Rf_rc_La(i, 2, :) = Rf_rc_La(i, 2, :) + fm*ff2
        		Rf_rc_Lb(i, 2) = Rf_rc_Lb(i, 2) + ff2		
				
				tmp_r = im*cross_r(k_hat(i, 1:3), fm)
				Rf_rc_k(i, 1, :) = Rf_rc_k(i, 1, :) - k1*tmp_r*ff1!
				Rf_rc_k(i, 2, :) = Rf_rc_k(i, 2, :) - k2*tmp_r*ff2!
			end do	
		end do	
		return
	end subroutine Receivef_rc_KL
	
	Subroutine Radiationf_rd_KL(struct, m, tt, R_b, Tf_rd_La, Tf_rd_Lb)	
	! rc: receive
	! rd: radiation
	! tr: transfer function
	! R_oa :radius around the observation point
	! R_pb: radius around the source point
	! R_ab : distance between the two points
		implicit none
		integer, intent(in) :: m, tt !sampling point on a sphere				
		real(kind = dp), intent(in) :: R_b(3) 	
		
		real(kind = dp) :: a(100), b(100), w(100)
		real(kind = 8) :: pm1(3), pm2(3), vm_t(3), R_pb(3), R_p(3), fm(3)
		complex(kind = dp), dimension(N_sp, 2, 3), intent(out) :: Tf_rd_La!, Tf_rd_K
		complex(kind = dp), dimension(N_sp, 2), intent(out) :: Tf_rd_Lb  
		real(kind = dp) :: area, Len_m, dfm
		integer :: i, j, n, pos_neg
		
		complex(kind = dp) :: ff1, ff2
		
		type(Structure) :: struct
		
		call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
		pos_neg = -1*(((tt-1)*tt)-1)            
		dfm = pos_neg*Len_m!				
		call Quadrature_tri(ng, a, b, w)
		do i = 1, N_sp
		   Tf_rd_La(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			Tf_rd_Lb(i, 1:2) = (0.0, 0.0)
		   Tf_rd_La(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
			!Attention: 
			!Tr_rd_K = Tf_rd_La
			do j = 1, ng
				R_p = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinates
				fm = 0.5*(R_p-vm_t)
				R_pb = R_p - R_b
				ff1 = dfm*exp(im*k1*dot_product(k_hat(i, :), R_pb(:)))*w(j)
				ff2 = dfm*exp(im*k2*dot_product(k_hat(i, :), R_pb(:)))*w(j)
				Tf_rd_La(i, 1, :) = Tf_rd_La(i, 1, :) + fm*ff1
				Tf_rd_Lb(i, 1) = Tf_rd_Lb(i, 1) + ff1
				Tf_rd_La(i, 2, :) = Tf_rd_La(i, 2, :) + fm*ff2
				Tf_rd_Lb(i, 2) = Tf_rd_Lb(i, 2) + ff2				
			end do	
		end do
	end subroutine Radiationf_rd_KL	
	
	!function TL_km_arr(r_ab, Lm, N_sp, km)
 !
 !     integer, intent(in) :: Lm, N_sp     !Lm is the order limit of the polynomials, N is the sampling points on the sphere
 !     real(kind = dp), intent(in) :: r_ab(3)
 !     real(kind = dp) :: r_hat(3), x, pm_arr(N_sp, Lm)
 !     complex(kind = dp) :: r_h
 !     complex(kind = dp), intent(in) :: km
 !     real (kind = dp), dimension(Lm) :: pm, dummy 
 !     complex(kind = dp) :: TL_km_arr(N_sp)
 !     complex(kind = dp), dimension(Lm) :: hl 
 !     integer :: j, fnu
	!	
	!	
 !     fnu = 0        
 !     r_hat = r_ab/vec_len(r_ab)
 !     r_h =  vec_len(r_ab)*km  
 !     hl =  lib_math_hankel_spherical_2(r_h, fnu, Lm) 
 !     TL_km_arr(1:N_sp) = (0.0, 0.0)         
 !     do j = 1, N_sp
 !        x = dot_product(r_hat, k_hat(j, :))
 !        call lib_math_legendre_polynomial(x, fnu, Lm, pm, dummy)
 !        pm_arr(j, :) = pm(:)
 !     end do
	!	
 !     do j = 1, Lm
 !        TL_km_arr(:) = TL_km_arr(:) + hl(j)*(-im)**j*(2*(j-1) + 1)*pm_arr(:, j)*km/(4*PI)**2 ! 
	!	end do
 !     return        
	!end function TL_km_arr
	!
	!Subroutine Radiationf_rd_KL(struct, m, tt, R_b, Tf_rd_La, Tf_rd_Lb)	
	!! rc: receive
	!! rd: radiation
	!! tr: transfer function
	!! R_oa :radius around the observation point
	!! R_pb: radius around the source point
	!! R_ab : distance between the two points
	!	implicit none
	!	integer, intent(in) :: m, tt !sampling point on a sphere				
	!	real(kind = dp), intent(in) :: R_b(3) 	
	!	
	!	real(kind = dp) :: a(100), b(100), w(100)
	!	real(kind = 8) :: pm1(3), pm2(3), vm_t(3), R_pb(3), R_p(3), fm(3)
	!	complex(kind = dp), dimension(N_sp, 2, 3), intent(out) :: Tf_rd_La!, Tf_rd_K
	!	complex(kind = dp), dimension(N_sp, 2), intent(out) :: Tf_rd_Lb  
	!	real(kind = dp) :: area, Len_m, dfm
	!	integer :: i, j, n, pos_neg
	!	
	!	complex(kind = dp) :: ff1, ff2
	!	
	!	type(Structure) :: struct
	!	
	!	call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
	!	pos_neg = -1*(((tt-1)*tt)-1)            
	!	dfm = pos_neg*Len_m!				
	!	call Quadrature_tri(ng, a, b, w)
	!	do i = 1, N_sp
	!	   Tf_rd_La(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		Tf_rd_Lb(i, 1:2) = (0.0, 0.0)
	!	   Tf_rd_La(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		!Attention: 
	!		!Tr_rd_K = Tf_rd_La
	!		do j = 1, ng
	!			R_p = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinates
	!			fm = 0.5*(R_p-vm_t)
	!			R_pb = R_p - R_b
	!			ff1 = dfm*exp(im*k1*dot_product(k_hat(i, :), R_pb(:)))*w(j)
	!			ff2 = dfm*exp(im*k2*dot_product(k_hat(i, :), R_pb(:)))*w(j)
	!			Tf_rd_La(i, 1, :) = Tf_rd_La(i, 1, :) + fm*ff1
	!			Tf_rd_Lb(i, 1) = Tf_rd_Lb(i, 1) + ff1
	!			Tf_rd_La(i, 2, :) = Tf_rd_La(i, 2, :) + fm*ff2
	!			Tf_rd_Lb(i, 2) = Tf_rd_Lb(i, 2) + ff2				
	!		end do	
	!	end do
	!end subroutine Radiationf_rd_KL
	!
	!subroutine Receivef_rc_KL(struct, m, tt, R_a, Rf_rc_k, Rf_rc_La, Rf_rc_Lb)
	!! rc: receive
	!! rd: radiation
	!! tr: transfer function
	!! R_oa :radius around the observation point
	!! R_pb: radius around the source point
	!! R_ab : distance between the two points
	!	implicit none
	!	integer, intent(in) :: m, tt !sampling point on a sphere				
	!	real(kind = dp), intent(in) :: R_a(3) 
	!	
	!	real(kind = dp), dimension(100) :: w(100), a(100), b(100)
	!	real(kind = 8) :: pm1(3), pm2(3), vm_t(3), R_oa(3), R_o(3), fm(3)
	!	complex(kind = dp), dimension(N_sp, 2, 3), intent(out) :: Rf_rc_La, Rf_rc_k
	!	complex(kind = dp), dimension(N_sp, 2), intent(out) :: Rf_rc_Lb 
	!	complex(kind = 8) :: ff1, ff2, tmp, tmp_r(3)
	!	!complex(kind = 8), dimension(N_sp) :: ff1, ff2
	!	real(kind = 8) :: area, Len_m, dfm
	!	integer :: i, j, n, pos_neg
	!	
	!	type(Structure) :: struct
	!	
	!	call fn_parameter(struct, m, tt, pm1, pm2, vm_t, Len_m, area)
	!	pos_neg = -1*(((tt-1)*tt)-1)            
	!	dfm = pos_neg*Len_m!				
	!	call Quadrature_tri(ng, a, b, w)
	!	do i = 1, N_sp
	!		Rf_rc_K(i, 1, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		Rf_rc_La(i, 1,:) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		Rf_rc_K(i, 2, :) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		Rf_rc_La(i, 2,:) = (/(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
	!		
	!		Rf_rc_Lb(i, 1:2) = (0.0, 0.0)
	!		
	!		do j = 1, ng
	!			R_o = r_simplex(a(j), b(j), pm1, pm2, vm_t) !point r in the original coordinate
	!			fm = 0.5*(R_o-vm_t)
	!			R_oa = R_o - R_a
	!			tmp = dot_product(k_hat(i, :), R_oa(:))
	!			ff1 = dfm*exp(-im*k1*tmp)*w(j)!
	!			ff2 = dfm*exp(-im*k2*tmp)*w(j)!
	!			
	!			Rf_rc_La(i, 1, :) = Rf_rc_La(i, 1, :) + fm*ff1
	!			Rf_rc_Lb(i, 1) = Rf_rc_Lb(i, 1) + ff1
	!			Rf_rc_La(i, 2, :) = Rf_rc_La(i, 2, :) + fm*ff2
 !       		Rf_rc_Lb(i, 2) = Rf_rc_Lb(i, 2) + ff2		
	!			
	!			tmp_r = im*cross_r(k_hat(i, 1:3), fm)
	!			Rf_rc_k(i, 1, :) = Rf_rc_k(i, 1, :) - k1*tmp_r*ff1!
	!			Rf_rc_k(i, 2, :) = Rf_rc_k(i, 2, :) - k2*tmp_r*ff2!
	!		end do	
	!	end do	
	!	return
	!end subroutine Receivef_rc_KL
	
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
	!end subroutine FMM_Far_MVProduct_Tri_OMP
	
	!subroutine FMM_cube_Tri_KL(transa, m, m_pairs, n_cube, struct, struct_cube, d_c, x1, rr) !	
	!!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
	!	implicit none		
	!	integer, intent(in) :: m, n_cube, m_pairs
	!	type(Structure), intent(in) :: struct
	!	type(Structure_cube), intent(in) :: struct_cube
	!	complex(kind = dp) :: sum_a, sum_c, sum_b
	!	character(len = 1), intent(in) :: transa
	!	
	!	complex(kind = dp) :: sum_r(2), sum_t(2)
	!	complex(kind = dp),  intent(out) :: rr(2)
	!	complex(kind = dp) :: ff1, gg1, ff2, gg2, ss1, ss2
	!	
	!	complex(kind = dp), dimension(:, :, :), allocatable :: Rfla, Tfla
	!	complex(kind = dp), dimension(:, :), allocatable :: Rflb, Tflb
	!	complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2, x1
	!	complex(kind = dp), dimension(:, :, :), allocatable :: Rfk!, Tfk 
	!	real(kind = dp) ::  Ra(3), Rb(3), Rab(3) !k_hat(ngp*n_phi, 3), TRF(ngp*n_phi),
	!	real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
	!	!-------------for test use-----------------
	!	!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
	!	!       Rq
	!	!       o
	!	!       | 
	!	!	     |
	!	!	  Ra o-----------------------o Rb
	!	!                                |
	!	!                                |
	!	!                                o Rp
	!	!------------------------------------------
	!	
	!	integer :: j, pos_neg, Lm, m_cube, tt, ll, nn_edge, nd, n_edge_number
	!	intrinsic :: sqrt
	!	
	!	!N_sp = n_theta*n_phi
	!	
	!	allocate(Rfla(N_sp, 2, 3), Tfla(N_sp, 2, 3)) 
	!	allocate(TL_1(N_sp), TL_2(N_sp), Rflb(N_sp, 2), Tflb(N_sp, 2))		
	!	allocate(Rfk(N_sp, 2, 3)) 
	!	
	!	sum_r(1) = (0.0, 0.0)
	!	sum_r(2) = (0.0, 0.0)
	!	
	!	m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located			
	!	Ra = struct_cube%cubes(m_cube)%cube_position
	!	Rb = struct_cube%cubes(n_cube)%cube_position	
	!	Rab = Ra - Rb
	!	
	!	Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
	!	TL_1 = TL_km_arr(Rab, Lm, k1)
	!	
	!	Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
	!	TL_2 = TL_km_arr(Rab, Lm, k2)
	!	
	!	n_edge_number = size(struct_cube%cubes(n_cube)%edges_in_cube) !number of edges in cube n_cube		
	!	do nn_edge = 1, n_edge_number		
	!		nd = struct_cube%cubes(n_cube)%edges_in_cube(nn_edge) !edge indices in all the near cubes of cube n				
	!		sum_a = (0.0, 0.0)
	!		sum_b = (0.0, 0.0)
	!		sum_c = (0.0, 0.0)
	!		do tt = 1, 2
	!			call Receivef_rc_KL(struct, m, tt, Ra, Rfk, Rfla, Rflb)
	!			do ll = 1, 2				
	!				call Radiationf_rd_KL(struct, nd, ll, Rb, Tfla, Tflb)	!Tfk = Tf1a
	!				do j = 1, N_sp
	!				   ss1 = dot_product(conjg(Rfk(j, 1, :)), Tfla(j, 1, :))
	!					ss2 = dot_product(conjg(Rfk(j, 2, :)), Tfla(j, 2, :))
	!					
	!					ff1 = omega*dot_product(conjg(Rfla(j, 1, :)), Tfla(j, 1, :))
	!					gg1 = Rflb(j, 1)*Tflb(j, 1)/omega
	!					ff2 = omega*dot_product(conjg(Rfla(j, 2, :)), Tfla(j, 2, :))
	!					gg2 = Rflb(j, 2)*Tflb(j, 2)/omega
	!					
	!					sum_b = sum_b + (TL_1(j)*(my_1*ff1 - gg1/eps_1) + TL_2(j)*(my_2*ff2 - gg2/eps_2))*TRF(j)*im
	!					sum_c = sum_c + (TL_1(j)*(eps_1*ff1 - gg1/my_1) + TL_2(j)*(eps_2*ff2 - gg2/my_2))*TRF(j)*im		
	!					sum_a = sum_a + (TL_1(j)*ss1 + TL_2(j)*ss2)*TRF(j)
	!				end do				
	!			end do				
	!		end do
	!		sum_r(1) = sum_r(1) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs)
	!		sum_r(2) = sum_r(2) + sum_b*x1(nd) - sum_a*x1(nd + m_pairs) ! for the item with index m_pairs + m
	!		
	!		sum_t(1) = sum_t(1) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs)
	!		sum_t(2) = sum_t(2) + sum_a*x1(nd) + sum_c*x1(nd + m_pairs) ! for the item with index m_pairs + m			
	!	end do		
	!	if (transa == 'T') then
	!		rr = sum_t
	!	else if (transa == 'N') then
	!		rr = sum_r
	!	else 
	!		print*, 'Not a proper multiplication type!'
	!	end if
	!end subroutine FMM_cube_Tri_KL
	!
	!subroutine FMM_edge_Tri(m, n, struct, struct_cube, d_c, sum_a, sum_b, sum_c) !	
	!!  | Z_EJ( = L1 + L2) Z_EM (=-(K1 + K2)) | |J|  
	!!  | Z_HJ( = K1 + K2) Z_HM( = eps_1/mu_1*L1 + eps_2/mu_2*L2) | |M|
	!	implicit none		
	!	integer, intent(in) :: m, n
	!	type(Structure), intent(in) :: struct
	!	type(Structure_cube), intent(in) :: struct_cube
	!	complex(kind = dp), intent(out) :: sum_a, sum_c, sum_b
	!	complex(kind = dp) :: tmp_1, tmp_2
	!	complex(kind = dp), dimension(:, :, :), allocatable :: Rfla, Tfla 
	!	complex(kind = dp), dimension(:, :), allocatable :: Rflb, Tflb	
	!	complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2
	!	complex(kind = dp), dimension(:, :, :), allocatable :: Rfk
	!	complex(kind = dp) :: ff1, gg1, ff2, gg2, ss1, ss2
	!	real(kind = dp) ::  Ra(3), Rb(3), Rab(3) !k_hat(ngp*n_phi, 3), TRF(ngp*n_phi),
	!	real(kind = dp), intent(in) :: d_c ! sphere diameter of the first level		
	!	!real(kind = dp) :: dc		
	!	!-------------for test use-----------------
	!	!real(kind = 8) :: Ra(3), Rb(3), Rab(3)
	!	!       Rq
	!	!       o
	!	!       | 
	!	!	     |
	!	!	  Ra o-----------------------o Rb
	!	!                                |
	!	!                                |
	!	!                                o Rp
	!	!------------------------------------------
	!	
	!	integer :: j, pos_neg, N_sp, Lm, m_cube, n_cube, tt, ll
	!	intrinsic :: sqrt
	!	
	!	
	!	allocate(Rfla(N_sp, 2, 3), Tfla(N_sp, 2, 3)) 
	!	allocate(TL_1(N_sp), TL_2(N_sp), Rflb(N_sp, 2), Tflb(N_sp, 2))
	!	allocate(Rfk(N_sp, 2, 3))
	!	!allocate(xm(n_edge_number))
	!	
	!	sum_a = (0.0, 0.0)
	!	sum_b = (0.0, 0.0)
	!	sum_c = (0.0, 0.0)
	!	
	!	m_cube = struct_cube%pair_cubes(m)%cube_index	!Cube index where edge m is located	
	!	n_cube = struct_cube%pair_cubes(n)%cube_index
	!	Ra = struct_cube%cubes(m_cube)%cube_position
	!	Rb = struct_cube%cubes(n_cube)%cube_position	
	!	Rab = Ra - Rb
	!	
	!	Lm = abs(k1*d_c) + 10*(abs(k1*d_c))**(1/3)				
	!	TL_1 = TL_km_arr(Rab, Lm, k1)
	!	
	!	Lm = abs(k2*d_c) + 10*(abs(k2*d_c))**(1/3)
	!	TL_2 = TL_km_arr(Rab, Lm, k2)
	!	
	!	sum_a = (0.0, 0.0)
	!	sum_b = (0.0, 0.0)
	!	sum_c = (0.0, 0.0)
	!	do tt = 1, 2
	!		call Receivef_rc_KL(struct, m, tt, N_sp, Ra, Rfk, Rfla, Rflb)
	!		do ll = 1, 2				
	!			call Radiationf_rd_KL(struct, n, ll, N_sp, Rb, Tfla, Tflb)	!Tfk = Tf1a
	!			do j = 1, N_sp
	!				ss1 = dot_product(conjg(Rfk(j, 1, :)), Tfla(j, 1, :))
	!				ss2 = dot_product(conjg(Rfk(j, 2, :)), Tfla(j, 2, :))
	!					
	!				ff1 = omega*dot_product(conjg(Rfla(j, 1, :)), Tfla(j, 1, :))
	!				gg1 = Rflb(j, 1)*Tflb(j, 1)/omega
	!				ff2 = omega*dot_product(conjg(Rfla(j, 2, :)), Tfla(j, 2, :))
	!				gg2 = Rflb(j, 2)*Tflb(j, 2)/omega
	!					
	!				sum_b = sum_b + (TL_1(j)*(my_1*ff1 - gg1/eps_1) + TL_2(j)*(my_2*ff2 - gg2/eps_2))*TRF(j)*im
	!				sum_c = sum_c + (TL_1(j)*(eps_1*ff1 - gg1/my_1) + TL_2(j)*(eps_2*ff2 - gg2/my_2))*TRF(j)*im		
	!				sum_a = sum_a + (TL_1(j)*ss1 + TL_2(j)*ss2)*TRF(j)
	!			end do				
	!		end do				
	!	end do		
 !   end subroutine FMM_edge_Tri
	
	end module FMM_tri_Spherical