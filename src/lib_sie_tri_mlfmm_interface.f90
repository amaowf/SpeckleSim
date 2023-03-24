module lib_sie_tri_mlfmm_interface
	use libtree
	use libmlfmm 
	use lib_sie_math
	use lib_sie_constants
	use lib_sie_tri_data_container
	use lib_sie_type_function
	
	use lib_sie_tri_calculation_mod
	use lib_sie_tri_mlfmm
	use Gauss_Quadrature_Fast
   use lib_sie_mlfmm 
	
	implicit none
	
	private
	
	public :: lib_sie_ml_fmm_calculate_vector_b_tri
	public :: lib_sie_ml_fmm_get_vector_b_tri
	public :: lib_sie_ml_fmm_get_init_vector_x_tri
	public :: lib_sie_ml_fmm_precalculation_tri
	public :: lib_sie_ml_fmm_run_tri
	public :: lib_sie_ml_fmm_save_vector_x_tri
	public :: run_MLFMM_one_forward_tri
	public :: find_elements_in_field_zcsr
	public :: lib_sie_ml_fmm_preconditioner_left
	
	! Test functions
	!public :: Test_TL_at_loc
	!public :: test_lib_sie_get_TL_from_boxi4
	!public :: test_MLFMM_versus_Normal	
	!public :: test_single_element_impedance
	!public :: Test_lagrange_Interpolation_New
	
   type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles
	
   ! Tree parameters
   integer(kind=UINDEX_BYTES) :: m_tree_neighbourhood_size_k
   integer(kind=4) :: m_tree_s_opt
   integer(kind=1) :: m_tree_l_min   
	type(lib_tree_universal_index), dimension(:), allocatable :: element_uindex
	type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
	 
   integer(kind=UINDEX_BYTES), dimension(3) :: m_data_element_hierarchy_type_length

   logical :: m_final_sum_calc_y_hierarchy
   logical :: m_final_sum_calc_xy_hierarchy
   logical :: m_use_own_sum 
	
	contains
	
	subroutine run_MLFMM_one_forward_tri()
		implicit none
		call set_parameters_tri()
		call lib_sie_ml_fmm_precalculation_tri()
		call lib_sie_ml_fmm_run_tri()	
	end subroutine
	
	subroutine lib_sie_ml_fmm_precalculation_tri()
		implicit none
		
		integer :: m, ot		
		double precision, dimension(TREE_DIMENSIONS) :: pos_edge
		type(lib_ml_fmm_data) :: data_elements		
	
		call data_import_tri()	!do not forget to restore it!!!!
		
		call precalculation_fn(struc_tri)
		n_el = size(struc_tri%elements)
		print*, 'n_el=', n_el
		m_pairs = size(struc_tri%neighbours)	
		print*, 'm_pairs=', m_pairs
		allocate(data_elements%xy(m_pairs))	
   
		ot = 1		
		do m = 1, m_pairs
			call fn_edge_center(struc_tri, m, ot, pos_edge)
			data_elements%xy(m)%point_x%x=pos_edge
			data_elements%xy(m)%element_type = 1 ! value .ne. -1
			data_elements%xy(m)%hierarchy = HIERARCHY_XY
		end do
		
		call lib_sie_destructor_tri()		
		call lib_sie_constructor_tri(data_elements, tree_s_opt)
		call Truncation_number_calculation()
		call precalculation_element_coefficient_tri()	
		!
		m_tree_neighbourhood_size_k = lib_ml_fmm_get_neighbourhood_size_k()
		call precalculation_LT3_TL(m_tree_neighbourhood_size_k, m_tree_l_max)
		call precalculation_near_field_csr_tri()
		call precalculation_diagonal_element()
	end subroutine lib_sie_ml_fmm_precalculation_tri
	
	subroutine precalculation_diagonal_element()
		implicit none
		complex(dp), dimension(2,2) :: D_mat_tmp		
		integer :: m
		
		if (allocated(diagonal_element))then
			deallocate(diagonal_element)
		end if
		allocate(diagonal_element(m_pairs*2))
		
		do m = 1, m_pairs
			call find_elements_in_field_zcsr(m, m, D_mat_tmp)
			diagonal_element(m) = D_mat_tmp(1, 1)
			diagonal_element(m + m_pairs) = D_mat_tmp(2, 2)		
		end do
	end subroutine
	
	subroutine lib_sie_constructor_tri(data_elements, s_opt)
		implicit none		
		
		! dummy
		type(lib_ml_fmm_data), intent(inout) :: data_elements
		integer, intent(in) :: s_opt
		    
		! auxiliary
		type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures
		type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures      
    
		ml_fmm_operator_procedures = ml_fmm_type_operator_get_procedures(2)
		ml_fmm_procedures = lib_sie_ml_fmm_get_procedures_tri()
		m_ml_fmm_handles = ml_fmm_procedures
		
		if (data_structure_type .eq. 'bottom_up') then
			call lib_ml_fmm_constructor(data_elements, &
					ml_fmm_operator_procedures, &
					ml_fmm_procedures, &
					tree_s_opt=s_opt, &
					final_sum_calc_y_hierarchy = .false., &
					final_sum_calc_xy_hierarchy = .true., &
					use_own_sum=.true.)
			m_tree_l_max = lib_tree_get_level_max(tree_s_opt)		
			
		else if (data_structure_type .eq. 'top_down') then
			call lib_ml_fmm_constructor(data_elements, &
					ml_fmm_operator_procedures, &
					ml_fmm_procedures, tree_s_opt, &
					tree_bounding_box, tree_l_max_predefined, &
					final_sum_calc_y_hierarchy = .false., &
					final_sum_calc_xy_hierarchy = .true., &
					use_own_sum=.true.)	
			m_tree_l_max = tree_l_max_predefined
		else
			print*, 'Wrong data construction type!'
			call exit
		end if
   end subroutine
		
	subroutine precalculation_element_coefficient_tri()
		implicit none
		
		call find_box_index_for_element(element_uindex)		
		call PreCalculate_RT(truncation_number, struc_tri, preCalculated_RT)
	end subroutine	precalculation_element_coefficient_tri

	subroutine lib_sie_destructor_tri()
		implicit none
		call lib_ml_fmm_destructor()
	end subroutine
	
	!The run without GMRES
	subroutine lib_sie_ml_fmm_run_tri()			
		implicit none
		complex(dp), dimension(:), allocatable :: vector_x
		complex(dp), dimension(:), allocatable :: vector_b
		
		!dummy
		integer :: i
		
		if (allocated(vector_b)) then
			deallocate(vector_b)
		end if		
		allocate(vector_b(2*m_pairs)) 
		
		call lib_sie_ml_fmm_get_vector_b_tri(vector_b)			
		call lib_sie_ml_fmm_get_init_vector_x_tri(vector_x) 
		
		call lib_sie_ml_fmm_calculate_vector_b_tri(vector_x, vector_b)
		
		!open (unit=206, file = 'vector_b_1.txt', action="write",status = 'replace')				
		!do i = 1, m_pairs*2
		!	write (206, '(201(es19.12, tr5))') vector_b(i)
		!end do
		!close(206)
		
	end subroutine
	
	subroutine lib_sie_ml_fmm_get_init_vector_x_tri(vector_x)        
		implicit none
		
		! dummy		
		double complex, dimension(:), allocatable, intent(inout) :: vector_x		
				
		!temporary variable		
		complex(dp), dimension(:), allocatable :: vector_b	!
		
		complex(dp), dimension(:), allocatable :: vector_x_PMCHWT, vector_x_MCTF
		!
		character(len = 50) :: dummy_character
		logical :: x_initial_file
		integer :: dummy_number, m
		real(dp) :: x_initial_re, x_initial_im
		
		if (allocated(vector_x)) then
			deallocate(vector_x)
		end if		
		
		allocate(vector_x(2*m_pairs), vector_x_PMCHWT(2*m_pairs), vector_x_MCTF(2*m_pairs))
		allocate(vector_b(2*m_pairs))		
		
		x_initial_file =  .false. !.true. !
		
		if (x_initial_file) then		
		
			open(unit = 81, file = 'Result_SEH_tri_ICTF.txt', status = 'old', action='read') !							
				read(81, *) dummy_character, dummy_number
				read(81, *) dummy_character, dummy_number								
				do m = 1, 2*m_pairs					
					read(81, *) x_initial_re, x_initial_im
					vector_x(m) = cmplx(x_initial_re, x_initial_im)					
				end do
			close(81)	
			!open(unit = 81, file = 'Result_SEH_tri_MCTF.txt', status = 'old', action='read') !							
			!	read(81, *) dummy_character, dummy_number
			!	read(81, *) dummy_character, dummy_number								
			!	do m = 1, 2*m_pairs					
			!		read(81, *) x_initial_re, x_initial_im
			!		vector_x_MCTF(m) = cmplx(x_initial_re, x_initial_im)					
			!	end do
			!close(81)	
			!vector_x = (vector_x_MCTF + vector_x_PMCHWT)/2
		else
			vector_x(1:m_pairs) = x_initial
			vector_x(m_pairs + 1 : 2*m_pairs) =x_initial				
		end if
		
	end subroutine lib_sie_ml_fmm_get_init_vector_x_tri !		
	
	subroutine lib_sie_ml_fmm_save_vector_x_tri(vector_x)        
		implicit none
       
		! dummy		
		complex(dp), dimension(:), allocatable, intent(in) :: vector_x
		complex(dp), dimension(:), allocatable :: buffer_x
		integer :: m
		
		! more has to be done
		m = size(vector_x)		
		allocate(buffer_x(m))
		buffer_x = vector_x		
		call move_alloc(buffer_x, simulation_data%vector_x)
		
		!open (unit=206, file = 'vector_x.txt', action="write",status = 'replace')				
		!do i = 1, m_pairs*2
		!	write (206, '(201(es19.12, tr5))') vector_x(i)
		!end do
		!close(206)
		
	end subroutine lib_sie_ml_fmm_save_vector_x_tri !
	
	subroutine lib_sie_ml_fmm_get_vector_b_tri(vector_b)
		implicit none
   	! dummy
		double complex, dimension(:), allocatable, intent(inout) :: vector_b		
		
		if (allocated(vector_b)) then
			deallocate(vector_b)
		end if		
		allocate(vector_b(2*m_pairs)) 
	
		!integer :: i
		call get_vector_b_tri(vector_b)			
		
	end subroutine lib_sie_ml_fmm_get_vector_b_tri !	
	
    !
    !        ! Argument
    !        ! ----
    !        !   vector_x: double complex, dimension(:)
    !        !       vector x: M x = b
    !        !   vector_b: double complex, dimension(:)
    !        !       vector x:
	subroutine lib_sie_ml_fmm_calculate_vector_b_tri(vector_x, vector_b)    
      implicit none
      ! dummy
      double complex, dimension(:), allocatable, intent(in) :: vector_x
      double complex, dimension(:), allocatable, intent(inout) :: vector_b

      ! auxiliary
      integer :: i

      type(lib_ml_fmm_v), dimension(:), allocatable :: b
      type(lib_ml_fmm_v), dimension(:), allocatable :: x
      
			
      allocate(x(m_pairs)) 
		  
      do i = 1, m_pairs
         allocate(x(i)%C(2))
         x(i)%C(1)=vector_x(i)
         x(i)%C(2)=vector_x(i + m_pairs)
      end do
        
      call lib_ml_fmm_run(x, b)		
		!save vector_b for the next run
      do i = 1, m_pairs        
         vector_b(i) = b(i)%C(1)
			vector_b(i + m_pairs) = b(i)%C(2)
      end do
		
		!open (unit=206, file = 'vector_b_2.txt', action="write",status = 'replace')				
		!	do i = 1, 2*m_pairs
		!		write (206, '(201(es19.12, tr5))') vector_b(i)
		!	end do
		!close(206)
		
	end subroutine    !
	
	!vector_1: should be vector_x, vector_2: should be M^{-1}
	!The result should be another vector_b
	subroutine lib_sie_ml_fmm_preconditioner_left(vector_1, vector_2)
      implicit none
      ! dummy
      double complex, dimension(:), allocatable, intent(in) :: vector_1
      double complex, dimension(:), allocatable, intent(out) :: vector_2

      ! auxiliary
      integer :: i, m
		m = size(vector_1)		
		allocate(vector_2(m))
		
      do i = 1, m_pairs
			vector_2(i) = vector_1(i)/diagonal_element(i)
			vector_2(i + m_pairs) = vector_1(i + m_pairs)/diagonal_element(i + m_pairs)
		end do
		
	end subroutine    !
    !
	function lib_sie_ml_fmm_get_procedures_tri() result(handle)
      implicit none
      ! dummy
      type(lib_ml_fmm_procedure_handles) :: handle
  
      ! Upward Pass
      ! Step 1
      handle%get_c => lib_sie_ml_fmm_get_C ! 
      ! Step 2
      handle%get_translation_SS => lib_sie_ml_fmm_translation_SS ! 
      !
      !! Downward Pass
      !! Step 1
      handle%get_translation_SR => lib_sie_ml_fmm_translation_SR ! 
      !! Step 2
      handle%get_translation_RR => lib_sie_ml_fmm_translation_RR ! 
        
      ! Final Summation        
      handle%get_v_y_j => lib_sie_ml_fmm_get_v_y_j !
	end function
	 
	 !subroutine get_truncation_number(l, K_m)
		!integer(kind = 1), intent(in) :: l
		!integer :: K_m
		!
		!m = size(truncation_number(1, :))
		!do i = 1, m
		!
		!
	 !end subroutine
	 
    !
    !        ! todo: define procedures
    !        ! - get_c
    !        ! - get_translation_SS
    !        ! - get_translation_SR
    !        ! - get_translation_RR
    !        ! - get_v_y_j
    !
    !		  ! Argument
    !        ! ----
    !        !   x: type(lib_tree_spatial_point)
    !        !       centre of the box
    !        !   data_element: type(lib_tree_data_element), dimension(:)
    !        !       data element list of the box
    !        !   element_number: integer(kind=4), dimension(:)
    !        !       number of the element
    !        !       HINT: X, Y, and XY lists are concatenated
    !        !
    !        ! Returns
    !        ! ----
    !        !   C: type(lib_ml_fmm_coefficient)
    !        !       coefficient of the box
   
	function lib_sie_ml_fmm_get_c(x_c, data_element, element_number) result(C_1)  
		use lib_ml_fmm_data_container
      implicit none
      ! dummy
      type(lib_tree_spatial_point), intent(in) :: x_c
      type(lib_tree_data_element), dimension(:), intent(in) :: data_element !
      integer(kind=4), dimension(:), intent(in) :: element_number
      type(lib_ml_fmm_coefficient) :: C_1, C1_tmp
        
      ! auxiliary
      integer :: i, K_m, m_element_number, n_arr
      !logical :: near_field
		type(lib_ml_fmm_coefficient) :: Fn_scs!
      type(lib_tree_spatial_point) :: data_element_x_unscaled
		integer :: Ns, k !nn,
      integer(kind = 1) :: l 
		complex(dp) :: xe, xh
		  
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		l = m_tree_l_max
		K_m = truncation_number(m_tree_l_max)%n
		Ns = K_m*K_m*2
		!
		!!call Khat_and_TRF_ExtendedK(K_m)
		n_arr = 4
		call init_list_sie(C_1, n_arr, Ns+2)
		call init_list_sie(C1_tmp, n_arr, Ns+2)
		n_arr = 2
		call init_list_sie(Fn_scs, n_arr, Ns+2)
		
		!+++++++++++++		
      data_element_x_unscaled = lib_tree_get_unscaled_point(x_c)		
		do i = 1, size(element_number)			
			m_element_number = element_number(i)
			
			xe = m_ml_fmm_u(m_element_number)%C(1)
			xh = m_ml_fmm_u(m_element_number)%C(2)
			
			do k = 1, 2 !k1 or k2
				fn_scs%a_nm%item(k)%item(:) = conjg(preCalculated_RT(m_element_number)%a_nm%item(k)%item(:))
				fn_scs%b_nm%item(k)%item(:) = conjg(preCalculated_RT(m_element_number)%b_nm%item(k)%item(:))

				!LE, KE*xe
				C_1%a_nm%item(k)%item(:) = C_1%a_nm%item(k)%item(:) + fn_scs%a_nm%item(k)%item(:)*xe!
				C_1%b_nm%item(k)%item(:) = C_1%b_nm%item(k)%item(:) + fn_scs%b_nm%item(k)%item(:)*xe!
				
				!KH, LH*xh
				C_1%a_nm%item(k+2)%item(:) = C_1%a_nm%item(k+2)%item(:) + fn_scs%a_nm%item(k)%item(:)*xh
				C_1%b_nm%item(k+2)%item(:) = C_1%b_nm%item(k+2)%item(:) + fn_scs%b_nm%item(k)%item(:)*xh
			end do! 
		end do
		!call exit
		!call fx_fy(C_1)
    end function
  
    ! Translation: far-to-far (Singular-to-Singular)
    !
    ! Arguments
    ! ----
    !   B_i_1
    !       set of expansion coefficients
    !   x_1: spatial point
    !       origin of coordinate system 1
    !   x_2: spatial point
    !       origin of coordinate system 2
    !
    ! Returns
    ! ----
    !   B_i_2
    !       set of expansion coefficients
    !
	function lib_sie_ml_fmm_translation_SS(C_1, x_1, x_2) result(B_i_2)
		implicit none
		! input
		type(lib_ml_fmm_coefficient), intent(in) :: C_1
		type(lib_tree_spatial_point), intent(in) :: x_1 !local
		type(lib_tree_spatial_point), intent(in) :: x_2 !parent
		
		! output
		type(lib_ml_fmm_coefficient) :: B_i_2
		
		! dummy
		type(lib_ml_fmm_coefficient) :: C1_inter
		type(lib_tree_spatial_point) :: r !
		real(dp) :: tmp
		integer :: Ns, K_m, K_n
		integer :: i, pp, j
		integer(kind = 1) :: l
		complex(dp), dimension(2) :: tmp_k
		
      B_i_2%uindex%l = C_1%uindex%l - int(1, 1)
		l = C_1%uindex%l	
		
		!in Gibson, r= r_local - r_parent
		r = lib_tree_get_unscaled_point(x_1) - lib_tree_get_unscaled_point(x_2) 
		K_m = truncation_number(l)%n
		K_n = truncation_number(l-1)%n
		!print*, 'K_n', K_n
		!+++++++++		
		if (K_n .gt. K_m) then
			pp = 3
			Ns = 2*K_n*K_n

			call init_list_sie(B_i_2, 4, Ns+2)
			call Khat_and_TRF_ExtendedK(K_n)
			call Lagrange_Interpolation_New(pp, K_m, K_n, C_1, C1_inter)			
			!+++++++++++++++++++++++++++++++++++++++++++
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(im*k1*tmp), exp(im*k2*tmp)/)
				do j = 1, 2
					B_i_2%a_nm%item(j)%item(i) = C1_inter%a_nm%item(j)%item(i)*tmp_k(j)
					B_i_2%b_nm%item(j)%item(i) = C1_inter%b_nm%item(j)%item(i)*tmp_k(j)
					B_i_2%a_nm%item(j+2)%item(i) = C1_inter%a_nm%item(j+2)%item(i)*tmp_k(j)
					B_i_2%b_nm%item(j+2)%item(i) = C1_inter%b_nm%item(j+2)%item(i)*tmp_k(j)
				end do
			end do		
		else			
			Ns = 2*K_m*K_m
			call init_list_sie(B_i_2, 4, Ns+2)
		
			!call Khat_and_TRF_ExtendedK(K_m)
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(im*k1*tmp), exp(im*k2*tmp)/)!the signs !!!!Eq.(9.72) by Gibson
				do j = 1, 2					
					B_i_2%a_nm%item(j)%item(i) = C_1%a_nm%item(j)%item(i)*tmp_k(j)
					B_i_2%b_nm%item(j)%item(i) = C_1%b_nm%item(j)%item(i)*tmp_k(j)			! xe
					B_i_2%a_nm%item(j+2)%item(i) = C_1%a_nm%item(j+2)%item(i)*tmp_k(j)	! xh
					B_i_2%b_nm%item(j+2)%item(i) = C_1%b_nm%item(j+2)%item(i)*tmp_k(j)					
				end do
			end do
		end if
	end function
    !
    ! Translation: far-to-local (Singular-to-Regular)
    !
    ! Arguments
    ! ----
    !   B_i_1
    !       set of expansion coefficients
    !   x_1: spatial point
    !       origin of coordinate system 1
    !   x_2: spatial point
    !       origin of coordinate system 2
    !
    ! Returns
    ! ----
    !   A_i_1
    !       set of expansion coefficients
    !
    !!++ In the normal calculation procedure, x_1 and x_2 are from the same level l. 
	function lib_sie_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_1)		
		implicit none
		
		! input
		type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
		type(lib_tree_spatial_point), intent(in) :: x_1 
		type(lib_tree_spatial_point), intent(in) :: x_2 !neighbour
		
		! output
		type(lib_ml_fmm_coefficient) :: A_i_1
		! dummy
		integer :: Ns, K_m, j !K_n, 
		integer(kind = 1) :: l
		type(lib_tree_spatial_point) :: ra, rb, r_ab
		type(lib_ml_fmm_coefficient) :: TL_pre
		complex(dp), dimension(:, :), allocatable :: TL
				
		l = B_i_1%uindex%l
		A_i_1%uindex%l = l
		j = l
		
		K_m = truncation_number(l)%n		
		Ns = K_m*K_m*2
		
		ra = lib_tree_get_unscaled_point(x_1)
		rb = lib_tree_get_unscaled_point(x_2)
		
		r_ab = ra - rb
		TL_pre = LT3_get_TL_at_loc(ra%x, rb%x, j)
		!TL = TL_km_arr_II(r_ab%x, K_m)
		call init_list_sie(A_i_1, 4, Ns+2)	
		
		!Anterpolation, w will be multiplied at the highest sampling points		
		do j = 1, 2
			A_i_1%a_nm%item(j)%item(:) = B_i_1%a_nm%item(j)%item(:)*TRF_p(:)*TL_pre%a_nm%item(j)%item(:)
			A_i_1%b_nm%item(j)%item(:) = B_i_1%b_nm%item(j)%item(:)*TRF_p(:)*TL_pre%a_nm%item(j)%item(:)
			A_i_1%a_nm%item(j+2)%item(:) = B_i_1%a_nm%item(j+2)%item(:)*TRF_p(:)*TL_pre%a_nm%item(j)%item(:)
			A_i_1%b_nm%item(j+2)%item(:) = B_i_1%b_nm%item(j+2)%item(:)*TRF_p(:)*TL_pre%a_nm%item(j)%item(:)
		end do
	end function
    !
    !! Translation: local-to-local (Regular-to-Regular)
    !    !
    !    ! Arguments
    !    ! ----
    !    !   A_i_1
    !    !   x_1
    !    !   x_2
    !    !
    !    ! Returns
    !    ! ----
    !    !   A_i_2
    !    !
   function lib_sie_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
      implicit none
		
      ! input
      type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
      type(lib_tree_spatial_point), intent(in) :: x_1 !parent
      type(lib_tree_spatial_point), intent(in) :: x_2
		! output
      type(lib_ml_fmm_coefficient) :: A_i_2, A1_inter, A2_tmp
                
      type(lib_tree_spatial_point) :: r
      real(dp) :: tmp
		complex(dp) :: tmp_k(2)
      integer :: Ns, K_m, K_n
      integer :: i, pp, j!
      integer(kind = 1) :: l
         
		l = A_i_1%uindex%l
		A_i_2%uindex%l = l + int(1, 1)
		  
      K_m = truncation_number(l)%n
		K_n = truncation_number(l+1)%n
		r = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1) !
		if ((K_n .lt. K_m) .and. (K_n .ne. 0)) then
			Ns = K_m*K_m*2
			call Khat_and_TRF_ExtendedK(K_m)
			call init_list_sie(A2_tmp, 4, Ns+2)
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
				do j = 1, 2
					A2_tmp%a_nm%item(j)%item(i) = A_i_1%a_nm%item(j)%item(i)*tmp_k(j)
					A2_tmp%b_nm%item(j)%item(i) = A_i_1%b_nm%item(j)%item(i)*tmp_k(j)
					A2_tmp%a_nm%item(j+2)%item(i) = A_i_1%a_nm%item(j+2)%item(i)*tmp_k(j)
					A2_tmp%b_nm%item(j+2)%item(i) = A_i_1%b_nm%item(j+2)%item(i)*tmp_k(j)
				end do  
			end do			
			pp = 3 
			call fx_fy(A2_tmp)			
			call Lagrange_Interpolation_New(pp, K_m, K_n, A2_tmp, A_i_2)				
		else if (K_n .eq. K_m) then
			Ns = K_m*K_m*2
			call init_list_sie(A_i_2, 4, Ns+2)
			!call Khat_and_TRF_ExtendedK(K_m)
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
				do j = 1, 2
					A_i_2%a_nm%item(j)%item(i) = A_i_1%a_nm%item(j)%item(i)*tmp_k(j)
					A_i_2%b_nm%item(j)%item(i) = A_i_1%b_nm%item(j)%item(i)*tmp_k(j)
					A_i_2%a_nm%item(j+2)%item(i) = A_i_1%a_nm%item(j+2)%item(i)*tmp_k(j)
					A_i_2%b_nm%item(j+2)%item(i) = A_i_1%b_nm%item(j+2)%item(i)*tmp_k(j)
				end do
			end do
		end if
   end function
    
	 ! Final Summation 	
	 ! Argument
			! ----
			!   data_element: type(lib_tree_data_element), dimension(:)
			!       data element
			!   element_number_i: integer
			!       number of the i-th data element at the concatenated element data list.
			!       CONVENTION:
			!           Internal representation of the element data list
			!           from the 1-st to the N-th element.
			!           ------------------------------------
			!           |1    X    |     Y    |     XY    N|
			!           ------------------------------------
			!           X-, Y-, XY- hierarchy
			!   y_j: type(lib_tree_spatial_point), dimension(:)
			!       scaled point of the j-th data element
			!       HINT: unscale with lib_tree_get_unscaled_point
			!   element_number_j: integer
			!       number of the j-th data element at the concatenated element data list.
			!       CONVENTION:
			!           Internal representation of the element data list
			!           from the 1-st to the N-th element.
			!           ------------------------------------
			!           |1    X    |     Y    |     XY    N|
			!           ------------------------------------
			!           X-, Y-, XY- hierarchy
			!
			! Returns
			! ----
			!   rv: type(lib_ml_fmm_v)
			!       the result of calculation of u_i * phi_i(y_j)
			!
			! Reference:  Data_Structures_Optimal_Choice_of_Parameters_and_C, eq. 38
		function lib_sie_ml_fmm_get_v_y_j(data_element_y_j, element_number_j, data_element_e2, element_number_e2, D) result(rv)
			implicit none
			! dummy
			type(lib_tree_data_element),  intent(in) :: data_element_y_j
			integer(kind=CORRESPONDENCE_VECTOR_KIND), intent(in) :: element_number_j
			type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_element_e2
			integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable, intent(in) :: element_number_e2
			type(lib_ml_fmm_coefficient), intent(in) :: D
         
			! auxiliary 
			type(lib_ml_fmm_v) :: rv!
			type(lib_ml_fmm_v) :: rv_tmp            
			!type(lib_tree_spatial_point) :: data_element_xj_unscaled			
            
			integer :: j, i
			integer(kind = 1) :: l
			integer(kind=4) :: element_number_i
         
			complex(kind= dp) :: sum_near(2), sum_far(2)			
			type(lib_ml_fmm_coefficient) :: ff_rl!
			integer :: K_m, Ns
                        
			l = data_element_y_j%uindex%l			
			j = l !must be changed
			
			K_m = truncation_number(l)%n
			Ns = 2*K_m*K_m
			
			allocate(rv_tmp%C(2))
			allocate(rv%C(2))
			rv_tmp%C(1:2) = (0.0, 0.0)
			rv%C(1:2) = (0.0, 0.0)
					
			call init_list_sie(ff_rl, 2, Ns+2)			
			do i = 1, 2
				ff_rl%a_nm%item(i)%item = preCalculated_RT(element_number_j)%a_nm%item(i)%item
				ff_rl%b_nm%item(i)%item = preCalculated_RT(element_number_j)%b_nm%item(i)%item
			end do
			
			! Far field summation
			call MLFMM_tri_summation_New(K_m, ff_rl, D, sum_far)
		
			!Near field summation
			sum_near(1:2) = (0.0, 0.0)
			do i = 1, size(data_element_e2)
				element_number_i = element_number_e2(i)
				rv_tmp = lib_sie_ml_fmm_get_u_phi_i_j(element_number_j, element_number_i) 
				do j = 1, 2!k1 and k2
					sum_near(j) = sum_near(j) + rv_tmp%C(j)
				end do
			end do
			
			!Final summation
			rv%C(1) = sum_near(1) + sum_far(1)			
			rv%C(2) = sum_near(2) + sum_far(2)
			
			!if (element_number_j .le. 4) then	
			!	print*, 'element_number_j', element_number_j
			!	print*,  ' sum_far(1) =',  rv%C(1)
			!	print*,  ' sum_far(2) =',  rv%C(2)
			!end if
		end function lib_sie_ml_fmm_get_v_y_j
      
		! Near field matrix
		!++++++++++++++++
		function lib_sie_ml_fmm_get_u_phi_i_j(element_number_j, &
																		element_number_i) result(rv)
			use lib_tree_public
			use ml_fmm_type
			use lib_ml_fmm_data_container			
			use lib_sie_tri_calculation_mod
			implicit none
						     
			integer(kind=4), intent(in) :: element_number_i          
			integer(kind=4), intent(in) :: element_number_j
			type(lib_ml_fmm_v) :: rv, buffer_x
			
			! auxiliary
			complex(dp), dimension(2,2) :: D_mat_tmp
			
			allocate(rv%C(2), buffer_x%C(2))
			buffer_x = m_ml_fmm_u(element_number_i)
			  
			call find_elements_in_field_zcsr(element_number_j, element_number_i, D_mat_tmp)			
			rv%C(1) = D_mat_tmp(1, 1)*buffer_x%C(1) + D_mat_tmp(1, 2)*buffer_x%C(2)!
			rv%C(2) = D_mat_tmp(2, 1)*buffer_x%C(1) + D_mat_tmp(2, 2)*buffer_x%C(2)
			
		end function
		
		subroutine precalculation_near_field_csr_tri()
			implicit none
			type(lib_tree_data_element), dimension(:), allocatable :: data_element_e2
			integer, dimension(:), allocatable :: element_number_e2
			complex(dp), dimension(2,2) :: D_mat_tmp 			
			integer, dimension(:), allocatable :: n_arr_tmp			
		 
			real(dp) :: r_ref
			integer :: m, n, ngp, j, n_csr, m_csr, p, q!
			
			m_tree_neighbourhood_size_k = lib_ml_fmm_get_neighbourhood_size_k()
			m_csr = 0
			
			!!$omp parallel default (shared)  &
			!!$omp private (data_element_e2, element_number_e2) &
			!!$omp shared (m, m_csr, m_pairs, m_tree_neighbourhood_size_k)
			!!$omp do	
			do m = 1, m_pairs
				data_element_e2 = lib_tree_get_domain_e2(m_tree_neighbourhood_size_k, element_uindex(m), element_number_e2)
				m_csr = m_csr + size(element_number_e2)
				deallocate(data_element_e2)
				deallocate(element_number_e2)
			end do
			
			!!$omp end do
			!!$omp end parallel
			
			if (allocated(field_csr%a))then
				deallocate(field_csr%a)
			end if
			
			if (allocated(field_csr%b))then
				deallocate(field_csr%b)
			end if
			
			if (allocated(field_csr%c))then
				deallocate(field_csr%c)
			end if
			
			if (allocated(field_csr%d))then
				deallocate(field_csr%d)
			end if
			
			if (allocated(field_csr%ia))then
				deallocate(field_csr%ia)
			end if			

			if (allocated(field_csr%ja))then
				deallocate(field_csr%ja)
			end if			
		
			allocate(field_csr%a(m_csr), field_csr%b(m_csr), field_csr%c(m_csr), field_csr%d(m_csr))
			allocate(field_csr%ia(m_pairs + 1), field_csr%ja(m_csr))
			!
			field_csr%ia(1:m_pairs+1) = 0
			field_csr%ja(1:m_csr) = 0
  
			m_csr = 0
			ngp = 4
			field_csr%ia(1) = 1
			do m = 1, m_pairs	
				n_csr = 0				
				data_element_e2 = lib_tree_get_domain_e2(m_tree_neighbourhood_size_k, element_uindex(m), element_number_e2)
				n_csr = size(element_number_e2)
				allocate(n_arr_tmp(n_csr))
				n_arr_tmp = element_number_e2				
				call sort_edge_index(n_arr_tmp(1:n_csr), n_csr, n_arr_tmp(1:n_csr))				
				field_csr%ja(m_csr + 1 : m_csr + n_csr) = n_arr_tmp(1: n_csr)
				p = struc_tri%neighbours(m)%element(1)
				q = struc_tri%neighbours(m)%element(2)	
				r_ref = vec_len(struc_tri%midpoint(p)%point - struc_tri%midpoint(q)%point)
				
				do j = 1, n_csr					
					n = n_arr_tmp(j)
					D_mat_tmp = lib_sie_tri_get_impedance_edges(m, n, r_ref)						
					field_csr%a(m_csr + j) = D_mat_tmp(1, 1)!
					field_csr%b(m_csr + j) = D_mat_tmp(1, 2)!
					field_csr%c(m_csr + j) = D_mat_tmp(2, 1)!
					field_csr%d(m_csr + j) = D_mat_tmp(2, 2)!
				end do
				m_csr = m_csr + n_csr
				field_csr%ia(m+1) = field_csr%ia(m) + n_csr
				deallocate(n_arr_tmp)				
			end do !loop m	
		end subroutine precalculation_near_field_csr_tri		

		
	!	subroutine find_elements_in_field_csr(i, j, ia, ja, a_csr, x)
	!		implicit none
	!		integer, intent(in) :: i, j
	!		integer, dimension(:), intent(in) :: a_csr
	!	
	!		integer, intent(out) :: x
	!		integer, dimension(:), intent(in) :: ia, ja
	!		integer :: is, it, m		
 !
	!		x = 0
	!		is = ia(i)
	!		it = ia(i+1)-1
	!		do m = is, it
	!			if (ja(m) .eq. j) then
	!				x = a_csr(m)
	!			else				
	!			 cycle	
	!			end if		
	!		end do
	!end subroutine
	

		
	!++++++++++++++++++++++++++++++++++++++++++++++++
		!Test functions
	!++++++++++++++++++++++++++++++++++++++++++++++++	
                                                 
	function get_Km(uindex, d_c, kd, d0) result(k_m)	 
		type(lib_tree_universal_index), intent(in) :: uindex    
		integer(kind = 1) :: l
		complex(dp), intent(in) :: kd
		integer :: k_m, j
		real(dp) :: d_c, d_tmp !diagonal of the box at level l
		integer, intent(in) :: d0		  
		l = uindex%l
		j = l !must be changed        
		d_tmp = abs(d_c/2*kd)
		k_m = d_tmp+1.8*(d0)**(2/3)*(d_tmp)**(1/3)
	end function
			
	function test_lib_sie_ml_fmm_get_c(x_c, data_element, element_number) result(C_1)
		use lib_ml_fmm_data_container
      implicit none
      ! dummy
      type(lib_tree_spatial_point), intent(in) :: x_c
      type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_element !
      integer(kind=4), dimension(:), intent(in) :: element_number
      type(lib_ml_fmm_coefficient) :: C_1
        
      ! auxiliary
      integer :: i, K_m, m_element_number, n_arr
      !logical :: near_field
		type(lib_ml_fmm_coefficient) :: Fn_scs!
      type(lib_tree_spatial_point) :: data_element_x_unscaled
		
		integer :: nn  
      integer :: Ns, k
      integer(kind = 1) :: l 
		
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		l = m_tree_l_max
		K_m = truncation_number(m_tree_l_max)%n
      Ns = K_m*K_m*2
		call Khat_and_TRF_ExtendedK(K_m)
		
		n_arr = 4
		call init_list_sie(C_1, n_arr, Ns+2)
		n_arr = 2
		call init_list_sie(Fn_scs, n_arr, Ns+2)

      data_element_x_unscaled = lib_tree_get_unscaled_point(x_c)
		!if (truncation_number(l)%near_field)
		nn = data_element(1)%uindex%n
		!print*, 'nn', nn
		!print*, 'l', data_element(1)%uindex%l		
		!print*, 'uindex%n of C1', data_element(1)%uindex%n		
		!print*, 'size(element_number)', size(element_number)
		!if ( nn .eq. 201) then	
			do i = 1, size(element_number)			
				m_element_number = element_number(i)
				!print*, 'm_element_number', m_element_number				
				!print*, 'i=', i
				do k = 1, 2 !k1 or k2
					fn_scs%a_nm%item(k)%item(:) = conjg(preCalculated_RT(m_element_number)%a_nm%item(k)%item(:))
					fn_scs%b_nm%item(k)%item(:) = conjg(preCalculated_RT(m_element_number)%b_nm%item(k)%item(:))
					!LE, KE *xe
					C_1%a_nm%item(k)%item(:) = C_1%a_nm%item(k)%item(:) + fn_scs%a_nm%item(k)%item(:)*m_ml_fmm_u(m_element_number)%C(1)!
					C_1%b_nm%item(k)%item(:) = C_1%b_nm%item(k)%item(:) + fn_scs%b_nm%item(k)%item(:)*m_ml_fmm_u(m_element_number)%C(1)!
					!KH, LH *xh
					C_1%a_nm%item(k+2)%item(:) = C_1%a_nm%item(k+2)%item(:) + fn_scs%a_nm%item(k)%item(:)*m_ml_fmm_u(m_element_number)%C(2)!
					C_1%b_nm%item(k+2)%item(:) = C_1%b_nm%item(k+2)%item(:) + fn_scs%b_nm%item(k)%item(:)*m_ml_fmm_u(m_element_number)%C(2)!
				end do! 
			end do
		!end if
		call fx_fy(C_1)
	end function
	 
	subroutine test_MLFMM_versus_Normal()
		implicit none
		
		real(dp) :: a_arr(2)
		complex(dp), dimension(:), allocatable :: D_mat
		integer :: j
		
		!integer:: m_pairs
		!m_pairs = 768		
		allocate(D_mat(2*m_pairs))
		
		
		open(unit = 12, file = 'Dmat_normal_sum_ICTF.txt')
		
		do j = 1, m_pairs*2
			read(12, *) a_arr
			D_mat(j) = cmplx(a_arr(1), a_arr(2))
		end do
		close(12)		
				
		do j = 1, 4!m_pairs
			print*, 'j=', j
			print*, 'D_mat_e(j) =', D_mat(j)
			print*, 'D_mat_h(j) =', D_mat(j + m_pairs)
		end do 

	end subroutine
	
	function test_lib_sie_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
      implicit none
		
      ! input
      type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
      type(lib_tree_spatial_point), intent(in) :: x_1 !parent
      type(lib_tree_spatial_point), intent(in) :: x_2
		! output
      type(lib_ml_fmm_coefficient) :: A_i_2, A1_inter, A2_tmp
                
      type(lib_tree_spatial_point) :: r
      real(dp) :: tmp
		complex(dp) :: tmp_k(2)
      integer :: Ns, K_m, K_n
      integer :: i, pp, j!
      integer(kind = 1) :: l
         
		l = A_i_1%uindex%l
		A_i_2%uindex%l = l + int(1, 1)
		  
      K_m = truncation_number(l)%n
		K_n = truncation_number(l+1)%n
		r = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1) !
		if ((K_n .lt. K_m) .and. (K_n .ne. 0)) then
			Ns = K_m*K_m*2
			call Khat_and_TRF_ExtendedK(K_m)
			call init_list_sie(A2_tmp, 4, Ns+2)
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
				do j = 1, 2
					A2_tmp%a_nm%item(j)%item(i) = A_i_1%a_nm%item(j)%item(i)*tmp_k(j)
					A2_tmp%b_nm%item(j)%item(i) = A_i_1%b_nm%item(j)%item(i)*tmp_k(j)
					A2_tmp%a_nm%item(j+2)%item(i) = A_i_1%a_nm%item(j+2)%item(i)*tmp_k(j)
					A2_tmp%b_nm%item(j+2)%item(i) = A_i_1%b_nm%item(j+2)%item(i)*tmp_k(j)
				end do  
			end do			
				!***********
				 !A2_tmp and A_i_2 are the same when km=14
				print*, 'A2_tmp%a_nm%item(2)%item(2)', A2_tmp%a_nm%item(2)%item(2)
				
				open (unit=206, file = 'A2_interface_km17.txt', action="write",status = 'replace')
				do i = 1, K_m*K_m*2+2
					write (206, '(16(es19.12, tr5))') 			&
					real(A2_tmp%a_nm%item(1)%item(i)), imag(A2_tmp%a_nm%item(1)%item(i)), &
					real(A2_tmp%a_nm%item(2)%item(i)), imag(A2_tmp%a_nm%item(2)%item(i)), &
					real(A2_tmp%b_nm%item(1)%item(i)), imag(A2_tmp%b_nm%item(1)%item(i)), &
					real(A2_tmp%b_nm%item(2)%item(i)), imag(A2_tmp%b_nm%item(2)%item(i)), &
					real(A2_tmp%a_nm%item(3)%item(i)), imag(A2_tmp%a_nm%item(3)%item(i)), &
					real(A2_tmp%a_nm%item(4)%item(i)), imag(A2_tmp%a_nm%item(4)%item(i)), &
					real(A2_tmp%b_nm%item(3)%item(i)), imag(A2_tmp%b_nm%item(3)%item(i)), &
					real(A2_tmp%b_nm%item(4)%item(i)), imag(A2_tmp%b_nm%item(4)%item(i))			
				end do
				close (206)
				
				!***********
				pp = 6 
				call fx_fy(A2_tmp)			
				call Lagrange_Interpolation_New(pp, K_m, K_n, A2_tmp, A_i_2)	
				
				open (unit=206, file = 'A2_interface_kn14.txt', action="write",status = 'replace')
				do i = 1, K_n*K_n*2+2
					write (206, '(16(es19.12, tr5))') 			&
					real(A_i_2%a_nm%item(1)%item(i)), imag(A_i_2%a_nm%item(1)%item(i)), &
					real(A_i_2%a_nm%item(2)%item(i)), imag(A_i_2%a_nm%item(2)%item(i)), &
					real(A_i_2%b_nm%item(1)%item(i)), imag(A_i_2%b_nm%item(1)%item(i)), &
					real(A_i_2%b_nm%item(2)%item(i)), imag(A_i_2%b_nm%item(2)%item(i)), &
					real(A_i_2%a_nm%item(3)%item(i)), imag(A_i_2%a_nm%item(3)%item(i)), &
					real(A_i_2%a_nm%item(4)%item(i)), imag(A_i_2%a_nm%item(4)%item(i)), &
					real(A_i_2%b_nm%item(3)%item(i)), imag(A_i_2%b_nm%item(3)%item(i)), &
					real(A_i_2%b_nm%item(4)%item(i)), imag(A_i_2%b_nm%item(4)%item(i))			
				end do
				close (206)
				!print*, 'A_i_2%a_nm%item(2)%item(2)', A_i_2%a_nm%item(2)%item(2)
				call exit
			else if (K_n .eq. K_m) then
				Ns = K_m*K_m*2
				call init_list_sie(A_i_2, 4, Ns+2)
				call Khat_and_TRF_ExtendedK(K_m)
				do i = 1, Ns+2
					tmp = dot_product(K_hat_p(i, :), r%x)
					tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
					do j = 1, 2
						A_i_2%a_nm%item(j)%item(i) = A_i_1%a_nm%item(j)%item(i)*tmp_k(j)
						A_i_2%b_nm%item(j)%item(i) = A_i_1%b_nm%item(j)%item(i)*tmp_k(j)
						A_i_2%a_nm%item(j+2)%item(i) = A_i_1%a_nm%item(j+2)%item(i)*tmp_k(j)
						A_i_2%b_nm%item(j+2)%item(i) = A_i_1%b_nm%item(j+2)%item(i)*tmp_k(j)
					end do
				end do
				!print*, 'Interpolated onto K_n', K_n
				!print*, 'Interpolated onto K_m', K_m
				!
				!print*, 'A_i_2%a_nm%item(2)%item(2)', A_i_2%a_nm%item(2)%item(2)
				
				open (unit=206, file = 'A_i_2_Interface_km14.txt', action="write",status = 'replace')
				do i = 1, K_n*K_n*2+2
					write (206, '(16(es19.12, tr5))') 			&
					real(A_i_2%a_nm%item(1)%item(i)), imag(A_i_2%a_nm%item(1)%item(i)), &
					real(A_i_2%a_nm%item(2)%item(i)), imag(A_i_2%a_nm%item(2)%item(i)), &
					real(A_i_2%b_nm%item(1)%item(i)), imag(A_i_2%b_nm%item(1)%item(i)), &
					real(A_i_2%b_nm%item(2)%item(i)), imag(A_i_2%b_nm%item(2)%item(i)), &
					real(A_i_2%a_nm%item(3)%item(i)), imag(A_i_2%a_nm%item(3)%item(i)), &
					real(A_i_2%a_nm%item(4)%item(i)), imag(A_i_2%a_nm%item(4)%item(i)), &
					real(A_i_2%b_nm%item(3)%item(i)), imag(A_i_2%b_nm%item(3)%item(i)), &
					real(A_i_2%b_nm%item(4)%item(i)), imag(A_i_2%b_nm%item(4)%item(i))			
				end do
				close (206)
				call exit
			else
				print*, 'No machted Kn and Km'
			end if
   end function
	
	subroutine Test_TL_at_loc()
		integer :: K_m, j
		complex(dp), dimension(:, :), allocatable :: TL
		
		type(lib_tree_universal_index) :: uindex, uindex_tmp
		type(lib_tree_spatial_point) :: ra, rb, rab
		type(lib_ml_fmm_coefficient) :: TL_pre
		
		uindex%l = 2
		uindex%n = 5
		
		uindex_tmp = uindex
		print*, 'x1, uindex_tmp' , uindex_tmp
		ra = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex_tmp))
		K_m = 13	
		
		call Khat_and_TRF_ExtendedK(K_m)
		
		uindex%l = 2
		uindex%n = 13		
		uindex_tmp = uindex
		print*, 'x2, uindex_tmp' , uindex_tmp
		
		rb = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex_tmp))
		
		ra%x = (/-1.856986195689589E-007, -6.189953985631962E-008,  6.189953985631972E-008/)
		rb%x = (/-6.189953985631962E-008, -6.189953985631962E-008, -1.856986195689589E-007/)
		rab = ra -rb
		
		print*, 'rab in test', rab
		TL = TL_km_arr_II(rab%x, K_m)
		j = uindex_tmp%l
		TL_pre = LT3_get_TL_at_loc(ra%x, rb%x, j)
		
		print*, 'TL(15, 2) =', TL(15, 2)
		print*, 'TL_pre =', TL_pre%a_nm%item(2)%item(15)
	end subroutine
	
	subroutine test_single_element_impedance()
		call Test_PreCalculate_RT_02(struc_tri)
	end subroutine 
	
	subroutine test_lib_sie_get_TL_from_boxi4() 
		
		implicit none
		! dummy
		type(lib_tree_universal_index) :: uindex
			
		! auxiliary			
		type(lib_tree_universal_index), dimension(:), allocatable :: boxes_i4
		type(lib_tree_spatial_point) :: x_c
		!type(lib_tree_spatial_point) :: x_i4
		
		!complex(dp), dimension(:, :), allocatable :: TL
		double precision :: d_c
		integer :: Ns, Lm, n
		integer(kind = 1) :: l
		type(lib_tree_spatial_point) :: r_ab
		
		integer :: K_m, j !sd,

		l = 3
		uindex%l = l
		do n = 1, 500
		uindex%n = n
		m_tree_neighbourhood_size_k = lib_ml_fmm_get_neighbourhood_size_k()
		!allocate(boxes_i4, source=lib_tree_get_domain_i4(uindex, m_tree_neighbourhood_size_k))
		boxes_i4 = lib_tree_get_domain_i4(uindex, m_tree_neighbourhood_size_k)
		if (size(boxes_i4) .gt. 219) then
		print*, 'n=', n
		print*, 'size(boxes_i4)', size(boxes_i4)
	   end if
		
		end do
		x_c = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex))!lib_tree_get_unscaled_point
		!print*, 'unscaled x_c', lib_tree_get_centre_of_box(uindex)
		
		K_m = truncation_number(l)%n
		Lm = K_m		
		Ns = K_m*K_m*2
		j = uindex%l
		d_c = lib_tree_get_box_edge_length(int(j, 1))*lib_tree_scaling_D%x(1)
		
		!call Khat_and_TRF_ExtendedK(K_m)
		
		
		!do i = 1, size(boxes_i4)
		!	x_i4 =lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(boxes_i4(i)))!lib_tree_get_unscaled_point			
		!	print*, 'uindex=',  boxes_i4(i)			
		!	
		!	r_ab = x_c-x_i4
		!	if (boxes_i4(i)%n .eq. 21) then
		!	!print*, 'rab =', r_ab%x
		!	!print*, 'x_i4%x', x_i4%x
		!	!print*, 'unscaled x_i4', lib_tree_get_centre_of_box(boxes_i4(i))
		!	!print*, 'scaling factor', lib_tree_scaling_D%x(1)
		!	!print*, 'scaling factor', lib_tree_scaling_D%x(2)
		!	!print*, 'unscaled distance', lib_tree_get_centre_of_box(uindex)-lib_tree_get_centre_of_box(boxes_i4(i))
		!	end if 
			!TL = TL_km_arr_II(r_ab%x, K_m)
		!end do
				
		!l = uindex%l			
		!j = l
		!d_c = lib_tree_get_box_diagonal(j)
		!print*, 'uindex%l=', uindex%l
		!print*, 'uindex%n=', uindex%n
			
	end subroutine test_lib_sie_get_TL_from_boxi4
	
	subroutine Test_lagrange_Interpolation_New()
		implicit none
		call Test_lagrange_interpolation(struc_tri)
	end subroutine

	
	end module lib_sie_tri_mlfmm_interface

