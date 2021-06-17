    module lib_sie_ml_fmm_interface
    use libtree
    use libmlfmm
    use disc_mod
    use input_tri

    use lib_sie_data_container
    
    implicit none
	 
    private

    public :: lib_sie_constructor
    public :: Truncation_number_calculation
    public :: lib_sie_ml_fmm_calculate_vector_b
	 public :: preCalculation_element_coefficient
	 public :: test_MLFMM_versus_Normal
	 public :: test_find_box_index_for_element
	 
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles
    !type(lib_ml_fmm_hierarchy), dimension(:), allocatable :: m_ml_fmm_hierarchy

    ! Tree parameters
    integer(kind=UINDEX_BYTES) :: m_tree_neighbourhood_size_k
    integer(kind=4) :: m_tree_s_opt
    integer(kind=1) :: m_tree_l_min
    integer(kind=1) :: m_tree_l_max
	 
    integer(kind=UINDEX_BYTES), dimension(3) :: m_data_element_hierarchy_type_length

    logical :: m_final_sum_calc_y_hierarchy
    logical :: m_final_sum_calc_xy_hierarchy
    logical :: m_use_own_sum   
	 
	contains
	
	subroutine test_MLFMM_versus_Normal()
		implicit none
		
		real(kind = 8) :: a_arr(2)
		complex(kind = 8), dimension(:), allocatable :: D_mat
		complex(kind = 8) :: sum_rr(3)
		integer :: i, j, m
		intrinsic :: cmplx
		
		allocate(D_mat(2*m_pairs))

		open(unit = 12, file = 'D_mat.txt')
		
		do j = 1, 4
			read(12, *) a_arr
			D_mat(j) = cmplx(a_arr(1), a_arr(2))			
		end do
		close(12)
		print*, 'D_mat_e(1) =', D_mat(1)		
		print*, 'D_mat_e(2) =', D_mat(2)
		print*, 'D_mat_e(3) =', D_mat(3)
		print*, 'D_mat_e(4) =', D_mat(4)
	end subroutine
    !        
	subroutine lib_sie_constructor(data_elements, s_opt)
		use libmlfmm
		use lib_sie_mlfmm_tri
		use libtree
		use calculation_mod_tri
		use input_tri
		implicit none
		
		! dummy
		type(lib_ml_fmm_data), intent(inout) :: data_elements
		integer, intent(in) :: s_opt
    
		! auxiliary
		type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures
		type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures      
    
		ml_fmm_operator_procedures = ml_fmm_type_operator_get_procedures(2)
		ml_fmm_procedures = lib_sie_ml_fmm_get_procedures()
		m_ml_fmm_handles = ml_fmm_procedures
		call lib_ml_fmm_constructor(data_elements, &
            ml_fmm_operator_procedures, &
            ml_fmm_procedures, &
            tree_s_opt=s_opt, &
            final_sum_calc_y_hierarchy = .false., &
            final_sum_calc_xy_hierarchy = .true., &
            use_own_sum=.true.)
   end subroutine
	
	subroutine Truncation_number_calculation()
		use libmlfmm
		use lib_sie_mlfmm_tri
		use lib_sie_data_container
		
		use libtree
		use calculation_mod_tri		
						
		! auxiliary
		integer :: k, K_m, counter
		integer(kind = 1) :: l
      type(lib_tree_spatial_point) :: x_c, ra, rb, rab
		logical :: near_field
		real(kind = 8) :: d_c, error, r_max
		intrinsic sqrt
				
		error = 1.0e-4
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		k = m_tree_l_max
		
		allocate(truncation_number(k))
		truncation_number(:)%n = 0
		truncation_number(:)%near_field = .false.
		
		d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
		
		!Assuming that only at the lowest level, k2 is considered. 
		!R_min is the smallest effective distance to consider k2			
		if (imag(k2) .ne. 0.0) then
			k = m_tree_l_max-1
			call Truncation_number_near(error, d_c, r_max, near_field, K_m)
			truncation_number(m_tree_l_max)%n = K_m
			truncation_number(m_tree_l_max)%near_field = near_field
			do
				if (k .lt. 2) exit
				d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
				if (d_c*2 .lt. r_max) then 
					call Truncation_number_near(error, d_c, r_max, near_field, K_m)
					truncation_number(k)%n = K_m
					truncation_number(k)%near_field = near_field	
				else			
					call Truncation_number_far(error, d_c, truncation_number(k+1)%n, K_m) !
					truncation_number(k)%n = K_m
				end if
				k = k-1
			end do
			print*, 'truncation_number =', truncation_number
		else
			k = m_tree_l_max
			print*, 'k=', k
			do
				if (k .lt. 2) exit
				d_c = lib_tree_get_box_edge_length(k)*lib_tree_scaling_D%x(1)
				call Truncation_number_far(error, d_c, truncation_number(k)%n, K_m) !
				truncation_number(k)%n = K_m				
				k = k-1
			end do
		end if
		
		!truncation_number(m_tree_l_max)%n = truncation_number(m_tree_l_max)%n
		!truncation_number(m_tree_l_max-1)%n = truncation_number(m_tree_l_max-1)%n
		do k = 1, m_tree_l_max
			truncation_number(k)%n = 13!truncation_number(k)%n + 1
		end do
		print*, 'truncation_number =', truncation_number
	end subroutine Truncation_number_calculation
	
	subroutine preCalculation_element_coefficient()
		use lib_sie_mlfmm_tri
		call PreCalculate_RT(truncation_number, struc, preCalculated_RT)
	end subroutine
	
	subroutine test_find_box_index_for_element()
		use lib_tree
		use lib_sie_data_container
		use libmlfmm
		
		implicit none

		type(lib_tree_universal_index) :: uindex
		
      ! auxiliaray
      integer(kind=UINDEX_BYTES) :: number_of_boxes, m_element_number
		type(lib_tree_data_element), dimension(:), allocatable :: data_element
		type(lib_tree_spatial_point) :: x_c
		
		integer(kind=4), dimension(:), allocatable :: element_number		
      integer(kind=UINDEX_BYTES) :: i, j
      integer(kind=UINDEX_BYTES) :: list_index
      integer(kind=1) :: hierarchy_type
		
		allocate(struc%x_c(m_pairs))				
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		number_of_boxes = size(m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index)
		!print*, 'number_of boxes', number_of_boxes
		uindex%l = m_tree_l_max
		if (m_ml_fmm_hierarchy(m_tree_l_max)%is_hashed) then
            !!$OMP PARALLEL DO PRIVATE(i, hierarchy_type, C) &
            !!$OMP  FIRSTPRIVATE(uindex)
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
										struc%x_c(m_element_number)%point = x_c%x
									end do								
							end if				
                end if
            end do
            !!$OMP END PARALLEL DO
        else
            !!$OMP PARALLEL DO PRIVATE(i, list_index, hierarchy_type, C) &
            !!$OMP  FIRSTPRIVATE(uindex)
            do i=1, number_of_boxes
					uindex%n = i - int(1, 1)                
					data_element = lib_tree_get_domain_e1(uindex, element_number)
					if ((allocated (data_element)) &
							.and. (size(data_element) .gt. 0)) then
						do j = 1, size(element_number)
							m_element_number = element_number(j)
							x_c = lib_tree_get_unscaled_point(lib_tree_get_centre_of_box(uindex))
							struc%x_c(m_element_number)%point = x_c%x
						end do								
					end if
            end do
            !!$OMP END PARALLEL DO
        end if
	end subroutine
		
	subroutine lib_sie_destructor()
		implicit none
		call lib_ml_fmm_destructor()
	end subroutine
    !
    !        ! Argument
    !        ! ----
    !        !   vector_x: double complex, dimension(:)
    !        !       vector x: M x = b
    !        !   vector_b: double complex, dimension(:)
    !        !       vector x:
	subroutine lib_sie_ml_fmm_calculate_vector_b(vector_x, vector_b)
        use lib_ml_fmm
        use lib_sie_data_container
        use lib_ml_fmm_data_container
        
        implicit none
        ! dummy
        double complex, dimension(:), allocatable, intent(in) :: vector_x
        double complex, dimension(:), allocatable, intent(inout) :: vector_b
        double complex, dimension(:), allocatable :: error_b
        ! auxiliary
        integer :: i, j

        type(lib_ml_fmm_v), dimension(:), allocatable :: b
        type(lib_ml_fmm_v), dimension(:), allocatable :: x
        
        allocate(x(m_pairs)) 
        allocate(error_b(2*m_pairs)) 
        
        do i = 1, m_pairs
            allocate(x(i)%C(2))
            x(i)%C(1)=vector_x(i)
            x(i)%C(2)=vector_x(i + m_pairs)
        end do
        
        call lib_ml_fmm_run(x, b)
        !print*, 'b=', b(1)%C(1)
        !do i = 1, size(vector_x)/2
        !    error_b(i) = (b(i)%C(1) - vector_b(i))/vector_b(i)
        !    !print*, 'vector_b =', vector_b(i)
        !    !print*, 'b =', b(i)%C(1)
        !end do
        
        ! todo: post-process: b
        ! v: 
	end subroutine
    !
    !
	function lib_sie_ml_fmm_get_procedures() result(handle)
        implicit none
        ! dummy
        type(lib_ml_fmm_procedure_handles) :: handle
  
        ! Upward Pass
        ! Step 1
        handle%get_c => lib_sie_ml_fmm_get_C ! todo: replace null() with own function ..done
        ! Step 2
        handle%get_translation_SS => lib_sie_ml_fmm_translation_SS ! todo: replace null() with own function
        !
        !! Downward Pass
        !! Step 1
        handle%get_translation_SR => lib_sie_ml_fmm_translation_SR ! todo: replace null() with own function
        !! Step 2
        handle%get_translation_RR => lib_sie_ml_fmm_translation_RR ! todo: replace null() with own function
        
        ! Final Summation        
        handle%get_v_y_j => lib_sie_ml_fmm_get_v_y_j !null() ! todo: replace null() with own function
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
		use lib_tree_public
      use ml_fmm_type
      use lib_ml_fmm
      use lib_sie_mlfmm_tri
      use lib_sie_data_container
      use lib_ml_fmm_data_container
		use disc_mod
  
      implicit none
      ! dummy
      type(lib_tree_spatial_point), intent(in) :: x_c
      type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_element !
      integer(kind=4), dimension(:), intent(in) :: element_number
      type(lib_ml_fmm_coefficient) :: C_1
        
      ! auxiliary
      integer :: i, K_m, m_element_number, n_arr
      logical :: interpol, near_field
		type(lib_ml_fmm_coefficient) :: Fn_scs!
		
		integer :: nn
  
      type(lib_tree_spatial_point) :: data_element_x_unscaled
      integer :: Ns, k
      integer(kind = 1) :: l   
		
        
		m_tree_l_max = lib_tree_get_level_max(tree_s_opt)
		l = m_tree_l_max
		K_m = truncation_number(m_tree_l_max)%n
      Ns = K_m*K_m*2
		call Khat_and_TRF_legendre(K_m)
		
		n_arr = 4
		call init_list_sie(C_1, n_arr, Ns+2)
		n_arr = 2
		call init_list_sie(Fn_scs, n_arr, Ns+2)

      data_element_x_unscaled = lib_tree_get_unscaled_point(x_c)
		!if (truncation_number(l)%near_field)
		nn = data_element(1)%uindex%n
		if ( nn .eq. 80) then
			print*, 'uindex%n', data_element(1)%uindex%n		
			print*, 'element_number', element_number
			do i = 1, size(element_number)			
				m_element_number = element_number(i)
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
		end if
		!print*, 'C_1%a_nm%item(1)%item(12)', C_1%a_nm%item(1)%item(12)
		call fx_fy(C_1)
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
		use lib_tree_public
		use ml_fmm_type
		use lib_sie_mlfmm_tri
		use lib_sie_data_container
		use calculation_mod_tri
		
		implicit none
		! dummy
		type(lib_ml_fmm_coefficient), intent(in) :: C_1
		type(lib_tree_spatial_point), intent(in) :: x_1 !parent
		type(lib_tree_spatial_point), intent(in) :: x_2 !child
		type(lib_ml_fmm_coefficient) :: B_i_2, C1_inter
		
		     
		type(lib_tree_spatial_point) :: ra, rb, r
		real(kind = dp) :: tmp, d_c
		integer :: Ns, K_m, K_D, Ns_p, K_n
		integer :: i, m, pp, j
		integer(kind = 1) :: l
		complex(kind = 8), dimension(2) :: tmp_k
      
		type(lib_ml_fmm_coefficient) :: C1_tmp
		
      B_i_2%uindex%l = C_1%uindex%l - int(1, 1)		
		l = C_1%uindex%l
		j = l
		d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
		!in Gibson, r= r_child - r_parent
		r = lib_tree_get_unscaled_point(x_2) - lib_tree_get_unscaled_point(x_1) 
		K_m = truncation_number(l)%n
		K_n = truncation_number(l-1)%n
		
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
					B_i_2%a_nm%item(j)%item(i) = tmp_k(j)*C1_inter%a_nm%item(j)%item(i)
					B_i_2%b_nm%item(j)%item(i) = tmp_k(j)*C1_inter%b_nm%item(j)%item(i)
					B_i_2%a_nm%item(j+2)%item(i) = tmp_k(j)*C1_inter%a_nm%item(j+2)%item(i)
					B_i_2%b_nm%item(j+2)%item(i) = tmp_k(j)*C1_inter%b_nm%item(j+2)%item(i)
				end do
			end do
		else
			Ns = 2*K_m*K_m
			call init_list_sie(B_i_2, 4, Ns+2)
			do i = 1, Ns+2
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(im*k1*tmp), exp(im*k2*tmp)/)!the signs !!!!Eq.(9.72) by Gibson
				do j = 1, 2
					B_i_2%a_nm%item(j)%item(i) = tmp_k(j)*C_1%a_nm%item(j)%item(i)
					B_i_2%b_nm%item(j)%item(i) = tmp_k(j)*C_1%b_nm%item(j)%item(i) !xe
					B_i_2%a_nm%item(j+2)%item(i) = tmp_k(j)*C_1%a_nm%item(j+2)%item(i) !xh
					B_i_2%b_nm%item(j+2)%item(i) = tmp_k(j)*C_1%b_nm%item(j+2)%item(i)					
				end do
			end do			
		end if
		!print*, 'B_i_2%a_nm%item(1)%item(12)', B_i_2%a_nm%item(1)%item(12)
	end function
    !
    !! Translation: far-to-local (Singular-to-Regular)
    !!
    !! Arguments
    !! ----
    !!   B_i_1
    !!       set of expansion coefficients
    !!   x_1: spatial point
    !!       origin of coordinate system 1
    !!   x_2: spatial point
    !!       origin of coordinate system 2
    !!
    !! Returns
    !! ----
    !!   A_i_1
    !!       set of expansion coefficients
    !!
    !!++ In the normal calculation procedure, x_1 and x_2 are from the same level l. 
	function lib_sie_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_1)
		use lib_tree_public
		use ml_fmm_type
		use lib_sie_mlfmm_tri
		use calculation_mod_tri
		use input_tri
		
		implicit none
		! dummy
		type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
		type(lib_tree_spatial_point), intent(in) :: x_1
		type(lib_tree_spatial_point), intent(in) :: x_2
		type(lib_tree_spatial_point) :: ra, rb
		type(lib_ml_fmm_coefficient) :: A_i_1
    
		double precision :: d_c, r_tmp	!problem with it, how can one obtain in a fast way?
		integer :: Ns, Lm, K_m, K_n, i, j, counter 
		integer(kind = 1) :: l, int_zero
		type(lib_tree_spatial_point) :: r_ab
		complex(kind = dp), dimension(:, :), allocatable :: TL
		intrinsic :: sqrt
		
		l = B_i_1%uindex%l
		A_i_1%uindex%l= l
		j = l
		d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
		!ra is receive cluster center (x1) while rb the radiation cluster center (x2)
		ra = lib_tree_get_unscaled_point(x_1)
		rb = lib_tree_get_unscaled_point(x_2)
		r_ab = ra-rb
		K_m = truncation_number(l)%n		
		
		Lm = K_m		
		Ns = K_m*K_m*2
		r_ab = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1) !x_2 target?!
		call Khat_and_TRF_ExtendedK(K_m)
		TL = TL_km_arr_II(r_ab%x, K_m)
		!print*, 'r_ab%x', vec_len(r_ab%x)
		!print*, 'ratio of rab/d_c =', vec_len(r_ab%x)/d_c
		
		!Km and m are separately calculated
		!TL = TL_km_arr_III(r_ab%x, K_m, d_c)
		call init_list_sie(A_i_1, 4, Ns+2)
		!Anterpolation, w will be mutiplied at the highest sampling points
		do j = 1, 2
			A_i_1%a_nm%item(j)%item(:) = TL(:, j)*B_i_1%a_nm%item(j)%item(:)*TRF_p(:)! 
			A_i_1%b_nm%item(j)%item(:) = TL(:, j)*B_i_1%b_nm%item(j)%item(:)*TRF_p(:)!
			A_i_1%a_nm%item(j+2)%item(:) = TL(:, j)*B_i_1%a_nm%item(j+2)%item(:)*TRF_p(:)! 
			A_i_1%b_nm%item(j+2)%item(:) = TL(:, j)*B_i_1%b_nm%item(j+2)%item(:)*TRF_p(:)!
		end do
		!print*, 'A_i_1%a_nm%item(1)%item(12)', A_i_1%a_nm%item(1)%item(12)
		!print*, 'TL(12, 1)', TL(12, 1)
		!print*, 'TL(12, 2)', TL(12, 2)
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
      use lib_tree_public
      use ml_fmm_type
      use lib_sie_mlfmm_tri
      use calculation_mod_tri
      use lib_sie_mlfmm_tri
      implicit none
      ! dummy
      type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
      type(lib_tree_spatial_point), intent(in) :: x_1
      type(lib_tree_spatial_point), intent(in) :: x_2
      type(lib_ml_fmm_coefficient) :: A_i_2, A2_inter
                
      type(lib_tree_spatial_point) :: r
      real(kind = dp) :: tmp
		complex(kind = 8) :: tmp_k(2)
      integer :: Ns, K_m, K_n
      integer :: i, m, pp, j  
      integer(kind = 1) :: l
         
		l = A_i_1%uindex%l
		A_i_2%uindex%l = l + int(1, 1)
		  
      K_m = truncation_number(l)%n
		K_n = truncation_number(l+1)%n
	   
		!r = x_children - x_parent, according Gibson
		r = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1) !		
		Ns = K_n*K_n*2
			
		if (K_n .lt. K_m) then				
			pp = 3
			Ns = K_n*K_n*2
			call init_list_sie(A_i_2, 4, Ns+2)
			call Khat_and_TRF_ExtendedK(K_n)
			call Lagrange_Interpolation_New(pp, K_m, K_n, A_i_1, A2_inter)
			do i = 1, Ns+1
				tmp = dot_product(K_hat_p(i, :), r%x)
				tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
				do j = 1, 2
					A_i_2%a_nm%item(j)%item(i) = tmp_k(j)*A2_inter%a_nm%item(j)%item(i)
					A_i_2%b_nm%item(j)%item(i) = tmp_k(j)*A2_inter%b_nm%item(j)%item(i)
					A_i_2%a_nm%item(j+2)%item(i) = tmp_k(j)*A2_inter%a_nm%item(j+2)%item(i)
					A_i_2%b_nm%item(j+2)%item(i) = tmp_k(j)*A2_inter%b_nm%item(j+2)%item(i)
				end do
			end do
		else
			Ns = K_m*K_m*2
			call init_list_sie(A_i_2, 4, Ns+2)
			do i = 1, Ns
				tmp = dot_product(K_hat(i, :), r%x)
				tmp_k = (/exp(-im*k1*tmp), exp(-im*k2*tmp)/)
				do j = 1, 2
					A_i_2%a_nm%item(j)%item(i) = tmp_k(j)*A_i_1%a_nm%item(j)%item(i)
					A_i_2%b_nm%item(j)%item(i) = tmp_k(j)*A_i_1%b_nm%item(j)%item(i)
					A_i_2%a_nm%item(j+2)%item(i) = tmp_k(j)*A_i_1%a_nm%item(j+2)%item(i)
					A_i_2%b_nm%item(j+2)%item(i) = tmp_k(j)*A_i_1%b_nm%item(j+2)%item(i)						
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
			use lib_tree_public
			use ml_fmm_type
			use lib_sie_mlfmm_tri            
			use lib_ml_fmm_data_container
			use lib_sie_data_container
			use calculation_mod_tri
            
			implicit none
			! dummy
			type(lib_tree_data_element),  intent(inout) :: data_element_y_j
			integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_number_j
			type(lib_tree_data_element), dimension(:), allocatable, intent(inout) :: data_element_e2
			integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2
			type(lib_ml_fmm_coefficient), intent(in) :: D
            
			type(lib_ml_fmm_coefficient) :: buffer_D
			type(lib_ml_fmm_v) :: rv, buffer_x
			type(lib_ml_fmm_v) :: rv_2, rv_tmp
            
			type(lib_tree_spatial_point) :: data_element_xj_unscaled
            
			integer :: j, i, m
			integer(kind = 1) :: l
			integer(kind=4) :: element_number_i
            
			double precision :: d_c
			complex(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: sum_near(2), sum_far(2)
			type(lib_ml_fmm_coefficient) :: ff_rl, ff_rk!
			integer :: K_m, number_j
			integer :: Ns
			real(kind = dp), dimension(:), allocatable :: TRF_II
			real(kind = dp), dimension(:, :), allocatable :: k_hat_dummy
                        
			l = data_element_y_j%uindex%l			
			j = l !must be changed
			d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
			
			K_m = truncation_number(l)%n
			Ns = 2*K_m*K_m
			
			if (element_number_j .eq. 1) then
			
			allocate(rv_tmp%C(2))
			allocate(rv%C(2))
			
			!Far field summation. f_rl, f_rk are receive functions for L and K parts, as is in .
			buffer_x = m_ml_fmm_u(element_number_j) !i or j?
			data_element_xj_unscaled = lib_tree_get_unscaled_point(data_element_y_j%point_x)
					
			call init_list_sie(ff_rl, 2, Ns)			
			do i = 1, 2
				ff_rl%a_nm%item(i)%item = preCalculated_RT(element_number_j)%a_nm%item(i)%item
				ff_rl%b_nm%item(i)%item = preCalculated_RT(element_number_j)%b_nm%item(i)%item
			end do
			call MLFMM_tri_summation_New(K_m, ff_rl, D, sum_far)
			
			!Near field summation
			sum_near(1:2) = (0.0, 0.0)
			do i = 1, size(data_element_e2)
				element_number_i = element_number_e2(i)
				rv_tmp = lib_sie_ml_fmm_get_u_phi_i_j(element_number_j, element_number_i) 
				sum_near(1:2) = sum_near(1:2) + rv_tmp%C(1:2)
			end do
			
			!Final summation
			rv%C(1) = sum_near(1) + sum_far(1)
			rv%C(2) = sum_near(2) + sum_far(2)
			if (element_number_j .le. 4) then
				print*, 'element_number_j', element_number_j
				print*,  'Far rv%C(1) =', sum_far(1) 
				print*,  'Far rv%C(2) =', sum_far(2) 
			end if
		end if		
		end function
      
		! Near field matrix
		!++++++++++++++++
		function lib_sie_ml_fmm_get_u_phi_i_j(element_number_j, &
																		element_number_i) result(rv)
			use lib_tree_public
			use ml_fmm_type
			use lib_ml_fmm_data_container
			use lib_sie_data_container
			use calculation_mod_tri
			implicit none
			! dummy
          
			integer(kind=4), intent(in) :: element_number_i          
			integer(kind=4), intent(in) :: element_number_j
			type(lib_ml_fmm_v) :: buffer_x, rv
			! auxiliary
			integer :: p, q, s, t, ngp
			real(kind = dp) :: r_ob, r_ref
			complex(kind = dp) :: ssum_a, ssum_b, ssum_c, sum_a, sum_b, sum_c
			
			allocate(buffer_x%C(2), rv%C(2))            
			buffer_x = m_ml_fmm_u(element_number_i)            
			p = struc%neighbours(element_number_i)%elements(1)
			q = struc%neighbours(element_number_i)%elements(2)
			r_ref = vec_len(struc%midpoint(p)%point - struc%midpoint(q)%point)	
			sum_a = (0.0, 0.0)
			sum_b = sum_a
			sum_c = sum_a
			do t = 1, 2
				do s = 1, 2!
					p = struc%neighbours(element_number_i)%elements(t)					
					q = struc%neighbours(element_number_j)%elements(s)						
					r_ob = vec_len(struc%midpoint(p)%point - struc%midpoint(q)%point)					
					if (r_ob > r_ref) then
						ngp = 3				
						call Normal_Integration_new(element_number_j, element_number_i, t, s, ngp, struc, ssum_a, ssum_b, ssum_c)
					else
						ngp = 3	
						call Singular_Integration(element_number_j, element_number_i, t, s, ngp, struc, ssum_a, ssum_b, ssum_c)
					end if
					sum_a = sum_a + ssum_a !
					sum_b = sum_b + ssum_b
					sum_c = sum_c + ssum_c
				end do !loop s 
			end do  !loop t
			rv%C(1) = ssum_b*buffer_x%C(1) - ssum_a*buffer_x%C(2)
			rv%C(2) = ssum_a*buffer_x%C(1) + ssum_c*buffer_x%C(2)
			!print*, ' rv%C(1)', rv%C(1)
			!print*, ' rv%C(2)', rv%C(2)
		end function
  !  !+++++++++++++++++
                                                 
		function get_Km(uindex, d_c, kd, d0) result(k_m)	 
			type(lib_tree_universal_index), intent(in) :: uindex    
			integer(kind = 1) :: l
			complex(kind = 8), intent(in) :: kd
			integer :: k_m, j
			real(kind = dp) :: d_c, d_tmp !diagonal of the box at level l
			integer, intent(in) :: d0		  
			l = uindex%l
			j = l !must be changed        
			d_tmp = abs(d_c/2*kd)
			k_m = d_tmp+1.8*(d0)**(2/3)*(d_tmp)**(1/3)
		 end function
	
	end module lib_sie_ml_fmm_interface

