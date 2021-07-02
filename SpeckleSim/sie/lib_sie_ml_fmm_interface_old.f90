    module lib_sie_ml_fmm_interface
    use libtree
    use libmlfmm
    use disc_mod
    use input_tri
    use lib_sie_data_container
    
    
    implicit none

    private

    public :: lib_sie_constructor
    
    !public :: lib_sie_destructor
    !
    public :: lib_sie_ml_fmm_calculate_vector_b
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles

    type(lib_ml_fmm_hierarchy), dimension(:), allocatable :: m_ml_fmm_hierarchy

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
    !        !
    subroutine lib_sie_constructor(data_elements, s_opt)
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
    !
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
        
        !print*, 'size vector_x', size(vector_x)
        do i = 1, m_pairs
            allocate(x(i)%C(2))
            x(i)%C(1)=vector_x(i)
            x(i)%C(2)=vector_x(i + m_pairs)
        end do
        
        call lib_ml_fmm_run(x, b)
        
        do i = 1, size(vector_x)/2
            error_b(i) = (b(i)%C(1) - vector_b(i))/vector_b(i)
            !print*, 'vector_b =', vector_b(i)
            !print*, 'b =', b(i)%C(1)
        end do
        
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
   
    function lib_sie_ml_fmm_get_c(x_c, data_element, element_number) result(C_box)
        use lib_tree_public
        use ml_fmm_type
        use lib_ml_fmm
        use lib_sie_mlfmm_tri
        use lib_sie_data_container
        use lib_ml_fmm_data_container
  
        implicit none
        ! dummy
        type(lib_tree_spatial_point), intent(in) :: x_c
        type(lib_tree_data_element), dimension(:), allocatable, intent(in) :: data_element !
        integer(kind=4), dimension(:), intent(in) :: element_number
        type(lib_ml_fmm_coefficient) :: C_box
        
        ! auxiliary
        integer :: i, K_m, m_element_number
        logical :: interpol
  
        type(lib_tree_spatial_point) :: data_element_x_unscaled
        complex(kind = dp), dimension(:, :), allocatable :: ff_out, ff_tmp, ff_inter
        integer :: Ns, m, n, j
        integer(kind = 1) :: l
        double precision :: d_c
        
        K_m = get_Km(data_element(1)%uindex) 
        Ns = K_m*K_m*2+2
        allocate(ff_tmp(Ns, 4))
        ff_tmp(1:Ns, 1:4) = (0.0, 0.0)
        call init_list_sie(C_box, 4, Ns)                
        call Khat_and_TRF_ExtendedK(K_m) 
        
        data_element_x_unscaled = lib_tree_get_unscaled_point(x_c)        
        do m = 1, 4
            do i = 1, size(element_number)
                m_element_number = element_number(i)
                call FMM_aggregation_coefficient_element(K_m, m_element_number, struc, data_element_x_unscaled%x, ff_out)
                ff_tmp(:, m) = ff_tmp(:, m) + ff_out(:, m)*m_ml_fmm_u(m_element_number)%C(1)
            end do!        
            C_box%a_nm%item(m)%item(:) = ff_tmp(:, m)
        end do
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
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: C_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2        
        type(lib_ml_fmm_coefficient) :: B_i_2
       
        type(lib_tree_spatial_point) :: r
        real(kind = dp), dimension(:, :), allocatable :: k_hat        
        real(kind = dp), dimension(:), allocatable :: TRF_dummy
        real(kind = dp) :: tmp
        integer :: Ns, K_m, K_D, Ns_p, K_n
        integer :: i, m, pp, j
        integer(kind = 1) :: l
        
        complex(kind = dp), dimension(:, :), allocatable :: ff_inter
		complex(kind = dp), dimension(:, :), allocatable :: fn_scs
        
        K_m = get_km(C_1%uindex)
        B_i_2%uindex%l = C_1%uindex%l - int(1, 1)
        K_n = get_km(B_i_2%uindex)
        
        Ns = K_m*K_m*2
        K_D = K_n - K_m          
        pp = 3
        Ns_p = 2*K_n*K_n
        
        allocate(K_hat(Ns, 3), fn_scs(Ns+2, 4))  !k_hat_p should be used
        call init_list_sie(B_i_2, 4, Ns_p)
        
        do m = 1, 2
            fn_scs(1:Ns, m) = C_1%a_nm%item(m)%item(1:Ns)
            fn_scs(1:Ns, m+2) = C_1%a_nm%item(m+2)%item(1:Ns)
        end do 
        fn_scs(Ns+1, :) = fxy_inter(1, :)
        fn_scs(Ns+2, :) = fxy_inter(2, :)
        r = lib_tree_get_unscaled_point(x_2)- lib_tree_get_unscaled_point(x_1) 
        
        call Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter)      
        !print*, 'ff_inter in translation_SS 1 =', ff_inter(10, 1)
        !print*, 'ff_inter in translation_SS 2 =', ff_inter(10, 2)
        !print*, 'ff_inter in translation_SS 3 =', ff_inter(10, 3)
        
        call Khat_and_TRF_II(K_n, K_hat, TRF_dummy) 
        
        do i = 1, Ns
            tmp = dot_product(K_hat(i, :), r%x)
            do m = 1, 2
                B_i_2%a_nm%item(m)%item(i) = exp(im*k1*tmp)*ff_inter(i, m)
                B_i_2%a_nm%item(m+2)%item(i) = exp(im*k2*tmp)*ff_inter(i, m+2)
            end do            
        end do   
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
        type(lib_ml_fmm_coefficient) :: A_i_1
    
        double precision :: d_c, r_tmp	!problem with it, how can one obtain in a fast way?
        integer :: Ns, Lm, K_m, K_n, i, j
        integer(kind = 1) :: l
        type(lib_tree_spatial_point) :: r_ab
        complex(kind = dp), dimension(:), allocatable :: TL_1, TL_2
        intrinsic :: sqrt
        
        l = B_i_1%uindex%l
        j = l
        d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
        
        K_m = get_km(B_i_1%uindex)
        Ns = K_m*K_m*2
        call init_list_sie(A_i_1, 4, Ns)        
        allocate(TL_1(Ns), TL_2(Ns))
        r_ab = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1) !x_2 target?!
        r_tmp = vec_len(r_ab%x)
        !Lm = abs(k1*r_tmp) + 10*(abs(k1*r_tmp))**(1/3)
        
        TL_1 = TL_km_arr_III(r_ab%x, K_m, k1)
		TL_2 = TL_km_arr_III(r_ab%x, K_m, k2)                
        do i = 1, 2
            A_i_1%a_nm%item(i)%item(:) = TL_1(:)*B_i_1%a_nm%item(i)%item(:)! 
            A_i_1%a_nm%item(i+2)%item(:)= TL_2(:)*B_i_1%a_nm%item(i+2)%item(:)
        end do
        
        !print*, '   '        
        !print*, 'translation SR:', A_i_1%a_nm%item(1)%item(2) !TL_1(5)!
        !print*, 'translation SR:', A_i_1%a_nm%item(4)%item(2) !TL_2(5)!
        
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
    function lib_sie_ml_fmm_translation_RR(A_i_2, x_1, x_2) result(A_i_1)
        use lib_tree_public
        use ml_fmm_type
        use lib_sie_mlfmm_tri
        use calculation_mod_tri
        use lib_sie_mlfmm_tri
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: A_i_2
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: A_i_1
        
        complex(kind = dp), dimension(:,:), allocatable :: fn_scs, ff_inter
        real(kind = dp), dimension(:), allocatable :: TRF_II
        real(kind = dp), dimension(:, :), allocatable :: k_hat
        type(lib_tree_spatial_point) :: r
        real(kind = dp) :: tmp
        integer :: Ns, K_m, K_n, Ns_p, K_D !The sampling point should be different in different levels
        integer :: i, m, pp
        integer(kind = 1) :: l
        
        K_m = get_km(A_i_2%uindex)       
        A_i_1%uindex%l = A_i_2%uindex%l + int(1, 1)
        K_n = get_km(A_i_1%uindex)
        K_D = K_n - K_m          
        
        Ns = K_m*K_m*2
        Ns_p = 2*K_n*K_n
        pp = 3
        
        allocate(K_hat(Ns, 3), fn_scs(Ns+2, 4))  !k_hat_p should be used
        call init_list_sie(A_i_1, 4, Ns_p)        
        do m = 1, 2
            fn_scs(1:Ns, m) = A_i_2%a_nm%item(m)%item(1:Ns)
            fn_scs(1:Ns, m+2) = A_i_2%a_nm%item(m+2)%item(1:Ns)
        end do 
        fn_scs(Ns+1, :) = fxy_inter(1, :)
        fn_scs(Ns+2, :) = fxy_inter(2, :)
        
        call Lagrange_Interpolated(pp, K_m, K_D, Fn_scs, ff_inter)      
        
        call Khat_and_TRF_II(K_n, K_hat, TRF_II)
       
        r = lib_tree_get_unscaled_point(x_2)-lib_tree_get_unscaled_point(x_1)        
        do i = 1, Ns_p
            tmp = dot_product(K_hat(i, :), r%x)
            do m = 1, 2
                A_i_1%a_nm%item(m)%item(i) = exp(im*k1*tmp)*A_i_2%a_nm%item(m)%item(i)
                A_i_1%a_nm%item(m+2)%item(i) = exp(im*k2*tmp)*A_i_2%a_nm%item(m+2)%item(i)
            end do
        end do
        
    end function
    
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
            complex(kind = dp), dimension(:, :), allocatable :: ff_rl, ff_rk, ff_D, D_abc
            complex(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: sum_rr(3), sum_r(2)
            
            integer :: K_m, number_j
		    integer :: Ns
            real(kind = dp), dimension(:), allocatable :: TRF_II
		    real(kind = dp), dimension(:, :), allocatable :: k_hat_dummy
            
                        
            l = data_element_y_j%uindex%l
            !print*, 'size of data_element_y_j', size(data_element_y_j)
            
            j = l !must be changed
            d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1) 
            K_m = get_Km(data_element_y_j%uindex) 
            !print*, 'K_m=', K_m
            Ns = 2*K_m*K_m
          
            allocate(ff_D(Ns, 4))     
            allocate(rv_tmp%C(2))
            allocate(rv%C(2))
            allocate(buffer_x%C(2))
            
            buffer_x = m_ml_fmm_u(element_number_j)
            call Khat_and_TRF_II(K_m, K_hat_dummy, TRF)            
            do m = 1, 4
                ff_D(:, m) = D%a_nm%item(m)%item(:)
            end do
            
            data_element_xj_unscaled = lib_tree_get_unscaled_point(data_element_y_j%point_x)               
            call FMM_disaggregation_coefficient_RL(K_m, element_number_j, struc, data_element_xj_unscaled%x, d_c, ff_rl) !
            call FMM_disaggregation_coefficient_RK(K_m, element_number_j, struc, data_element_xj_unscaled%x, d_c, ff_rk)
            
            !print*, 'ff_rl', ff_rl(10, 1)
            !print*, 'ff_rk', ff_rk(10, 1)
            
            sum_r(1:2) = (0.0, 0.0)
            do i = 1, size(data_element_e2)
                element_number_i = element_number_e2(i)
                rv_tmp = lib_sie_ml_fmm_get_u_phi_i_j(element_number_i, element_number_j)
                sum_r(1:2) = sum_r(1:2) + rv_tmp%C(1:2)
            end do
            !
            call  MLFMM_tri_Dabc_calculation(ff_rk, ff_rl, ff_D, D_abc)
            do i = 1, 3
                sum_rr(i) = sum(D_abc(i, :)*TRF(:))                    
            end do
            sum_r(1) = sum_r(1) + sum_rr(2)*buffer_x%C(1) - sum_rr(1)*buffer_x%C(2)
            sum_r(2) = sum_r(2) + sum_rr(1)*buffer_x%C(1) + sum_rr(3)*buffer_x%C(2)
            !print*, 'sum_r', sum_r(1)
            !print*, 'sum_r', sum_r(2)
        end function
        
        !++++++++++++++++
        function lib_sie_ml_fmm_get_u_phi_i_j(element_number_i, &
                                                 element_number_j) result(rv)
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
			r_ref = vec_len(struc%midpoints(p)%point - struc%midpoints(q)%point)	
            sum_a = (0.0, 0.0)
            sum_b = sum_a
            sum_c = sum_a
            do t = 1, 2
				do s = 1, 2!
					p = struc%neighbours(element_number_i)%elements(t)					
					q = struc%neighbours(element_number_j)%elements(s)						
					r_ob = vec_len(struc%midpoints(p)%point - struc%midpoints(q)%point)					
					if (r_ob > r_ref) then
						ngp = 3				
						call Normal_Integration_new(element_number_i, element_number_j, t, s, ngp, struc, ssum_a, ssum_b, ssum_c)								
					else
						ngp = 3	
						call Singular_Integration(element_number_i, element_number_j, t, s, ngp, struc, ssum_a, ssum_b, ssum_c)
					end if
					sum_a = sum_a + ssum_a !
					sum_b = sum_b + ssum_b
					sum_c = sum_c + ssum_c
				end do !loop s 
            end do  !loop t
            
            rv%C(1) = ssum_b*buffer_x%C(1) - ssum_a*buffer_x%C(2)
            rv%C(2) = ssum_a*buffer_x%C(1) + ssum_c*buffer_x%C(2)
    end function
  !  !+++++++++++++++++
                                                 
    function get_Km(uindex) result(K_m)
        use lib_tree_public
        use ml_fmm_type
        use lib_ml_fmm_data_container
        use lib_sie_data_container       

        type(lib_tree_universal_index), intent(in) :: uindex    
        integer(kind = 1) :: l
        integer :: K_m, j
        real(kind = dp) :: d_c
        l = uindex%l
        j = l !must be changed
        d_c = lib_tree_get_box_diagonal(j)*lib_tree_scaling_D%x(1)
        k_m = abs(k2*d_c)+1.8*(d0)**(2/3)*(abs(k2)*d_c)**(1/3)
    end function
    
    end module lib_sie_ml_fmm_interface

