program main
    
	use input_tri	   

	use calculation_mod_tri
	use disc_mod	
	use time_stamp
	use omp_lib
	use libmath
	use libtree
	use libmlfmm
    use lib_sie_ml_fmm_interface
    use lib_sie_data_container
    
	implicit none
    
    
	
	integer :: n_el, m, n, mr
	integer :: ngp, m_cube, nn_edge, n_cube, nd, n_cube_number, n_edge_number, M_edge
	
	real(kind = dp) :: D_cube
	double precision, dimension(TREE_DIMENSIONS) :: pos_edge
	!type(Structure) :: struct	
	
	type(lib_ml_fmm_data) :: data_elements
	
	type(Structure_cube) :: struct_cube
	complex(kind = dp), dimension(:), allocatable :: b1, x_estimate, b1_tmp
	complex(kind = dp) :: sum_EH(2), sum_a, sum_b, sum_c
	complex(kind = dp), dimension(:,:), allocatable :: D_mat, Dmat_f, Dmat_n
	integer, dimension(:), allocatable :: ipiv
	integer :: info, ot
	
	!++++++++++++++++++++++++++++	
	type(ml_fmm_type_operator_procedures) :: operator_procedures
	type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures
    integer(kind=4) :: tree_s_opt
 
    logical :: final_sum_calc_y_hierarchy
    logical :: final_sum_calc_xy_hierarchy
	
    type(lib_tree_universal_index) :: uindex
	logical :: ignore_box
    type(lib_ml_fmm_coefficient) :: C_i
    
    double complex, dimension(:), allocatable :: vector_x 
    double complex, dimension(:), allocatable :: vector_b
    double precision, dimension(:, :), allocatable :: x_tmp
    intrinsic :: cmplx
	!+++++++++++
	
    call Data_Import()	
    n_el = size(struc%elements)
	print*, 'n_el=', n_el
	m_pairs = size(struc%neighbours)	
	print*, 'm_pairs=', m_pairs
	allocate(data_elements%xy(m_pairs))
    allocate(vector_x(2*m_pairs), x_tmp(2*m_pairs, 2))
    allocate(vector_b(2*m_pairs)) 
    
    ! calculate the rhs
    do m = 1, m_pairs
		ngp = 6		
		if (Illumination == 'Gauss' .or. Illumination == 'gauss') then            
			call incident_Gauss_tri(struc, m, ngp, Sum_EH)
		else if ((illumination == 'Plane') .or. (illumination == 'plane')) then
			call incident_plane_tri(struc, m, ngp, Sum_EH) !calculate separately		
		else if ((illumination == 'FB') .or. (illumination == 'fb')) then    
			call incident_FocusedGauss_tri(struc, m, ngp, Sum_EH) !calculate separately	
      else 
			print*, 'Not a proper illumination type'            
			call exit
      end if
		vector_b(m) = Sum_EH(1)
		vector_b(m + m_pairs) = Sum_EH(2)	
		vector_x(m) = (1.0, 0.0)
		vector_x(m + m_pairs) = (1.0, 0.0)!Initial guess
    end do
    !+++++++ vector_x for test+++++++++++++
    
    open(unit = 12, file = 'Result_SEH_tri.txt')
        do m = 1, m_pairs*2
            read(12, *) x_tmp(m, :)
    end do
    close(12)
    
    vector_x(:) = cmplx(x_tmp(:, 1), x_tmp(:, 2))
    print*, 'vector_x(12)=', vector_x(12)
    
    !+++++++++++++++    
    ot = 1
	do m = 1, m_pairs
		call fn_edge_center(struc, m, ot, pos_edge)
		data_elements%xy(m)%point_x%x=pos_edge
        data_elements%XY(m)%element_type = 1 ! value .ne. -1
        data_elements%XY(m)%hierarchy = HIERARCHY_XY
    end do
	tree_s_opt = 4
    
    call lib_sie_constructor(data_elements, tree_s_opt)
    call lib_sie_ml_fmm_calculate_vector_b(vector_x, vector_b)
   
	!call lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box) result(C_i)
	!lib_ml_fmm_get_C_i_from_elements_at_box!
    
    !subroutine test(a)
    !    integer, dimension(:), allocatable, intent(inout)  :: a
    !    integer :: b(2), c(3)
    !    b(1) = 2
    !    b(2) = 1
    !    allocate(a, source = b)
    !end 
end 