module lib_sie_data_container
    use ml_fmm_type
    
    implicit none
    
   	integer, parameter :: ng = 3		
    integer, parameter :: d0 = 4 !corret number of digits
	
	double precision, dimension(:, :), allocatable :: k_hat
	double precision, dimension(:), allocatable :: TRF
    integer :: m_pairs	!number of edge elements
    
	type list_list_cmplx_vector
        type(list_list_cmplx), dimension(3) :: vector
    end type
    
    double precision, dimension(:,:), allocatable :: k_hat_p
    double precision, dimension(:,:), allocatable :: transfer_theta
    double precision, dimension(:,:), allocatable :: transfer_phi
    double precision, dimension(:), allocatable :: TRF_p
    double precision, dimension(2, 4) :: fxy_inter
    
end module lib_sie_data_container