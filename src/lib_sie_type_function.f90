module lib_sie_type_function
	use lib_sie_constants
	use lib_sie_lookup_table	
	use libmath	
	
	implicit none
	
	type :: point
      real(dp) :: point(3)
	end type point
		
	type :: fn_rwg
		integer :: sign
		type(point), dimension(3) :: corners
		real(dp) :: length !
		real(dp) :: area
		real(dp) :: df
		real(dp) :: mid(3)
		real(dp) :: nq(3)
	end type fn_rwg
	
	type :: edges_fn
		type(fn_rwg) :: fn
	end type edges_fn
	
	type :: evaluation_r_media
		real(dp) :: point(3)
		complex(dp) :: eps_r
		integer :: media !1 or 2: 1 for free space and 2 for medium
		integer :: ngp
	end type
	
	type :: edges_in_structure
		type(edges_fn), dimension(2) :: fn_ot
	end type edges_in_structure
	
	type :: edge_tri
		real(dp) :: p_s(3) !start point
		real(dp) :: p_t(3) !terminate point
	end type edge_tri
	
	type :: vector_r
		real(dp) :: vector(3)
	end type vector_r

	type :: array_c
		complex, dimension(:), allocatable :: Array
	end type array_c

	type :: matrix_c
		complex, dimension(10,10) :: matrix
	end type matrix_c
	
	!******************++
	type :: Element_tri
		integer, dimension(:), allocatable :: corners !index of each corner point				
	end type Element_tri

	
   type :: Pair
		integer :: corner(2) !corner indices
		integer :: element(2) !element indices		
   end type Pair
	
	type csr_format
		complex(dp), dimension(:), allocatable :: a
		complex(dp), dimension(:), allocatable :: b
		complex(dp), dimension(:), allocatable :: c
		complex(dp), dimension(:), allocatable :: d
		integer, dimension(:), allocatable :: ia
		integer, dimension(:), allocatable :: ja
	end type csr_format

	type preset_type
		character(len = 15) :: illumination
		character(len = 15) :: formulation
		character(len = 15) :: evaluation
		character(len = 15) :: object
		character(len = 15) :: calculation	
	end type preset_type
	
   type :: Vector_c
		complex(dp) :: vector(3) !changed into point
	end type Vector_c
	
	type structure_tri
      type(Point), dimension(:), allocatable :: points
      type(Element_tri), dimension(:), allocatable :: elements
      type(Pair), dimension(:), allocatable :: neighbours !should be changed into pairs
      type(Point), dimension(:), allocatable :: midpoint !midpoint of element
		type(Point), dimension(:), allocatable :: midpoint_edge
		type(Point), dimension(:), allocatable :: x_c !box center at level lmax
   end type structure_tri
	
	type tri_surface_parameter
		real(dp), dimension(3) :: D !dimension of the surface, D(1:3) = (/Dx, Dy, Dz/)
		type(point) :: centroid		
		complex(dp) :: eps_r2
		complex(dp) :: eps_r1
	end type tri_surface_parameter
	
	type tri_sphere_parameter
		real(dp) :: D !sphere diameter for sphere_type
		integer :: n_disc !discretization number
		type(point) :: centroid
		complex(dp) :: eps_r2
		complex(dp) :: eps_r1
	end type tri_sphere_parameter
	
	type quad_surface_parameter
		real(dp), dimension(3) :: D !dimension of the surface, D(1:3) = (/Dx, Dy, Dz/)
		integer :: np  !element number for surface
		integer :: nt  !element number for surface
		type(point) :: centroid		
		complex(dp) :: eps_r2
		complex(dp) :: eps_r1
	end type quad_surface_parameter
	
	type quad_sphere_parameter
		real(dp) :: D
		integer :: np
		integer :: nt
		type(point) :: centroid
		complex(dp) :: eps_r2
		complex(dp) :: eps_r1
	end type quad_sphere_parameter
	
			
	type edges_parameter
		integer, dimension(5) :: edge_p
	end type
	
	type element_edges
		type(edges_parameter), dimension(:), allocatable :: nr_edges
	end type
	
	type truncation_hl
		integer :: n
		logical :: near_field		
	end type
	
	type box_size_at_l
		real(dp) :: dc
		integer(kind = 1) :: lmax
	end type
	
	type lib_sie_illumination_parameter
		real(dp) :: lambda
		real(dp) :: theta_in
		real(dp) :: phi_in        
        real(dp) :: phi_max
		real(dp) :: k_in(3)
		real(dp) :: E_in(3) !amplitude
		integer :: Nt ! for conical illumination, sampling points along theta
		integer :: Np ! for conical illumination, sampling points along phi
		character(len = 5) :: pol
	end type lib_sie_illumination_parameter
	
	type w_matrix
		real(dp), dimension(:, :), allocatable :: w_tp
	end type w_matrix
	
	type Lagrange_inter
		type(w_matrix), dimension(:), allocatable :: w_coeff
	end type Lagrange_inter
	
	!i4point - kd - khat
	type lib_sie_mlfmm_TL_coefficient	
      type(list_list_list_cmplx), dimension(:), allocatable :: TL ! dimension reserved for level 
		type(LT3_t), dimension(:), allocatable :: my_lt					! dimension reserved for level 
   end type 
	
	type result_inner_JMCFIE
		complex(dp) :: IK(3)
		complex(dp) :: IT_a(3)
		complex(dp) :: IT_bn(3)
		complex(dp) :: IT_bt
	end type
	
	type result_element_JMCFIE
		complex(dp) :: Kn
		complex(dp) :: Kt
		complex(dp) :: Tn
		complex(dp) :: Tt
	end type

	type lib_sie_parameter_objective
		real(dp) :: M !magnification
		real(dp) :: NA
		real(dp) :: r_ph !radius of the pinhole		
	end type lib_sie_parameter_objective
	!********
	!For quadrilateral elements
	type :: Edge_Adjacent !Pair is changed into "Edge_Adjacent"
		integer :: CommonEdge_Nodes(3) ! The nodes of the overlapped edge
		integer :: CommonEdge_Elements(8) ! indices of the two elements and their node indices
	end type Edge_Adjacent 
    
	type :: Vertex_Adjacent!Pair_point !used for common corner
		integer :: CommonCorner(1)
		integer :: CommonCorner_Elements(6) !indices of the two elements 
	end type Vertex_Adjacent
	
	type :: Element_quad
		integer, dimension(:), allocatable :: vertices  !8 vertex index of each element 
		type(Edge_Adjacent), dimension(:), allocatable :: edge_overlap !
		type(Vertex_Adjacent), dimension(:), allocatable :: corner_overlap !	
	end type Element_quad
	
	type structure_quad
		type(Point), dimension(:), allocatable :: node_vector!
		type(Element_quad), dimension(:), allocatable :: elements	
		type(Point), dimension(:), allocatable :: midpoint !midpoint of element
		type(Point), dimension(:), allocatable :: x_c !box center at level lmax
	end type structure_quad
	
	type lib_sie_evaluation_point_type
		type(point) :: coordinate
		type(vector_c) :: e_field
		type(vector_c) :: h_field
		integer :: ngp
	end type

	 type lib_sie_evaluation_parameter_type			
		real(dp) :: dim_a(2) !min and max dimension range
		real(dp) :: dim_b(2)	!min and max dimension range		
		integer :: N_dim(2) !sampling points in the two dimensions
		integer(kind = 1) :: tot_field
		integer(kind = 1) :: ngp
		real(dp) :: dim_c !position of the third dimension
	 end type
	
	type lib_sie_simulation_parameter_type
        !double precision :: refractive_index_medium
        !type(lib_sie_illumination_parameter) :: illumination
        complex(dp), dimension(:), allocatable :: vector_x
    end type lib_sie_simulation_parameter_type	 	 
	 
end module
	