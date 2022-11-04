module lib_ml_fmm_type_operator
    implicit none

    private

    ! ----- public functions -----
    public :: lib_ml_fmm_type_operator_constructor

    public :: operator (*)
    public :: operator (+)
    public :: operator (-)
    public :: operator (.CoR.)

    public :: operator (.eq.)
    public :: operator (.ne.)

    public :: lib_ml_fmm_type_operator_set_coefficient_zero
!    public :: lib_ml_fmm_type_operator_allocate_coefficient_list
!    public :: lib_ml_fmm_type_operator_deallocate_coefficient_list
    public :: lib_ml_fmm_type_operator_set_coefficient
    public :: lib_ml_fmm_type_operator_get_coefficient

    public :: lib_ml_fmm_get_B_i
    public :: lib_ml_fmm_phi_i_j
    public :: lib_ml_fmm_translation_RR
    public :: lib_ml_fmm_translation_SR
    public :: lib_ml_fmm_translation_SS

    ! ---- public type definitions -----
    public :: ml_fmm_type_operator_procedures
    public :: lib_ml_fmm_procedure_handles
!    public :: ml_fmm_coefficient_add_operator

    ! ----- operator -----
    interface operator (+)
        module procedure lib_ml_fmm_m_coefficient_add
        !module procedure m_coefficient_add
!        module procedure m_procedures%m_coefficient_add
        module procedure lib_ml_fmm_v_operator_add
        module procedure lib_ml_fmm_v_operator_add_0d
    end interface

    interface operator (-)
        module procedure lib_ml_fmm_v_operator_sub
        module procedure lib_ml_fmm_v_operator_sub_0d
    end interface

    interface operator (*)
        module procedure m_u_dot_coefficient
    end interface

    interface operator (.eq.)
        module procedure m_coefficient_eq
    end interface

    interface operator (.ne.)
        module procedure m_coefficient_ne
    end interface

    ! Sum_q=0_p-1[C_q R_q(y_j − x_∗)] = C o R(y_j − x_∗)
    !                                     ^ cor-operator
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    interface operator (.CoR.)
        module procedure m_cor
    end interface

    ! ----- interfaces -----
    interface
        function ml_fmm_coefficient_add_operator(lhs,rhs) result (rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_coefficient_add_operator

        function ml_fmm_u_dot_coefficient_operator(lhs,rhs) result (rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_coefficient), intent(in) :: rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_u_dot_coefficient_operator

        function ml_fmm_cor_operator(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_tree_spatial_point), intent(in) :: rhs
            type(lib_ml_fmm_v) :: rv
        end function ml_fmm_cor_operator

        function ml_fmm_coefficient_eq(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_coefficient), intent(in) :: rhs
            logical :: rv
        end function

        function ml_fmm_coefficient_ne(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_coefficient), intent(in) :: rhs
            logical :: rv
        end function

        subroutine ml_fmm_coefficient_set_zero(coefficient)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
        end subroutine

!        subroutine ml_fmm_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!        end subroutine
!
!        subroutine ml_fmm_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            type(lib_ml_fmm_coefficient_list_list), intent(inout) :: coefficient_list
!        end subroutine

        subroutine ml_fmm_set_coefficient(coefficient, uindex, hierarchy)
            use ml_fmm_type
            use lib_tree_public
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_ml_fmm_coefficient), intent(in) :: coefficient
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
        end subroutine

        function ml_fmm_get_coefficient(uindex, hierarchy) result(coefficient)
            use ml_fmm_type
            use lib_tree_public
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_ml_fmm_coefficient) :: coefficient
        end function

!        ! Basis function: A
!        function lib_ml_fmm_get_A(x) result(A)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_A) :: A
!        end function

!        ! Basis function: A_i
!        function lib_ml_fmm_get_A_i(x, data_element) result(A_i)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_data_element) :: data_element
!            type(lib_ml_fmm_coefficient) :: A_i
!        end function

!        ! Basis function: B
!        function lib_ml_fmm_get_B(x) result(B)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_B) :: B
!        end function

        ! Basis function: B_i
        function lib_ml_fmm_get_B_i(x, data_element) result(B_i)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element) :: data_element
            type(lib_ml_fmm_coefficient) :: B_i
        end function

!        ! Basis function: S
!        function lib_ml_fmm_get_S(x) result(S)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_S) :: S
!        end function
!
!        ! Basis function: R
!        function lib_ml_fmm_get_R(x) result(R)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_R) :: R
!        end function

!        ! Local expansion (inner or Regular expansion)
!        !
!        ! Restriction
!        ! ----
!        ! Example calculation:
!        !       phi_i(y) = A_i (x_∗) o R(y − x_∗)     (8)
!        !
!        !   "Here the series is valid in the domain |y − x_∗| .le. r_c |x_i − x_∗| (see Fig. 2),
!        !   where 0 < r c < 1 is some real number."
!        !
!        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   i
!        !   x_*
!        !   y
!        !
!        ! Returns
!        ! ----
!        !   phi_i
!        !
!        !
!        function lib_ml_fmm_expansion_R(i, x, y) result(phi_i)
!            use lib_tree_type
!            implicit none
!            ! dummy
!            integer, intent(in) :: i
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_spatial_point), intent(in) :: y
!            integer :: phi_i                   ! todo: define type
!        end function lib_ml_fmm_expansion_R

!        ! Far field expansion (outer, Singular or multipole expansion)
!        !
!        ! Restriction
!        ! ----
!        ! Example calculation:
!        !   "Any function phi_i(y) has a complementary expansion valid outside
!        !   a d-dimensional sphere centered at y = x_∗ with radius R_c |x_i − x_∗| :
!        !       phi_i(y) = B_i (x_∗) o S(y − x_∗), |y - x_*| .ge. R_c |x_i - x_*|,     (10)
!        !
!        !   where R_c > 1 is a real number similar to r_c".
!        !   "Even though for many physical fields, such as the Green’s function
!        !   for Laplace’s equation, the function S(y − x_∗) is singular at y = x_∗ ,
!        !   this condition is not necessary. In particular we can have S = R."
!        !
!        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   i
!        !   x_*
!        !   y
!        !
!        ! Returns
!        ! ----
!        !   phi_i
!        !
!        !
!        function lib_ml_fmm_expansion_S(i, x, y) result(phi_i_j)
!            use lib_tree_type
!            implicit none
!            ! dummy
!            integer, intent(in) :: i
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_tree_spatial_point), intent(in) :: y
!            integer :: phi_i_j                   ! todo: define type
!        end function lib_ml_fmm_expansion_S

        function lib_ml_fmm_phi_i_j(data_element_i, y_j) result(rv)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_data_element), intent(inout) :: data_element_i
            type(lib_tree_spatial_point), intent(inout) :: y_j
            type(lib_ml_fmm_v) :: rv

        end function

        ! Translation: local-to-local (Regular-to-Regular)
        !
        ! Arguments
        ! ----
        !   A_i_1
        !   x_1
        !   x_2
        !
        ! Returns
        ! ----
        !   A_i_2
        !
        function lib_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
        end function

        ! Translation: far-to-local (Singular-to-Regular)
        !
        ! Arguments
        ! ----
        !   B_i_1
        !       set of expansion coefficients
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_1: spatial point
        !       origin of coordinate system 2
        !
        ! Returns
        ! ----
        !   A_i_2
        !       set of expansion coefficients
        !
        function lib_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
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
        function lib_ml_fmm_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
            use lib_tree_public
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: B_i_2
        end function
    end interface

    type lib_ml_fmm_procedure_handles
!        procedure(lib_ml_fmm_get_A_i), pointer, nopass :: get_A_i => null()
        procedure(lib_ml_fmm_get_B_i), pointer, nopass :: get_B_i => null()
!        procedure(lib_ml_fmm_get_S), pointer, nopass :: get_S => null()
!        procedure(lib_ml_fmm_get_R), pointer, nopass :: get_R => null()
!        procedure(lib_ml_fmm_expansion_R), pointer, nopass :: expansion_R => null()
!        procedure(lib_ml_fmm_expansion_S), pointer, nopass :: expansion_S => null()
        procedure(lib_ml_fmm_phi_i_j), pointer, nopass :: get_phi_i_j => null()
        procedure(lib_ml_fmm_translation_RR), pointer, nopass :: get_translation_RR => null()
        procedure(lib_ml_fmm_translation_SR), pointer, nopass :: get_translation_SR => null()
        procedure(lib_ml_fmm_translation_SS), pointer, nopass :: get_translation_SS => null()
    end type lib_ml_fmm_procedure_handles

    type ml_fmm_type_operator_procedures
        procedure(ml_fmm_coefficient_add_operator), pointer, nopass :: coefficient_add => null()
        procedure(ml_fmm_u_dot_coefficient_operator), pointer, nopass :: u_dot_coefficient => null()
        procedure(ml_fmm_cor_operator), pointer, nopass :: cor => null()
        procedure(ml_fmm_coefficient_set_zero), pointer, nopass :: coefficient_set_zero => null()
!        procedure(ml_fmm_allocate_coefficient_list), pointer, nopass :: allocate_coefficient_list => null()
        procedure(ml_fmm_set_coefficient), pointer, nopass :: set_coefficient => null()
        procedure(ml_fmm_get_coefficient), pointer, nopass :: get_coefficient => null()
!        procedure(ml_fmm_deallocate_coefficient_list), pointer, nopass :: deallocate_coefficient_list => null()
        procedure(ml_fmm_coefficient_eq), pointer, nopass :: coefficient_eq => null()
        procedure(ml_fmm_coefficient_eq), pointer, nopass :: coefficient_ne => null()
    end type

    ! ----- member procedures -----
    procedure(ml_fmm_coefficient_add_operator), pointer :: m_coefficient_add => null()
    procedure(ml_fmm_u_dot_coefficient_operator), pointer :: m_u_dot_coefficient => null()
    procedure(ml_fmm_cor_operator), pointer :: m_cor => null()
    procedure(ml_fmm_coefficient_set_zero), pointer :: m_coefficient_set_zero => null()
!    procedure(ml_fmm_allocate_coefficient_list), pointer :: m_allocate_coefficient_list => null()
    procedure(ml_fmm_set_coefficient), pointer :: m_set_coefficient => null()
    procedure(ml_fmm_get_coefficient), pointer :: m_get_coefficient => null()
!    procedure(ml_fmm_deallocate_coefficient_list), pointer :: m_deallocate_coefficient_list => null()
    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_eq => null()
    procedure(ml_fmm_coefficient_eq), pointer :: m_coefficient_ne => null()

    contains

    subroutine lib_ml_fmm_type_operator_constructor(operator_procedures)
!    subroutine lib_ml_fmm_type_operator_constructor(coefficient_add, u_dot_coefficient, cor, &
!                                                    set_coefficient_zero)!, &
!                                                    allocate_coefficient_list, &
!                                                    set_coefficient, get_coefficient, &
!                                                    deallocate_coefficient_list)
        use ml_fmm_type
        implicit none
        ! dummy
!        procedure(ml_fmm_coefficient_add_operator) :: coefficient_add
!        procedure(ml_fmm_u_dot_coefficient_operator) :: u_dot_coefficient
!        procedure(ml_fmm_cor_operator) :: cor
!        procedure(ml_fmm_coefficient_set_zero) :: set_coefficient_zero
!        procedure(ml_fmm_allocate_coefficient_list) :: allocate_coefficient_list
!        procedure(ml_fmm_set_coefficient) :: set_coefficient
!        procedure(ml_fmm_get_coefficient) :: get_coefficient
!        procedure(ml_fmm_deallocate_coefficient_list) :: deallocate_coefficient_list
!
!        m_coefficient_add => coefficient_add
!        m_u_dot_coefficient => u_dot_coefficient
!        m_cor => cor
!        m_coefficient_set_zero => set_coefficient_zero
!        m_allocate_coefficient_list => allocate_coefficient_list
!        m_set_coefficient => set_coefficient
!        m_get_coefficient => get_coefficient
!        m_deallocate_coefficient_list => deallocate_coefficient_list


        type(ml_fmm_type_operator_procedures) :: operator_procedures
        m_coefficient_add => operator_procedures%coefficient_add
        m_u_dot_coefficient => operator_procedures%u_dot_coefficient
        m_cor => operator_procedures%cor
        m_coefficient_set_zero => operator_procedures%coefficient_set_zero
!        m_allocate_coefficient_list => operator_procedures%allocate_coefficient_list
        m_set_coefficient => operator_procedures%set_coefficient
        m_get_coefficient => operator_procedures%get_coefficient
!        m_deallocate_coefficient_list => operator_procedures%deallocate_coefficient_list

        m_coefficient_eq => operator_procedures%coefficient_eq
        m_coefficient_ne => operator_procedures%coefficient_ne

    end subroutine lib_ml_fmm_type_operator_constructor

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (associated (m_coefficient_set_zero)) then
            call m_coefficient_set_zero(coefficient)
        else
            print *, "lib_ml_fmm_type_operator_set_coefficient_zero:  ERROR"
            print *, "  m_coefficient_set_zero is not associated"
        end if
    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero

!    subroutine lib_ml_fmm_type_operator_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_allocate_coefficient_list) ) then
!                call m_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            else
!                print *, "lib_ml_fmm_type_operator_allocate_coefficient_list:  ERROR"
!                print *, "  m_allocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_allocate_coefficient_list
!
!    subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_deallocate_coefficient_list) ) then
!                call m_deallocate_coefficient_list(coefficient_list)
!            else
!                print *, "lib_ml_fmm_type_operator_deallocate_coefficient_list:  ERROR"
!                print *, "  m_deallocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list

    subroutine lib_ml_fmm_type_operator_set_coefficient(coefficient, uindex, hierarchy)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(inout) :: uindex
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy

            if ( associated(m_set_coefficient)) then
                call m_set_coefficient(coefficient, uindex, hierarchy)
            else
                print *, "lib_ml_fmm_type_operator_set_coefficient:  ERROR"
                print *, "  m_set_coefficient is not associated"
            end if
        end subroutine

        function lib_ml_fmm_type_operator_get_coefficient(uindex, hierarchy) result(coefficient)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(inout) :: uindex
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_ml_fmm_coefficient) :: coefficient

            if ( associated(m_set_coefficient)) then
                coefficient = m_get_coefficient(uindex, hierarchy)
            else
                print *, "lib_ml_fmm_type_operator_get_coefficient:  ERROR"
                print *, "  m_get_coefficient is not associated"
            end if
        end function

        function lib_ml_fmm_v_operator_add(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

            ! auxilary
            integer :: i

            do i=1, size(lhs)
                rv(i)%dummy = lhs(i)%dummy + rhs(i)%dummy
            end do
        end function
        
        function lib_ml_fmm_m_coefficient_add(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv

            
            rv = m_coefficient_add(lhs, rhs)
            !do i=1, size(lhs)
            !    rv(i)%dummy = lhs(i)%dummy + rhs(i)%dummy
            !end do
        end function

        function lib_ml_fmm_v_operator_add_0d(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_v), intent(in) :: rhs
            type (lib_ml_fmm_v) :: rv

!            allocate(rv%dummy, source = lhs%dummy + rhs%dummy)
            rv%dummy = lhs%dummy + rhs%dummy
        end function

        function lib_ml_fmm_v_operator_sub(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
            type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
            type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

            ! auxilary
            integer :: i

            do i=1, size(lhs)
                rv(i)%dummy = lhs(i)%dummy - rhs(i)%dummy
            end do
        end function

        function lib_ml_fmm_v_operator_sub_0d(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_v), intent(in) :: rhs
            type (lib_ml_fmm_v) :: rv

!            allocate(rv%dummy, source = lhs%dummy - rhs%dummy)
            rv%dummy = lhs%dummy - rhs%dummy
        end function

end module lib_ml_fmm_type_operator
