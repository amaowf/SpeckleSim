!    Copyright (C) 2020  Max Daiber-Huppert <max_daiber-huppert@gmx.de>
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
! Created on Thu Jan 30 13:07:51 2020
!
! @author: Max Daiber-Huppert
!

module lib_math_solver_GMRES
!    use libmath
    use zPackgmres
    implicit none

    private

    public :: lib_math_solver_gmres_get_parameter_std_values
    public :: lib_math_solver_gmres_run

    public :: solver_gmres_parameter_type
    public :: solver_gmres_callback_type

    public :: GMRES_test
    public :: GMRES_test_with_callback

    ! parameter
    integer, parameter, public :: GMRES_ORTHOGONALIZATION_SCHEME_MGS = 0
    integer, parameter, public :: GMRES_ORTHOGONALIZATION_SCHEME_IMGS = 1
    integer, parameter, public :: GMRES_ORTHOGONALIZATION_SCHEME_CGS = 2
    integer, parameter, public :: GMRES_ORTHOGONALIZATION_SCHEME_ICGS = 3

    integer, parameter, public :: GMRES_PRECONDITIONING_NO = 0
    integer, parameter, public :: GMRES_PRECONDITIONING_LEFT = 1
    integer, parameter, public :: GMRES_PRECONDITIONING_RIGHT = 2
    integer, parameter, public :: GMRES_PRECONDITIONING_DOUBLE_SIDE = 3
    integer, parameter, public :: GMRES_PRECONDITIONING_ERROR = 4

    type solver_gmres_parameter_type
        integer :: max_iterations
        integer :: restart
        logical :: use_initial_guess
        double precision :: convergence_tolerance
        integer :: orthogonalization_scheme
        logical :: use_recurence_formula_at_restart
        logical :: residual_calc_explicitly
        integer :: no_of_elements_vector_x

        integer :: preconditioning

        double precision :: backward_error
    end type solver_gmres_parameter_type

    ! Internal Callback Function
    ! ----
    !   Reference: libmath_external/algorithm_842/zPackgmres.f90
    !
    !   internal%backward_error:
    !       subroutine zgmres_backward_error(iteration, backward_error_arnoldi, backward_error_true)
    !           implicit none
    !           ! dummy
    !           integer, intent(in) :: iteration
    !           double precision, intent(in) :: backward_error_arnoldi
    !           double precision, intent(in) :: backward_error_true
    !       end subroutine zgmres_backward_error
    type solver_gmres_callback_type
        procedure(lib_math_solver_gmres_get_vector_x_init), pointer, nopass :: get_vector_x_init => null()
        procedure(lib_math_solver_gmres_get_vector_b), pointer, nopass :: get_vector_b => null()
        procedure(lib_math_solver_gmres_calculate_vector_b), pointer, nopass :: calc_vector_b => null()
        procedure(lib_math_solver_gmres_save_vector_x), pointer, nopass :: save_vector_x => null()
        procedure(lib_math_solver_gmres_preconditioner_left), pointer, nopass :: preconditioner_left => null()
        procedure(lib_math_solver_gmres_preconditioner_right), pointer, nopass :: preconditioner_right => null()

        type(zgmres_callback_type) :: internal
    end type solver_gmres_callback_type

    interface
        subroutine lib_math_solver_gmres_get_vector_x_init(vector)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(inout) :: vector
        end subroutine lib_math_solver_gmres_get_vector_x_init

        subroutine lib_math_solver_gmres_get_vector_b(vector)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(inout) :: vector
        end subroutine lib_math_solver_gmres_get_vector_b

        subroutine lib_math_solver_gmres_calculate_vector_b(vector_x, vector_b)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: vector_x
            double complex, dimension(:), allocatable, intent(inout) :: vector_b
        end subroutine lib_math_solver_gmres_calculate_vector_b

        subroutine lib_math_solver_gmres_save_vector_x(vector)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: vector
        end subroutine lib_math_solver_gmres_save_vector_x

        subroutine lib_math_solver_gmres_preconditioner_left(vector_1, vector_2)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: vector_1
            double complex, dimension(:), allocatable, intent(out) :: vector_2
        end subroutine

        subroutine lib_math_solver_gmres_preconditioner_right(vector_1, vector_2)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: vector_1
            double complex, dimension(:), allocatable, intent(out) :: vector_2
        end subroutine

    end interface

    contains

        ! Returns
        ! ----
        !   rv: type(solver_gmres_parameter_type)
        !       default values of the GMRES solver
        function lib_math_solver_gmres_get_parameter_std_values() result(rv)
            implicit none
            ! dummy
            type(solver_gmres_parameter_type) :: rv

            rv%max_iterations = 1000
            rv%restart = 10
            rv%use_initial_guess = .false.
            rv%convergence_tolerance = 1D-15
            rv%orthogonalization_scheme = GMRES_ORTHOGONALIZATION_SCHEME_MGS
            rv%use_recurence_formula_at_restart = .false.
            rv%residual_calc_explicitly = .true.
            rv%preconditioning = GMRES_PRECONDITIONING_NO

        end function lib_math_solver_gmres_get_parameter_std_values

        ! Argument
        ! ----
        !   parameter: type(solver_gmres_parameter_type)
        !       parameter of the GMRES solver
        !   simulation_data: type(lib_mie_simulation_parameter_type) @ lib_mie_ms_data_container
        !       dataset of the simulation
        !   save_solution: logical, optional (std: .true.)
        !       true: save solution x into simulation_data
        subroutine lib_math_solver_gmres_run(gmres_parameter, gmres_callback, save_solution)
            implicit none
            ! dummy
            type(solver_gmres_parameter_type), intent(inout) :: gmres_parameter
            type(solver_gmres_callback_type), intent(in) :: gmres_callback
            logical, intent(in), optional :: save_solution

            ! parameter
            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)

            ! auxiliaray
            integer :: n
            integer :: m
            integer :: lda
            integer :: lwork
            integer :: ldstrt
            integer revcom, colx, coly, colz, nbscal
            integer, dimension(5) :: irc
            integer, dimension(8) :: icntl
            integer, dimension(3) :: info

            integer nout

            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            integer :: m_use_initial_guess
            integer :: m_residual_calc

            ! auxiliary
            logical :: m_save_solution

            double complex, dimension(:), allocatable :: vector
            double complex, dimension(:), allocatable :: vector_2

            ldstrt = gmres_parameter%restart

            if (gmres_parameter%use_initial_guess) then
                m_use_initial_guess = 1
            else
                m_use_initial_guess = 0
            end if

            if (gmres_parameter%residual_calc_explicitly) then
                m_residual_calc = 1
            else
                m_residual_calc = 0
            end if

            m_save_solution = .true.
            if (present(save_solution)) m_save_solution = save_solution

            nout = 6

            lda = gmres_parameter%no_of_elements_vector_x

            ! Reference: A Set of GMRES Routines for Real and Complex Arithmetics on High Performance Computers,
            !            Valerie Frayss LucGiraud Serge Gratton Julien Langou, p. 9
            if (gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_MGS &
                .or. gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_IMGS) then
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 2
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + 2
                end if
            else
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + ldstrt + 1
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + ldstrt + 1
                end if
            end if

            allocate(work(lwork))

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------
            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = gmres_parameter%convergence_tolerance
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = gmres_parameter%max_iterations

            ! preconditioner location
            icntl(4) = gmres_parameter%preconditioning
            ! orthogonalization scheme
            icntl(5)=gmres_parameter%orthogonalization_scheme
            ! initial guess
            icntl(6) = m_use_initial_guess
            ! residual calculation strategy at restart
            icntl(8) = m_residual_calc

            ! set initial guess x_0
            if (gmres_parameter%use_initial_guess) then
                call gmres_callback%get_vector_x_init(vector)
                work(1:lda) = vector
            else
                work(1:lda) = dcmplx(0,0)
            end if

            ! Initialise the right hand side b
            call gmres_callback%get_vector_b(vector)
            work(lda+1:2*lda) = vector
            n = lda
            m = ldstrt

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------
            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo, &
                                  gmres_callback%internal)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
!                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
!                               ZERO,work(colz),1)
                    vector_2 = work(colx:colx+lda-1)
                    call gmres_callback%calc_vector_b(vector_2, vector)
                    work(colz:colz+lda-1) = vector

                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    if (associated(gmres_callback%preconditioner_left)) then
                        vector_2 = work(colx:colx+lda-1)
                        call gmres_callback%preconditioner_left(vector_2, vector)
                        work(colz:colz+lda-1) = vector
								
                        cycle
                    else
                        print *, "lib_math_solver_gmres_run_without_ml_fmm: ERROR"
                        print *, "  no preconditioning (left), but required"

                        exit
                    end if

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning

                    if (associated(gmres_callback%preconditioner_right)) then
                        vector_2 = work(colx:colx+lda-1)
                        call gmres_callback%preconditioner_right(vector_2, vector)
                        work(colz:colz+lda-1) = vector

                        cycle
                    else
                        print *, "lib_math_solver_gmres_run_without_ml_fmm: ERROR"
                        print *, "  no preconditioning (right), but required"

                        exit
                    end if

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                gmres_parameter%backward_error = rinfo(2)

                exit

            end do

            if (m_save_solution) then
                vector = work(1:lda)
                call gmres_callback%save_vector_x(vector)
            end if

        end subroutine lib_math_solver_gmres_run

        subroutine GMRES_test()
            implicit none
            ! dummy

            ! parameter
            integer, parameter :: lda = 1000
            integer, parameter :: ldstrt = 60

            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)


            ! auxiliary
            integer i, j, n, m
            integer :: lwork
            integer revcom, colx, coly, colz, nbscal
            integer irc(5), icntl(8), info(3)

            integer nout

            double complex, dimension(:,:), allocatable :: a
            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            ! ------------------------------------------------------------
            !  Generate the test matrix a and set the right-hand side
            !  in positions (n+1) to 2n of the array work.
            !  The right-hand side is chosen such that the exact solution
            !  is the vector of all ones.
            ! ------------------------------------------------------------

            write(*,*) '***********************************************'
            write(*,*) 'This code is an example of use of GMRES'
            write(*,*) 'in double precision complex arithmetic'
            write(*,*) 'Results are written standard output'
            write(*,*) '***********************************************'
            write(*,*)
!            write(*,*) 'Matrix size < ', lda
            n = 100
            write(*,*) 'Matrix size = ', n
            if (n.gt.lda) then
              write(*,*) 'You are asking for a too large matrix'
              return
            end if

            lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1

            allocate(a(lda, lda))
            allocate(work(lwork))

            write(*,*) "test"

            do j = 1,n
              do i = 1,n
                a(i,j) = ZERO
              end do
            end do

            do i = 1,n
              a(i,i) = dcmplx(4.d0, 0.d0)
            end do
            do i = 1,n-1
              a(i,i+1) = dcmplx(-2.d0, 1.d0)
              a(i+1,i) = dcmplx(-2.d0, 1.d0)
            end do
            nout = 1
            ! -------------------------------
            !  Choose the restart parameter
            ! -------------------------------

!            write(*,*) 'Restart  <', ldstrt
            m = ldstrt
            write(*,*) 'Restart  =', m

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------

            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = 1.d-15
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = 1000

            ! preconditioner location
            icntl(4) = 1
            ! orthogonalization scheme
            icntl(5)=0
            ! initial guess
            icntl(6) = 0
            ! residual calculation strategy at restart
            icntl(8) = 0
            ! Initialise the right hand side
            do j = 1,n
                work(j) = ONE
            end do
            call ZGEMV('N',n,n,ONE,A,lda,work(1),1,ZERO,work(n+1),1)
            do j = 1,n
                work(j) = ONE/2.0
            end do

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------

            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
                               ZERO,work(colz),1)
                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','L','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','U','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                exit

            end do
    end subroutine

    subroutine zgmres_backward_error_test(iteration, backward_error_arnoldi, backward_error_true)
        implicit none
        ! dummy
        integer, intent(in) :: iteration
        double precision, intent(in) :: backward_error_arnoldi
        double precision, intent(in) :: backward_error_true


        print *, "zgmres_backward_error_test: "
        print *, "  iteration: ", iteration
        print *, "  Arnoldi backward error: ", backward_error_arnoldi
        print *, "  True backward error", backward_error_true

    end subroutine zgmres_backward_error_test

    subroutine GMRES_test_with_callback()
            implicit none
            ! dummy

            ! parameter
            integer, parameter :: lda = 1000
            integer, parameter :: ldstrt = 60

            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)


            ! auxiliary
            integer i, j, n, m
            integer :: lwork
            integer revcom, colx, coly, colz, nbscal
            integer irc(5), icntl(8), info(3)

            integer nout

            double complex, dimension(:,:), allocatable :: a
            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            type(solver_gmres_callback_type) :: callback

            callback%internal%backward_error => zgmres_backward_error_test

            ! ------------------------------------------------------------
            !  Generate the test matrix a and set the right-hand side
            !  in positions (n+1) to 2n of the array work.
            !  The right-hand side is chosen such that the exact solution
            !  is the vector of all ones.
            ! ------------------------------------------------------------

            write(*,*) '***********************************************'
            write(*,*) 'This code is an example of use of GMRES'
            write(*,*) 'in double precision complex arithmetic'
            write(*,*) 'Results are written standard output'
            write(*,*) '***********************************************'
            write(*,*)
!            write(*,*) 'Matrix size < ', lda
            n = 100
            write(*,*) 'Matrix size = ', n
            if (n.gt.lda) then
              write(*,*) 'You are asking for a too large matrix'
              return
            end if

            lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1

            allocate(a(lda, lda))
            allocate(work(lwork))

            write(*,*) "test"

            do j = 1,n
              do i = 1,n
                a(i,j) = ZERO
              end do
            end do

            do i = 1,n
              a(i,i) = dcmplx(4.d0, 0.d0)
            end do
            do i = 1,n-1
              a(i,i+1) = dcmplx(-2.d0, 1.d0)
              a(i+1,i) = dcmplx(-2.d0, 1.d0)
            end do
            nout = 1
            ! -------------------------------
            !  Choose the restart parameter
            ! -------------------------------

!            write(*,*) 'Restart  <', ldstrt
            m = ldstrt
            write(*,*) 'Restart  =', m

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------

            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = 1.d-15
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = 1000

            ! preconditioner location
            icntl(4) = 1
            ! orthogonalization scheme
            icntl(5)=0
            ! initial guess
            icntl(6) = 0
            ! residual calculation strategy at restart
            icntl(8) = 0
            ! Initialise the right hand side
            do j = 1,n
                work(j) = ONE
            end do
            call ZGEMV('N',n,n,ONE,A,lda,work(1),1,ZERO,work(n+1),1)
            do j = 1,n
                work(j) = ONE/2.0
            end do

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------

            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo, callback%internal)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
                               ZERO,work(colz),1)
                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','L','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','U','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                exit

            end do
    end subroutine GMRES_test_with_callback

end module lib_math_solver_GMRES
