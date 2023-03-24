module solver_conjugate_gradient_method
    implicit none

    private

    public :: dm_bcg
    public :: dm_bcg_2

    ! parameter
    integer, parameter :: mvulswitch = 10

    interface MVMUL
        module procedure MVMULc
        module procedure MVMULz
    end interface

    interface TMVMUL
        module procedure TMVMULc
        module procedure TMVMULz
    end interface

    interface
        ! Calculates the matrix vector product
        !
        ! Formula: A x = b
        !          A^T x_t = b_t
        !
        ! Argument
        ! ----
        !   x: double complex, dimension(:)
        !       x vector
        !   x_t: double complex, dimension(:), optional
        !       transposed matrix multiplied by x_t or x if x_t is not present.
        !   calc_b: logical, optional (std: .true.)
        !       flag to calculate vector b = A x
        !   calc_b_t: logical, optional (std: .false.)
        !       flag to calculate vector b_t = A^T x[,_t]
        !
        ! Retruns
        ! ----
        !   b: double complex, dimension(:)
        !       result of the product A x
        !   b_t: double complex, dimension(:)
        !       result of the product A^T x[,_t]
        !
        subroutine bcg_calculate_matrix_vector_product(x, b, b_t, x_t, calc_b, calc_b_t)
            implicit none
            ! dummy
            double complex, dimension(:), allocatable, intent(in) :: x
            double complex, dimension(:), allocatable, intent(inout) :: b
            double complex, dimension(:), allocatable, intent(inout) :: b_t

            double complex, dimension(:), allocatable, intent(in), optional :: x_t
            logical, intent(in), optional :: calc_b
            logical, intent(in), optional :: calc_b_t
        end subroutine bcg_calculate_matrix_vector_product
    end interface

    contains

        subroutine io_error(str)
            implicit none
            ! dummy
            character(len=*) :: str

            print *, "io ERROR"
            print *, "  ", str
        end subroutine

        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        !
        ! Changelog
        !   16.09.2019: make double (Max Daiber-Huppert)
        !
        ! CG Algo
        double complex function dm_bcg(x,genau,maxit,sp,b)
        implicit none
        double precision,intent(in)::genau
        integer,intent(in)::maxit
        double complex,dimension(:,:),allocatable,intent(inout)::sp
        double complex,dimension(:),allocatable,intent(inout)::x,b

          double complex,dimension(:),allocatable ::x1,r1,rr1,p1,pp1,r2,rr2,p2,pp2,x2,tmp
          double precision res
          double complex alpha1,beta1
          integer i,dim,status

          dim=size(sp,1)
        allocate(x1(dim),r1(dim),rr1(dim),p1(dim),pp1(dim),r2(dim),rr2(dim),p2(dim),pp2(dim),x2(dim),tmp(dim), STAT=status)
          IF (status.ne.0) call io_error( "*** memory error 139 ***")
          x1=0;r1=0;rr1=0;p1=0;pp1=0;r2=0;rr2=0;p2=0;pp2=0;x2=0

          call MVMUL(tmp, sp,x1)
          r1=b-tmp
          rr1=r1
          p1=r1
          pp1=rr1
          res=1.

          do i=1,maxit
            if (res<=genau) exit
            call MVMUL(tmp, sp,p1)
            alpha1=dot_product(conjg(rr1),r1)/dot_product(conjg(pp1),tmp)
            r2=r1-alpha1*tmp
            call TMVMUL(tmp, sp,pp1)
            rr2=rr1-alpha1*tmp
            beta1=dot_product(conjg(rr2),r2)/dot_product(conjg(rr1),r1)
            p2=r2+beta1*p1
            pp2=rr2+beta1*pp1
            x2=x1+alpha1*p1
            x1=x2
            r1=r2
            rr1=rr2
            p1=p2
            pp1=pp2
            res=sqrt(sum(abs(r1**2))) !VektorNorm
            !call io_printrsr("Step ",real(i),"  RES ",res)
          end do
            x=x1
            dm_bcg=dcmplx(i,res)
        end function

        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        !
        ! Changelog
        !   16.09.2019: make double (Max Daiber-Huppert)
        !
        ! CG Algo
        double complex function dm_bcg_2(x,genau,maxit,b, callback_mvmul, use_init_guess)
        implicit none
        double precision,intent(in)::genau
        integer,intent(in)::maxit
        double complex,dimension(:),allocatable,intent(inout)::x,b
        procedure(bcg_calculate_matrix_vector_product), pointer :: callback_mvmul
        logical, optional :: use_init_guess

          double complex,dimension(:),allocatable ::x1,r1,rr1,p1,pp1,r2,rr2,p2,pp2,x2,tmp,tmp2
          double precision res
          double complex alpha1,beta1
          integer i,dim,status

          dim=size(b,1)
        allocate(x1(dim),r1(dim),rr1(dim),p1(dim),pp1(dim),r2(dim),rr2(dim),p2(dim),pp2(dim),x2(dim),tmp(dim), STAT=status)
          IF (status.ne.0) call io_error( "*** memory error 139 ***")
          x1=0;r1=0;rr1=0;p1=0;pp1=0;r2=0;rr2=0;p2=0;pp2=0;x2=0

          if (present(use_init_guess)) then
              if (use_init_guess) then
!                  call MVMUL(tmp, sp,x1)
                  call callback_mvmul(x, tmp, tmp2)
              end if
          end if
          r1=b-tmp
          rr1=r1
          p1=r1
          pp1=rr1
          res=1.

          allocate(tmp2(dim))
          do i=1,maxit
            if (res<=genau) exit
!            call MVMUL(tmp, sp,p1)
            call callback_mvmul(p1, tmp, tmp2, pp1, calc_b_t = .true.)
            alpha1=dot_product(conjg(rr1),r1)/dot_product(conjg(pp1),tmp)
            r2=r1-alpha1*tmp
!            call TMVMUL(tmp, sp,pp1)
            rr2=rr1-alpha1*tmp2
            beta1=dot_product(conjg(rr2),r2)/dot_product(conjg(rr1),r1)
            p2=r2+beta1*p1
            pp2=rr2+beta1*pp1
            x2=x1+alpha1*p1
            x1=x2
            r1=r2
            rr1=rr2
            p1=p2
            pp1=pp2
            res=sqrt(sum(abs(r1**2))) !VektorNorm
            !call io_printrsr("Step ",real(i),"  RES ",res)
          end do
            x=x1
            dm_bcg_2=dcmplx(i,res)
        end function


        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        subroutine MVMULc(v1,a,v2)
        implicit none
        complex,dimension(:),allocatable,intent(inout):: v1,v2
        complex,dimension(:,:),allocatable,intent(inout):: a
          integer n,m,lda,i
          m=size(a,1)
          if (m<mvulswitch) then
            forall(i=1:m) v1(i)=dot_product(conjg(a(i,:)),v2)
          else
            n=size(a,2);lda=m
            !call cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            call cgemv('n',m,n,cmplx(1.,0.), a,lda, v2, 1, cmplx(0.,0.), v1,1)
          end if
        end subroutine

        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        subroutine TMVMULc(v1,a,v2)
        implicit none
        complex,dimension(:),allocatable,intent(inout):: v1,v2
        complex,dimension(:,:),allocatable,intent(inout):: a
          integer n,m,lda,i
          m=size(a,1)
          if (m<mvulswitch) then
            forall(i=1:m) v1(i)=dot_product(conjg(a(:,i)),v2)
          else
            n=size(a,2);lda=m
            call cgemv('t',m,n,cmplx(1.,0.), a,lda, v2, 1, cmplx(0.,0.), v1,1)
          end if
        end subroutine

        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        ! Changelog
        !   16.09.2019: make double (Max Daiber-Huppert)
        !
        subroutine MVMULz(v1,a,v2)
        implicit none
        double complex,dimension(:),allocatable,intent(inout):: v1,v2
        double complex,dimension(:,:),allocatable,intent(inout):: a
          integer n,m,lda,i
          m=size(a,1)
          if (m<mvulswitch) then
            forall(i=1:m) v1(i)=dot_product(conjg(a(i,:)),v2)
          else
            n=size(a,2);lda=m
            !call cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            call zgemv('n',m,n,dcmplx(1.,0.), a,lda, v2, 1, dcmplx(0.,0.), v1,1)
          end if
        end subroutine

        ! author: Karsten Frenner
        ! date: 11.09.2019
        !
        ! Changelog
        !   16.09.2019: make double (Max Daiber-Huppert)
        !
        subroutine TMVMULz(v1,a,v2)
        implicit none
        double complex,dimension(:),allocatable,intent(inout):: v1,v2
        double complex,dimension(:,:),allocatable,intent(inout):: a
          integer n,m,lda,i
          m=size(a,1)
          if (m<mvulswitch) then
            forall(i=1:m) v1(i)=dot_product(conjg(a(:,i)),v2)
          else
            n=size(a,2);lda=m
            call zgemv('t',m,n,dcmplx(1.,0.), a,lda, v2, 1, dcmplx(0.,0.), v1,1)
          end if
        end subroutine

end module solver_conjugate_gradient_method
