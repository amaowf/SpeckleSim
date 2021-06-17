module Test_ZGESV
    contains

    subroutine Test_ZGESV_test()
      INTEGER          N, NRHS
      PARAMETER        ( N = 4, NRHS = 1)
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )
!
!     .. Local Scalars ..
      INTEGER          INFO
!
!     .. Local Arrays ..
      INTEGER          IPIV( N )
      COMPLEX*16       A( LDA, N ), B( LDB, NRHS ), B_new( LDB, NRHS ),  A_new( LDA, N )
      !DATA             A/&
      !( 1.23,-5.50),(-2.14,-1.12),(-4.30,-7.10),( 1.27, 7.29),&
      !( 7.91,-5.38),(-9.92,-0.79),(-6.47, 2.52),( 8.90, 6.92),&
      !(-9.80,-4.86),(-9.18,-1.12),(-6.51,-2.67),(-8.82, 1.25),&
      !(-7.32, 7.57),( 1.37, 0.43),(-5.86, 7.38),( 5.41, 5.37)&
      !                 /
      DATA             A/&
      ( 1.23,0.0),(-2.14, 0.0),(-4.30, 0.0),( 1.27, 0.0),&
      ( 7.91,0.0),(-9.92, 0.0),(-6.47, 0.0),( 8.90, 0.0),&
      (-9.80,0.0),(-9.18, 0.0),(-6.51, 0.0),(-8.82, 0.0),&
      (-7.32, 0.0),( 1.37, 0.0),(-5.86,0.0),( 5.41, 0.0)&
                       /
    

      !DATA             B/&
      !( 8.33,-7.32),(-6.18,-4.80),(-5.71,-2.80),(-1.60, 3.08),&
      !(-6.11,-3.81),( 0.14,-7.71),( 1.41, 3.40),( 8.54,-4.05)&
      !                 /

      DATA             B/&
      !( 8.33,-7.32),(-6.18,-4.80),(-5.71,-2.80),(-1.60, 3.08)&
      ( 7.33,0.0),(-6.18,0.0),(-5.71,0.0),(-1.60, 0.0)&
                       /

!
!     .. External Subroutines ..
!      EXTERNAL         ZGESV
!      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR
!
!     .. Executable Statements ..
      WRITE(*,*)'ZGESV Example Program Results'
      A_new = A

!     Solve the equations A*X = B.
      !CALL PRINT_MATRIX( 'Matrix A', N, N, A, LDA ) 
      CALL ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!     Check for the exact singularity.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF

      CALL PRINT_MATRIX( 'Solution', N, NRHS, B, LDB )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test the reverse process
    !print*, 'B =', B
    !print*, 'A(1, 1)', A(1, 1)
    !print*, 'A(1, 2)', A(1, 2)
    !print*, 'A(1, 3)', A(1, 3)
    do i = 1, N
    B_new(i, 1) = (0.0, 0.0)
    do j = 1, N
       print*, B_New(i, 1)
       B_new(i, 1) = B_new(i, 1) + A_new(i, j)*B(j, 1)
    end do
    end do

    CALL PRINT_MATRIX( 'Reverse result', N, NRHS, B_new, LDB )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print details of LU factorization.

          CALL PRINT_MATRIX( 'Details of LU factorization', N, N, A, LDA )

    !    Print pivot indices.
    !
          !CALL PRINT_INT_VECTOR( 'Pivot indices', N, IPIV )
          print*, IPIV
!          STOP
!    END

    contains

        SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
              CHARACTER*(*)    DESC
              INTEGER          M, N, LDA
              COMPLEX*16       A( LDA, * )
        !
              INTEGER          I, J
        !
              WRITE(*,*)
              WRITE(*,*) DESC
              DO I = 1, M
                 WRITE(*,9998) ( A( I, J ), J = 1, N )
              END DO
        !
         9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
              RETURN
              END SUBROUTINE PRINT_MATRIX
        !
        !     Auxiliary routine: printing a vector of integers.
        !
              SUBROUTINE PRINT_INT_VECTOR( DESC, N, A )
              CHARACTER*(*)    DESC
              INTEGER          N
              INTEGER          A( N )
        !
              INTEGER          I
        !
              WRITE(*,*)
              WRITE(*,*) DESC
              WRITE(*,9999) ( A( I ), I = 1, N )
        !
         9999 FORMAT( 11(:,1X, I16) )
              RETURN
              END SUBROUTINE PRINT_INT_VECTOR

      end subroutine Test_ZGESV_test

end module Test_ZGESV
