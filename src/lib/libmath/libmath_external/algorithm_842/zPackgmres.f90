module zPackgmres

!*
!*  Copyright (C) CERFACS 1998
!*
!*  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED
!*  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALL
!*  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD A
!*  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS.
!*
!*  1 - This software program provided in source code format ("the " Sou
!*  Code ") and any associated documentation (the " Documentation ") are
!*  licensed, not sold, to you.
!*
!*  2 - CERFACS grants you a personal, non-exclusive, non-transferable a
!*  royalty-free right to use, copy or modify the Source Code and
!*  Documentation, provided that you agree to comply with the terms and
!*  restrictions of this agreement. You may modify the Source Code and
!*  Documentation to make source code derivative works, object code
!*  derivative works and/or documentation derivative Works (called "    
!*  Derivative Works "). The Source Code, Documentation and Derivative  
!*  Works (called " Licensed Software ") may be used by you for personal
!*  and non-commercial use only. " non-commercial use " means uses that
!*  not or will not result in the sale, lease or rental of the Licensed
!*  Software and/or the use of the Licensed Software in any commercial
!*  product or service. CERFACS reserves all rights not expressly grante
!*  to you. No other licenses are granted or implied.
!*
!*  3 - The Source Code and Documentation are and will remain the sole
!*  property of CERFACS. The Source Code and Documentation are copyright
!*  works. You agree to treat any modification or derivative work of the
!*  Licensed Software as if it were part of the Licensed Software itself
!*  In return for this license, you grant CERFACS a non-exclusive perpet
!*  paid-up royalty-free license to make, sell, have made, copy, distrib
!*  and make derivative works of any modification or derivative work you
!*  make of the Licensed Software.
!*
!*  4- The licensee shall acknowledge the contribution of the Source Cod
!*  (using the reference [1]) in any publication of material dependent u
!*  upon the use of the Source Code. The licensee shall use reasonable
!*  endeavours to notify the authors of the package of this publication.
!*
!*  [1] V. Frayss�, L. Giraud, S. Gratton, and J. Langou, A set of GMR
!*    routines for real and complex arithmetics on high performance
!*    computers, CERFACS Technical Report TR/PA/03/3, public domain soft
!*    available on www.cerfacs/algor/Softs, 2003
!*
!*  5- CERFACS has no obligation to support the Licensed Software it is
!*  providing under this license.
!*
!*  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO
!*  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE
!*  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF
!*  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE
!*  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD
!*  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL
!*  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR
!*  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR
!*  USE OF THE LICENSED SOFTWARE.
!*
!*  6- For information regarding a commercial license for the Source Cod
!*  and Documentation, please contact Mrs Campassens (campasse@cerfacs.f
!*
!*  7- This license is effective until terminated. You may terminate thi
!*  license at any time by destroying the Licensed Software.
!*
!*    I agree all the terms and conditions of the above license agreemen
!*
!
!      type zgmres_callback_type
!       procedure(zgmres_be), pointer, nopass :: backward_error => null(
!      end type zgmres_callback_type
!
!      interface
!        ! Argument
!        ! ----
!        !   it: integer
!        !       iteration number
!        !   bea: double precision
!        !       Arnoldi backward error
!        !   be: double precision
!        !       true backward error
!        subroutine zgmres_be(it, bea, be)
!            implicit none
!            ! dummy
!            integer, intent(out) :: iteration
!            double precision, intent(out) :: backward_error_arnoldi
!            double precision, intent(out) :: backward_error_true
!        end subroutine zgmres_backward_error
!      end interface

    type zgmres_callback_type
        procedure(zgmres_backward_error), pointer, nopass :: backward_error => null()
    end type zgmres_callback_type

    interface
        subroutine zgmres_backward_error(iteration, backward_error_arnoldi, backward_error_true)
            implicit none
            ! dummy
            integer, intent(in) :: iteration
            double precision, intent(in) :: backward_error_arnoldi
            double precision, intent(in) :: backward_error_true
        end subroutine zgmres_backward_error
    end interface

contains


      SUBROUTINE drive_zgmres (n, nloc, m, lwork, work, irc, icntl,     &
      cntl, info, rinfo, callback)
!
!  Purpose
!  =======
!    drive_zgmres is the driver routine for solving the linear system
!  Ax = b using the *  Generalized Minimal Residual iterative method
!  with preconditioning.
!  This solver is implemented with a reverse communication scheme: contr
!  is returned to the user for computing the (preconditioned)
!  matrix-vector product.
!  See the User's Guide for an example of use.
!
!
! Written : June 1996
! Authors : Luc Giraud, Serge Gratton, V. Fraysse
!             Parallel Algorithms - CERFACS
!
! Updated : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS
!
! Updated : June 1998
! Authors : Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS
! Purpose : Make clear that the warning and error messages come from the
!           zgmres modules.
!
! Updated : December 2002 - L. Giraud, J.Langou
! Purpose : Add the capability to avoid explicit residual calculation at
!
!
!  Arguments
!  =========
!
!  n      (input) INTEGER.
!          On entry, the dimension of the problem.
!          Unchanged on exit.
!
!  nloc   (input) INTEGER.
!          On entry, the dimension of the local problem.
!          In a parallel distributed envirionment, this corresponds
!          to the size of the subset of entries of the right hand side
!          and solution allocated to the calling process.
!          Unchanged on exit.
!
!
!  m      (input) INTEGER
!          Restart parameter, <= N. This parameter controls the amount
!          of memory required for matrix H (see WORK and H).
!          Unchanged on exit.
!
!  lwork  (input) INTEGER
!          size of the workspace
!          lwork >= m*m + m*(n+5) + 5*n+1, if icntl(8) = 1
!          lwork >= m*m + m*(n+5) + 6*n+1, if icntl(8) = 0
!
!  work   (workspace) double precision/double complex array, length lwor
!          work contains the required vector and matrices stored in the
!          following order :
!            x  (n,1)       : computed solution.
!            b  (n,1)       : right hand side.
!            r0 (n,1)       : vector workspace.
!            w  (n,1)       : vector workspace.
!            V  (n,m)       : Krylov basis.
!            H  (m+1,m+1)   : Hessenberg matrix (full storage).
!            yCurrent (m,1) : solution of the current LS
!            xCurrent (n,1) : current iterate
!            rotSin (m,1)   : Sine of the Givens rotation
!            rotCos (m,1)   : Cosine of the Givens rotation
!
!  irc     (input/output) INTEGER array. length 5
!            irc(1) : REVCOM   used for reverse communication
!                             (type of external operation)
!            irc(2) : COLX     used for reverse communication
!            irc(3) : COLY     used for reverse communication
!            irc(4) : COLZ     used for reverse communication
!            irc(5) : NBSCAL   used for reverse communication
!
!  icntl   (input) INTEGER array. length 7
!            icntl(1) : stdout for error messages
!            icntl(2) : stdout for warnings
!            icntl(3) : stdout for convergence history
!            icntl(4) : 0 - no preconditioning
!                       1 - left preconditioning
!                       2 - right preconditioning
!                       3 - double side preconditioning
!                       4 - error, default set in Init
!            icntl(5) : 0 - modified Gram-Schmidt
!                       1 - iterative modified Gram-Schmidt
!                       2 - classical Gram-Schmidt
!                       3 - iterative classical Gram-Schmidt
!            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                       1 - user supplied initial guess
!            icntl(7) : maximum number of iterations
!            icntl(8) : 0 - use recurence formula at restart
!                       1 - default compute the true residual at each re
!
!  cntl    (input) double precision array, length 5
!            cntl(1) : tolerance for convergence
!            cntl(2) : scaling factor for normwise perturbation on A
!            cntl(3) : scaling factor for normwise perturbation on b
!            cntl(4) : scaling factor for normwise perturbation on the
!                      preconditioned matrix
!            cntl(5) : scaling factor for normwise perturbation on
!                      preconditioned right hand side
!
!  info    (output) INTEGER array, length 3
!            info(1) :  0 - normal exit
!                      -1 - n < 1
!                      -2 - m < 1
!                      -3 - lwork too small
!                      -4 - convergence not achieved after icntl(7) iter
!                      -5 - precondition type not set by user
!            info(2) : if info(1)=0 - number of iteration to converge
!                      if info(1)=-3 - minimum workspace size necessary
!            info(3) : optimal size for the workspace
!
!  rinfo   (output) double precision array, length 2
!            if info(1)=0
!              rinfo(1) : backward error for the preconditioned system
!              rinfo(2) : backward error for the unpreconditioned system
!
! Input variables
! ---------------
      INTEGER n, nloc, lwork, icntl ( * )
      DOUBLEPRECISION cntl ( * )
      DOUBLEPRECISION sA, sb, sPA, sPb
      type(zgmres_callback_type), intent(in), optional :: callback
! Output variables
! ----------------
      INTEGER info ( * )
      DOUBLEPRECISION rinfo ( * )
! Input/Output variables
! ----------------------
      INTEGER m, irc ( * )
      DOuble complex work ( * )
! Local variables
! ---------------
      INTEGER xptr, bptr, wptr, r0ptr, Vptr, Hptr
      INTEGER yCurrent, rotSin, rotCos, xCurrent
      INTEGER sizeWrk, newRestart
      INTEGER iwarn, ierr, ihist, compRsd
      REAL rn
      DOUBLEPRECISION DZRO
      PARAMETER (DZRO = 0.0d0)
!
      INTEGER icheck
      DATA icheck / 0 /
      SAVE icheck
!
      INTRINSIC ifix
      INTRINSIC float
!
!       Executable statements :
!
      ierr = icntl (1)
      iwarn = icntl (2)
      ihist = icntl (3)
      compRsd = icntl (8)
!
      IF (ierr.lt.0) ierr = 6
!
      IF (compRsd.eq.1) then
         sizeWrk = m * m + m * (nloc + 5) + 5 * nloc + 1
      ELSE
         sizeWrk = m * m + m * (nloc + 5) + 6 * nloc + 1
      ENDIF
!
      IF (icheck.eq.0) then
! Check the value of the arguments
         IF ( (n.lt.1) .or. (nloc.lt.1) ) then
            WRITE (ierr, * )
            WRITE (ierr, * ) ' ERROR GMRES : '
      WRITE (ierr,  * ) '     N < 1 '
            WRITE (ierr, * )
            info (1) = - 1
            irc (1) = 0
            RETURN
         ENDIF
         IF (m.lt.1) then
            WRITE (ierr, * )
            WRITE (ierr, * ) ' ERROR GMRES :'
      WRITE (ierr,  * ) '     M < 1 '
            WRITE (ierr, * )
            info (1) = - 2
            irc (1) = 0
            RETURN
         ENDIF
         IF ( (icntl (4) .ne.0) .and. (icntl (4) .ne.1) .and. (icntl (4)&
         .ne.2) .and. (icntl (4) .ne.3) ) then
            WRITE (ierr, * )
            WRITE (ierr, * ) ' ERROR GMRES : '
      WRITE (ierr,  * ) '     Undefined preconditioner '
            WRITE (ierr, * )
            info (1) = - 5
            irc (1) = 0
            RETURN
         ENDIF
!
         IF ( (icntl (5) .lt.0) .or. (icntl (5) .gt.3) ) then
            icntl (5) = 0
            IF (iwarn.ne.0) then
               WRITE (iwarn, * )
      WRITE (iwarn,  * ) ' WARNING  GMRES : '
      WRITE (iwarn,  * ) '       Undefined orthogonalisation '
      WRITE (iwarn,  * ) '       Default MGS '
               WRITE (iwarn, * )
            ENDIF
         ENDIF
         IF ( (icntl (6) .ne.0) .and. (icntl (6) .ne.1) ) then
            icntl (6) = 0
            IF (iwarn.ne.0) then
               WRITE (iwarn, * )
               WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       Undefined intial guess '
      WRITE (iwarn,  * ) '       Default x0 = 0 '
               WRITE (iwarn, * )
            ENDIF
         ENDIF
         IF (icntl (7) .le.0) then
            icntl (7) = n
            IF (iwarn.ne.0) then
               WRITE (iwarn, * )
               WRITE (iwarn, * ) ' WARNING GMRES :'
      WRITE (iwarn,  * ) '       Negative max number of iterations'
      WRITE (iwarn,  * ) '       Default N '
               WRITE (iwarn, * )
            ENDIF
         ENDIF
         IF ( (icntl (8) .ne.0) .and. (icntl (8) .ne.1) ) then
            icntl (8) = 1
            WRITE (iwarn, * )
            WRITE (iwarn, * ) ' WARNING GMRES :'
      WRITE (iwarn,  * ) '       Undefined strategy for the residual'
      WRITE (iwarn,  * ) '       at restart'
      WRITE (iwarn,  * ) '       Default 1 '
            WRITE (iwarn, * )
         ENDIF
! Check if the restart parameter is correct and if the size of the
!  workspace is big enough for the restart.
! If not try to fix correctly the parameters
!
         IF ( (m.gt.n) .or. (lwork.lt.sizeWrk) ) then
            IF (m.gt.n) then
               m = n
               IF (iwarn.ne.0) then
                  WRITE (iwarn, * )
                  WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       Parameter M bigger than N'
      WRITE (iwarn,  * ) '       New value for M ', m
                  WRITE (iwarn, * )
               ENDIF
               IF (compRsd.eq.1) then
                  sizeWrk = m * m + m * (nloc + 5) + 5 * nloc + 1
               ELSE
                  sizeWrk = m * m + m * (nloc + 5) + 6 * nloc + 1
               ENDIF
            ENDIF
!
            IF ( (lwork.lt.sizeWrk) .and. (n.eq.nloc) ) then
! Compute the maximum size of the m according to the memory space
               rn = float (n)
               newRestart = ifix ( ( - 5.0 - rn + sqrt ( (rn + 5.0) **2 &
               - 4.0 * (5.0 * rn + 1.0 - float (lwork) ) ) ) / 2.0)
               IF (compRsd.eq.0) then
                  newRestart = newRestart - 1
               ENDIF
               IF (newRestart.gt.0) then
                  m = newRestart
                  IF (iwarn.ne.0) then
                     WRITE (iwarn, * )
                     WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       Workspace too small for M'
      WRITE (iwarn,  * ) '       New value for M ', m
                     WRITE (iwarn, * )
                  ENDIF
! Update the value of the optimal size
                  IF (compRsd.eq.1) then
                     sizeWrk = m * m + m * (nloc + 5) + 5 * nloc + 1
                  ELSE
                     sizeWrk = m * m + m * (nloc + 5) + 6 * nloc + 1
                  ENDIF
               ELSE
                  WRITE (ierr, * )
                  WRITE (ierr, * ) ' ERROR GMRES : '
      WRITE (ierr,  * ) '     Not enough space for the problem'
      WRITE (ierr,  * ) '     the space does not permit any m'
                  WRITE (ierr, * )
                  info (1) = - 3
                  irc (1) = 0
                  RETURN
               ENDIF
            ENDIF
            IF ( (lwork.lt.sizeWrk) .and. (n.ne.nloc) ) then
               WRITE (ierr, * )
               WRITE (ierr, * ) ' ERROR GMRES : '
      WRITE (ierr,  * ) '     Not enough space for the problem'
               WRITE (ierr, * )
               info (1) = - 3
               irc (1) = 0
               RETURN
            ENDIF
         ENDIF
!
         IF (iwarn.ne.0) then
            WRITE (iwarn, * )
            WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       For M = ', m, ' optimal value '
      WRITE (iwarn,  * ) '       for LWORK =  ', sizeWrk
            WRITE (iwarn, * )
         ENDIF
!
         info (3) = sizeWrk
         icheck = 1
!
! save the parameters the the history file
!
         IF (ihist.ne.0) then
            WRITE (ihist, '(10x,A39)') 'CONVERGENCE HISTORY FOR GMRES'
            WRITE (ihist, * )
            WRITE (ihist, '(A30,I2)') 'Errors are displayed in unit: ', &
            ierr
            IF (iwarn.eq.0) then
               WRITE (ihist, '(A27)') 'Warnings are not displayed:'
            ELSE
      WRITE (ihist, '(A32,I2)') 'Warnings are displayed in unit: ', iwarn
            ENDIF
            WRITE (ihist, '(A13,I7)') 'Matrix size: ', n
            WRITE (ihist, '(A19,I7)') 'Local matrix size: ', nloc
            WRITE (ihist, '(A9,I7)') 'Restart: ', m
            IF (icntl (4) .eq.0) then
               WRITE (ihist, '(A18)') 'No preconditioning'
            ELSEIF (icntl (4) .eq.1) then
               WRITE (ihist, '(A20)') 'Left preconditioning'
            ELSEIF (icntl (4) .eq.2) then
               WRITE (ihist, '(A21)') 'Right preconditioning'
            ELSEIF (icntl (4) .eq.3) then
               WRITE (ihist, '(A30)') 'Left and right preconditioning'
            ENDIF
            IF (icntl (5) .eq.0) then
               WRITE (ihist, '(A21)') 'Modified Gram-Schmidt'
            ELSEIF (icntl (5) .eq.1) then
               WRITE (ihist, '(A31)') 'Iterative modified Gram-Schmidt'
            ELSEIF (icntl (5) .eq.2) then
               WRITE (ihist, '(A22)') 'Classical Gram-Schmidt'
            ELSE
      WRITE (ihist, '(A32)') 'Iterative classical Gram-Schmidt'
            ENDIF
            IF (icntl (6) .eq.0) then
               WRITE (ihist, '(A29)') 'Default initial guess x_0 = 0'
            ELSE
               WRITE (ihist, '(A27)') 'User supplied initial guess'
            ENDIF
            IF (icntl (8) .eq.1) then
      WRITE (ihist, '(A33)') 'True residual computed at restart'
            ELSE
               WRITE (ihist, '(A30)') 'Recurrence residual at restart'
            ENDIF
            WRITE (ihist, '(A30,I5)') 'Maximum number of iterations: ', &
            icntl (7)
            WRITE (ihist, '(A27,E8.2)') 'Tolerance for convergence: ',  &
            cntl (1)
!
      WRITE (ihist, '(A53)') 'Backward error on the unpreconditioned sys&
     &tem Ax = b:'
            sA = cntl (2)
            sb = cntl (3)
            IF ( (sA.eq.DZRO) .and. (sb.eq.DZRO) ) then
      WRITE (ihist, '(A39)') '    the residual is normalised by ||b||'
            ELSE
               WRITE (ihist, 1) sA, sb
    1 FORMAT('    the residual is normalised by         ',E8.2,         &
     &       ' * ||x|| + ',E8.2)
            ENDIF
            sPA = cntl (4)
            sPb = cntl (5)
            WRITE (ihist, 2)
    2 FORMAT('Backward error on the preconditioned system',             &
     &        ' (P1)A(P2)y = (P1)b:')
            IF ( (sPA.eq.DZRO) .and. (sPb.eq.DZRO) ) then
               WRITE (ihist, 3)
    3 FORMAT('    the preconditioned residual is normalised ',          &
     &       'by ||(P1)b||')
            ELSE
               WRITE (ihist, 4) sPA, sPb
    4 FORMAT('    the preconditioned residual is normalised by ', E8.2, &
     &        ' * ||(P2)y|| + ',E8.2)
            ENDIF
!
            WRITE (ihist, '(A31,I7)') 'Optimal size for the workspace:',&
            info (3)
            WRITE (ihist, * )
            WRITE (ihist, 5)
    5 FORMAT('Convergence history: b.e. on the preconditioned system')
            WRITE (ihist, 6)
    6 FORMAT(' Iteration   Arnoldi b.e.    True b.e.')
         ENDIF
!
      ENDIF
! setup some pointers on the workspace
      xptr = 1
      bptr = xptr + nloc
      r0ptr = bptr + nloc
      wptr = r0ptr + nloc
      Vptr = wptr + nloc
      IF (compRsd.eq.1) then
         Hptr = Vptr + m * nloc
      ELSE
         Hptr = Vptr + (m + 1) * nloc
      ENDIF
      yCurrent = Hptr + (m + 1) * (m + 1)
      xCurrent = yCurrent + m
      rotSin = xCurrent + nloc
      rotCos = rotSin + m
!
      CALL zgmres (nloc, m, work (bptr), work (xptr), work (Hptr),      &
      work (wptr), work (r0ptr), work (Vptr), work (yCurrent), work (   &
      xCurrent), work (rotSin), work (rotCos), irc, icntl, cntl, info,  &
      rinfo, callback)
!
      IF (irc (1) .eq.0) then
         icheck = 0
      ENDIF
!
      RETURN
      END SUBROUTINE drive_zgmres
!
      SUBROUTINE zgmres (n, m, b, x, H, w, r0, V, yCurrent, xCurrent,   &
      rotSin, rotCos, irc, icntl, cntl, info, rinfo, callback)
!
!
!  Purpose
!  =======
!  zgmres solves the linear system Ax = b using the
!  Generalized Minimal Residual iterative method
!
! When preconditioning is used we solve :
!     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
!     x = M_2^{-1} y
!
!   Convergence test based on the normwise backward error for
!  the preconditioned system
!
! Written : June 1996
! Authors : Luc Giraud, Serge Gratton, V. Fraysse
!             Parallel Algorithms - CERFACS
!
! Updated : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS
!
! Updated : March 1998
! Purpose : Pb with F90 on DEC ws
!           cure : remove "ZDSCAL" when used to initialize vectors to ze
!
! Updated : May 1998
! Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot produc
!           cure : w(1) <--  r0'r0
!
! Updated : June 1998
! Purpose : Make clear that the warning and error messages come from the
!           zgmres modules.
!
! Updated : February 2001 - L. Giraud
! Purpose : In complex version, initializations to zero performed  in co
!           arithmetic to avoid implicit conversion by the compiler.
!
! Updated : July 2001 - L. Giraud, J. Langou
! Purpose : Avoid to compute the approximate solution at each step of
!           the Krylov space construction when spA is zero.
!
! Updated : November 2002 - S. Gratton
! Purpose : Use Givens rotations conform to the classical definition.
!           No impact one the convergence history.
!
! Updated : November 2002 - L. Giraud
! Purpose : Properly handle the situation when the convergence is obtain
!           exactly at the "IterMax" iteration
!
! Updated : December 2002 - L. Giraud, J.Langou
! Purpose : Add the capability to avoid explicit residual calculation at
!
! Updated : January  2003 - L. Giraud, S. Gratton
! Purpose : Use Givens rotations from BLAS.
!
! Updated : March    2003 - L. Giraud
! Purpose : Set back retlbl to zero, if initial guess is solution
!           or right-hand side is zero
!
!  Arguments
!  =========
!
!  n       (input) INTEGER.
!           On entry, the dimension of the problem.
!           Unchanged on exit.
!
!  m        (input) INTEGER
!           Restart parameter, <= N. This parameter controls the amount
!           of memory required for matrix H (see WORK and H).
!           Unchanged on exit.
!
!  b        (input) double precision/double complex
!           Right hand side of the linear system.
!
!  x        (output) double precision/double complex
!           Computed solution of the linear system.
!
!  H        (workspace)  double precision/double complex
!           Hessenberg matrix built within dgmres
!
!  w        (workspace)  double precision/double complex
!           Vector used as temporary storage
!
!  r0       (workspace)  double precision/double complex
!           Vector used as temporary storage
!
!  V        (workspace)  double precision/double complex
!           Basis computed by the Arnoldi's procedure.
!
!  yCurrent (workspace) double precision/double complex
!           solution of the current LS
!
!  xCurrent (workspace) double precision/double complex
!           current iterate
!
!  rotSin   (workspace) double precision/double complex
!           Sine of the Givens rotation
!
!  rotCos   (workspace) double precision
!           Cosine of the Givens rotation
!
!  irc      (input/output) INTEGER array. length 3
!             irc(1) : REVCOM   used for reverse communication
!                              (type of external operation)
!             irc(2) : COLX     used for reverse communication
!             irc(3) : COLY     used for reverse communication
!             irc(4) : COLZ     used for reverse communication
!             irc(5) : NBSCAL   used for reverse communication
!
!  icntl    (input) INTEGER array. length 7
!             icntl(1) : stdout for error messages
!             icntl(2) : stdout for warnings
!             icntl(3) : stdout for convergence history
!             icntl(4) : 0 - no preconditioning
!                        1 - left preconditioning
!                        2 - right preconditioning
!                        3 - double side preconditioning
!                        4 - error, default set in Init
!             icntl(5) : 0 - modified Gram-Schmidt
!                        1 - iterative modified Gram-Schmidt
!                        2 - classical Gram-Schmidt
!                        3 - iterative classical Gram-Schmidt
!             icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                        1 - user supplied initial guess
!             icntl(7) : maximum number of iterations
!             icntl(8) : 1 - default compute the true residual at each r
!                        0 - use recurence formula at restart
!
!  cntl     (input) double precision array, length 5
!             cntl(1) : tolerance for convergence
!             cntl(2) : scaling factor for normwise perturbation on A
!             cntl(3) : scaling factor for normwise perturbation on b
!             cntl(4) : scaling factor for normwise perturbation on the
!                       preconditioned matrix
!             cntl(5) : scaling factor for normwise perturbation on
!                       preconditioned right hand side
!
!  info     (output) INTEGER array, length 2
!             info(1) :  0 - normal exit
!                       -1 - n < 1
!                       -2 - m < 1
!                       -3 - lwork too small
!                       -4 - convergence not achieved after icntl(7) ite
!                       -5 - precondition type not set by user
!             info(2) : if info(1)=0 - number of iteration to converge
!                       if info(1)=-3 - minimum workspace size necessary
!             info(3) : optimal size for the workspace
!
! rinfo     (output) double precision array, length 2
!             if info(1)=0
!               rinfo(1) : backward error for the preconditioned system
!               rinfo(2) : backward error for the unpreconditioned syste
!
! Input variables
! ---------------
      INTEGER n, m, icntl ( * )
      double complex b ( * )
      DOUBLEPRECISION cntl ( * )
      type(zgmres_callback_type), intent(in), optional :: callback
!
! Output variables
! ----------------
      INTEGER info ( * )
      DOUBLEPRECISION rinfo ( * )
!
! Input/Output variables
! ----------------------
      INTEGER irc ( * )
      double complex x ( * ), H (m + 1, * ), w ( * ), r0 ( * ), V (n,    &
      * ), yCurrent ( * )
      double complex xCurrent ( * ), rotSin ( * )
      double complex rotCos ( * )
!
! Local variables
! ---------------
      INTEGER j, jH, iterOut, nOrtho, iterMax, initGuess, iOrthog
      INTEGER xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
      INTEGER typePrec, leftPrec, rightPrec, dblePrec, noPrec
      INTEGER iwarn, ihist
      INTEGER compRsd
      DOUBLEPRECISION beta, bn, sA, sb, sPA, sPb, bea, be, temp
      DOUBLEPRECISION dloo, dnormw, dnormx, dnormres, trueNormRes
      double complex dVi, aux
      double complex auxHjj, auxHjp1j
!
      PARAMETER (noPrec = 0, leftPrec = 1)
      PARAMETER (rightPrec = 2, dblePrec = 3)
!
      double complex ZERO, ONE
      PARAMETER (ZERO = (0.0d0, 0.0d0), ONE = (1.0d0, 0.0d0) )
      DOUBLEPRECISION DZRO, DONE
      PARAMETER (DZRO = 0.0d0, DONE = 1.0d0)
!
!
! External functions
! ------------------
      DOUBLEPRECISION dznrm2
      EXTERNAL dznrm2
!
! Reverse communication variables
! -------------------------------
      INTEGER retlbl
      DATA retlbl / 0 /
      INTEGER matvec, precondLeft, precondRight, prosca
      PARAMETER (matvec = 1, precondLeft = 2, precondRight = 3, prosca =&
      4)
!
! Saved variables
! ---------------
      SAVE iterOut, jH, beta, bn, dnormres, retlbl, j
      SAVE sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be
      SAVE dloo, nOrtho, compRsd
!
! Intrinsic function
! ------------------
      INTRINSIC abs, dsqrt, dreal, dcmplx
!
!       Executable statements
!
! setup some pointers on the workspace
      xptr = 1
      bptr = xptr + n
      r0ptr = bptr + n
      wptr = r0ptr + n
      Vptr = wptr + n
      IF (icntl (8) .eq.1) then
         Hptr = Vptr + m * n
      ELSE
         Hptr = Vptr + (m + 1) * n
      ENDIF
      yptr = Hptr + (m + 1) * (m + 1)
      xcuptr = yptr + m
!
      iwarn = icntl (2)
      ihist = icntl (3)
      typePrec = icntl (4)
      iOrthog = icntl (5)
      initGuess = icntl (6)
      iterMax = icntl (7)
!
      IF (retlbl.eq.0) then
         compRsd = icntl (8)
      ENDIF
!
      IF (retlbl.ne.0) then
         IF (retlbl.eq.5) then
            GOTO 5
         ELSEIF (retlbl.eq.6) then
            GOTO 6
         ELSEIF (retlbl.eq.8) then
            GOTO 8
         ELSEIF (retlbl.eq.11) then
            GOTO 11
         ELSEIF (retlbl.eq.16) then
            GOTO 16
         ELSEIF (retlbl.eq.18) then
            GOTO 18
         ELSEIF (retlbl.eq.21) then
            GOTO 21
         ELSEIF (retlbl.eq.26) then
            GOTO 26
         ELSEIF (retlbl.eq.31) then
            GOTO 31
         ELSEIF (retlbl.eq.32) then
            GOTO 32
         ELSEIF (retlbl.eq.33) then
            GOTO 33
         ELSEIF (retlbl.eq.34) then
            GOTO 34
         ELSEIF (retlbl.eq.36) then
            GOTO 36
         ELSEIF (retlbl.eq.37) then
            GOTO 37
         ELSEIF (retlbl.eq.38) then
            GOTO 38
         ELSEIF (retlbl.eq.41) then
            GOTO 41
         ELSEIF (retlbl.eq.43) then
            GOTO 43
         ELSEIF (retlbl.eq.46) then
            GOTO 46
         ELSEIF (retlbl.eq.48) then
            GOTO 48
         ELSEIF (retlbl.eq.51) then
            GOTO 51
         ELSEIF (retlbl.eq.52) then
            GOTO 52
         ELSEIF (retlbl.eq.61) then
            GOTO 61
         ELSEIF (retlbl.eq.66) then
            GOTO 66
         ELSEIF (retlbl.eq.68) then
            GOTO 68
         ENDIF
      ENDIF
!
!
! intialization of various variables
!
      iterOut = 0
      beta = DZRO
!
      IF (initGuess.eq.0) then
         DO j = 1, n
         x (j) = ZERO
         enddo
      ENDIF
!
!        bn = dznrm2(n,b,1)
!
      irc (1) = prosca
      irc (2) = bptr
      irc (3) = bptr
      irc (4) = r0ptr
      irc (5) = 1
      retlbl = 5
      RETURN
    5 CONTINUE
      bn = dsqrt (dreal (r0 (1) ) )
!
      IF (bn.eq.DZRO) then
         DO j = 1, n
         x (j) = ZERO
         enddo
         IF (iwarn.ne.0) then
            WRITE (iwarn, * )
            WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       Null right hand side'
      WRITE (iwarn,  * ) '       solution set to zero'
            WRITE (iwarn, * )
         ENDIF
         info (1) = 0
         info (2) = 0
         rinfo (1) = DZRO
         rinfo (2) = DZRO
         irc (1) = 0
         retlbl = 0
         RETURN
      ENDIF
!
! Compute the scaling factor for the backward error on the
!  unpreconditioned sytem
!
      sA = cntl (2)
      sb = cntl (3)
      IF ( (sA.eq.DZRO) .and. (sb.eq.DZRO) ) then
         sb = bn
      ENDIF
! Compute the scaling factor for the backward error on the
!  preconditioned sytem
!
      sPA = cntl (4)
      sPb = cntl (5)
      IF ( (sPA.eq.DZRO) .and. (sPb.eq.DZRO) ) then
         IF ( (typePrec.eq.noPrec) .or. (typePrec.eq.rightPrec) ) then
            sPb = bn
         ELSE
            irc (1) = precondLeft
            irc (2) = bptr
            irc (4) = r0ptr
            retlbl = 6
            RETURN
         ENDIF
      ENDIF
    6 CONTINUE
      IF ( (sPA.eq.DZRO) .and. (sPb.eq.DZRO) ) then
         IF ( (typePrec.eq.dblePrec) .or. (typePrec.eq.leftPrec) ) then
!
!           sPb = dznrm2(n,r0,1)
!
            irc (1) = prosca
            irc (2) = r0ptr
            irc (3) = r0ptr
            irc (4) = wptr
            irc (5) = 1
            retlbl = 8
            RETURN
         ENDIF
      ENDIF
    8 CONTINUE
      IF ( (sPA.eq.DZRO) .and. (sPb.eq.DZRO) ) then
         IF ( (typePrec.eq.dblePrec) .or. (typePrec.eq.leftPrec) ) then
            sPb = dsqrt (dreal (w (1) ) )
!
         ENDIF
      ENDIF
!
!
! Compute the first residual
!           Y = AX : r0 <-- A x
!
! The residual is computed only if the initial guess is not zero
!
      IF (initGuess.ne.0) then
         irc (1) = matvec
         irc (2) = xptr
         irc (4) = r0ptr
         retlbl = 11
         RETURN
      ENDIF
   11 CONTINUE
      IF (initGuess.ne.0) then
         DO j = 1, n
         r0 (j) = b (j) - r0 (j)
         enddo
      ELSE
         CALL zcopy (n, b, 1, r0, 1)
      ENDIF
!
! Compute the preconditioned residual if necessary
!      M_1Y = X : w <-- M_1^{-1} r0
!
      IF ( (typePrec.eq.noPrec) .or. (typePrec.eq.rightPrec) ) then
         CALL zcopy (n, r0, 1, w, 1)
      ELSE
         irc (1) = precondLeft
         irc (2) = r0ptr
         irc (4) = wptr
         retlbl = 16
         RETURN
      ENDIF
   16 CONTINUE
!
!
!       beta = dznrm2(n,w,1)
!
!
      irc (1) = prosca
      irc (2) = wptr
      irc (3) = wptr
      irc (4) = r0ptr
      irc (5) = 1
      retlbl = 18
      RETURN
   18 CONTINUE
      beta = dsqrt (dreal (r0 (1) ) )
!
      IF (beta.eq.DZRO) then
!  The residual is exactly zero : x is the exact solution
         info (1) = 0
         info (2) = 0
         rinfo (1) = DZRO
         rinfo (2) = DZRO
         irc (1) = 0
         retlbl = 0
         IF (iwarn.ne.0) then
            WRITE (iwarn, * )
            WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       Intial residual is zero'
      WRITE (iwarn,  * ) '       initial guess is solution'
            WRITE (iwarn, * )
         ENDIF
         RETURN
      ENDIF
!
      aux = ONE / beta
      DO j = 1, n
      V (j, 1) = ZERO
      enddo
      CALL zaxpy (n, aux, w, 1, V (1, 1), 1)
!
!       Most outer loop : zgmres iteration
!
!       REPEAT
    7 CONTINUE
!
!
      H (1, m + 1) = beta
      DO j = 1, m
      H (j + 1, m + 1) = ZERO
      enddo
!
!        Construction of the hessenberg matrix WORK and of the orthogona
!        basis V such that AV=VH
!
      jH = 1
   10 CONTINUE
! Remark : this  do loop has been written with a while do
!          because the
!               " do jH=1,restart "
!         fails with the reverse communication.
!      do  jH=1,restart
!
!
! Compute the preconditioned residual if necessary
!
      IF ( (typePrec.eq.rightPrec) .or. (typePrec.eq.dblePrec) ) then
!
!           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH)
!
         irc (1) = precondRight
         irc (2) = vptr + (jH - 1) * n
         irc (4) = wptr
         retlbl = 21
         RETURN
      ELSE
         CALL zcopy (n, V (1, jH), 1, w, 1)
      ENDIF
   21 CONTINUE
!
!           Y = AX : r0 <-- A w
!
      irc (1) = matvec
      irc (2) = wptr
      irc (4) = r0ptr
      retlbl = 26
      RETURN
   26 CONTINUE
!
!      MY = X : w <-- M_1^{-1} r0
!
      IF ( (typePrec.eq.noPrec) .or. (typePrec.eq.rightPrec) ) then
         CALL zcopy (n, r0, 1, w, 1)
      ELSE
         irc (1) = precondLeft
         irc (2) = r0ptr
         irc (4) = wptr
         retlbl = 31
         RETURN
      ENDIF
   31 CONTINUE
!
! Orthogonalization using either MGS or IMGS
!
! initialize the Hessenberg matrix to zero in order to be able to use
!     IMGS as orthogonalization procedure.
      DO j = 1, jH
      H (j, jH) = ZERO
      enddo
      nOrtho = 0
   19 CONTINUE
      nOrtho = nOrtho + 1
      dloo = DZRO
!
      IF ( (iOrthog.eq.0) .or. (iOrthog.eq.1) ) then
! MGS
!
!           do j=1,jH
!
         j = 1
!           REPEAT
      ENDIF
   23 CONTINUE
      IF ( (iOrthog.eq.0) .or. (iOrthog.eq.1) ) then
!
!             dVi     = zdotc(n,V(1,j),1,w,1)
!
         irc (1) = prosca
         irc (2) = vptr + (j - 1) * n
         irc (3) = wptr
         irc (4) = r0ptr
         irc (5) = 1
         retlbl = 32
         RETURN
      ENDIF
   32 CONTINUE
      IF ( (iOrthog.eq.0) .or. (iOrthog.eq.1) ) then
         dVi = r0 (1)
         H (j, jH) = H (j, jH) + dVi
         dloo = dloo + abs (dVi) **2
         aux = - ONE * dVi
         CALL zaxpy (n, aux, V (1, j), 1, w, 1)
         j = j + 1
         IF (j.le.jH) goto 23
!          enddo_j
      ELSE
! CGS
! Gathered dot product calculation
!
!           call zgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)
!
         irc (1) = prosca
         irc (2) = vptr
         irc (3) = wptr
         irc (4) = r0ptr
         irc (5) = jH
         retlbl = 34
         RETURN
      ENDIF
   34 CONTINUE
      IF ( (iOrthog.eq.2) .or. (iOrthog.eq.3) ) then
!
         CALL zaxpy (jH, ONE, r0, 1, H (1, jH), 1)
         CALL zgemv ('N', n, jH, - ONE, V (1, 1) , n, r0, 1, ONE, w, 1)
         dloo = dznrm2 (jH, r0, 1) **2
      ENDIF
!
!         dnormw = dznrm2(n,w,1)
!
      irc (1) = prosca
      irc (2) = wptr
      irc (3) = wptr
      irc (4) = r0ptr
      irc (5) = 1
      retlbl = 33
      RETURN
   33 CONTINUE
      dnormw = dsqrt (dreal (r0 (1) ) )
!
      IF ( (iOrthog.eq.1) .or. (iOrthog.eq.3) ) then
! IMGS / CGS orthogonalisation
         dloo = dsqrt (dloo)
! check the orthogonalization quality
         IF ( (dnormw.le.dloo) .and. (nOrtho.lt.3) ) then
            GOTO 19
         ENDIF
      ENDIF
!
      H (jH + 1, jH) = dnormw
      IF ( (jH.lt.m) .or. (icntl (8) .eq.0) ) then
         aux = ONE / dnormw
         DO j = 1, n
         V (j, jH + 1) = ZERO
         enddo
         CALL zaxpy (n, aux, w, 1, V (1, jH + 1), 1)
      ENDIF
! Apply previous Givens rotations to the new column of H
      DO j = 1, jH - 1
      CALL zrot (1, H (j, jH), 1, H (j + 1, jH), 1, dreal (rotCos (j) ),&
      rotSin (j) )
      enddo
      auxHjj = H (jH, jH)
      auxHjp1j = H (jH + 1, jH)
      CALL zrotg (auxHjj, auxHjp1j, temp, rotSin (jH) )
      rotCos (jH) = dcmplx (temp, DZRO)
! Apply current rotation to the rhs of the least squares problem
      CALL zrot (1, H (jH, m + 1), 1, H (jH + 1, m + 1), 1, dreal (     &
      rotCos (jH) ), rotSin (jH) )
!
! zabs(H(jH+1,m+1)) is the residual computed using the least squares
!          solver
! Complete the QR factorisation of the Hessenberg matrix by apply the cu
! rotation to the last entry of the collumn
      CALL zrot (1, H (jH, jH), 1, H (jH + 1, jH), 1, dreal (rotCos (jH)&
      ), rotSin (jH) )
      H (jH + 1, jH) = ZERO
!
! Get the Least square residual
!
      dnormres = abs (H (jH + 1, m + 1) )
      IF (sPa.ne.DZRO) then
!
! Compute the solution of the current linear least squares problem
!
         CALL zcopy (jH, H (1, m + 1), 1, yCurrent, 1)
         CALL ztrsv ('U', 'N', 'N', jH, H, m + 1, yCurrent, 1)
!
! Compute the value of the new iterate
!
         CALL zgemv ('N', n, jH, ONE, v, n, yCurrent, 1, ZERO, xCurrent,&
         1)
!
         IF ( (typePrec.eq.rightPrec) .or. (typePrec.eq.dblePrec) )     &
         then
!
!         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!
            irc (1) = precondRight
            irc (2) = xcuptr
            irc (4) = r0ptr
            retlbl = 36
            RETURN
         ELSE
            CALL zcopy (n, xCurrent, 1, r0, 1)
         ENDIF
      ENDIF
   36 CONTINUE
!
!
      IF (sPa.ne.DZRO) then
! Update the current solution
         CALL zcopy (n, x, 1, xCurrent, 1)
         CALL zaxpy (n, ONE, r0, 1, xCurrent, 1)
!
!         dnormx = dznrm2(n,xCurrent,1)
!
         irc (1) = prosca
         irc (2) = xcuptr
         irc (3) = xcuptr
         irc (4) = r0ptr
         irc (5) = 1
         retlbl = 38
         RETURN
      ELSE
         dnormx = DONE
      ENDIF
   38 CONTINUE
      IF (sPa.ne.DZRO) then
         dnormx = dsqrt (dreal (r0 (1) ) )
      ENDIF
!
      bea = dnormres / (sPA * dnormx + sPb)
!
! Check the convergence based on the Arnoldi Backward error for the
! preconditioned system
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
!
! The Arnoldi Backward error indicates that zgmres might have converge
! enforce the calculation of the true residual at next restart
         compRsd = 1
!
!  If the update of X has not yet been performed
         IF (sPA.eq.DZRO) then
!
! Compute the solution of the current linear least squares problem
!
            CALL zcopy (jH, H (1, m + 1), 1, yCurrent, 1)
            CALL ztrsv ('U', 'N', 'N', jH, H, m + 1, yCurrent, 1)
!
! Compute the value of the new iterate
!
            CALL zgemv ('N', n, jH, ONE, v, n, yCurrent, 1, ZERO,       &
            xCurrent, 1)
!
            IF ( (typePrec.eq.rightPrec) .or. (typePrec.eq.dblePrec) )  &
            then
!
!         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!
               irc (1) = precondRight
               irc (2) = xcuptr
               irc (4) = r0ptr
               retlbl = 37
               RETURN
            ELSE
               CALL zcopy (n, xCurrent, 1, r0, 1)
            ENDIF
         ENDIF
      ENDIF
   37 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
         IF (sPA.eq.DZRO) then
! Update the current solution
            CALL zcopy (n, x, 1, xCurrent, 1)
            CALL zaxpy (n, ONE, r0, 1, xCurrent, 1)
         ENDIF
!
         CALL zcopy (n, xCurrent, 1, r0, 1)
! Compute the true residual, the Arnoldi one may be unaccurate
!
!           Y = AX : w  <-- A r0
!
         irc (1) = matvec
         irc (2) = r0ptr
         irc (4) = wptr
         retlbl = 41
         RETURN
      ENDIF
   41 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
!
         DO j = 1, n
         w (j) = b (j) - w (j)
         enddo
! Compute the norm of the unpreconditioned residual
!
!        trueNormRes = dznrm2(n,w,1)
!
         irc (1) = prosca
         irc (2) = wptr
         irc (3) = wptr
         irc (4) = r0ptr
         irc (5) = 1
         retlbl = 43
         RETURN
      ENDIF
   43 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
         trueNormRes = dsqrt (dreal (r0 (1) ) )
!
         IF ( (typePrec.eq.leftPrec) .or. (typePrec.eq.dblePrec) ) then
!
!      MY = X : r0 <-- M_1^{-1} w
!
            irc (1) = precondLeft
            irc (2) = wptr
            irc (4) = r0ptr
            retlbl = 46
            RETURN
         ELSE
            CALL zcopy (n, w, 1, r0, 1)
         ENDIF
      ENDIF
   46 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
!
!        dnormres = dznrm2(n,r0,1)
!
         irc (1) = prosca
         irc (2) = r0ptr
         irc (3) = r0ptr
         irc (4) = wptr
         irc (5) = 1
         retlbl = 48
         RETURN
      ENDIF
   48 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
         dnormres = dsqrt (dreal (w (1) ) )
!
         be = dnormres / (sPA * dnormx + sPb)
! Save the backward error on a file if convergence history requested
         IF (ihist.ne.0) then
            WRITE (ihist, 1000) iterOut * m + jH, bea, be
 1000 FORMAT(I5,11x,E8.2,7x,E8.2)
         ENDIF

        if (present(callback)) then
            if (associated(callback%backward_error)) then
                call callback%backward_error(iterOut*m+jH, bea, be)
            end if
        end if
!
      ENDIF
!
!
! Check again the convergence
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
         IF ( (be.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) )    &
         then
! The convergence has been achieved, we restore the solution in x
! and compute the two backward errors.
            CALL zcopy (n, xCurrent, 1, x, 1)
!
            IF (sA.ne.DZRO) then
!
!            dnormx = dznrm2(n,x,1)
!
               irc (1) = prosca
               irc (2) = xptr
               irc (3) = xptr
               irc (4) = r0ptr
               irc (5) = 1
               retlbl = 51
               RETURN
            ENDIF
         ENDIF
      ENDIF
   51 CONTINUE
      IF ( (bea.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) ) then
         IF ( (be.le.cntl (1) ) .or. (iterOut * m + jH.ge.iterMax) )    &
         then
            IF (sA.ne.DZRO) then
               dnormx = dsqrt (dreal (r0 (1) ) )
!
            ELSE
               dnormx = DONE
            ENDIF
! Return the backward errors
            rinfo (1) = be
            rinfo (2) = trueNormRes / (sA * dnormx + sb)
            IF (be.le.cntl (1) ) then
               info (1) = 0
               IF (ihist.ne.0) then
                  WRITE (ihist, * )
                  WRITE (ihist, '(A20)') 'Convergence achieved'
               ENDIF
            ELSEIF (be.gt.cntl (1) ) then
               IF (iwarn.ne.0) then
                  WRITE (iwarn, * )
                  WRITE (iwarn, * ) ' WARNING GMRES : '
      WRITE (iwarn,  * ) '       No convergence after '
                  WRITE (iwarn, * ) iterOut * m + jH, ' iterations '
                  WRITE (iwarn, * )
               ENDIF
               IF (ihist.ne.0) then
                  WRITE (ihist, * )
                  WRITE (ihist, * ) ' WARNING GMRES :'
      WRITE (ihist,  * ) '       No convergence after '
                  WRITE (ihist, * ) iterOut * m + jH, ' iterations '
                  WRITE (ihist, * )
               ENDIF
               info (1) = - 4
            ENDIF
            IF (ihist.ne.0) then
               WRITE (ihist, 1010) rinfo (1)
               WRITE (ihist, 1011) rinfo (2)
 1010 FORMAT('B.E. on the preconditioned system:   ',E8.2)
 1011 FORMAT('B.E. on the unpreconditioned system: ',E8.2)
            ENDIF
            info (2) = iterOut * m + jH
            IF (ihist.ne.0) then
               WRITE (ihist, '(A10,I2)') 'info(1) = ', info (1)
      WRITE (ihist, '(A32,I5)') 'Number of iterations (info(2)): ', info&
     & (2)
            ENDIF
            irc (1) = 0
            retlbl = 0
            RETURN
         ENDIF
      ELSE
! Save the backward error on a file if convergence history requested
         IF (ihist.ne.0) then
            WRITE (ihist, 1001) iterOut * m + jH, bea
 1001 FORMAT(I5,11x,E8.2,7x,'--')
         ENDIF

         if (present(callback)) then
            if (associated(callback%backward_error)) then
                call callback%backward_error(iterOut*m+jH, bea, be)
            end if
         end if
!
      ENDIF
!
      jH = jH + 1
      IF (jH.le.m) then
         GOTO 10
      ENDIF
!
      iterOut = iterOut + 1
!
! we have completed the Krylov space construction, we restart if
! we have not yet exceeded the maximum number of iterations allowed.
!
      IF ( (sPa.eq.DZRO) .and. (bea.gt.cntl (1) ) ) then
!
! Compute the solution of the current linear least squares problem
!
         jH = jH - 1
         CALL zcopy (jH, H (1, m + 1), 1, yCurrent, 1)
         CALL ztrsv ('U', 'N', 'N', jH, H, m + 1, yCurrent, 1)
!
! Compute the value of the new iterate
!
         CALL zgemv ('N', n, jH, ONE, v, n, yCurrent, 1, ZERO, xCurrent,&
         1)
!
         IF ( (typePrec.eq.rightPrec) .or. (typePrec.eq.dblePrec) )     &
         then
!
!         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!
            irc (1) = precondRight
            irc (2) = xcuptr
            irc (4) = r0ptr
            retlbl = 52
            RETURN
         ELSE
            CALL zcopy (n, xCurrent, 1, r0, 1)
         ENDIF
      ENDIF
   52 CONTINUE
      IF ( (sPa.eq.DZRO) .and. (bea.gt.cntl (1) ) ) then
! Update the current solution
         CALL zcopy (n, x, 1, xCurrent, 1)
         CALL zaxpy (n, ONE, r0, 1, xCurrent, 1)
      ENDIF
!
      CALL zcopy (n, xCurrent, 1, x, 1)
!
      IF (compRsd.eq.1) then
!
! Compute the true residual
!
         CALL zcopy (n, x, 1, w, 1)
         irc (1) = matvec
         irc (2) = wptr
         irc (4) = r0ptr
         retlbl = 61
         RETURN
      ENDIF
   61 CONTINUE
      IF (compRsd.eq.1) then
         DO j = 1, n
         r0 (j) = b (j) - r0 (j)
         enddo
!
! Precondition the new residual if necessary
!
         IF ( (typePrec.eq.leftPrec) .or. (typePrec.eq.dblePrec) ) then
!
!      MY = X : w <-- M_1^{-1} r0
!
            irc (1) = precondLeft
            irc (2) = r0ptr
            irc (4) = wptr
            retlbl = 66
            RETURN
         ELSE
            CALL zcopy (n, r0, 1, w, 1)
         ENDIF
      ENDIF
   66 CONTINUE
!
!           beta = dznrm2(n,w,1)
!
      IF (compRsd.eq.1) then
         irc (1) = prosca
         irc (2) = wptr
         irc (3) = wptr
         irc (4) = r0ptr
         irc (5) = 1
         retlbl = 68
         RETURN
      ENDIF
   68 CONTINUE
      IF (compRsd.eq.1) then
         beta = dsqrt (dreal (r0 (1) ) )
!
      ELSE
! Use recurrence to approximate the residual at restart
         beta = abs (H (m + 1, m + 1) )
! Apply the Givens rotation is the reverse order
         DO j = m, 1, - 1
         H (j, m + 1) = ZERO
         temp = dreal (rotCos (j) )
         CALL zrot (1, H (j, m + 1), 1, H (j + 1, m + 1), 1, temp,      &
         - rotSin (j) )
         enddo
!
! On applique les vecteurs V
!
         CALL zgemv ('N', n, m + 1, ONE, v, n, H (1, m + 1) , 1, ZERO,  &
         w, 1)
!
      ENDIF
      DO j = 1, n
      V (j, 1) = ZERO
      enddo
      aux = ONE / beta
      CALL zaxpy (n, aux, w, 1, V (1, 1), 1)
!
      GOTO 7
!
      END SUBROUTINE zgmres
!
!
      SUBROUTINE init_zgmres (icntl, cntl)
!
!  Purpose
!  =======
!    Set default values for the parameters defining the characteristics
! of the Gmres algorithm.
!  See the User's Guide for an example of use.
!
!
! Written : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS
!
!
!  Arguments
!  =========
!
! icntl    (input) INTEGER array. length 6
!            icntl(1) : stdout for error messages
!            icntl(2) : stdout for warnings
!            icntl(3) : stdout for convergence history
!            icntl(4) : 0 - no preconditioning
!                       1 - left preconditioning
!                       2 - right preconditioning
!                       3 - double side preconditioning
!                       4 - error, default set in Init
!            icntl(5) : 0 - modified Gram-Schmidt
!                       1 - iterative modified Gram-Schmidt
!                       2 - classical Gram-Schmidt
!                       3 - iterative classical Gram-Schmidt
!            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                       1 - user supplied initial guess
!            icntl(7) : maximum number of iterations
!            icntl(8) : 1 - default compute the true residual at each re
!                       0 - use recurence formaula at restart
!
! cntl     (input) double precision array, length 5
!            cntl(1) : tolerance for convergence
!            cntl(2) : scaling factor for normwise perturbation on A
!            cntl(3) : scaling factor for normwise perturbation on b
!            cntl(4) : scaling factor for normwise perturbation on the
!                      preconditioned matrix
!            cntl(5) : scaling factor for normwise perturbation on
!                      preconditioned right hand side
!
! Output variables
! ----------------
      INTEGER icntl ( * )
      DOUBLEPRECISION cntl ( * )
!
      icntl (1) = 6
      icntl (2) = 6
      icntl (3) = 0
      icntl (4) = 4
      icntl (5) = 0
      icntl (6) = 0
      icntl (7) = - 1
      icntl (8) = 1
!
      cntl (1) = 1.0d-5
      cntl (2) = 0.0d0
      cntl (3) = 0.0d0
      cntl (4) = 0.0d0
      cntl (5) = 0.0d0
!
      RETURN
      END SUBROUTINE init_zgmres

end module zPackgmres
