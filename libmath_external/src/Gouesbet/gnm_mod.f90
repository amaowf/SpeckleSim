! Programme de calcul des gnm
!
!    issu de gnm_test. for de Gerard GREHAN
!    Version prise au 7 avril 2010
MODULE gnm_mod
    INTEGER, PARAMETER :: MAX_M = 10
    INTEGER, PARAMETER :: MAX_N    = 1000
    INTEGER, PARAMETER :: DBL = 8 ! Double prcision
    COMPLEX (KIND=DBL), DIMENSION(MAX_N, -MAX_M:MAX_M) :: gtm, gte
    ! enables logging messages
    LOGICAL :: writelog = .TRUE.
    integer :: unitlog = 13

    contains

        !
        ! subroutine gnmf !
        ! Computes g(n,m) !
        ! References :
        !
        ! History :
        !
        ! translated in Fortran 95 7/04/2010 !
        ! Parameters:
        !     x0ad, y0ad, z0ad : DOUBLE PRECISION,
        !         Dimensionless coordinates of the beam waist center viewed
        !         from the particle.
        !     nfing : INTEGER, maximum n of g(n,m) to compute
        !     S : DOUBLE PRECISION,
        !         beam shape factor, s=w0/l=1/k.w0
        !         W0 : beam waist radius
        !         K : wave-number.
        !     esj : DOUBLE PRECISION,
        !        accuracy of g(n,m).
        !         The program stops when the term to add is smaller than esj.
        !
        ! Results :
        !     gtm : COMPLEX, g(n,m) for the tm wave.
        !     gte : COMPLEX, g(n,m) for the te wave.
        !
        ! THIS SUBROUTINE DOES NOT CALL ANY OTHER SUBROUTINE.
        !
        SUBROUTINE gnmf(nfing,s,x0ad,y0ad,z0ad,esj)
            IMPLICIT NONE
            INTEGER :: nfing
            REAL (KIND=DBL) :: s,x0ad,y0ad,z0ad,esj
            ! local variables
            INTEGER :: n,j,m,mm
            REAL (KIND=DBL) :: en,ej,em,r1
            COMPLEX (KIND=DBL) :: zi,a,b,c,d,am,bm
            COMPLEX (KIND=DBL) :: iqronw0,iqw2,xyn,xyp,iqbar,ht0
            COMPLEX (KIND=DBL), DIMENSION(-MAX_M:MAX_M) :: htm, hte
            IF(writelog)WRITE(unitlog,*) 'In gnmf ',nfing

            zi = CMPLX(0.0_DBL,1.0_DBL, DBL)
            xyn= CMPLX(x0ad, -y0ad, DBL)
            xyp= CMPLX(x0ad, y0ad, DBL)
            iqbar = 1.0_DBL/(1.0_DBL + 2.0_DBL*zi*z0ad)
            !
            DO n=1,nfing
                en = REAL(n,KIND=DBL)
                iqronw0 = iqbar*(en + 0.5_DBL)*S
                iqw2 = iqronw0*iqronw0
                !
                ! Compute htm(m) and hte(m)
                !
                ht0 = iqronw0
                !
                ! Compute htn(0) and hte(0)
                !
                a = iqronw0
                !
                ! j loop
                !
                j = 0
                DO
                    j = j+1
                    ej= REAL(j*(j+1),KIND=DBL)
                    a = a*iqw2*xyn*xyp/ej
                    ht0 = ht0+a
                    IF( ABS(a)<=esj ) EXIT
                ENDDO
                !
                ! m>0 and m<0
                !
                a = CMPLX(1.0_DBL,0.0_DBL,KIND=DBL)
                b = a
                c = a
                mm = MIN(n,MAX_M)
                DO m=1,mm
                    EM=REAL(m,KIND=DBL)
                    htm(m) = a*b
                    hte(m) = htm(m)
                    htm(-m)= a*c
                    hte(-m)= -htm(-m)
                    a = a*iqronw0
                    b = b*xyn/em
                    c = c*xyp/em
                ENDDO
                !
                ! Compute htm(m) and hte(m)
                !
                am = iqronw0
                bm = iqronw0
                d = iqw2*xyn*xyp
                !
                DO m=1,mm
                    am = am*iqronw0*xyn/REAL(m,KIND=DBL)
                    htm(m) = htm(m)+am*(xyp+xyn/REAL(m+1,KIND=DBL))
                    hte(m) = hte(m)+am*(xyp-xyn/REAL(m+1,KIND=DBL))
                    bm = bm*iqronw0*xyp/REAL(m,KIND=DBL)
                    htm(-m)= htm(-m)+bm*(xyp/REAL(m+1,KIND=DBL)+xyn)
                    hte(-m)= hte(-m)+bm*(xyp/REAL(m+1,KIND=DBL)-xyn)
                    a = am
                    b = bm
                    j = m
                    ! J loop
                    DO
                        j = j+1
                        c = d/REAL((j)*(j-m))
                        a = a*c
                        htm(m) = htm(m)+a*(xyp/REAL(j-m+1,KIND=DBL) + xyn/REAL(j+1,KIND=DBL))
                        hte(m) = hte(m)+a*(xyp/REAL(j-m+1,KIND=DBL) - xyn/REAL(j+1,KIND=DBL))
                        b = b*c
                        htm(-m)= htm(-m)+b*(xyp/REAL(j+1,KIND=DBL) + xyn/REAL(j-m+1,KIND=DBL))
                        hte(-m)= hte(-m)+b*(xyp/REAL(j+1,KIND=DBL) - xyn/REAL(j-m+1,KIND=DBL))
                        IF((ABS(a)<=esj).AND.(ABS(b)<=esj)) EXIT
                    ENDDO
                ENDDO
                !
                ! Compute gtm(n,m) and gte(n,m)
                !
                r1 = 2.0_DBL/(2.0_DBL*en+1.0_DBL)
                !
                a = zi*z0ad/(s*s)-iqbar*(x0ad*x0ad + y0ad*y0ad)
                a = 0.5_DBL*iqbar*EXP(a-iqronw0*iqronw0/iqbar)
                !
                gtm(n,0) = r1*en*(en+1.0_DBL)*a*htm(0)*zi
                gte(n,0) = r1*en*(en+1.0_DBL)*a*hte(0)
                c = a
                DO m=1,mm
                    gtm(n,m) = c*htm(m)
                    gte(n,m) = c*hte(m)*(-zi)
                    gtm(n,-m)= c*htm(-m)
                    gte(n,-m)= c*hte(-m)*(-zi)
                    c = c*r1*(-zi)
                ENDDO
            ENDDO
        END SUBROUTINE gnmf
END MODULE gnm_mod
