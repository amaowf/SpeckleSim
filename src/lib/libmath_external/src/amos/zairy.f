﻿      SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
C***BEGIN PROLOGUE  ZAIRY
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
C         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
C         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
C         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
C         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
C         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
C
C         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
C         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
C         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
C         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI)
C           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             AI=AI(Z)                ON ID=0 OR
C                             AI=DAI(Z)/DZ            ON ID=1
C                        = 2  RETURNS
C                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
C                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
C                             ZTA=(2/3)*Z*CSQRT(Z)
C
C         OUTPUT     AIR,AII ARE DOUBLE PRECISION
C           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
C                    KODE
C           NZ     - UNDERFLOW INDICATOR
C                    NZ= 0   , NORMAL RETURN
C                    NZ= 1   , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN
C                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
C                            TOO LARGE ON KODE=1
C                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
C         FUNCTIONS BY
C
C            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
C                           C=1.0/(PI*SQRT(3.0))
C                            ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
C         MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 1986
C
C***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,I1MACH,D1MACH
C***END PROLOGUE  ZAIRY
C     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB,
     * ZABS_amos
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CYR(1), CYI(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667D-01,
     * 3.55028053887817240D-01,2.58819403792806799D-01,
     * 1.83776298473930683D-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS_amos(ZR,ZI)
      TOL = DMAX1(D1MACH(4),1.0D-18)
      FID = DBLE(FLOAT(ID))
      IF (AZ.GT.1.0D0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR CABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 170
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = DMIN1(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0D0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR CABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      ALAZ = DLOG(AZ)
C--------------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA=DMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -DABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0) GO TO 90
      IF (ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 100
      AA = -AA + 0.25D0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 270
  100 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     * ELIM, ALIM)
      IF (NN.LT.0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      IF (KODE.EQ.2) GO TO 120
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 120
      AA = -AA - 0.25D0*ALAZ
      IFLAG = 2
      SFAC = 1.0D0/TOL
      IF (AA.LT.(-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     * ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      IF (IFLAG.NE.0) GO TO 150
      IF (ID.EQ.1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      IF (ID.EQ.1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0D+3*D1MACH(1)
      S1R = ZEROR
      S1I = ZEROI
      IF (ID.EQ.1) GO TO 190
      IF (AZ.LE.AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0D0
      AA = DSQRT(AA)
      IF (AZ.LE.AA) GO TO 200
      S1R = 0.5D0*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      IF(NN.EQ.(-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
