************************************************************************
* ******************************************************************** *
* *                                                                  * *
* *   MM      MM     AA     TTTTTTTT  HH    HH             QQQQQQ    * *
* *   MMM    MMM    A  A       TT     HH    HH            QQ    QQ   * *
* *   MM M  M MM   A    A      TT     HH    HH            QQ    QQ   * *
* *   MM  MM  MM  AAAAAAAA     TT     HHHHHHHH            QQ    QQ   * *
* *   MM      MM  AA    AA     TT     HH    HH            QQ    QQ   * *
* *   MM      MM  AA    AA     TT     HH    HH            QQ  Q QQ   * *
* *   MM      MM  AA    AA     TT     HH    HH  =========  QQQQQQ    * *
* *                                                             QQ   * *
* *           General mathematical procedures: QUADRATURES           * *
* *                                                                  * *
* ******************************************************************** *
*                                                                      *
*  ##################################################################  *
*  #                                                                #  *
*  #     Vadim Naumov, INFN, Sezione di Firenze and LTP, ISU        #  *
*  #                                                                #  *
*  # Version of August 15, 1997  for MS 32-bit FORTRAN PowerStation #  *
*  #                                                                #  *
*  ##################################################################  *
*                                                                      *
* ******************************************************************** *
* *                                                                  * *
* *   MM      MM            ll          -------------------------    * *
* *   MMM    MMM            ll         | ADAPTIVE QUADRATURE FOR |   * *
* *   MM M  M MM  uu    uu  ll         | MULTIPLE INTEGRALS OVER |   * *
* *   MM  MM  MM  uu    uu  ll         |   RECTANGULAR REGIONS   |   * *
* *   MM      MM  uu    uu  ll         | by Genz-Malik algorithm |   * *
* *   MM      MM  uu    uu  ll         |    (revized version)    |   * *
* *   MM      MM    uuuuuu   llllll     -------------------------    * *
* *                                                                  * *
* *                                                                  * *
* *   ------------------------------------------------------------   * *
* *                                                                  * *
* *                       T H E    S O U R C E                       * *
* *                                                                  * *
* *    --------                                            ------    * *
* *   | RADMUL |          CERN Program Library            | D120 |   * *
* *    --------                                            ------    * *
* *                                                                  * *
* *   Authors: A.C. Genz & A.A. Malik        Library: MATHLIB        * *
* *   Submitter: K.S. Kolbig                 Submitted: 15.11.1995   * *
* *                                                                  * *
* *   References:                                                    * *
* *                                                                  * *
* *   [1] A.C.Genz & A.A.Malik, Remarks on algorithm 006: An adap-   * *
* *       tive algorithm for numerical integration over N-dimensi-   * *
* *       onal rectangular region, J. Comput. Appl. Math. 6 (1980)   * *
* *       295-302.                                                   * *
* *   [2] A. van Doren & L. de Ridder,  An adaptive algorithm  for   * *
* *       numerical integration over N-dimensional cube, J.Comput.   * *
* *       Appl. Math. 2 (1976) 207-217.                              * *
* *                                                                  * *
* *   ------------------------------------------------------------   * *
* *   E N T R I E S:                                                 * *
* *   MulSet(Fun,Res,RelErr,MinCal,IntDim,*)                         * *
* *   MulTot(Fun,Res,RelErr,MinCal,*)                                * *
* *   MulInt(Fun,Res,*)                                              * *
* *   MulInf                                                         * *
* *   ------------------------------------------------------------   * *
* *                                                                  * *
* ******************************************************************** *
************************************************************************
      SUBROUTINE MuLSet(Fun,Res,RelErr,MinCal,IntDim,*)
************************************************************************

        IMPLICIT REAL (A-H,O-Z),INTEGER (I-N)

            SAVE

              LOGICAL(2) Flag

*     THE FOLLOWING PARAMETERS  MUST  BE  SPECIFIED  BEFORE  USING!
         INTEGER, PARAMETER ::
     ,            Nlog  =  10,
     ,            Length= 100000000,                                     > (2*MaxDim+3)(1+k)/2, k > 0
     ,            MaxDim=   3                                            > 1 [MAX(2,IntDim)]
            REAL, PARAMETER ::
     ,            Zero  =    0.0,
     ,            One   =    1.0,
     ,            Half  =    1.0/2,
     ,            cVN   =  epsilon(one)*1.0d+1,
     ,            V2    =  245.0/486,
     ,            V4    =   25.0/729,
     ,            W2    =  980.0/6561,
     ,            W4    =  200.0/19683,
     ,            C0    =    4.0/19683,
     ,            C1    =    1.0/729,
     ,            C3    =    5.0/1458,
     ,            C5    = 6859.0/19683,
c    ,            Y2    =    0.3585685828003181,                         3/SQRT(70)
c    ,            Y4    =    0.9486832980505138,                         3/SQRT(10)
c    ,            Y5    =    0.6882472016116853                          3/SQRT(19)
     ,            Y2    =    0.358568582800318091990645153907938d+00,    3/SQRT(70)
     ,            Y4    =    0.948683298050513799599668063329816d+00,    3/SQRT(10)
     ,            Y5    =    0.688247201611685297721628734293623d+00     3/SQRT(19)

      ALLOCATABLE STORE(:)                                               storage array
      DIMENSION CTR(MaxDim),WTH(MaxDim),VTH(MaxDim),X(MaxDim)

         COMMON     /MulLim/Xlow(MaxDim),Xupp(MaxDim)                    MuL integration limits
*     the arrays are lower limits/upper limits of integration
      EQUIVALENCE (CTR1,CTR(1)),(WTH1,WTH(1)),(WTH2,WTH(2))

*        THIS  FRAGMENT  MAY BE REMOVED AS  A UNIT  WHEN DEBUGGED
         IF (IntDim.LT.1.OR.IntDim.GT.MaxDim.OR.MinCal.LT.1) THEN
           WRITE(Nlog,*) ' '
           WRITE(Nlog,*) '   ------------------------ '
           WRITE(Nlog,*) '  | MISTAKE AT CALL MULSET |'
           WRITE(Nlog,*) '   ------------------------ '
           WRITE(Nlog,*) ' '
           WRITE(Nlog,*) '    MaxDim =',MaxDim
           WRITE(Nlog,*) '    IntDim =',IntDim
           WRITE(Nlog,*) '    MinCal =',MinCal
           WRITE(Nlog,*) ' '
           STOP          '    MISTAKE AT CALL MULSET  '
      endIF
         EPS=RelErr                                                         "relative integration accuracy"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  COMMENTS # The authors of the method credited the parameter EPS as  c
c             being the "specified relative accuracy of integration".  c
c             Unfortunately, this is not the case,  especially in the  c
c             original REAL*4 version for VAX. Generally this parame-  c
c             ter (together with parameters Length and MinCal) provi-  c
c             des a control upon the accuracy and a means for minimi-  c
c             zation of CPU time and memory.  A good rule of thumb is  c
c             "THE REAL INTEGRATION ACCURACY RISES AS EPS DECREASES".  c
c             Usually, for "not-too-terrible" integrands the real ac-  c
c             curacy is BETTER than EPS but from a certain (dependent  c
c             of the integral) value,  further decrease of EPS leaves  c
c             the result invariant. It is therefore strongly recomme-  c
c             nded that  FOR EVERY NEW INTEGRAND AND/OR EVERY NEW SET  c
c             OF THE LIMITS OF INTEGRATION, THE INTEGRATION PROCEDURE  c
c             BE REPEATED WITH DIFFERENT EPS.                          c
c                                                                      c
c           # Use entry 'MulInt' for any series of integrals with the  c
c             common parameters RelErr, MinCal and IntDim.             c
c                                                                      c
c           # Use entry 'MulTot' if, with the same IntDim, the RelErr  c
c             and/or MinCal must be changed.                           c
c                                                                      c
c           # When the parameters  Length  [the length of array STORE  c
c             of working storage],  MinCal [the minimum number of in-  c
c             tegrand Fun calls requested to be allowed] and RelErr =  c
c             EPS are fitted closely,  the VAX REAL*16 version of the  c
c             routine gives very good results  (accuracy-time-memory)  c
c             for integrands with variation up to (and including) 150  c
c             ORDERS OF MAGNITUDE!                                     c
c                                                                      c
c           # Maximum multiplicity is formally unbounded but actually  c
c             it is limited by your computer. For instance, ALPHA-VAX  c
c             usually has no troubles with 20-fold integrals.          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        Mult=IntDim                  ! Multiplicity of integral.
           N=MAX(2,Mult)             ! "Pseudo-multiplicity" at Mult=1.
        Jump=MIN(2,Mult)             ! Switch 'Mult = 1 <--> Mult > 1'.
          N0=2*N+3                   !
          N1=N0-1                    !
          N2=N0-2                    ! Ancillari constants.
      N0dupl=2*N0                    !
      Cardin=2**N                    !
          V1=(729+(50*N-950)*N)*C1   ! V1,...,V4  and  W1,...,W5 are the
          V3=(53-20*N)*C3            ! weights for two integration rules
          W1=(3206+(100*N-2280)*N)   !
     *       *C0                     ! B'F and BF of degree of exactness
          W3=(455-100*N)*C0          ! 5 and 7, respectively, normalized
          W5=C5/Cardin               ! to the factor 2**N (see Ref.[1])
      NRcall=2**N+2*N*(N+1)+1        ! Number of rule B'F calls.
      MinMax=NRcall*                 !
     *       INT(2*One*Length/N0-1)  ! MIN(MaxPTS) for the given Length.
      MinPTS=MinCal                  ! The minimum number of  Fun calls;
                                     ! simultaneously, MinPTS prescribes
      MaxPTS=MinMax+MinPTS           ! the maximum number of  Fun calls.
      MinFin=MaxPTS                  ! Initial value for 'MulInf'
      MaxFin=MinPTS                  ! Initial value for 'MulInf'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  If parameters MaxPTS or Length proved to be inadequate to the req-  c
c  uested accuracy RelErr, the routine prints an error message and is  c
c  interrupted with alternative RETURN.                                c
c                                                                      c
c  MinFin and MaxFin contain on exit  the minimum and maximum numbers  c
c  of integrand evaluations performed within a given series of integ-  c
c  rations. To print these numbers it will suffice to insert operator  c
c  'CALL MULINF' in any point after that but before the next series.   c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (Mult.EQ.1) NRcall=0
      Lim1=MaxPTS-2*NRcall
      Lim2=Length-N0

      RETURN

*     ==================================================================
      ENTRY MulTot(Fun,Res,RelErr,MinCal,*)
*     ==================================================================
         EPS=RelErr
        Lim1=Lim1-MinPTS+MinCal
      MinPTS=MinCal
      MaxPTS=MinMax+MinPTS

*     ==================================================================
      ENTRY MulInt(Fun,Res,*)
*     ==================================================================
        Flag=.FALSE.
      NFcall=0
       NsubR=N0
       LsubR=N0
      Output=Zero
      AbsErr=Zero
      ALLOCATE(STORE(Length))
      IF (Mult.GT.1) GOTO 2
      CTR1=(Xupp(1)+Xlow(1))*Half
      WTH1= Xupp(1)-CTR1
      WTH2= Half
    1 VolRGN=4*WTH1*WTH2
      X(1)=CTR1
      Sum1=Fun(X)
       V=WTH1*Y2
      X(1)=CTR1-V
      Sum2=Fun(X)
      X(1)=CTR1+V
      Sum2=Sum2+Fun(X)
       V=WTH1*Y4
      X(1)=CTR1-V
      Sum4=Fun(X)
      X(1)=CTR1+V
      Sum4=Sum4+Fun(X)
          IF (7*Sum2.EQ.Sum4+12*Sum1) THEN
            Ndiv=2
                                      ELSE
            Ndiv=1
      endIF
         V=2*Sum1
      Sum2=V+Sum2
      Sum3=V+Sum4
      Sum4=2*Sum4
         IF (WTH1.EQ.Zero) THEN
           NFcall=NFcall+5
           IF (WTH2.GT.Zero) THEN
             Sum5=Sum1
                             ELSE
             Sum5=V
        endIF
           GOTO 13
      endIF
       V=WTH1*Y5
      X(1)=CTR1-V
      Sum5=Fun(X)
         IF (WTH1.GT.Zero) THEN
           X(1)=CTR1+V
           NFcall=NFcall+7
           IF (WTH2.GT.Zero) THEN
             Sum5=2*(Sum5+Fun(X))
                             ELSE
             Sum5=Sum5+Fun(X)
        endIF
                           ELSE
           IF (WTH2.GT.Zero) THEN
             X(1)=CTR1+V
             Sum5=2*Sum5+Fun(X)
             NFcall=NFcall+7
                             ELSE
             NFcall=NFcall+6
        endIF
      endIF
      GOTO 13
    2 DO 3 k=1,N
      CTR(k)=(Xupp(k)+Xlow(k))*Half
    3 WTH(k)= Xupp(k)-CTR(k)
    4 VolRGN=Cardin
      DO 5 k=1,N
      VolRGN=VolRGN*WTH(k)
    5 X(k)=CTR(k)
      Sum1=Fun(X)
      Sum2=Zero
      Sum3=Zero
      DifMax=Zero
      DO 6 k=1,N
      V=CTR(k)
      VTH(k)=WTH(k)*Y2
      X(k)=V-VTH(k)
      S2=Fun(X)
      X(k)=V+VTH(k)
      S2=S2+Fun(X)
      Sum2=Sum2+S2
      VTH(k)=WTH(k)*Y4
      X(k)=V-VTH(k)
      S3=Fun(X)
      X(k)=V+VTH(k)
      S3=S3+Fun(X)
      Sum3=Sum3+S3
      Dif=ABS(7*S2-S3-12*Sum1)
      DifMax=MAX(Dif,DifMax)
      IF (DifMax.EQ.Dif) Ndiv=k
    6 X(k)=V
      Sum4=Zero
      DO 9 j0=2,N
      j=j0-1
      DO 8 k=j0,N
      DO 7 j1=1,2
      VTH(j)=    -VTH(j)
      X(j)=CTR(j)+VTH(j)
      DO 7 k1=1,2
      VTH(k)=    -VTH(k)
      X(k)=CTR(k)+VTH(k)
    7 Sum4=Sum4+Fun(X)
    8 X(k)=CTR(k)
    9 X(j)=CTR(j)
      Sum5=Zero
      DO 10 j=1,N
      VTH(j)=    -WTH(j)*Y5
   10 X(j)=CTR(j)+VTH(j)
   11 Sum5=Sum5+Fun(X)
      DO 12 j=1,N
      VTH(j)=    -VTH(j)
      X(j)=CTR(j)+VTH(j)
   12 IF (VTH(j).GT.Zero) GOTO 11
      NFcall=NFcall+NRcall
   13 ErrRGN=VolRGN*(V1*Sum1+V2*Sum2+V3*Sum3+V4*Sum4)
      ValRGN=VolRGN*(W1*Sum1+W2*Sum2+W3*Sum3+W4*Sum4+W5*Sum5)
      ErrRGN=ABS(ErrRGN-ValRGN)
      Output=    Output+ValRGN
      AbsErr=    AbsErr+ErrRGN
         IF (Flag) THEN
   14      IsubR=2*NsubR
           IF (IsubR.GT.LsubR) GOTO 16
           IF (IsubR.LT.LsubR) THEN
             I=IsubR+N0
             IF (STORE(IsubR).LT.STORE(I)) IsubR=I
        endIF
           IF (ErrRGN.GE.STORE(IsubR)) GOTO 16
           DO k=0,N1
             STORE(NsubR-k)=STORE(IsubR-k)
        endDO
         NsubR=IsubR
         GOTO 14
      endIF
   15 IsubR=(NsubR/N0dupl)*N0
         IF (IsubR.GE.N0) THEN
           IF (ErrRGN.GT.STORE(IsubR)) THEN
c          IF ((IsubR.GE.N0).AND.(ErrRGN.GT.STORE(IsubR))) THEN
             DO k=0,N1
               STORE(NsubR-k)=STORE(IsubR-k)
          endDO
             NsubR=IsubR
             GOTO 15
        endIF
      endIF
   16 I=NsubR-2
      DO 17 j=1,N
      IsubR=I-(j+j)
      STORE(IsubR  )=WTH(j)
   17 STORE(IsubR+1)=CTR(j)
      STORE(NsubR  )=ErrRGN
      STORE(NsubR-1)=ValRGN
      STORE(I)      =Ndiv
         IF (Flag) THEN
           Flag=.FALSE.
           CTR(Nd)=CTR(Nd)+2*WTH(Nd)
           LsubR=LsubR+N0
           NsubR=LsubR
           GOTO (1,4) Jump
      endIF
         IF ((AbsErr.LE.Eps*ABS(Output).OR.ABS(Output).LE.cVN)
     #       .AND.NFcall.GE.MinPTS) THEN
           DEALLOCATE(STORE)
           MinFin=MIN(MinFin,NFcall)
           MaxFin=MAX(MaxFin,NFcall)
           Res=Output                                                    Tutto va bene
           RETURN
      endIF
         IF (NFcall.GT.Lim1) THEN
           WRITE(*   ,101) Mult,MaxPTS                                   Error message (screen)
           WRITE(Nlog,101) Mult,MaxPTS                                   Error message (LogFile)
           GOTO 18
      endIF
         IF (LsubR.GT.Lim2) THEN
           WRITE(*   ,102) Mult,Length                                   Error message (screen)
           WRITE(Nlog,102) Mult,Length                                   Error message (LogFile)
           GOTO 18
      endIF
      Flag=.TRUE.
      NsubR=N0
         DO j=1,N
           IsubR=N2-(j+j)
           WTH(j)=    STORE(IsubR  )
           CTR(j)=    STORE(IsubR+1)
      endDO
      AbsErr=AbsErr-STORE(N0)
      Output=Output-STORE(N1)
      Nd    =       STORE(N2)
      WTH(Nd)=   Half*WTH(Nd)
      CTR(Nd)=CTR(Nd)-WTH(Nd)
      GOTO (1,4) Jump
   18 CONTINUE
      DEALLOCATE(STORE)
         IF (Output.NE.Zero) THEN
         WRITE(Nlog,103) AbsErr/ABS(Output),NFcall
                             ELSE
         WRITE(Nlog,104) AbsErr,NFcall
      endIF
      RETURN 1

*     ==================================================================
      ENTRY MulInf
*     ==================================================================

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  This entry together with all starred operators (see above)  may be  c
c  removed if the output (see FORMATs) is not warranted. However that  c
c  information is helpful in the debug stage, specifically for reduc-  c
c  tion of CPU time and volume occupied in the main memory.            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      WRITE(*   ,105) Mult,EPS,MaxPTS,MinPTS,Length,MaxFin,MinFin
      WRITE(Nlog,105) Mult,EPS,MaxPTS,MinPTS,Length,MaxFin,MinFin
      MinFin=MaxPTS
      MaxFin=MinPTS
      RETURN

  101 FORMAT(/ '  MULT =',I3/'  THE MAXIMUM NUMBER OF INTEGRAND CALLS,',
     #     I15,', IS TOO SMALL'/'  FOR THE REQUIRED ACCURACY:')
  102 FORMAT(/ '  MULT =',I3/'  THE STORAGE ARRAY LENGTH,',
     #     I10,', IS TOO SMALL'/'  FOR THE REQUIRED ACCURACY:')
  103 FORMAT('  RELATIVE ERROR =',1pd9.2,
     #                            ', NUMBER OF INTEGRAND CALLS =',I14/)
  104 FORMAT('  RESULT = 0! ABSOLUTE ERROR =',1pd9.2,
     #                            ', NUMBER OF INTEGRAND CALLS =',I14/)
  105 FORMAT(/16x,'MULINFORM (multiplicity of integral is',I4,')'/
     #         3x,69('-')/'  | RelErr =',1PD10.3,' | MaxCal =',I15,
     #                     ' | MinCal =',    I12,' |'/
     #                    '  | Length =',    I10,' | MaxFin =',I15,
     #                     ' | MinFin =',    I12,' |'/3x,69('-')/)
      END SUBROUTINE MulSet
c  =======================================================================
c  |             MULINFORM (multiplicity of integral is   2)             |
c  |---------------------------------------------------------------------|
c  | RelErr = 1.000D-05 | MaxCal =        9714359 | MinCal =         100 |
c  | Length =   2000000 | MaxFin =           6341 | MinFin =        1377 |
c  =======================================================================