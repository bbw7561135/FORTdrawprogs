************************************************************************
* ******************************************************************** *
* *                                                                  * *
* *   MM      MM     AA     TTTTTTTT  HH    HH             SSSSSS    * *
* *   MMM    MMM    A  A       TT     HH    HH            SS    SS   * *
* *   MM M  M MM   A    A      TT     HH    HH            SS         * *
* *   MM  MM  MM  AAAAAAAA     TT     HHHHHHHH              SSSS     * *
* *   MM      MM  AA    AA     TT     HH    HH                  SS   * *
* *   MM      MM  AA    AA     TT     HH    HH            SS    SS   * *
* *   MM      MM  AA    AA     TT     HH    HH  =========  SSSSS     * *
* *                                                                  * *
* *       General mathematical procedures: SPLINE FUNCTIONS          * *
* *                                                                  * *
* *       =======  ONE-DIMENSIONAL SPLINE FUNCTIONS  ======          * *
* *                                                                  * *
* ******************************************************************** *
*  ##################################################################  *
*  #                                                                #  *
*  #     Vadim Naumov, INFN, Sezione di Firenze and LTP, ISU        #  *
*  #                                                                #  *
*  # Version of August 15, 1997  for MS 32-bit FORTRAN PowerStation #  *
*  #                                                                #  *
*  ##################################################################  *
************************************************************************
************************************************************************
      SUBROUTINE Coeff1(Mult,Mode,Tie,Eps,Issue,NX,Xmin,Xmax,F,C,Quiz,L)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)
            LOGICAL(2) Tie,Quiz
                  REAL Inter1

         INTEGER, PARAMETER ::
     #            NXmax=999999,                                          User's setting (NXmax > 2)
     #            Nlog =    10                                           LogFile (must be opened in MAIN)
            REAL, PARAMETER ::
     #            Zero=0, One=1, Const=1.0

         DIMENSION F(*),C(*)
         ALLOCATABLE G(:)                                                storage array
*     ---------------------------------------------------------------- *
*     ATTENTION !   IF it is impossible to use ALLOCATABLE  dimen-     *
*                   sion G(:) THEN comment lines 39, 119, 151  AND     *
*                   uncomment lines 45, 46, 47, 48.                    *
*                                                                      *
*     DIMENSION G(NXmax+2)                                               Storage array
*     INCLUDE 'E:/FORTRAN/Sources/Shared/Include/Heap.fi'                V.A. Naumov
*     INCLUDE 'E:/FORTRAN/Sources/Shared/Include/Heap.fi'                S.I. Sinegovsky
*     INCLUDE 'G:/FORTRAN/Sources/Shared/Include/Heap.fi'                K. Kuzmin
*     EQUIVALENCE (G(1),Heap(1))                                         Memory saving
*     ---------------------------------------------------------------- *
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c   ------------------------------------------------------------------ c
c     D E S C R I P T I O N    O F    D U M M Y    A R G U M E N T S   c
c   ------------------------------------------------------------------ c
c                                                                      c
c   Mult   is the "installation code" for installation of the entries  c
c          in function 'Inter1' and subroutine 'IRoots':               c
c                                                                      c
c          # Mult=0 for entries 'Sp1' and 'rSp1' (see 'Inter1');       c
c          # Mult=1 for entry 'Root' (see 'Inter1');                   c
c          # Mult>1 (=Kmax) for entry 'Roots' (see 'IRoots').          c
c                                                                      c
c   Mode   is the "interpolation code" which sets the method of inte-  c
c          rpolation:                                                  c
c                                                                      c
c          # Mode=1 for conventional interpolation,                    c
c          # Mode=2 for logarithmic interpolation,                     c
c          # Mode=3 for the "shifted" logarithmic interpolation        c
c          # Mode=4 for interpolation with a function 'Transf'.        c
c                                                                      c
c   The specification of the function 'Transf' may be realized by di-  c
c   fferent ways; in particular, by declaration of the desirable tra-  c
c   nsformation before the 1st imperative statement, as in the example c
c                                                                      c
      Transf(Func)=SIN(Func)    ! Example of using the Mode=3 regime   c
c                                                                      c
c   Of course,  the same must be done for the  inverse transformation  c
c   in the function 'Inter1' [see example there].                      c
c                                                                      c
c   -----------------------------------------------------------------  c
c   ATTENTION! There is no limit to the volume of supplement interpo-  c
c              lation modes but concurrently the corresponding change  c
c              must be made in the function 'Inter1'.                  c
c   -----------------------------------------------------------------  c
c                                                                      c
c   Tie    is the "procedure code"  for interpolation and/or solution  c
c          to transcendental equations:                                c
c                                                                      c
c          # Tie=T with a stringent clamping of tips;                  c
c          # Tie=F without the clamping.                               c
c                                                                      c
c   Eps    is the const to manage the clamping.                        c
c                                                                      c
c   -----------------------------------------------------------------  c
c   ATTENTION! For entry 'Roots' of the subroutine 'IRoots' the para-  c
c              meters Tie and Eps play slightly another part.          c
c   -----------------------------------------------------------------  c
c                                                                      c
c   Issue  numerates different sets of interpolation ranges;           c
c   NX     is the number of the interpolation nodes [NX.GE.3];         c
c   Xmin   is the lower boundary of the interpolation range;           c
c   Xmax   is the upper boundary of the interpolation range;           c
c   F(*)   is the input array for the function values evaluated on an  c
c          evenly divided scale [*.GE.NX];                             c
c   C(*)   is the output array for coefficients of local one-dimensi-  c
c          onal B-spline [*.GE.NX+2];                                  c
c   Quiz   manages a check upon the accuracy of interpolation  in the  c
c          reference points (interpolation nodes);                     c
c   L      switchs calculating the spline coefficients (when L=0, the  c
c          routine doesn't calculate the coefficients but does install c
c          the parameters for the functions 'Inter1' or 'IRoots').     c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*     ---------------------------------------------------------------- *
*     Installation of functions 'Inter1/IRoots':
*     ---------------------------------------------------------------- *

         ALLOCATE(G(NXmax+2))

         DO m=1,NX
           IF (Mode.EQ.1) G(m)=          F(m)
           IF (Mode.EQ.2) G(m)=      LOG(F(m))
           IF (Mode.EQ.3) G(m)=LOG(Const+F(m))
           IF (Mode.EQ.4) G(m)=   Transf(F(m))
      endDO
         IF (Tie.AND.Mult.EQ.0) THEN
           G(NX+1)=F(1)
           G(NX+2)=F(NX)
      endIF
         IF (Mult.LE.1) THEN
           discrM=Inter1(Issue,G,Xmax,Mult,Tie,Eps,Mode,NX,Xmin)
                        ELSE
           CALL IRoots(G,Xmax,Mult,Tie,Eps,NX,Xmin)
      endIF
         IF (L.EQ.0) RETURN
         StepX=G(NX+1)
*     ---------------------------------------------------------------- *
*     Calculation of spline coefficients:
*     ---------------------------------------------------------------- *
         C(1)=1.4375*G(1)-1.3750*G(2)+0.4375*G(3)
         C(2)=0.4375*G(1)+0.1250*G(2)-0.0625*G(3)
         DO m=3,NX
           C(m)=0.6250*G(m-1)-0.0625*(G(m-2)+G(m))
      endDO
         C(NX+1)=0.4375*G(NX)+0.1250*G(NX-1)-0.0625*G(NX-2)
         C(NX+2)=1.4375*G(NX)-1.3750*G(NX-1)+0.4375*G(NX-2)
         DEALLOCATE(G)
*     ---------------------------------------------------------------- *
*     Checking up the interpolation accuracy (for 'Sp1')
*     ---------------------------------------------------------------- *
         IF (Quiz.AND.Mult.LE.1) THEN
           discrA=Zero
           discrQ=Zero
           Null=0
           DO m=1,NX
             IF (F(m).NE.Zero) THEN
               X=Xmin+StepX*(m-1)
               Z=One-Sp1(Issue,C,X)/F(m)
               X=ABS(Z)
               IF (X.GT.discrM) THEN
                 discrM=X
                 m1=m
            endIF
               discrA=discrA+Z
               discrQ=discrQ+X**2
                                ELSE                                     F(m)=0
               Null=Null+1
          endIF
        endDO
           IF (Tie.AND.Mult.LE.1) THEN
             WRITE(*,1) Issue,Mode,' Fixed',NX
                                  ELSE                                   Unfixed tips
             WRITE(*,1) Issue,Mode,'Untied',NX
        endIF
           IF (Null.LT.NX) THEN
             discrA=discrA/(NX-Null)
             discrQ=SQRT(discrQ)/(NX-Null)
             WRITE(*,2) discrM,m1,Null,discrA,discrQ
                           ELSE                                          Null=NX
             WRITE(*,3)
        endIF
      endIF
         RETURN

    1 FORMAT(/3x,73('_')/             '  | ',71('-'),' |'/
     #'  | Test of Inter1 [Issue=',I4,', Mode=',I1,', ',A6,' Tips]',
     #'  in',I5,' reference points |'/'  | ',71('-'),' |')
    2 FORMAT('  |   Maximum Discrepancy =',1pd9.1,'   [node =',I6,
     #',   Null =',I6,']', 7x,'|'/'  |  Averaged Discrepancy =',1pd9.1,
     #40x,'|'/'  | Quadratic Discrepancy =',1pd9.1,40x,'|'/
     #                                               '  |',73('_'),'|'/)
    3 FORMAT('  | Input function is null equation',41x,'|'/
     #                                               '  |',73('_'),'|'/)
      END SUBROUTINE Coeff1

************************************************************************
      FUNCTION Inter1(Issue,C,X,Mult,Tie0,Eps,Mode0,NX0,Xmin0)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)
         LOGICAL(2) Tie,Tie0
                  REAL Inter1,Invers

      SAVE Tie,Mode,NX,Delta,StepX,Xmin,Xmax,Gmin,Gmax,TipL,TipR
      SAVE m,A1,A2,A3,Nstop

      PARAMETER (Zero=0, One=1, Two=2, Half=One/2, Quart=One/4)
      PARAMETER (Const  = 1.0)

      PARAMETER (MaxIss =10000)                                          User's setting (0 < MaxIss < 100)
      PARAMETER (MaxNst =   10)                                          User's setting (MaxNst > 0)
      PARAMETER (Nlog   =   10)                                          LogFile (must be opened in MAIN)

************************************************************************
*                                                                      *
*     ATTENTION! Maximum value of Issue is equal to 99 and the ma-     *
*                ximum NX0 is 9999.  If however there is a need in     *
*                more than 99 issues(that is almost beyond belief)     *
*                or in NX more  than 9999 (that is quite possible)     *
*                one must to reconstruct slightly the FORMAT under     *
*                the label 109; this is an easy matter.                *
*                                                                      *
*     ATTENTION! Variable  X  plays 2  different roles:  for entry     *
*                'Root' it is the right part of equation F(Root)=X     *
*                and  for the  rest entries it  is the argument of     *
*                function F. Clearly X must fall in the range Xmin     *
*                to Xmax but for reasonably smooth  functions  one     *
*                can take risks  to use  X  for the extended range     *
*                Xmin - StepX/2 to Xmax + StepX/2.  A tool for ex-     *
*                tending the interpolation range is the pa rameter     *
*                Eps.                                                  *
*                                                                      *
************************************************************************

      DIMENSION  Tie(MaxIss),Mode(MaxIss),  NX(MaxIss),C(*)
      DIMENSION Xmin(MaxIss),Gmin(MaxIss),TipL(MaxIss),Delta(MaxIss)
      DIMENSION Xmax(MaxIss),Gmax(MaxIss),TipR(MaxIss),StepX(MaxIss)

      Invers(Func)=ASIN(Func)                                            Example of inverse function

       Tie (Issue)=Tie0                                                  "procedure code" (see 'Coeff1')
       Mode(Issue)=Mode0                                                 "interpolation code" (see 'Coeff1')
       Xmin(Issue)=Xmin0                                                 left end point of interpolation range
       Xmax(Issue)=X                                                     right end point of interpolation range
       NX  (Issue)=NX0                                                   number of the interpolation nodes
                 S=(X-Xmin0)/(NX0-1)
                                                                         =(Xmax(Issue)-Xmin(Issue))/(NX(Issue)-1)
      StepX(Issue)=S                                                     interpolation step
      Delta(Issue)=Eps*S                                                 constant to manage the "clamping"
             Nstop=0                                                     a constant to manage the emergency exit
*     ---------------------------------------------------------------- *
*     Installation of tips and extrema. Check up of monotony.
*     ---------------------------------------------------------------- *
************************************************************************
*                                                                      *
*     When (in entry Root) X is brought into the  Delta-neighbour-     *
*     hood of Gmin or Gmax,  the corresponding boundary  point for     *
*     the  interpolation range is  taken as the root of considered     *
*     equation.                                                        *
*                                                                      *
*     When (in entries Sp1/rSp1) X is brought into the Delta-neig-     *
*     hbourhood of Xmin or Xmax, the value of F is taken to be the     *
*     corresponding tip of F.                                          *
*                                                                      *
************************************************************************
         IF (Tie0) THEN
           TipL(Issue)=C(NX0+1)                                          Left  tip of F/LOG(F)/Transf(F)
           TipR(Issue)=C(NX0+2)                                          Right tip of F/LOG(F)/Transf(F)
           IF (Mult.EQ.1) THEN                                           preparations for entry 'Root'
             IF (C(1).LE.C(NX0)) THEN
               A1=C(1)
               A2=C(NX0)
               A0=A2
               DO m=2,NX0-1
                 A0=MIN(A0,C(m))
                 IF (A0.LE.A1) Nstop=m
            endDO
               IF (Nstop.GT.0) THEN
                 D=Xmin0+S*(Nstop-1)
                 WRITE(Nlog,101)
                 WRITE(Nlog,102) A0,D
                 WRITE(Nlog,104) A1,A2,Xmin0,X
                 STOP 'CRASH LANDING FROM INTER1 at Mult=1'
            endIF
                               ELSE   ! C(1) > C(NX0)
               A1=C(NX0)
               A2=C(1)
               A0=A1
               DO m=2,NX0-1
                 A0=MAX(A0,C(m))
                 IF (A0.GE.A2) Nstop=m
            endDO
               IF (Nstop.GT.0) THEN
                 D=Xmin0+S*(Nstop-1)
                 WRITE(Nlog,101)
                 WRITE(Nlog,103) A0,D
                 WRITE(Nlog,104) A2,A1,Xmin0,X
                 STOP 'CRASH LANDING FROM INTER1 at Mult=1'
            endIF
          endIF
             Gmin(Issue)=A1                                              min value of F/LOG(F)/Transf(F)
             Gmax(Issue)=A2                                              max value of F/LOG(F)/Transf(F)
        endIF
      endIF
         Inter1=Zero                                                     sending discrM to 'Coeff1'
         C(NX0+1)=S                                                      sending StepX  to 'Coeff1'

         RETURN
************************************************************************
*                                                                      *
*     ------------------------------------------------------------     *
*            D E S C R I P T I O N    O F    E N T R I E S             *
*     ------------------------------------------------------------     *
*                                                                      *
*     ENTRY 'Root' to find the roots of following equations            *
*                                                                      *
*                  #         F(Root)  = X (if Mode0=1)                 *
*                  #     LOG(F(Root)) = X (if Mode0=2)                 *
*                  #  Transf(F(Root)) = X (if Mode0=3)                 *
*                                                                      *
*     ATTENTION! Entry 'Root'  can be used ONLY  for DIMINISHED or     *
*                INCREASING  functions (or corresponding transfor-     *
*                mations). An  extension to  the case of ARBITRARY     *
*                functions [when the above equations may have more     *
*                than 1 root] is the subroutine 'Roots'.               *
*     ------------------------------------------------------------     *
*                                                                      *
*     ENTRY 'Sp1'  to find the value of function F in the point X.     *
*                                                                      *
*     ENTRY 'rSp1' to find the value of  some  ANOTHER function  f     *
*            interpolated just as F in the SAME point X.  The call     *
*            of this  subfunction  must be placed  DIRECTLY BEHIND     *
*            the call of 'Sp1'.                                        *
*                                                                      *
*     ------------------------------------------------------------     *
*     ATTENTION! Contrary  to entry  'Sp1', the entry  'rSp1' does     *
*                not tie the spline to the tips.                       *
*     ------------------------------------------------------------     *
*                                                                      *
************************************************************************

*     ==================================================================
      ENTRY Root(Issue,C,X)
*     ==================================================================
         IF (Tie(Issue)) THEN
           D=Delta(Issue)
           A1=X-Gmin(Issue)
           A2=Gmax(Issue)-X
           IF (A1.LT.-D.OR.A2.LT.-D) THEN                                ERROR MESSAGE
             Root=Xmin(Issue)
             WRITE(Nlog,105) X
             GOTO 4
                                     ELSE                                ALL CORRECT
             IF (A1.LE.D) THEN
               Root=Xmin(Issue)
               RETURN
          endIF
             IF (A2.LE.D) THEN
               Root=Xmax(Issue)
               RETURN
          endIF
        endIF
      endIF
         DO 3 m=1,NX(Issue)
           A0=    C(m  )
           A1=Two*C(m+1)
           A2=    C(m+2)
           A3=A0-A1+A2
           IF (A3.NE.Zero)     GOTO 1
           S=A2-A0
           IF (S .EQ.Zero)     GOTO 3
           Z=(X-Quart*A3-A1)/S
                               GOTO 2
    1      S=Half*(A0-A2)
           A1=S**2+A3*(X-Quart*A3-A1)
           IF (A1.LT.Zero)     GOTO 3
           A1=SQRT(A1)
           Z=(S-A1)/A3
           IF (ABS(Z).LE.Half) GOTO 2
           Z=(S+A1)/A3
           IF (ABS(Z).GT.Half) GOTO 3
    2      S=Z+m-One
           IF (NINT(Z).EQ.0) THEN
             Root=StepX(Issue)*S+Xmin(Issue)
             RETURN
        endIF
    3 endDO
         WRITE(Nlog,106) X
    4    WRITE(Nlog,107) Gmin(Issue)-D,Gmax(Issue)+D
         WRITE(Nlog,108) Xmin(Issue)-D,Xmax(Issue)+D
         WRITE(Nlog,109) D,Issue,NX(Issue)
         Nstop=Nstop+1
         IF (Nstop.GT.MaxNst) STOP 'CRASH LANDING FROM ROOT'
         RETURN

*     ==================================================================
      ENTRY Sp1(Issue,C,X)
*     ==================================================================
         D =Delta (Issue)
         A1=X-Xmin(Issue)
         A2=  Xmax(Issue)-X
         IF (A1.LT.-D.OR.A2.LT.-D) THEN
           Sp1=TipL(Issue)
           WRITE(Nlog,105) X
           WRITE(Nlog,108) Xmin(Issue)-D,Xmax(Issue)+D
           WRITE(Nlog,109) D,Issue,NX(Issue)
           Nstop=Nstop+1
           IF (Nstop.GT.MaxNst) STOP 'CRASH LANDING FROM Sp1'
           RETURN
                                   ELSE                                  ALL CORRECT
           IF (Tie(Issue)) THEN
             IF (A1.LE.D) THEN
               Sp1=TipL(Issue)
               RETURN
          endIF
             IF (A2.LE.D) THEN
               Sp1=TipR(Issue)
               RETURN
          endIF
        endIF
      endIF

      A0=A1/StepX(Issue)
       m=NINT(A0)                                                        =INT(A0+Half)
      A0=A0-m
       S=A0**2+Quart
      A1=S-A0
      A2=Two*(One-S)
      A3=S+A0
      GOTO (10,20,30,40) Mode(Issue)
   10 Sp1=       A1*C(m+1)+A2*C(m+2)+A3*C(m+3)
      RETURN
   20 Sp1=EXP   (A1*C(m+1)+A2*C(m+2)+A3*C(m+3))
      RETURN
   30 Sp1=EXP   (A1*C(m+1)+A2*C(m+2)+A3*C(m+3))-Const
      RETURN
   40 Sp1=Invers(A1*C(m+1)+A2*C(m+2)+A3*C(m+3))
      RETURN

*     ==================================================================
      ENTRY rSp1(Issue,C)
*     ==================================================================
      GOTO (100,200,300,400) Mode(Issue)
  100 rSp1=       A1*C(m+1)+A2*C(m+2)+A3*C(m+3)
      RETURN
  200 rSp1=EXP   (A1*C(m+1)+A2*C(m+2)+A3*C(m+3))
      RETURN
  300 rSp1=EXP   (A1*C(m+1)+A2*C(m+2)+A3*C(m+3))-Const
      RETURN
  400 rSp1=Invers(A1*C(m+1)+A2*C(m+2)+A3*C(m+3))
      RETURN

  101 FORMAT(/'  Function is neither diminished not increasing.',
     #        '  Try to use subroutine ROOTS.'/)
  102 FORMAT( '  Real minimum:',1PD22.14,' (x =',1PD22.14,')')
  103 FORMAT( '  Real maximum:',1PD22.14,' (x =',1PD22.14,')')
  104 FORMAT( '          Tips:',1PD22.14,'  and',1PD22.14/
     #        '    Boundaries:',1PD22.14,'  and',1PD22.14/75('-')/)
  105 FORMAT(/'  Mistake! X =',1PD20.13,' lies OUTSIDE the range')
  106 FORMAT(/'  Mistake! Root escaped detection despite the fact that'/
     #                 '  X =',1PD20.13,' lies INSIDE the range')
  107 FORMAT( '  Gmin-Delta =',1PD20.13,' to Gmax+Delta =',1PD20.13)
  108 FORMAT( '  Xmin-Delta =',1PD20.13,' to Xmax+Delta =',1PD20.13)
  109 FORMAT( '  Delta =',1PD8.2,'  [Issue =',I3,', NX =',I5,']'/)
      END FUNCTION Inter1

************************************************************************
      SUBROUTINE IRoots(C,X,Mult,Tie0,Eps,NX0,Xmin0)
************************************************************************
*                                                                      *
*   This subroutine is completely similar to the entry 'Root'  of the  *
*   function 'Inter1' but it can be used for arbitrary functions.      *
*                                                                      *
*   Tie and Eps are used here to extend the interpolation range;       *
*                                                                      *
*   Kmax is the maximum number of roots expected (Kmax < MaxMax), and  *
*   constants from common domain are being used to manage the output.  *
*                                                                      *
*   In contrast with entry 'Root', parameters Gmin and Gmax are being  *
*   used here for somewhat another aims  (in particular,  to speed up  *
*   the calculations) and hence logical Tie also plays another role.   *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)
               LOGICAL(2) Tie,Tie0,outP,outW

         SAVE Kmax,Tie,NX,Xmin,Xmax,Delta,StepX,Gmin,Gmax

         PARAMETER (MaxMax = 99)                                         User's setting (0 < MaxMax < 100)

         PARAMETER (Zero=0, One=1, Two=2, Half=One/2, Quart=One/4)
         DIMENSION C(*)
         ALLOCATABLE R(:)
*     ---------------------------------------------------------------- *
*        ATTENTION ! If it is impossible to use  ALLOCATABLE  dimensi- *
*                    on R(:), comment  lines 509, 552, 602 and  uncom- *
*                    ment lines 515, 516, 517, 518.                    *
*                                                                      *
*        DIMENSION R(MaxMax)
*        INCLUDE 'E:/FORTRAN/Sources/Shared/Include/Heap.fi'             V.A. Naumov
*        INCLUDE 'G:/FORTRAN/Sources/Shared/Include/Heap.fi'             K. Kuzmin
*        EQUIVALENCE (R(1),Heap(1))                                      memory saving
*     ---------------------------------------------------------------- *
         COMMON /R_Info/ Nlog,Nout,outP,outW

         Kmax=Mult
         IF (MaxMax.LT.Kmax) THEN
           WRITE(Nlog,101) MaxMax,Kmax
           STOP 'Fault! Increase MaxMax.'
      endIF
         Tie=Tie0
         NX=NX0
         Xmin=Xmin0
         Xmax=X
         Delta=Eps                                                       /NX
         StepX=(Xmax-Xmin)/(NX-1)
         Gmin=C(1)
         Gmax=Gmin
         DO m=2,NX
           Gmin=MIN(Gmin,C(m))
           Gmax=MAX(Gmax,C(m))
      endDO
         C(1)=StepX
         RETURN

*     ==================================================================
      ENTRY Roots(C,X)
*     ==================================================================
         IF (Tie) THEN
           IF (ABS(X-Gmin).LE.Delta) Right=Gmin+Delta
           IF (ABS(X-Gmax).LE.Delta) Right=Gmax-Delta
                  ELSE
           Right=X
      endIF

         ALLOCATE (R(MaxMax))

         N=0
         DO 4 m=1,NX
           A0=  C(m  )
           A1=Two*C(m+1)
           A2=  C(m+2)
           A3=A0-A1+A2
           IF (A3.NE.Zero)                 GOTO 1
           B=A2-A0
           IF ( B.EQ.Zero)                 GOTO 4
           D=(Right-Quart*A3-A1)/B
                                           GOTO 2
    1      B=Half*(A0-A2)
           A1=B**2+A3*(Right-Quart*A3-A1)
           IF (A1.LT.Zero)                 GOTO 4
           A1=SQRT(A1)
           D=(B-A1)/A3
           IF (ABS(D).LE.Half)             GOTO 2
           D=(B+A1)/A3
           IF (ABS(D).GT.Half)             GOTO 4
    2      B=D+m-1
           IF (NINT(D).NE.0)               GOTO 4
           Root0=StepX*B+Xmin
           IF (N.LT.1)                     GOTO 3
           DO i=1,N
             IF (ABS(R(i)-Root0).LE.Delta) GOTO 4
        endDO
    3      IF (N.LT.Kmax) THEN
             N=N+1
             R(N)=Root0
                          ELSE                                           N = Kmax
                                           GOTO 5
        endIF
    4 endDO
    5    IF (N.EQ.Kmax) THEN                                             Number of roots exceeds Kmax
           IF (outP) WRITE(Nlog,102) Kmax
           IF (outW) WRITE(Nout,102) Kmax
      endIF
         IF (N.EQ.0) THEN                                                Roots escaped detection
           IF (outP) WRITE(Nlog,103) NX,Delta,X,Xmin,Xmax
           IF (outW) WRITE(Nout,103) NX,Delta,X,Xmin,Xmax
                     ELSE                                                Roots are founded
           IF (outP) WRITE(Nlog,104) X,N
           IF (outW) WRITE(Nout,104) X,N
           DO i=1,N
             IF (outP) WRITE(Nlog,105) i,R(i)
             IF (outW) WRITE(Nout,105) i,R(i)
        endDO
      endIF
         DEALLOCATE(R)

         RETURN

  101 FORMAT(/'  Mistake: MAX(Kmax) =',I6,', Kmax =',I6//)
  102 FORMAT(/'  Number of roots exceeds Kmax =',I6/)
  103 FORMAT(/'  Roots escaped detection.  NX =',I6,', Delta =',1PD10.4/
     #       /'  X =',1PD21.14,' Xmin=',1PD21.15,' Xmax=',1PD21.14/)
  104 FORMAT(/'  X =',1PD21.14,', K =',I6/)
  105 FORMAT(I8,'  Root =',1PD21.14)
      END SUBROUTINE IRoots