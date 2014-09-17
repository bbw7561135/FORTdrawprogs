************************************************************************
      FUNCTION DZEROX(A0,B0,EPS,MAXF,F,MODE)
************************************************************************
*                                                                      *
*     This FUNCTION returns the zero of a function of one real va-     *
*     riable. Mathlib, C200.  The program is based on algorithm by     *
*     J.C.P. Bus  and T.J. Dekker, "Two efficient algorithms  with     *
*     guaranteed  convergence  for finding  a zero of a function,"     *
*     ACM  Trans. Math. Soft. 1 (1975) 330-345.                        *
*                                                                      *
*     MODE 1 for algorithm "M", MODE 2 for algorithm "R".              *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z), INTEGER (I,J,K,L,M,N)

         LOGICAL(2)  LMT(2)

         REAL,PARAMETER::
     #        Z1  = 1.0d+00,
     #        HALF= Z1/2

         REAL IM1(2), IM2(2)

         DATA IM1 /2,3/, IM2 /-1,3/

         IF (MODE.ne.1 .and. MODE.ne.2) THEN
           C = 0.0d+00
           WRITE(6,101) MODE
           STOP 'STOP 1. FUNCTION DZEROX'
      endIF
         FA= F(B0)
         FB= F(A0)
         IF (FA*FB.gt.0.0d+00) THEN
           C= 0.0d+00
           WRITE(6,102) A0, B0
           STOP 'STOP 2. FUNCTION DZEROX'
      endIF
         ATL   = abs(EPS)
         B     = A0
         A     = B0
         LMT(2)=.true.
         MF    = 2
    1    C     = A
         FC    = FA
    2    IE    = 0
    3    IF (abs(FC).lt.abs(FB)) THEN
           IF (C.ne.A) THEN
             D = A
             FD= FA
        endIF
           A = B
           B = C
           C = A
           FA= FB
           FB= FC
           FC= FA
      endIF
         TOL= ATL*(1+abs(C))
         H  = HALF*(C+B)
         HB = H-B
         IF (abs(HB).gt.TOL) THEN
           IF (IE.gt.IM1(MODE)) THEN
             W     = HB
                                ELSE
             TOL   = TOL*sign(Z1,HB)
             P     = (B-A)*FB
             LMT(1)= IE .le. 1
             IF (LMT(MODE)) THEN
               Q     = FA-FB
               LMT(2)=.FALSE.
                            ELSE
               FDB= (FD-FB)/(D-B)
               FDA= (FD-FA)/(D-A)
               P  = FDA*P
               Q  = FDB*FA-FDA*FB
          endIF
             IF (P.lt.0.0d+00) THEN
               P=-P
               Q=-Q
          endIF
             IF (IE.eq.IM2(MODE)) P = P + P
             IF (P .eq.0.0d+00 .or. P.le.Q*TOL) THEN
               W= TOL
         ELSEIF (P.lt.HB*Q                    ) THEN
               W= P/Q
                                                ELSE
               W= HB
          endIF
        endIF
           D = A
           A = B
           FD= FA
           FA= FB
           B = B+W
           MF= MF+1
           IF (MF.gt.MAXF) THEN 
             WRITE(6,103) 
             STOP 'STOP 3. FUNCTION DZEROX'
        endIF
           FB= F(B)
           IF (FB.eq.0 .or. sign(Z1,FC).eq.sign(Z1,FB)) GOTO 1
           IF (W.eq.HB) GOTO 2
           IE= IE+1
           GOTO 3
      endIF
         DZEROX= C
   99    CONTINUE

         RETURN
  100 FORMAT('ERROR in Mathlib C200 subroutine:')
  101 FORMAT('MODE = ',I3,' ILLEGAL')
  102 FORMAT('F(A) AND F(B) HAVE THE SAME SIGN, A = ',1P,D15.8,
     #       ', B = ',D15.8)
 103  FORMAT('TOO MANY FUNCTION CALLS')
      END FUNCTION DZEROX