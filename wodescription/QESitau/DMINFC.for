************************************************************************
      SUBROUTINE DMINFC(F,A,B,EPS,DELTA,X,Y,LLM)
************************************************************************
*                                                                      *
*     This SUBROUTINE returns the minimum of a function of one va-     *
*     riable. Mathlib, D503                                            *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z), INTEGER (I,J,K,L,M,N)

         LOGICAL(2) LLM,LLT,LGE

         REAL,PARAMETER::
     #        W5= 2.236067977499790d+00,
     #        HV= (3-W5)/2,
     #        HW= (W5-1)/2,
     #        R1= 1.0d+00,
     #        HF= R1/2

         N=-1
         IF (A.ne.B) N=nint(2.08*log(abs((A-B)/EPS)))
         C=A
         D=B
         IF (A.gt.B) THEN
           C=B
           D=A
      endIF
         LLT=.TRUE.
         LGE=.TRUE.
    1    H  =D-C
         IF (N.lt.0) THEN
           X  =HF*(C+D)
           Y  =F(X)
           LLM=abs(X-A).gt.DELTA .and. abs(X-B).gt.DELTA
           RETURN
      endIF
         IF (LLT) THEN
           V =C+HV*H
           FV=F(V)
      endIF
         IF (LGE) THEN
           W =C+HW*H
           FW=F(W)
      endIF
         IF (FV.lt.FW) THEN
           LLT=.TRUE.
           LGE=.FALSE.
           D  =W
           W  =V
           FW =FV
                       ELSE
           LLT=.FALSE.
           LGE=.TRUE.
           C  =V
           V  =W
           FV =FW
      endIF
         N=N-1
         GOTO 1
      
         RETURN
      END SUBROUTINE DMINFC