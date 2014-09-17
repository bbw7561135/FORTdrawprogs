************************************************************************
      MODULE Routines
************************************************************************
*                                                                      *
*     Shared SUBROUTINEs and FUNCTIONs                                 *
*                                                                      *
*                                      ITEP, Moscow, Russia 2005/05/28 *
************************************************************************

         USE PhysMathConstants

         CONTAINS

************************************************************************
      FUNCTION Momentum(E,m)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         x=(m/E)**2
         IF (x.ge.one) THEN; Momentum=Precision
                       ELSE; Momentum=E*sqrt(one-x)
      endIF

         RETURN
      END FUNCTION Momentum

************************************************************************
      FUNCTION sine(cosine)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         IF (cosine.ge.+one) THEN; sine=zero; cosine= one; RETURN; endIF
         IF (cosine.le.-one) THEN; sine=zero; cosine=-one; RETURN; endIF
         sine=sqrt(one-cosine**2)

         RETURN
      END FUNCTION sine

************************************************************************
      FUNCTION Tf(x,Q2)
************************************************************************
*                                                                      *
*     ---------------------------------------------------------------- * ORIGINAL CODE
*        IF (Q2.lt.1.5d-01) THEN
*          Tf=(one+12*0.15/((1+0.15)*(1+64*x**2)))/log(0.15/0.04)
*                           ELSE
*          Tf=(one+12*Q2  /((1+Q2  )*(1+64*x**2)))/log(Q2  /0.04)
*     endIF
*     ---------------------------------------------------------------- *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         IF (Q2.lt.1.5d-01) THEN
           Tf=(one+12*0.15/((1+0.15)*(1+64*x**2)))/log(0.15/0.04)
                            ELSE
           Tf=(one+12*Q2  /((1+Q2  )*(1+64*x**2)))/log(Q2  /0.04)
      endIF

         RETURN
      END FUNCTION Tf

************************************************************************
      FUNCTION Rg(x,Q2,r1,r2,r3)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         Rg=r1*Tf(x,Q2)+r2/Q2+r3/(Q2**2+9.0d-02)

         RETURN
      END FUNCTION Rg

************************************************************************
      FUNCTION Rf(x,Q2,Q2_0,r1,r2,r3)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         Rf=2*Q2_0*Q2/(Q2**2+Q2_0**2)*Rg(x,Q2_0,r1,r2,r3)

         RETURN
      END FUNCTION Rf
************************************************************************

      END MODULE Routines