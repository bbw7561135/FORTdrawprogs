************************************************************************
      FUNCTION FunMul(var)
************************************************************************
*                                                                      *
*                                                                      *
*                                    ITEP, Moscow, Russia, 2010/02/10  *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         COMMON        /m_ini/m_ini,mm_ini                               Mass of target nucleon
         COMMON         /E_nu/E_nu                                       Energy of initial neutrino
         COMMON       /W2_lim/dW2,W2_min                                 W^2 kinematic limits

         DIMENSION var(*)

         FunMul= one
         RETURN

*     ================================================================ *
      ENTRY MuL_Q2d2sDISCC_dxdy(var)
*     ================================================================ *
         W2= dW2*var(1)+W2_min
         CALL Q2DIS_lim(E_nu,W2,Q2_min,Q2_max)
         IF (Q2_min.ge.Q2_max) THEN
           MuL_Q2d2sDISCC_dxdy=zero
                               ELSE
           Q2= (Q2_max-Q2_min)*var(2)+Q2_min                             Square of mometum transfer (Q^2=-q^2)
           x = Q2/(Q2+W2-mm_ini)
           CALL xDIS_lim(E_nu,x_min,x_max)
           IF (x.le.x_min .or. x.ge.x_max) THEN
             MuL_Q2d2sDISCC_dxdy= zero
                                           ELSE
             y = Q2/(2*m_ini*x*E_nu)
             CALL yDIS_lim(E_nu,x,y_min,y_max)
             IF (y.le.y_min .or. y.ge.y_max) THEN
               MuL_Q2d2sDISCC_dxdy= zero
                                             ELSE
               MuL_Q2d2sDISCC_dxdy= (Q2_max-Q2_min)*
     #         max(Q2*d2sDISCC_dxdy(E_nu,x,y),Precision)/y
          endIF
        endIF
      endIF
         RETURN

*     ================================================================ *
      ENTRY MuL_d2sDISCC_dxdy(var)
*     ================================================================ *
         W2= dW2*var(1)+W2_min
         CALL Q2DIS_lim(E_nu,W2,Q2_min,Q2_max)
         IF (Q2_min.ge.Q2_max) THEN
           MuL_d2sDISCC_dxdy= zero
                               ELSE
           Q2= (Q2_max-Q2_min)*var(2)+Q2_min                             Square of mometum transfer (Q^2=-q^2)
           x = Q2/(Q2+W2-mm_ini)
           CALL xDIS_lim(E_nu,x_min,x_max)
           IF (x.le.x_min .or. x.ge.x_max) THEN
             MuL_d2sDISCC_dxdy= zero
                                           ELSE
             y= Q2/(2*m_ini*x*E_nu)
             CALL yDIS_lim(E_nu,x,y_min,y_max)
             IF (y.le.y_min .or. y.ge.y_max) THEN
               MuL_d2sDISCC_dxdy= zero
                                             ELSE
               MuL_d2sDISCC_dxdy= (Q2_max-Q2_min)*
     #         max(d2sDISCC_dxdy(E_nu,x,y),Precision)/y
          endIF
        endIF
      endIF
         RETURN

*     ================================================================ *
      ENTRY MuL_d2sDISNC_dxdy(var)
*     ================================================================ *
         W2= dW2*var(1)+W2_min
         CALL Q2DIS_lim(E_nu,W2,Q2_min,Q2_max)
         IF (Q2_min .ge. Q2_max) THEN
           MuL_d2sDISNC_dxdy= Precision
                                 ELSE
           Q2= (Q2_max-Q2_min)*var(2)+Q2_min                             Square of mometum transfer (Q^2=-q^2)
           x = Q2/(Q2+W2-mm_ini)
           CALL xDIS_lim(E_nu,x_min,x_max)
           IF (x.le.x_min .or. x.ge.x_max) THEN
             MuL_d2sDISNC_dxdy= Precision
                                           ELSE
             y= Q2/(2*m_ini*x*E_nu)
             CALL yDIS_lim(E_nu,x,y_min,y_max)
             IF (y.le.y_min .or. y.ge.y_max) THEN
               MuL_d2sDISNC_dxdy= Precision
                                             ELSE
               MuL_d2sDISNC_dxdy= (Q2_max-Q2_min)*
     #         max(d2sDISNC_dxdy(E_nu,x,y),Precision)/y
          endIF
        endIF
      endIF
         RETURN

*     ================================================================ *

      END FUNCTION FunMuL