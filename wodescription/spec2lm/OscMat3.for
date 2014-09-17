************************************************************************
      SUBROUTINE OscMat3
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012        *
*                                BLTP JINR, Dubna, Russia, 2013/02/25  *
************************************************************************
*                                                                      *
*     This SUBROUTINE returns the neutrino mixing matrix V and         * 
*     PREM coefficients                                                *
*                                                                      *
************************************************************************

         USE OscMatParameters

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
         SAVE
                 REAL
     #                th12,th23,th13,
     #                c12,c13,c23,s12,s13,s23,CPV,
     #                PREM(4,0:KN)
              COMPLEX
     #                V(3,3)
              INTEGER
     #                j,Nr,n_NT
     
               COMMON /PREM/    PREM
               COMMON /Vmatrix/ V
               COMMON /Numbers/ Nr
               COMMON /n_NT/    n_NT

         th12=asin(sqrt(s2th12))
         th23=asin(sqrt(s2th23))
         th13=asin(sqrt(s2th13))

         c12=cos(th12)
         c13=cos(th13)
         c23=cos(th23)

         s12=sin(th12)
         s13=sin(th13)
         s23=sin(th23)

         CPV=n_NT*CPVdel

         V(1,1)= c12*c13
         V(1,2)= s12*c13
         V(1,3)= s13*exp(-i*CPV)

         V(2,1)=-s12*c23-c12*s23*s13*exp(i*CPV)
         V(2,2)= c12*c23-s12*s23*s13*exp(i*CPV)
         V(2,3)= s23*c13

         V(3,1)= s12*s23-c12*c23*s13*exp(i*CPV)
         V(3,2)=-c12*s23-s12*c23*s13*exp(i*CPV)
         V(3,3)= c23*c13
         
         RETURN
         
      END SUBROUTINE OscMat3
