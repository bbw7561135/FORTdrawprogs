************************************************************************
      SUBROUTINE OscMat3
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012        *
*                                BLTP JINR, Dubna, Russia, 2013/06/05  *
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
     #                PREM(4,0:KN),
     #                rad(0:KN)
              COMPLEX
     #                V(3,3)
              INTEGER
     #                n_NT,j,LL,UL,NDet,Nr,beta(0:KN)
     
               COMMON /PREM/    PREM,beta
               COMMON /rad/     rad
               COMMON /Vmatrix/ V
               COMMON /Numbers/ NDet,Nr
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

         LL=0
         UL=N1
         DO j=LL,UL
           PREM(1,j)= 1.30885d+01                                       
           PREM(2,j)= zero                                              
           PREM(3,j)=-8.8381d+00                                        
           PREM(4,j)= zero
           beta(j)=0
      endDO
         LL=UL+1
         UL=UL+N2
         DO j=LL,UL
           PREM(1,j)= 1.25815d+01
           PREM(2,j)=-1.26380d+00
           PREM(3,j)=-3.64260d+00
           PREM(4,j)=-5.52810d+00
           beta(j)=0
      endDO
         LL=UL+1
         UL=UL+N3
         DO j=LL,UL
           PREM(1,j)= 7.95650d+00
           PREM(2,j)=-6.47610d+00
           PREM(3,j)= 5.52830d+00
           PREM(4,j)=-3.08070d+00
           beta(j)=0
      endDO
         LL=UL+1
         UL=UL+N4
         DO j=LL,UL
           PREM(1,j)= zero
           PREM(2,j)= zero
           PREM(3,j)= zero
           PREM(4,j)= zero
           beta(j)=1
      endDO
         LL=UL+1
         UL=UL+N5
         DO j=LL,UL
           PREM(1,j)= 2.90d+00
           PREM(2,j)= zero
           PREM(3,j)= zero
           PREM(4,j)= zero
           beta(j)=0
      endDO
         LL=UL+1
         UL=UL+N6
         DO j=LL,UL
           PREM(1,j)= 2.60d+00
           PREM(2,j)= zero
           PREM(3,j)= zero
           PREM(4,j)= zero
           beta(j)=0
      endDO
         LL=UL+1
         UL=UL+N7
         DO j=LL,UL
           PREM(1,j)= 1.02d+00
           PREM(2,j)= zero
           PREM(3,j)= zero
           PREM(4,j)= zero
           beta(j)=0
      endDO

         CALL SubRadii(rad)
         
         CALL Dichotomy(RadDet,NDet)
         rad(NDet)=RadDet
         
         RETURN
         
      END SUBROUTINE OscMat3
