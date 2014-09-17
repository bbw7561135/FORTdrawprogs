************************************************************************
      FUNCTION oscvac3(n_alpha,n_beta,E,L)
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/08/19  *
************************************************************************
*                                                                      *
*     This FUNCTION returns the probability of vacuum neutrino         *
*     oscillation from one lepton flavor to another for any neutrino   *
*     energy E (GeV) and distance L (km)                               *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
         SAVE

                 REAL,PARAMETER::
      INCLUDE 'OscParFogli.f'
!      INCLUDE 'OscParValle.f'
     #                dm32  =-(dm21+dm13),
     #                dm12  =-dm21,
     #                dm31  =-dm13,
     #                dm23  =-dm32
                 REAL
     #                dm(3,3)/0,dm21,dm31,dm12,0,dm32,dm13,dm23,0/,      !WTF?!
     #                L0(3,3),
     #                o12(3,3)/9*0/,
     #                o23(3,3)/9*0/,
     #                o31(3,3)/9*0/,
     #                CPV
              COMPLEX
     #                v(3,3),me(3,3,3,3),
     #                g(3,3)/1,3*0,1,3*0,1/,
     #                gp(3,3),
     #                ex(3,3)
              INTEGER
     #                a,b,j,k,n_NT

               COMMON /n_NT/n_NT

         CPV=n_NT*CPVdel

         g(3,3)=exp((0,1)*CPV)                                           Dirac's matrix
         gp=transpose(conjg(g))
         
         th12=asin(sqrt(s2th12))
         th23=asin(sqrt(s2th23))
         th13=asin(sqrt(s2th13))
         
         o12(1,1)= cos(th12); o12(1,2)= sin(th12)
         o12(2,1)=-sin(th12); o12(2,2)= cos(th12)                        Solar matrix
         o12(3,3)=1

         o23(1,1)=1
         o23(2,2)= cos(th23); o23(2,3)= sin(th23)                        Atmospheric matrix
         o23(3,2)=-sin(th23); o23(3,3)= cos(th23)

         o31(1,1)= cos(th13); o31(1,3)= sin(th13)
         o31(2,2)=1                                                      Reactor matrix
         o31(3,1)=-sin(th13); o31(3,3)= cos(th13)

         v=matmul(o23,matmul(g,matmul(o31,matmul(gp,o12))))              PMNS-matrix
         
         DO a=1,3
           DO b=1,3
             DO j=1,3
               DO k=1,3
               me(a,b,j,k)=v(a,j)*v(b,k)*conjg(v(a,k)*v(b,j))
            endDO
          endDO
        endDO
      endDO
         
         oscvac3=one
         RETURN
*     ==================================================================
      ENTRY Pab(n_alpha,n_beta,E,L)
*     ==================================================================
         a=n_alpha
         b=n_beta
         L0=dm/(2*E)/hbarc                                               1/L_ij!/(2*pi)
         ex=(0,1)*L*L0*1.0d-13                                          !*2*pi
         
         s=0
         DO j=1,3
           DO k=1,3
             s=s+me(a,b,j,k)*exp(ex(j,k))
        endDO
      endDO

         Pab=s
         RETURN

      END FUNCTION oscvac3
