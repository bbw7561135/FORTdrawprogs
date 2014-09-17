************************************************************************
      FUNCTION oscvac3(E,L)
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/08/21  *
************************************************************************
*                                                                      *
*     This FUNCTION returns the probability of vacuum neutrino         *
*     oscillation from one lepton flavor to another for any neutrino   *
*     energy E (GeV) and distance L (km)                               *
*                                                                      *
************************************************************************

         USE OscMatNeeds, ONLY: zero,one,i,hbarc
         USE OscParameters

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
         SAVE

                 REAL
     #                L0(3,3),
     #                CPV
              COMPLEX
     #                v(3,3),me(3,3,3,3),
     #                g(3,3)/one,3*zero,one,3*zero,one/,
     #                gp(3,3),
     #                ex(3,3)
              INTEGER
     #                alpha,beta,a,b,j,k,n_NT

               COMMON /n_NT/n_NT

         CALL OscPar

         CPV=n_NT*CPVdel

         g(3,3)=exp(i*CPV)                                               Dirac's matrix
         gp=transpose(conjg(g))

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
      ENTRY Pab(alpha,beta,E,L)
*     ==================================================================
         L0=dm/(2*E)/hbarc                                               1/L_ij!/(2*pi)
         ex=i*L*L0*1.0d-13                                               !*2*pi
         
         s=0
         DO j=1,3
           DO k=1,3
             s=s+me(alpha,beta,j,k)*exp(ex(j,k))
        endDO
      endDO

         Pab=s
         RETURN

      END FUNCTION oscvac3
