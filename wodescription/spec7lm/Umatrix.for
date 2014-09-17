************************************************************************
      SUBROUTINE Umatrix(E,Psi,U)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2013/02/25  *
************************************************************************
*                                                                      *
*     This SUBROUTINE calculates matrix U                              *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                rho(0:KN),
     #                W(3),q(3),
     #                Omega(3),SqOmega,Omegas,
     #                ModH(3),SqH(3),H2,
     #                Fi(3),Fis,     
     #                KubV,Angle,
     #                E(3,0:KN),
     #                Psi(3,3,0:KN),
     #                Xi(3,3),Xis(3),
     #                U(3,3,0:KN)
              COMPLEX
     #                H(3),Bk(3,3)
              INTEGER
     #                a,b,j,k,m,n,
     #                l,NDet,Nr,n_NT
     
               COMMON /UmtrxDblPr/    W,Fi,Fis,ModH
               COMMON /UmtrxDblCmplx/ H
               COMMON /Numbers/       NDet,Nr
               COMMON /n_NT/          n_NT

         CALL EarthRho(Rho)

!         n_NT= 1

         DO l=KN,Nr,-1
           q(1)=-n_NT*Coeff*Rho(l)*VarDens
           q(2)=-half*q(1)
           q(3)= q(2)

           Omega=W-q

           SqOmega=Omega(1)**2+Omega(2)**2+Omega(3)**2
           SqH=H*conjg(H)
           H2=sum(SqH)
           Omegas=sqrt(sixth*SqOmega+third*H2)
           
           KubV=half*(product(Omega)+two*cos(Fis)*product(ModH)-
     #          sum(Omega*SqH))

           Angle=third*acos(KubV/Omegas**3)
           
           E(1,l)=two*Omegas*cos(Angle+third*twopi)
           E(2,l)=two*Omegas*cos(Angle-third*twopi)
           E(3,l)=two*Omegas*cos(Angle)
           
           DO n=1,3
             Bk(1,n)=H(1)*(E(n,l)-Omega(1))+conjg(H(2)*H(3))
             Bk(2,n)=H(2)*(E(n,l)-Omega(2))+conjg(H(3)*H(1))
             Bk(3,n)=H(3)*(E(n,l)-Omega(3))+conjg(H(1)*H(2))
        endDO

           DO a=1,3
             DO n=1,3
               Psi(a,n,l)=ArgZ(Bk(a,n))
          endDO
        endDO
           
           DO n=1,3
             Xi(1,n)=E(n,l)**2+E(n,l)*Omega(1)+Omega(2)*Omega(3)-SqH(1)
             Xi(2,n)=E(n,l)**2+E(n,l)*Omega(2)+Omega(3)*Omega(1)-SqH(2)  !sum(Omega)=sum(W)-sum(q)=0-0=0
             Xi(3,n)=E(n,l)**2+E(n,l)*Omega(3)+Omega(1)*Omega(2)-SqH(3)    
        endDO
        
           Xis=sum(Xi,1)
           
           DO a=1,3
             DO n=1,3
               U(a,n,l)=sqrt(max((Xi(a,n)/Xis(n)),zero))
          endDO
        endDO
        
      endDO
     
      RETURN
         
      END SUBROUTINE Umatrix
