************************************************************************
      SUBROUTINE Pmatrix(EneV,Theta)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2013/06/05  *
************************************************************************
*                                                                      *
*     This SUBROUTINE calculates neutrino oscillation probabilities P  * 
*     for neutrino energy EneV (eV) and zenith angle Theta (pi/2,pi)   *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                delta(3),
     #                W(3),Fi(3),Fis,ModH(3),
     #                R,Theta,EneV,
     #                MaxD(3,3),
     #                Length(0:KN),
     #                E(3,0:KN),Psi(3,3,0:KN),U(3,3,0:KN),
     #                DynPh(3),
     #                DeltaI(3,3),
     #                P(3,3)
              COMPLEX
     #                V(3,3),H(3),
     #                expDP(3),Y(3,3),
     #                Matrix(3,3),Smatrix(3,3)
              INTEGER
     #                a,b,j,k,m,n,
     #                Nr,l,k1,k2
     
               COMMON /UmtrxDblCmplx/ H
               COMMON /UmtrxDblPr/    W,Fi,Fis,ModH
               COMMON /radii/         R
               COMMON /Numbers/       Nr
               COMMON /MaxD/          MaxD
               COMMON /Vmatrix/       V
               COMMON /P/             P
         
         delta(1)=sixth*(dm13-dm21)/EneV
         delta(2)=sixth*(dm21-dm32)/EneV
         delta(3)=sixth*(dm32-dm13)/EneV
         
         DO a=1,3
           W(a)=zero
           DO j=1,3
           W(a)=W(a)+V(a,j)*conjg(V(a,j))*delta(j)
        endDO
      endDO
      
         DO a=1,3
           H(a)=zero
      endDO
         DO j=1,3
           H(1)=H(1)+V(2,j)*conjg(V(3,j))*delta(j)
           H(2)=H(2)+V(3,j)*conjg(V(1,j))*delta(j)
           H(3)=H(3)+V(1,j)*conjg(V(2,j))*delta(j)
      endDO
         
         DO a=1,3
           Fi(a)=ArgZ(H(a))
      endDO

         Fis=sum(Fi)

         ModH=sqrt(H*conjg(H))
         
         R=RadDet*sin(Theta)
         
          IF (R.LE.Precision) THEN
            Nr=0
      ELSEIF (((Rcr-R).GE.zero).and.(R.GT.zero)) THEN
            Nr=1
      ELSEIF (((R_Earth-R).GE.zero).and.((R-Rcr).GT.zero)) THEN
            Nr=2
       endIF
   
         CALL SegLengths(Length)  
         CALL Umatrix(E,Psi,U)

         Smatrix(1,1)=( one,zero)
         Smatrix(1,2)=(zero,zero)
         Smatrix(1,3)=(zero,zero)

         Smatrix(2,1)=(zero,zero)
         Smatrix(2,2)=( one,zero)
         Smatrix(2,3)=(zero,zero)

         Smatrix(3,1)=(zero,zero)
         Smatrix(3,2)=(zero,zero)
         Smatrix(3,3)=( one,zero)
      
         IF (UTestSw) THEN 
           DO k1=1,3
             DO k2=1,3
               MaxD(k1,k2)=zero
          endDO
        endDO
      endIF

         DO l=KN,Nr,-1
           DO n=1,3
             DynPh(n)=E(n,l)*Length(l)/hbarc
        endDO
           expDP=exp(-i*DynPh)
           
           DO a=1,3
             DO n=1,3
               Y(a,n)=U(a,n,l)*exp(-i*Psi(a,n,l))
          endDO
        endDO
           
           DO a=1,3
             DO b=1,3
               Matrix(a,b)=(zero,zero)
          endDO
        endDO
        
           DO n=1,3
             DO a=1,3
               Matrix(a,a)=Matrix(a,a)+U(a,n,l)*U(a,n,l)*expDP(n)
          endDO

             Matrix(1,2)=Matrix(1,2)+Y(1,n)*Y(2,n)*expDP(n)
             Matrix(2,3)=Matrix(2,3)+Y(2,n)*Y(3,n)*expDP(n)
             Matrix(3,1)=Matrix(3,1)+Y(3,n)*Y(1,n)*expDP(n)
     
             Matrix(2,1)=Matrix(2,1)+conjg(Y(2,n)*Y(1,n))*expDP(n)
             Matrix(3,2)=Matrix(3,2)+conjg(Y(3,n)*Y(2,n))*expDP(n)
             Matrix(1,3)=Matrix(1,3)+conjg(Y(1,n)*Y(3,n))*expDP(n)
        endDO
     
           Smatrix=matmul(Smatrix,Matrix)
                   
           IF (UTestSw) THEN
             CALL UTest(Smatrix,DeltaI)

             DO k1=1,3
               DO k2=1,3
                 IF (abs(DeltaI(k1,k2)).GE.MaxD(k1,k2)) THEN 
                   MaxD(k1,k2)=abs(DeltaI(k1,k2))
              endIF
            endDO
          endDO
        endIF
      endDO
      
         DO l=Nr,NDet
           DO n=1,3
             DynPh(n)=E(n,l)*Length(l)/hbarc
        endDO
           expDP=exp(-i*DynPh)
           
           DO a=1,3
             DO n=1,3
               Y(a,n)=U(a,n,l)*exp(-i*Psi(a,n,l))
          endDO
        endDO
           
           DO a=1,3
             DO b=1,3
               Matrix(a,b)=(zero,zero)
          endDO
        endDO
        
           DO n=1,3
             DO a=1,3
               Matrix(a,a)=Matrix(a,a)+U(a,n,l)*U(a,n,l)*expDP(n)
          endDO

             Matrix(1,2)=Matrix(1,2)+Y(1,n)*Y(2,n)*expDP(n)
             Matrix(2,3)=Matrix(2,3)+Y(2,n)*Y(3,n)*expDP(n)
             Matrix(3,1)=Matrix(3,1)+Y(3,n)*Y(1,n)*expDP(n)
     
             Matrix(2,1)=Matrix(2,1)+conjg(Y(2,n)*Y(1,n))*expDP(n)
             Matrix(3,2)=Matrix(3,2)+conjg(Y(3,n)*Y(2,n))*expDP(n)
             Matrix(1,3)=Matrix(1,3)+conjg(Y(1,n)*Y(3,n))*expDP(n)
        endDO
     
           Smatrix=matmul(Smatrix,Matrix)

           IF (UTestSw) THEN
             CALL UTest(Smatrix,DeltaI)

             DO k1=1,3
               DO k2=1,3
                 IF (abs(DeltaI(k1,k2)).GE.MaxD(k1,k2)) THEN 
                   MaxD(k1,k2)=abs(DeltaI(k1,k2))
              endIF
            endDO
          endDO
        endIF
      endDO
         
         P=Smatrix*conjg(Smatrix)

      RETURN
         
      END SUBROUTINE Pmatrix
