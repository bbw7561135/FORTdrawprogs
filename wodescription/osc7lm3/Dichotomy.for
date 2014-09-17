************************************************************************
      SUBROUTINE Dichotomy(Rd,N)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012        *
*                                BLTP JINR, Dubna, Russia, 2013/02/25  *
************************************************************************
*                                                                      *
*     This SUBROUTINE finds number N of the Earth sublayer which       * 
*     contained preassigned radius R                                   *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                rad(0:KN),Rd
              INTEGER
     #                j,Xmx,Xmn,X,N,ResX
     
               COMMON /rad/rad
               
         Xmx=KN
         Xmn=0
     
         IF (Rd.EQ.zero) THEN
           N=0
           RETURN
      endIF
         
         DO j=Xmn,Xmx
           X=Xmx-Xmn
           IF (X.EQ.1) THEN
             N=Xmx
                       ELSE  
             ResX=mod(X,2)
             IF (ResX.EQ.0) THEN                                         X is even
               X=X/2
               IF (Rd.EQ.rad(Xmn+X)) THEN
                 N=Xmn+X
                                     ELSE
                 IF (Rd.LT.rad(Xmn+X)) THEN
                   Xmx=Xmn+X
                                       ELSE 
                   Xmn=Xmn+X
              endIF
            endIF
                            ELSE
               X=(X+1)/2                                                 X is odd
               IF (Rd.EQ.rad(Xmn+X)) THEN
                 N=Xmn+X
                                     ELSE
                 IF (Rd.LT.rad(Xmn+X)) THEN
                   Xmx=Xmn+X
                                       ELSE
                   Xmn=Xmn+X
              endIF
            endIF
          endIF
        endIF
      endDO
     
      RETURN
         
      END SUBROUTINE Dichotomy
