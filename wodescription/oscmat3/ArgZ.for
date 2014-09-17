************************************************************************
      FUNCTION ArgZ(Z)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012        *
*                                BLTP JINR, Dubna, Russia, 2013/02/24  *
************************************************************************
*                                                                      *
*     This FUNCTION returns the argument of the complex number Z       *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                ArgZ,ReZ,ImZ,t
             COMPLEX
     #                Z
     
         ReZ=real(Z)
         ImZ=aimag(Z)
         t  =atan(abs(ImZ/ReZ))
         
           IF (ReZ<zero) THEN
             IF (ImZ<zero)    THEN
               ArgZ= t+pi
         ELSEIF (ImZ.EQ.zero) THEN
               ArgZ= pi
                              ELSE
               ArgZ=-t+pi
          endIF
       ELSEIF (ReZ.EQ.zero)   THEN
             IF (ImZ<zero)    THEN
               ArgZ= three*halfpi
         ELSEIF (ImZ.EQ.zero) THEN
               STOP 'Z=0 in ArgZ'
                              ELSE
               ArgZ= halfpi
          endIF
                         ELSE
             IF (ImZ<zero)    THEN
               ArgZ=-t+twopi
         ELSEIF (ImZ.EQ.zero) THEN
               ArgZ= zero
                              ELSE
               ArgZ= t
          endIF
        endIF
      
      RETURN
         
      END FUNCTION ArgZ
