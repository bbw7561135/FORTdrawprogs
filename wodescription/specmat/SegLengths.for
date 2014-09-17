************************************************************************
      SUBROUTINE SegLengths(Length)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2012/11/15  *
************************************************************************
*                                                                      *
*     This SUBROUTINE calculates lengths of the segments which         *
*     sublayer radii cut on the neutrino path chord                    *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                Length(0:KN),
     #                rad(0:KN),R
              INTEGER
     #                j,NDet,Nr
     
               COMMON /radii/  R
               COMMON /rad/    rad
               COMMON /Numbers/NDet,Nr
               
         DO j=KN,Nr+1,-1
           Length(j)=sqrt(rad(j)**2-R**2)-sqrt(rad(j-1)**2-R**2)
      endDO
      
      Length(Nr)=sqrt(rad(Nr)**2-R**2)
     
      RETURN
         
      END SUBROUTINE SegLengths
