************************************************************************
      SUBROUTINE EarthRho(Rho)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2012/11/15  *
************************************************************************
*                                                                      *
*     This SUBROUTINE calculates density of the Earth matter in each   *
*     sublayer by PREM                                                 *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                PREM(4,0:KN),
     #                Rho(0:KN),RadMid,R,
     #                rad(0:KN)
              INTEGER
     #                j,NDet,Nr,beta(0:KN)
     
               COMMON /PREM/    PREM,beta
               COMMON /radii/   R
               COMMON /rad/     rad
               COMMON /Numbers/ NDet,Nr
         
         DO j=KN,Nr+1,-1
           RadMid=sqrt((half*(sqrt(rad(j)**2-R**2)+
     #            sqrt(rad(j-1)**2-R**2)))**2+R**2)

           Rho(j)=sixth*(PREM(1,j)+PREM(2,j)*(rad(j)/R_Earth)+
     #                             PREM(3,j)*(rad(j)/R_Earth)**2+
     #                             PREM(4,j)*(rad(j)/R_Earth)**3+
     
     #         beta(j)*(B2+B12/(1.0d0+exp((rad(j)/R_Earth-x0 )/dx)))+
     
     #                   PREM(1,j)+PREM(2,j)*(rad(j-1)/R_Earth)+
     #                             PREM(3,j)*(rad(j-1)/R_Earth)**2+
     #                             PREM(4,j)*(rad(j-1)/R_Earth)**3+
     
     #         beta(j)*(B2+B12/(1.0d0+exp((rad(j-1)/R_Earth-x0 )/dx)))+
     
     #             four*(PREM(1,j)+PREM(2,j)*(RadMid/R_Earth)+
     #                             PREM(3,j)*(RadMid/R_Earth)**2+
     #                             PREM(4,j)*(RadMid/R_Earth)**3+
     
     #        beta(j)*(B2+B12/(1.0d0+exp((RadMid/R_Earth-x0 )/dx)))))
      endDO
     
         RadMid=sqrt((half*sqrt(rad(Nr)**2-R**2))**2+R**2)

         Rho(Nr)=sixth*(PREM(1,Nr)+PREM(2,Nr)*(rad(Nr)/R_Earth)+
     #                        PREM(3,Nr)*(rad(Nr)/R_Earth)**2+
     #                        PREM(4,Nr)*(rad(Nr)/R_Earth)**3+
     
     #         beta(Nr)*(B2+B12/(1.0d0+exp((rad(Nr)/R_Earth-x0 )/dx)))+
     
     #                        PREM(1,Nr)+PREM(2,Nr)*(R/R_Earth)+
     #                        PREM(3,Nr)*(R/R_Earth)**2+
     #                        PREM(4,Nr)*(R/R_Earth)**3+
     
     #           beta(Nr)*(B2+B12/(1.0d0+exp((R/R_Earth-x0 )/dx)))+
     
     #             four*(PREM(1,j)+PREM(2,j)*(RadMid/R_Earth)+
     #                             PREM(3,j)*(RadMid/R_Earth)**2+
     #                             PREM(4,j)*(RadMid/R_Earth)**3+
     
     #          beta(j)*(B2+B12/(1.0d0+exp((RadMid/R_Earth-x0 )/dx)))))
     
      RETURN
         
      END SUBROUTINE EarthRho
