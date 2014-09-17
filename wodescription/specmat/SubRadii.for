************************************************************************
      SUBROUTINE SubRadii(rad)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2013/02/25  *
************************************************************************
*                                                                      *
*     This SUBROUTINE calculates radii the of Earth sublayers          *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nl = 10                                            Number of Earth layers by PREM
                 REAL
     #                radius(0:Nl)/zero,
     #                            1.2215d+06,
     #                            3.4800d+06,
     #                            5.7010d+06,
     #                            5.7710d+06,
     #                            5.9710d+06,
     #                            6.1510d+06,
     #                            6.3466d+06,
     #                            6.3560d+06,
     #                            6.3680d+06,
     #                            6.3710d+06/,
     #                rad(0:KN),step
              INTEGER
     #                Nsl(Nl)/N1,N2,N3,N4,N5,N6,N7,N8,N9,N10/,
     #                j,n_Nl,n_sl
     
         n_sl=0
         rad(0)=zero
         
         DO n_Nl=1,Nl
           step=(radius(n_Nl)-radius(n_Nl-1))/Nsl(n_Nl)
           DO j=1,Nsl(n_Nl)
             n_sl=n_sl+1
             rad(n_sl)=radius(n_Nl-1)+j*step
        endDO
      endDO
     
         RETURN
         
      END SUBROUTINE SubRadii
