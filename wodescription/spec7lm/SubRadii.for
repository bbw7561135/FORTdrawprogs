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
     #                Nl = 7                                            Number of Earth layers by PREM
                 REAL
     #                radius(0:Nl)/0.0000d0,
     #						       1221.5d3,
     #						       3480.0d3,
     #						       5701.0d3,
     #				    	       6346.6d3,
     #						       6356.0d3,
     #						       6368.0d3,
     #						       6371.0d3/,
     #                rad(0:KN),step
              INTEGER
     #                Nsl(Nl)/N1,N2,N3,N4,N5,N6,N7/,
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
