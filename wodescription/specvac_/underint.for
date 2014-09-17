************************************************************************
      FUNCTION underint(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
         SAVE
         
                 REAL,PARAMETER::
     #                Re=6.37122d+03,                                    Earth radius!*1.0d+05
     #                Ra=1.5d+01                                         Atmosphere depth!*1.0d+05

         COMMON    /N/N
         COMMON /n_fl/n_fl
         COMMON /n_NT/n_NT
         COMMON/C_lim/C_min,deltaC
         COMMON /E_nu/E_nu
         
         R2 =(Re+Ra)**2
         Re2=Re**2

         underint=one
         
         RETURN
         
*     ==================================================================
      ENTRY uif(var)
*     ==================================================================
         
         cosin=deltaC*var+C_min
         root=sqrt(R2-Re2*(1-cosin**2))
         bite=Re*cosin
         eL=root-bite
         SELECTCASE(n_fl)
               CASE(1)
           IF (n_NT.EQ.1) THEN
             uif=Pab(1,1,E_nu,eL)*dFAen_dE(E_nu,cosin,N)+
     #           Pab(2,1,E_nu,eL)*dFAmn_dE(E_nu,cosin,N)              !*deltaC
                          ELSE
             uif=Pab(1,1,E_nu,eL)*dFAea_dE(E_nu,cosin,N)+
     #           Pab(2,1,E_nu,eL)*dFAma_dE(E_nu,cosin,N)              !*deltaC
        endIF
               CASE(2)
           IF (n_NT.EQ.1) THEN
             uif=Pab(2,2,E_nu,eL)*dFAmn_dE(E_nu,cosin,N)+
     #           Pab(1,2,E_nu,eL)*dFAen_dE(E_nu,cosin,N)              !*deltaC
                          ELSE
             uif=Pab(2,2,E_nu,eL)*dFAma_dE(E_nu,cosin,N)+
     #           Pab(1,2,E_nu,eL)*dFAea_dE(E_nu,cosin,N)              !*deltaC
        endIF
               CASE(3)
           IF (n_NT.EQ.1) THEN
             uif=Pab(1,3,E_nu,eL)*dFAen_dE(E_nu,cosin,N)+
     #           Pab(2,3,E_nu,eL)*dFAmn_dE(E_nu,cosin,N)              !*deltaC
                          ELSE
             uif=Pab(1,3,E_nu,eL)*dFAea_dE(E_nu,cosin,N)+
     #           Pab(2,3,E_nu,eL)*dFAma_dE(E_nu,cosin,N)              !*deltaC
        endIF
      endSELECT
         
         RETURN

      END FUNCTION underint