************************************************************************
      FUNCTION GeM_FV_SM(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: pi

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON   /n_SM_TMD/n_SM_TMD                                     Switch for target nucleon momentum distribution
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon

         p=var*P_FeMAX
         GeM_FV_SM=4*pi*p**2*rho_SM(n_SM_TMD,p)*P_FeMAX

         RETURN
      END FUNCTION GeM_FV_SM