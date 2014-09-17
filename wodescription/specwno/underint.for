************************************************************************
      FUNCTION underint(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
                 REAL
     #                cosin

         COMMON    /N/N
         COMMON/C_lim/C_min,deltaC
         COMMON /E_nu/E_nu

         underint=one
         
         RETURN
         
*     ==================================================================
      ENTRY uifen(var)
*     ==================================================================
         
         cosin=deltaC*var+C_min
         uifen=dFAen_dE(E_nu,cosin,N)                                   !*deltaC
         
         RETURN

*     ==================================================================
      ENTRY uifea(var)
*     ==================================================================
         
         cosin=deltaC*var+C_min
         uifea=dFAea_dE(E_nu,cosin,N)                                   !*deltaC
         
         RETURN

*     ==================================================================
      ENTRY uifmn(var)
*     ==================================================================
         
         cosin=deltaC*var+C_min
         uifmn=dFAmn_dE(E_nu,cosin,N)                                   !*deltaC
         
         RETURN

*     ==================================================================
      ENTRY uifma(var)
*     ==================================================================
         
         cosin=deltaC*var+C_min
         uifma=dFAma_dE(E_nu,cosin,N)                                   !*deltaC
         
         RETURN

      END FUNCTION underint