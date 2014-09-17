************************************************************************
      SUBROUTINE AN_Honda11
************************************************************************
*                                                                      *
*                         HONDA-2011 FLUX                              *
*                                                                      *
*                Energy range is 0,1 GeV to 10 TeV.                    *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [ 1] M. Honda, T. Kajita, K. Kasahara, S. Midorikawa.            *
*          "Improvement of low energy atmospheric neutrino flux        *
*          calculation using the JAM nuclear interaction model,"       *
*          Phys. Rev. D 83 (2011) 123001 [astro-ph/1102.2688].         *
*                                                                      *
*     In source the values of  atmospheric neutrino  fluxes  shown     *
*     (arrays F_en, F_ea, F_mn, F_ma) as dF/dE 1/(m^2 sec sr GeV).     *                                                    *
*     Zenith angle bins cos = -1.0-(-0.9), ..., 0.0-0.1, 0.9-1.0.      *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: one

         IMPLICIT REAL (A-H,K-M,O-Z), INTEGER (I-J,N)
                  
         SAVE

         LOGICAL(2),PARAMETER::
     #              Quiz    = .FALSE.
            INTEGER,PARAMETER::
     #              NC      =   20,
     #              Np      =   12,
     #              NE      =  101,
     #              NCoef   =(NC+2)*(NE+2),
     #              Issue_en=    1,
     #              Issue_ea=    2,
     #              Issue_mn=    3,
     #              Issue_ma=    4,
     #              Mode    =    2,                                            Interpolation mode
     #              L       =    1
               REAL,PARAMETER::
     #              C_min   =-9.5d-01,
     #              C_max   = 9.5d-01,
     #              E_min   = 1.0d-01,
     #              E_max   = 1.0d+04

         DIMENSION E(NE),C(NC),
     #             F_en(NE,NC),C_en(NCoef),
     #             F_ea(NE,NC),C_ea(NCoef),
     #             F_mn(NE,NC),C_mn(NCoef),
     #             F_ma(NE,NC),C_ma(NCoef)

         DIMENSION phi(Np),Faz_en(NE,Np,NC),
     #                     Faz_ea(NE,Np,NC),
     #                     Faz_mn(NE,Np,NC),
     #                     Faz_ma(NE,Np,NC)

         OPEN (Ndat00,FILE=datACN//'AN_Honda11_az.data',STATUS='old',
     #                                                  ACTION='read')   f90
         DO n_NC=0,NC-1
           DO n_Np=1,Np
             READ(Ndat00,*) C(NC-n_NC),phi(n_Np)
             DO n_NE=1,NE
               READ(Ndat00,*) E(n_NE),
     #                        Faz_mn(n_NE,n_Np,NC-n_NC),
     #                        Faz_ma(n_NE,n_Np,NC-n_NC),
     #                        Faz_en(n_NE,n_Np,NC-n_NC),
     #                        Faz_ea(n_NE,n_Np,NC-n_NC)
          endDO
        endDO
      endDO
         CLOSE (Ndat00)
         PRINT *,' File AN_Honda11_az.data was read '

         F_en  = sum(Faz_en,2)/Np
         F_ea  = sum(Faz_ea,2)/Np
         F_mn  = sum(Faz_mn,2)/Np
         F_ma  = sum(Faz_ma,2)/Np

         lgE_min= log10(E_min)
         lgE_max= log10(E_max)
         steplgE= (lgE_min-lgE_max)/(NE-1)
         stepC  = (C_max-C_min)/(NC-1)
         DO n_NE=1,NE
           E(n_NE)= ten**(lgE_min+(n_NE-1)*steplgE)
           DO n_NC=1,NC
             C(n_NC)=C_min+(n_NC-1)*stepC
        endDO
      endDO
        
         CALL Coeff2(Mode,Issue_en,NE,NC,lnE_min,C_min,
     #                                lnE_max,C_max,F_en,C_en,Quiz,L)
         CALL Coeff2(Mode,Issue_ea,NE,NC,lnE_min,C_min,
     #                                lnE_max,C_max,F_ea,C_ea,Quiz,L)
         CALL Coeff2(Mode,Issue_mn,NE,NC,lnE_min,C_min,
     #                                lnE_max,C_max,F_mn,C_mn,Quiz,L)
         CALL Coeff2(Mode,Issue_ma,NE,NC,lnE_min,C_min,
     #                                lnE_max,C_max,F_ma,C_ma,Quiz,L)
      
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_en(F0,E0,C0)
*     ==================================================================
         F0=Sp2(Issue_en,C_en,log10(E0),C0)
         RETURN
*     ==================================================================
      ENTRY AN_Honda11_ea(F0,E0,C0)
*     ==================================================================
         F0=Sp2(Issue_ea,C_ea,log10(E0),C0)
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_mn(F0,E0,C0)
*     ==================================================================
         F0=Sp2(Issue_mn,C_mn,log10(E0),C0)
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_ma(F0,E0,C0)
*     ==================================================================
         F0=Sp2(Issue_ma,C_ma,log10(E0),C0)
         RETURN
*     ==================================================================

      END SUBROUTINE AN_Honda11
