************************************************************************
      SUBROUTINE AN_Honda11
************************************************************************
*                                                                      *
*                         HONDA-2011 FLUX                              *
*                      in 1/(cm^2 sec sr GeV)                          *
*                Energy range is 0,1 GeV to 10 TeV.                    *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [ 1] M. Honda,  T. Kajita,  K. Kasahara, and  S. Midorikawa,     *
*          "Improvement  of low energy  atmospheric neutrino  flux     *
*          calculation  using the JAM nuclear  interaction model,"     *
*          Phys. Rev. D83 (2011) 123001 [arXiv:astro-ph/1102.2688].    *
*                                                                      *
*     In source the values of  atmospheric neutrino  fluxes  shown     *
*     (arrays F_en, F_ea, F_mn, F_ma) as dF/dE 1/(m^2 sec sr GeV).     *
*     Zenith angle bins cos = -1.0-(-0.9), ..., 0.0-0.1, 0.9-1.0.      *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: ten

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

         LOGICAL(2),PARAMETER::
     #              TEST   =.TRUE.,                                      Test of spline in reference points
     #              Spectra=.TRUE.                                       Test of spline for energy spectra
            INTEGER,PARAMETER::
     #              Np     =  12,
     #              NE     = 101,
     #              NC     =  20,
     #              NE_test= 101,
     #              NC_test=  20

         DIMENSION E(NE),C(NC),
     #             F_en(NE,NC),C_en(NE,NC),SF_en(NC_test),
     #             F_ea(NE,NC),C_ea(NE,NC),SF_ea(NC_test),
     #             F_mn(NE,NC),C_mn(NE,NC),SF_mn(NC_test),
     #             F_ma(NE,NC),C_ma(NE,NC),SF_ma(NC_test)

         DIMENSION phi(Np),Faz_en(NE,Np,NC),
     #                     Faz_ea(NE,Np,NC),
     #                     Faz_mn(NE,Np,NC),
     #                     Faz_ma(NE,Np,NC)

         OPEN (Ndat00,FILE=DatACN//'AN_Honda11.data',STATUS='old',
     #                                               ACTION='read')      f90
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
         PRINT *,' File AN_Honda11.data was read '
         
         Faz_mn=Faz_mn*1.0d-04
         Faz_ma=Faz_ma*1.0d-04
         Faz_en=Faz_en*1.0d-04
         Faz_ea=Faz_ea*1.0d-04
         
         F_en  = sum(Faz_en,2)/Np
         F_ea  = sum(Faz_ea,2)/Np
         F_mn  = sum(Faz_mn,2)/Np
         F_ma  = sum(Faz_ma,2)/Np

         lgE_ini= log10(1.0d-01)
         lgE_fin= log10(1.0d+04)
         steplgE= (lgE_fin-lgE_ini)/(NE-1)

         cos_ini= -9.5d-01
         cos_fin=  9.5d-01
         stepcos= (cos_fin-cos_ini)/(NC-1)

!         DO n_NE=1,NE
!           E(n_NE)= ten**(lgE_ini+(n_NE-1)*steplgE)
!           DO n_NC=1,NC
!             C(n_NC)= cos_ini+(n_NC-1)*stepcos
!        endDO
!      endDO

         CALL Splie2_mod(E,C,F_en,NE,NC,C_en,TEST)
         CALL Splie2_mod(E,C,F_ea,NE,NC,C_ea,TEST)
         CALL Splie2_mod(E,C,F_mn,NE,NC,C_mn,TEST)
         CALL Splie2_mod(E,C,F_ma,NE,NC,C_ma,TEST)

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN
           OPEN (Ndat01,FILE=OutACN//'Honda/FE_Honda_en.dat')
           OPEN (Ndat02,FILE=OutACN//'Honda/FE_Honda_ea.dat')
           OPEN (Ndat03,FILE=OutACN//'Honda/FE_Honda_mn.dat')
           OPEN (Ndat04,FILE=OutACN//'Honda/FE_Honda_ma.dat')
           DO n_NE=1,NE
             WRITE(Ndat01,1) E(n_NE),(F_en(n_NE,n_NC),n_NC=1,NC)
             WRITE(Ndat02,1) E(n_NE),(F_ea(n_NE,n_NC),n_NC=1,NC)
             WRITE(Ndat03,1) E(n_NE),(F_mn(n_NE,n_NC),n_NC=1,NC)
             WRITE(Ndat04,1) E(n_NE),(F_ma(n_NE,n_NC),n_NC=1,NC)
        endDO
           CLOSE(Ndat01)
           CLOSE(Ndat02)
           CLOSE(Ndat03)
           CLOSE(Ndat04)

           OPEN (Ndat01,FILE=OutACN//'Honda/SE_Honda_en.dat')
           OPEN (Ndat02,FILE=OutACN//'Honda/SE_Honda_ea.dat')
           OPEN (Ndat03,FILE=OutACN//'Honda/SE_Honda_mn.dat')
           OPEN (Ndat04,FILE=OutACN//'Honda/SE_Honda_ma.dat')
           steplgE= (lgE_fin-lgE_ini)/(NE_test-1)
           stepcos= (cos_fin-cos_ini)/(NC_test-1)

           DO n_NE=1,NE_test
             lgEnergy= lgE_ini+(n_NE-1)*steplgE
             Energy  = ten**lgEnergy
             DO n_NC=1,NC_test
               Cosine     = cos_ini+(n_NC-1)*stepcos

               SF_en(n_NC)= Splin2_mod
     #                     (E,C,F_en,C_en,NE,NC,Energy,Cosine)
               SF_ea(n_NC)= Splin2_mod
     #                     (E,C,F_ea,C_ea,NE,NC,Energy,Cosine)
               SF_mn(n_NC)= Splin2_mod
     #                     (E,C,F_mn,C_mn,NE,NC,Energy,Cosine)
               SF_ma(n_NC)= Splin2_mod
     #                     (E,C,F_ma,C_ma,NE,NC,Energy,Cosine)
          endDO
             WRITE(Ndat01,2) Energy,(SF_en(n_NC),n_NC=1,NC_test)
             WRITE(Ndat02,2) Energy,(SF_ea(n_NC),n_NC=1,NC_test)
             WRITE(Ndat03,2) Energy,(SF_mn(n_NC),n_NC=1,NC_test)
             WRITE(Ndat04,2) Energy,(SF_ma(n_NC),n_NC=1,NC_test)
        endDO
           CLOSE(Ndat01)
           CLOSE(Ndat02)
           CLOSE(Ndat03)
           CLOSE(Ndat04)
      endIF
*        ============================================================= *

    1 FORMAT(1PE10.3,20(1PE12.4))
    2 FORMAT(1PE10.3,20(1PE12.4))

         RETURN

*     ==================================================================
      ENTRY AN_Honda11_en(F0,E0,C0)
*     ==================================================================
         F0=Splin2_mod(E,C,F_en,C_en,NE,NC,E0,C0)
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_ea(F0,E0,C0)
*     ==================================================================
         F0=Splin2_mod(E,C,F_ea,C_ea,NE,NC,E0,C0)
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_mn(F0,E0,C0)
*     ==================================================================
         F0=Splin2_mod(E,C,F_mn,C_mn,NE,NC,E0,C0)
         RETURN

*     ==================================================================
      ENTRY AN_Honda11_ma(F0,E0,C0)
*     ==================================================================
         F0=Splin2_mod(E,C,F_ma,C_ma,NE,NC,E0,C0)
         RETURN
*     ==================================================================

      END SUBROUTINE AN_Honda11
