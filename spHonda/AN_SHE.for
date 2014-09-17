************************************************************************
      SUBROUTINE PowerLaw_en(E1,E2,F1,F2,Gamma_en,F0_en)
************************************************************************
*                                                                      *
*     Atmospheric neutrinos and antineutrinos fluxes is polinomial     *
*     approximation F0*Energy^(-Gamma) with employment "GZK"cutoff     *
*                                                                      *
*             Energy range is 10^9 GeV to 3 10^10 GeV.                 *
*                                                                      *
*     AN_SHE (Super Higt Energy) consists is two particls:             *
*                                                                      *
*     SUBROUTINE PowerLaw_en with entries for muon neutrinos, ele-     *
*     ctron (muon) antineutrinos, initials  the  parameters Gamma,     *
*     F0 for each neutrino tipe.                                       *
*                                                                      *
*     FUNCTION AN_Spec_SHE  returns  values AN fluxes with  emplo-     *
*     yment "GZK" cutoff.                                              *
*                                                                      *
************************************************************************

         IMPLICIT NONE
                  
         SAVE

         REAL E1,F1,E2,F2,Gamma_en,Gamma_ea,Gamma_mn,Gamma_ma,
     #        F0_en,F0_ea,F0_mn,F0_ma

         REAL,PARAMETER::
     #        Zero=0, G_min=3

         IF (E1.EQ.E2.OR.E1*F2.EQ.Zero)
     #   STOP 'Mistake in SUBR. PowerLaw_en'
         Gamma_en=log(F1/F2)/log(E2/E1)
         IF (Gamma_en.LE.G_min)
     #   STOP 'Warning from SUBROUTINE PowerLaw_en: a rum Gamma_en!'

         F0_en=F1*E1**Gamma_en
         RETURN

*     ==================================================================
      ENTRY PowerLaw_ea(E1,E2,F1,F2,Gamma_ea,F0_ea)
*     ==================================================================
         IF (E1.EQ.E2.OR.E1*F2.EQ.Zero)
     #   STOP 'Mistake in SUBR. PowerLaw_ea'
         Gamma_ea=log(F1/F2)/log(E2/E1)
         IF (Gamma_ea.LE.G_min)
     #   STOP 'Warning from SUBROUTINE PowerLaw_ea: a rum Gamma_ea!'

         F0_ea=F1*E1**Gamma_ea
         RETURN

*     ==================================================================
      ENTRY PowerLaw_mn(E1,E2,F1,F2,Gamma_mn,F0_mn)
*     ==================================================================
         IF (E1.EQ.E2.OR.E1*F2.EQ.Zero)
     #   STOP 'Mistake in SUBR. PowerLaw_mn'
         Gamma_mn=log(F1/F2)/log(E2/E1)
         IF (Gamma_mn.LE.G_min)
     #   STOP 'Warning from SUBROUTINE PowerLaw_mn: a rum Gamma_mn!'

         F0_mn=F1*E1**Gamma_mn
         RETURN

*     ==================================================================
      ENTRY PowerLaw_ma(E1,E2,F1,F2,Gamma_ma,F0_ma)
*     ==================================================================
         IF (E1.EQ.E2.OR.E1*F2.EQ.Zero)
     #   STOP 'Mistake in SUBR. PowerLaw_ma'
         Gamma_ma=log(F1/F2)/log(E2/E1)
         IF (Gamma_ma.LE.G_min)
     #   STOP 'Warning from SUBROUTINE PowerLaw_ma: a rum Gamma_ma!'
         F0_ma=F1*E1**Gamma_ma
         RETURN
*     ==================================================================

      END SUBROUTINE PowerLaw_en

************************************************************************
      FUNCTION AN_Spec_SHE(E)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one

         IMPLICIT NONE

         REAL AN_Spec_SHE,AN_Spec_cut,Cut,E

         SAVE cut

         REAL,PARAMETER:: E_cut=3.0d+10                                  "GZK" cutoff energy

         cut=asin(one)/E_cut
         AN_Spec_SHE=cut
         RETURN

*     ==================================================================
      ENTRY AN_Spec_cut(E)
*     ==================================================================
         AN_Spec_cut=one/(one+tan(cut*E))
         RETURN
*     ==================================================================

      END FUNCTION AN_Spec_SHE