************************************************************************
      SUBROUTINE PowerLaw(E1,E2,F1,F2,gamm,F00)
************************************************************************
*                                                                      *
*     Atmospheric neutrino and antineutrino fluxes is polinomial       *
*     approximation F00*Energy^(-gamm) with employment "GZK" cutoff    *
*                                                                      *
*             Energy range is 10^9 GeV to 3*10^10 GeV.                 *
*                                                                      *
*     AN_SHE (Super High Energy) consists of two parts:                *
*                                                                      *
*     SUBROUTINE PowerLaw initials the parameters gamm, F00 for each   *
*     neutrino type.                                                   *
*                                                                      *
*     FUNCTION AN_Spec_SHE returns values of AN fluxes with employment *
*     "GZK" cutoff.                                                    *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

            INTEGER,PARAMETER::
     #              Nfl = 2
               REAL,PARAMETER::
     #              gamm_min=3
               REAL
     #              E1,F1,E2,F2,
     #              gamm,
     #              F00

         IF (E1.EQ.E2.OR.E1*F2.EQ.zero)
     #   STOP 'Mistake in SUBROUTINE PowerLaw'
         gamm=log(F1/F2)/log(E2/E1)
         IF (gamm.LE.gamm_min)
     #   STOP 'Warning from SUBROUTINE PowerLaw: a rum gamm!'

         F00=F1*E1**gamm
         RETURN

      END SUBROUTINE PowerLaw

************************************************************************
      FUNCTION AN_Spec_SHE(E)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one

         IMPLICIT NONE

         SAVE cut

               REAL,PARAMETER::
     #              E_cut=3.0d+10                                        "GZK"-cutoff energy
               REAL
     #              AN_Spec_SHE,AN_Spec_cut,cut,E

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