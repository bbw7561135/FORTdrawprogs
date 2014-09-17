************************************************************************
      FUNCTION dFAN_dE(P,C,N)
************************************************************************
*                                                                      *
*     This FUNCTION returns the values of atmospheric neutrino fluxes  *
*     dF/dE [1/(cm^2 s sr GeV)] for any E_\nu from range 50 MeV to     *
*     1 EeV and any zenith angle.                                      *
*                                                                      *
*            Data neutrino energy range is 50 MeV to 1 EeV:            *
*                                                                      *
*    0     1    10   10^2  10^3  10^4  10^5  10^6  10^7  10^8  10^9    *
*    ||----|-----|-----|-----|-----|-----|-----|-----|-----|-----|     *
*    0.05                                                              *
*                      |--------------HE_ISU13-------------|           *
*     |---------------------------CORT---------------------------|     *
*      |----------Honda11----------|                                   *
*                                                                      *
*     cosine range                                                     *
* |-----------------|                                                  *
* |HE_ISU13|Honda11 |                                                  *
* |-----------------|                                                  *
* |   --   |        |                                                  *
* |   --   | -0.95  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.85  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.75  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.65  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.55  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.45  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.35  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.25  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.15  |                                                  *
* |   --   |        |                                                  *
* |   --   | -0.05  |                                                  *
* |  0.00  |        |                                                  *
* |        |  0.05  |                                                  *
* |  0.10  |        |                                                  *
* |        |  0.15  |                                                  *
* |  0.20  |        |                                                  *
* |        |  0.25  |                                                  *
* |  0.30  |        |                                                  *
* |        |  0.35  |                                                  *
* |  0.40  |        |                                                  *
* |        |  0.45  |                                                  *
* |  0.50  |        |                                                  *
* |        |  0.55  |                                                  *
* |  0.60  |        |                                                  *
* |        |  0.65  |                                                  *
* |  0.70  |        |                                                  *
* |        |  0.75  |                                                  *
* |  0.80  |        |                                                  *
* |        |  0.85  |                                                  *
* |  0.90  |        |                                                  *
* |        |  0.95  |                                                  *
* |  1.00  |        |                                                  *
* |-----------------|                                                  *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one,E_cut

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER(3),PARAMETER::
     #                SA   ='max'                                        Level of Solar activity
              INTEGER,PARAMETER::
     #                Exp  =5,                                           Neutrino experiment number
     #                PNmod=2                                            Prompt neutrino production model
                 REAL,PARAMETER::
     #                HE_ISU_ini=1.00d+03, HE_ISU_fin=1.0d+08

         PRINT *,' ATMOSPHERIC NEUTRINO MODEL : '
         SELECTCASE(N)
               CASE(1)
                       PRINT *,' AN_Honda11 + AN_HE_ISU + AN_SHE '
                       CALL AN_Honda11
                       CALL AN_HE_ISU13
               CASE(2)
                       PRINT *,' CORTout '
!                      CALL AN_CORTout(Exp,SA,PNmod)
      endSELECT
         dFACN_dE=one
         RETURN

*     ==================================================================
      ENTRY dFnu_dE(nflavor,ntype,P,C,N)
*     ==================================================================
          IF (P.gt.E_cut) P=E_cut
          GOTO (101,102,103) N

  101     IF (P< HE_ISU_ini) THEN
            CALL AN_Honda(nflavor,ntype,dFA_dE,P,C)
      ELSEIF (P<=HE_ISU_fin) THEN
            CALL AN_HE_ISU(nflavor,dFA_dE,P,C)
                             ELSE 
            CALL AN_SHE_ISU(nflavor,dFA_dE,P,C)
       endIF
          RETURN

  102     STOP 'CORT' 
  103     STOP 'THIS CASE IS UNDER CONSTRUCTION'

      END FUNCTION dFAN_dE