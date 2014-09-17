************************************************************************
      FUNCTION dFACN_dE(P,C,N)
************************************************************************
*                                                                      *
*     This FUNCTION returns the values of atmospheric electron and     *
*     muon (anti)neutrino fluxes dF/dE [1/(cm^2 s sr GeV)] for any     *
*     any E_\nu from range 50 MeV to 10^9 GeV and zenith angle.        *
*                                                                      *
*            Data neutrino energy range is 50 MeV to 1 EeV:            *
*                                                                      *
*    0     1   10    10^2  10^3  10^4  10^5  10^6  10^7  10^8  10^9    *
*    |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|     *
*     |----|-----|---|-------|-----|-----|-----|-----|-----------|     *
*    0.05  1   10   70     10^3  10^4  10^5  10^6  10^7        10^9    *
*                            |--------------HE_ISU---------------|     *
*     |---------------------------CORT---------------------------|     *
*      |----------Honda11----------|                                   *
*                                                                      *
*                                                                      *
*     cosine range                                                     *
* |-----------------|                                                  *
* | HE_ISU |Honda11 |                                                  *
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
* |  0.05  |  0.05  |                                                  *
* |  0.10  |        |                                                  *
* |  0.15  |  0.15  |                                                  *
* |  0.20  |        |                                                  *
* |  0.25  |  0.25  |                                                  *
* |  0.30  |        |                                                  *
* |  0.35  |  0.35  |                                                  *
* |  0.40  |        |                                                  *
* |  0.45  |  0.45  |                                                  *
* |  0.50  |        |                                                  *
* |  0.55  |  0.55  |                                                  *
* |  0.60  |        |                                                  *
* |  0.65  |  0.65  |                                                  *
* |  0.70  |        |                                                  *
* |  0.75  |  0.75  |                                                  *
* |  0.80  |        |                                                  *
* |  0.85  |  0.85  |                                                  *
* |  0.90  |        |                                                  *
* |  0.95  |  0.95  |                                                  *
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
     #                HE_ISU_ini=1.00d+06, HE_ISU_fin=1.0d+09,           HE_ISU_ini=1.00d+03
     #                Honda11_ini=1.0d-01,Honda11_fin=1.0d+04,
     #                E07=5.0d+03, E08=5.0d+03, E11=5.0d+00,E12=5.0d+00,
     #                E13=5.0d+00, E14=5.0d+00, E15=5.0d+00,E16=5.0d+00

         PRINT *,' ATMOSPHERIC NEUTRINO MODEL : '
         SELECTCASE(N)
               CASE(1);PRINT *,' AN_Honda11 + AN_HE_ISU + AN_SHE '
                       CALL AN_Honda11                                   0.1  - 10^4 (100 MeV - 10 TeV)
                       CALL AN_HE_ISU13                                  10^3 - 10^9 (  1 TeV -  1 EeV)
*              CASE(2);PRINT *,' CORTout '
!*                      CALL AN_CORTout(Exp,SA,PNmod)                     0.05 - 10^9 (50 MeV -   1 EeV)    
      endSELECT
         dFACN_dE=one
         RETURN

*     ==================================================================
      ENTRY dFAen_dE(P,C,N)
*     ==================================================================
          IF (P.gt.E_cut) P=E_cut
          GOTO (101,102,103) N

  101     IF (P< HE_ISU_ini) THEN; CALL AN_Honda11_en(dFAen_dE,P,C)
      ELSEIF (P<=HE_ISU_fin) THEN; CALL AN_HE_ISU_e(dFAen_dE,P,C)
                             ELSE; CALL AN_HE_ISU_C_e(dFAen_dE,P,C)
       endIF; RETURN
  102     STOP 'CORT' 
  103     STOP 'THIS CASE IS UNDER CONSTRUCTION'
          RETURN

*     ==================================================================
      ENTRY dFAea_dE(P,C,N)
*     ==================================================================
          IF (P.gt.E_cut) P=E_cut
          GOTO (201,202,203) N

  201     IF (P< HE_ISU_ini) THEN; CALL AN_Honda11_ea(dFAea_dE,P,C)
      ELSEIF (P<=HE_ISU_fin) THEN; CALL AN_HE_ISU_e(dFAea_dE,P,C)
                             ELSE; CALL AN_HE_ISU_C_e(dFAea_dE,P,C)
       endIF; RETURN
  202     STOP 'CORT' 
  203     STOP 'THIS CASE IS UNDER CONSTRUCTION'
          RETURN

*     ==================================================================
      ENTRY dFAmn_dE(P,C,N)
*     ==================================================================
          IF (P.gt.E_cut) P=E_cut
          GOTO (301,302,303) N

  301     IF (P< HE_ISU_ini) THEN; CALL AN_Honda11_mn(dFAmn_dE,P,C)
      ELSEIF (P<=HE_ISU_fin) THEN; CALL AN_HE_ISU_m(dFAmn_dE,P,C)
                             ELSE; CALL AN_HE_ISU_C_m(dFAmn_dE,P,C)
       endIF; RETURN
  302     STOP 'CORT' 
  303     STOP 'THIS CASE IS UNDER CONSTRUCTION'
          RETURN

*     ==================================================================
      ENTRY dFAma_dE(P,C,N)
*     ==================================================================
          IF (P.gt.E_cut) P=E_cut
          GOTO (401,402,403) N

  401     IF (P< HE_ISU_ini) THEN; CALL AN_Honda11_ma(dFAma_dE,P,C)
      ELSEIF (P<=HE_ISU_fin) THEN; CALL AN_HE_ISU_m(dFAma_dE,P,C)
                             ELSE; CALL AN_HE_ISU_C_m(dFAma_dE,P,C)
       endIF; RETURN
  402     STOP 'CORT' 
  403     STOP 'THIS CASE IS UNDER CONSTRUCTION'
          RETURN
*     ==================================================================

      END FUNCTION dFACN_dE
