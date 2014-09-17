************************************************************************
      FUNCTION dFANom_dE(P)
************************************************************************
*                                                                      *
*     This FUNCTION returns the values of atmospheric electron and     *
*     muon (anti)neutrino fluxes dF/dE [1/(cm^2 s sr GeV)] modulated   *
*     neutrino oscillations and averaged over the cosine of            *
*     zenith angle for any E_\nu from range 50 MeV to 10^9 GeV.        *
*                                                                      *
*            Data neutrino energy range is 50 MeV to 1 EeV:            *
*                                                                      *
*    0     1   10    10^2  10^3  10^4  10^5  10^6  10^7  10^8  10^9    *
*    |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|     *
*     |----|-----|---|-------|-----|-----|-----|-----|-----------|     *
*    0.05  1   10   70     10^3  10^4  10^5  10^6  10^7        10^9    *
*                            |--------------ISU_HE---------------|     *
*     |---------------------------CORT---------------------------|     *
*      |----------Honda11----------|                                   *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE
         
           LOGICAL(2),PARAMETER::
     #                Tie  =.TRUE.,
     #                Quiz =.TRUE.
              INTEGER,PARAMETER::
     #                Nfen = 301,
     #                Nfea = 302,              
     #                Nfmn = 303,
     #                Nfma = 304,
     #                Nftn = 305,
     #                Nfta = 306,
     #                Mult =   0,
     #                Mode =   1,
     #                Issen=   1,
     #                Issea=   2,
     #                Issmn=   3,
     #                Issma=   4,
     #                Isstn=   3,
     #                Issta=   4,
     #                NP   = 100,
     #                L    =   1
                 REAL,PARAMETER::
     #                P_min = 1.0d-01,
     #                P_max = 1.0d+04,
     #                Eps   = 1.0d-12
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*57
     #                fnen,fnea,fnmn,fnma,fntn,fnta
         CHARACTER*3
     #                Sp,
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                SA,
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/


         DIMENSION Fen(NP),Fea(NP),Fmn(NP),Fma(NP),Ftn(NP),Fta(NP),
     #             Cen(NP+2),Cea(NP+2),Cmn(NP+2),Cma(NP+2),
     #             Ctn(NP+2),Cta(NP+2)

         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON     /n_DM/n_DM                                           Name of the Earth density model
         COMMON     /n_hi/n_hi                                           Switch for neutrino mass hierarchy

         SELECTCASE(N)
               CASE(1)
                     WRITE(*,*) ' AN_Honda11 + AN_ISU_HE + AN_SHE '
                     Sp='H11'
                     WRITE(*,*) ' Maximal solar activity '
                     SA='x'
               CASE(2)
                     WRITE(*,*) ' CORTout '
                     Sp='CRT'
                     WRITE(*,*) ' Minimal solar activity '
                     SA='n'
      endSELECT
         fnen=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(1)//
     #                                                      NTn(1)//ext
         fnea=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(1)//
     #                                                      NTn(2)//ext
         fnmn=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(2)//
     #                                                      NTn(1)//ext
         fnma=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(2)//
     #                                                      NTn(2)//ext
         fntn=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(3)//
     #                                                      NTn(1)//ext
         fnta=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(3)//
     #                                                      NTn(2)//ext
         OPEN(Nfen,FILE=fnen)
         OPEN(Nfea,FILE=fnea)
         OPEN(Nfmn,FILE=fnmn)
         OPEN(Nfma,FILE=fnma)
         OPEN(Nftn,FILE=fntn)
         OPEN(Nfta,FILE=fnta)

         DO n_NP=1,NP
           READ(Nfen,101) E,Fen(n_NP)
           READ(Nfea,101) E,Fea(n_NP)
           READ(Nfmn,101) E,Fmn(n_NP)
           READ(Nfma,101) E,Fma(n_NP)
           READ(Nftn,101) E,Ftn(n_NP)
           READ(Nfta,101) E,Fta(n_NP)
      endDO
      
         CLOSE(Nfen)
         CLOSE(Nfea)
         CLOSE(Nfmn)
         CLOSE(Nfma)
         CLOSE(Nftn)
         CLOSE(Nfta)
         
         lgP_min = log10(P_min)
         lgP_max = log10(P_max)
         CALL Coeff1(Mult,Mode,Tie,Eps,Issen,NP,lgP_min,lgP_max,
     #               Fen,Cen,Quiz,L)
         CALL Coeff1(Mult,Mode,Tie,Eps,Issea,NP,lgP_min,lgP_max,
     #               Fea,Cea,Quiz,L)
         CALL Coeff1(Mult,Mode,Tie,Eps,Issmn,NP,lgP_min,lgP_max,
     #               Fmn,Cmn,Quiz,L)
         CALL Coeff1(Mult,Mode,Tie,Eps,Issma,NP,lgP_min,lgP_max,
     #               Fma,Cma,Quiz,L)
         CALL Coeff1(Mult,Mode,Tie,Eps,Isstn,NP,lgP_min,lgP_max,
     #               Ftn,Ctn,Quiz,L)
         CALL Coeff1(Mult,Mode,Tie,Eps,Issta,NP,lgP_min,lgP_max,
     #               Fta,Cta,Quiz,L)

         dFANom_dE=one
         RETURN
*     ==================================================================
      ENTRY dFANomen_dE(P)
*     ==================================================================
         dFANomen_dE=Sp1(Issen,Cen,log10(P))
         RETURN
         
*     ==================================================================
      ENTRY dFANomea_dE(P)
*     ==================================================================
         dFANomea_dE=Sp1(Issea,Cea,log10(P))
         RETURN
         
*     ==================================================================
      ENTRY dFANommn_dE(P)
*     ==================================================================
         dFANommn_dE=Sp1(Issmn,Cmn,log10(P))
         RETURN
         
*     ==================================================================
      ENTRY dFANomma_dE(P)
*     ==================================================================
         dFANomma_dE=Sp1(Issma,Cma,log10(P))
         RETURN

*     ==================================================================
      ENTRY dFANomtn_dE(P)
*     ==================================================================
         dFANomtn_dE=Sp1(Isstn,Ctn,log10(P))
         RETURN
         
*     ==================================================================
      ENTRY dFANomta_dE(P)
*     ==================================================================
         dFANomta_dE=Sp1(Issta,Cta,log10(P))
         RETURN
         
  101 FORMAT(2(1PE16.8))
      END FUNCTION dFANom_dE
