************************************************************************
      SUBROUTINE AN_HE_ISU13
************************************************************************
*                                                                      *
*                         HE_ISU-2013 FLUX                             *
*                      in 1/(cm^2 sec sr GeV)                          *
*                  Energy range is 100 GeV to 100 EeV.                 *
*                                                                      *
*     HE_ISU (Hilas-Gaisser+QGS) fluxes of atmospheric neutrinos and   *
*     antineutrinos as a function of energy and zenith angle.          *
*                                                                      *
*     The source data (arrays F_en, F_ea, F_mn and F_ma) are dF/dE     *
*     in units of 1/(cm^2 s sr GeV).  The  energy reference points     *
*     are calculated on the equidistant over log10(E) grid and the     *
*     reference points for cos(theta)  [where theta  is the zenith     *
*     angle] are calculated on equidistant grid with step = 0.1.       *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: zero,one,ten,E_cut

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

           LOGICAL(2),PARAMETER::
     #                Test   =.TRUE.,                                    Test of spline in reference points
     #                Spectra=.TRUE.,                                    Test of spline for energy spectra
     #                ZenDist=.TRUE.,                                    Test of spline for Z-A distributions
     #                Ratio  =.TRUE.                                     Test of spline for nu/antinu ratios
              INTEGER,PARAMETER::
     #                NE_test=  60,
     #                NC_test=  10,
     #                IncrE  =  10,
     #                K      =   3,
     #                NC     =  11,
     #                NE     =  61,
     #                NE_SHE =  11,
     #                Nfl    =   2
              INTEGER
     #                Ndat(0:2)/200,201,202/
                 REAL,PARAMETER::
     #                E_ini  = 1.0d+02,
     #                E_fin  = 1.0d+08
                 REAL
     #                gamm(Nfl),F00(Nfl),
     #                C(NC),E(NE),logE(NE),
     #                T(NC+1,NE+1),F(Nfl,NC,NE),
     #                lnF_e(NC,NE),lnF_m(NC,NE),
     #                lnF2e(NC,NE),lnF2m(NC,NE),
     #                SF(Nfl,NE),SF0(Nfl,NE),Rs(Nfl,NE),
     #                E_SHE(NE_SHE),logE_SHE(NE_SHE),                    DIMENSIONs for SHE
     #                F_SHE(Nfl,NC,NE_SHE),
     #                lnF_SHE_e(NC,NE_SHE),lnF_SHE_m(NC,NE_SHE),
     #                lnF2SHE_e(NC,NE_SHE),lnF2SHE_m(NC,NE_SHE)
         CHARACTER(*),PARAMETER::
     #                dir='HE_ISU/',
     #                ext='.dat'
         CHARACTER*1
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/

         OPEN(Ndat(0),FILE=datACN//'AN_HE_ISU13.data')
         READ(Ndat(0),*) T
         DO n_NC=1,NC
           C(n_NC)=T(NC+2-n_NC,1)
      endDO
         DO n_NE=1,NE
           E(n_NE)=T(1,n_NE+1)
           DO n_NC=1,NC
             F(1,n_NC,n_NE)=T(NC+2-n_NC,n_NE+1)
        endDO
      endDO
         READ(Ndat(0),*) T
         DO n_NE=1,NE
           DO n_NC=1,NC
             F(2,n_NC,n_NE)=T(NC+2-n_NC,n_NE+1)
        endDO
      endDO
         CLOSE(Ndat(0))
         PRINT *,' File AN_HE_ISU13.data was read '

         lgEmin  =log10(E_ini)
         lgEmax  =log10(E_fin)
         steplgE =(lgEmax-lgEmin)/(NE-1)
         steplgEt=(lgEmax-lgEmin)/NE_test
         lgEcut  =log10(E_cut)
         stepC   =one/(NC-1)
         logE    =log10(E)
         DO n_NC=1,NC
           DO n_NE=1,NE
             lnF_e(n_NC,n_NE)=log(F(1,n_NC,n_NE))
             lnF_m(n_NC,n_NE)=log(F(2,n_NC,n_NE))
        endDO
      endDO

         CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_e,NC,NE,lnF2e,
     #                                                            Test)
         CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_m,NC,NE,lnF2m,
     #                                                            Test)

         steplgE_SHE=(lgEcut-lgEmax)/(NE_SHE-1)
         DO n_NE=1,NE_SHE
           logE_SHE(n_NE)=lgEmax+steplgE_SHE*(n_NE-1)
           E_SHE(n_NE)=ten**logE_SHE(n_NE)
      endDO

         DO n_NC=1,NC
           DO n_NE=1,NE_SHE
             DO n_Nfl=1,Nfl
               CALL PowerLaw(E(NE-1),E(NE),F(n_Nfl,n_NC,NE-1),
     #                         F(n_Nfl,n_NC,NE),gamm(n_Nfl),F00(n_Nfl))
               F_SHE(n_Nfl,n_NC,n_NE)=AN_Spec_cut(E_SHE(n_NE))*
     #                           F00(n_Nfl)*E_SHE(n_NE)**(-gamm(n_Nfl))
          endDO
             lnF_SHE_e(n_NC,n_NE)=log(F_SHE(1,n_NC,n_NE))
             lnF_SHE_m(n_NC,n_NE)=log(F_SHE(2,n_NC,n_NE))
        endDO
      endDO

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN                                               Test for energy spectra
           WRITE(*,*) 'Emax and output energies [GeV] for AN_HE_ISU:'
           WRITE(*,3)  E_max,(E(n_NE),n_NE=1,NE,IncrE)
           DO n_Nfl=1,Nfl
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'FE_'//fln(n_Nfl)//ext)
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             Scale =Energy**K
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),1) Energy,(Scale*F(n_Nfl,n_NC,n_NE),
     #                                                       n_NC=1,NC)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'SE_'//fln(n_Nfl)//ext)
        endDO
           DO n_NE=0,NE_test
             logEt =lgEmin+steplgEt*n_NE
             Energy=ten**logEt
             Scale =Energy**K
             DO n_NC=1,NC
               Cosine=C(n_NC)
               SF(1,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_e,lnF2e,NC,NE,Cosine,logEt))
               SF(2,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_m,lnF2m,NC,NE,Cosine,logEt))
          endDO
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),1) Energy,(Scale*SF(n_Nfl,n_NC),
     #                                                       n_NC=1,NC)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
        endDO
      endIF
*        ============================================================= *
*        TEST FOR ZENITH-ANGLE DISTRIBUTIONS                           *
*        ============================================================= *
         IF (ZenDist) THEN                                               Test for zenith-angle distributions
           DO n_Nfl=1,Nfl
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'FZ_'//fln(n_Nfl)//ext)
        endDO
           DO n_NC=1,NC
             Cosine=C(n_NC)
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),3) Cosine,(F(n_Nfl,n_NC,n_NE)/
     #                                F(n_Nfl,NC,n_NE),n_NE=1,NE,IncrE)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'SZ_'//fln(n_Nfl)//ext)
        endDO
           stepCt=1.0/NC_test
           DO n_NE=1,NE,IncrE                                            Calculating the vertical fluxes
             logEt=logE(n_NE)
             SF0(1,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_e,
     #                                          lnF2e,NC,NE,one,logEt))
             SF0(2,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_m,
     #                                          lnF2m,NC,NE,one,logEt))
        endDO
           DO n_NC=0,NC_test
             Cosine=stepCt*n_NC
             DO n_NE=1,NE,IncrE
               logEt=logE(n_NE)
               SF(1,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_e,
     #                                       lnF2e,NC,NE,Cosine,logEt))
               SF(2,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_m,
     #                                       lnF2m,NC,NE,Cosine,logEt))
          endDO
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),3) Cosine,(SF(n_Nfl,n_NE)/
     #                                 SF0(n_Nfl,n_NE),n_NE=1,NE,IncrE)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
        endDO
      endIF

*        ============================================================= *
*        TEST FOR NEUTRINO/ANTINEUTRINO RATIO                          *
*        ============================================================= *
         IF (Ratio) THEN                                                 Test for neutrino/antineutrino ratio
           DO n_Nfl=1,Nfl
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'FR_'//fln(n_Nfl)//ext)
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),1) Energy,(F(n_Nfl,n_NC,n_NE)/
     #                                    F(n_Nfl,n_NC,n_NE),n_NC=1,NC)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
             OPEN(Ndat(n_Nfl),FILE=outACN//dir//'SR_'//fln(n_Nfl)//ext)
        endDO
           DO n_NE=0,NE_test
             logEt=lgEmin+steplgEt*n_NE
             Energy=ten**logEt
             DO n_NC=1,NC
               Cosine=C(n_NC)
               Rs(1,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_e,lnF2e,NC,NE,Cosine,logEt)-
     #                        Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_e,lnF2e,NC,NE,Cosine,logEt))
               Rs(2,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_m,lnF2m,NC,NE,Cosine,logEt)-
     #                        Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_m,lnF2m,NC,NE,Cosine,logEt))
          endDO
             DO n_Nfl=1,Nfl
               WRITE(Ndat(n_Nfl),1) Energy,(Rs(n_Nfl,n_NC),n_NC=1,NC)
          endDO
        endDO
           DO n_Nfl=1,Nfl
             CLOSE(Ndat(n_Nfl))
        endDO
      endIF
*        ============================================================= *

         RETURN
    1 FORMAT(1PE9.3,21(1PE11.4))
    2 FORMAT(1PE9.3, 8(1PE11.4))
    3 FORMAT(1PE9.3,12(1PE11.4))
    4 FORMAT(12(1PE11.4))

*     ==================================================================
      ENTRY AN_HE_ISU_e(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_e,lnF2e,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_HE_ISU_m(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_m,lnF2m,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

************************************************************************
*                                                                      *
*     Entries  for approximation  AN  fluxies  polinomial function     *
*     F0*Energy^(-Gamma) and employment "GZK" cutoff. (AN_SHE.for)     *
*                                                                      *
************************************************************************
*     ==================================================================
      ENTRY AN_HE_ISU_C_e(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_e,
     #                          lnF2SHE_a,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_HE_ISU_C_m(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_m,
     #                          lnF2SHE_m,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN

  100 FORMAT((1PE10.8)$)

      END SUBROUTINE AN_HE_ISU13