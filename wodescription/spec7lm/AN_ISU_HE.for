************************************************************************
      SUBROUTINE AN_ISU_HE
************************************************************************
*                                                                      *
*     ISU_HE fluxes of atmospheric  neutrinos and antineutrinos as     *
*     a function of energy and zenith angle.                           *
*                                                                      *
*                  Energy range is 1 TeV to 1 EeV.                     *
*                                                                      *
*     The source data (arrays F_en, F_ea, F_mn and F_ma) are dF/dE     *
*     in units of 1/(cm^2 s sr GeV).  The  energy reference points     *
*     are calculated on the equidistant over log10(E) grid and the     *
*     reference points for cos(theta)  [where theta  is the zenith     *
*     angle] are calculated on equidistant grid with  step = 0.05.     *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: one,ten,E_cut

         IMPLICIT REAL (A-H,K-Z), INTEGER (I,J)
                  
         SAVE

         LOGICAL(2),PARAMETER::
     #              Test   =.TRUE.,                                      Test of spline in reference points
     #              Spectra=.TRUE.,                                      Test of spline for energy spectra
     #              ZenDist=.TRUE.,                                      Test of spline for Z-A distributions
     #              Ratio  =.TRUE.,                                      Test of spline for nu/antinu ratios
     #              Vgrid  =.FALSE.                                      Output on the grid of Volkova
            INTEGER,PARAMETER::
     #              NE_test=150,
     #              NC_test=100,
     #              IncrE  =  9,
     #              K      =  3,
     #              NC     = 21,
     #              NE     =101,
     #              NE_SHE = 11
               REAL,PARAMETER::
     #              E_min  =1.0d+03,
     #              E_max  =1.0d+09,
     #              Rescale=1.0d-04

         REAL F_en1(7,NE),F_en2(7,NE),F_en3(7,NE),
     #        F_ea1(7,NE),F_ea2(7,NE),F_ea3(7,NE),
     #        F_mn1(7,NE),F_mn2(7,NE),F_mn3(7,NE),
     #        F_ma1(7,NE),F_ma2(7,NE),F_ma3(7,NE)

         REAL SF_en(NE),SF_ea(NE),SF_mn(NE),SF_ma(NE)

         REAL C(NC),E(NE),logE(NE),C_V(8),E_V(15),
     #        F_en(NC,NE),F_ea(NC,NE),Rs_e(NE),
     #        F_mn(NC,NE),F_ma(NC,NE),Rs_m(NE),
     #        lnF_en(NC,NE),lnF2en(NC,NE),SF0en(NE),
     #        lnF_ea(NC,NE),lnF2ea(NC,NE),SF0ea(NE),
     #        lnF_mn(NC,NE),lnF2mn(NC,NE),SF0mn(NE),
     #        lnF_ma(NC,NE),lnF2ma(NC,NE),SF0ma(NE)

*        DIMENSIONs for "super high energies range"
         REAL E_SHE(NE_SHE),logE_SHE(NE_SHE),
     #        F_SHE_en(NC,NE_SHE),lnF_SHE_en(NC,NE_SHE) ,
     #                            lnF2_SHE_en(NC,NE_SHE),
     #        F_SHE_ea(NC,NE_SHE),lnF_SHE_ea(NC,NE_SHE) ,
     #                            lnF2_SHE_ea(NC,NE_SHE),
     #        F_SHE_mn(NC,NE_SHE),lnF_SHE_mn(NC,NE_SHE) ,
     #                            lnF2_SHE_mn(NC,NE_SHE),
     #        F_SHE_ma(NC,NE_SHE),lnF_SHE_ma(NC,NE_SHE) ,
     #                            lnF2_SHE_ma(NC,NE_SHE)

         EQUIVALENCE (SF_en(1),Rs_e(1)),(SF_mn(1),Rs_m(1))

         COMMON    /Volkova/C_V,E_V

         NAMELIST /AN_ISU_HE_data/F_en1,F_en2,F_en3,F_ea1,F_ea2,F_ea3,
     #                            F_mn1,F_mn2,F_mn3,F_ma1,F_ma2,F_ma3
         OPEN (Ndat00,FILE=datACN//'AN_ISU_HE.data',STATUS='old',
     #                                              ACTION='read')       f90
         READ (Ndat00,AN_ISU_HE_data); CLOSE (Ndat00)
         PRINT *,' File AN_ISU_HE.data was read '

         lgEmin     =log10(E_min)
         lgEmax     =log10(E_max)
         lgEcut     =log10(E_cut)
         steplgE    =(lgEmax-lgEmin)/(NE-1)
         steplgE_SHE=(lgEcut-lgEmax)/(NE_SHE-1)
         stepC      =1.0/(NC-1)

         DO i=1,NE
           logE(i)=lgEmin+(i-1)*SteplgE
           E   (i)=ten**logE(i)
           DO j=1,7
             F_en(22-j,i)=F_en1(j,i)*Rescale
             F_en(15-j,i)=F_en2(j,i)*Rescale
             F_en( 8-j,i)=F_en3(j,i)*Rescale
             F_ea(22-j,i)=F_ea1(j,i)*Rescale
             F_ea(15-j,i)=F_ea2(j,i)*Rescale
             F_ea( 8-j,i)=F_ea3(j,i)*Rescale
             F_mn(22-j,i)=F_mn1(j,i)*Rescale
             F_mn(15-j,i)=F_mn2(j,i)*Rescale
             F_mn( 8-j,i)=F_mn3(j,i)*Rescale
             F_ma(22-j,i)=F_ma1(j,i)*Rescale
             F_ma(15-j,i)=F_ma2(j,i)*Rescale
             F_ma( 8-j,i)=F_ma3(j,i)*Rescale
        endDO
           DO j=1,NC
             lnF_en(j,i)=log(F_en(j,i))
             lnF_ea(j,i)=log(F_ea(j,i))
             lnF_mn(j,i)=log(F_mn(j,i))
             lnF_ma(j,i)=log(F_ma(j,i))
        endDO
      endDO
         DO j=1,NC
           C(j)=stepC*(j-1)
      endDO
         DO i=1,NE_SHE
           logE_SHE(i)=lgEmax+stepLgE_SHE*(i-1)
           E_SHE(i)=ten**logE_SHE(i)
      endDO
         DO j=1,NC
           DO i=1,NE_SHE
             CALL PowerLaw_en(E(NE-1),E(NE),F_en(j,NE-1),F_en(j,NE),
     #                                                  Gamma_en,F0_en)
             F_SHE_en(j,i)=AN_Spec_cut(E_SHE(i))*
     #                                      F0_en*E_SHE(i)**(-Gamma_en)
             CALL PowerLaw_ea(E(NE-1),E(NE),F_ea(j,NE-1),F_ea(j,NE),
     #                                                  Gamma_ea,F0_ea)
             F_SHE_ea(j,i)=AN_Spec_cut(E_SHE(i))*
     #                                      F0_ea*E_SHE(i)**(-Gamma_ea)
             CALL PowerLaw_mn(E(NE-1),E(NE),F_mn(j,NE-1),F_mn(j,NE),
     #                                                  Gamma_mn,F0_mn)
             F_SHE_mn(j,i)=AN_Spec_cut(E_SHE(i))*
     #                                      F0_mn*E_SHE(i)**(-Gamma_mn)
             CALL PowerLaw_ma(E(NE-1),E(NE),F_ma(j,NE-1),F_ma(j,NE),
     #                                                  Gamma_ma,F0_ma)
             F_SHE_ma(j,i)=AN_Spec_cut(E_SHE(i))*
     #                                      F0_ma*E_SHE(i)**(-Gamma_ma)
        endDO
      endDO

         lnF_SHE_en=log(F_SHE_en)
         lnF_SHE_ea=log(F_SHE_ea)
         lnF_SHE_mn=log(F_SHE_mn)
         lnF_SHE_ma=log(F_SHE_ma)

      CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_en,NC,NE,lnF2en,Test)
      CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_ea,NC,NE,lnF2ea,Test)
      CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_mn,NC,NE,lnF2mn,Test)
      CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF_ma,NC,NE,lnF2ma,Test)

         steplgEt=(lgEmax-lgEmin)/NE_test

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN                                               Test for energy spectra
           WRITE(*,*) 'Emax and output energies [GeV] for AN_ISU_HE:'
           WRITE(*,3)  E_max,(E(i),i=1,NE,IncrE)
           OPEN (Ndat01,FILE=outACN//'ISU_HE/FE_ISU_HE_en.dat')
           OPEN (Ndat02,FILE=outACN//'ISU_HE/FE_ISU_HE_ea.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/FE_ISU_HE_mn.dat')
           OPEN (Ndat04,FILE=outACN//'ISU_HE/FE_ISU_HE_ma.dat')
           DO i=1,NE
             Energy=E(i)
             Scale =Energy**K
             WRITE(Ndat01,1) Energy, (Scale*F_en(j,i),j=1,NC)
             WRITE(Ndat02,1) Energy, (Scale*F_ea(j,i),j=1,NC)
             WRITE(Ndat03,1) Energy, (Scale*F_mn(j,i),j=1,NC)
             WRITE(Ndat04,1) Energy, (Scale*F_ma(j,i),j=1,NC)
        endDO
           OPEN (Ndat01,FILE=outACN//'ISU_HE/SE_ISU_HE_en.dat')
           OPEN (Ndat02,FILE=outACN//'ISU_HE/SE_ISU_HE_ea.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/SE_ISU_HE_mn.dat')
           OPEN (Ndat04,FILE=outACN//'ISU_HE/SE_ISU_HE_ma.dat')
           DO i=0,NE_test
             lgE   =lgEmin+steplgEt*i
             Energy=ten**lgE
             Scale =Energy**K
             DO j=1,NC
               Cosine  =C(j)
               SF_en(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_en,lnF2en,NC,NE,Cosine,lgE))
               SF_ea(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_ea,lnF2ea,NC,NE,Cosine,lgE))
               SF_mn(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_mn,lnF2mn,NC,NE,Cosine,lgE))
               SF_ma(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_ma,lnF2ma,NC,NE,Cosine,lgE))
          endDO
             WRITE(Ndat01,1) Energy,(Scale*SF_en(j),j=1,NC)
             WRITE(Ndat02,1) Energy,(Scale*SF_ea(j),j=1,NC)
             WRITE(Ndat03,1) Energy,(Scale*SF_mn(j),j=1,NC)
             WRITE(Ndat04,1) Energy,(Scale*SF_ma(j),j=1,NC)
        endDO
      endIF
*        ============================================================= *
*        TEST FOR ZENITH-ANGLE DISTRIBUTIONS                           *
*        ============================================================= *
         IF (ZenDist) THEN                                               Test for zenith-angle distributions
           OPEN (Ndat01,FILE=outACN//'ISU_HE/FZ_ISU_HE_en.dat')
           OPEN (Ndat02,FILE=outACN//'ISU_HE/FZ_ISU_HE_ea.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/FZ_ISU_HE_mn.dat')
           OPEN (Ndat04,FILE=outACN//'ISU_HE/FZ_ISU_HE_ma.dat')
           DO j=1,NC
             Cosine=C(j)
             WRITE(Ndat01,3) Cosine,(F_en(j,i)/F_en(NC,i),i=1,NE,IncrE)
             WRITE(Ndat02,3) Cosine,(F_ea(j,i)/F_ea(NC,i),i=1,NE,IncrE)
             WRITE(Ndat03,3) Cosine,(F_mn(j,i)/F_mn(NC,i),i=1,NE,IncrE)
             WRITE(Ndat04,3) Cosine,(F_ma(j,i)/F_ma(NC,i),i=1,NE,IncrE)
        endDO
           OPEN (Ndat01,FILE=outACN//'ISU_HE/SZ_ISU_HE_en.dat')
           OPEN (Ndat02,FILE=outACN//'ISU_HE/SZ_ISU_HE_ea.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/SZ_ISU_HE_mn.dat')
           OPEN (Ndat04,FILE=outACN//'ISU_HE/SZ_ISU_HE_ma.dat')
           stepCt=1.0/NC_test
           DO i=1,NE,IncrE                                               Calculating the vertical fluxes
             lgE=logE(i)
             SF0en(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                              lnF_en,lnF2en,NC,NE,one,lgE))
             SF0ea(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                              lnF_ea,lnF2ea,NC,NE,one,lgE))
             SF0mn(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                              lnF_mn,lnF2mn,NC,NE,one,lgE))
             SF0ma(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                              lnF_ma,lnF2ma,NC,NE,one,lgE))
        endDO
           DO j=0,NC_test-1
             Cosine=stepCt*j
             DO i=1,NE,IncrE
               lgE=logE(i)
               SF_en(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_en,lnF2en,NC,NE,Cosine,lgE))
               SF_ea(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_ea,lnF2ea,NC,NE,Cosine,lgE))
               SF_mn(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_mn,lnF2mn,NC,NE,Cosine,lgE))
               SF_ma(i)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                lnF_ma,lnF2ma,NC,NE,Cosine,lgE))
          endDO
             WRITE(Ndat01,3) Cosine,(SF_en(i)/SF0en(i),i=1,NE,IncrE)
             WRITE(Ndat02,3) Cosine,(SF_ea(i)/SF0ea(i),i=1,NE,IncrE)
             WRITE(Ndat03,3) Cosine,(SF_mn(i)/SF0mn(i),i=1,NE,IncrE)
             WRITE(Ndat04,3) Cosine,(SF_ma(i)/SF0ma(i),i=1,NE,IncrE)
        endDO
      endIF
*        ============================================================= *
*        TEST FOR NEUTRINO/ANTINEUTRINO RATIO                          *
*        ============================================================= *
         IF (Ratio) THEN                                                 Test for neutrino/antineutrino ratio
           OPEN (Ndat01,FILE=outACN//'ISU_HE/FR_ISU_HE_e.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/FR_ISU_HE_m.dat')
           DO i=1,NE
             Energy=E(i)
             WRITE(Ndat01,1) Energy,(F_en(j,i)/F_ea(j,i),j=1,NC)
             WRITE(Ndat03,1) Energy,(F_mn(j,i)/F_ma(j,i),j=1,NC)
        endDO
           OPEN (Ndat01,FILE=outACN//'ISU_HE/SR_ISU_HE_e.dat')
           OPEN (Ndat03,FILE=outACN//'ISU_HE/SR_ISU_HE_m.dat')
           IF (Vgrid) THEN                                               Output on the grid of Volkova
             DO i=1,15
               Energy=E_V(i)
               lgE=log10(Energy)
               DO j=1,8
                 Cosine=C_V(j)
                 Rs_e(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_en,lnF2en,NC,NE,Cosine,lgE)-
     #                       Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_ea,lnF2ea,NC,NE,Cosine,lgE))
                 Rs_m(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_mn,lnF2mn,NC,NE,Cosine,lgE)-
     #                       Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_ma,lnF2ma,NC,NE,Cosine,lgE))
            endDO
               WRITE(Ndat01,2) Energy,(Rs_e(j),j=1,8)
               WRITE(Ndat03,2) Energy,(Rs_m(j),j=1,8)
          endDO
                      ELSE                                               "Standard" output (ISU_HE grid)
             DO i=0,NE_test
               lgE=lgEmin+steplgEt*i
               Energy=ten**lgE
               DO j=1,NC
                 Cosine=C(j)
                 Rs_e(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_en,lnF2en,NC,NE,Cosine,lgE)-
     #                       Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_ea,lnF2ea,NC,NE,Cosine,lgE))
                 Rs_m(j)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_mn,lnF2mn,NC,NE,Cosine,lgE)-
     #                       Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                                 lnF_ma,lnF2ma,NC,NE,Cosine,lgE))
            endDO
               WRITE(Ndat01,1) Energy,(Rs_e(j),j=1,NC)
               WRITE(Ndat03,1) Energy,(Rs_m(j),j=1,NC)
          endDO
        endIF
      endIF
*        ============================================================= *

         RETURN
    1 FORMAT(1PE9.3,21(1PE11.4))
    2 FORMAT(1PE9.3, 8(1PE11.4))
    3 FORMAT(1PE9.3,12(1PE11.4))

*     ==================================================================
      ENTRY AN_ISU_HE_en (F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_en,lnF2en,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_HE_ea (F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_ea,lnF2ea,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_HE_mn (F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_mn,lnF2mn,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_HE_ma (F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,lnF_ma,lnF2ma,NC,NE,
     #                                               abs(C0),log10(E0)))
         RETURN

************************************************************************
*                                                                      *
*     Entries  for approximation  AN  fluxies  polinomial function     *
*     F0*Energy^(-Gamma) and employment "GZK" cutoff. (AN_SHE.for)     *
*                                                                      *
************************************************************************
*     ==================================================================
      ENTRY AN_ISU_C_en(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_en,
     #                    lnF2_SHE_en,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_C_ea(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_ea,
     #                    lnF2_SHE_ea,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_C_mn(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_mn,
     #                    lnF2_SHE_mn,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN

*     ==================================================================
      ENTRY AN_ISU_C_ma(F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,lnF_SHE_ma,
     #                    lnF2_SHE_ma,NC,NE_SHE,abs(C0),log10(E0)))
         RETURN
*     ==================================================================

      END SUBROUTINE AN_ISU_HE
