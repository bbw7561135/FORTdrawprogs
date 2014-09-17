************************************************************************
      SUBROUTINE AN_HE_ISU13
************************************************************************
*                                                                      *
*                         HE_ISU-2013 FLUX                             *
*                      in 1/(cm^2 sec sr GeV)                          *
*                Energy range is 100 GeV to 100 PeV.                   *
*                                                                      *
*     Primary cosmic ray spectra: Hilas-Gaisser.                       * 
*     Hardron interaction model: QGSJET II-03.                         *
*                                                                      *
*     Fluxes of conventional high-energy atmospheric neutrinos and     *
*     antineutrinos (total) as a function of energy and zenith angle.  *
*                                                                      *
*     The source data array F is dF/dE in units of 1/(cm^2 s sr GeV).  *
*     The energy reference points are calculated on the equidistant    *
*     over log10(E) grid and the reference points for cos(theta)       *
*     (where theta is the zenith angle) are calculated on equidistant  *
*     grid from 0 to 1 with step = 0.1.                                *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: zero,one,ten,E_cut

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

           LOGICAL(2),PARAMETER::
     #                Test   =.FALSE.,                                    Test of spline in reference points
     #                Spectra=.FALSE.,                                    Test of spline for energy spectra
     #                ZenDist=.FALSE.,                                    Test of spline for Z-A distributions
     #                Ratio  =.FALSE.                                     Test of spline for nu/antinu ratios
              INTEGER,PARAMETER::
     #                NE_test=  61,
     #                NC_test=  11,
     #                IncrE  =  10,
     #                K      =   3,
     #                NE     =  61,
     #                NE_SHE =  11,
     #                NC     =  11,
     #                Nfl    =   2,
     #                Nt     =   2
              INTEGER
     #                Ndat(Nfl)/201,202/
                 REAL,PARAMETER::
     #                E_min  = 1.0d+02,
     #                E_max  = 1.0d+08,
     #                C_min  = zero,
     #                C_max  = one
                 REAL
     #                gamm,F00,
     #                C(NC),E(NE),
     #                T(NC+1,NE+1),F(Nfl,NC,NE),
     #                lnF(Nfl,NC,NE),ClnF(Nfl,NC,NE),
     #                E_SHE(NE_SHE),logE_SHE(NE_SHE),
     #                F_SHE(Nfl,NC,NE_SHE),
     #                lnF_SHE(Nfl,NC,NE_SHE),ClnF_SHE(Nfl,NC,NE_SHE),
     #                SF(Nfl,NE_test),SF0(Nfl,NE_test),
     #                Rs(Nfl,NC_test),Factor
         CHARACTER(*),PARAMETER::
     #                dir='HEISU/',
     #                ext='.dat'
         CHARACTER*1
     #                fln(Nfl)/'e','m'/,
     #                NTn(Nt)/'n','a'/

         OPEN(Ndat00,FILE=datACN//'AN_HE_ISU13.data')
         READ(Ndat00,*) T
         DO n_NC=1,NC
           C(n_NC)=T(NC+2-n_NC,1)
      endDO
         DO n_NE=1,NE
           E(n_NE)=T(1,n_NE+1)
           DO n_NC=1,NC
             F(1,n_NC,n_NE)=T(NC+2-n_NC,n_NE+1)
        endDO
      endDO
         READ(Ndat00,*) T
         DO n_NE=1,NE
           DO n_NC=1,NC
             F(2,n_NC,n_NE)=T(NC+2-n_NC,n_NE+1)
        endDO
      endDO
         CLOSE(Ndat00)
         PRINT *,' File AN_HE_ISU13.data was read '

         lgEmin  =log10(E_min)
         lgEmax  =log10(E_max)
         steplgE =(lgEmax-lgEmin)/(NE-1)
         steplgEt=(lgEmax-lgEmin)/(NE_test-1)
         stepC   =(C_max-C_min)/(NC-1)
         stepCt  =(C_max-C_min)/(NC_test-1)

         lnF=log(F)
         DO n_fl=1,Nfl
           CALL Splie2_ED(zero,lgEmin,stepC,steplgE,lnF(n_fl,:,:),NC,NE,
     #                                              ClnF(n_fl,:,:),Test)
      endDO

         steplgE_SHE=(log10(E_cut)-lgEmax)/(NE_SHE-1)
         DO n_NE=1,NE_SHE
           logE_SHE(n_NE)=lgEmax+steplgE_SHE*(n_NE-1)
           E_SHE(n_NE)=ten**logE_SHE(n_NE)
      endDO
         DO n_NC=1,NC
           DO n_fl=1,Nfl
             CALL PowerLaw(E(NE-1),E(NE),F(n_fl,n_NC,NE-1),
     #                                         F(n_fl,n_NC,NE),gamm,F00)
             DO n_NE=1,NE_SHE
               F_SHE(n_fl,n_NC,n_NE)=AN_SHE_cut(E_SHE(n_NE))*F00*
     #                                              E_SHE(n_NE)**(-gamm)
          endDO
        endDO
      endDO
         lnF_SHE=log(F_SHE)

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN                                               Test for energy spectra
           DO n_fl=1,Nfl
             OPEN(Ndat(n_fl),FILE=outACN//dir//'FE_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             Factor=Energy**K
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),1) Energy,(Factor*F(n_fl,n_NC,n_NE),
     #                                                        n_NC=1,NC)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
             OPEN(Ndat(n_fl),FILE=outACN//dir//'SE_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=ten**lgEt
             Factor=Energy**K
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO n_fl=1,Nfl
                 SF(n_fl,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                  lnF(n_fl,:,:),ClnF(n_fl,:,:),NC,NE,Cosine,lgEt))
            endDO
          endDO
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),3) Energy,(Factor*SF(n_fl,n_NC),
     #                                                   n_NC=1,NC_test)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
        endDO
      endIF
*        ============================================================= *
*        TEST FOR ZENITH-ANGLE DISTRIBUTIONS                           *
*        ============================================================= *
         IF (ZenDist) THEN                                               Test for zenith-angle distributions
           DO n_fl=1,Nfl
             OPEN(Ndat(n_fl),FILE=outACN//dir//'FZ_'//fln(n_fl)//ext)
        endDO
           DO n_NC=1,NC
             Cosine=C(n_NC)
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),2) Cosine,(F(n_fl,n_NC,n_NE)/
     #                                  F(n_fl,NC,n_NE),n_NE=1,NE,IncrE)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
             OPEN(Ndat(n_fl),FILE=outACN//dir//'SZ_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE_test,IncrE                                       Calculating the vertical fluxes
             lgEt=lgEmin+(n_NE-1)*steplgEt
             DO n_fl=1,Nfl
               SF0(n_fl,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                     lnF(n_fl,:,:),ClnF(n_fl,:,:),NC,NE,one,lgEt))
          endDO
        endDO
           DO n_NC=1,NC_test
             Cosine=C_min+(n_NC-1)*stepCt
             DO n_NE=1,NE_test,IncrE
               lgEt=lgEmin+(n_NE-1)*steplgEt
               DO n_fl=1,Nfl
                 SF(n_fl,n_NE)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #                  lnF(n_fl,:,:),ClnF(n_fl,:,:),NC,NE,Cosine,lgEt))
            endDO
          endDO
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),4) Cosine,(SF(n_fl,n_NE)/SF0(n_fl,n_NE),
     #                                                  n_NE=1,NE,IncrE)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
        endDO
      endIF

*        ============================================================= *
*        TEST FOR NEUTRINO/ANTINEUTRINO RATIO                          *
*        ============================================================= *
         IF (Ratio) THEN                                                 Test for neutrino/antineutrino ratio
           DO n_fl=1,Nfl
             OPEN(Ndat(n_fl),FILE=outACN//dir//'FR_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),1) Energy,(F(n_fl,n_NC,n_NE)/
     #                                      F(n_fl,n_NC,n_NE),n_NC=1,NC)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
             OPEN(Ndat(n_fl),FILE=outACN//dir//'SR_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=ten**lgEt
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO n_fl=1,Nfl
                 Rs(n_fl,n_NC)=exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #           lnF(n_fl,:,:),ClnF(n_fl,:,:),NC,NE,Cosine,lgEt)-
     #           Splin2_ED(zero,lgEmin,stepC,steplgE,
     #           lnF(n_fl,:,:),ClnF(n_fl,:,:),NC,NE,Cosine,lgEt))
            endDO
          endDO
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl),3) Energy,(Rs(n_fl,n_NC),n_NC=1,NC)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl))
        endDO
      endIF
*        ============================================================= *

         RETURN

    1 FORMAT(1PE9.3,11(1PE12.4))
    2 FORMAT(1PE9.3,7(1PE12.4))
    3 FORMAT(1PE9.3,11(1PE12.4))
    4 FORMAT(1PE9.3,7(1PE12.4))


*     ==================================================================
      ENTRY AN_HE_ISU(nflavor,F0,E0,C0)
*     ==================================================================
         F0=0.5*exp(Splin2_ED(zero,lgEmin,stepC,steplgE,
     #   lnF(nflavor,:,:),ClnF(nflavor,:,:),NC,NE,abs(C0),log10(E0)))
         RETURN

************************************************************************
*                                                                      *
*     Entry for extrapolation of atmospheric neutrino fluxes with the  *
*     power function F00*Energy^(-gamm) and employment GZK cutoff      *
*     (AN_SHE.for)                                                     *
*                                                                      *
************************************************************************
*     ==================================================================
      ENTRY AN_SHE_ISU(nflavor,F0,E0,C0)
*     ==================================================================
         F0=0.5*exp(Splin2_ED(zero,lgEmax,stepC,steplgE_SHE,
     #   lnF_SHE(nflavor,:,:),ClnF_SHE(nflavor,:,:),NC,NE_SHE,abs(C0),
     #   log10(E0)))
         RETURN

  100 FORMAT((1PE10.8)$)

      END SUBROUTINE AN_HE_ISU13