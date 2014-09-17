cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     This file includes the subroutine 'CORTout' and a few support    c
c     routines for calculating the atmospheric neutrino energy and     c
c     zenith-angle distributions near the Earth surface.               c
c                                                                      c
c     -----------------------------------------------------------      c
c    |                                                           |     c
c    |                     by Vadim Naumov                       |     c
c    |                                                           |     c
c    |       Dipartimento di Fisica,  Universita di Ferrara      |     c
c    |                and INFN, Sezione di Ferrara               |     c
c    |             FORTRAN 77 version (March 2, 2002)            |     c
c    |                                                           |     c
c     -----------------------------------------------------------      c
c                                                                      c
c     The input data for the routines are calculated in cooperation    c
c     with Kostya Kuz'min, Tanya Sinegovskaya and Sergey Sinegovsky    c
c     (Physics Department of Irkutsk State University).                c
c                                                                      c
c     The physical models and mathematical methods were described      c
c     in papers [1-7].                                                 c
c                                                                      c
c     Examples of using the routines are given in file 'Test77.for'    c
c                                                                      c
c     ---------------------------------------------------------------- c
c     REFERENCES                                                       c
c     ---------------------------------------------------------------- c
c                                                                      c
c     [1] G. Fiorentini, V.A. Naumov, and F.L. Villante, "Atmospheric  c
c         neutrino flux supported by recent muon experiments", Phys.   c
c         Lett. B 510 (2001) 173 (hep-ph/0103322).                     c
c     [2] G. Fiorentini, V.A. Naumov, and F.L. Villante, "Atmospheric  c
c         neutrino flux and muon data", in: Proceedings of the 27th    c
c         International Cosmic Ray Conference, Hamburg, 2001, Vol.3,   c
c         p.1218 (hep-ph/0106014).                                     c
c     [3] V.A. Naumov, "Atmospheric muons and neutrinos", in: Proc.    c
c         of the 2nd Workshop on Methodical Aspects of Underwater/     c
c         Underice Neutrino Telescopes, Hamburg, August 15-16, 2001    c
c         [in press] (hep-ph/0201310).                                 c
c     [4] E.V. Bugaev and V.A. Naumov, "Cosmic-ray muons and neutrinos c
c         at low and intermediate energies", Yad. Fiz. 45 (1987) 1380  c
c         [Sov. J. Nucl. Phys. 45 (1987) 857].                         c
c     [5] E.V. Bugaev and V.A. Naumov, "On the interpretation of the   c
c         Kamiokande neutrino experiment", Phys. Lett. B 232 (1989)    c
c         391.                                                         c
c     [6] V.A. Naumov, T.S. Sinegovskaya, and S.I. Sinegovsky,         c
c         "The Kl3 form factors and atmospheric neutrino flavor ratio  c
c         at high energies", Il Nuovo Cimento 111A, No.2 (1998) 129    c
c         [see also e-print hep-ph/9802410].                           c
c     [7] E.V. Bugaev, V.A. Naumov, S.I. Sinegovsky, and E.S.          c
c         Zaslavskaya, "Prompt leptons in cosmic rays", Il Nuovo       c
c         Cimento 12C, No.1 (1989) 41.                                 c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

************************************************************************
      SUBROUTINE CORTout(ExpV,SA,PNmod)
************************************************************************

      INCLUDE 'Include.f'

      SAVE

      INTEGER*4 ExpV     ! Neutrino experiment number (1 to 10)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     Possible values for 'ExpV':                                       c
c     -------------------------                                        c
c      1 - BUST, Baksan Valley, Kabardino-Balkaria, Russia             c
c      2 - Frejus, Frejus tunnel, Alps, France                         c
c      3 - HPW (Harvard-Purdue-Wisconsin), Park City, Utah, USA        c
c      4 - IMB (Irvine-Michigan-Brookhaven), Cleveland, Ohio, USA      c
c      5 - Kamiokande & Super-Kamiokande, Kamioka Mine, Japan          c
c      6 - KGF (Kolar Gold Fields), India                              c
c      7 - LNGS (MACRO, LVD, ICARUS, etc.), Gran Sasso, Italy          c
c      8 - NUSEX, Mont Blanc, France                                   c
c      9 - SNO, Sudbury, Canada                                        c
c     10 - SOUDAN 2, Mine State Park, Minnesota, USA                   c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CHARACTER*3 SA    ! Level of solar activity ('Min','Mid','Max')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     Possible values for 'SA':                                        c
c     ------------------------                                         c
c     'Min' - Minimum level of solar activity [default],               c
c     'Mid' - Medium  level of solar activity,                         c
c     'max' - Maximum level of solar activity,                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER*4 PNmod   ! Prompt Neutrino (PN) production Model (0,1,2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     Possible values for 'PNmod':                                     c
c     ---------------------------                                      c
c     0 - No PN contribution (only conventional neutrinos) [default],  c
c     1 - Recombination Quark-Parton Model (RQPM),                     c
c     2 - Quark-Gluon String Model (QGSM).                             c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
cc    ---------------------------------------------------------------- c
cc    ATTENTION! The following 5 variables must be specified by user:  c
cc    ---------------------------------------------------------------- c
c                                                                      c
      CHARACTER*41 Tab/'/home/redponick/phd/fortran/Input/NuData/'/
      CHARACTER*43 Out/'/home/redponick/phd/fortran/Output/'/
      LOGICAL*2    Test/.False./          ! Default value is .False.
      LOGICAL*2    Res /.False./          ! Default value is .False.
      LOGICAL*2    B   /.True. /          ! Default value is .True.
c                                                                      c
c     1. 'Tab' is the path to the AN flux tables.                      c
c     2. 'Out' is the path to the file Test.out for a test output      c
c        (if Test=.True.).                                             c
c     3. 'Test' is the control variable for the test output            c
c        (if Test=.True., the file 'Test.out' will be created and      c
c        filled).                                                      c
c     4. 'Res' is the control variable for testing the accuracy of     c
c        spline approximation of the AN fluxes in reference points     c
c        (when Res = .True. the worst residuals are printed into       c
c        the screen for each neutrino type).                           c
c     5. 'B' is the key to switch between two interpolation methods.   c
c        For interpolating the AN fluxes through the whole energy      c
c        range, the routine may use local parabolic B-spline           c
c        (when B=.True.) or natural cubic spline (when B=.False.).     c
c        Note that B-spline works much faster (it is recommended)      c
c        while the cubic spline is a little bit more accurate          c
c        ("exact" for the reference points) and thus it may be         c
c        useful for a control.                                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      PARAMETER        ! Do not change these parameters!
     ,(
     , NE    =400,     ! Number of energies in the (Emin,Emax) range
     , KE    =50,      ! Number of energies for the low-energy range
     , NZ    =41,      ! Number of zenith angles in the (0,pi) range
     , KZ    =21,      ! Number of zenith angles for every semisphere
     , Emin  =5.0D-2,  ! Minimum neutrino energy [GeV]
     , Emax  =3.0D10,  ! Maximum neutrino energy [GeV]
     , Emid  =7.0D+1,  ! Maximum energy for the low-energy range [GeV]
     , E0    =4.0D+1,  ! Binding boundary [GeV]
     , E1    =1.0D4,   ! Bound for the PN fitting formulas
     , EPmin =1.0D3,   ! Minimum energy [GeV] for prompt neutrinos
     , StLgEP=0.1D0,   ! Step over lg(E) for prompt neutrinos
     , NEP   =76,      ! Number of energies in the (EPmin,EPmax) range
     , NZP   =KZ       ! Number of zenith angles in (0,pi) range for PN
     ,)

      PARAMETER (NC1=(KE+2)*(NZ+2))
      PARAMETER (NC2=(NE+2)*(NZ+2))
      PARAMETER (NCP=(NEP+2)*(NZP+2))

      DIMENSION D1(KE+4,NZ+4),D2(NE+4,NZ+4),DPN(NEP+4,NZP+4),
     ,          F1e1(KE,NZ),C1e1(NC1),F2e1(NE,NZ),C2e1(NC2),D2e1(NE,NZ),
     ,          F1e2(KE,NZ),C1e2(NC1),F2e2(NE,NZ),C2e2(NC2),D2e2(NE,NZ),
     ,          F1m1(KE,NZ),C1m1(NC1),F2m1(NE,NZ),C2m1(NC2),D2m1(NE,NZ),
     ,          F1m2(KE,NZ),C1m2(NC1),F2m2(NE,NZ),C2m2(NC2),D2m2(NE,NZ),
     ,          FPNe1(NEP,NZP),CPNe1(NCP),FPNe2(NEP,NZP),CPNe2(NCP)

      COMMON /BSpl/BXmin,BYmin,BStepX,BStepY
      COMMON /QSpl/QXmin,QYmin,QStepX,QStepY
      COMMON /InOut/Enu,cosZA,Fnu_e1,Fnu_e2,Fnu_m1,Fnu_m2

      EQUIVALENCE (D1(1,1),D2(1,1),DPN(1,1))
      EQUIVALENCE (F1e1(1,1),F2e1(1,1)),(F1e2(1,1),F2e2(1,1))
      EQUIVALENCE (F1m1(1,1),F2m1(1,1)),(F1m2(1,1),F2m2(1,1))

         IF (SA.NE.'Min') THEN
         SA='Min'
         PRINT *, '                                        '
         PRINT *, '  Sorry. Data for the requested level of'
         PRINT *, '  solar activity  are not available yet.'
         PRINT *, '  Default value (SA = Min) will be used.'
         PRINT *, '                                        '
      endIF

         IF (PNmod.NE.0.AND.PNmod.NE.1.AND.PNmod.NE.2) THEN
         PNmod=0
         PRINT *,'                                                    '
         PRINT *,' Wrong PNmod. Default value (PNmod=0) will be used. '
         PRINT *,'                                                    '
      endIF

      LgEmin=LOG10(Emin)
      LgEmax=LOG10(Emid)
      StLgE=(LgEmax-LgEmin)/(KE-1)
      LnEmin=LOG(Emin)
      LnEmax=LOG(Emax)
      StLnE=(LnEmax-LnEmin)/(NE-1)

      BStepY=One/(KZ-1)
      QStepY=Two/(NZ-1)
      BYmin=-One
      QYmin=-One

      Elim1=Emin-Ten**(StLgE/Ten)
      Elim2=Emax+EXP(StLnE/Ten)

      PRINT *, '  '
      GO TO(1,2,3,4,5,6,7,8,9,10) ExpV
      STOP ' Wrong experiment number'
    1 PRINT *, ' Lab 1 (BUST, Russia)'
      OPEN (1,FILE=Tab//SA//'/Z/BNO_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/BNO_D.dat',STATUS='OLD')
      GOTO 20
    2 PRINT *, ' Lab 2 (Frejus, France)'
      OPEN (1,FILE=Tab//SA//'/Z/Fre_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/Fre_D.dat',STATUS='OLD')
      GOTO 20
    3 PRINT *, ' Lab 3 (HPW, USA)'
      OPEN (1,FILE=Tab//SA//'/Z/HPW_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/HPW_D.dat',STATUS='OLD')
      GOTO 20
    4 PRINT *, ' Lab 4 (IMB, USA)'
      OPEN (1,FILE=Tab//SA//'/Z/IMB_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/IMB_D.dat',STATUS='OLD')
      GOTO 20
    5 PRINT *, ' Lab 5 (Kamiokande/Super-Kamiokande, Japan)'
      OPEN (1,FILE=Tab//SA//'/Z/Kam_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/Kam_D.dat',STATUS='OLD')
      GOTO 20
    6 PRINT *, ' Lab 6 (KGF, India)'
      OPEN (1,FILE=Tab//SA//'/Z/KGF_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/KGF_D.dat',STATUS='OLD')
      GOTO 20
    7 PRINT *, ' Lab 7 (LNGS, Italy)'
      OPEN (1,FILE=Tab//SA//'/Z/LGS_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/LGS_D.dat',STATUS='OLD')
      GOTO 20
    8 PRINT *, ' Lab 8 (NUSEX, France)'
      OPEN (1,FILE=Tab//SA//'/Z/MtB_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/MtB_D.dat',STATUS='OLD')
      GOTO 20
    9 PRINT *, ' Lab 9 (SNO, Canada)'
      OPEN (1,FILE=Tab//SA//'/Z/SNO_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/SNO_D.dat',STATUS='OLD')
      GOTO 20
   10 PRINT *, ' Lab 10 (SOUDAN 2, USA)'
      OPEN (1,FILE=Tab//SA//'/Z/SOU_U.dat',STATUS='OLD')
      OPEN (2,FILE=Tab//SA//'/Z/SOU_D.dat',STATUS='OLD')

   20 CONTINUE          ! Open full tables for conventional AN fluxes

      OPEN (11,FILE=Tab//'CN/e1.dat',STATUS='OLD')
      OPEN (12,FILE=Tab//'CN/e2.dat',STATUS='OLD')
      OPEN (13,FILE=Tab//'CN/m1.dat',STATUS='OLD')
      OPEN (14,FILE=Tab//'CN/m2.dat',STATUS='OLD')

c---- Low-energy range (Emin to Emid) ----------------------------------

         DO k=KZ,1,-1
           DO i=1,KE
           READ (1,100) F1m1(i,k),F1m2(i,k),F1e1(i,k),F1e2(i,k)
         endDO
      endDO
      CLOSE(1)
         DO k=KZ,NZ
           DO i=1,KE
           READ (2,100) F1m1(i,k),F1m2(i,k),F1e1(i,k),F1e2(i,k)
         endDO
      endDO
      CLOSE(2)

         IF (Test) THEN
         OPEN (1,FILE=Out//'Test.out')
            DO j=1,NZ
            C=BYmin+BStepY*(j-1)
               DO i=1,KE
               WRITE(1,200) Ten**(LgEmin+StLgE*(i-1)),C,
     ,                      F1m1(i,j),F1m2(i,j),F1e1(i,j),F1e2(i,j)
            endDO
         endDO
         CLOSE(1)
      endIF

      BXmin=LgEmin
      BStepX=StLgE
      CALL BSPL2(KE,NZ,F1e1,D1,C1e1,Res)
      CALL BSPL2(KE,NZ,F1e2,D1,C1e2,Res)
      CALL BSPL2(KE,NZ,F1m1,D1,C1m1,Res)
      CALL BSPL2(KE,NZ,F1m2,D1,C1m2,Res)

c---- All-energy range (Emin to Emax) ----------------------------------

         DO i=1,NE
         READ(11,300) (F2e1(i,k),k=1,NZ)
         READ(12,300) (F2e2(i,k),k=1,NZ)
         READ(13,300) (F2m1(i,k),k=1,NZ)
         READ(14,300) (F2m2(i,k),k=1,NZ)
      endDO
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)

         IF (B) THEN                           ! Use parabolic B-spline
         BXmin=LnEmin
         BStepX=StLnE
         CALL BSPL2(NE,NZ,F2e1,D2,C2e1,Res)
         CALL BSPL2(NE,NZ,F2e2,D2,C2e2,Res)
         CALL BSPL2(NE,NZ,F2m1,D2,C2m1,Res)
         CALL BSPL2(NE,NZ,F2m2,D2,C2m2,Res)
                ELSE                           ! Use cubic spline
            DO k=1,NZ
               DO i=1,NE
               F2e1(i,k)=LOG(F2e1(i,k))
               F2e2(i,k)=LOG(F2e2(i,k))
               F2m1(i,k)=LOG(F2m1(i,k))
               F2m2(i,k)=LOG(F2m2(i,k))
            endDO
         endDO
         QXmin=LnEmin
         QStepX=StLnE
         CALL QSpl2(F2e1,D2e1,NE,NZ,Res)
         CALL QSpl2(F2e2,D2e2,NE,NZ,Res)
         CALL QSpl2(F2m1,D2m1,NE,NZ,Res)
         CALL QSpl2(F2m2,D2m2,NE,NZ,Res)
      endIF

      IF (PNmod.EQ.0) RETURN

      LgE1=LOG10(EPmin)
      LgE2=(NEP-1)*StLgEP
      EPmax=Ten**LgE2

         IF (PNmod.EQ.1) THEN             ! Tables for PN fluxes (RQPM)
         OPEN (11,FILE=Tab//'PN/RQPM_N.dat',STATUS='OLD')
         OPEN (12,FILE=Tab//'PN/RQPM_A.dat',STATUS='OLD')
      endIF
         IF (PNmod.EQ.2) THEN             ! Tables for PN fluxes (QGSM)
         OPEN (11,FILE=Tab//'PN/QGSM_N.dat',STATUS='OLD')
         OPEN (12,FILE=Tab//'PN/QGSM_A.dat',STATUS='OLD')
      endIF
         DO i=1,NEP
         READ(11,400) (FPNe1(i,k),k=1,NZP)
         READ(12,400) (FPNe2(i,k),k=1,NZP)
      endDO
      CLOSE(11)
      CLOSE(12)

      BXmin=LgE1
      BYmin=Zero
      BStepX=StLgEP

      CALL BSPL2(NEP,NZP,FPNe1,DPN,CPNe1,Res)
      CALL BSPL2(NEP,NZP,FPNe2,DPN,CPNe2,Res)

      RETURN

*     ==================================================================
      ENTRY Flux(PNmod)
*     ==================================================================

         IF (cosZA.LT.-1.0001D0.OR.cosZA.GT.1.0001D0) THEN
         PRINT *,' Cos(Zenith Angle)',cosZA
         PRINT *,' is out of range'
         STOP
      endIF
         IF (Enu.LT.Elim1.OR.Enu.GT.Elim2) THEN
         PRINT *,' Neutrino energy',Enu,' GeV'
         PRINT *,' is out of range'
         STOP
      endIF

         BYmin=-One

         IF (Enu.LE.E0) THEN
         Q=One/Enu**2
         BXmin=LgEmin
         BStepX=StLgE
         Fnu_e1= Spl2(C1e1,LOG10(Enu),CosZA,KE,NZ)*Q
         Fnu_e2=rSpl2(C1e2)*Q
         Fnu_m1=rSpl2(C1m1)*Q
         Fnu_m2=rSpl2(C1m2)*Q
                        ELSE
         Q=CutOff(Enu)/Enu**3
         LnE=LOG(Enu)
            IF (B) THEN                        ! Use parabolic B-spline
            BXmin=LnEmin
            BStepX=StLnE
            Fnu_e1= Spl2(C2e1,LnE,CosZA,NE,NZ)*Q
            Fnu_e2=rSpl2(C2e2)*Q
            Fnu_m1=rSpl2(C2m1)*Q
            Fnu_m2=rSpl2(C2m2)*Q
                   ELSE                        ! Use cubic spline
            QXmin=LnEmin
            QStepX=StLnE
            Fnu_e1=EXP(SplQ2(F2e1,D2e1,NE,NZ,LnE,CosZA))*Q
            Fnu_e2=EXP(SplQ2(F2e2,D2e2,NE,NZ,LnE,CosZA))*Q
            Fnu_m1=EXP(SplQ2(F2m1,D2m1,NE,NZ,LnE,CosZA))*Q
            Fnu_m2=EXP(SplQ2(F2m2,D2m2,NE,NZ,LnE,CosZA))*Q
         endIF
      endIF

      IF (PNmod.EQ.0) RETURN

         IF (Enu.LT.E1) THEN
         CALL Fit_LE(PNmod,Enu,F1,F2)
                        ELSE
         Q=CutOff(Enu)/Enu**3
         BXmin=LgE1
         BStepX=StLgEP
         BYmin=Zero
         F1= Spl2(CPNe1,LOG10(Enu),ABS(CosZA),NEP,NZP)*Q
         F2=rSpl2(CPNe2)*Q
      endIF

      Fnu_e1=Fnu_e1+F1
      Fnu_e2=Fnu_e2+F2
      Fnu_m1=Fnu_m1+F1
      Fnu_m2=Fnu_m2+F2

      RETURN

  100 FORMAT(4(1pE11.4))
  200 FORMAT(1pE12.4,0pF8.4,4(1pE12.4))
  300 FORMAT(41(1pE10.3))                      ! 41=NZ
  400 FORMAT(21(1pE10.3))                      ! 21=NZP

      END

************************************************************************
      SUBROUTINE Fit_LE(PNmod,Enu,F1,F2)
************************************************************************

      INCLUDE 'Include.f'

      INTEGER*4 PNmod                      ! PN production Model (1,2)

c                 0: nu+nubar        1: nu              2: nubar
      PARAMETER(F10=1.0207027D-13, F11=4.9022820D-14, F12=5.3109261D-14,
     ,          E10=7.5960543D+04, E11=7.5470209D+04, E12=7.6379462D+04,
     ,          g10=1.95264526453, g11=1.95622562256, g12=1.94950495050,
     ,          a10=0.17470624922, a11=0.17584851090, a12=0.17373967777,
     ,          F20=3.3631239D-14, F21=2.3284651D-14, F22=1.0986783D-14,
     ,          E20=7.1233662D+04, E21=6.7292334D+04, E22=7.7192928D+04,
     ,          g20=1.99904990500, g21=1.96684668467, g22=2.05121512151,
     ,          a20=0.19418425218, a21=0.18390531416, a22=0.21067494637)

      s=Step(Enu)                          ! Low-energy cutoff
      GOTO(1,2) PNmod
    1 x=E11/Enu
      y=x**g11
      F1=F11*y*x*s/(One+y)**a11            ! RQPM: nu flux
      x=E12/Enu
      y=x**g12
      F2=F12*y*x*s/(One+y)**a12            ! RQPM: nubar flux
      RETURN
    2 x=E21/Enu
      y=x**g21
      F1=F21*y*x*s/(One+y)**a21            ! QGSM: nu flux
      x=E22/Enu
      y=x**g22
      F2=F22*y*x*s/(One+y)**a22            ! QGSM: nubar flux
      RETURN

      END

************************************************************************
      FUNCTION CutOff(Enu)                 ! Model for GZK cutoff
************************************************************************

      INCLUDE 'Include.f'

      PARAMETER (Ecut=3.0D10, Ethr=3.0D1, a=pi/(2*Ecut), b=pi*Ethr/2)

         IF (Enu.GT.Ecut) THEN
         CutOff=Zero
                          ELSE
         CutOff=One/(One+ABS(TAN(a*Enu)))
      endIF
      RETURN

*     ==================================================================
      ENTRY Step(Enu)                      ! Charm production threshold
*     ==================================================================

         IF (Enu.LT.Ethr) THEN
         Step=Zero
                          ELSE
         Step=One/(One+ABS(TAN((b/Enu)**2)))
      endIF
      RETURN

      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INCLUDE 'Extra/StandardFit.for'      ! Debug mode
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
