************************************************************************
      PROGRAM qesitau
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfilof   = 100,
     #                Nfilobn  = 101,
     #                Nfiloba  = 102,
     #                Nfiloe   = 103,
     #                NP_lep   = 100,
     #                NE_nu    = 100,
     #                MinCal   = 100
                 REAL,PARAMETER::
     #                Xlow     = zero,
     #                Xupp     = one,
     #                RelErr   = 1.0d-13,
     #                f1       = c2C/(8*pi)*mm_I*G_Fermi**2*hc2,         Coefficient for section
     #                f2       = 4*pi,                                   Coefficient for neutrino flux     
     #                f3       = 1.055d+39,                              Coefficient for number of events
     #                factorf  = 2*f1*f2*f3*mm_W**2,                     (m_W is from W-boson propagator)
     #                factorb  = 8*f2*f3*1.00d-38,                       (section is multiplied by 1.00d+38)
     #                E_nu_min = 1.0d-01,                                Minimal energy given by AN spectrum
     #                E_nu_max = 1.0d+03,                                Maximal energy given by AN spectrum
     #                P_lep_min= 9.0d-02,
     #                P_lep_max= 5.0d+01,
     #                Bground  = zero!5.0d-01
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*60
     #                namfof,namfobn,namfoba,namfoe
         CHARACTER*4
     #                MAn(5)/'0.90','0.95','1.10','1.20','1.35'/
         CHARACTER*3
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/
                 REAL
     #                ValP(NP_lep),ValgP(NP_lep),
     #                ValE(NE_nu),ValgE(NE_nu),
     #                Exe(NP_lep,NE_nu),Exf(NP_lep,NE_nu),
     #                Exba(NP_lep,NE_nu),Exbn(NP_lep,NE_nu)

         COMMON     /n_MA/n_MA                                           Switch for MA_QES
         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON     /n_DM/n_DM                                           Name of the Earth density model
         COMMON     /n_hi/n_hi                                           Switch for neutrino mass hierarchy
         COMMON     /n_NT/n_NT                                           Switch for neutrino type
         COMMON     /n_fl/n_fl                                           Switch fot lepton flavor
         COMMON /n_FF_QES/n_FF_QES                                       Switch for model of nucleon form factors in QES reactions
         COMMON    /P_lep/P_lep,E_lep                                    Charged lepton momentum
         COMMON    /x_lim/x_ini,deltax                                   Limits (for neutrino energy)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON /E_nu_thr/E_nu_thr,P_cher,O_lep,P_kth                    Neutrino energy threshold, Cherenkov threshold of lepton momentum
         
         EXTERNAL fui,fuifea,fuifen,fuifma,fuifmn,fuifta,fuiftn,
     #                fuibea,fuiben,fuibma,fuibmn,fuibta,fuibtn

!         DO n_DM=1,4
         n_DM    = 2
         SELECTCASE(n_DM)
               CASE(0)
                     WRITE(*,*) ' with no oscillation '
               CASE(1)
                     WRITE(*,*) ' vacuum case '
               CASE(2)
                     WRITE(*,*) ' 2LEM '
               CASE(3)
                     WRITE(*,*) ' 7LEM '
               CASE(4)
                     WRITE(*,*) ' PREM '
      endSELECT
!         DO n_hi=1,2
         n_hi    = 1
         SELECTCASE(n_hi)
               CASE(1)
                     WRITE(*,*) ' normal hierarchy '
               CASE(2)
                     WRITE(*,*) ' inverse hierarchy '
      endSELECT

!         DO N=1,2
         N       = 1
         set=dFANom_dE(one)

         n_FF_QES= 8                                                     (Bodek,Avvakumov,Bradford&Budd form factor)
         CALL NucQESFF(one,one,one,one,one,one,one,
     #                     one,one,one,one,one,one)

         DO n_fl=1,3
!         n_fl    = 3
         SELECTCASE(n_fl)
               CASE(1)
                     WRITE(*,*) ' electron '
               CASE(2)
                     WRITE(*,*) ' muon '
               CASE(3)
                     WRITE(*,*) ' tau-lepton '
      endSELECT
!         DO n_MA=1,5
         n_MA    = 2
         WRITE(*,*) ' M_A = ',MAn(n_MA)
         
         buSM=dsQESCC_dQ2_SM_set(zero,zero,MA_QES)
         bufN=dsQESCC_dQ2_fN_set(zero,zero,MA_QES)
         buFP=dsQESCC_dQ2_FP_set(zero,zero,MA_QES)

         P_lep_ini=P_lep_min
         P_lep_fin=P_lep_max
         lgP_lep_ini= log10(P_lep_ini)
         lgP_lep_fin= log10(P_lep_fin)
         steplgP_lep= (lgP_lep_fin-lgP_lep_ini)/(NP_lep-1)
         lgE_nu_min= log10(E_nu_min)
         lgE_nu_max= log10(E_nu_max)
         steplgE_nu= (lgE_nu_max-lgE_nu_min)/(NE_nu-1)

         DO n_NP_lep=1,NP_lep
           ValgP(n_NP_lep)=lgP_lep_ini+(n_NP_lep-1)*steplgP_lep
      endDO
         DO n_NE_nu=1,NE_nu
           ValgE(n_NE_nu)=lgE_nu_min+(n_NE_nu-1)*steplgE_nu
      endDO
         ValP=ten**ValgP
         ValE=ten**ValgE

         n_NT=-1
         CALL setEds

         namfof=Out//'QESitau/f'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'a'//
     #                                                    MAn(n_MA)//ext
         OPEN(Nfilof,FILE=namfof)
         WRITE(Nfilof,104) zero,ValgE
         DO n_NP_lep=1,NP_lep
           P_lep=ValP(n_NP_lep)

           E_lep= sqrt(P_lep**2+mm_lep)
           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
!           IF (P_lep.GT.O_lep) THEN
!             E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
!             IF (P_lep.LT.P_kth) THEN
!               E_nu_fin=min(E_nu_upp,E_nu_max)
!                                 ELSE
!               E_nu_fin=E_nu_max
!          endIF
!                               ELSE
!             E_nu_ini=1.0
!             E_nu_fin=1.0
!        endIF
           E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
           IF (P_lep.LE.abs(P_kth)) THEN
             IF (P_kth.LT.zero) THEN
               E_nu_fin=E_nu_ini
                                ELSE
               E_nu_fin=min(E_nu_upp,E_nu_max)
          endIF
                                    ELSE
             E_nu_fin=E_nu_max
        endIF
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-

           DO n_NE_nu=1,NE_nu
             E_nu=ValE(n_NE_nu)
             x=P_lep/E_nu
             IF ((x.LT.x_ini).OR.(x.GT.x_fin)) THEN
               Exf(n_NP_lep,n_NE_nu)=Bground
                                                ELSE
               SELECTCASE(n_fl)
                     CASE(1)
                           Exf(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuifea(x)
                     CASE(2)
                           Exf(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuifma(x)
                     CASE(3)
                           Exf(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuifta(x)
            endSELECT
          endIF
        endDO
             Exf(n_NP_lep,:)=2*m_ini*factorf*deltax*Exf(n_NP_lep,:)
             WRITE(Nfilof,104) ValgP(n_NP_lep),log10(Exf(n_NP_lep,:))
      endDO
         CLOSE(Nfilof)

         n_NT= 1
         CALL setEds

         namfobn=Out//'QESitau/b'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'n'//
     #                                                    MAn(n_MA)//ext
         OPEN(Nfilobn,FILE=namfobn)
         WRITE(Nfilobn,104) zero,ValgE
         DO n_NP_lep=1,NP_lep
           P_lep= ValP(n_NP_lep)
           E_lep= sqrt(P_lep**2+mm_lep)

           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
!           IF (P_lep.GT.O_lep) THEN
!             E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
!             IF (P_lep.LT.P_kth) THEN
!               E_nu_fin=min(E_nu_upp,E_nu_max)
!                                 ELSE
!               E_nu_fin=E_nu_max
!          endIF
!                               ELSE
!             E_nu_ini=1.0
!             E_nu_fin=1.0
!        endIF
           E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
           IF (P_lep.LE.abs(P_kth)) THEN
             IF (P_kth.LT.zero) THEN
               E_nu_fin=E_nu_ini
                                ELSE
               E_nu_fin=min(E_nu_upp,E_nu_max)
          endIF
                                    ELSE
             E_nu_fin=E_nu_max
        endIF
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-

           DO n_NE_nu=1,NE_nu
             E_nu= ValE(n_NE_nu)
             x=P_lep/E_nu
             IF ((x.LT.x_ini).OR.(x.GT.x_fin)) THEN
               Exbn(n_NP_lep,n_NE_nu)=Bground
                                               ELSE
               SELECTCASE(n_fl)
                     CASE(1)
                           Exbn(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuiben(x)
                     CASE(2)
                           Exbn(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuibmn(x)
                     CASE(3)
                           Exbn(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuibtn(x)
            endSELECT
          endIF
        endDO
             Exbn(n_NP_lep,:)=2*m_ini*factorb*deltax*Exbn(n_NP_lep,:)
             WRITE(Nfilobn,104) ValgP(n_NP_lep),log10(Exbn(n_NP_lep,:))
      endDO
         CLOSE(Nfilobn)

         n_NT=-1
         CALL setEds

         namfoba=Out//'QESitau/b'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'a'//
     #                                                    MAn(n_MA)//ext
         OPEN(Nfiloba,FILE=namfoba)
         WRITE(Nfiloba,104) zero,ValgE
         DO n_NP_lep=1,NP_lep
           P_lep= ValP(n_NP_lep)
           E_lep= sqrt(P_lep**2+mm_lep)

           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
!           IF (P_lep.GT.O_lep) THEN
!             E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
!             IF (P_lep.LT.P_kth) THEN
!               E_nu_fin=min(E_nu_upp,E_nu_max)
!                                 ELSE
!               E_nu_fin=E_nu_max
!          endIF
!                               ELSE
!             E_nu_ini=1.0
!             E_nu_fin=1.0
!        endIF
           E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
           IF (P_lep.LE.abs(P_kth)) THEN
             IF (P_kth.LT.zero) THEN
               E_nu_fin=E_nu_ini
                                ELSE
               E_nu_fin=min(E_nu_upp,E_nu_max)
          endIF
                                    ELSE
             E_nu_fin=E_nu_max
        endIF
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-

           DO n_NE_nu=1,NE_nu
             E_nu= ValE(n_NE_nu)
             x=P_lep/E_nu
             IF ((x.LT.x_ini).OR.(x.GT.x_fin)) THEN
               Exba(n_NP_lep,n_NE_nu)=Bground
                                               ELSE
               SELECTCASE(n_fl)
                     CASE(1)
                           Exba(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuibea(x)
                     CASE(2)
                           Exba(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuibma(x)
                     CASE(3)
                           Exba(n_NP_lep,n_NE_nu)=P_lep/E_lep*fuibta(x)
            endSELECT
          endIF
        endDO
             Exba(n_NP_lep,:)=2*m_ini*factorb*deltax*Exba(n_NP_lep,:)
             WRITE(Nfiloba,104) ValgP(n_NP_lep),log10(Exba(n_NP_lep,:))
      endDO
         CLOSE(Nfiloba)

         namfoe=Out//'QESitau/e'//DM(n_DM)//hin(n_hi)//fln(n_fl)//
     #                                                    MAn(n_MA)//ext
         OPEN(Nfiloe,FILE=namfoe)
         WRITE(Nfiloe,104) zero,ValgE
         DO n_NP_lep=1,NP_lep
           Exe(n_NP_lep,:)=Exf(n_NP_lep,:)+Exbn(n_NP_lep,:)+
     #                                                  Exba(n_NP_lep,:)
           WRITE(Nfiloe,104) ValgP(n_NP_lep),log10(Exe(n_NP_lep,:))
      endDO
         CLOSE(Nfiloe)

!      endDO
      endDO
!      endDO
!      endDO
!      endDO
         
         STOP 'THE END OF PROGRAM qesitau'

   99    STOP 'ERROR WITH GeMSet. PROGRAM qesitau'
  100    STOP 'ERROR WITH GeMInt. PROGRAM qesitau'

  101 FORMAT(A,$)
  103 FORMAT(I4,$)
  104 FORMAT(101(1PE16.8))
      END PROGRAM qesitau
