************************************************************************
      PROGRAM qesemat
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
     #                P_lep_max= 5.0d+01
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*63
     #                namfof,namfobn,namfoba,namfoe
         CHARACTER*4
     #                MAn(5)/'0.90','0.95','1.10','1.20','1.35'/
         CHARACTER*3
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                nc(2)/'c','v'/,
     #                nb(3)/'l','c','r'/,
     #                nl(5)/'0','n','m','x','1'/,
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/
                 REAL
     #                ValP(NP_lep),
     #                Intf(NP_lep),Intbn(NP_lep),Intba(NP_lep)

         COMMON     /n_TT/n_TT                                           Switch for nuclear target type
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
         COMMON      /n_l/n_l
         COMMON      /n_b/n_b
         COMMON   /n_corv/n_corv
         
         EXTERNAL fui,fuifea,fuifen,fuifma,fuifmn,fuifta,fuiftn,
     #                fuibea,fuiben,fuibma,fuibmn,fuibta,fuibtn
         
         CALL GeMSet(fui,one,Xlow,Xupp,RelErr,MinCal,*99)

         DO n_corv=1,2
!         n_corv    = 1

         DO n_b=1,3
!         n_b    = 2
!         DO n_l=1,5
         n_l    = 3

!         DO n_DM=1,4
         n_DM    = 4
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
         DO n_hi=1,2
!         n_hi    = 1
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

!         DO n_fl=1,3
         n_fl    = 3
         SELECTCASE(n_fl)
               CASE(1)
                     WRITE(*,*) ' electron '
               CASE(2)
                     WRITE(*,*) ' muon '
               CASE(3)
                     WRITE(*,*) ' tau-lepton '
      endSELECT
         DO n_MA=1,5
!         n_MA    = 2
         WRITE(*,*) ' M_A = ',MAn(n_MA)
         
         buSM=dsQESCC_dQ2_SM_set(zero,zero,MA_QES)
         bufN=dsQESCC_dQ2_fN_set(zero,zero,MA_QES)
         buFP=dsQESCC_dQ2_FP_set(zero,zero,MA_QES)

         P_lep_ini=P_lep_min
         P_lep_fin=P_lep_max
         lgP_lep_ini= log10(P_lep_ini)
         lgP_lep_fin= log10(P_lep_fin)
         steplgP_lep= (lgP_lep_fin-lgP_lep_ini)/(NP_lep-1)

         n_NT=-1
         CALL setEds

         namfof=Out//'QESnewP/F'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'a'//
     #                 MAn(n_MA)//'_'//nl(n_l)//nb(n_b)//nc(n_corv)//ext
         OPEN(Nfilof,FILE=namfof)
         DO n_NP_lep=1,NP_lep
           P_lep= 10**(lgP_lep_ini+(n_NP_lep-1)*steplgP_lep)
           ValP(n_NP_lep)=P_lep

           E_lep= sqrt(P_lep**2+mm_lep)
           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
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

         SELECTCASE(n_fl)
               CASE(1)
                     CALL GeMInt(fuifea,Res,Xlow,Xupp,*100)
               CASE(2)
                     CALL GeMInt(fuifma,Res,Xlow,Xupp,*100)
               CASE(3)
                     CALL GeMInt(fuifta,Res,Xlow,Xupp,*100)
      endSELECT
           Ya   = 2*m_ini*P_lep/E_lep
           Intf(n_NP_lep)=factorf*deltax*Res*Ya
           WRITE(Nfilof,102) P_lep,Intf(n_NP_lep)
      endDO
         CLOSE(Nfilof)
         !CALL GeMInf

         n_NT= 1
         CALL setEds

         namfobn=Out//'QESnewP/B'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'n'//
     #                 MAn(n_MA)//'_'//nl(n_l)//nb(n_b)//nc(n_corv)//ext
         OPEN(Nfilobn,FILE=namfobn)
         DO n_NP_lep=1,NP_lep
           P_lep= ValP(n_NP_lep)
           E_lep= sqrt(P_lep**2+mm_lep)

           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
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

         SELECTCASE(n_fl)
               CASE(1)
                     CALL GeMInt(fuiben,Resn,Xlow,Xupp,*100)
               CASE(2)
                     CALL GeMInt(fuibmn,Resn,Xlow,Xupp,*100)
               CASE(3)
                     CALL GeMInt(fuibtn,Resn,Xlow,Xupp,*100)
      endSELECT
           Ya   = 2*m_ini*P_lep/E_lep
           Intbn(n_NP_lep)=factorb*deltax*Resn*Ya
           WRITE(Nfilobn,102) P_lep,Intbn(n_NP_lep)
      endDO
         CLOSE(Nfilobn)
         !CALL GeMInf

         n_NT=-1
         CALL setEds

         namfoba=Out//'QESnewP/B'//DM(n_DM)//hin(n_hi)//fln(n_fl)//'a'//
     #                 MAn(n_MA)//'_'//nl(n_l)//nb(n_b)//nc(n_corv)//ext
         OPEN(Nfiloba,FILE=namfoba)
         DO n_NP_lep=1,NP_lep
           P_lep= ValP(n_NP_lep)
           E_lep= sqrt(P_lep**2+mm_lep)

           EpP_lep= E_lep+P_lep
           EmO_lep= E_lep-O_lep
           E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
           E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
           E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
        if(P_lep.GT.abs(P_kth))then
        E_nu_fin=E_nu_max                                              !+Inf actually
        else
        if(P_kth.LT.zero)then
            E_nu_fin=E_nu_ini
        else
            E_nu_fin=min(E_nu_upp,E_nu_max)
        endif
        endif
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-

         SELECTCASE(n_fl)
               CASE(1)
                     CALL GeMInt(fuibea,Resa,Xlow,Xupp,*100)
               CASE(2)
                     CALL GeMInt(fuibma,Resa,Xlow,Xupp,*100)
               CASE(3)
                     CALL GeMInt(fuibta,Resa,Xlow,Xupp,*100)
      endSELECT
           Ya   = 2*m_ini*P_lep/E_lep
           Intba(n_NP_lep)=factorb*deltax*Resa*Ya
           WRITE(Nfiloba,102) P_lep,Intba(n_NP_lep)
      endDO
         CLOSE(Nfiloba)
         !CALL GeMInf

         namfoe=Out//'QESnewP/E'//DM(n_DM)//hin(n_hi)//fln(n_fl)//
     #                 MAn(n_MA)//'_'//nl(n_l)//nb(n_b)//nc(n_corv)//ext
         OPEN(Nfiloe,FILE=namfoe)
         DO n_NP_lep=1,NP_lep
           P_lep= ValP(n_NP_lep)
           WRITE(Nfiloe,102) P_lep,Intf(n_NP_lep)+
     #                Intbn(n_NP_lep)+Intba(n_NP_lep)
      endDO
         CLOSE(Nfiloe)

      endDO
!      endDO
      endDO
!      endDO
!      endDO
!      endDO
      endDO
      endDO
         
         STOP 'THE END OF PROGRAM qesemat'

   99    STOP 'ERROR WITH GeMSet. PROGRAM qesemat'
  100    STOP 'ERROR WITH GeMInt. PROGRAM qesemat'

  101 FORMAT(A,$)
  102 FORMAT(2(1PE16.8))
  103 FORMAT(I4,$)
      END PROGRAM qesemat
