************************************************************************
      PROGRAM qessect
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfilofn  = 100,
     #                Nfilofa  = 101,
     #                Nfilobn  = 102,
     #                Nfiloba  = 103,
     #                NE_nu    = 100,
     #                MinCal   = 100
                 REAL,PARAMETER::
     #                Xlow     = zero,
     #                Xupp     = one,
     #                RelErr   = 1.0d-13,
     #                f1       = c2C/(8*pi)*mm_I*G_Fermi**2*hc2,         Coefficient for section
     #                factorf  = f1*mm_W**2*1.00d+38,                    (m_W is from W-boson propagator)
     #                factorb  = one,
     #                E_nu_min = 1.0d-02,
     #                E_nu_max = 1.0d+03
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*60
     #                namfofn,namfofa,namfobn,namfoba
         CHARACTER*4
     #                MAn(5)/'0.90','0.95','1.10','1.20','1.35'/
         CHARACTER*1
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/
                 REAL
     #                ValE(NE_nu),
     #                Intfn(NE_nu),Intfa(NE_nu),
     #                Intbn(NE_nu),Intba(NE_nu)

         COMMON     /n_MA/n_MA                                           Switch for MA_QES
         COMMON     /n_NT/n_NT                                           Switch for neutrino type
         COMMON     /n_fl/n_fl                                           Switch fot lepton flavor
         COMMON /n_FF_QES/n_FF_QES                                       Switch for model of nucleon form factors in QES reactions
         COMMON    /x_lim/x_ini,deltax                                   Limits (for Q2)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON /E_nu_thr/E_nu_thr,P_cher,O_lep,P_kth                    Neutrino energy threshold, Cherenkov threshold of lepton momentum
         COMMON     /E_nu/E_nu
         
         EXTERNAL fui,fuifea,fuifen,fuifma,fuifmn,fuifta,fuiftn,
     #                fuibea,fuiben,fuibma,fuibmn,fuibta,fuibtn
         
         CALL GeMSet(fui,one,Xlow,Xupp,RelErr,MinCal,*99)

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

         CALL setEds
         E_nu_ini=E_nu_min
         E_nu_fin=E_nu_max
         lgE_nu_ini= log10(E_nu_ini)
         lgE_nu_fin= log10(E_nu_fin)
         steplgE_nu= (lgE_nu_fin-lgE_nu_ini)/(NE_nu-1)

         n_NT= 1
         CALL setEds

         namfofn=Out//'section/F'//fln(n_fl)//'n'//MAn(n_MA)//ext
         OPEN(Nfilofn,FILE=namfofn)
         DO n_NE_nu=1,NE_nu
           E_nu= 10**(lgE_nu_ini+(n_NE_nu-1)*steplgE_nu)
           ValE(n_NE_nu)=E_nu

           IF (E_nu.LT.E_nu_thr) THEN
             Intfn(n_NE_nu)=zero
                                 ELSE
             CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
             x_ini = Q2_min
             x_fin = Q2_max
             deltax= x_fin-x_ini                                        !-

             SELECTCASE(n_fl)
                   CASE(1)
                         CALL GeMInt(fuifen,Res,Xlow,Xupp,*100)
                   CASE(2)
                         CALL GeMInt(fuifmn,Res,Xlow,Xupp,*100)
                   CASE(3)
                         CALL GeMInt(fuiftn,Res,Xlow,Xupp,*100)
          endSELECT
             Intfn(n_NE_nu)=factorf*deltax*Res
        endIF
           WRITE(Nfilofn,102) E_nu,Intfn(n_NE_nu)
      endDO
         CLOSE(Nfilofn)
         CALL GeMInf

         n_NT=-1
         CALL setEds

         namfofa=Out//'section/F'//fln(n_fl)//'a'//MAn(n_MA)//ext
         OPEN(Nfilofa,FILE=namfofa)
         DO n_NE_nu=1,NE_nu
           E_nu=ValE(n_NE_nu)

           IF (E_nu.LT.E_nu_thr) THEN
             Intfa(n_NE_nu)=zero
                                 ELSE
             CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
             x_ini = Q2_min
             x_fin = Q2_max
             deltax= x_fin-x_ini                                        !-

             SELECTCASE(n_fl)
                   CASE(1)
                         CALL GeMInt(fuifea,Res,Xlow,Xupp,*100)
                   CASE(2)
                         CALL GeMInt(fuifma,Res,Xlow,Xupp,*100)
                   CASE(3)
                         CALL GeMInt(fuifta,Res,Xlow,Xupp,*100)
          endSELECT
             Intfa(n_NE_nu)=factorf*deltax*Res
        endIF
           WRITE(Nfilofa,102) E_nu,Intfa(n_NE_nu)
      endDO
         CLOSE(Nfilofa)
         CALL GeMInf

         n_NT= 1
         CALL setEds

         namfobn=Out//'section/B'//fln(n_fl)//'n'//MAn(n_MA)//ext
         OPEN(Nfilobn,FILE=namfobn)
         DO n_NE_nu=1,NE_nu
           E_nu=ValE(n_NE_nu)

           IF (E_nu.LT.E_nu_thr) THEN
             Intbn(n_NE_nu)=zero
                                 ELSE
             CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
             x_ini = Q2_min
             x_fin = Q2_max
             deltax= x_fin-x_ini                                        !-

             SELECTCASE(n_fl)
                   CASE(1)
                         CALL GeMInt(fuiben,Res,Xlow,Xupp,*100)
                   CASE(2)
                         CALL GeMInt(fuibmn,Res,Xlow,Xupp,*100)
                   CASE(3)
                         CALL GeMInt(fuibtn,Res,Xlow,Xupp,*100)
          endSELECT
             Intbn(n_NE_nu)=factorb*deltax*Res
        endIF
           WRITE(Nfilobn,102) E_nu,Intbn(n_NE_nu)
      endDO
         CLOSE(Nfilobn)
         CALL GeMInf

         n_NT=-1
         CALL setEds

         namfoba=Out//'section/B'//fln(n_fl)//'a'//MAn(n_MA)//ext
         OPEN(Nfiloba,FILE=namfoba)
         DO n_NE_nu=1,NE_nu
           E_nu=ValE(n_NE_nu)

           IF (E_nu.LT.E_nu_thr) THEN
             Intba(n_NE_nu)=zero
                                 ELSE
             CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
             x_ini = Q2_min
             x_fin = Q2_max
             deltax= x_fin-x_ini                                        !-

             SELECTCASE(n_fl)
                   CASE(1)
                         CALL GeMInt(fuibea,Res,Xlow,Xupp,*100)
                   CASE(2)
                         CALL GeMInt(fuibma,Res,Xlow,Xupp,*100)
                   CASE(3)
                         CALL GeMInt(fuibta,Res,Xlow,Xupp,*100)
          endSELECT
             Intba(n_NE_nu)=factorb*deltax*Res
        endIF
           WRITE(Nfiloba,102) E_nu,Intba(n_NE_nu)
      endDO
         CLOSE(Nfiloba)
         CALL GeMInf

      endDO
!      endDO
         
         STOP 'THE END OF PROGRAM qessect'

   99    STOP 'ERROR WITH GeMSet. PROGRAM qessect'
  100    STOP 'ERROR WITH GeMInt. PROGRAM qessect'

  101 FORMAT(A,$)
  102 FORMAT(2(1PE16.8))
  103 FORMAT(I4,$)
      END PROGRAM qessect
