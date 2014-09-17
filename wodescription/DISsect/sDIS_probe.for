************************************************************************
      SUBROUTINE sDIS_probe
************************************************************************
*                                                                      *
*                                                                      *
*     NOTES                                                            *
*                                                                      *
*     Switch for model of DIS structure functions:                     *
*     0 - E.A. Paschos and J.Y. Yu, Ref.[];                            *
*     1 - This work, SF scheme "2-1";                                  *
*     2 - A. Bodek and U.K. Yang,   Ref.[];                            *
*     3 - This work, SF scheme "1-2";                                  *
*     4 - U.K. Yang, SF scheme "1-2";                                  *
*     5 - This work, QCD scheme;                                       *
*     6 - This work, temporary scheme.                                 *
*     Switch for model of F_L(x,Q^2) function:                         *
*     0 - F_L(x,Q^2) = 0, Ref.[];                                      *
*     1 - QCD calculation, Ref.[];                                     *
*     2 - Phenomenologycal calculation, Ref.[]. (default set).         *
*     Switch for model of R function:                                  *
*     1 - R(x,Q^2) = 0, Ref.[];                                        *
*     2 - Liang + R_a SLAC E143, Ref.[] (default set);                 *
*     4 - Whitlow et al. + smoothing, Ref.[];                          *
*     12 - R_a SLAC E143, Ref.[].                                      *
*     Switch for model of modification of R function:                  *
*     0 - no modification,                                             *
*     1 - "R_e --> R_\nu" accordin Yang's Ph.D., Ref.[].               *
*     Switch for lepton polarization type:                             *
*     0 - no polarization (default set), Ref.[];                       *
*     1 - "correct" or "uncorrect" polarization, Ref.[].               *
*     Switch for type of Q^2 fixing:                                   *
*     0 - frozen (default se),                                         *
*     1 - cutoff.                                                      *
*                                                                      *
*     REFEENCES                                                        *
*                                                                      *
*     [1]
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         CHARACTER(4) AG_PDF,R
         CHARACTER(3) AG_DIS
         CHARACTER(2) RT_DIS
         CHARACTER(1) Qc_DIS,BF_DIS

         CHARACTER(*),PARAMETER::
     #                PATH1='/home/redponick/phd/fortran/Output/DISout/'
           LOGICAL(2),PARAMETER::
     #                CC           =.TRUE.,                              CC cross sections
     #                NC           =.FALSE.,                             NC cross sections
     #                Tables_sDISCC=.TRUE.,
     #                Tables_qDISCC=.FALSE.,
     #                Tables_sDISNC=.FALSE.,
     #                Tables_qDISNC=.FALSE.
              INTEGER,PARAMETER::
     #                NE_nu        = 100,                                Number of points for neutrino energy
     #                MinCal       = 120                                 Number of integrand calls
                 REAL,PARAMETER::
     #                RelErr       = 1.000d-04,                          Relative integration accuracy [see MuL.for]
     #                m_fin_ini    = 1.400d+00,                          Initial W^DIS_cut
     #                E_nu_ini     = ((m_fin_ini+m_e)**2-mm_n)/(2*m_n),  Initial neutrino energy
     #                E_nu_fin     = 1.000d+11,                          Final   neutrino energy
     #                factor       = G_Fermi**2*hc2*1.0d+38/pi           Global factor for the cross sections

         COMMON       /n_LP/n_LP                                         Switch for lepton polarization type
         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON       /n_TT/n_TT                                         Switch for nuclear target type
         COMMON   /n_AG_DIS/n_AG_DIS                                     Switch for model of DIS structure functions
         COMMON   /n_FL_DIS/n_FL_DIS                                     Switch for model of function F_L
         COMMON   /n_RT_DIS/n_RT_DIS                                     Switch for model of function R
         COMMON   /n_Rc_DIS/n_Rc_DIS                                     Switch for model of modification of function R
         COMMON   /n_Qc_DIS/n_Qc_DIS                                     Switch for type of PDF limitation
         COMMON   /n_BF_DIS/n_BF_DIS                                     Switch for type of "bend-factor" for DIS structure functions
         COMMON     /PDFLIB/NGROUP,NSET                                  Parameters for PDFLIB setup
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON       /E_nu/E_nu                                         Neutrino energy
         COMMON     /Q2_DIS/Q2_DIS                                       Minimal Q^2_DIS value for PFDs
         COMMON     /A0_DIS/A0_DIS                                       Parameter in "bend-factor" of DIS structure functions
         COMMON     /B0_DIS/B0_DIS                                       Parameter in "bend-factor" of DIS structure functions
         COMMON     /C0_DIS/C0_DIS                                       Parameter in "bend-factor" of DIS structure functions
         COMMON     /W2_lim/dW2,W2_min                                   W^2 kinematic limits
         COMMON     /MulLim/Xlow(3),Xupp(3)                              MuL integration limits

         DIMENSION E_nu_out(NE_nu),sDIS(NE_nu),
     #                             qDIS(NE_nu)

         EXTERNAL FunsGeM,GeM_d2sDISCC_dQ2dnu,GeM_Q2d2sDISCC_dQ2dnu,
     #                    GeM_d2sDISNC_dQ2dnu

         n_LP    = 0                                                     Switch for lepton polarization type
         n_AG_DIS= 1                                                     Switch for model of DIS structure functions 
         n_FL_DIS= 2                                                     Switch for model of F_L(x,Q^2) function
         n_RT_DIS= 2                                                     Switch for model of R function
         n_Rc_DIS= 1                                                     Switch for model of modification of R function
         n_Qc_DIS= 0                                                     Switch for type of Q^2 fixing

         m_fin = 1.400d+00
         mm_fin= m_fin**2

         PRINT *, 'en + I', ((m_fin_ini+m_e  )**2-mm_I)/(2*m_I)
         PRINT *, 'en + n', ((m_fin_ini+m_e  )**2-mm_n)/(2*m_n)
         PRINT *, 'en + p', ((m_fin_ini+m_e  )**2-mm_p)/(2*m_p)
         PRINT *, 'mn + I', ((m_fin_ini+m_mu )**2-mm_I)/(2*m_I)
         PRINT *, 'mn + n', ((m_fin_ini+m_mu )**2-mm_n)/(2*m_n)
         PRINT *, 'mn + p', ((m_fin_ini+m_mu )**2-mm_p)/(2*m_p)
         PRINT *, 'tn + I', ((m_fin_ini+m_tau)**2-mm_I)/(2*m_I)
         PRINT *, 'tn + n', ((m_fin_ini+m_tau)**2-mm_n)/(2*m_n)
         PRINT *, 'tn + p', ((m_fin_ini+m_tau)**2-mm_p)/(2*m_p)

         SELECTCASE(n_AG_DIS)
               CASE(       0);AG_DIS='P&Y'                               E.A. Paschos and J.Y. Yu, Ref.[2]
               CASE(       1);AG_DIS='KLN'                               This work, SF scheme "2-1" 
               CASE(       2);AG_DIS='B&Y'                               A. Bodek and U.K. Yang,   Ref.[3]
               CASE(       3);AG_DIS='K12'                               This work, SF scheme "1-2"
               CASE(       4);AG_DIS='Y12'                               U.K. Yang, SF scheme "1-2"
               CASE(       5);AG_DIS='QCD'; set=FQCD_L_set(one,one)       This work, QCD scheme
               CASE(       6);AG_DIS='tmp'                               This work, temporary scheme
      endSELECT
         SELECTCASE(n_Qc_DIS)
               CASE(       0);Qc_DIS='f'                                 "frozen"
               CASE(       1);Qc_DIS='c'                                 "cutoff"
      endSELECT
         SELECTCASE(n_RT_DIS)
               CASE(       0);RT_DIS='00'                                R=0
               CASE(       1);RT_DIS='01'                                QCD
               CASE(       2);RT_DIS='02'                                Liang + R_a
               CASE(       3);RT_DIS='03'                                Whitlow et al.
               CASE(       4);RT_DIS='04'                                Whitlow et al.   + smoothing
               CASE(       5);RT_DIS='05'                                Whitlow et al.   + JLab
               CASE(       6);RT_DIS='06'                                Bartelski et al.
               CASE(       7);RT_DIS='07'                                Bartelski et al. + smoothing
               CASE(       8);RT_DIS='08'                                Bartelski et al. + JLab
               CASE(       9);RT_DIS='09'                                Alekhin
               CASE(      10);RT_DIS='10'                                Alekhin          + smoothing
               CASE(      11);RT_DIS='11'                                Alekhin          + JLab
               CASE(      12);RT_DIS='12'                                R_a E-143
               CASE(      13);RT_DIS='13'                                R_b E-143
               CASE(      14);RT_DIS='14'                                R_c E-143
               CASE(      15);RT_DIS='15'                                R1998 E-143
               CASE(      16);RT_DIS='16'                                R_a E-143        + JLab
               CASE(      17);RT_DIS='17'                                R_b E-143        + JLab
               CASE(      18);RT_DIS='18'                                R_c E-143        + JLab
               CASE(      19);RT_DIS='19'                                R1998 E-143      + JLab
      endSELECT

         Xlow=zero; Xupp=one                                             Intergation limits
         CALL MuLSet(FunsGeM,Res,RelErr,MinCal,2,*100)                   Intergation setup
         set=R_set(one,one,one)

         lgE_nu_ini=log10(E_nu_ini)
         lgE_nu_fin=log10(E_nu_fin)
         steplgE_nu=(lgE_nu_fin-lgE_nu_ini)/(NE_nu-1)
*        ------------------------------------------------------------- *
         IF (CC) THEN
*        ------------------------------------------------------------- *
           DO n_FM_DIS=3,3
             SELECTCASE(n_FM_DIS)
*                  --------------------------------------------------- *
*
*                  --------------------------------------------------- *
                   CASE(       1);AG_PDF='3.00'; NGROUP  =3; NSET=  0    MRST'04
                                  BF_DIS='1';    n_BF_DIS=1;
                                  Q2_DIS= 5.000d-01
                                  A0_DIS= 0.000d+00
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       2);AG_PDF='5.14'; NGROUP  =5; NSET= 14    GRV'98
                                  BF_DIS='1';    n_BF_DIS=1
                                  Q2_DIS= 6.100d-01
                                  A0_DIS= 0.000d+00
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       3);AG_PDF='6.01'; NGROUP  =6; NSET=  1    CTEQ6M
                                  BF_DIS='1';    n_BF_DIS=1;
                                  Q2_DIS= 8.100d-01
                                  A0_DIS= 0.000d+00
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       4);AG_PDF='6.02'; NGROUP  =6; NSET=  2    CTEQ6D
                                  BF_DIS='1';    n_BF_DIS=1;
                                  Q2_DIS= 5.000d-02
                                  A0_DIS= 0.000d+00
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       5);AG_PDF='6.5M'; NGROUP  =6; NSET=300    CTEQ6.5M
                                  BF_DIS='1';    n_BF_DIS=1;
                                  Q2_DIS= 4.600d-01
                                  A0_DIS= 0.000d+00
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00
*                  --------------------------------------------------- *
*
*                  --------------------------------------------------- *
                   CASE(       6);AG_PDF='3.00'; NGROUP  =3; NSET=  0    MRST'04
                                  BF_DIS='2';    n_BF_DIS=2
                                  Q2_DIS= 3.100d-01
                                  A0_DIS= 1.290d-01
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       7);AG_PDF='5.14'; NGROUP  =5; NSET= 14    GRV'98
                                  BF_DIS='2';    n_BF_DIS=2
                                  Q2_DIS= 2.000d-01
                                  A0_DIS= 8.500d-02
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00

                   CASE(       8);AG_PDF='6.5M'; NGROUP  =6; NSET=300    CTEQ6.5M
                                  BF_DIS='2';    n_BF_DIS=2
                                  Q2_DIS= 4.700d-01
                                  A0_DIS= 8.500d-02
                                  B0_DIS= 0.000d+00
                                  C0_DIS= 0.000d+00
*                  --------------------------------------------------- *
*
*                  --------------------------------------------------- *
                   CASE(       9);AG_PDF='5.14'; NGROUP  =5; NSET= 14    GRV'98
                                  BF_DIS='4';    n_BF_DIS=4;
                                  Q2_DIS= 1.700d-01
                                  A0_DIS= 4.030d-01
                                  B0_DIS= 3.500d-01
                                  C0_DIS= 1.407d+00

                   CASE(      10);AG_PDF='6.5M'; NGROUP  =6; NSET=300    CTEQ6.5M
                                  BF_DIS='4';    n_BF_DIS=4;
                                  Q2_DIS= 4.600d-01
                                  A0_DIS= 1.330d+00
                                  B0_DIS= 2.147d+00
                                  C0_DIS= 3.460d-01
*                  --------------------------------------------------- *
          endSELECT
             PRINT 202, Q2_DIS, A0_DIS, B0_DIS, C0_DIS
             CALL NuclStructFuns(one,one,one,one,one,one,one,one,one)
             DO n=5,5
               SELECTCASE( n)
                     CASE( 1);PRINT  1; R='en_n'; n_NT=+1; n_TT=2
                              m_ini= m_n; m_lep= m_e                     en + n --> e^-   + X
                     CASE( 2);PRINT  2; R='en_p'; n_NT=+1; n_TT=1
                              m_ini= m_p; m_lep= m_e                     en + p --> e^-   + X
                     CASE( 3);PRINT  3; R='ea_n'; n_NT=-1; n_TT=2
                              m_ini= m_n; m_lep= m_e                     ea + n --> e^+   + X
                     CASE( 4);PRINT  4; R='ea_p'; n_NT=-1; n_TT=1
                              m_ini= m_p; m_lep= m_e                     ea + p --> e^+   + X

                     CASE( 5);PRINT  5; R='mn_n'; n_NT=+1; n_TT=2
                              m_ini= m_n; m_lep= m_mu                    mn + n --> mu^-  + X
                     CASE( 6);PRINT  6; R='mn_p'; n_NT=+1; n_TT=1
                              m_ini= m_p; m_lep= m_mu                    mn + p --> mu^-  + X
                     CASE( 7);PRINT  7; R='ma_n'; n_NT=-1; n_TT=2
                              m_ini= m_n; m_lep= m_mu                    ma + n --> mu^+  + X
                     CASE( 8);PRINT  8; R='ma_p'; n_NT=-1; n_TT=1
                              m_ini= m_p; m_lep= m_mu                    ma + p --> mu^+  + X

                     CASE( 9);PRINT  9; R='tn_n'; n_NT=+1; n_TT=2
                              m_ini= m_n; m_lep= m_tau                   tn + n --> tau^- + X
                     CASE(10);PRINT 10; R='tn_p'; n_NT=+1; n_TT=1
                              m_ini= m_p; m_lep= m_tau                   tn + p --> tau^- + X
                     CASE(11);PRINT 11; R='ta_n'; n_NT=-1; n_TT=2
                              m_ini= m_n; m_lep= m_tau                   ta + n --> tau^+ + X
                     CASE(12);PRINT 12; R='ta_p'; n_NT=-1; n_TT=1
                              m_ini= m_p; m_lep= m_tau                   ta + p --> tau^+ + X
            endSELECT
               mm_ini= m_ini**2
               mm_lep= m_lep**2
                 E_nu_thr= ((m_fin+m_lep)**2-mm_ini)/(2*m_ini)
                 PRINT 204, n_FM_DIS,n,m_fin
                 DO n_E=1,NE_nu
                   E_nu= ten**(lgE_nu_ini+(n_E-1)*steplgE_nu)
                   E_nu_out(n_E)=E_nu
                   IF (E_nu.le.E_nu_thr) THEN
                     Res= zero
                                         ELSE
                     N_thr= (log10(E_nu_thr)-lgE_nu_ini)/steplgE_nu+1
                     IF (n_E.eq.N_thr .or. n_E.eq.N_thr+1) THEN
                       Res= zero
                                                           ELSE
                       CALL W2DIS_lim(E_nu,W2_min,W2_max)
                       dW2= W2_max-W2_min
                       f  = dW2*factor/(4*m_ini*E_nu**3)
                       IF (dW2.lt.zero) THEN
                         Res= zero
                                        ELSE
                         IF (Tables_sDISCC) THEN
                           CALL MuLInt(GeM_d2sDISCC_dQ2dnu,  Res,*101)
                      endIF
                         IF (Tables_qDISCC) THEN
                           CALL MuLInt(GeM_Q2d2sDISCC_dQ2dnu,Res,*101)
                      endIF
                    endIF
                  endIF
                endIF
                   IF (Tables_sDISCC) THEN
                     sDIS(n_E)= f*Res
                     PRINT 203, n_E, m_fin,E_nu,sDIS(n_E)
                endIF
                   IF (Tables_qDISCC) THEN
                     qDIS(n_E)= f*Res
                     PRINT 203, n_E, m_fin,E_nu,qDIS(n_E)
                endIF
              endDO
               IF (Tables_sDISCC) THEN
                 OPEN (Ndat01,FILE=PATH1//AG_DIS//'/'//AG_PDF//
     #           '/sDIS_'//R//'_'//Qc_DIS//
     #           '_'//BF_DIS//'_FIT_'//RT_DIS//'_1.40.dat')
                 DO n_E=1,NE_nu
                   WRITE(Ndat01,201)
     #             E_nu_out(n_E), sDIS(n_E)
              endDO
                 CLOSE(Ndat01)
            endIF
               IF (Tables_qDISCC) THEN
                 OPEN (Ndat01,FILE=PATH1//AG_DIS//'/'//AG_PDF//
     #           '/qDIS_'//R//'_'//Qc_DIS//'_'//BF_DIS//
     #           '_FIT_'//RT_DIS//'_1.40.dat')
                 DO n_E=1,NE_nu
                   WRITE(Ndat01,201)
     #             E_nu_out(n_E), qDIS(n_E)
              endDO
                 CLOSE(Ndat01)
            endIF
!               CALL RunTime(.FALSE.,.TRUE.)
          endDO
        endDO
      endIF
*        ------------------------------------------------------------- *
         IF (NC) THEN
*        ------------------------------------------------------------- *
           DO n_FM_DIS=4,4
             SELECTCASE(n_FM_DIS)
*                  --------------------------------------------------- *
*
*                  --------------------------------------------------- *
                   CASE(       1); AG_PDF='3.00'; NGROUP  =3; NSET=  0   MRST'04
                                   BF_DIS='1';    n_BF_DIS=1
                                   Q2_DIS= 5.000d-01
                                   A0_DIS= 0.000d+00
                                   B0_DIS= 0.000d+00
                                   C0_DIS= 0.000d+00

                   CASE(       2); AG_PDF='5.14'; NGROUP  =5; NSET= 14   GRV'98
                                   BF_DIS='1';    n_BF_DIS=1
                                   Q2_DIS= 6.100d-01
                                   A0_DIS= 0.000d+00
                                   B0_DIS= 0.000d+00
                                   C0_DIS= 0.000d+00

                   CASE(       3); AG_PDF='6.02'; NGROUP  =6; NSET=  2   CTEQ6D
                                   BF_DIS='1';    n_BF_DIS=1
                                   Q2_DIS= 2.300d-01
                                   A0_DIS= 0.000d+00
                                   B0_DIS= 0.000d+00
                                   C0_DIS= 0.000d+00

                   CASE(       4); AG_PDF='6.5M'; NGROUP  =6; NSET=300   CTEQ6.5M
                                   BF_DIS='1';    n_BF_DIS=1
                                   Q2_DIS= 4.600d-01
                                   A0_DIS= 0.000d+00
                                   B0_DIS= 0.000d+00
                                   C0_DIS= 0.000d+00
*                  --------------------------------------------------- *
          endSELECT
             PRINT 202, Q2_DIS, A0_DIS, B0_DIS, C0_DIS
             CALL NuclStructFuns(one,one,one,one,one,one,one,one,one)
             DO n=13,16
               SELECTCASE( n)
                     CASE(13);PRINT 13; R='nu_n'; n_NT=+1; n_TT=2        nu + n --> nu + X
                              m_ini= m_n; m_lep= zero
                     CASE(14);PRINT 14; R='nu_p'; n_NT=+1; n_TT=1        nu + p --> nu + X
                              m_ini= m_p; m_lep= zero
                     CASE(15);PRINT 15; R='an_n'; n_NT=-1; n_TT=2        an + n --> nu + X
                              m_ini= m_n; m_lep= zero
                     CASE(16);PRINT 16; R='an_p'; n_NT=-1; n_TT=1        an + p --> nu + X
                              m_ini= m_p; m_lep= zero
            endSELECT
               mm_ini= m_ini**2
               mm_lep= m_lep**2
                 E_nu_thr= ((m_fin+m_lep)**2-mm_ini)/(2*m_ini)
                 PRINT 204, n_FM_DIS,n,m_fin
                 DO n_E=1,NE_nu
                   E_nu= ten**(lgE_nu_ini+(n_E-1)*steplgE_nu)
                   E_nu_out(n_E)=E_nu
                   IF (E_nu.le.E_nu_thr) THEN
                     Res= zero
                                         ELSE
                     CALL W2DIS_lim(E_nu,W2_min,W2_max)
                     dW2= W2_max-W2_min
                     f  = dW2*factor/(4*m_ini*E_nu**3)
                     IF (dW2.lt.zero) THEN
                       Res= zero
                                      ELSE
                       CALL MuLInt(GeM_d2sDISNC_dQ2dnu,Res,*101)
                  endIF
                endIF
                   sDIS(n_E)= f*Res
                   PRINT 203, n_E, m_fin,E_nu,sDIS(n_E)
              endDO
               OPEN (Ndat01,FILE=PATH1//AG_DIS//'/'//AG_PDF//
     #         '/sDIS_'//R//'_'//Qc_DIS//'_'//BF_DIS//
     #         '_FIT_1.40.dat')
               DO n_E=1,NE_nu
                 WRITE(Ndat01,201)
     #           E_nu_out(n_E), sDIS(n_E)
            endDO
               CLOSE(Ndat01)
!               CALL RunTime(.FALSE.,.TRUE.)
          endDO
        endDO
      endIF
*        ------------------------------------------------------------- *
         CALL MulInf

         RETURN
  100    STOP 'ERROR IN MuLSet. SUBROUTINE sDIS_probe'
  101    STOP 'ERROR IN MuLInt. SUBROUTINE sDIS_probe'

  201 FORMAT(1PE11.5,1PE12.5)
  202 FORMAT(2x,'Q2_DIS =',F5.2,', A0_DIS =',F6.3,', B0_DIS =',F6.3,
     #', C0_DIS =',F6.3)
  203 FORMAT(2x,I3,2x,'m_fin =',F5.2,', E_nu =',1PE12.5,
     #', sDIS =',1PE12.5)
  204 FORMAT(2x,'n_FM_DIS =',I2,', n =',I2,', m_fin =',F5.2)

    1 FORMAT(2x,' 1',2x,'en + n --> e^-   + X')
    2 FORMAT(2x,' 2',2x,'en + p --> e^-   + X')
    3 FORMAT(2x,' 3',2x,'ea + n --> e^+   + X')
    4 FORMAT(2x,' 4',2x,'ea + p --> e^+   + X')
    5 FORMAT(2x,' 5',2x,'mn + n --> mu^-  + X')
    6 FORMAT(2x,' 6',2x,'mn + p --> mu^-  + X')
    7 FORMAT(2x,' 7',2x,'ma + n --> mu^+  + X')
    8 FORMAT(2x,' 8',2x,'ma + p --> mu^+  + X')
    9 FORMAT(2x,' 9',2x,'tn + n --> tau^- + X')
   10 FORMAT(2x,'10',2x,'tn + p --> tau^- + X')
   11 FORMAT(2x,'11',2x,'ta + n --> tau^+ + X')
   12 FORMAT(2x,'12',2x,'ta + p --> tau^+ + X')
   13 FORMAT(2x,'13',2x,'nu + n --> nu    + X')
   14 FORMAT(2x,'14',2x,'nu + p --> nu    + X')
   15 FORMAT(2x,'15',2x,'an + n --> nu    + X')
   16 FORMAT(2x,'16',2x,'an + p --> nu    + X')
      END SUBROUTINE sDIS_probe
