************************************************************************
      FUNCTION dsQESCC_dQ2_SM_set(E_nu,Q2,MA_QES)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER(2) R,T

         SAVE

         LOGICAL(2),PARAMETER::
     #              OLDSCHEME =.FALSE.
            INTEGER,PARAMETER::
     #              MinCal_GeM= 100,
     #              MinCal_MuL= 100
               REAL,PARAMETER::
     #              RelErr_GeM= 1.0d-05,
     #              RelErr_MuL= 5.0d-05,
     #              factor    = G_Fermi**2*hc2*1.0d+38/pi

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON       /n_PT/n_PT                                         Switch for lepton polarization type
         COMMON       /n_TT/n_TT                                         Switch for nuclear target type
         COMMON   /n_AP_QES/n_AP                                         Switch for model of axial form factor in QES reactions
         COMMON   /n_MS_QES/n_MS                                         Switch for value of axial mass in Sehgal's model
         COMMON   /n_GE_QES/n_GE                                         Switch for parametrization of Sachs electric form factor of neutron
         COMMON   /n_AG_QES/n_AG                                         Switch for model of QES reactions
         COMMON   /n_MC_QES/n_MC                                         Switch for correction of masses of initial and final nucleons
         COMMON   /n_SM_TMD/n_SM_TMD                                     Switch for target nucleon momentum distribution
         COMMON   /n_SM_PEF/n_SM_PEF                                     Switch for Pauli exclusion factor for outgoing nucleon
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus
         COMMON         /Q2/Q2_tmp                                       Square of mometum transfer (Q^2=-q^2)
         COMMON       /E_nu/E_nu_tmp                                     Neutrino energy
         COMMON     /MA_QES/MA_QES_tmp                                   Mass of axial-vector in QES CC reactions
         COMMON     /MV_QES/MV_QES                                       Mass of isovector in QES reactions
         COMMON     /MM_QES/MM_QES                                       Mass of monopole axial-vector in QES reactions
         COMMON     /MS_QES/MS_QES                                       Mass of scalar in QES reactions
         COMMON     /MT_QES/MT_QES                                       Mass of tensor in QES reactions
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /P_Fermi/P_Fermi                                      Fermi momentum of target nucleon
         COMMON    /T_Fermi/T_Fermi                                      Effective Fermi temperature of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy
         COMMON      /FV_SM/FV_SM                                        Normalization factor for nuclear volume
         COMMON       /xi_V/xi_V                                         Normalization of vector form factor
         COMMON       /xi_M/xi_M                                         Normalization of monopole form factor
         COMMON       /xi_S/xi_S                                         Normalization of scalar form factor
         COMMON       /xi_A/xi_A                                         Normalization of axial form factor
         COMMON       /xi_P/xi_P                                         Normalization of pseudoscalar form factor
         COMMON       /xi_T/xi_T                                         Normalization of tensor form factor
         COMMON      /phi_T/phi_T                                        Phase of tensor form factor
         COMMON      /phi_S/phi_S                                        Phase of scalar form factor
         COMMON     /MulLim/Xlow(3),Xupp(3)                              MuL integration limits

         EXTERNAL MuL_dsQESCC_dQ2_SM,GeM_FV_SM

         n_SM_TMD= 1                                                     Switch for target nucleon momentum distribution
         n_SM_PEF= 2                                                     Switch for Pauli exclusion factor for outgoing nucleon

         CALL NucQESFF(one,one,one,one,one,one,one,
     #                     one,one,one,one,one,one)

         Xlow= zero
         Xupp= one
         CALL MuLSet(MuL_dsQESCC_dQ2_SM,R1,RelErr_MuL,MinCal_MuL,2,*101)
         CALL GeMSet(GeM_FV_SM,R1,zero,one,RelErr_GeM,MinCal_GeM,  *102)

         DO n_R=1,6
           SELECTCASE(n_R)
                 CASE(  1)
                 R='en'; n_NT=+1; m_ini=m_n; m_fin=m_p; m_lep=m_e
                 CASE(  2)
                 R='ea'; n_NT=-1; m_ini=m_p; m_fin=m_n; m_lep=m_e
                 CASE(  3)
                 R='mn'; n_NT=+1; m_ini=m_n; m_fin=m_p; m_lep=m_mu
                 CASE(  4)
                 R='ma'; n_NT=-1; m_ini=m_p; m_fin=m_n; m_lep=m_mu
                 CASE(  5)
                 R='tn'; n_NT=+1; m_ini=m_n; m_fin=m_p; m_lep=m_tau
                 CASE(  6)
                 R='ta'; n_NT=-1; m_ini=m_p; m_fin=m_n; m_lep=m_tau
        endSELECT
           DO n_N=1,18
             SELECTCASE(n_N)
*                  --------------------------------------------------- *
                   CASE(  1); T='Li'; Z=  3; A=  7                       Li (LITHIUM)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0151; P_Fermi=0.1690
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0151; P_Fermi=0.1690
                endIF
*                  --------------------------------------------------- *
                   CASE(  2); T='_C'; Z=  6; A= 12                       C  (CARBON)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0256; P_Fermi=0.2210
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0257; P_Fermi=0.2210
                endIF
*                  --------------------------------------------------- *
                   CASE(  3); T='_O'; Z=  8; A= 16                       O  (OXYGEN)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0266; P_Fermi=0.2250
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0266; P_Fermi=0.2250
                endIF
*                  --------------------------------------------------- *
                   CASE(  4); T='Mg'; Z= 12; A= 24                       Mg (MAGNESIUM)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0289; P_Fermi=0.2350
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0290; P_Fermi=0.2350
                endIF
*                  --------------------------------------------------- *
                   CASE(  5); T='Si'; Z= 14; A= 28                       Si (SILICON)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0299; P_Fermi=0.2390
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0300; P_Fermi=0.2390
                endIF
*                  --------------------------------------------------- *
                   CASE(  6); T='Ca'; Z= 20; A= 40                       Ca (CALCIUM)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0329; P_Fermi=0.2510
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0330; P_Fermi=0.2510
                endIF
*                  --------------------------------------------------- *
                   CASE(  7); T='Fe'; Z= 26; A= 56                       Fe (IRON)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0361; P_Fermi=0.2630
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0330; P_Fermi=0.2510
                endIF 
*                  --------------------------------------------------- *
                   CASE(  8); T='Ni'; Z= 28; A= 59                       Ni (NICKEL)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0361; P_Fermi=0.2630
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0346; P_Fermi=0.2570
                endIF 
*                  --------------------------------------------------- *
                   CASE(  9); T='Ar'; Z= 18; A= 40                       Ar (ARGON)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0350; P_Fermi=0.2590
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0307; P_Fermi=0.2420
                endIF 
*                  --------------------------------------------------- *
                   CASE( 10); T='Sn'; Z= 50; A=119                       Sn (TIN)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0391; P_Fermi=0.2740
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0300; P_Fermi=0.2390
                endIF
*                  --------------------------------------------------- *
                   CASE( 11); T='Pb'; Z= 82; A=207                       Pb (LEAD)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0417; P_Fermi=0.2830
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0315; P_Fermi=0.2450
                endIF
*                  --------------------------------------------------- *
                   CASE( 12); T='_N'; Z=  7; A=14                        N  (NITROGEN)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0261;P_Fermi=0.2230
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0262;P_Fermi=0.2230
                endIF
*                  --------------------------------------------------- *
                   CASE( 13); T='_F'; Z=  9; A=19                        F  (FLUORINE)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0283;P_Fermi=0.2325
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0284;P_Fermi=0.2325
                endIF
*                  --------------------------------------------------- *
                   CASE( 14); T='Ne'; Z= 10; A=20                        Ne (NEON)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0278;P_Fermi=0.2300
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0278;P_Fermi=0.2300
                endIF
*                  --------------------------------------------------- *
                   CASE( 15); T='Al'; Z= 13; A=27                        Al (ALUMINUM)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0294;P_Fermi=0.2390
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0295;P_Fermi=0.2390
                endIF
*                  --------------------------------------------------- *
                   CASE( 16); T='Cl'; Z= 17; A=36                        Cl (CHLORINE)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0314;P_Fermi=0.2450
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0315;P_Fermi=0.2450
                endIF
*                  --------------------------------------------------- *
                   CASE( 17); T='Cu'; Z= 29; A=64                        Cu (COPPER)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0361;P_Fermi=0.2630
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0346;P_Fermi=0.2570
                endIF
*                  --------------------------------------------------- *
                   CASE( 18); T='Br'; Z= 35; A=80                        Br (BROMINE)
*                  --------------------------------------------------- *
                   IF (n_NT.eq.+1) THEN; E_nuBIN=0.0381;P_Fermi=0.2703
               ELSEIF (n_NT.eq.-1) THEN; E_nuBIN=0.0315;P_Fermi=0.2450
                endIF
*                  --------------------------------------------------- *
          endSELECT
             P_FeMAX= P_Fermi
             T_Fermi= 1.0d-02
             m_tar  = MassNucleus(Z,A)
             IF (n_NT.eq.+1) THEN; m_rnu= MassNucleus(Z,  A-1)           NEUTRINO
         ELSEIF (n_NT.eq.-1) THEN; m_rnu= MassNucleus(Z-1,A-1)           ANTINEUTRINO
          endIF
             mm_ini= m_ini**2
             mm_fin= m_fin**2
             mm_lep= m_lep**2
             mm_tar= m_tar**2
             mm_rnu= m_rnu**2
             CALL GeMInt(GeM_FV_SM,FV_SM,zero,one,*103)
             CALL E_nu_th_SM(E_nu_thr)
*            --------------------------------------------------------- *
             IF (T.eq.'Li') THEN                                         Li (LITHIUM)
*            --------------------------------------------------------- *
               m_tar_Li=m_tar; mm_tar_Li=mm_tar
               IF (R.eq.'en') THEN;FV_en_Li=FV_SM;E_en_Li_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Li=FV_SM;E_ea_Li_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Li=FV_SM;E_mn_Li_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Li=FV_SM;E_ma_Li_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Li=FV_SM;E_tn_Li_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Li=FV_SM;E_ta_Li_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Li  =m_rnu;   mm_rnu_nu_Li =mm_rnu
                 E_nuBIN_nu_Li=E_nuBIN; P_Fermi_nu_Li=P_Fermi
                 P_FeMAX_nu_Li=P_FeMAX; T_Fermi_nu_Li=T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Li  =m_rnu;   mm_rnu_an_Li =mm_rnu
                 E_nuBIN_an_Li=E_nuBIN; P_Fermi_an_Li=P_Fermi
                 P_FeMAX_an_Li=P_FeMAX; T_Fermi_an_Li=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'_C') THEN                                         C  (CARBON)
*            --------------------------------------------------------- *
               m_tar__C=m_tar; mm_tar__C=mm_tar
               IF (R.eq.'en') THEN;FV_en__C=FV_SM;E_en__C_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea__C=FV_SM;E_ea__C_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn__C=FV_SM;E_mn__C_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma__C=FV_SM;E_ma__C_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn__C=FV_SM;E_tn__C_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta__C=FV_SM;E_ta__C_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu__C  =m_rnu;   mm_rnu_nu__C =mm_rnu
                 E_nuBIN_nu__C=E_nuBIN; P_Fermi_nu__C=P_Fermi
                 P_FeMAX_nu__C=P_FeMAX; T_Fermi_nu__C=T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an__C  =m_rnu;   mm_rnu_an__C =mm_rnu
                 E_nuBIN_an__C=E_nuBIN; P_Fermi_an__C=P_Fermi
                 P_FeMAX_an__C=P_FeMAX; T_Fermi_an__C=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'_O') THEN                                         O  (OXYGEN)
*            --------------------------------------------------------- *
               m_tar__O=m_tar; mm_tar__O=mm_tar
               IF (R.eq.'en') THEN;FV_en__O=FV_SM;E_en__O_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea__O=FV_SM;E_ea__O_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn__O=FV_SM;E_mn__O_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma__O=FV_SM;E_ma__O_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn__O=FV_SM;E_tn__O_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta__O=FV_SM;E_ta__O_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu__O  =m_rnu;   mm_rnu_nu__O =mm_rnu
                 E_nuBIN_nu__O=E_nuBIN; P_Fermi_nu__O=P_Fermi
                 P_FeMAX_nu__O=P_FeMAX; T_Fermi_nu__O=T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an__O  =m_rnu;   mm_rnu_an__O =mm_rnu
                 E_nuBIN_an__O=E_nuBIN; P_Fermi_an__O=P_Fermi
                 P_FeMAX_an__O=P_FeMAX; T_Fermi_an__O=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Mg') THEN                                         Mg (MAGNESIUM)
*            --------------------------------------------------------- *
               m_tar_Mg=m_tar; mm_tar_Mg=mm_tar
               IF (R.eq.'en') THEN;FV_en_Mg=FV_SM;E_en_Mg_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Mg=FV_SM;E_ea_Mg_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Mg=FV_SM;E_mn_Mg_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Mg=FV_SM;E_ma_Mg_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Mg=FV_SM;E_tn_Mg_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Mg=FV_SM;E_ta_Mg_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Mg  =m_rnu;   mm_rnu_nu_Mg  =mm_rnu
                 E_nuBIN_nu_Mg=E_nuBIN; P_Fermi_nu_Mg =P_Fermi
                 P_FeMAX_nu_Mg=P_FeMAX; T_Fermi_nu_Mg =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Mg  =m_rnu;   mm_rnu_an_Mg  =mm_rnu
                 E_nuBIN_an_Mg=E_nuBIN; P_Fermi_an_Mg =P_Fermi
                 P_FeMAX_an_Mg=P_FeMAX; T_Fermi_an_Mg =T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Si') THEN                                         Si (SILICON)
*            --------------------------------------------------------- *
               m_tar_Si=m_tar; mm_tar_Si=mm_tar
               IF (R.eq.'en') THEN;FV_en_Si=FV_SM;E_en_Si_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Si=FV_SM;E_ea_Si_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Si=FV_SM;E_mn_Si_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Si=FV_SM;E_ma_Si_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Si=FV_SM;E_tn_Si_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Si=FV_SM;E_ta_Si_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Si  =m_rnu;   mm_rnu_nu_Si  =mm_rnu
                 E_nuBIN_nu_Si=E_nuBIN; P_Fermi_nu_Si =P_Fermi
                 P_FeMAX_nu_Si=P_FeMAX; T_Fermi_nu_Si =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Si  =m_rnu;   mm_rnu_an_Si  =mm_rnu
                 E_nuBIN_an_Si=E_nuBIN; P_Fermi_an_Si =P_Fermi
                 P_FeMAX_an_Si=P_FeMAX; T_Fermi_an_Si =T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Ca') THEN                                         Ca (CALCIUM)
*            --------------------------------------------------------- *
               m_tar_Ca=m_tar; mm_tar_Ca=mm_tar
               IF (R.eq.'en') THEN;FV_en_Ca=FV_SM;E_en_Ca_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Ca=FV_SM;E_ea_Ca_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Ca=FV_SM;E_mn_Ca_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Ca=FV_SM;E_ma_Ca_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Ca=FV_SM;E_tn_Ca_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Ca=FV_SM;E_ta_Ca_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Ca  =m_rnu;   mm_rnu_nu_Ca  =mm_rnu
                 E_nuBIN_nu_Ca=E_nuBIN; P_Fermi_nu_Ca =P_Fermi
                 P_FeMAX_nu_Ca=P_FeMAX; T_Fermi_nu_Ca =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Ca  =m_rnu;   mm_rnu_an_Ca  =mm_rnu
                 E_nuBIN_an_Ca=E_nuBIN; P_Fermi_an_Ca =P_Fermi
                 P_FeMAX_an_Ca=P_FeMAX; T_Fermi_an_Ca =T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Fe') THEN                                         Fe (IRON)
*            --------------------------------------------------------- *
               m_tar_Fe=m_tar; mm_tar_Fe=mm_tar
               IF (R.eq.'en') THEN;FV_en_Fe=FV_SM;E_en_Fe_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Fe=FV_SM;E_ea_Fe_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Fe=FV_SM;E_mn_Fe_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Fe=FV_SM;E_ma_Fe_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Fe=FV_SM;E_tn_Fe_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Fe=FV_SM;E_ta_Fe_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Fe  =m_rnu;   mm_rnu_nu_Fe  =mm_rnu
                 E_nuBIN_nu_Fe=E_nuBIN; P_Fermi_nu_Fe =P_Fermi
                 P_FeMAX_nu_Fe=P_FeMAX; T_Fermi_nu_Fe =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Fe  =m_rnu;   mm_rnu_an_Fe  =mm_rnu
                 E_nuBIN_an_Fe=E_nuBIN; P_Fermi_an_Fe=P_Fermi
                 P_FeMAX_an_Fe=P_FeMAX; T_Fermi_an_Fe=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Ni') THEN                                         Ni (NICKEL)
*            --------------------------------------------------------- *
               m_tar_Ni=m_tar; mm_tar_Ni=mm_tar
               IF (R.eq.'en') THEN;FV_en_Ni=FV_SM;E_en_Ni_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Ni=FV_SM;E_ea_Ni_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Ni=FV_SM;E_mn_Ni_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Ni=FV_SM;E_ma_Ni_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Ni=FV_SM;E_tn_Ni_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Ni=FV_SM;E_ta_Ni_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Ni  =m_rnu;   mm_rnu_nu_Ni  =mm_rnu
                 E_nuBIN_nu_Ni=E_nuBIN; P_Fermi_nu_Ni =P_Fermi
                 P_FeMAX_nu_Ni=P_FeMAX; T_Fermi_nu_Ni =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Ni  =m_rnu;   mm_rnu_an_Ni  =mm_rnu
                 E_nuBIN_an_Ni=E_nuBIN; P_Fermi_an_Ni =P_Fermi
                 P_FeMAX_an_Ni=P_FeMAX; T_Fermi_an_Ni =T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Ar') THEN                                         Ar (ARGON)
*            --------------------------------------------------------- *
               m_tar_Ar=m_tar; mm_tar_Ar=mm_tar
               IF (R.eq.'en') THEN;FV_en_Ar=FV_SM;E_en_Ar_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Ar=FV_SM;E_ea_Ar_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Ar=FV_SM;E_mn_Ar_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Ar=FV_SM;E_ma_Ar_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Ar=FV_SM;E_tn_Ar_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Ar=FV_SM;E_ta_Ar_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Ar  =m_rnu;   mm_rnu_nu_Ar  =mm_rnu
                 E_nuBIN_nu_Ar=E_nuBIN; P_Fermi_nu_Ar =P_Fermi
                 P_FeMAX_nu_Ar=P_FeMAX; T_Fermi_nu_Ar =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Ar  =m_rnu;   mm_rnu_an_Ar  =mm_rnu
                 E_nuBIN_an_Ar=E_nuBIN; P_Fermi_an_Ar=P_Fermi
                 P_FeMAX_an_Ar=P_FeMAX; T_Fermi_an_Ar=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Sn') THEN                                         Sn (TIN)
*            --------------------------------------------------------- *
               m_tar_Sn=m_tar; mm_tar_Sn=mm_tar
               IF (R.eq.'en') THEN;FV_en_Sn=FV_SM;E_en_Sn_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Sn=FV_SM;E_ea_Sn_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Sn=FV_SM;E_mn_Sn_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Sn=FV_SM;E_ma_Sn_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Sn=FV_SM;E_tn_Sn_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Sn=FV_SM;E_ta_Sn_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Sn=m_rnu;     mm_rnu_nu_Sn  =mm_rnu
                 E_nuBIN_nu_Sn=E_nuBIN; P_Fermi_nu_Sn =P_Fermi
                 P_FeMAX_nu_Sn=P_FeMAX; T_Fermi_nu_Sn =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Sn=m_rnu;     mm_rnu_an_Sn  =mm_rnu
                 E_nuBIN_an_Sn=E_nuBIN; P_Fermi_an_Sn=P_Fermi
                 P_FeMAX_an_Sn=P_FeMAX; T_Fermi_an_Sn=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Pb') THEN                                         Pb (LEAD)
*            --------------------------------------------------------- *
               m_tar_Pb=m_tar; mm_tar_Pb=mm_tar
               IF (R.eq.'en') THEN;FV_en_Pb=FV_SM;E_en_Pb_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Pb=FV_SM;E_ea_Pb_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Pb=FV_SM;E_mn_Pb_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Pb=FV_SM;E_ma_Pb_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Pb=FV_SM;E_tn_Pb_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Pb=FV_SM;E_ta_Pb_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Pb  =m_rnu;   mm_rnu_nu_Pb  =mm_rnu
                 E_nuBIN_nu_Pb=E_nuBIN; P_Fermi_nu_Pb =P_Fermi
                 P_FeMAX_nu_Pb=P_FeMAX; T_Fermi_nu_Pb =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Pb  =m_rnu;   mm_rnu_an_Pb  =mm_rnu
                 E_nuBIN_an_Pb=E_nuBIN; P_Fermi_an_Pb=P_Fermi
                 P_FeMAX_an_Pb=P_FeMAX; T_Fermi_an_Pb=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'_N') THEN                                         N  (NITROGEN)
*            --------------------------------------------------------- *
               m_tar__N=m_tar; mm_tar__N=mm_tar
               IF (R.eq.'en') THEN;FV_en__N=FV_SM;E_en__N_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea__N=FV_SM;E_ea__N_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn__N=FV_SM;E_mn__N_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma__N=FV_SM;E_ma__N_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn__N=FV_SM;E_tn__N_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta__N=FV_SM;E_ta__N_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu__N  =m_rnu;   mm_rnu_nu__N  =mm_rnu
                 E_nuBIN_nu__N=E_nuBIN; P_Fermi_nu__N =P_Fermi
                 P_FeMAX_nu__N=P_FeMAX; T_Fermi_nu__N =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an__N  =m_rnu;   mm_rnu_an__N  =mm_rnu
                 E_nuBIN_an__N=E_nuBIN; P_Fermi_an__N=P_Fermi
                 P_FeMAX_an__N=P_FeMAX; T_Fermi_an__N=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'_F') THEN                                         F  (FLUORINE)
*            --------------------------------------------------------- *
               m_tar__F=m_tar; mm_tar__F=mm_tar
               IF (R.eq.'en') THEN;FV_en__F=FV_SM;E_en__F_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea__F=FV_SM;E_ea__F_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn__F=FV_SM;E_mn__F_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma__F=FV_SM;E_ma__F_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn__F=FV_SM;E_tn__F_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta__F=FV_SM;E_ta__F_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu__F  =m_rnu;   mm_rnu_nu__F  =mm_rnu
                 E_nuBIN_nu__F=E_nuBIN; P_Fermi_nu__F =P_Fermi
                 P_FeMAX_nu__F=P_FeMAX; T_Fermi_nu__F =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an__F  =m_rnu;   mm_rnu_an__F  =mm_rnu
                 E_nuBIN_an__F=E_nuBIN; P_Fermi_an__F=P_Fermi
                 P_FeMAX_an__F=P_FeMAX; T_Fermi_an__F=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Ne') THEN                                         Ne (NEON)
*            --------------------------------------------------------- *
               m_tar_Ne=m_tar; mm_tar_Ne=mm_tar
               IF (R.eq.'en') THEN;FV_en_Ne=FV_SM;E_en_Ne_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Ne=FV_SM;E_ea_Ne_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Ne=FV_SM;E_mn_Ne_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Ne=FV_SM;E_ma_Ne_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Ne=FV_SM;E_tn_Ne_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Ne=FV_SM;E_ta_Ne_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Ne  =m_rnu;   mm_rnu_nu_Ne  =mm_rnu
                 E_nuBIN_nu_Ne=E_nuBIN; P_Fermi_nu_Ne =P_Fermi
                 P_FeMAX_nu_Ne=P_FeMAX; T_Fermi_nu_Ne =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Ne  =m_rnu;   mm_rnu_an_Ne  =mm_rnu
                 E_nuBIN_an_Ne=E_nuBIN; P_Fermi_an_Ne=P_Fermi
                 P_FeMAX_an_Ne=P_FeMAX; T_Fermi_an_Ne=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Al') THEN                                         Al (ALUMINIUM)
*            --------------------------------------------------------- *
               m_tar_Al=m_tar; mm_tar_Al=mm_tar
               IF (R.eq.'en') THEN;FV_en_Al=FV_SM;E_en_Al_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Al=FV_SM;E_ea_Al_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Al=FV_SM;E_mn_Al_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Al=FV_SM;E_ma_Al_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Al=FV_SM;E_tn_Al_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Al=FV_SM;E_ta_Al_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Al  =m_rnu;   mm_rnu_nu_Al  =mm_rnu
                 E_nuBIN_nu_Al=E_nuBIN; P_Fermi_nu_Al =P_Fermi
                 P_FeMAX_nu_Al=P_FeMAX; T_Fermi_nu_Al =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Al  =m_rnu;   mm_rnu_an_Al  =mm_rnu
                 E_nuBIN_an_Al=E_nuBIN; P_Fermi_an_Al=P_Fermi
                 P_FeMAX_an_Al=P_FeMAX; T_Fermi_an_Al=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Cl') THEN                                         Cl (CHLORINE)
*            --------------------------------------------------------- *
               m_tar_Cl=m_tar; mm_tar_Cl=mm_tar
               IF (R.eq.'en') THEN;FV_en_Cl=FV_SM;E_en_Cl_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Cl=FV_SM;E_ea_Cl_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Cl=FV_SM;E_mn_Cl_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Cl=FV_SM;E_ma_Cl_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Cl=FV_SM;E_tn_Cl_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Cl=FV_SM;E_ta_Cl_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Cl  =m_rnu;   mm_rnu_nu_Cl  =mm_rnu
                 E_nuBIN_nu_Cl=E_nuBIN; P_Fermi_nu_Cl =P_Fermi
                 P_FeMAX_nu_Cl=P_FeMAX; T_Fermi_nu_Cl =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Cl  =m_rnu;   mm_rnu_an_Cl  =mm_rnu
                 E_nuBIN_an_Cl=E_nuBIN; P_Fermi_an_Cl=P_Fermi
                 P_FeMAX_an_Cl=P_FeMAX; T_Fermi_an_Cl=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Cu') THEN                                         Cu (COPPER)
*            --------------------------------------------------------- *
               m_tar_Cu=m_tar; mm_tar_Cu=mm_tar
               IF (R.eq.'en') THEN;FV_en_Cu=FV_SM;E_en_Cu_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Cu=FV_SM;E_ea_Cu_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Cu=FV_SM;E_mn_Cu_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Cu=FV_SM;E_ma_Cu_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Cu=FV_SM;E_tn_Cu_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Cu=FV_SM;E_ta_Cu_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Cu  =m_rnu;   mm_rnu_nu_Cu  =mm_rnu
                 E_nuBIN_nu_Cu=E_nuBIN; P_Fermi_nu_Cu =P_Fermi
                 P_FeMAX_nu_Cu=P_FeMAX; T_Fermi_nu_Cu =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Cu  =m_rnu;   mm_rnu_an_Cu  =mm_rnu
                 E_nuBIN_an_Cu=E_nuBIN; P_Fermi_an_Cu=P_Fermi
                 P_FeMAX_an_Cu=P_FeMAX; T_Fermi_an_Cu=T_Fermi
            endIF
*            --------------------------------------------------------- *
         ELSEIF (T.eq.'Br') THEN                                         Br (BROMINE)
*            --------------------------------------------------------- *
               m_tar_Br=m_tar; mm_tar_Br=mm_tar
               IF (R.eq.'en') THEN;FV_en_Br=FV_SM;E_en_Br_thr=E_nu_thr
           ELSEIF (R.eq.'ea') THEN;FV_ea_Br=FV_SM;E_ea_Br_thr=E_nu_thr
           ELSEIF (R.eq.'mn') THEN;FV_mn_Br=FV_SM;E_mn_Br_thr=E_nu_thr
           ELSEIF (R.eq.'ma') THEN;FV_ma_Br=FV_SM;E_ma_Br_thr=E_nu_thr
           ELSEIF (R.eq.'tn') THEN;FV_tn_Br=FV_SM;E_tn_Br_thr=E_nu_thr
           ELSEIF (R.eq.'ta') THEN;FV_ta_Br=FV_SM;E_ta_Br_thr=E_nu_thr
            endIF
               IF (R.eq.'en' .or. R.eq.'mn' .or. R.eq.'tn') THEN
                 m_rnu_nu_Br  =m_rnu;   mm_rnu_nu_Br  =mm_rnu
                 E_nuBIN_nu_Br=E_nuBIN; P_Fermi_nu_Br =P_Fermi
                 P_FeMAX_nu_Br=P_FeMAX; T_Fermi_nu_Br =T_Fermi
           ELSEIF (R.eq.'ea' .or. R.eq.'ma' .or. R.eq.'ta') THEN
                 m_rnu_an_Br  =m_rnu;   mm_rnu_an_Br  =mm_rnu
                 E_nuBIN_an_Br=E_nuBIN; P_Fermi_an_Br=P_Fermi
                 P_FeMAX_an_Br=P_FeMAX; T_Fermi_an_Br=T_Fermi
            endIF
*            --------------------------------------------------------- *
          endIF
        endDO
      endDO

         E_en_thr= ((m_p+m_e  )**2-mm_n)/(2*m_n)
         E_mn_thr= ((m_p+m_mu )**2-mm_n)/(2*m_n)
         E_tn_thr= ((m_p+m_tau)**2-mm_n)/(2*m_n)
         E_ea_thr= ((m_n+m_e  )**2-mm_p)/(2*m_p)
         E_ma_thr= ((m_n+m_mu )**2-mm_p)/(2*m_p)
         E_ta_thr= ((m_n+m_tau)**2-mm_p)/(2*m_p)

         dsQESCC_dQ2_SM_set= one

         RETURN
  101    STOP 'ERROR IN MuLSet. FUNCTION dsQESCC_dQ2_SM_set'
  102    STOP 'ERROR IN GeMSet. FUNCTION dsQESCC_dQ2_SM_set'
  103    STOP 'ERROR IN MuLInt. FUNCTION dsQESCC_dQ2_SM_set'
  104    STOP 'ERROR IN GeMInt. FUNCTION dsQESCC_dQ2_SM_set'

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_en(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_n;   mm_ini    = mm_n
         m_fin     = m_p;   mm_fin    = mm_p
         m_lep     = m_e;   mm_lep    = mm_e
         n_NT      =+1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES
         SELECTCASE(n_TT)
*              ------------------------------------------------------- *
               CASE(   0) ! FREE NEUTRON
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_en= dsQESCC_dQ2_fN_en(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE(  10) ! O (OXIGEN)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_en__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 dsQESCC_dQ2_SM_en=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   dsQESCC_dQ2_SM_en=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_nu__O;   mm_rnu =mm_rnu_nu__O
                   E_nuBIN=E_nuBIN_nu__O; P_Fermi=P_Fermi_nu__O
                   P_FeMAX=P_FeMAX_nu__O; T_Fermi=T_Fermi_nu__O
                   FV_SM  =FV_en__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   dsQESCC_dQ2_SM_en=factor*S__O
              endIF
            endIF
*              ------------------------------------------------------- *
               CASE(  23) ! BNL 1985 SCINTILLATOR                        60%(C_17H_36)-40%(C_9H_12)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_en__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                                     ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_en__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
               dsQESCC_dQ2_SM_en= S__C
*              ------------------------------------------------------- *
      endSELECT
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_mn(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_n;   mm_ini    = mm_n
         m_fin     = m_p;   mm_fin    = mm_p
         m_lep     = m_mu;  mm_lep    = mm_mu
         n_NT      =+1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES
         SELECTCASE(n_TT)
*              ------------------------------------------------------- *
               CASE(   0) ! FREE NEUTRON
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_mn= dsQESCC_dQ2_fN_mn(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE(   1) ! Al (ALUMINUM)
*              ------------------------------------------------------- *
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg= Precision
                                                       ELSE
                     m_tar  = m_tar_Mg;      mm_tar = mm_tar_Mg
                     m_rnu  = m_rnu_nu_Mg;   mm_rnu = mm_rnu_nu_Mg
                     E_nuBIN= E_nuBIN_nu_Mg; P_Fermi= P_Fermi_nu_Mg
                     P_FeMAX= P_FeMAX_nu_Mg; T_Fermi= T_Fermi_nu_Mg
                     FV_SM  = FV_mn_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg= factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_mn_Si_thr) THEN                           Si (SILICON),   Z=14
                   S_Si= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Si= Precision
                                                       ELSE
                     m_tar  = m_tar_Si;      mm_tar = mm_tar_Si
                     m_rnu  = m_rnu_nu_Si;   mm_rnu = mm_rnu_nu_Si
                     E_nuBIN= E_nuBIN_nu_Si; P_Fermi= P_Fermi_nu_Si
                     P_FeMAX= P_FeMAX_nu_Si; T_Fermi= T_Fermi_nu_Si
                     FV_SM  = FV_mn_Si
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Si,*104)
                     S_Si= factor*S_Si
                endIF
              endIF
                 S_Al= (S_Mg+S_Si)/2                                     Al (ALUMINIUM), Z=13
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn_Al_thr) THEN                           Al (ALUMINIUM), Z=13
                   S_Al= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Al= Precision
                                                       ELSE
                     m_tar  = m_tar_Al;      mm_tar = mm_tar_Al
                     m_rnu  = m_rnu_nu_Al;   mm_rnu = mm_rnu_nu_Al
                     E_nuBIN= E_nuBIN_nu_Al; P_Fermi= P_Fermi_nu_Al
                     P_FeMAX= P_FeMAX_nu_Al; T_Fermi= T_Fermi_nu_Al
                     FV_SM  = FV_mn_Al
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Al,*104)
                     S_Al= factor*S_Al
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               dsQESCC_dQ2_SM_mn= S_Al
*              ------------------------------------------------------- *
               CASE(   2) ! C_3H_8 (PROPANE)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                        ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
               S_C3H_8= (3*6*S__C)/(3*6)

               dsQESCC_dQ2_SM_mn= S_C3H_8
*              ------------------------------------------------------- *
               CASE(   3) ! C_3H_8-CF_3Br (PROPANE-BROMINE FREON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                        ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__O_thr) THEN                           O  (OXYGEN),    Z= 8
                   S__O= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__O= Precision
                                          ELSE
                     m_tar  = m_tar__O;      mm_tar = mm_tar__O
                     m_rnu  = m_rnu_nu__O;   mm_rnu = mm_rnu_nu__O
                     E_nuBIN= E_nuBIN_nu__O; P_Fermi= P_Fermi_nu__O
                     P_FeMAX= P_FeMAX_nu__O; T_Fermi= T_Fermi_nu__O
                     FV_SM  = FV_mn__O
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                     S__O= factor*S__O
                endIF
              endIF
                 IF (E_nu.le.E_mn_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg= Precision
                                                       ELSE
                     m_tar  = m_tar_Mg;      mm_tar = mm_tar_Mg
                     m_rnu  = m_rnu_nu_Mg;   mm_rnu = mm_rnu_nu_Mg
                     E_nuBIN= E_nuBIN_nu_Mg; P_Fermi= P_Fermi_nu_Mg
                     P_FeMAX= P_FeMAX_nu_Mg; T_Fermi= T_Fermi_nu_Mg
                     FV_SM  = FV_mn_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg= factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_mn_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni= Precision
                                                       ELSE
                     m_tar  = m_tar_Ni;      mm_tar = mm_tar_Ni
                     m_rnu  = m_rnu_nu_Ni;   mm_rnu = mm_rnu_nu_Ni
                     E_nuBIN= E_nuBIN_nu_Ni; P_Fermi= P_Fermi_nu_Ni
                     P_FeMAX= P_FeMAX_nu_Ni; T_Fermi= T_Fermi_nu_Ni
                     FV_SM  = FV_mn_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni= factor*S_Ni
                endIF
              endIF
                 IF (E_nu.le.E_mn_Sn_thr) THEN                           Sn (TIN),       Z=50
                   S_Sn= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Sn= Precision
                                          ELSE
                     m_tar  = m_tar_Sn;      mm_tar = mm_tar_Sn
                     m_rnu  = m_rnu_nu_Sn;   mm_rnu = mm_rnu_nu_Sn
                     E_nuBIN= E_nuBIN_nu_Sn; P_Fermi= P_Fermi_nu_Sn
                     P_FeMAX= P_FeMAX_nu_Sn; T_Fermi= T_Fermi_nu_Sn
                     FV_SM  = FV_mn_Sn
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Sn,*104)
                     S_Sn= factor*S_Sn
                endIF
              endIF
                 S__F    =(   S__O+ 3*S_Mg)/ 4                           F  (FLUORINE), Z= 9
                 S_Br    =( 7*S_Ni+14*S_Sn)/21                           Br (BROMINE),  Z=35
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__F_thr) THEN                           F  (FLUORINE), Z= 9
                   S__F= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__F= Precision
                                          ELSE
                     m_tar  = m_tar__F;      mm_tar = mm_tar__F
                     m_rnu  = m_rnu_nu__F;   mm_rnu = mm_rnu_nu__F
                     E_nuBIN= E_nuBIN_nu__F; P_Fermi= P_Fermi_nu__F
                     P_FeMAX= P_FeMAX_nu__F; T_Fermi= T_Fermi_nu__F
                     FV_SM  = FV_mn__F
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__F,*104)
                     S__F= factor*S__F
                endIF
              endIF
                 IF (E_nu.le.E_mn_Br_thr) THEN                           Br (BROMINE),  Z=35
                   S_Br= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Br= Precision
                                          ELSE
                     m_tar  = m_tar_Br;      mm_tar = mm_tar_Br
                     m_rnu  = m_rnu_nu_Br;   mm_rnu = mm_rnu_nu_Br
                     E_nuBIN= E_nuBIN_nu_Br; P_Fermi= P_Fermi_nu_Br
                     P_FeMAX= P_FeMAX_nu_Br; T_Fermi= T_Fermi_nu_Br
                     FV_SM  = FV_mn_Br
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Br,*104)
                     S_Br= factor*S_Br
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               S_C_3H_8= ( 3*6*S__C                  )/(3*6)             C_3H_8 (PROPANE)
               S_CF_3Br= (   6*S__C+3*10*S__F+45*S_Br)/(6+3*10+45)       F_3Br  (BROMINE FREON)

               dsQESCC_dQ2_SM_mn= (8.7d-01*S_C_3H_8+1.3d-01*S_CF_3Br)
*              ------------------------------------------------------- *
               CASE(   4) ! C  (CARBON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                                     ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
               dsQESCC_dQ2_SM_mn= S__C
*              ------------------------------------------------------- *
               CASE(   5) ! CF_3Br (BROMINE FREON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                                     ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__O_thr) THEN                           O  (OXYGEN),    Z= 8
                   S__O= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                    S__O= Precision
                                                        ELSE
                     m_tar  =m_tar__O;      mm_tar =mm_tar__O
                     m_rnu  =m_rnu_nu__O;   mm_rnu =mm_rnu_nu__O
                     E_nuBIN=E_nuBIN_nu__O; P_Fermi=P_Fermi_nu__O
                     P_FeMAX=P_FeMAX_nu__O; T_Fermi=T_Fermi_nu__O
                     FV_SM  =FV_mn__O
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                     S__O=factor*S__O
                endIF
              endIF
                 IF (E_nu.le.E_mn_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_nu_Mg;   mm_rnu =mm_rnu_nu_Mg
                     E_nuBIN=E_nuBIN_nu_Mg; P_Fermi=P_Fermi_nu_Mg
                     P_FeMAX=P_FeMAX_nu_Mg; T_Fermi=T_Fermi_nu_Mg
                     FV_SM  =FV_mn_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_mn_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni=Precision
                                                       ELSE
                     m_tar  =m_tar_Ni;      mm_tar =mm_tar_Ni
                     m_rnu  =m_rnu_nu_Ni;   mm_rnu =mm_rnu_nu_Ni
                     E_nuBIN=E_nuBIN_nu_Ni; P_Fermi=P_Fermi_nu_Ni
                     P_FeMAX=P_FeMAX_nu_Ni; T_Fermi=T_Fermi_nu_Ni
                     FV_SM  =FV_mn_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni=factor*S_Ni
                endIF
              endIF
                 IF (E_nu.le.E_mn_Sn_thr) THEN                           Sn (TIN),       Z=50
                   S_Sn=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Sn=Precision
                                                       ELSE
                     m_tar  =m_tar_Sn;      mm_tar =mm_tar_Sn
                     m_rnu  =m_rnu_nu_Sn;   mm_rnu =mm_rnu_nu_Sn
                     E_nuBIN=E_nuBIN_nu_Sn; P_Fermi=P_Fermi_nu_Sn
                     P_FeMAX=P_FeMAX_nu_Sn; T_Fermi=T_Fermi_nu_Sn
                     FV_SM  =FV_mn_Sn
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Sn,*104)
                     S_Sn=factor*S_Sn
                endIF
              endIF
                 S__F=(  S__O+ 3*S_Mg)/ 4                                F  (FLUORINE),  Z= 9
                 S_Br=(7*S_Ni+14*S_Sn)/21                                Br (BROMINE),   Z=35
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__F_thr) THEN                           F  (FLUORINE), Z= 9
                   S__F=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__F=Precision
                                          ELSE
                     m_tar  =m_tar__F;      mm_tar =mm_tar__F
                     m_rnu  =m_rnu_nu__F;   mm_rnu =mm_rnu_nu__F
                     E_nuBIN=E_nuBIN_nu__F; P_Fermi=P_Fermi_nu__F
                     P_FeMAX=P_FeMAX_nu__F; T_Fermi=T_Fermi_nu__F
                     FV_SM  =FV_mn__F
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__F,*104)
                     S__F=factor*S__F
                endIF
              endIF
                 IF (E_nu.le.E_mn_Br_thr) THEN                           Br (BROMINE),  Z=35
                   S_Br=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Br=Precision
                                          ELSE
                     m_tar  =m_tar_Br;      mm_tar =mm_tar_Br
                     m_rnu  =m_rnu_nu_Br;   mm_rnu =mm_rnu_nu_Br
                     E_nuBIN=E_nuBIN_nu_Br; P_Fermi=P_Fermi_nu_Br
                     P_FeMAX=P_FeMAX_nu_Br; T_Fermi=T_Fermi_nu_Br
                     FV_SM  =FV_mn_Br
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Br,*104)
                     S_Br=factor*S_Br
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               S_CF_3Br=(6*S__C+3*10*S__F+45*S_Br)/(6+3*10+45)

               dsQESCC_dQ2_SM_mn=S_CF_3Br
*              ------------------------------------------------------- *
               CASE(   6) ! D_2 (DEUTERIUM)
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_mn= dsQESCC_dQ2_FP_mn(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE( 7:8) ! Fe (IRON, STEEL)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn_Fe_thr) THEN                             Fe (IRON),      Z=26
                 S_Fe=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Fe=Precision
                                                     ELSE
                   m_tar  =m_tar_Fe;      mm_tar =mm_tar_Fe
                   m_rnu  =m_rnu_nu_Fe;   mm_rnu =mm_rnu_nu_Fe
                   E_nuBIN=E_nuBIN_nu_Fe; P_Fermi=P_Fermi_nu_Fe
                   P_FeMAX=P_FeMAX_nu_Fe; T_Fermi=T_Fermi_nu_Fe
                   FV_SM  =FV_mn_Fe
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Fe,*104)
                   S_Fe=factor*S_Fe
              endIF
            endIF
               dsQESCC_dQ2_SM_mn=S_Fe
*              ------------------------------------------------------- *
               CASE(   9) ! Ne-H_2 (64% NEOUN-HYDROGEN MIXTURE)
*              ------------------------------------------------------- *
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__O_thr) THEN                           O  (OXYGEN),    Z= 8
                   S__O= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__O= Precision
                                                       ELSE
                     m_tar  = m_tar__O;      mm_tar = mm_tar__O
                     m_rnu  = m_rnu_nu__O;   mm_rnu = mm_rnu_nu__O
                     E_nuBIN= E_nuBIN_nu__O; P_Fermi= P_Fermi_nu__O
                     P_FeMAX= P_FeMAX_nu__O; T_Fermi= T_Fermi_nu__O
                     FV_SM  = FV_mn__O
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                     S__O= factor*S__O
                endIF
              endIF
                 IF (E_nu.le.E_mn_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg= Precision
                                                       ELSE
                     m_tar  = m_tar_Mg;      mm_tar = mm_tar_Mg
                     m_rnu  = m_rnu_nu_Mg;   mm_rnu = mm_rnu_nu_Mg
                     E_nuBIN= E_nuBIN_nu_Mg; P_Fermi= P_Fermi_nu_Mg
                     P_FeMAX= P_FeMAX_nu_Mg; T_Fermi= T_Fermi_nu_Mg
                     FV_SM  = FV_mn_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg= factor*S_Mg
                endIF
              endIF
                 S_Ne= (2*S__O+2*S_Mg)/4                                 Ne (NEON),      Z=10
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn_Ne_thr) THEN                           Ne (NEON),      Z=10
                   S_Ne= Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ne= Precision
                                                       ELSE
                     m_tar  = m_tar_Ne;      mm_tar = mm_tar_Ne
                     m_rnu  = m_rnu_nu_Ne;   mm_rnu = mm_rnu_nu_Ne
                     E_nuBIN= E_nuBIN_nu_Ne; P_Fermi= P_Fermi_nu_Ne
                     P_FeMAX= P_FeMAX_nu_Ne; T_Fermi= T_Fermi_nu_Ne
                     FV_SM  = FV_mn_Ne
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ne,*104)
                     S_Ne= factor*S_Ne
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               dsQESCC_dQ2_SM_mn= 6.4d-01*S_Ne
*              ------------------------------------------------------- *
               CASE(  10) ! O (OXIGEN)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 dsQESCC_dQ2_SM_mn=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   dsQESCC_dQ2_SM_mn=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_nu__O;   mm_rnu =mm_rnu_nu__O
                   E_nuBIN=E_nuBIN_nu__O; P_Fermi=P_Fermi_nu__O
                   P_FeMAX=P_FeMAX_nu__O; T_Fermi=T_Fermi_nu__O
                   FV_SM  =FV_mn__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   dsQESCC_dQ2_SM_mn=factor*S__O
              endIF
            endIF
*              ------------------------------------------------------- *
               CASE(  11) ! Ar (ARGON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn_Ar_thr) THEN                             Ar (ARGON),     Z=18
                 S_Ar=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Ar=Precision
                                                     ELSE
                   m_tar  =m_tar_Ar;      mm_tar =mm_tar_Ar
                   m_rnu  =m_rnu_nu_Ar;   mm_rnu =mm_rnu_nu_Ar
                   E_nuBIN=E_nuBIN_nu_Ar; P_Fermi=P_Fermi_nu_Ar
                   P_FeMAX=P_FeMAX_nu_Ar; T_Fermi=T_Fermi_nu_Ar
                   FV_SM  =FV_mn_Ar
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ar,*104)
                   dsQESCC_dQ2_SM_mn= factor*S_Ar
              endIF
            endIF
*              ------------------------------------------------------- *
c              CASE(  12) ! N  (NITROGEN)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  13) ! Cl (CHLORINE)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  14) ! Si (SILICON)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  15) ! Cu (COPPER)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
               CASE(  16) ! CH_2 (MINERAL OIL)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                                     ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_nu__C;   mm_rnu =mm_rnu_nu__C
                   E_nuBIN=E_nuBIN_nu__C; P_Fermi=P_Fermi_nu__C
                   P_FeMAX=P_FeMAX_nu__C; T_Fermi=T_Fermi_nu__C
                   FV_SM  =FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               dsQESCC_dQ2_SM_mn=S__C
*              ------------------------------------------------------- *
               CASE(  17) ! NOMAD MIXTURE
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                                     ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_nu__C;   mm_rnu =mm_rnu_nu__C
                   E_nuBIN=E_nuBIN_nu__C; P_Fermi=P_Fermi_nu__C
                   P_FeMAX=P_FeMAX_nu__C; T_Fermi=T_Fermi_nu__C
                   FV_SM  =FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               IF (E_nu.le.E_mn__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 S__O=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__O=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_nu__O;   mm_rnu =mm_rnu_nu__O
                   E_nuBIN=E_nuBIN_nu__O; P_Fermi=P_Fermi_nu__O
                   P_FeMAX=P_FeMAX_nu__O; T_Fermi=T_Fermi_nu__O
                   FV_SM  =FV_mn__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   S__O=factor*S__O
              endIF
            endIF
               IF (E_nu.le.E_mn_Ar_thr) THEN                             Ar (ARGON),     Z=18
                 S_Ar=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Ar=Precision
                                                     ELSE
                   m_tar  =m_tar_Ar;      mm_tar =mm_tar_Ar
                   m_rnu  =m_rnu_nu_Ar;   mm_rnu =mm_rnu_nu_Ar
                   E_nuBIN=E_nuBIN_nu_Ar; P_Fermi=P_Fermi_nu_Ar
                   P_FeMAX=P_FeMAX_nu_Ar; T_Fermi=T_Fermi_nu_Ar
                   FV_SM  =FV_mn_Ar
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ar,*104)
                   S_Ar=factor*S_Ar
              endIF
            endIF
               S__H=0                                                    H  (HYDROGEN),  Z= 1
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn_Si_thr) THEN                           Si (SILICON),   Z=14
                   S_Si=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Si=Precision
                                                       ELSE
                     m_tar  =m_tar_Si;      mm_tar =mm_tar_Si
                     m_rnu  =m_rnu_nu_Si;   mm_rnu =mm_rnu_nu_Si
                     E_nuBIN=E_nuBIN_nu_Si; P_Fermi=P_Fermi_nu_Si
                     P_FeMAX=P_FeMAX_nu_Si; T_Fermi=T_Fermi_nu_Si
                     FV_SM  =FV_mn_Si
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Si,*104)
                     S_Si=factor*S_Si
                endIF
              endIF
                 IF (E_nu.le.E_mn_Ca_thr) THEN                           Ca (CALCIUM),   Z=20
                   S_Ca=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ca=Precision
                                                       ELSE
                     m_tar  =m_tar_Ca;      mm_tar =mm_tar_Ca
                     m_rnu  =m_rnu_nu_Ca;   mm_rnu =mm_rnu_nu_Ca
                     E_nuBIN=E_nuBIN_nu_Ca; P_Fermi=P_Fermi_nu_Ca
                     P_FeMAX=P_FeMAX_nu_Ca; T_Fermi=T_Fermi_nu_Ca
                     FV_SM  =FV_mn_Ca
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ca,*104)
                     S_Ca=factor*S_Ca
                endIF
              endIF
                 IF (E_nu.le.E_mn_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_nu_Mg;   mm_rnu =mm_rnu_nu_Mg
                     E_nuBIN=E_nuBIN_nu_Mg; P_Fermi=P_Fermi_nu_Mg
                     P_FeMAX=P_FeMAX_nu_Mg; T_Fermi=T_Fermi_nu_Mg
                     FV_SM  =FV_mn_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_mn_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni=Precision
                                                       ELSE
                     m_tar  =m_tar_Ni;      mm_tar =mm_tar_Ni
                     m_rnu  =m_rnu_nu_Ni;   mm_rnu =mm_rnu_nu_Ni
                     E_nuBIN=E_nuBIN_nu_Ni; P_Fermi=P_Fermi_nu_Ni
                     P_FeMAX=P_FeMAX_nu_Ni; T_Fermi=T_Fermi_nu_Ni
                     FV_SM  =FV_mn_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni=factor*S_Ni
                endIF
              endIF
                 S__N=(  S__C+  S__O)/2                                  N  (NITROGEN),  Z= 7
                 S_Al=(  S_Mg+  S_Si)/2                                  Al (ALUMINUM),  Z=13
                 S_Cl=(3*S_Si+3*S_Ca)/6                                  Cl (CHLORINE),  Z=17
                 S_Cu=S_Ni                                               Cu (COPPER),    Z=29
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_mn__N_thr) THEN                           N  (NITROGEN),  Z= 7
                   S__N=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__N=Precision
                                                       ELSE
                     m_tar  =m_tar__N;      mm_tar =mm_tar__N
                     m_rnu  =m_rnu_nu__N;   mm_rnu =mm_rnu_nu__N
                     E_nuBIN=E_nuBIN_nu__N; P_Fermi=P_Fermi_nu__N
                     P_FeMAX=P_FeMAX_nu__N; T_Fermi=T_Fermi_nu__N
                     FV_SM  =FV_mn__N
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__N,*104)
                     S__N=factor*S__N
                endIF
              endIF
                 IF (E_nu.le.E_mn_Al_thr) THEN                           Al (ALUMINIUM),  Z=13
                   S_Al=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Al=Precision
                                                       ELSE
                     m_tar  =m_tar_Al;      mm_tar =mm_tar_Al
                     m_rnu  =m_rnu_nu_Al;   mm_rnu =mm_rnu_nu_Al
                     E_nuBIN=E_nuBIN_nu_Al; P_Fermi=P_Fermi_nu_Al
                     P_FeMAX=P_FeMAX_nu_Al; T_Fermi=T_Fermi_nu_Al
                     FV_SM  =FV_mn_Al
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Al,*104)
                     S_Al=factor*S_Al
                endIF
              endIF
                 IF (E_nu.le.E_mn_Cl_thr) THEN                           Cl (CHLORINE),  Z=17
                   S_Cl=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Cl=Precision
                                                       ELSE
                     m_tar  =m_tar_Cl;      mm_tar =mm_tar_Cl
                     m_rnu  =m_rnu_nu_Cl;   mm_rnu =mm_rnu_nu_Cl
                     E_nuBIN=E_nuBIN_nu_Cl; P_Fermi=P_Fermi_nu_Cl
                     P_FeMAX=P_FeMAX_nu_Cl; T_Fermi=T_Fermi_nu_Cl
                     FV_SM  =FV_mn_Cl
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Cl,*104)
                     S_Cl=factor*S_Cl
                endIF
              endIF
                 IF (E_nu.le.E_mn_Cu_thr) THEN                           Cu (COPPER),    Z=29
                   S_Cu=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Cu=Precision
                                                       ELSE
                     m_tar  =m_tar_Cu;      mm_tar =mm_tar_Cu
                     m_rnu  =m_rnu_nu_Cu;   mm_rnu =mm_rnu_nu_Cu
                     E_nuBIN=E_nuBIN_nu_Cu; P_Fermi=P_Fermi_nu_Cu
                     P_FeMAX=P_FeMAX_nu_Cu; T_Fermi=T_Fermi_nu_Cu
                     FV_SM  =FV_mn_Cu
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Cu,*104)
                     S_Cu=factor*S_Cu
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               C__H= 5.14d-02; N__H= 0                                   H  (HYDROGEN),  Z= 1
               C__C=64.30d-02; N__C= 6                                   C  (CARBON),    Z= 6
               C__N= 5.92d-00; N__N= 7                                   N  (NITROGEN),  Z= 7
               C__O=22.13d-02; N__O= 8                                   O  (OXIGEN),    Z= 8
               C_Al= 1.71d-00; N_Al=14                                   Al (ALUMINUM),  Z=13
               C_Si= 0.27d-00; N_Si=14                                   Si (SILICON),   Z=14
               C_Cl= 0.30d-00; N_Cl=18                                   Cl (CHLORINE),  Z=17
               C_Ar= 0.19d-00; N_Ar=22                                   Ar (ARGON),     Z=18
               C_Cu= 0.03d-00; N_Cu=35                                   Cu (COPPER),    Z=29

               dsQESCC_dQ2_SM_mn=(C__C*N__C*S__C+
     #                              C__H*N__H*S__H+
     #                              C__O*N__O*S__O+
     #                              C__N*N__N*S__N+
     #                              C_Cl*N_Cl*S_Cl+
     #                              C_Al*N_Al*S_Al+
     #                              C_Si*N_Si*S_Si+
     #                              C_Ar*N_Ar*S_Ar+
     #                              C_Cu*N_Cu*S_Cu )/
     #                             (C__C*N__C     +
     #                              C__H*N__H     +
     #                              C__O*N__O     +
     #                              C__N*N__N     +
     #                              C_Cl*N_Cl     +
     #                              C_Al*N_Al     +
     #                              C_Si*N_Si     +
     #                              C_Ar*N_Ar     +
     #                              C_Cu*N_Cu      )
*              ------------------------------------------------------- *
               CASE(  23) ! BNL 1985 SCINTILLATOR                        60%(C_17H_36)-40%(C_9H_12)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_mn__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C= Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C= Precision
                                                     ELSE
                   m_tar  = m_tar__C;      mm_tar = mm_tar__C
                   m_rnu  = m_rnu_nu__C;   mm_rnu = mm_rnu_nu__C
                   E_nuBIN= E_nuBIN_nu__C; P_Fermi= P_Fermi_nu__C
                   P_FeMAX= P_FeMAX_nu__C; T_Fermi= T_Fermi_nu__C
                   FV_SM  = FV_mn__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C= factor*S__C
              endIF
            endIF
               dsQESCC_dQ2_SM_mn= S__C
*              ------------------------------------------------------- *
      endSELECT
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_tn(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_n;   mm_ini    = mm_n
         m_fin     = m_p;   mm_fin    = mm_p
         m_lep     = m_tau; mm_lep    = mm_tau
         n_NT      =+1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES

         dsQESCC_dQ2_SM_tn= zero
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_ea(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_p;   mm_ini    = mm_p
         m_fin     = m_n;   mm_fin    = mm_n
         m_lep     = m_e;   mm_lep    = mm_e
         n_NT      =-1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES
         SELECTCASE(n_TT)
*              ------------------------------------------------------- *
               CASE(   0) ! FREE NEUTRON
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_ea= dsQESCC_dQ2_fN_ea(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE(  10) ! O (OXIGEN)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_ea__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 dsQESCC_dQ2_SM_ea=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   dsQESCC_dQ2_SM_ea=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                   E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                   P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                   FV_SM  =FV_ea__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   dsQESCC_dQ2_SM_ea=factor*S__O
              endIF
            endIF
      endSELECT
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_ma(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_p;   mm_ini    = mm_p
         m_fin     = m_n;   mm_fin    = mm_n
         m_lep     = m_mu;  mm_lep    = mm_mu
         n_NT      =-1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES
         SELECTCASE(n_TT)
*              ------------------------------------------------------- *
               CASE(   0) ! FREE PROTON (HYDROGEN)
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_ma=dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE(   1) ! Al (ALUMINUM)
*              ------------------------------------------------------- *
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_an_Mg;   mm_rnu =mm_rnu_an_Mg
                     E_nuBIN=E_nuBIN_an_Mg; P_Fermi=P_Fermi_an_Mg
                     P_FeMAX=P_FeMAX_an_Mg; T_Fermi=T_Fermi_an_Mg
                     FV_SM  =FV_ma_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_ma_Si_thr) THEN                           Si (SILICON),   Z=14
                   S_Si=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Si=Precision
                                                       ELSE
                     m_tar  =m_tar_Si;      mm_tar =mm_tar_Si
                     m_rnu  =m_rnu_an_Si;   mm_rnu =mm_rnu_an_Si
                     E_nuBIN=E_nuBIN_an_Si; P_Fermi=P_Fermi_an_Si
                     P_FeMAX=P_FeMAX_an_Si; T_Fermi=T_Fermi_an_Si
                     FV_SM  =FV_ma_Si
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Si,*104)
                     S_Si=factor*S_Si
                endIF
              endIF
                 S_Al=(S_Mg+S_Si)/2                                      Al (ALUMINIUM), Z=13
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma_Al_thr) THEN                           Al (ALUMINIUM), Z=13
                   S_Al=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Al=Precision
                                                       ELSE
                     m_tar  =m_tar_Al;      mm_tar =mm_tar_Al
                     m_rnu  =m_rnu_an_Al;   mm_rnu =mm_rnu_an_Al
                     E_nuBIN=E_nuBIN_an_Al; P_Fermi=P_Fermi_an_Al
                     P_FeMAX=P_FeMAX_an_Al; T_Fermi=T_Fermi_an_Al
                     FV_SM  =FV_ma_Al
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Al,*104)
                     S_Al=factor*S_Al
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               dsQESCC_dQ2_SM_ma=S_Al
*              ------------------------------------------------------- *
               CASE(   2) ! C_3H_8 (PROPANE)
*              ------------------------------------------------------- *
               S__p=dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)                    FREE proton
               S__H=S__p                                                 H  (HYDROGEN),  Z= 1
               IF (E_nu.le.E_ma__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                        ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_an__C;   mm_rnu =mm_rnu_an__C
                   E_nuBIN=E_nuBIN_an__C; P_Fermi=P_Fermi_an__C
                   P_FeMAX=P_FeMAX_an__C; T_Fermi=T_Fermi_an__C
                   FV_SM  =FV_ma__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               S_C_3H_8=(3*6*S__C+8*S__H)/(3*6+8)
               dsQESCC_dQ2_SM_ma=S_C_3H_8
*              ------------------------------------------------------- *
               CASE(   3) ! C_3H_8-CF_3Br (PROPANE-FREON MIXTURE)
*              ------------------------------------------------------- *
               S__p=dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)                    FREE proton
               S__H=S__p                                                 H  (HYDROGEN),  Z= 1
               IF (E_nu.le.E_ma__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                        ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_an__C;   mm_rnu =mm_rnu_an__C
                   E_nuBIN=E_nuBIN_an__C; P_Fermi=P_Fermi_an__C
                   P_FeMAX=P_FeMAX_an__C; T_Fermi=T_Fermi_an__C
                   FV_SM  =FV_ma__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               IF (E_nu.le.E_ma__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 S__O=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__O=Precision
                                        ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                   E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                   P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                   FV_SM  =FV_ma__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   S__O=factor*S__O
              endIF
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_an_Mg;   mm_rnu =mm_rnu_an_Mg
                     E_nuBIN=E_nuBIN_an_Mg; P_Fermi=P_Fermi_an_Mg
                     P_FeMAX=P_FeMAX_an_Mg; T_Fermi=T_Fermi_an_Mg
                     FV_SM  =FV_ma_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_ma_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni=Precision
                                                       ELSE
                     m_tar  =m_tar_Ni;      mm_tar =mm_tar_Ni
                     m_rnu  =m_rnu_an_Ni;   mm_rnu =mm_rnu_an_Ni
                     E_nuBIN=E_nuBIN_an_Ni; P_Fermi=P_Fermi_an_Ni
                     P_FeMAX=P_FeMAX_an_Ni; T_Fermi=T_Fermi_an_Ni
                     FV_SM  =FV_ma_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni=factor*S_Ni
                endIF
              endIF
                 IF (E_nu.le.E_ma_Sn_thr) THEN                           Sn (TIN),       Z=50
                   S_Sn=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Sn=Precision
                                          ELSE
                     m_tar  =m_tar_Sn;      mm_tar =mm_tar_Sn
                     m_rnu  =m_rnu_an_Sn;   mm_rnu =mm_rnu_an_Sn
                     E_nuBIN=E_nuBIN_an_Sn; P_Fermi=P_Fermi_an_Sn
                     P_FeMAX=P_FeMAX_an_Sn; T_Fermi=T_Fermi_an_Sn
                     FV_SM  =FV_ma_Sn
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Sn,*104)
                     S_Sn=factor*S_Sn
                endIF
              endIF
                 S__F=(  S__O+ 3*S_Mg)/ 4                                F  (FLUORINE),  Z= 9
                 S_Br=(7*S_Ni+14*S_Sn)/21                                Br (BROMINE),   Z=35
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma__F_thr) THEN                           F  (FLUORINE),  Z= 9
                   S__F=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__F=Precision
                                                       ELSE
                     m_tar  =m_tar__F;      mm_tar =mm_tar__F
                     m_rnu  =m_rnu_an__F;   mm_rnu =mm_rnu_an__F
                     E_nuBIN=E_nuBIN_an__F; P_Fermi=P_Fermi_an__F
                     P_FeMAX=P_FeMAX_an__F; T_Fermi=T_Fermi_an__F
                     FV_SM  =FV_ma__F
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__F,*104)
                     S__F=factor*S__F
                endIF
              endIF
                 IF (E_nu.le.E_ma_Br_thr) THEN                           Br (BROMINE),   Z=35
                   S_Br=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Br=Precision
                                                       ELSE
                     m_tar  =m_tar_Br;      mm_tar =mm_tar_Br
                     m_rnu  =m_rnu_an_Br;   mm_rnu =mm_rnu_an_Br
                     E_nuBIN=E_nuBIN_an_Br; P_Fermi=P_Fermi_an_Br
                     P_FeMAX=P_FeMAX_an_Br; T_Fermi=T_Fermi_an_Br
                     FV_SM  =FV_ma_Br
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Br,*104)
                     S_Br=factor*S_Br
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               S_C_3H_8=(3*6*S__C+8*  S__H        )/(3*6+8)              C_3H_8 (PROPANE)
               S_CF_3Br=(  6*S__C+3*9*S__F+35*S_Br)/(6+3*9+35)           F_3Br  (BROMINE FREON)

               dsQESCC_dQ2_SM_ma=(0.87*S_C_3H_8+0.13*S_CF_3Br)
*              ------------------------------------------------------- *
               CASE(   4) ! C (CARBON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_ma__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                                     ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_an__C;   mm_rnu =mm_rnu_an__C
                   E_nuBIN=E_nuBIN_an__C; P_Fermi=P_Fermi_an__C
                   P_FeMAX=P_FeMAX_an__C; T_Fermi=T_Fermi_an__C
                   FV_SM  =FV_ma__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               dsQESCC_dQ2_SM_ma=S__C
*              ------------------------------------------------------- *
               CASE(   5) ! CF_3Br (BROMINE FREON)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_ma__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                                     ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_an__C;   mm_rnu =mm_rnu_an__C
                   E_nuBIN=E_nuBIN_an__C; P_Fermi=P_Fermi_an__C
                   P_FeMAX=P_FeMAX_an__C; T_Fermi=T_Fermi_an__C
                   FV_SM  =FV_ma__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma__O_thr) THEN                           O  (OXYGEN),    Z= 8
                   S__O=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__O=Precision
                                                       ELSE
                     m_tar  =m_tar__O;      mm_tar =mm_tar__O
                     m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                     E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                     P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                     FV_SM  =FV_ma__O
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                     S__O=factor*S__O
                endIF
              endIF
                 IF (E_nu.le.E_ma_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_an_Mg;   mm_rnu =mm_rnu_an_Mg
                     E_nuBIN=E_nuBIN_an_Mg; P_Fermi=P_Fermi_an_Mg
                     P_FeMAX=P_FeMAX_an_Mg; T_Fermi=T_Fermi_an_Mg
                     FV_SM  =FV_ma_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_ma_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni=Precision
                                                       ELSE
                     m_tar  =m_tar_Ni;      mm_tar =mm_tar_Ni
                     m_rnu  =m_rnu_an_Ni;   mm_rnu =mm_rnu_an_Ni
                     E_nuBIN=E_nuBIN_an_Ni; P_Fermi=P_Fermi_an_Ni
                     P_FeMAX=P_FeMAX_an_Ni; T_Fermi=T_Fermi_an_Ni
                     FV_SM  =FV_ma_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni=factor*S_Ni
                endIF
              endIF
                 IF (E_nu.le.E_ma_Sn_thr) THEN                           Sn (TIN),       Z=50
                   S_Sn=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Sn=Precision
                                                       ELSE
                     m_tar  =m_tar_Sn;      mm_tar =mm_tar_Sn
                     m_rnu  =m_rnu_an_Sn;   mm_rnu =mm_rnu_an_Sn
                     E_nuBIN=E_nuBIN_an_Sn; P_Fermi=P_Fermi_an_Sn
                     P_FeMAX=P_FeMAX_an_Sn; T_Fermi=T_Fermi_an_Sn
                     FV_SM  =FV_ma_Sn
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Sn,*104)
                     S_Sn=factor*S_Sn
                endIF
              endIF
                 S__F=(  S__O+ 3*S_Mg)/ 4                                F  (FLUORINE),  Z= 9
                 S_Br=(7*S_Ni+14*S_Sn)/21                                Br (BROMINE),   Z=35
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma__F_thr) THEN                           F  (FLUORINE),  Z= 9
                   S__F=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__F=Precision
                                                       ELSE
                     m_tar  =m_tar__F;      mm_tar =mm_tar__F
                     m_rnu  =m_rnu_an__F;   mm_rnu =mm_rnu_an__F
                     E_nuBIN=E_nuBIN_an__F; P_Fermi=P_Fermi_an__F
                     P_FeMAX=P_FeMAX_an__F; T_Fermi=T_Fermi_an__F
                     FV_SM  =FV_ma__F
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__F,*104)
                     S__F=factor*S__F
                endIF
              endIF
                 IF (E_nu.le.E_ma_Br_thr) THEN                           Br (BROMINE),   Z=35
                   S_Br=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Br=Precision
                                                       ELSE
                     m_tar  =m_tar_Br;      mm_tar =mm_tar_Br
                     m_rnu  =m_rnu_an_Br;   mm_rnu =mm_rnu_an_Br
                     E_nuBIN=E_nuBIN_an_Br; P_Fermi=P_Fermi_an_Br
                     P_FeMAX=P_FeMAX_an_Br; T_Fermi=T_Fermi_an_Br
                     FV_SM  =FV_ma_Br
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Br,*104)
                     S_Br=factor*S_Br
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               S_CF_3Br=(6*S__C+3*9*S__F+35*S_Br)/(6+3*9+35)

               dsQESCC_dQ2_SM_ma=S_CF_3Br
*              ------------------------------------------------------- *
               CASE(   6) ! D_2 (DEUTERIUM)
*              ------------------------------------------------------- *
               dsQESCC_dQ2_SM_ma= dsQESCC_dQ2_FP_ma(E_nu,Q2,MA_QES)
*              ------------------------------------------------------- *
               CASE( 7:8) ! Fe (IRON), STEEL
*              ------------------------------------------------------- *
               IF (E_nu.le.E_ma_Fe_thr) THEN                             Fe (IRON),      Z=26
                 S_Fe=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Fe=Precision
                                                     ELSE
                   m_tar  =m_tar_Fe;      mm_tar =mm_tar_Fe
                   m_rnu  =m_rnu_an_Fe;   mm_rnu =mm_rnu_an_Fe
                   E_nuBIN=E_nuBIN_an_Fe; P_Fermi=P_Fermi_an_Fe
                   P_FeMAX=P_FeMAX_an_Fe; T_Fermi=T_Fermi_an_Fe
                   FV_SM  =FV_ma_Fe
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Fe,*104)
                   S_Fe=factor*S_Fe
              endIF
            endIF
               dsQESCC_dQ2_SM_ma=S_Fe
*              ------------------------------------------------------- *
               CASE(   9) ! Ne-H_2 (NEOUN-HYDROGEN MIXTURE)
*              ------------------------------------------------------- *
               S__p=dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)                    FREE proton
               S__H=S__p                                                 H  (HYDROGEN),  Z= 1
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma__O_thr) THEN                           O  (OXYGEN),    Z= 8
                   S__O=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__O=Precision
                                                       ELSE
                     m_tar  =m_tar__O;      mm_tar =mm_tar__O
                     m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                     E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                     P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                     FV_SM  =FV_ma__O
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                     S__O=factor*S__O
                endIF
              endIF
                 IF (E_nu.le.E_ma_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_an_Mg;   mm_rnu =mm_rnu_an_Mg
                     E_nuBIN=E_nuBIN_an_Mg; P_Fermi=P_Fermi_an_Mg
                     P_FeMAX=P_FeMAX_an_Mg; T_Fermi=T_Fermi_an_Mg
                     FV_SM  =FV_ma_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 S_Ne=(S__O+S_Mg)/2                                      Ne (NEON),      Z=10
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma_Ne_thr) THEN                           Ne (NEON),      Z=10
                   S_Ne=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ne=Precision
                                                       ELSE
                     m_tar  =m_tar_Ne;      mm_tar =mm_tar_Ne
                     m_rnu  =m_rnu_an_Ne;   mm_rnu =mm_rnu_an_Ne
                     E_nuBIN=E_nuBIN_an_Ne; P_Fermi=P_Fermi_an_Ne
                     P_FeMAX=P_FeMAX_an_Ne; T_Fermi=T_Fermi_an_Ne
                     FV_SM  =FV_ma_Ne
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ne,*104)
                     S_Ne=factor*S_Ne
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               C_Ne=0.64; N_Ne=10                                        Ne (NEON),      Z=10
               C__H=0.36; N__H= 1                                        H  (HYDROGEN),  Z= 1
               S_Ne_H=(C_Ne*N_Ne*S_Ne+C__H*N__H*S__H)/
     #                (C_Ne*N_Ne     +C__H*N__H     )

               dsQESCC_dQ2_SM_ma=S_Ne_H
*              ------------------------------------------------------- *
               CASE(  10) ! O (OXIGEN)
*              ------------------------------------------------------- *
               IF (E_nu.le.E_ma__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 dsQESCC_dQ2_SM_ma=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   dsQESCC_dQ2_SM_ma=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                   E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                   P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                   FV_SM  =FV_ma__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   dsQESCC_dQ2_SM_ma=factor*S__O
              endIF
            endIF
*              ------------------------------------------------------- *
c              CASE(  12) ! N  (NITROGEN)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  13) ! Cl (CHLORINE)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  14) ! Si (SILICON)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  15) ! Cu (COPPER)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
c              CASE(  16) ! CH_2 (MINERAL OIL)
*              ------------------------------------------------------- *
*              ------------------------------------------------------- *
               CASE(  17) ! NOMAD MIXTURE
*              ------------------------------------------------------- *
               S__p=dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)                    FREE proton
               S__H=S__p                                                 H  (HYDROGEN),  Z= 1
               IF (E_nu.le.E_ma__C_thr) THEN                             C  (CARBON),    Z= 6
                 S__C=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__C=Precision
                                                     ELSE
                   m_tar  =m_tar__C;      mm_tar =mm_tar__C
                   m_rnu  =m_rnu_an__C;   mm_rnu =mm_rnu_an__C
                   E_nuBIN=E_nuBIN_an__C; P_Fermi=P_Fermi_an__C
                   P_FeMAX=P_FeMAX_an__C; T_Fermi=T_Fermi_an__C
                   FV_SM  =FV_ma__C
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__C,*104)
                   S__C=factor*S__C
              endIF
            endIF
               IF (E_nu.le.E_ma__O_thr) THEN                             O  (OXYGEN),    Z= 8
                 S__O=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S__O=Precision
                                                     ELSE
                   m_tar  =m_tar__O;      mm_tar =mm_tar__O
                   m_rnu  =m_rnu_an__O;   mm_rnu =mm_rnu_an__O
                   E_nuBIN=E_nuBIN_an__O; P_Fermi=P_Fermi_an__O
                   P_FeMAX=P_FeMAX_an__O; T_Fermi=T_Fermi_an__O
                   FV_SM  =FV_ma__O
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__O,*104)
                   S__O=factor*S__O
              endIF
            endIF
               IF (E_nu.le.E_ma_Ar_thr) THEN                             Ni (NICKEL),    Z=28
                 S_Ar=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Ar=Precision
                                                     ELSE
                   m_tar  =m_tar_Ar;      mm_tar =mm_tar_Ar
                   m_rnu  =m_rnu_an_Ar;   mm_rnu =mm_rnu_an_Ar
                   E_nuBIN=E_nuBIN_an_Ar; P_Fermi=P_Fermi_an_Ar
                   P_FeMAX=P_FeMAX_an_Ar; T_Fermi=T_Fermi_an_Ar
                   FV_SM  =FV_ma_Ar
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ar,*104)
                   S_Ar=factor*S_Ar
              endIF
            endIF
               IF (E_nu.le.E_ma_Si_thr) THEN                             Si (SILICON),   Z=14
                 S_Si=Precision
                                        ELSE
                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                   S_Si=Precision
                                                     ELSE
                   m_tar  =m_tar_Si;      mm_tar =mm_tar_Si
                   m_rnu  =m_rnu_an_Si;   mm_rnu =mm_rnu_an_Si
                   E_nuBIN=E_nuBIN_an_Si; P_Fermi=P_Fermi_an_Si
                   P_FeMAX=P_FeMAX_an_Si; T_Fermi=T_Fermi_an_Si
                   FV_SM  =FV_ma_Si
                   CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Si,*104)
                   S_Si=factor*S_Si
              endIF
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               IF (OLDSCHEME) THEN
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma_Si_thr) THEN                           Si (SILICON),   Z=14
                   S_Si=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Si=Precision
                                                       ELSE
                     m_tar  =m_tar_Si;      mm_tar =mm_tar_Si
                     m_rnu  =m_rnu_an_Si;   mm_rnu =mm_rnu_an_Si
                     E_nuBIN=E_nuBIN_an_Si; P_Fermi=P_Fermi_an_Si
                     P_FeMAX=P_FeMAX_an_Si; T_Fermi=T_Fermi_an_Si
                     FV_SM  =FV_ma_Si
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Si,*104)
                     S_Si=factor*S_Si
                endIF
              endIF
                 IF (E_nu.le.E_ma_Ca_thr) THEN                           Ca (CALCIUM),   Z=20
                   S_Ca=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ca=Precision
                                                       ELSE
                     m_tar  =m_tar_Ca;      mm_tar =mm_tar_Ca
                     m_rnu  =m_rnu_an_Ca;   mm_rnu =mm_rnu_an_Ca
                     E_nuBIN=E_nuBIN_an_Ca; P_Fermi=P_Fermi_an_Ca
                     P_FeMAX=P_FeMAX_an_Ca; T_Fermi=T_Fermi_an_Ca
                     FV_SM  =FV_ma_Ca
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ca,*104)
                     S_Ca=factor*S_Ca
                endIF
              endIF
                 IF (E_nu.le.E_ma_Mg_thr) THEN                           Mg (MAGNESIUM), Z=12
                   S_Mg=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Mg=Precision
                                                       ELSE
                     m_tar  =m_tar_Mg;      mm_tar =mm_tar_Mg
                     m_rnu  =m_rnu_an_Mg;   mm_rnu =mm_rnu_an_Mg
                     E_nuBIN=E_nuBIN_an_Mg; P_Fermi=P_Fermi_an_Mg
                     P_FeMAX=P_FeMAX_an_Mg; T_Fermi=T_Fermi_an_Mg
                     FV_SM  =FV_ma_Mg
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Mg,*104)
                     S_Mg=factor*S_Mg
                endIF
              endIF
                 IF (E_nu.le.E_ma_Ni_thr) THEN                           Ni (NICKEL),    Z=28
                   S_Ni=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Ni=Precision
                                                       ELSE
                     m_tar  =m_tar_Ni;      mm_tar =mm_tar_Ni
                     m_rnu  =m_rnu_an_Ni;   mm_rnu =mm_rnu_an_Ni
                     E_nuBIN=E_nuBIN_an_Ni; P_Fermi=P_Fermi_an_Ni
                     P_FeMAX=P_FeMAX_an_Ni; T_Fermi=T_Fermi_an_Ni
                     FV_SM  =FV_ma_Ni
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Ni,*104)
                     S_Ni=factor*S_Ni
                endIF
              endIF
                 S__N=(  S__C+  S__O)/2                                  N  (NITROGEN),  Z= 7
                 S_Al=(  S_Mg+  S_Si)/2                                  Al (ALUMINUM),  Z=13
                 S_Cl=(3*S_Si+3*S_Ca)/6                                  Cl (CHLORINE),  Z=17
                 S_Cu=S_Ni                                               Cu (COPPER),    Z=29
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                              ELSE
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
                 IF (E_nu.le.E_ma__N_thr) THEN                           N  (NITROGEN),  Z= 7
                   S__N=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S__N=Precision
                                                       ELSE
                     m_tar  =m_tar__N;      mm_tar =mm_tar__N
                     m_rnu  =m_rnu_an__N;   mm_rnu =mm_rnu_an__N
                     E_nuBIN=E_nuBIN_an__N; P_Fermi=P_Fermi_an__N
                     P_FeMAX=P_FeMAX_an__N; T_Fermi=T_Fermi_an__N
                     FV_SM  =FV_ma__N
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S__N,*104)
                     S__N=factor*S__N
                endIF
              endIF
                 IF (E_nu.le.E_ma_Al_thr) THEN                           Al (ALUMINUM),  Z=13
                   S_Al=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Al=Precision
                                                       ELSE
                     m_tar  =m_tar_Al;      mm_tar =mm_tar_Al
                     m_rnu  =m_rnu_an_Al;   mm_rnu =mm_rnu_an_Al
                     E_nuBIN=E_nuBIN_an_Al; P_Fermi=P_Fermi_an_Al
                     P_FeMAX=P_FeMAX_an_Al; T_Fermi=T_Fermi_an_Al
                     FV_SM  =FV_ma_Al
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Al,*104)
                     S_Al=factor*S_Al
                endIF
              endIF
                 IF (E_nu.le.E_ma_Cl_thr) THEN                           Cl (CHLORINE),  Z=17
                   S_Cl=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Cl=Precision
                                                       ELSE
                     m_tar  =m_tar_Cl;      mm_tar =mm_tar_Cl
                     m_rnu  =m_rnu_an_Cl;   mm_rnu =mm_rnu_an_Cl
                     E_nuBIN=E_nuBIN_an_Cl; P_Fermi=P_Fermi_an_Cl
                     P_FeMAX=P_FeMAX_an_Cl; T_Fermi=T_Fermi_an_Cl
                     FV_SM  =FV_ma_Cl
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Cl,*104)
                     S_Cl=factor*S_Cl
                endIF
              endIF
                 IF (E_nu.le.E_ma_Cu_thr) THEN                           Cu (COPPER),    Z=29
                   S_Cu=Precision
                                          ELSE
                   CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
                   IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
                     S_Cu=Precision
                                                       ELSE
                     m_tar  =m_tar_Cu;      mm_tar =mm_tar_Cu
                     m_rnu  =m_rnu_an_Cu;   mm_rnu =mm_rnu_an_Cu
                     E_nuBIN=E_nuBIN_an_Cu; P_Fermi=P_Fermi_an_Cu
                     P_FeMAX=P_FeMAX_an_Cu; T_Fermi=T_Fermi_an_Cu
                     FV_SM  =FV_ma_Cu
                     CALL MuLInt(MuL_dsQESCC_dQ2_SM,S_Cu,*104)
                     S_Cu=factor*S_Cu
                endIF
              endIF
*                - - - - - - - - - - - - - - - - - - - - - - - - - - - *
            endIF
*              - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
               C__H= 5.14d-02; N__H= 1                                   H  (HYDROGEN),  Z= 1
               C__C=64.30d-02; N__C= 6                                   C  (CARBON),    Z= 6
               C__N= 5.92d-00; N__N= 7                                   N  (NITROGEN),  Z= 7
               C__O=22.13d-02; N__O= 8                                   O  (OXYGEN),    Z= 8
               C_Al= 1.71d-00; N_Al=13                                   Al (ALUMINUM),  Z=13
               C_Si= 0.27d-00; N_Si=14                                   Si (SILICON),   Z=14
               C_Cl= 0.30d-00; N_Cl=17                                   Cl (CHLORINE),  Z=17
               C_Ar= 0.19d-00; N_Ar=18                                   Ar (ARGON),     Z=18
               C_Cu= 0.03d-00; N_Cu=28                                   Cu (COPPER),    Z=29

               dsQESCC_dQ2_SM_ma=(C__C*N__C*S__C+
     #                              C__H*N__H*S__H+
     #                              C__O*N__O*S__O+
     #                              C__N*N__N*S__N+
     #                              C_Cl*N_Cl*S_Cl+
     #                              C_Al*N_Al*S_Al+
     #                              C_Si*N_Si*S_Si+
     #                              C_Ar*N_Ar*S_Ar+
     #                              C_Cu*N_Cu*S_Cu )/
     #                             (C__C*N__C     +
     #                              C__H*N__H     +
     #                              C__O*N__O     +
     #                              C__N*N__N     +
     #                              C_Cl*N_Cl     +
     #                              C_Al*N_Al     +
     #                              C_Si*N_Si     +
     #                              C_Ar*N_Ar     +
     #                              C_Cu*N_Cu      )
*              ------------------------------------------------------- *
      endSELECT
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_SM_ta(E_nu,Q2,MA_QES)
*     ==================================================================
         m_ini     = m_p;   mm_ini    = mm_p
         m_fin     = m_n;   mm_fin    = mm_n
         m_lep     = m_tau; mm_lep    = mm_e
         n_NT      =-1
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES

         dsQESCC_dQ2_SM_ta= zero
         RETURN

*     ==================================================================

      END FUNCTION dsQESCC_dQ2_SM_set
