************************************************************************
      FUNCTION dsQESCC_dQ2_fN_set(E_nu,Q2,MA_QES)
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         REAL,PARAMETER:: factor= G_Fermi**2*c2C*mm_I/(8*pi)*
     #                            hc2*1.0d+38

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON     /MA_QES/MA_QES_tmp                                   Mass of axial-vector in QES CC reactions

         E_en_thr= ((m_p+m_e  )**2-mm_n)/(2*m_n)
         E_mn_thr= ((m_p+m_mu )**2-mm_n)/(2*m_n)
         E_tn_thr= ((m_p+m_tau)**2-mm_n)/(2*m_n)

         E_ea_thr= ((m_n+m_e  )**2-mm_p)/(2*m_p)
         E_ma_thr= ((m_n+m_mu )**2-mm_p)/(2*m_p)
         E_ta_thr= ((m_n+m_tau)**2-mm_p)/(2*m_p)

         dsQESCC_dQ2_fN_set= one
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_en(E_nu,Q2,MA_QES)                            en + n --> e^-   + p
*     ==================================================================
         n_NT      =+1                                                   Switch for neutrino type
         m_ini     = m_n;   mm_ini= mm_n                                 Mass of target nucleon
         m_fin     = m_p;   mm_fin= mm_p                                 Mass of final hadron or hadron system
         m_lep     = m_e;   mm_lep= mm_e
         MA_QES_tmp= MA_QES
         CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
         IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
           dsQESCC_dQ2_fN_en= Precision
                                             ELSE
           dsQESCC_dQ2_fN_en= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_ea(E_nu,Q2,MA_QES)                            ea + p --> e^+   + n
*     ==================================================================
         n_NT      =-1                                                   Switch for neutrino type
         m_ini     = m_p;   mm_ini= mm_p                                 Mass of target nucleon
         m_fin     = m_n;   mm_fin= mm_n                                 Mass of final hadron or hadron system
         m_lep     = m_e;   mm_lep= mm_e
         MA_QES_tmp= MA_QES
         IF (E_nu.le.E_ea_thr) THEN
           dsQESCC_dQ2_fN_ea=Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_fN_ea= Precision
                                               ELSE
             dsQESCC_dQ2_fN_ea= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_mn(E_nu,Q2,MA_QES)                            mn + n --> mu^-  + p
*     ==================================================================
         n_NT      =+1                                                   Switch for neutrino type
         m_ini     = m_n;   mm_ini= mm_n                                 Mass of target nucleon
         m_fin     = m_p;   mm_fin= mm_p                                 Mass of final hadron or hadron system
         m_lep     = m_mu;  mm_lep= mm_mu
         MA_QES_tmp= MA_QES
         IF (E_nu.le.E_mn_thr) THEN
           dsQESCC_dQ2_fN_mn= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_fN_mn= Precision
                                               ELSE
             dsQESCC_dQ2_fN_mn= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)                            ma + p --> mu^+  n
*     ==================================================================
         n_NT      =-1                                                   Switch for neutrino type
         m_ini     = m_p;   mm_ini= mm_p                                 Mass of target nucleon
         m_fin     = m_n;   mm_fin= mm_n                                 Mass of final hadron or hadron system
         m_lep     = m_mu;  mm_lep= mm_mu
         MA_QES_tmp= MA_QES
         IF (E_nu.le.E_ma_thr) THEN
           dsQESCC_dQ2_fN_ma= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_fN_ma= Precision
                                               ELSE
             dsQESCC_dQ2_fN_ma= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_tn(E_nu,Q2,MA_QES)                            tn + n --> tau^-  + p
*     ==================================================================
         n_NT      =+1                                                   Switch for neutrino type
         m_ini     = m_n;   mm_ini= mm_n                                 Mass of target nucleon
         m_fin     = m_p;   mm_fin= mm_p                                 Mass of final hadron or hadron system
         m_lep     = m_tau; mm_lep= mm_tau
         MA_QES_tmp= MA_QES
         IF (E_nu.le.E_tn_thr) THEN
           dsQESCC_dQ2_fN_tn= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_fN_tn= Precision
                                               ELSE
             dsQESCC_dQ2_fN_tn= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN_ta(E_nu,Q2,MA_QES)                            ta + p --> tau^+  n
*     ==================================================================
         n_NT      =-1                                                   Switch for neutrino type
         m_ini     = m_p;   mm_ini= mm_p                                 Mass of target nucleon
         m_fin     = m_n;   mm_fin= mm_n                                 Mass of final hadron or hadron system
         m_lep     = m_tau; mm_lep= mm_tau
         MA_QES_tmp= MA_QES
         IF (E_nu.le.E_ta_thr) THEN
           dsQESCC_dQ2_fN_ta= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_fN_ta= Precision
                                               ELSE
             dsQESCC_dQ2_fN_ta= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
        endIF
      endIF
         RETURN
*     ==================================================================

      END FUNCTION dsQESCC_dQ2_fN_set
