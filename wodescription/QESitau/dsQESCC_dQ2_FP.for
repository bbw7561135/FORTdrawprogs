************************************************************************
      FUNCTION dsQESCC_dQ2_FP_set(E_nu,Q2,MA_QES)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         REAL,PARAMETER:: factor_D2= G_Fermi**2*c2C*mm_I/(8*pi)*
     #                               hc2*1.0d+38

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

         dsQESCC_dQ2_FP_set= one
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_FP_mn(E_nu,Q2,MA_QES)                            mn + n --> mu^-  p
*     ==================================================================
         n_NT      =+1
         MA_QES_tmp= MA_QES
         m_ini     = m_n;   mm_ini= mm_n
         m_fin     = m_p;   mm_fin= mm_p
         m_lep     = m_mu;  mm_lep= mm_mu
         IF (E_nu.le.E_mn_thr) THEN
           dsQESCC_dQ2_SM_mn= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_FP_mn= Precision
                                               ELSE
             dsQESCC_dQ2_FP_mn= factor_D2*dsQESCC_dQ2(E_nu,Q2)*
     #                            FactorPauli_D2(Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_FP_ma(E_nu,Q2,MA_QES)                            ma + p --> mu^+  n
*     ==================================================================
         n_NT      =-1
         MA_QES_tmp= MA_QES
         m_ini     = m_p;   mm_ini= mm_p
         m_fin     = m_n;   mm_fin= mm_n
         m_lep     = m_mu;  mm_lep= mm_mu
         IF (E_nu.le.E_ma_thr) THEN
           dsQESCC_dQ2_SM_ma= Precision
                                 ELSE
           CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
           IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
             dsQESCC_dQ2_FP_ma= Precision
                                               ELSE
             dsQESCC_dQ2_FP_ma= factor_D2*dsQESCC_dQ2(E_nu,Q2)*
     #                            FactorPauli_D2(Q2)/E_nu**2
        endIF
      endIF
         RETURN

*     ==================================================================

      END FUNCTION dsQESCC_dQ2_FP_set
