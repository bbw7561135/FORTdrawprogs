************************************************************************
      FUNCTION fui(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         EXTERNAL MA_QES_EFF

         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON    /P_lep/P_lep,E_lep                                    Charged lepton momentum
         COMMON    /x_lim/x_ini,deltax                                   Limits (for neutrino energy)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON   /MA_QES/MA_QES                                         Mass of axial-vector in QES CC reactions


         fui= one
         RETURN

*     ==================================================================
      ENTRY fuifen(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifen= dFANomen_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN
                  
*     ==================================================================
      ENTRY fuifea(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifea= dFANomea_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN

*     ==================================================================
      ENTRY fuifmn(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifmn= dFANommn_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN
                  
*     ==================================================================
      ENTRY fuifma(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifma= dFANomma_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN

*     ==================================================================
      ENTRY fuiftn(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuiftn= dFANomtn_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN
                  
*     ==================================================================
      ENTRY fuifta(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifta= dFANomta_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)  !-
         RETURN

*     ==================================================================
      ENTRY fuiben(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuiben= dFANomen_dE(E_nu)*dsQESCC_dQ2_SM_en(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN
                  
*     ==================================================================
      ENTRY fuibea(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibea= dFANomea_dE(E_nu)*dsQESCC_dQ2_SM_ea(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN

*     ==================================================================
      ENTRY fuibmn(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibmn= dFANommn_dE(E_nu)*dsQESCC_dQ2_SM_mn(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN
                  
*     ==================================================================
      ENTRY fuibma(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibma= dFANomma_dE(E_nu)*dsQESCC_dQ2_SM_ma(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN

*     ==================================================================
      ENTRY fuibtn(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibtn= dFANomtn_dE(E_nu)*dsQESCC_dQ2_SM_tn(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN
                  
*     ==================================================================
      ENTRY fuibta(var)
*     ==================================================================
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibta= dFANomta_dE(E_nu)*dsQESCC_dQ2_SM_ta(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
         RETURN

      END FUNCTION fui