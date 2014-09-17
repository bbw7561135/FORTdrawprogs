************************************************************************
      FUNCTION fui(x)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         EXTERNAL MA_QES_EFF

         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON    /P_lep/P_lep,E_lep                                    Charged lepton momentum
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON   /MA_QES/MA_QES                                         Mass of axial-vector in QES CC reactions


         fui= one
         RETURN

*     ==================================================================
      ENTRY fuifen(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifen= dFANomen_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuifen= one
         RETURN
                  
*     ==================================================================
      ENTRY fuifea(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifea= dFANomea_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuifea= one
         RETURN

*     ==================================================================
      ENTRY fuifmn(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifmn= dFANommn_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuifmn= one
         RETURN
                  
*     ==================================================================
      ENTRY fuifma(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifma= dFANomma_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuifma= one
         RETURN

*     ==================================================================
      ENTRY fuiftn(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuiftn= dFANomtn_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuiftn= one
         RETURN
                  
*     ==================================================================
      ENTRY fuifta(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifta= dFANomta_dE(E_nu)*dsQESCC_dQ2(E_nu,Q2)/(den**2*P_lep)
!         fuifta= one
         RETURN

*     ==================================================================
      ENTRY fuiben(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuiben= dFANomen_dE(E_nu)*dsQESCC_dQ2_SM_en(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuiben= one
         RETURN
                  
*     ==================================================================
      ENTRY fuibea(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibea= dFANomea_dE(E_nu)*dsQESCC_dQ2_SM_ea(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuibea= one
         RETURN

*     ==================================================================
      ENTRY fuibmn(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibmn= dFANommn_dE(E_nu)*dsQESCC_dQ2_SM_mn(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuibmn= one
         RETURN
                  
*     ==================================================================
      ENTRY fuibma(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibma= dFANomma_dE(E_nu)*dsQESCC_dQ2_SM_ma(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuibma= one
         RETURN

*     ==================================================================
      ENTRY fuibtn(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibtn= dFANomtn_dE(E_nu)*dsQESCC_dQ2_SM_tn(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuibtn= one
         RETURN
                  
*     ==================================================================
      ENTRY fuibta(x)
*     ==================================================================
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         fuibta= dFANomta_dE(E_nu)*dsQESCC_dQ2_SM_ta(E_nu,Q2,
     #                                  MA_QES_EFF(E_nu))*E_nu**2/P_lep !-
!         fuibta= one
         RETURN

      END FUNCTION fui