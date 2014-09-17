************************************************************************
      FUNCTION fui(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         EXTERNAL MA_QES_EFF

         COMMON    /x_lim/x_ini,deltax                                   Limits (for Q2)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON   /MA_QES/MA_QES                                         Mass of axial-vector in QES CC reactions
         COMMON     /E_nu/E_nu

         fui= one
         RETURN

*     ==================================================================
      ENTRY fuifen(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifen= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuifea(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifea= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN

*     ==================================================================
      ENTRY fuifmn(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifmn= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuifma(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifma= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN

*     ==================================================================
      ENTRY fuiftn(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuiftn= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuifta(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
         fuifta= dsQESCC_dQ2(E_nu,Q2)/(den*E_nu)**2                     !-!
         RETURN

*     ==================================================================
      ENTRY fuiben(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuiben= dsQESCC_dQ2_SM_en(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuibea(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuibea= dsQESCC_dQ2_SM_ea(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN

*     ==================================================================
      ENTRY fuibmn(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuibmn= dsQESCC_dQ2_SM_mn(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuibma(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuibma= dsQESCC_dQ2_SM_ma(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN

*     ==================================================================
      ENTRY fuibtn(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuibtn= dsQESCC_dQ2_SM_tn(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN
                  
*     ==================================================================
      ENTRY fuibta(var)
*     ==================================================================
         x=deltax*var+x_ini
         Q2=x
         E_lep=E_nu-0.5*Q2/m_ini
         P_lep=sqrt(E_lep**2-mm_lep)
         fuibta= dsQESCC_dQ2_SM_ta(E_nu,Q2,MA_QES_EFF(E_nu))            !-!
         RETURN

      END FUNCTION fui