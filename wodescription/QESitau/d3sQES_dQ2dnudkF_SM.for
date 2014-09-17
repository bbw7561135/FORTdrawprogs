************************************************************************
      FUNCTION d3sQES_dQ2dnudkF_SM(E_nu,Q2,v,kF)
************************************************************************
*                                                                      *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] R.A. Smith and E.J. Moniz, "Neutrino reactions on nucle-     *
*         ar  targets,"  Nucl. Phys. B 43 (1972) 605-622;  Erratum     *
*         Nucl. Phys. B 101 (1975) 547.                                *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one,half,quarter,
     #                                m_I,mm_I,mm_W,c2C,pi

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON   /n_SM_TMD/n_SM_TMD                                     Switch for target nucleon momentum distribution
         COMMON   /n_SM_PEF/n_SM_PEF                                     Switch for Pauli exclusion factor for outgoing nucleon
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy
         COMMON      /FV_SM/FV_SM                                        Normalization factor for nuclear volume

*        ------------------------------------------------------------- *
*        QES FORM FACTORS                                              *
*        ------------------------------------------------------------- *
         CALL FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
         F_V      = ReF_V
         F_M      = ReF_M
         F_A      = ReF_A
         F_P      = ReF_P
         F_T      = ReF_T
         F_S      = ReF_S
*        ------------------------------------------------------------- *
*        SMITH-MONIZ STRUCTURE FUNCTIONS                               *
*        ------------------------------------------------------------- *
         qv       = sqrt(v**2+Q2)                                        |\vec{q}|
         E_p      = sqrt(mm_ini+kF**2)-E_nuBIN

         cosT_p   = ((v-E_nuBIN)*(2*E_p+v+E_nuBIN)-qv**2+mm_ini-mm_fin)/
     #              (2*kF*qv)                                            \cos\theta_p
         pF       = sqrt(kF**2+2*kF*qv*cosT_p+qv**2)
         cosT_k   = (v+(Q2+mm_lep)/(2*E_nu))/qv

         b2_flux  = (E_p-kF*cosT_k*cosT_p)**2
         c2_flux  = kF**2*(one-cosT_p**2)*(one-cosT_k**2)

         one_v_rel= m_tar*kF*pi/(FV_SM*qv*sqrt(b2_flux-c2_flux))

         rho_kF   = rho_SM(n_SM_TMD,kF)
         rho_pF   = rho_SM(n_SM_PEF,pF)

         t        = Q2/(4*mm_I)

         T_1      = F_A**2*(one+t)+t*(F_V+F_M)**2                        Ref.[1], \tilde{T}_1
         T_2      = F_A**2+F_V**2+t*F_M**2                               Ref.[1], \tilde{T}_2
         T_3      =-2*F_A*(F_V+F_M)                                      Ref.[1], \tilde{T}_8
         T_4      =-half*F_V*F_M-F_A*F_P+t*F_P**2-quarter*(one-t)*F_M**2 Ref.[1], \tilde{T}_\alpha
         T_5      = F_V**2+t*F_M**2+F_A**2

         a1       = one
         a2       = (kF/m_I)**2
         a3       = (kF*cosT_p/m_I)**2
         a4       = (E_p/m_I)**2
         a5       = E_p*kF*cosT_p/mm_I
         a6       = kF*cosT_p/m_I
         a7       = E_p/m_I

         W_1      = a1*T_1+half*(a2-a3)*T_2                              Ref.[1], W_1
         W_2      = (Q2/(2*qv**2)*(a2-a3)+a4-2*v/qv*a5+(v/qv)**2*a3)*T_2 Ref.[1], W_2
         W_3      = m_tar/m_I*(a7-(v/qv)*a6)*T_3                         Ref.[1], W_8
         W_4      = mm_tar/qv**2*(3*a3-a2)*half*T_2+mm_tar/mm_I*a1*T_4+  Ref.[1], W_\alpha
     #              mm_tar/(m_I*qv)*a6*T_5
         W_5      = m_tar*(a7-v*a6/qv)*T_5/m_I+
     #              (m_tar*v*(a2-3*a3)/qv**2+2*m_tar*a5/qv)*T_2
*        ------------------------------------------------------------- *
*        d2\sigma_dE_\ell d\cos\theta_\ell                             *
*        ------------------------------------------------------------- *
         E_lep      = E_nu-v
         P_lep      = sqrt(E_lep**2-mm_lep)
         cosT_lep   = (E_lep-(Q2+mm_lep)/(2*E_nu))/P_lep
         d2s_dEdcosT= (E_lep-P_lep*cosT_lep)/m_tar*
     #                (W_1+mm_lep/(2*mm_tar)*W_4)+
     #                (E_lep+P_lep*cosT_lep)/(2*m_tar)*W_2+
     #                n_NT*W_3*
     #                ((E_nu+E_lep)*(E_lep-P_lep*cosT_lep)-mm_lep)/
     #                 (2*mm_tar)-
     #                mm_lep/(2*mm_tar)*W_5
*        ------------------------------------------------------------- *
         d3sQES_dQ2dnudkF_SM= (c2C*mm_I)*
     #                        one_v_rel*rho_kF*(one-rho_pF)*
     #                        d2s_dEdcosT*(mm_W/(mm_W+Q2))**2/E_nu

         RETURN
      END FUNCTION d3sQES_dQ2dnudkF_SM