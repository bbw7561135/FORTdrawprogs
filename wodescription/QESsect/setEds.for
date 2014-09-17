************************************************************************
      SUBROUTINE setEds
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL,PARAMETER::
     #                CoefMA = 1.00d+00
                 REAL
     #                MA_cen(5)/0.90,0.95,1.10,1.20,1.35/

         COMMON   /E_nu_thr/E_nu_thr,P_cher,O_lep,P_kth                  Neutrino energy threshold, Cherenkov threshold of lepton momentum
         COMMON       /n_MA/n_MA                                         Switch for MA_QES
         COMMON     /MA_QES/MA_QES                                       Mass of axial-vector in QES CC reactions
         COMMON     /MA_ELS/MA_ELS                                       Mass of axial-vector in ELS NC reactions
         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON       /n_fl/n_fl                                         Switch for neutrino flavor
         COMMON       /n_PT/n_PT                                         Switch for lepton polarization type
         COMMON   /n_AG_QES/n_AG                                         Switch for model of QES reactions
         COMMON   /n_MC_QES/n_MC                                         Switch for correction of masses of initial and final nucleons
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON   /n_AP_QES/n_AP                                         Switch for model of axial form factor in QES reactions
         COMMON   /n_MS_QES/n_MS                                         Switch for value of axial mass in Sehgal's model
         COMMON   /n_GE_QES/n_GE                                         Switch for parametrization of Sachs electric form factor of neutron
         COMMON     /MV_QES/MV_QES                                       Mass of isovector in QES reactions
         COMMON     /MM_QES/MM_QES                                       Mass of monopole axial-vector in QES reactions
         COMMON     /MS_QES/MS_QES                                       Mass of scalar in QES reactions
         COMMON     /MT_QES/MT_QES                                       Mass of tensor in QES reactions
         COMMON       /xi_V/xi_V                                         Normalization of vector form factor
         COMMON       /xi_M/xi_M                                         Normalization of monopole form factor
         COMMON       /xi_S/xi_S                                         Normalization of scalar form factor
         COMMON       /xi_A/xi_A                                         Normalization of axial form factor
         COMMON       /xi_P/xi_P                                         Normalization of pseudoscalar form factor
         COMMON       /xi_T/xi_T                                         Normalization of tensor form factor
         COMMON      /phi_T/phi_T                                        Phase of tensor form factor
         COMMON      /phi_S/phi_S                                        Phase of scalar form factor

         MA_QES  = CoefMA*MA_cen(n_MA)

         MA_ELS  = MA_QES
         n_PT    = 0
         n_AG    = 1
         n_MC    = 1
         n_AP    = 1
         n_MS    = 1
         n_GE    = 5
         MV_QES  = 8.4d-01
         MM_QES  = 8.0d-01
         MS_QES  = one
         MT_QES  = 1.5d+00
         xi_V    = one
         xi_M    = one
         xi_S    = zero
         xi_A    = one
         xi_P    = one
         xi_T    = zero
         phi_T   = zero
         phi_T   = phi_T*dtr
         phi_S   = zero
         phi_S   = phi_S*dtr

         IF (n_NT.EQ.1) THEN
           m_ini  = m_n; mm_ini = mm_n
           m_fin  = m_p; mm_fin = mm_p
                        ELSE
           m_ini  = m_p; mm_ini = mm_p
           m_fin  = m_n; mm_fin = mm_n
      endIF
         SELECTCASE(n_fl)
*              ------------------------------------------------------- *
               CASE(1)                                                   electron (anti)neutrino
*              ------------------------------------------------------- *
               m_lep = m_e; mm_lep= mm_e
               IF (n_NT.eq.1) THEN
                 E_nu_thr= Precision
                              ELSE
                 E_nu_thr=half*((m_fin+m_lep)**2-mm_ini)/m_ini
            endIF
*              ------------------------------------------------------- *
               CASE(2)                                                   muon (anti)neutrino
*              ------------------------------------------------------- *
               m_lep = m_mu
               mm_lep= mm_mu
               E_nu_thr=half*((m_fin+m_lep)**2-mm_ini)/m_ini
*              ------------------------------------------------------- *
               CASE(3)                                                   tau (anti)neutrino
*              ------------------------------------------------------- *
               m_lep = m_tau
               mm_lep= mm_tau
               E_nu_thr=half*((m_fin+m_lep)**2-mm_ini)/m_ini+5*Precision
*              ------------------------------------------------------- *
      endSELECT

         c_med=one/n_w
         gf_cm=one/sqrt(one-c_med**2)
         P_cher=gf_cm*m_lep*c_med
         O_lep= half*mm_lep/m_I
         P_kth= half*m_I-O_lep

         RETURN
      END SUBROUTINE setEds
