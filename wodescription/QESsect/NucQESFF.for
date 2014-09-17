************************************************************************
      SUBROUTINE NucQESFF(Q2,
     #                        ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                        ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
************************************************************************
*                                                                      *
*     This SUBROUTINE  returns  the vector  and axial-vector  form     *
*     factors of the nucleon using several phenomenological models     *
*     and parametrizations of the experimental data.                   *
*                                                                      *
*     There are two models for axial form factor:  standard dipole     *
*     form [n_AP=1] and fit by L.M. Sehgal [n_AP=2]. The parameter     *
*     n_MS  switchs the  model of  the parameter m_A in  this fit.     *
*     We use the following notations  for the mesons: "w0" denotes     *
*     \omega(782), "r1" - \rho (1450), "w1" - \omega(1420), "f0" -     *
*     \phi(1020).  Here "Q2"  denotes Q^2 (> 0).  Symbols "Re" and     *
*     "Im" denote real and imaginary parts of the form factors.        *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [ 1] C.H. Llewellyn Smith,  "Neutrino reactions at accelera-     *
*          tor energies," Phys. Rept. 3 (1972) 261-379  [SLAC-PUB-     *
*          958 (1971)].                                                *
*     [ 2] E.L. Lomon, "Effect of recent  $R_p$ and $R_n$ measure-     *
*          ments  on extended  Gari-Krumpelmann model fits to nuc-     *
*          leon  electromagnetic  form  factors,"  Phys. Rev. C 66     *
*          (2002) 045501 [arXiv: nucl-th/0203081].                     *
*     [ 3] N.J. Baker  et al., "Quasielastic neutrino  scattering:     *
*          A measurement of the weak nucleon axial-vector form fa-     *
*          ctor," Phys. Rev. D23 (1981) 2499-2505.                     *
*     [ 4] L.M. Sehgal,  in Proc. of the European Phys. Soc. Int.      *
*          Conf.in High Energy Physics (1979), p.98.                   *
*     [ 5] S.K. Singh and E. Oset, "Quasielastic neutrino (antine-     *
*          utrino) reactions in nuclei and the  axial vector  form     *
*          factor of the nucleon," Nucl. Phys.A 542(1992) 587-615.     *
*     [ 6] T. Kitagaki et al., "Study  of  neutrino $\nu d\to\mu^-     *
*          pp_s$ $\nu d\to \mu^-\Delat^++(1232)n_s$  using the BNL     *
*          7-foot  deuterium-filled bubble chamber,"  Phys. Rev. D     *
*          42 (1990) 1331-1338.                                        *
*     [ 7] G.M. Radecky et al., "Study of  single  pion production     *
*          by weak charged currents in low-energy neutrino $\nu d$     *
*          interactions," Phys. Rev. D 25 (1982) 1161-1173.            *
*     [ 8] L.A. Ahrens et al., "A study  of the  axial-vector form     *
*          factor  and second-class  currents in antineutrino qua-     *
*          sielastic scattering," Phys. Lett. B202 (1988) 284-288.     *
*     [ 9] H. Budd, A. Bodek  and J. Arrington,  "Modeling  quasi-     *
*          elastic form factors  for electron and neutrino scatte-     *
*          ring", arXiv: hep-ex/0308005.                               *
*     [10] G. Warren et al.  (The Jefferson Lab E93-026 Collabora-     *
*          tion), "Measurement of the electric form factor  of the     *
*          neutron at $Q^2=0.5$ and $1.0$ (GeV/c)$^2$," Phys. Rev.     *
*          Lett. 92 (2004) 042301 [arXiv: nucl-ex/0308021].            *
*     [11] S. Galster et al.,  "Elastic electron-deuteron  scatte-     *
*          ring and  the electric neutron form  factor at four mo-     *
*          mentum  transfers 5 fm$^{-2}<q^2<14$ fm$^{-2}$,"  Nucl.     *
*          Phys. B 32 (1971) 221-237.                                  *
*     [12] E.L. Lomon, "Effect  of revised $R_n$  measurements  on     *
*          extended  Gari-Kruempelmann model fits to nucleon elec-     *
*          tromagnetic form factors," arXiv: nucl-th/0609020.          *
*     [13] A. Bodek, S. Avvakumov, R. Bradford and H. Budd,  "Dua-     *
*          lity constrained  parameterization  of vector and axial     *
*          nucleon  form  factors," Eur. Phys. J. C 53 (2008) 349-     *
*          354 [arXiv: 0708.1946[hep-ex]].                             *
*     [14] G. Vereshkov and O. Lalakulich,  "Logarithmic correcti-     *
*          ons and soft photon phenomenology in  the  multipole mo-    *
*          del of  the nucleon form  factors,"  Eur. Phys. J. A 34     *
*    :     (2007) 223-236 [arXiv: 0705.1476[hep-ph]].                  *
*     [15] A. Bradford, A. Bodek, H. Budd and J. Arrington, "A new     *
*          parametrization  of the nucleon  elastic form factors,"     *
*          Nucl. Phys. B (Proc. Suppl.) 159 (2006) 127-132 [arXiv:     *
*          hep-ex/0602017].                                            *
*     [16] J.J. Kelly, "Simple  parametrization  of  nucleon  form     *
*          factors," Phys. Rev. C 70 (2004) 068202.                    *
*     [17] K.M. Graczyk, P. Plonski, R. Sulej, "Neural network pa-     *
*          rameterizations  of electromagnetic  nucleon form  fac-     *
*          tors," JHEP 53 (2010) 1009 [arXiv: 1006.0342[hep-ph]].      *
*     [18] C. Crawford et al., "The role of mesons in the electro-     *
*          magnetic form factors of the nucleon," arXiv: 1003.0903     *
*          [nucl-th] (Submitted to Phys. Rev. C).                      *
*     [19] A. Bodek, H.S. Budd and M.E. Christy,  "Neutrino quasi-     *
*          elastic  scattering on  nuclear targets," Eur. Phys. J.     *
*          C 71 (2011) 1726 [arVix: 1106.0340[hep-ph]].                *
*                                                                      *
*                                         Angarsk, Russia, 2011/07/14  *
*                                    ITEP, Moscow, Russia, 2011/10/05  *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         REAL,PARAMETER::
     #        Gp0_M = one+k_p,                                           G^p_M  (Q^2=0)
     #        Gp0_E = one,                                               G^p_E  (Q^2=0)
     #        Gn0_M = k_n,                                               G^n_M  (Q^2=0)
     #        Gn0_E = zero,                                              G^n_E  (Q^2=0)
     #        G0_M  = Gp0_M-Gn0_M,                                       G_M    (Q^2=0)=1+k_p-k_n
     #        G0_E  = Gp0_E-Gn0_E,                                       G_E    (Q^2=0)=1
     #        FCC_A0=-1.2695d+00,                                        F^CC_A (Q^2=0)
     #        M_BNL7= 1.37d+00,                                          Axial mass according to Ref.[6]
     #        mu    = k_p-mu_n,
     #        x1    = zero,
     #        x2    = 1.0d+00/six,
     #        x3    = 2.0d+00/six,
     #        x4    = 3.0d+00/six,
     #        x5    = 4.0d+00/six,
     #        x6    = 5.0d+00/six,
     #        x7    = one,
     #        w01GnE= 1.019704d+01,
     #        w02GnE= 2.368120d+00,
     #        w03GnE=-1.144266d+00,
     #        w04GnE=-4.274101d+00,
     #        w05GnE= 8.149924d-01,
     #        w06GnE= 2.985524d+00,
     #        w07GnE=-7.864434d-01,
     #        w01GnM= 3.196460d+00,
     #        w02GnM= 2.565681d+00,
     #        w03GnM= 6.441526d+00,
     #        w04GnM=-2.004055d+00,
     #        w05GnM=-2.972361d-01,
     #        w06GnM= 3.606737d+00,
     #        w07GnM=-3.135199d+00,
     #        w08GnM= 2.995230d-01,
     #        w09GnM= 1.261638d+00,
     #        w10GnM= 2.647470d+00,
     #        w01GpE= 3.930227d+00,
     #        w02GpE= 1.108384d-01,
     #        w03GpE=-5.325479d+00,
     #        w04GpE=-2.846154d+00,
     #        w05GpE=-2.071328d-01,
     #        w06GpE= 8.742101d-01,
     #        w07GpE= 4.283194d-01,
     #        w08GpE= 2.568322d+00,
     #        w09GpE= 2.577635d+00,
     #        w10GpE=-1.185632d+00,
     #        w01GpM=-2.862682d+00,
     #        w02GpM=-1.560675d+00,
     #        w03GpM= 2.321148d+00,
     #        w04GpM= 1.283189d-01,
     #        w05GpM=-2.803566d-01,
     #        w06GpM= 2.794296d+00,
     #        w07GpM= 1.726774d+00,
     #        w08GpM= 8.610830d-01,
     #        w09GpM= 4.184286d-01,
     #        w10GpM=-1.526676d-01

         COMMON   /n_RT_QES/n_RT_QES                                     Switch for QES reaction type
         COMMON   /n_FF_QES/n_FF_QES                                     Switch for model of nucleon form factors in QES reactions
         COMMON   /n_AP_QES/n_AP                                         Switch for model of axial form factor in QES reactions
         COMMON   /n_MS_QES/n_MS                                         Switch for value of axial mass in Sehgal's model
         COMMON   /n_GE_QES/n_GE                                         Switch for parametrization of Sachs electric form factor of neutron
         COMMON     /MA_QES/MA_QES                                       Mass of axial-vector in QES CC reactions
         COMMON     /MA_ELS/MA_ELS                                       Mass of axial-vector in ELS NC reactions
         COMMON     /MV_QES/MV_QES                                       Mass of isovector in QES reactions
         COMMON     /MM_QES/MM_QES                                       Mass of monopole axial-vector in QES reactions
         COMMON     /MS_QES/MS_QES                                       Mass of scalar in QES reactions
         COMMON     /MT_QES/MT_QES                                       Mass of tensor in QES reactions
         COMMON     /gs_ELS/gs_ELS                                       Parameter for strange Fs_A form factor
         COMMON     /ms_ELS/ms_ELS                                       Parameter for strange Gs_M form factor
         COMMON     /rs_ELS/rs_ELS                                       Parameter for strange Gs_E form factor
         COMMON       /Gp_M/Gp_M                                         Sachs magnetic form factor of proton
         COMMON       /Gp_E/Gp_E                                         Sachs electric form factor of proton
         COMMON       /Gn_M/Gn_M                                         Sachs magnetic form factor of neutron
         COMMON       /Gn_E/Gn_E                                         Sachs electric form factor of neutron
         COMMON        /G_D/G_D                                          Standard dipole form
         COMMON        /G_M/G_M                                          Sachs magnetic form factor of nucleon
         COMMON        /G_E/G_E                                          Sachs electric form factor of nucleon
         COMMON       /xi_V/xi_V                                         Normalization of vector form factor
         COMMON       /xi_M/xi_M                                         Normalization of monopole form factor
         COMMON       /xi_S/xi_S                                         Normalization of scalar form factor
         COMMON       /xi_A/xi_A                                         Normalization of axial form factor
         COMMON       /xi_P/xi_P                                         Normalization of pseudoscalar form factor
         COMMON       /xi_T/xi_T                                         Normalization of tensor form factor
         COMMON      /phi_T/phi_T                                        Phase of tensor form factor
         COMMON      /phi_S/phi_S                                        Phase of scalar form factor
         COMMON      /ReFsA/ReFsA                                        Real part of strange axial-vector form factor of nucleon
         COMMON       /Gs_M/Gs_M                                         Strange Sachs magnetic form factor of nucleon
         COMMON       /Gs_E/Gs_E                                         Strange Sachs electric form factor of nucleon

          P6(x,a02,a04,a06,a08,a10,a12)=
     #    one+x*(a02+x*(a04+x*(a06+x*(a08+x*(a10+x*a12)))))

          AN(x,c1,c2,c3,c4,c5,c6,c7)=
     #    c1*       (x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7)+
     #    c2*(x-x1)*       (x-x3)*(x-x4)*(x-x5)*(x-x6)*(x-x7)+
     #    c3*(x-x1)*(x-x2)*       (x-x4)*(x-x5)*(x-x6)*(x-x7)+
     #    c4*(x-x1)*(x-x2)*(x-x3)*       (x-x5)*(x-x6)*(x-x7)+
     #    c5*(x-x1)*(x-x2)*(x-x3)*(x-x4)*       (x-x6)*(x-x7)+
     #    c6*(x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5)*       (x-x7)+
     #    c7*(x-x1)*(x-x2)*(x-x3)*(x-x4)*(x-x5)*(x-x6)

          KG(x,a0,a1,a2,b1,b2,b3,b4)=
     #    (a0+x*(a1+x*a2))/(one+x*(b1+x*(b2+x*(b3+x*b4))))

          f_GPS(x)=one/(one+exp(-x))

          PRINT *,' MODEL FOR THE NUCLEON FORM FACTORS IS'
*         ------------------------------------------------------------ *
          IF (n_FF_QES.eq. 0) THEN;   PRINT *, ' FIXED FF'
*         ------------------------------------------------------------ *
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 1) THEN;   PRINT *, ' DIPOLE model'
*         ------------------------------------------------------------ *
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 2) THEN;   PRINT *, ' BBA03p fit'
*         ------------------------------------------------------------ *
            GEp02= 3.253d+00; GMp02= 3.104d+00; GMn02= 3.043d+00
            GEp04= 1.422d+00; GMp04= 1.428d+00; GMn04= 8.548d-01
            GEp06= 8.582d-02; GMp06= 1.112d-01; GMn06= 6.806d-01
            GEp08= 3.318d-01; GMp08=-6.981d-03; GMn08=-1.287d-01
            GEp10=-9.371d-02; GMp10= 3.705d-04; GMn10= 8.912d-03
            GEp12= 1.076d-02; GMp12=-7.063d-06; GMn12= 0.000d+00
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 3) THEN;   PRINT *, ' BBA_03 fit'
*         ------------------------------------------------------------ *
            GEp02= 3.226d+00; GMp02= 3.188d+00; GMn02= 3.043d+00
            GEp04= 1.508d+00; GMp04= 1.354d+00; GMn04= 8.548d-01
            GEp06=-3.773d-01; GMp06= 1.511d-01; GMn06= 6.806d-01
            GEp08= 6.109d-01; GMp08=-1.135d-02; GMn08=-1.287d-01
            GEp10=-1.853d-01; GMp10= 5.330d-04; GMn10= 8.912d-03
            GEp12= 1.596d-02; GMp12=-9.005d-06; GMn12= 0.000d+00
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 4) THEN;  PRINT *, ' GKex(02L) fit'
*         ------------------------------------------------------------ *
            k_r1  = 5.30380d+00                                          \kappa_{\rho'}
            k_f0  = 1.30037d+01                                          \kappa_\phi
            k_w0  =-2.85850d+00                                          \kappa_\omega
            k_w1  = 1.82284d+01                                          \kappa_{\omega'}
            gf_r1 = 6.08000d-02                                          g_{\rho'}/f_{\rho'}
            gf_f0 =-1.85200d-01                                          g_\phi/f_\phi
            gf_w0 = 6.89600d-01                                          g_\omega/f_\omega
            gf_w1 = 2.34600d-01                                          g_{\omega'}/f_{\omega'}
            LL_1  = 9.44100d-01**2                                       \Lambda^2_1
            LL_2  = 2.82680d+00**2                                       \Lambda^2_2
            LL_D  = 1.23500d+00**2                                       \Lambda^2_D
            mu2_f0= 6.84800d-01**2                                       \mu_\phi
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 5) THEN;  PRINT *, ' GKex(02S) fit'
*         ------------------------------------------------------------ *
            k_r1  = 6.81900d+00                                          \kappa_{\rho'}
            k_f0  = 7.01720d+00                                          \kappa_\phi
            k_w0  = 8.76200d-01                                          \kappa_\omega
            k_w1  = 1.49160d+00                                          \kappa_{\omega'}
            gf_r1 = 4.01000d-02                                          g_{\rho'}/f_{\rho'}
            gf_f0 =-1.67600d-01                                          g_\phi/f_\phi
            gf_w0 = 6.73900d-01                                          g_\omega/f_\omega
            gf_w1 = 2.55200d-01                                          g_{\omega'}/f_{\omega'}
            LL_1  = 9.40700d-01**2                                       \Lambda^2_1
            LL_2  = 2.78910d+00**2                                       \Lambda^2_2
            LL_D  = 1.21110d+00**2                                       \Lambda^2_D
            mu2_f0= 8.54400d-01**2                                       \mu_\phi
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 6) THEN;  PRINT *, ' GKex(05) fit'
*         ------------------------------------------------------------ *
            k_r1  = 1.20000d+01                                          \kappa_{\rho'}
            k_f0  = 1.00000d-02                                          \kappa_\phi
            k_w0  = 4.02700d-01                                          \kappa_\omega
            k_w1  =-2.97300d+00                                          \kappa_{\omega'}
            gf_r1 = 7.20890d-03                                          g_{\rho'}/f_{\rho'}
            gf_f0 =-1.71100d-01                                          g_\phi/f_\phi
            gf_w0 = 7.02100d-01                                          g_\omega/f_\omega
            gf_w1 = 1.64000d-01                                          g_{\omega'}/f_{\omega'}
            LL_1  = 9.30880d-01**2                                       \Lambda^2_1
            LL_2  = 2.61150d+00**2                                       \Lambda^2_2
            LL_D  = 1.18100d+00**2                                       \Lambda^2_D
            mu2_f0= 2.00000d-01**2                                       \mu_\phi
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 7) THEN;  PRINT *, ' BBBA05 fit'
*         ------------------------------------------------------------ *
            a0Mp= 1.00d+00;a0Ep= 1.00d+00;a0Mn= 1.00d+00;a0En= 0.00d+00
            a1Mp= 1.50d-01;a1Ep=-5.78d-02;a1Mn= 1.81d+00;a1En= 1.25d+00
            a2Mp= 0.00d+00;a2Ep= 0.00d+00;a2Mn= 0.00d+00;a2En= 1.30d+00
            b1Mp= 1.11d+01;b1Ep= 1.11d+01;b1Mn= 1.41d+01;b1En=-9.86d+00
            b2Mp= 1.96d+01;b2Ep= 1.36d+01;b2Mn= 2.07d+01;b2En= 3.05d+02
            b3Mp= 7.54d+00;b3Ep= 3.30d+01;b3Mn= 6.87d+01;b3En=-7.58d+02
            b4Mp= 0.00d+00;b4Ep= 0.00d+00;b4Mn= 0.00d+00;b4En= 8.02d+02
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 8) THEN;  PRINT *, ' BBBA07_25 fit'
*         ------------------------------------------------------------ *
*           A^p_M(dipole)  A^p_E(dipole)  A^n_M(25)      A^p_E(25)     *
            p1Mp= 1.0000;p1Ep=   1.0000;p1Mn=   1.0000;p1En=   1.0000
            p2Mp= 1.0011;p2Ep=   0.9927;p2Mn=   0.9958;p2En=   1.1011
            p3Mp= 0.9992;p3Ep=   0.9898;p3Mn=   0.9877;p3En=   1.1392
            p4Mp= 0.9974;p4Ep=   0.9975;p4Mn=   1.0193;p4En=   1.0203
            p5Mp= 1.0010;p5Ep=   0.9812;p5Mn=   1.0350;p5En=   1.1093
            p6Mp= 1.0003;p6Ep=   0.9340;p6Mn=   0.9164;p6En=   1.5429
            p7Mp= 1.0000;p7Ep=   1.0000;p7Mn=   0.7300;p7En=   0.9706
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq. 9) THEN;  PRINT *, ' BBBA07_43 fit'
*         ------------------------------------------------------------ *
*           A^p_M(dipole)  A^p_E(dipole)  A^n_M(43)      A^p_E(43)     *
            p1Mp= 1.0000;p1Ep=   1.0000;p1Mn=   1.0000;p1En=   1.0000
            p2Mp= 1.0011;p2Ep=   0.9927;p2Mn=   0.9959;p2En=   1.1019
            p3Mp= 0.9992;p3Ep=   0.9898;p3Mn=   0.9851;p3En=   1.1387
            p4Mp= 0.9974;p4Ep=   0.9975;p4Mn=   1.0187;p4En=   1.0234
            p5Mp= 1.0010;p5Ep=   0.9812;p5Mn=   1.0307;p5En=   1.1046
            p6Mp= 1.0003;p6Ep=   0.9340;p6Mn=   0.9080;p6En=   1.5395
            p7Mp= 1.0000;p7Ep=   1.0000;p7Mn=   0.9557;p7En=   1.2708
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq.10) THEN; PRINT *, ' VL 2007 model'
*         ------------------------------------------------------------ *

*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq.11) THEN; PRINT *, ' GPS 2010 fit'
*         ------------------------------------------------------------ *

*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq.12) THEN; PRINT *, ' GKex 2010 fit'            Ref.[18]
*         ------------------------------------------------------------ *
            mm_r0 = mm_rho
            mm_r01= (m_rho-0.03465)**2                                   m^2_\rho-\deta_1
            mm_r02= (m_rho-0.04374)**2                                   m^2_\rho-\deta_2
            k_r0  =  5.51564                                             \kappa_{\rho}
            k_r1  = 12.0000                                              \kappa_{\rho'}
            k_f0  =  0.0100                                              \kappa_\phi
            k_w0  =  0.4027                                              \kappa_\omega
            k_w1  =- 2.973                                               \kappa_{\omega'}
            gf_r0 =  0.5596                                              g_{\rho}/f_{\rho}
            gf_r1 =  0.0072089                                           g_{\rho'}/f_{\rho'}
            gf_w0 =  0.7021                                              g_\omega/f_\omega
            gf_w1 =  0.164                                               g_{\omega'}/f_{\omega'}
            gf_f0 =- 0.1711                                              g_\phi/f_\phi
            LL_1  =  0.93088**2                                          \Lambda^2_1
            LL_2  =  2.6115**2                                           \Lambda^2_2
            LL_D  =  1.1810**2                                           \Lambda^2_D
            LL_QCD=  0.1500**2                                           \Lambda^2_QCD
            mu2_f0=  0.2000**2                                           \mu_\phi
            Q2_V1 =  0.3176                                              Q^2_{V1}
            Q2_M1 =  0.1422                                              Q^2_{M1}
            a_1   =  0.0781808                                           \alpha_1
            a_2   =  0.0632907                                           \alpha_2
            a_3   =  one-a_1                                             1-\alpha_1
            a_4   =  one-a_2                                             1-\alpha_2
            k_v   =  3.7060                                              \kappa_v
            k_s   =- 0.1200                                              \kappa_s
            Cs_V  =  one-     gf_w0-     gf_w1
            Cs_M  =  k_s-k_w0*gf_w0-k_w1*gf_w1-k_f0*gf_f0
            Cv_V  =  one-     gf_r0-     gf_r1                           1-g_{\rho}/f_{\rho}-g_{\rho'}/f_{\rho'}
            Cv_M  =  k_v-k_r0*gf_r0-k_r1*gf_r1                           1-\kappa_{\rho}g_{\rho}/f_{\rho}-\kappa_{\rho'}g_{\rho'}/f_{\rho'}
            lnLL  =  log(LL_D/LL_QCD)
*         ------------------------------------------------------------ *
      ELSEIF (n_FF_QES.eq.13) THEN; PRINT *, ' BBC 2011'                 Ref.[19]
*         ------------------------------------------------------------ *
*           A^p_M(dipole)  A^p_E(dipole)  A^n_M(25)      A^p_E(25)     *
            p1Mp=   1.0000;p1Ep=   1.0000;p1Mn=   1.0000;p1En=   1.0000
            p2Mp=   0.9916;p2Ep=   0.9839;p2Mn=   0.9958;p2En=   1.1011
            p3Mp=   0.9771;p3Ep=   0.9632;p3Mn=   0.9877;p3En=   1.1392
            p4Mp=   0.9801;p4Ep=   0.9748;p4Mn=   1.0193;p4En=   1.0203
            p5Mp=   1.0321;p5Ep=   0.9136;p5Mn=   1.0350;p5En=   1.1093
            p6Mp=   1.0429;p6Ep=   0.5447;p6Mn=   0.9164;p6En=   1.5429
            p7Mp=   0.5084;p7Ep=  -0.2682;p7Mn=   0.7300;p7En=   0.9706

            A_TE_C= 6.000d+00; B_TE_C=-one/3.400d-01
*         ------------------------------------------------------------ *
       endIF
          IF (n_FF_QES.ge. 4 .and. n_FF_QES.le. 6) THEN
            A_V   = 1.0317d+00
            A_M   = 5.7824d+00
            B_V   = 8.7500d-02
            B_M   = 3.9070d-01
            Q2_V1 = 3.1760d-01                                           Q^2_{V1}
            Q2_M1 = 1.4220d-01                                           Q^2_{M1}
            Q2_V2 = 5.4960d-01                                           Q^2_{V2}
            Q2_M2 = 5.3620d-01                                           Q^2_{M2}
            k_v   = 3.7060d+00                                           \kappa_v
            k_s   =-1.2000d-01                                           \kappa_s
            LL_QCD= 1.5000d-01**2                                        \Lambda^2_QCD
            C     = half                                                 N/2
            Cs_V  = one-gf_w0-gf_w1
            Cs_M  = k_s-k_w0*gf_w0-k_w1*gf_w1-k_f0*gf_f0
            Cv_V  = one-1.1192*C-gf_r1
            Cv_M  = k_v-6.1731*C-k_r1*gf_r1
            lnLL  = log(LL_D/LL_QCD)
       endIF
          IF (n_FF_QES.eq. 8 .or. n_FF_QES.eq. 9 .or.
     #        n_FF_QES.eq.13                         ) THEN
            a1Mp= 1.717d-01; a1Ep=-2.400d-01
            b1Mp= 1.126d+01; b1Ep= 1.098d+01
            b2Mp= 1.932d+01; b2Ep= 1.282d+01
            b3Mp= 8.330d+00; b3Ep= 2.197d+01

            d1=(x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)*(x1-x7)
            d2=(x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)*(x2-x7)
            d3=(x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)*(x3-x6)*(x3-x7)
            d4=(x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)*(x4-x6)*(x4-x7)
            d5=(x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)*(x5-x6)*(x5-x7)
            d6=(x6-x1)*(x6-x2)*(x6-x3)*(x6-x4)*(x6-x5)*(x6-x7)
            d7=(x7-x1)*(x7-x2)*(x7-x3)*(x7-x4)*(x7-x5)*(x7-x6)

            c1Mp=p1Mp/d1; c1Ep=p1Ep/d1; c1Mn=p1Mn/d1; c1En=p1En/d1
            c2Mp=p2Mp/d2; c2Ep=p2Ep/d2; c2Mn=p2Mn/d2; c2En=p2En/d2
            c3Mp=p3Mp/d3; c3Ep=p3Ep/d3; c3Mn=p3Mn/d3; c3En=p3En/d3
            c4Mp=p4Mp/d4; c4Ep=p4Ep/d4; c4Mn=p4Mn/d4; c4En=p4En/d4
            c5Mp=p5Mp/d5; c5Ep=p5Ep/d5; c5Mn=p5Mn/d5; c5En=p5En/d5
            c6Mp=p6Mp/d6; c6Ep=p6Ep/d6; c6Mn=p6Mn/d6; c6En=p6En/d6
            c7Mp=p7Mp/d7; c7Ep=p7Ep/d7; c7Mn=p7Mn/d7; c7En=p7En/d7
       endIF
          IF (n_AP.eq.3) THEN
            SELECTCASE(n_MS)                                             Switch for value of axial mass in Sehgal's model
                  CASE(   1);mm_A= mm_a1                                 Ref.[5]
                  CASE(   2);mm_A= 2*mm_rho                              Ref.[7]
                  CASE(   3);mm_A= M_BNL7**2                             Ref.[6]
         endSELECT
       endIF
         SELECTCASE(n_GE)                                                Switch for parametrization of Sachs electric form factor of neutron
               CASE(   1);a= 1.00d+00; b= 5.60d+00                       Galster et al., Ref.[11]
               CASE(   2);a= 1.25d+00; b= 1.83d+01                       Platchikov et al., Ref.[]
               CASE(   3);a= 1.00d+00; b= 3.40d+00                       Herberg et al., Ref.[]
               CASE(   4);a= 9.42d-01; b= 4.65d+00                       Krutov and Troitsky, Ref.[]
               CASE(   5);a= 9.42d-01; b= 4.61d+00                       Warren et al., Ref.[10]
               CASE(   6);a= 8.95d-01; b= 3.69d+00                       Budd, Bodek and Arrington, Ref.[9]
      endSELECT

         RETURN

*     ==================================================================
      ENTRY FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #              ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
*     ==================================================================
         MMV_QES= MV_QES**2                                              Mass of isovector in QES reactions
         MMS_QES= MS_QES**2                                              Mass of scalar in QES reactions
         MMT_QES= MT_QES**2                                              Mass of tensor in QES reactions
         MMA_QES= MA_QES**2
         G_D    = (MMV_QES/(MMV_QES+Q2))**2
         SELECTCASE(n_FF_QES)
*              ------------------------------------------------------- *
               CASE(       0)                                            Fixed form factors
*              ------------------------------------------------------- *
c              G_M = G0_M
c              G_E = G0_E
               G_E = (MMV_QES/(MMV_QES-mm_e))**2
               G_M = (one+mu)*G_E
*              ------------------------------------------------------- *
               CASE(       1)                                            Standard dipole formula, Ref.[1]
*              ------------------------------------------------------- *
               G_E = (MMV_QES/(MMV_QES+Q2))**2
               G_M = (one+mu)*G_E
*              ------------------------------------------------------- *
               CASE(     2,3)                                            BBA (2003), Ref[9]
*              ------------------------------------------------------- *
               Gn_M= Gn0_M/P6(Q2,GMn02,GMn04,GMn06,GMn08,GMn10,GMn12)
               Gp_E= Gp0_E/P6(Q2,GEp02,GEp04,GEp06,GEp08,GEp10,GEp12)
               IF (Q2 .ge. 2.0d+01) THEN
c                Gp_M= Gp0_M/P6(Q2,GMp02,GMp04,GMp06,GMp08,GMp10,GMp12)  BBA
                 Gp_M= mu_p*G_D/(0.304*Q2-2.500)**0.222                  KLN
                                    ELSE
                 Gp_M= Gp0_M/P6(Q2,GMp02,GMp04,GMp06,GMp08,GMp10,GMp12)
            endIF
               t   = Q2/(4*mm_n)
               Gn_E=-mu_n*(a*t/(one+b*t))*G_D
               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(     4:6)                                            Gari-Kruempelmann-Lomon model, Ref.[2]
*              ------------------------------------------------------- *
               IF (Q2.le.zero) THEN
                 G_E = (MMV_QES/(MMV_QES+Q2))**2
                 G_M = (one+mu)*G_E
                 GOTO 1000
                               ELSE
                 Q2_t= Q2*log((LL_D+Q2)/LL_QCD)/lnLL
                 FD_V= LL_D*LL_2/((LL_D+Q2_t)*(LL_2+Q2_t))
                 Fa_V= LL_1*LL_2/((LL_1+Q2_t)*(LL_2+Q2_t))
                 FD_M= LL_D/(LL_D+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Fa_M= LL_1/(LL_1+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Ff_V= Fa_V*(Q2/(LL_1+Q2))**1.5
                 Ff_M= Fa_M*(((mu2_f0+Q2)/mu2_f0)/((LL_1+Q2)/LL_1))**1.5

                 Fs_1= (     (gf_w0*mm_w0/(mm_w0+Q2))+
     #                       (gf_w1*mm_w1/(mm_w1+Q2)))*Fa_V+
     #                 (     (gf_f0*mm_f0/(mm_f0+Q2)))*Ff_V+Cs_V*FD_V
                 Fs_2= (k_w0*(gf_w0*mm_w0/(mm_w0+Q2))+
     #                  k_w1*(gf_w1*mm_w1/(mm_w1+Q2)))*Fa_M+
     #                 (k_f0*(gf_f0*mm_f0/(mm_f0+Q2)))*Ff_M+Cs_M*FD_M
                 Fv_1= (C*(A_V+B_V*(Q2_V1/(Q2_V1+Q2))**2)*Q2_V2/
     #                 (Q2_V2+Q2)+
     #                 (gf_r1*mm_r1/(mm_r1+Q2)))*Fa_V+Cv_V*FD_V
                 Fv_2= (C*(A_M+B_M*(Q2_M1/(Q2_M1+Q2))   )*Q2_M2/
     #                 (Q2_M2+Q2)+
     #                  k_r1*(gf_r1*mm_r1/(mm_r1+Q2)))*Fa_M+Cv_M*FD_M

                 Fp_1= half*(Fs_1+Fv_1)
                 Fp_2= half*(Fs_2+Fv_2)
                 Fn_1= half*(Fs_1-Fv_1)
                 Fn_2= half*(Fs_2-Fv_2)
                 Gp_M= Fp_1+Fp_2
                 Gp_E= Fp_1-Fp_2*Q2/(4*mm_I)
                 Gn_M= Fn_1+Fn_2
                 Gn_E= Fn_1-Fn_2*Q2/(4*mm_I)
            endIF
               G_M= Gp_M-Gn_M
               G_E= Gp_E-Gn_E
c              PRINT *, ' Q^2  =', Q2,' GeV^2'
c              PRINT *, ' Fp_1 =', Fp_1, ' (must be 1)'
c              PRINT *, ' Fn_1 =', Fn_1, ' (must be 0)'
c              PRINT *, ' Fp_2 =', Fp_2, ' (must be ',k_p,')'
c              PRINT *, ' Fn_2 =', Fn_2, ' (must be ',k_n,')'
c              PRINT *, ' Gp_M =', Gp_M, ' (must be ',1+k_p,')'
c              PRINT *, ' Gp_E =', Gp_E, ' (must be 1)'
c              PRINT *, ' Gn_M =', Gn_M, ' (must be ',k_n,')'
c              PRINT *, ' Gn_E =', Gn_E, ' (must be 0)'; STOP
c              ------------------------------------------------------- *
               CASE(       7)                                            BBBA (2005), Ref.[15,16]
*              ------------------------------------------------------- *
               tI  = Q2/(4*mm_I)

               Gp_M= KG(tI,a0Mp,a1Mp,a2Mp,b1Mp,b2Mp,b3Mp,b4Mp)*mu_p
               Gp_E= KG(tI,a0Ep,a1Ep,a2Ep,b1Ep,b2Ep,b3Ep,b4Ep)
               Gn_M= KG(tI,a0Mn,a1Mn,a2Mn,b1Mn,b2Mn,b3Mn,b4Mn)*mu_n
               Gn_E= KG(tI,a0En,a1En,a2En,b1En,b2En,b3En,b4En)

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(     8:9)                                            BBBA_{25,43} (2007), Ref.[13]
*              ------------------------------------------------------- *
               Q2  = max(Q2,Precision)
               tp  = Q2/(4*mm_p)
               tn  = Q2/(4*mm_n)
               xp  = two/(one+sqrt(one+one/tp))
               xn  = two/(one+sqrt(one+one/tn))

               GMp = (one+a1Mp*tp)/(one+tp*(b1Mp+tp*(b2Mp+b3Mp*tp)))
               GEp = (one+a1Ep*tp)/(one+tp*(b1Ep+tp*(b2Ep+b3Ep*tp)))

               AMp = AN(xp,c1Mp,c2Mp,c3Mp,c4Mp,c5Mp,c6Mp,c7Mp)
               AEp = AN(xp,c1Ep,c2Ep,c3Ep,c4Ep,c5Ep,c6Ep,c7Ep)
               AMn = AN(xn,c1Mn,c2Mn,c3Mn,c4Mn,c5Mn,c6Mn,c7Mn)
               AEn = AN(xn,c1En,c2En,c3En,c4En,c5En,c6En,c7En)

               Gp_M= AMp*GMp*mu_p
               Gp_E= AEp*GEp
               Gn_M= AMn*Gp_M*mu_n/mu_p
               Gn_E= AEn*Gp_E*1.7*tn/(one+3.3*tn)

*              PRINT *, ' Q^2  =', Q2,   ' GeV^2'
*              PRINT *, ' Gp_M =', Gp_M, ' (must be ',1+k_p,')'
*              PRINT *, ' Gp_E =', Gp_E, ' (must be 1)'
*              PRINT *, ' Gn_M =', Gn_M, ' (must be ',k_n,')'
*              PRINT *, ' Gn_E =', Gn_E, ' (must be 0)'; STOP

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(      10)                                            Vereshkov-Lalakulich (2007), Ref.[14]
*              ------------------------------------------------------- *
               t   = Q2/(4*mm_I)
               z   = log(one+Q2/2.660d-02)
               p   = 1.100d+01-6.700d-01*
     #               (2.000d+00+Q2/(9.000d-02+Q2)+Q2/(9.000d+00+Q2))
               p_1 = (3.555d+00/p+2.000d+00)*5.000d-01
               p_2 = (2.667d+00/p)*5.000d-01
               p_1p= 1.160d-03*log(one+Q2/1.640d-05)**2

               Fs_1= (9.230d-01/(6.120d-01+Q2)-
     #                1.314d+00/(2.031d+00+Q2)+
     #                3.910d-01/(2.789d+00+Q2))/
     #               (1.000d+00+1.380d-02*z**2)**p_1
               Fs_2= (3.880d-01/(2.031d+00+Q2)-
     #                1.350d-01/(6.120d-01+Q2)-
     #                2.530d-01/(2.789d+00+Q2))/
     #               (1.000d+00-z*(5.680d-01-8.130d-02*z))**p_2
               Fv_1= (1.6018d+01/(2.147d+00+Q2)-
     #                7.930d-01 /(6.020d-01+Q2)-
     #                1.5226d+01/(2.958d+00+Q2))/
     #               (1.000d+00+9.490d-02*z**2)**p_1
               Fv_2= (3.893d+00/(6.020d-01+Q2)-
     #                1.1295d+01/(2.147d+00+Q2)+
     #                7.402d+00/(2.958d+00+Q2))/
     #               (one-z*(1.180d-01-3.260d-01*z))**p_2

               Fp_1= half*(Fs_1+Fv_1)
               Fp_2= half*(Fs_2+Fv_2)
               Fn_1= half*(Fs_1-Fv_1)
               Fn_2= half*(Fs_2-Fv_2)

               Gp_M= sqrt((Fp_1+  Fp_2)**2-
     #               2*(one-exp(-p_1p))*  Fp_1*Fp_2)
               Gp_E= sqrt((Fp_1-t*Fp_2)**2+
     #               2*(one-exp(-p_1p))*t*Fp_1*Fp_2)
               Gn_M= Fn_1+  Fn_2
               Gn_E= Fn_1-t*Fn_2

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E

*              PRINT *, ' Q^2  =', Q2,   ' GeV^2'
*              PRINT *, ' Fp_1 =', Fp_1, ' (must be 1)'
*              PRINT *, ' Fn_1 =', Fn_1, ' (must be 0)'
*              PRINT *, ' Fp_2 =', Fp_2, ' (must be ',k_p,')'
*              PRINT *, ' Fn_2 =', Fn_2, ' (must be ',k_n,')'
*              PRINT *, ' Gp_M =', Gp_M, ' (must be ',1+k_p,')'
*              PRINT *, ' Gp_E =', Gp_E, ' (must be 1)'
*              PRINT *, ' Gn_M =', Gn_M, ' (must be ',k_n,')'
*              PRINT *, ' Gn_E =', Gn_E, ' (must be 0)'; STOP
*              ------------------------------------------------------- *
               CASE(      11)                                            Graczyk et al., Ref.[17]
*              ------------------------------------------------------- *
               Gp_M= (w07GpM*f_GPS(Q2*w01GpM+w02GpM)+
     #                w08GpM*f_GPS(Q2*w03GpM+w04GpM)+
     #                w09GpM*f_GPS(Q2*w05GpM+w06GpM)+w10GpM)*mu_p*G_D
               Gp_E= (w07GpE*f_GPS(Q2*w01GpE+w02GpE)+
     #                w08GpE*f_GPS(Q2*w03GpE+w04GpE)+
     #                w09GpE*f_GPS(Q2*w05GpE+w06GpE)+w10GpE)*G_D
               Gn_M= (w07GnM*f_GPS(Q2*w01GnM+w02GnM)+
     #                w08GnM*f_GPS(Q2*w03GnM+w04GnM)+
     #                w09GnM*f_GPS(Q2*w05GnM+w06GnM)+w10GnM)*mu_n*G_D
               Gn_E=  w05GnE*f_GPS(Q2*w01GnE+w02GnE)+
     #                w06GnE*f_GPS(Q2*w03GnE+w04GnE)+w07GnE

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(      12)                                            C. Crawford et al., Ref.[12]
*              ------------------------------------------------------- *
               IF (Q2.le.zero) THEN
                 G_E = (MMV_QES/(MMV_QES+Q2))**2
                 G_M = (one+mu)*G_E
                 GOTO 1000
                               ELSE
                 Q2_t= Q2*log((LL_D+Q2)/LL_QCD)/lnLL
                 FD_V= LL_D*LL_2/((LL_D+Q2_t)*(LL_2+Q2_t))
                 Fa_V= LL_1*LL_2/((LL_1+Q2_t)*(LL_2+Q2_t))
                 FD_M= LL_D/(LL_D+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Fa_M= LL_1/(LL_1+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Ff_V= Fa_V*(Q2/(LL_1+Q2))**1.5
                 Ff_M= Fa_M*(((mu2_f0+Q2)/mu2_f0)/((LL_1+Q2)/LL_1))**1.5

                 Fs_1= (     (gf_w0*mm_w0/(mm_w0+Q2))+
     #                       (gf_w1*mm_w1/(mm_w1+Q2)))*Fa_V+
     #                 (     (gf_f0*mm_f0/(mm_f0+Q2)))*Ff_V+Cs_V*FD_V
                 Fs_2= (k_w0*(gf_w0*mm_w0/(mm_w0+Q2))+
     #                  k_w1*(gf_w1*mm_w1/(mm_w1+Q2)))*Fa_M+
     #                 (k_f0*(gf_f0*mm_f0/(mm_f0+Q2)))*Ff_M+Cs_M*FD_M

                 Fv_1= ((gf_r0*(mm_r01/(mm_r01+Q2)))*
     #                  (a_3+a_1*(Q2_V1/(Q2_V1+Q2))**2)+
     #                   (gf_r1*(mm_r1/(mm_r1+Q2))) )*Fa_V+Cv_V*FD_V
     #                
                 Fv_2= ( k_r0*(gf_r0*(mm_r02/(mm_r02+Q2)))*
     #                  (a_4+a_2*(Q2_M1/(Q2_M1+Q2))**2)+
     #                   k_r1*(gf_r1*(mm_r1/(mm_r1+Q2))))*Fa_M+Cv_M*FD_M
     
                 Fp_1= half*(Fs_1+Fv_1)
                 Fp_2= half*(Fs_2+Fv_2)
                 Fn_1= half*(Fs_1-Fv_1)
                 Fn_2= half*(Fs_2-Fv_2)
                 Gp_M= Fp_1+Fp_2
                 Gp_E= Fp_1-Fp_2*Q2/(4*mm_I)
                 Gn_M= Fn_1+Fn_2
                 Gn_E= Fn_1-Fn_2*Q2/(4*mm_I)
            endIF
               G_M= Gp_M-Gn_M
               G_E= Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(      13)                                            A. Bodek et al. (2011), Ref.[19]
*              ------------------------------------------------------- *
               Q2  = max(Q2,Precision)
               tp  = Q2/(4*mm_p)
               tn  = Q2/(4*mm_n)
               xp  = two/(one+sqrt(one+one/tp))
               xn  = two/(one+sqrt(one+one/tn))

               AMp = AN(xp,c1Mp,c2Mp,c3Mp,c4Mp,c5Mp,c6Mp,c7Mp)           A^p_M(dipole)
               AEp = AN(xp,c1Ep,c2Ep,c3Ep,c4Ep,c5Ep,c6Ep,c7Ep)           A^p_E(dipole)
               AMn = AN(xn,c1Mn,c2Mn,c3Mn,c4Mn,c5Mn,c6Mn,c7Mn)           A^n_M(25)
               AEn = AN(xn,c1En,c2En,c3En,c4En,c5En,c6En,c7En)           A^n_E(25)

               Gp_M= AMp*G_D*mu_p*sqrt(one+A_TE_C*Q2*exp(B_TE_C*Q2))     G^p_M
               Gp_E= AEp*G_D                                             G^p_E
               Gn_M= AMn*Gp_M*mu_n/mu_p                                  G^n_M
               Gn_E= AEn*Gp_E*1.7*tn/(one+3.3*tn)                        G^n_E

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
      endSELECT
 1000    ReF_V= xi_V*(G_E+Q2*G_M/(4*mm_I))/(one+Q2/(4*mm_I))
         ReF_M= xi_M*(G_M-   G_E         )/(one+Q2/(4*mm_I))
         SELECTCASE(n_AP)                                                Switch for model of axial form factor in QES reactions
*              ------------------------------------------------------- *
               CASE(   1)                                                Standard dipole formula, Ref.[1]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*(MMA_QES/(MMA_QES+Q2))**2
               ImF_A= zero
*              ------------------------------------------------------- *
               CASE(   2)                                                Monopole formula, Ref.[9]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*MM_QES/(MM_QES+Q2)
               ImF_A= zero
*              ------------------------------------------------------- *
               CASE(   3)                                                L.M. Sehgal's modification, Ref.[4,5]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*
     #                exp(max(-Q2*(one+Q2/(4*mm_I)),-1.0d+02))/
     #                        (one+Q2/mm_A)
               ImF_A= zero
*              ------------------------------------------------------- *
      endSELECT
         ReF_P= xi_P*2*ReF_A*mm_I/(Q2+mm_pi)
 
         ReF_T= xi_T*(MMT_QES/(MMT_QES+Q2))**2*cos(phi_T)*FCC_A0
         ReF_S= xi_S*(MMS_QES/(MMS_QES+Q2))**2*cos(phi_S)                Ref.[8]

         ImF_V= zero
         ImF_M= zero
         ImF_A= zero
         ImF_P= zero
         ImF_T= xi_T*(MMT_QES/(MMT_QES+Q2))**2*sin(phi_T)*FCC_A0
         ImF_S= xi_S*(MMS_QES/(MMS_QES+Q2))**2*sin(phi_S)                Ref.[8]

         RETURN

*     ==================================================================
      ENTRY FFNC(Q2,ReF_V,ReF_M,ReF_A,ReF_T,
     #              ImF_V,ImF_M,ImF_A,ImF_T)
*     ==================================================================
         MMV_QES= MV_QES**2                                              Mass of isovector in QES reactions
         MMT_QES= MT_QES**2                                              Mass of tensor in QES reactions
         MMA_ELS= MA_ELS**2
         G_D    = (MMV_QES/(MMV_QES+Q2))**2
         SELECTCASE(n_FF_QES)
*              ------------------------------------------------------- *
               CASE(     4:6)                                            Gari-Kruempelmann-Lomon model, Ref.[2]
*              ------------------------------------------------------- *
               IF (Q2.le.zero) THEN
                 G_E = (MMV_QES/(MMV_QES+Q2))**2
                 G_M = (one+mu)*G_E
                 GOTO 2000
                               ELSE
                 Q2_t= Q2*log((LL_D+Q2)/LL_QCD)/lnLL
                 FD_V= LL_D*LL_2/((LL_D+Q2_t)*(LL_2+Q2_t))
                 Fa_V= LL_1*LL_2/((LL_1+Q2_t)*(LL_2+Q2_t))
                 FD_M= LL_D/(LL_D+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Fa_M= LL_1/(LL_1+Q2_t)*(LL_2/(LL_2+Q2_t))**2
                 Ff_V= Fa_V*(Q2/(LL_1+Q2))**1.5
                 Ff_M= Fa_M*(((mu2_f0+Q2)/mu2_f0)/((LL_1+Q2)/LL_1))**1.5
                 Fs_1= (     (gf_w0*mm_w0/(mm_w0+Q2))+
     #                       (gf_w1*mm_w1/(mm_w1+Q2)))*Fa_V+
     #                 (     (gf_f0*mm_f0/(mm_f0+Q2)))*Ff_V+Cs_V*FD_V
                 Fs_2= (k_w0*(gf_w0*mm_w0/(mm_w0+Q2))+
     #                  k_w1*(gf_w1*mm_w1/(mm_w1+Q2)))*Fa_M+
     #                 (k_f0*(gf_f0*mm_f0/(mm_f0+Q2)))*Ff_M+Cs_M*FD_M
                 Fv_1= (C*(A_V+B_V*(Q2_V1/(Q2_V1+Q2))**2)*Q2_V2/
     #                 (Q2_V2+Q2)+
     #                 (gf_r1*mm_r1/(mm_r1+Q2)))*Fa_V+Cv_V*FD_V
                 Fv_2= (C*(A_M+B_M*(Q2_M1/(Q2_M1+Q2))   )*Q2_M2/
     #                 (Q2_M2+Q2)+
     #                  k_r1*(gf_r1*mm_r1/(mm_r1+Q2)))*Fa_M+Cv_M*FD_M
                 Fp_1= half*(Fs_1+Fv_1)
                 Fp_2= half*(Fs_2+Fv_2)
                 Fn_1= half*(Fs_1-Fv_1)
                 Fn_2= half*(Fs_2-Fv_2)
                 Gp_M= Fp_1+Fp_2
                 Gp_E= Fp_1-Fp_2*(Q2/(4*mm_I))
                 Gn_M= Fn_1+Fn_2
                 Gn_E= Fn_1-Fn_2*(Q2/(4*mm_I))
            endIF
               G_M= Gp_M-Gn_M
               G_E= Gp_E-Gn_E

c              PRINT *, ' Q^2  =', Q2,' GeV^2'
c              PRINT *, ' Fp_1 =', Fp_1, ' (must be 1)'
c              PRINT *, ' Fn_1 =', Fn_1, ' (must be 0)'
c              PRINT *, ' Fp_2 =', Fp_2, ' (must be ',k_p,')'
c              PRINT *, ' Fn_2 =', Fn_2, ' (must be ',k_n,')'
c              PRINT *, ' Gp_M =', Gp_M, ' (must be ',1+k_p,')'
c              PRINT *, ' Gp_E =', Gp_E, ' (must be 1)'
c              PRINT *, ' Gn_M =', Gn_M, ' (must be ',k_n,')'
c              PRINT *, ' Gn_E =', Gn_E, ' (must be 0)'; STOP
*              ------------------------------------------------------- *
               CASE(     8:9)                                            BBBA_{25,43} (2007), Ref.[13]
*              ------------------------------------------------------- *
               Q2  = max(Q2,Precision)
               tp  = Q2/(4*mm_p)
               tn  = Q2/(4*mm_n)
               xp  = two/(one+sqrt(one+one/tp))
               xn  = two/(one+sqrt(one+one/tn))

               GMp = (one+a1Mp*tp)/(one+tp*(b1Mp+tp*(b2Mp+b3Mp*tp)))
               GEp = (one+a1Ep*tp)/(one+tp*(b1Ep+tp*(b2Ep+b3Ep*tp)))

               AMp = AN(xp,c1Mp,c2Mp,c3Mp,c4Mp,c5Mp,c6Mp,c7Mp)
               AEp = AN(xp,c1Ep,c2Ep,c3Ep,c4Ep,c5Ep,c6Ep,c7Ep)
               AMn = AN(xn,c1Mn,c2Mn,c3Mn,c4Mn,c5Mn,c6Mn,c7Mn)
               AEn = AN(xn,c1En,c2En,c3En,c4En,c5En,c6En,c7En)

               Gp_M= AMp*GMp*mu_p
               Gp_E= AEp*GEp
               Gn_M= AMn*Gp_M*mu_n/mu_p
               Gn_E= AEn*Gp_E*1.7*tn/(one+3.3*tn)

               G_M = Gp_M-Gn_M
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
               CASE(      13)                                            A. Bodek et al. (2011), Ref.[19]
*              ------------------------------------------------------- *
               Q2  = max(Q2,Precision)
               tp  = Q2/(4*mm_p)
               tn  = Q2/(4*mm_n)
               xp  = two/(one+sqrt(one+one/tp))
               xn  = two/(one+sqrt(one+one/tn))

               AMp = AN(xp,c1Mp,c2Mp,c3Mp,c4Mp,c5Mp,c6Mp,c7Mp)           A^p_M(dipole)
               AEp = AN(xp,c1Ep,c2Ep,c3Ep,c4Ep,c5Ep,c6Ep,c7Ep)           A^p_E(dipole)
               AMn = AN(xn,c1Mn,c2Mn,c3Mn,c4Mn,c5Mn,c6Mn,c7Mn)           A^n_M(25)
               AEn = AN(xn,c1En,c2En,c3En,c4En,c5En,c6En,c7En)           A^n_E(25)

               TEF = sqrt(one+A_TE_C*Q2*exp(B_TE_C*Q2))

               Gp_M= AMp*G_D*mu_p*TEF                                    G^p_M
               Gp_E= AEp*G_D                                             G^p_E
               Gn_M= AMn*Gp_M*mu_n/mu_p*TEF                              G^n_M
               Gn_E= AEn*Gp_E*1.7*tn/(one+3.3*tn)                        G^n_E

               G_M = Gp_M-Gn_M                                           *TEF
               G_E = Gp_E-Gn_E
*              ------------------------------------------------------- *
      endSELECT
 2000    ReF_V= xi_V*(G_E+Q2*G_M/(4*mm_I))/(one+Q2/(4*mm_I))
         ReF_M= xi_M*(G_M-   G_E         )/(one+Q2/(4*mm_I))
         SELECTCASE(n_AP)                                                Switch for model of axial form factor in QES reactions
*              ------------------------------------------------------- *
               CASE(   1)                                                Standard dipole formula, Ref.[1]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*(MMA_ELS/(MMA_ELS+Q2))**2
               ImF_A= zero
*              ------------------------------------------------------- *
               CASE(   2)                                                Monopole formula, Ref.[9]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*MM_QES/(MM_QES+Q2)
               ImF_A= zero
*              ------------------------------------------------------- *
               CASE(   3)                                                L.M. Sehgal's Mmodification, Ref.[4,5]
*              ------------------------------------------------------- *
               ReF_A= xi_A*FCC_A0*
     #                exp(max(-Q2*(one+Q2/(4*mm_I)),-1.0d+02))/
     #                        (one+Q2/mm_A)
               ImF_A= zero
*              ------------------------------------------------------- *
      endSELECT
         tp  = Q2/(4*mm_p)
         tn  = Q2/(4*mm_n)
         Fp_1= (Gp_E+tp*Gp_M)/(one+tp)
         Fn_1= (Gn_E+tn*Gn_M)/(one+tp)
         Fp_2= (Gp_M-   Gp_E)/(one+tp)
         Fn_2= (Gn_M-   Gn_E)/(one+tn)
         FV_1= half*(Fp_1-Fn_1)
         FV_2= half*(Fp_2-Fn_2)
         SELECTCASE(n_RT_QES)
*              ------------------------------------------------------- *
               CASE(       1)                                            nu + p, an + p
*              ------------------------------------------------------- *
               Gs_M = ms_ELS/(one+Q2/MMV_QES)**2
               Gs_E = rs_ELS*tp/(one+Q2/MMV_QES)**2

               ReFsV= (Gs_E+tp*Gs_M)/(one+tp)
               ReFsM= (Gs_M-   Gs_E)/(one+tp)
               ReFsA= gs_ELS/(one+Q2/MMA_ELS)**2

               ReF_V= FV_1-2*sin2W*Fp_1-half*ReFsV
               ReF_M= FV_2-2*sin2W*Fp_2-half*ReFsM
               ReF_A=        half*ReF_A+half*ReFsA
*              ------------------------------------------------------- *
               CASE(       2)                                            nu + n, an + n
*              ------------------------------------------------------- *
               Gs_M = ms_ELS/(one+Q2/MMV_QES)**2
               Gs_E = rs_ELS*tn/(one+Q2/MMV_QES)**2

               ReFsV= (Gs_E+tn*Gs_M)/(one+tp)
               ReFsM= (Gs_M-   Gs_E)/(one+tp)
               ReFsA= gs_ELS/(one+Q2/MMA_ELS)**2

               ReF_V=-FV_1-2*sin2W*Fn_1-half*ReFsV
               ReF_M=-FV_2-2*sin2W*Fn_2-half*ReFsM
               ReF_A=       -half*ReF_A+half*ReFsA
*              ------------------------------------------------------- *
      endSELECT
         ReF_T= zero

         ImF_V= zero
         ImF_M= zero
         ImF_A= zero
         ImF_T= zero

c        IF (gs_ELS.ne.zero) print *, ' gs_ELS =', gs_ELS
c        IF (ms_ELS.ne.zero) print *, ' ms_ELS =', ms_ELS
c        IF (rs_ELS.ne.zero) print *, ' rs_ELS =', rs_ELS

         RETURN
*     ==================================================================

      END SUBROUTINE NucQESFF
