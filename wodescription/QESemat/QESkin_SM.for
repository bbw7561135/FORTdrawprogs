************************************************************************
      FUNCTION Q2lim1_SM(Q2)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: half

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON       /E_nu/E_nu                                         Neutrino energy
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy

         s     = 2*E_nu*m_tar+mm_tar
         nu_max= (s-mm_tar-mm_lep*(s-mm_tar)/
     #           (Q2+mm_lep)-mm_tar*(Q2+mm_lep)/(s-mm_tar))/(2*m_tar)
         E     = sqrt(P_FeMAX**2+mm_ini)
         b     = (E-E_nuBIN)**2-P_FeMAX**2
         a     = (Q2+mm_fin-b)*half
         nu_1  = (a*(E-E_nuBIN)-P_FeMAX*sqrt(a**2+Q2*b))/b

         Q2lim1_SM=nu_1-nu_max

         RETURN
      END FUNCTION Q2lim1_SM

************************************************************************
      FUNCTION Q2lim2_SM(Q2)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: half

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy

         nu_min= ((m_rnu+m_fin)**2+Q2-mm_tar)/(2*m_tar)
         E     = sqrt(P_FeMAX**2+mm_ini)
         b     = (E-E_nuBIN)**2-P_FeMAX**2
         a     = (Q2+mm_fin-b)*half
         nu_2  = (a*(E-E_nuBIN)+P_FeMAX*sqrt(a**2+Q2*b))/b

         Q2lim2_SM=nu_min-nu_2

         RETURN
      END FUNCTION Q2lim2_SM

************************************************************************
      SUBROUTINE nuQES_SM_lim(E_nu,Q2,nu_min,nu_max)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: half

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy

         s     = 2*E_nu*m_tar+mm_tar
         nu_min= ((m_rnu+m_fin)**2+Q2-mm_tar)/(2*m_tar)
         nu_max= (s-mm_tar-mm_lep*(s-mm_tar)/(Q2+mm_lep)-
     #           mm_tar*(Q2+mm_lep)/(s-mm_tar))/(2*m_tar)
         E     = sqrt(P_FeMAX**2+mm_ini)
         b     = (E-E_nuBIN)**2-P_FeMAX**2
         a     = (Q2+mm_fin-b)*half
         tmp1  = a*(E-E_nuBIN)
         tmp2  = P_FeMAX*sqrt(a**2+Q2*b)
         nu_1  = (tmp1-tmp2)/b
         nu_2  = (tmp1+tmp2)/b
         nu_min= max(nu_min,nu_1)
         nu_max= min(nu_max,nu_2)

         RETURN
      END SUBROUTINE nuQES_SM_lim

************************************************************************
      SUBROUTINE E_nu_th_SM(E_min)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero

         IMPLICIT REAL (A-Z)

         INTEGER,PARAMETER::
     #           MFC= 1000
            REAL,PARAMETER::
     #           EPS= 1.0d-08

         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus

         E_min= ((m_lep+m_rnu+m_fin)**2-mm_tar)/(2*m_tar)
         Enu_2= 5.0d+00                                                  ???
         IF (QEL_EnuMin_SM(E_min).gt.zero) THEN
           Enu_rf= DZEROX(E_min,Enu_2,EPS,MFC,QEL_EnuMin_SM,1)
                                           ELSE
           Enu_rf=-1.0d+01                                               ???
      endIF
         E_min= max(E_min,Enu_rf)
         IF (E_min.lt.zero) THEN
           E_min= zero
           print *, 'E_min =', E_min
      endIF

         RETURN
      END SUBROUTINE E_nu_th_SM

************************************************************************
      FUNCTION QEL_EnuMin_SM(E_nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one,half,Precision

         IMPLICIT REAL (A-Z)

               LOGICAL(2) LLM
            
            INTEGER,PARAMETER::
     #              MFC  = 1000
               REAL,PARAMETER::
     #              EPS  = 1.0d-06,
     #              Delta= 1.0d-14

         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus
         COMMON       /E_nu/Enu_in                                       Neutrino energy

         EXTERNAL Q2lim1_SM

         Enu_in = E_nu
         s      = 2*E_nu*m_tar+mm_tar
         W2     = (m_rnu+m_fin)**2
         c      = half*(W2+mm_lep-mm_tar*(W2-mm_lep)/s)
         sqrtD  = sqrt(max(Precision,LambdaFUNCTION(one,mm_lep/s,W2/s)))
         tmp    = half*(s-mm_tar)
         Q2_lim1= tmp*(one-sqrtD)-c
         Q2_lim2= tmp*(one+sqrtD)-c
      
         CALL
     #   DMINFC(Q2lim1_SM,Q2_lim1,Q2_lim2,EPS,Delta,Q2_0,F_MIN,LLM)

         QEL_EnuMin_SM=F_MIN

         RETURN
      END FUNCTION QEL_EnuMin_SM

************************************************************************
      SUBROUTINE Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one,half

         IMPLICIT REAL (A-Z)

         LOGICAL(2) LLM

         INTEGER,PARAMETER::
     #           MFC  = 1000
            REAL,PARAMETER::
     #           Eps  = 1.0d-08,
     #           Delta= 1.0d-14

         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus

         EXTERNAL Q2lim1_SM,Q2lim2_SM

         s      = 2*E_nu*m_tar+mm_tar
         W2     = (m_rnu+m_fin)**2
         c      = half*(W2+mm_lep-mm_tar*(W2-mm_lep)/s)
         sqrtD  = sqrt(LambdaFUNCTION(one,mm_lep/s,W2/s))
         tmp    = half*(s-mm_tar)
         Q2_min = tmp*(one-sqrtD)-c
         Q2_max = tmp*(one+sqrtD)-c

         CALL DMINFC(Q2lim1_SM,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM)

         IF (F_MIN.gt.zero) THEN
           PRINT 100, E_nu, Q2_min, Q2_max, Q2_0, F_MIN
           STOP 'STOP. SUBROUTINE Q2QES_SM_lim'
      endIF
         IF (Q2lim1_SM(Q2_min).gt.zero) THEN
           Q2_RF = DZEROX(Q2_min,Q2_0,EPS,MFC,Q2lim1_SM,1)
           Q2_min= max(Q2_min,Q2_RF)
      endIF
         IF (Q2lim1_SM(Q2_max).gt.zero) THEN
           Q2_RF = DZEROX(Q2_0,Q2_max,Eps,MFC,Q2lim1_SM,1)
           Q2_max= min(Q2_max,Q2_RF)
      endIF
         IF (Q2lim2_SM(Q2_min).gt.zero) THEN
           IF (Q2lim2_SM(Q2_max).gt.zero) THEN
             PRINT 200, E_nu
             STOP 'STOP. SUBROUTINE Q2QES_SM_lim'
                                          ELSE
             Q2_RF = DZEROX(Q2_min,Q2_max,Eps,MFC,Q2lim2_SM,1)
             Q2_min= max(Q2_min,Q2_RF)
        endIF
      endIF
         Q2_min= max(Q2_min,zero)

         RETURN
  100 FORMAT('Error 1: no overlapped area for energy ',E12.5
     #     ,/,'Q2_min = ',E12.5,' Q2_max = ',E12.5
     #     ,/,'Q2_0   = ',E12.5,' F_min  = ',E12.5)
  200 FORMAT('Error 2: no overlapped area for energy ',E12.5)
      END SUBROUTINE Q2QES_SM_lim

************************************************************************
      SUBROUTINE kFQES_SM_lim(Q2,nu,kF_min,kF_max)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy

         qv    = sqrt(nu**2+Q2)
         c_f   = (nu-E_nuBIN)/qv
         d_f   = (E_nuBIN**2-2*nu*E_nuBIN-Q2+mm_ini-mm_fin)/(2*qv*m_ini)
         Ef_min= m_ini*(c_f*d_f+sqrt(one-c_f**2+d_f**2))/(one-c_f**2)
         kF_min= sqrt(Ef_min**2-mm_ini)
         kF_max= P_FeMAX

         RETURN
      END SUBROUTINE kFQES_SM_lim