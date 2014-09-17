************************************************************************
      SUBROUTINE yQES_lim(E_nu,y_min,y_max)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system

         s      = mm_ini+2*m_ini*E_nu
         sqrt_s = sqrt(s)
         Ecm_lep= (s+mm_lep-mm_fin)/(2*sqrt_s)
         Pcm_lep= sqrt(Ecm_lep**2-mm_lep)
         y_min  = one-(Ecm_lep*(one+m_ini/E_nu)+Pcm_lep)/sqrt_s
         y_max  = one-(Ecm_lep*(one+m_ini/E_nu)-Pcm_lep)/sqrt_s

         RETURN
      END SUBROUTINE yQES_lim

************************************************************************
      SUBROUTINE Q2QES_lim(E_nu,Q2_min,Q2_max)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one,Precision

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system

         s      = mm_ini+2*m_ini*E_nu
         sqrt_s = sqrt(s)
         Ecm_lep= (s+mm_lep-mm_fin)/(2*sqrt_s)
         Pcm_lep= sqrt(Ecm_lep**2-mm_lep)
         y_min  = one-(Ecm_lep*(one+m_ini/E_nu)+Pcm_lep)/sqrt_s
         y_max  = one-(Ecm_lep*(one+m_ini/E_nu)-Pcm_lep)/sqrt_s
         Q2_min = 2*m_ini*E_nu*y_min+mm_ini-mm_fin
         Q2_max = 2*m_ini*E_nu*y_max+mm_ini-mm_fin

         RETURN
      END SUBROUTINE Q2QES_lim

************************************************************************
      FUNCTION zetaQES(E_nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system

         s= m_ini*(2*E_nu+m_ini)
         zetaQES= (m_ini*sqrt((s+mm_lep-mm_fin)**2-4*mm_lep*s))/
     #            (m_lep*(s-mm_ini))

         RETURN
      END FUNCTION zetaQES