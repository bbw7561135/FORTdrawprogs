************************************************************************
      FUNCTION d2sDISNC_dxdy(E_nu,x,y)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one,half,mm_Z

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon

         Q2=(2*E_nu)*(m_ini*(x*y))

         CALL SFNC(E_nu,x,Q2,F1,F2,F3)

         A1= (x*y)*y
         A2= one-y-(m_ini*(x*y))/(2*E_nu)
         A3= (x*y)*(one-y*half)

         d2sDISNC_dxdy= E_nu*(F1*A1+F2*A2+n_NT*F3*A3)/(one+Q2/mm_Z)**2

         RETURN
      END FUNCTION d2sDISNC_dxdy