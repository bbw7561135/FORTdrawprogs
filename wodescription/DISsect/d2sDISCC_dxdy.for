************************************************************************
      FUNCTION d2sDISCC_dxdy(E_nu,x,y)
************************************************************************
*                                                                      *
*     This FUNCTION returns the charged current double differenti-     *
*     alneutrino-nucleon cross section  $d^2\sigam^CC_DIS/(dx dy)$     *
*     for  (anti)neutrino-nucleon  deep inelastic scattering. Here     *
*     Fi are parton structure functions of target nucleon,  Ai are     *
*     dynamic factors  and Ei are from the  $~(q^\mu q^\nu)/M^2_W$     *
*     part of the massive boson  propagator according to Ref. [1].     *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] S. Kretzer,  M.H. Reno.  "Tau  neutrino  deep  inelastic     *
*         charged  current  interactions,"  Phys. Rev. D 66 (2002)     *
*         113007 [arXiv: hep-ph/0208187].                              *
*     [2] K.S. Kuzmin, "Neutrino scattering off nucleons and pola-     *
*         rization of charged leptons in  quasielastic reactions,"     *
*         Ph.D. Thesis, JINR, Dubna, 2009/04/01  (Ph.D. Thesis ad-     *
*         visor V.A. Naumov, BLTP JINR).                               *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one,half,mm_W

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON       /n_LP/n_LP                                         Switch for lepton polarization type
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton

         Q2= 2*m_ini*x*y*E_nu
         CALL SFCC(E_nu,x,Q2,F1,F2,F3,F4,F5,F6)

         E1= 1                                                           1+(m_lep/m_W)**2*(1+((Q2/mm_W)/2))
         E2= 1                                                           1+(m_lep/m_W)**2*y*(1+y*(Q2+mm_lep)/(4*mm_W))*(1-y-(Q2+mm_lep)/(2*E_nu)**2)
         E4= 1                                                           1+(Q2/mm_W)*(1+((Q2/mm_W)/2))*2
         E5= 1                                                           1+(Q2/mm_W)+(y/2)*(1+(Q2/mm_W))*((m_lep/m_W)**2+(Q2/mm_W))

         a= mm_lep/(2*m_ini*E_nu)
         IF (n_LP.eq.0) THEN
           GOTO 1
                        ELSE
           E_lep= (one-y)*E_nu
           IF (E_lep**2-mm_lep.gt.zero) THEN
             P_lep= sqrt(E_lep**2-mm_lep)
                                        ELSE
             P_lep= E_lep*(one-mm_lep/(2*E_lep**2))
        endIF
           p= (x*y)+a*(one-(2*E_nu)/(E_lep+P_lep)       )
           m= (x*y)+a*(one-(2*E_nu)*(E_lep+P_lep)/mm_lep)
           IF (n_NT.eq.1) THEN
             IF (n_LP.eq.1) THEN; GOTO 2
                            ELSE; GOTO 3
          endIF
                          ELSE
             IF (n_LP.eq.1) THEN; GOTO 3
                            ELSE; GOTO 2
          endIF
        endIF
      endIF
*        ------------------------------------------------------------- *
*        NO POLARIZATION                                               *
*        ------------------------------------------------------------- *
    1    f = E_nu
         A1= y*((x*y)+a)
         A2= one-y-m_ini*(x*y)/(2*E_nu)-(m_lep/(2*E_nu))**2
         A3= y*(x*(one-y*half)-a*half)
         A4= a*((x*y)+a)
         A5=-a
         GOTO 4
*        ------------------------------------------------------------- *
*        "CORRECT" POLARIZATION                                        *
*        ------------------------------------------------------------- *
    2    f = E_nu*(E_lep+P_lep)/(2*P_lep)
         A1= p*y
         A2= (one-p*m_ini/(2*P_lep))*P_lep/E_nu
         A3= p*(E_nu+P_lep)/(2*E_nu)
         A4=-m*(mm_lep/ (E_lep+P_lep))**2/((2*E_nu)*m_ini)
         A5= m* mm_lep/((E_lep+P_lep)    * (2*E_nu)      )
         GOTO 4
*        ------------------------------------------------------------- *
*        "UNCORRECT" POLARIZATION                                      *
*        ------------------------------------------------------------- *
    3    f = E_nu*mm_lep/(2*(E_lep+P_lep)*P_lep)
         A1=-m*y
         A2= (one+m*m_ini/(2*P_lep))*P_lep/E_nu
         A3=-m*(E_nu -P_lep)   / (2*E_nu)
         A4= p*(E_lep+P_lep)**2/((2*E_nu)*m_ini)
         A5=-p*(E_lep+P_lep)   / (2*E_nu)
*        ------------------------------------------------------------- *

    4    d2sDISCC_dxdy=
     #   f*(F1*A1*E1+F2*A2*E2+n_NT*F3*A3+F4*A4*E4+F5*A5*E5)/
     #     (one+(Q2/mm_W))**2

         RETURN
      END FUNCTION d2sDISCC_dxdy