************************************************************************
      FUNCTION FactorPauli(A,Z,Q2)
************************************************************************
*                                                                      *
*     For deuterium the  parametrization of Factor Pauli are taken     *
*     from  paper by  S.K. Singh, Nucl. Phys. B 36 (1972) 419-435.     *
*                                                                      *
************************************************************************

         IMPLICIT REAL(8) (A-Z)

         REAL(8),PARAMETER::
     #           A1=0.38455,
     #           A2=0.18452,
     #           t1=0.03648,
     #           t2=0.00920

         COMMON      /m_ini/m_ini,mm_ini                                 mass of a initial nucleon  [GeV, GeV^2]
         COMMON    /P_Fermi/P_Fermi                                      Fermi momentum

*        ------------------------------------------------------------- *
         IF (A.eq.2.0d+00 .and. Z.eq.1.0d+00) THEN                       DEUTERIUM 
*        ------------------------------------------------------------- *
           FactorPauli=1-A1*exp(-Q2/t1)-A2*exp(-Q2/t2)
*        ------------------------------------------------------------- *
                                              ELSE
*        ------------------------------------------------------------- *
           Q2=max(0.0d+00,Q2)
           QV=sqrt(Q2**2+4*mm_ini*Q2)/(2*m_ini)
           xf=QV/(2*P_Fermi)
           u =(2*(A-Z)/A)**(1.0/3)
           v =(2*   Z /A)**(1.0/3)
           IF ((2*xf).le.(u-v)) THEN 
             D=Z
       ELSEIF ((2*xf).gt.(u-v) .and. (2*xf).le.(u+v)) THEN
             D=0.5*A*(1-0.75*xf*(u**2+v**2)+(xf**3)*0.5-
     #         3*((u**2-v**2)**2)/(32*xf))
       ELSEIF ((2*xf).gt.(u+v)) THEN
             D=0.0d+00
        endIF
           FactorPauli=1-D/(A-Z)
*        ------------------------------------------------------------- *
      endIF

         RETURN
      END FUNCTION FactorPauli

************************************************************************
      FUNCTION FactorPauli_D2(Q2)
************************************************************************

         IMPLICIT REAL(8) (A-Z)

         REAL(8),PARAMETER::
     #           A1=0.38455,
     #           A2=0.18452,
     #           t1=0.03648,
     #           t2=0.00920

         FactorPauli_D2=1-A1*exp(-Q2/t1)-A2*exp(-Q2/t2)

         RETURN
      END FUNCTION FactorPauli_D2