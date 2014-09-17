************************************************************************
      FUNCTION FunGeM(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE Routines

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         COMMON            /x/x                                          Bjorken scaling variable x
         COMMON           /Q2/Q2                                         Square of mometum transfer (Q^2=-q^2)
         COMMON         /E_nu/E_nu                                       Neutrino energy

         FunGeM= one
         RETURN

*     ================================================================ *
      ENTRY GeM_FQCD_L(var)
*     ================================================================ *
         CALL SFCC(E_nu,var,Q2,F1,F2,F3,F4,F5,F6)
         CALL PDF (E_nu,var,Q2,Uq,Ua,Dq,Da,Sq,Sa,Cq,Ca,G,A)
         GeM_FQCD_L= A*(8*F2/3+16*(var-x)*G)/var
         RETURN

      END FUNCTION FunGeM