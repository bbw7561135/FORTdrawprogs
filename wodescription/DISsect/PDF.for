************************************************************************
      SUBROUTINE PDFsetup(E_nu,x,Q2,Uq,Ua,Dq,Da,Sq,Sa,Cq,Ca)
************************************************************************
*                                                                      *
*                                                                      *
*                                    ITEP, Moscow, Russia, 2005/10/21  *
*                                     JINR, Dubna, Russia, 2006/03/03  *
************************************************************************

         USE PhysMathConstants, ONLY: zero

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         REAL,PARAMETER:: alphaCTEQ6_s=0.114

         COMMON       /n_TT/n_TT                                         Switch for nuclear target type
         COMMON   /n_Qc_DIS/n_Qc_DIS                                     Switch for type of PDF limitation
         COMMON     /PDFLIB/NGROUP,NSET                                  Parameters for PDFLIB setup
         COMMON     /Q2_DIS/Q2_DIS                                       Minimal Q^2_DIS value for PFDs

         DIMENSION pdfs(0:9)
c        DIMENSION dpdfs(0:9,23),delpdf(0:9)

         SELECTCASE(NGROUP)
               CASE(     3)
               SELECTCASE(NSET)
                     CASE(  99); CALL PDFlib8_04set
            endSELECT
               CASE(     5); CALL PDFlib8_04set
               CASE(     6); CALL SETCTQ6(NSET)
      endSELECT
         RETURN

*     ==================================================================
      ENTRY PDF(E_nu,x,Q2,Uq,Ua,Dq,Da,Sq,Sa,Cq,Ca,G,A)
*     ==================================================================
         E_nu_SAVE= E_nu
         x_SAVE   = x
         Q2_save  = Q2

         IF (n_Qc_DIS.eq.0 .and. Q2.lt.Q2_DIS) Q2=Q2_DIS

         CALL xDIS_lim(E_nu,x_min,x_max)
         IF (x.le.x_min .or. x.ge.x_max) THEN
           Uq= zero; Ua= zero; Dq= zero; Da= zero
           Sq= zero; Sa= zero; Cq= zero; Ca= zero; G = zero; A = zero
                                         ELSE
           SELECTCASE(NGROUP)
*                ----------------------------------------------------- *
                 CASE(     3)                                            MRST
*                ----------------------------------------------------- *
                 SELECTCASE(NSET)
                       CASE(   0)                                        MRST 2004 QED /02.11.2004/
                       SELECTCASE(n_TT)
*                            - - - - - - - - - - - - - - - - - - - - - *
                             CASE(   1)
*                            - - - - - - - - - - - - - - - - - - - - - *
                             CALL MRST1(x,Q2,Uv,Dv,Us,Ds,S,C,B,G,P)
*                            - - - - - - - - - - - - - - - - - - - - - *
                             CASE(   2)
*                            - - - - - - - - - - - - - - - - - - - - - *
                             CALL MRST2(x,Q2,Dv,Uv,Ds,Us,S,C,B,G,P)
*                            - - - - - - - - - - - - - - - - - - - - - *
                    endSELECT
*                      - - - - - - - - - - - - - - - - - - - - - - - - *
                       CASE(  99)                                        MRST 1998
*                      - - - - - - - - - - - - - - - - - - - - - - - - *
                       CALL STRUCTM(x,sqrt(Q2),Uv,Dv,Us,Ds,S,C,B,T,G)
*                      - - - - - - - - - - - - - - - - - - - - - - - - *
              endSELECT
                 IF (Uv.lt.zero) Uv= zero
                 IF (Dv.lt.zero) Dv= zero
                 IF (Us.lt.zero) Us= zero
                 IF (Ds.lt.zero) Ds= zero
                 IF (S .lt.zero) S = zero
                 IF (C .lt.zero) C = zero
                 IF (B .lt.zero) B = zero

                 Ua= Us/x; Uq= Uv/x+Us/x
                 Da= Ds/x; Dq= Dv/x+Ds/x
                 Sq= S /x; Sa= Sq
                 Cq= C /x; Ca= Cq
                 G = G /x; A = alphaCTEQ6_s
*                ----------------------------------------------------- *
                 CASE(     5)                                            GRV98
*                ----------------------------------------------------- *
                 CALL STRUCTM(x,sqrt(Q2),Uv,Dv,Us,Ds,S,C,B,T,G)
                 Ua= Us/x; Uq= Uv/x+Us/x
                 Da= Ds/x; Dq= Dv/x+Ds/x
                 Sq= S /x; Sa= Sq
                 Cq= C /x; Ca= Cq
                 G = G /x; A = ALPHAS2(Q)
*                ----------------------------------------------------- *
                 CASE(     6)                                            CTEQ6/12.12.2004/, CTEQ6.5M/12.12.2006/
*                ----------------------------------------------------- *
                 Uq= CTQ6PDF( 1,x,(sqrt(Q2)))
                 Ua= CTQ6PDF(-1,x,(sqrt(Q2)))
                 Dq= CTQ6PDF( 2,x,(sqrt(Q2)))
                 Da= CTQ6PDF(-2,x,(sqrt(Q2)))
                 Sq= CTQ6PDF( 3,x,(sqrt(Q2)))
                 Sa= CTQ6PDF(-3,x,(sqrt(Q2)))
                 Cq= CTQ6PDF( 4,x,(sqrt(Q2)))
                 Ca= CTQ6PDF(-4,x,(sqrt(Q2)))
                 G = CTQ6PDF( 0,x,(sqrt(Q2)))
                 A = alphaCTEQ6_s
*                ----------------------------------------------------- *
                 CASE(     7)                                            ALEKHIN'2006
*                ----------------------------------------------------- *
c                CALL ALEKHIN06(x,Q2,    pdfs,dpdfs,npdf,npar)
                 CALL ALEKHIN06(x,Q2,pdfs,          npdf,npar)

                 Uq= pdfs(1)/x; Ua= pdfs(4)/x
                 Dq= pdfs(2)/x; Da= pdfs(6)/x
                 Sq= pdfs(5)/x; Sa= Sq
                 Cq= pdfs(7)/x; Ca= Cq
                 Bq= pdfs(8)/x; Ba= Bq
                 Tq= pdfs(9)/x; Ta= Tq
                 G = pdfs(3)/x; A = pdfs(0)/x
*                ----------------------------------------------------- *
        endSELECT
      endIF

         E_nu= E_nu_SAVE
         x   = x_SAVE
         Q2  = Q2_SAVE

         RETURN
*     ==================================================================

      END SUBROUTINE PDFsetup