************************************************************************
      FUNCTION MA_QES_EFF(E_nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL,PARAMETER::
     #                MA0  = 9.64d-01,
     #                E0   = 3.56d-01,
     #                a0   = 0.846d+00,
     #                a1   = 0.10583186025394,
     #                b1   = 0.36563826143285,
     #                c1   = 0.91005905942944,
     #                a3   = 0.00245436835458,
     #                b3   = 0.45861634912994,
     #                c3   = 1.03324251687687
                 REAL
     #                MA_cen(5)/0.90,0.964,1.10,1.20,1.35/

         COMMON       /n_MA/n_MA                                         Switch for MA0
         COMMON        /n_b/n_b
         COMMON     /CoefMA/CoefMA
         COMMON     /n_corv/n_corv

               IF (n_corv.EQ.2) THEN
                 SELECTCASE(n_b)
                       CASE(1)
                             MA_QES_EFF=c1+b1/(E_nu+a1)
                       CASE(2)
                             MA_QES_EFF=MA0*(1+(E0/E_nu)**a0)
                       CASE(3)
                             MA_QES_EFF=c3+b3/(E_nu+a3)
              endSELECT
                           ELSE
                 MA_QES_EFF = CoefMA*MA_cen(n_MA)
            endIF
        
         RETURN

      END FUNCTION MA_QES_EFF