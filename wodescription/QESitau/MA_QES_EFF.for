************************************************************************
      FUNCTION MA_QES_EFF(E_nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL,PARAMETER::
     #                MA02 = 9.73d-01,
!     #                MA02 = 9.64d-01,
     #                E0  = 3.56d-01,
     #                a   = 0.846d+00,
     #                ErrC= 1.0

         COMMON       /n_MA/n_MA                                         Switch for MA0

         SELECTCASE(n_MA)
               CASE(1)
                      MA0=0.90d+00
               CASE(2)
                      MA0=MA02
               CASE(3)
                      MA0=1.10d+00
               CASE(4)
                      MA0=1.20d+00
               CASE(5)
                      MA0=1.35d+00
      endSELECT

         MA_QES_EFF=ErrC*MA0
!         MA_QES_EFF=ErrC*MA0*(1+(E0/E_nu)**a)
        
         RETURN

      END FUNCTION MA_QES_EFF