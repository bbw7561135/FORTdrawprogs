************************************************************************
      FUNCTION MA_QES_EFF(E_nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL,PARAMETER::
!     #                MA0= 9.780d-01,
     #                E0 = 4.590d-01,
     #                a  = 1.283d+00

         COMMON       /n_MA/n_MA                                         Switch for MA0

         SELECTCASE(n_MA)
               CASE(1)
                      MA0=0.90d+00
               CASE(2)
                      !MA0=0.95d+00
                      MA0=0.978d+00
               CASE(3)
                      MA0=1.10d+00
               CASE(4)
                      MA0=1.20d+00
               CASE(5)
                      MA0=1.35d+00
      endSELECT

!         MA_QES_EFF=MA0
         MA_QES_EFF=MA0*(1+(E0/E_nu)**a)
!         MA_QES_EFF=0.93*MA0*(1+(E0/E_nu)**a)
!         MA_QES_EFF=1.07*MA0*(1+(E0/E_nu)**a)
        
         RETURN

      END FUNCTION MA_QES_EFF