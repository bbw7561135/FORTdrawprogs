************************************************************************
      SUBROUTINE UTest(U,DeltaI)
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2012/11/15  *
************************************************************************
*                                                                      *
*     This SUBROUTINE tests unitarity of U-matrix                      *
*                                                                      *
************************************************************************

         USE OscMatParameters
         
         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

                 REAL
     #                Imtrx(3,3),DeltaI(3,3)
              COMPLEX
     #                U(3,3)
     
         Imtrx=matmul(transpose(conjg(U)),U)
         
         DeltaI(1,1)=Imtrx(1,1)-one
         DeltaI(1,2)=Imtrx(1,2)
         DeltaI(1,3)=Imtrx(1,3)

         DeltaI(2,1)=Imtrx(2,1)
         DeltaI(2,2)=Imtrx(2,2)-one
         DeltaI(2,3)=Imtrx(2,3)

         DeltaI(3,1)=Imtrx(3,1)
         DeltaI(3,2)=Imtrx(3,2)
         DeltaI(3,3)=Imtrx(3,3)-one
     
      RETURN
         
      END SUBROUTINE UTest
