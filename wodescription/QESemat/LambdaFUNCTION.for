************************************************************************
      FUNCTION LambdaFUNCTION(a,b,c)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         LambdaFUNCTION=a**2+b**2+c**2-2*(a*b+b*c+a*c)

         RETURN
      END FUNCTION LambdaFUNCTION