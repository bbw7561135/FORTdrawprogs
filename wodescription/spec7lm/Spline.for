cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     This file includes several routines for spline interpolation.    c
c                                                                      c
c     -----------------------------------------------------------      c
c    |                     by Vadim Naumov                       |     c
c    |                                                           |     c
c    |       Dipartimento di Fisica,  Universita di Ferrara      |     c
c    |                and INFN, Sezione di Ferrara               |     c
c    |             FORTRAN 77 version of May 10, 2001            |     c
c     -----------------------------------------------------------      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

************************************************************************
      SUBROUTINE BSPL2(NX,NY,F,D,Coef,Res)
************************************************************************

      INCLUDE 'Include.f'
      LOGICAL*2 Res
      DIMENSION F(NX,NY),D(NX+4,NY+4),Coef(*)
         COMMON /BSpl/Xmin,Ymin,StepX,StepY
      DO 1 j=1,NY
      J2=j+2
      DO 1 i=1,NX
      I2=i+2
    1 D(I2,J2)=3.90625d-3*LOG(F(i,j))
      J1=NY+1
      J2=J1+1   ! One can comment this line for some compilers (verify)
      J3=J2+1
      J4=J3+1
      DO 2 i=3,I2
      A=D(i,3)
      B=D(i,4)
      D(i,2)=Three*(A-B)+D(i,5)
      D(i,1)=Three*(D(i,2)-A)+B
      A=D(i,J1)
      B=D(i,J2)
      D(i,J3)=Three*(B-A)+D(i,NY)
    2 D(i,J4)=Three*(D(i,J3)-B)+A
      I1=NX+1
      I2=I1+1   ! One can comment this line for some compilers (verify)
      I3=I2+1
      I4=I3+1
      DO 3 j=1,J4
      A=D(3,j)
      B=D(4,j)
      D(2,j)=Three*(A-B)+D(5,j)
      D(1,j)=Three*(D(2,j)-A)+B
      A=D(I1,j)
      B=D(I2,j)
      D(I3,j)=Three*(B-A)+D(NX,j)
    3 D(I4,j)=Three*(D(I3,j)-B)+A
      DO 4 j=1,J2
      J3=j+1
      J4=j+2
      m=(j-1)*I2
      DO 4 i=1,I2
      I3=i+1
      I4=i+2
    4 Coef(m+i)= D(I4,j)+D(i,J4)+D(I4,J4)+D(I3,J3) *Cent
     -         -(D(I3,j)+D(i,J3)+D(I4,J3)+D(I3,J4))*Ten+D(i,j)
         IF (Res) THEN
         A=Zero
            DO j=1,NY
            Y=Ymin+StepY*(j-1)
               DO i=1,NX
               X=Xmin+StepX*(i-1)
               B=One-Spl2(Coef,X,Y,NX,NY)/F(i,j)
               IF(ABS(B).GT.ABS(A)) A=B
            endDO
         endDO
         WRITE(*,5) Cent*A
      endIF
      RETURN
    5 FORMAT('  Worst Residual is',F6.2,'%')

      END

************************************************************************
      FUNCTION Spl2(Coef,X,Y,NX,NY)
************************************************************************

      INCLUDE 'Include.f'
      SAVE A1,A2,A3,B1,B2,B3,M1,M2,M3
      DIMENSION Coef(*)
      COMMON /BSpl/Xmin,Ymin,StepX,StepY
      A3=(X-Xmin)/StepX
      B3=(Y-Ymin)/StepY
      M1=NINT(A3)
      M2=NINT(B3)
      IF (M1.LT.0.OR.M1.GE.NX.OR.M2.LT.0.OR.M2.GE.NY)
     .    WRITE(*,1) X,M1,Y,M2
      M3=NX+2
      A3=A3-M1
      B3=B3-M2
      M1=M2*M3+M1+1
      M2=M1+M3
      M3=M2+M3
      A2=A3**2+Quad
      B2=B3**2+Quad
      A1=A2-A3
      B1=B2-B3
      A3=A2+A3
      B3=B2+B3
      A2=Two-(A2+A2)
      B2=Two-(B2+B2)
       Spl2 = EXP((A1*Coef(M1)+A2*Coef(M1+1)+A3*Coef(M1+2))*B1
     +            +(A1*Coef(M2)+A2*Coef(M2+1)+A3*Coef(M2+2))*B2
     +            +(A1*Coef(M3)+A2*Coef(M3+1)+A3*Coef(M3+2))*B3)
      RETURN

*     ==================================================================
      ENTRY rSpl2(Coef)
*     ==================================================================

      rSpl2 = EXP((A1*Coef(M1)+A2*Coef(M1+1)+A3*Coef(M1+2))*B1
     +            +(A1*Coef(M2)+A2*Coef(M2+1)+A3*Coef(M2+2))*B2
     +            +(A1*Coef(M3)+A2*Coef(M3+1)+A3*Coef(M3+2))*B3)
      RETURN

    1 FORMAT(' Mistake in Spl2:'/'  X=',D22.15,' MX=',I4/
     ,                           '  Y=',D22.15,' MY=',I4)

      END

************************************************************************
      SUBROUTINE QSpl1(F,N,s,Yp1,Ypn,D)
************************************************************************

      INCLUDE 'Include.f'
      DIMENSION F(*),D(*),Ftmp(Nmax)

         IF (Yp1.GT.VeryBig) THEN
         D(1)=Zero
         Ftmp(1)=Zero
                             ELSE
         D(1)=-Half
         Ftmp(1)=(Three/s)*((F(2)-F(1))/s-Yp1)
      endIF
      a=Six/s**2
      DO 1 i=2,N-1
      Q=D(i-1)+Four
      D(i)=-One/Q
    1 Ftmp(i)=(a*(F(i+1)-Two*F(i)+F(i-1))-Ftmp(i-1))/Q
         IF (Ypn.GT.VeryBig) THEN
         D(N)=Zero
                             ELSE
         D(N)=(a*(Ypn*s-F(N)+F(N-1))-Ftmp(N-1))/(D(N-1)+Two)
      endIF
      DO 2 k=N-1,1,-1
    2 D(k)=D(k)*D(k+1)+Ftmp(k)
      RETURN

*     ==================================================================
      ENTRY SplQ1N(F,N,s,D)                         ! "Natural" spline
*     ==================================================================

      D(1)=Zero
      Ftmp(1)=Zero
      a=Six/s**2
      DO 3 i=2,N-1
      Q=D(i-1)+Four
      D(i)=-One/Q
    3 Ftmp(i)=(a*(F(i+1)-Two*F(i)+F(i-1))-Ftmp(i-1))/Q
      D(N)=Zero
      DO 4 k=N-1,1,-1
    4 D(k)=D(k)*D(k+1)+Ftmp(k)
      RETURN

      END

************************************************************************
      FUNCTION SplQ1(delta,s,F,D)
************************************************************************

      INCLUDE 'Include.f'
      SAVE
      DIMENSION D(*),F(*)

      c=delta/s
      k1=INT(c)+1
      k2=k1+1
      c=c+One
      a=DFLOAT(k2)-c
      b=c-DFLOAT(k1)
      c=s**2*OS
      fa=a*(a**2*c-c)
      fb=b*(b**2*c-c)
      SplQ1=a*F(k1)+b*F(k2)+fa*D(k1)+fb*D(k2)
      RETURN

*     ==================================================================
      ENTRY rSplQ2(F,D)
*     ==================================================================

      rSplQ2=a*F(k1)+b*F(k2)+fa*D(k1)+fb*D(k2)
      RETURN

      END

************************************************************************
      SUBROUTINE QSpl2(F,D,M,N,Res)
************************************************************************

      INCLUDE 'Include.f'
      LOGICAL*2 Res
      DIMENSION D(M,N),F(M,N),Dtmp(Nmax),Ftmp(Nmax)
      COMMON /QSpl/Xmin,Ymin,stepX,stepY
      DO j=1,M
        DO 1 k=1,N
    1   Ftmp(k)=F(j,k)
      CALL SplQ1N(Ftmp,N,stepY,Dtmp)
        DO 2 k=1,N
    2   D(j,k)=Dtmp(k)
      endDO
      IF(.NOT.Res) RETURN
c     Test calculation of the maximum discrepancy in reference points
      Discr_max=Zero
      i=0
      X0=Xmin-stepX
      Y0=Ymin-stepY
      DO 3 j=1,M
      X=X0+stepX*j
      DO 3 k=1,N
      Y=Y0+stepY*k
      z=F(j,k)
         IF(z.EQ.Zero) THEN
         i=i+1
         Discr=        SplQ2(F,D,M,N,X,Y)
                       ELSE
         Discr=ABS(One-SplQ2(F,D,M,N,X,Y)/z)
      endIF
    3 IF (Discr.GT.Discr_max) Discr_max=Discr
      PRINT 4, Discr_max
      IF (i.GT.0) PRINT *,' Function vanishs in',i,' reference points'
      RETURN
    4 FORMAT(' Maximum discrepancy in reference points =',1pD10.3)

      END

************************************************************************
      FUNCTION SplQ2(F,D,M,N,X,Y)
************************************************************************

      INCLUDE 'Include.f'
      DIMENSION D(M,N),F(M,N),Dtmp(Nmax),Ftmp(Nmax),Ztmp(Nmax)
      COMMON /QSpl/Xmin,Ymin,stepX,stepY
      DO 2 j=1,M
        DO 1 k=1,N
        Ftmp(k)=F(j,k)
    1   Dtmp(k)=D(j,k)
    2   Ztmp(j)=SplQ1(Y-Ymin,stepY,Ftmp,Dtmp)
           CALL SplQ1N(Ztmp,M,stepX,Dtmp)
          SplQ2=SplQ1(X-Xmin,stepX,Ztmp,Dtmp)
      RETURN

      END
