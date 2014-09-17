************************************************************************
*                                                                      *
*         A PACKAGE FOR THE MRST 2004-QED PARTON DISTRIBUTIONS         *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] A.D. Martin, R.G. Roberts, W.J. Stirling, and R.S. Thor-     *
*         ne, "Parton  distributions incorporating QED contributi-     *
*         ons,"  Eur.  Phys. J. C39 (2005) 155-161  [arXiv: hep-ph     *
*         /0411040].                                                    *
*                                                                      *
*     There are 2 pdf sets: proton distributions (mode=1), neutron     *
*     distributions (mode=2).  Note the  extra argument "phot" for     *
*     the photon  distribution in  the proton/neutron.  As always,     *
*     the quantity returned is x (f,x,Q^2). This  SUBROUTINE  uses     *
*     an improved interpolation procedure for extracting values of     *
*     the PDF's from the grid.                                         *
*                                                                      *
************************************************************************
      SUBROUTINE MRST04QED(x,q,mode,
     #                     upv,dnv,usea,dsea,str,chm,bot,glu,phot)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         DATA xmin,xmax,qsqmin,qsqmax/1.0d-05,1.0d+00,1.25d+00,1.0d+07/

          q2=q*q
          IF (q2.lt.qsqmin .or. q2.gt.qsqmax) PRINT 99, q2
          IF (x .lt.xmin   .or.  x.gt.  xmax) PRINT 98, x
          IF (mode.eq.1) THEN
            CALL mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
      ELSEIF (mode.eq.2) THEN
            CALL mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
       endIF

         RETURN
   99 FORMAT(' WARNING: Q^2 VALUE IS OUT OF RANGE ','q2= ',e10.5)
   98 FORMAT(' WARNING: X   VALUE IS OUT OF RANGE ','x = ',e10.5)
      END SUBROUTINE MRST04QED

************************************************************************
      SUBROUTINE mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE InpOutUnits

         IMPLICIT REAL (A-H,O-Z)

         PARAMETER(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)

         REAL f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     #        f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc),
     #        fb(nx,nqb)
         REAL qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     #        cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),
     #        cc8(nx,nq,4,4),cc9(nx,nq,4,4),ccc(nx,nqc,4,4),
     #        ccb(nx,nqb,4,4)
         REAL xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)

         DATA xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     #           1d-4,2d-4,4d-4,6d-4,8d-4,
     #           1d-3,2d-3,4d-3,6d-3,8d-3,
     #           1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     #           .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     #           .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     #           .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     #           .8d0,.9d0,1d0/
         DATA qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     #           1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     #           1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     #           1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     #           1.8d6,3.2d6,5.6d6,1d7/
c        DATA xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/               KK: This variables have not been used
         DATA init/0/

         SAVE

         xsave =x
         q2save=qsq
         IF (init.ne.0) GOTO 10
         OPEN(unit=33,file=IniPP//'PDFLIB-like/MRST/'//                  KK
     #   'MRST 2004 QED/MRST-02.11.2004-p.dat',status='old')             KK
         DO 20 n=1,nx-1
         DO 20 m=1,nq
           READ(33,50) f1(n,m), ! uval
     #                 f2(n,m), ! val
     #                 f3(n,m), ! glue
     #                 f4(n,m), ! usea
     #                 f5(n,m), ! chm
     #                 f7(n,m), ! str
     #                 f6(n,m), ! btm
     #                 f8(n,m), ! dsea
     #                 f9(n,m)  ! photon
   20 conTINUE
         DO 40 m=1,nq
           f1(nx,m)=0.0d+00
           f2(nx,m)=0.0d+00
           f3(nx,m)=0.0d+00
           f4(nx,m)=0.0d+00
           f5(nx,m)=0.0d+00
           f6(nx,m)=0.0d+00
           f7(nx,m)=0.0d+00
           f8(nx,m)=0.0d+00
           f9(nx,m)=0.0d+00
   40 conTINUE
         DO n=1,nx
           xxl(n)=log(xx(n))
      endDO
         DO m=1,nq
           qql(m)=log(qq(m))
      endDO
         CALL jeppe1(nx,nq,xxl,qql,f1,cc1)
         CALL jeppe1(nx,nq,xxl,qql,f2,cc2)
         CALL jeppe1(nx,nq,xxl,qql,f3,cc3)
         CALL jeppe1(nx,nq,xxl,qql,f4,cc4)
         CALL jeppe1(nx,nq,xxl,qql,f6,cc6)
         CALL jeppe1(nx,nq,xxl,qql,f8,cc8)
         CALL jeppe1(nx,nq,xxl,qql,f9,cc9)
         emc2= 2.045
         emb2=18.500
         DO 44 m=1,nqc
           qqlc(m)=qql(m+nqc0)
         DO 44 n=1,nx
           fc(n,m)=f5(n,m+nqc0)
   44 conTINUE
         qqlc(1)=log(emc2)
         CALL jeppe1(nx,nqc,xxl,qqlc,fc,ccc)
         DO 45 m=1,nqb
           qqlb(m)=qql(m+nqb0)
         DO 45 n=1,nx
           fb(n,m)=f7(n,m+nqb0)
   45 conTINUE
         qqlb(1)=log(emb2)
         CALL jeppe1(nx,nqb,xxl,qqlb,fb,ccb)
         init=1

   10    CONTINUE

         xlog  =log(x)
         qsqlog=log(qsq)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)
         chm=0.0d+00
         IF (qsq.gt.emc2) THEN
           CALL jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endIF
         bot=0.0d+00
         if (qsq.gt.emb2) THEN
           CALL jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endIF
         x  =xsave
         qsq=q2save

         RETURN
   50 FORMAT(9f10.5)
      END SUBROUTINE mrst1

************************************************************************
      SUBROUTINE mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,phot)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE InpOutUnits

         IMPLICIT REAL (A-H,O-Z)

         PARAMETER(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)

         REAL f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     #        f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc),
     #        fb(nx,nqb)
         REAL qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     #        cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),
     #        cc8(nx,nq,4,4),cc9(nx,nq,4,4),ccc(nx,nqc,4,4),
     #        ccb(nx,nqb,4,4)
         REAL xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)

         DATA xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     #           1d-4,2d-4,4d-4,6d-4,8d-4,
     #           1d-3,2d-3,4d-3,6d-3,8d-3,
     #           1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     #           .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     #           .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     #           .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     #           .8d0,.9d0,1d0/
         DATA qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     #           1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     #           1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     #           1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     #           1.8d6,3.2d6,5.6d6,1d7/
c        DATA xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/               KK: This variables have not been used
         DATA init/0/

         SAVE

         xsave =x
         q2save=qsq
         IF (init.ne.0) GOTO 10
         OPEN(unit=34,file=IniPP//'PDFLIB-like/MRST/'//                  KK
     #   'MRST 2004 QED/MRST-02.11.2004-n.dat',status='old')             KK
         DO 20 n=1,nx-1
         DO 20 m=1,nq
         READ(34,50) f1(n,m), ! uval
     #               f2(n,m), ! val
     #               f3(n,m), ! glue
     #               f4(n,m), ! usea
     #               f5(n,m), ! chm
     #               f7(n,m), ! str
     #               f6(n,m), ! btm
     #               f8(n,m), ! dsea
     #               f9(n,m)  ! photon
   20 conTINUE
         DO 40 m=1,nq
           f1(nx,m)=0.0d+00
           f2(nx,m)=0.0d+00
           f3(nx,m)=0.0d+00
           f4(nx,m)=0.0d+00
           f5(nx,m)=0.0d+00
           f6(nx,m)=0.0d+00
           f7(nx,m)=0.0d+00
           f8(nx,m)=0.0d+00
           f9(nx,m)=0.0d+00
   40 conTINUE
         DO n=1,nx
           xxl(n)=log(xx(n))
      endDO
         DO m=1,nq
           qql(m)=log(qq(m))
      endDO
         CALL jeppe1(nx,nq,xxl,qql,f1,cc1)
         CALL jeppe1(nx,nq,xxl,qql,f2,cc2)
         CALL jeppe1(nx,nq,xxl,qql,f3,cc3)
         CALL jeppe1(nx,nq,xxl,qql,f4,cc4)
         CALL jeppe1(nx,nq,xxl,qql,f6,cc6)
         CALL jeppe1(nx,nq,xxl,qql,f8,cc8)
         CALL jeppe1(nx,nq,xxl,qql,f9,cc9)
         emc2=2.045
         emb2=18.5
         DO 44 m=1,nqc
           qqlc(m)=qql(m+nqc0)
         DO 44 n=1,nx
           fc(n,m)=f5(n,m+nqc0)
   44 conTINUE
         qqlc(1)=log(emc2)
         CALL jeppe1(nx,nqc,xxl,qqlc,fc,ccc)
         DO 45 m=1,nqb
           qqlb(m)=qql(m+nqb0)
         DO 45 n=1,nx
           fb(n,m)=f7(n,m+nqb0)
   45 conTINUE
         qqlb(1)=log(emb2)
         CALL jeppe1(nx,nqb,xxl,qqlb,fb,ccb)
         init=1

   10    CONTINUE

         xlog  =log(x)
         qsqlog=log(qsq)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
         CALL jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)
         chm=0.0d+00
         IF (qsq.gt.emc2) THEN
           CALL jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endIF
         bot=0.0d+00
         IF (qsq.gt.emb2) THEN
           CALL jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endIF
         x  =xsave
         qsq=q2save

         RETURN
   50 FORMAT(9f10.5)
      END SUBROUTINE mrst2

************************************************************************
      SUBROUTINE jeppe1(nx,my,xx,yy,ff,cc)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         DIMENSION xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my),
     #             ff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),
     #             cl(16),cc(nx,my,4,4),iwt(16,16)
c    #             wt(16,16)                                             KK: This variable has not been used

         DATA iwt/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     #             0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     #            -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
     #             2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
     #             0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     #             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     #             0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
     #             0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
     #            -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     #             0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
     #             9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
     #            -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
     #             2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     #             0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
     #            -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
     #             4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1/

         DO 42 m=1,my
           dx=xx(2)-xx(1)
           ff1(1,m)=(ff(2,m)-ff(1,m))/dx
           dx=xx(nx)-xx(nx-1)
           ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
           DO 41 n=2,nx-1
             ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     #                         ff(n+1,m))
   41   conTINUE
   42 conTINUE
         DO 44 n=1,nx
           dy=yy(2)-yy(1)
           ff2(n,1)=(ff(n,2)-ff(n,1))/dy
           dy=yy(my)-yy(my-1)
           ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
           DO 43 m=2,my-1
             ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     #                         ff(n,m+1))
   43   conTINUE
   44 conTINUE
         DO 46 m=1,my
           dx=xx(2)-xx(1)
           ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
           dx=xx(nx)-xx(nx-1)
           ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
           DO 45 n=2,nx-1
             ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),
     #                          ff2(n,m),ff2(n+1,m))
   45   conTINUE
   46 conTINUE
         DO 53 n=1,nx-1
           DO 52 m=1,my-1
             d1    =xx(n+1)-xx(n)
             d2    =yy(m+1)-yy(m)
             d1d2  =d1*d2
             yy0(1)=ff(n,m)
             yy0(2)=ff(n+1,m)
             yy0(3)=ff(n+1,m+1)
             yy0(4)=ff(n,m+1)
             yy1(1)=ff1(n,m)
             yy1(2)=ff1(n+1,m)
             yy1(3)=ff1(n+1,m+1)
             yy1(4)=ff1(n,m+1)
             yy2(1)=ff2(n,m)
             yy2(2)=ff2(n+1,m)
             yy2(3)=ff2(n+1,m+1)
             yy2(4)=ff2(n,m+1)
             yy12(1)=ff12(n,m)
             yy12(2)=ff12(n+1,m)
             yy12(3)=ff12(n+1,m+1)
             yy12(4)=ff12(n,m+1)
             DO 47 k=1,4
               z(k)   =yy0(k)
               z(k+4) =yy1(k)*d1
               z(k+8) =yy2(k)*d2
               z(k+12)=yy12(k)*d1d2
   47     conTINUE
             DO 49 l=1,16
               xxd=0.0d+00
               DO 48 k=1,16
                 xxd=xxd+iwt(k,l)*z(k)
   48       conTINUE
               cl(l)=xxd
   49     conTINUE
             l=0
             DO 51 k=1,4
               DO 50 j=1,4
                 l=l+1
                 cc(n,m,k,j)=cl(l)
   50       conTINUE
   51     conTINUE
   52   conTINUE
   53 conTINUE

         RETURN
      END SUBROUTINE jeppe1

************************************************************************
      SUBROUTINE jeppe2(x,y,nx,my,xx,yy,cc,z)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z), INTEGER (I,J,K,L,M,N)

         DIMENSION xx(nx),yy(my),cc(nx,my,4,4)

         n=locx(xx,nx,x)
         m=locx(yy,my,y)
         t=(x-xx(n))/(xx(n+1)-xx(n))
         u=(y-yy(m))/(yy(m+1)-yy(m))
         z=0.0d+00
         DO 1 l=4,1,-1
           z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u+cc(n,m,l,2))*u+
     #             cc(n,m,l,1)
    1 conTINUE

         RETURN
      END SUBROUTINE jeppe2

************************************************************************
      INTEGER FUNCTION locx(xx,nx,x)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z), INTEGER (I,J,K,L,M,N)

         DIMENSION xx(nx)

         IF (x.le.xx(1) ) THEN
           locx=1
           RETURN
      endIF
         IF (x.ge.xx(nx)) THEN 
           locx=nx-1
           RETURN
      endIF
         ju=nx+1
         jl=0
    1    IF ((ju-jl).le.1) GOTO 2
           jm=(ju+jl)*5.0d-01
         IF (x.ge.xx(jm) ) THEN
           jl=jm
                           ELSE
           ju=jm
      endIF
         GOTO 1
    2    locx=jl

         RETURN
      END FUNCTION locx

************************************************************************
      REAL FUNCTION polderiv(x1,x2,x3,y1,y2,y3)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z), INTEGER (I,J,K,L,M,N)

         polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*(y2-y3))+
     #            x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))

         RETURN
      END FUNCTION polderiv
