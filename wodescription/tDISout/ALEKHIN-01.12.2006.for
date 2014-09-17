************************************************************************
c     SUBROUTINE ALEKHIN06(xb,q2,pdfs,dpdfs,npdf,npar)
      SUBROUTINE ALEKHIN06(xb,q2,pdfs,      npdf,npar)
************************************************************************
*                                                                      *
*     This is a code for the  NNLO parton distributions in the va-     *
*     riable-flavor-number (VFN) schem with account of their expe-     *
*     rimental  and theoretical  uncertainties.  The Q^2  range is     *
*     0.8 < Q^2 < 2.0d+8 GeV^2, the x range is 1.0d-7 < x < 1.0d+0     *
*     (for the values of  PDFs and strong coupling constant at Q^2     *
*     < 0.8 GeV^2  (x < 1.0d-7) their  values  at  Q^2 = 0.8 GeV^2     *
*     (x=1.0d-7) are returned).                                        *
*                                                                      *
*     The output array  PDFS contains  fitted values of the strong     *
*     coupling  constant  and the  parton distributions  at  given     *
*     x  and  Q:  PDFS(0) - \alpha_s;  PDFS(1) - valence u-quarks;     *
*     PDFS(2) - valence d-quarks; PDFS(3) - gluons;  PDFS(4) - sea     *
*     u-quarks; PDFS(5) - s-quarks; PDFS(6) -sea d-quarks; PDFS(7)     *
*     - c-quarks; PDFS(8) - b-quarks; PDFS(9) - t-quarks.  npdf is     *
*     the  number of  PDFs returned  (npdf=9 for the  VFN scheme).     *
*     The output array dpdfs(0:npdf,npar) contains  derivatives of     *
*     \alpha_s and the PDFs on the fitted parameters with the num-     *
*     ber of the parameters returned in npar.  With the  derivati-     *
*     ves of \alpha_s  included one can take into account the cor-     *
*     relations  of the  fitted PDFs  with  \alpha_s as  well. All     *
*     derivatives are transformed  to the orthonormal basis of ei-     *
*     genvectors of the  parameters error matrix.  For this reason     *
*     the variation of the PDFs in the  derivatives directions can     *
*     be  performed independently.  For example the  dispersion of     *
*     the i-th PDF can be stored in DELPDF using the code.             *
*                                                                      *
*     -------------------------------------                            *
*        DELPDF=0.                                                     *
*        DO k=1,npar                                                   *
*          DELPDF=DELPDF+dpdfs(i,k)**2                                 *
*     endDO                                                            *
*     -------------------------------------                            *
*     and its random value can be stored in RPDF using the code        *
*     -------------------------------------                            *
*        RPDF=pdfs(i)                                                  *
*        DO k=1,npar                                                   *
*          s=0.0                                                       *
*          DO l=1,96                                                   *
*            s=s+(2*rndm(xxx)-1)/sqrt(32.0)                            *
*       endDO                                                          *
*          RPDF=RPDF+s*dpdfs(i,k)                                      *
*     endDO                                                            *
*     -------------------------------------                            *
*                                                                      *
*     REFERENCE: S. Alekhin, K. Melnikov and F. Petriello,  "Fixed     *
*     target Drell-Yan data and NNLO QCD fits  of parton distribu-     *
*     tion functions," Phys. Rev. D74 (2006) 054033                    *
*     [hep-ph/0606237].                                                *
*                                                                      *
*     COMMENTS: Sergey.Alekhin@ihep.ru                                 *
*                                                                      *
*                       Initial version: Dec 2006                      *
*                                                                      *
************************************************************************

         USE InpOutUnits

         IMPLICIT NONE

         CHARACTER(1) pdford
         CHARACTER(3) pdfschem
           INTEGER k,i,n,m,kx,nxbb,npdf,npar,KORD,kschem,kords,
     #             kschems,nport
              REAL x,qsq,dels,delx,x1,delx1,xlog1,xd,b,aa,ss,f0,fp,
     #             fm,xb,q2,df0,dfp,dfm,xmin,xmax,qsqmin,qsqmax

         INTEGER,PARAMETER::
     #           nxb =99,
     #           nq  =20,
     #           np  = 9,
     #           nvar=23

         INTEGER nexp(0:np)
            REAL f(nxb,nq+1,0:np),
     #           xx(nxb),
     #           fsp(nxb),
     #           bs (nxb),
     #           cs (nxb),
     #           ds (nxb),
     #           bsp(nxb,nq+1,0:np),
     #           csp(nxb,nq+1,0:np),
     #           dsp(nxb,nq+1,0:np),
     #           bspd(nvar,nxb,nq+1,0:np),
     #           cspd(nvar,nxb,nq+1,0:np),
     #           dspd(nvar,nxb,nq+1,0:np),
     #           pdfs(0:np),
     #           dpdfs(0:np,nvar),
     #           df(nvar,0:np,nxb,nq+1)

         DIMENSION pdfschem(0:1),pdford(3)

         DATA nport                  /1/
         DATA pdford                 /'1','2','3'/
         DATA pdfschem               /'ffn','vfn'/
         DATA nexp                   /0,3,4,5,5,5,5,5,5,5/
         DATA xmin,xmax,qsqmin,qsqmax/1.0d-07,1.0d+00,0.8d+00,2.0d+08/
         DATA KORDS,KSCHEMS          /-1,-1/
         DATA KORD,KSCHEM            /3,1/

         SAVE kords,kschems,f,df,dels,delx,x1,delx1,xlog1,nxbb,xx

         IF (kschem.eq.0) THEN
           npdf=6
                          ELSE
           npdf=9
      endIF
         npar=nvar
*        ------------------------------------------------------------- *
*        RESET ARRAYS                                                  *
*        ------------------------------------------------------------- *
         DO i=0,npdf
           pdfs(i)=0
           DO k=1,npar
             dpdfs(i,k)=0
        endDO
      endDO
         IF ((kords.eq.kord).and.(kschems.eq.kschem)) GOTO 10
         kords  =kord
         kschems=kschem
c KK     dels   =(dlog(dlog(qsqmax/0.04d+00))-
c KK #            dlog(dlog(qsqmin/0.04d+00)))/dble(nq-1)
         dels   =(log(log(qsqmax*25))-                                   KK
     #            log(log(qsqmin*25)))/dble(nq-1)                        KK

         nxbb   =nxb/2
c KK     x1     =0.3d+00
         x1     =3.0d-01                                                 KK
         xlog1  = log(x1)
         delx   =(log(x1)-log(xmin))/dble(nxbb-1)
         DELX1  =(1-x1)**2/dble(nxbb+1)
*        ------------------------------------------------------------- *
*        X GRID                                                        *
*        ------------------------------------------------------------- *
         DO kx=1,nxbb
           xx(kx)=exp(log(xmin)+delx*dble(kx-1))
      endDO
         DO kx=nxbb+1,nxb-1
           xx(kx)=1-sqrt(abs((1-x1)**2-delx1*dble(kx-nxbb)))
      endDO
         xx(nxb)=1
*        ------------------------------------------------------------- *
*        READ INPUT TABLES                                             *
*        ------------------------------------------------------------- *
c KK     PRINT *,'***** Reading PDFs from tables *****'                  KK
         OPEN(unit=nport,status='old',
     #        FILE=IniPP//'PDFLIB-like/ALEKHIN/2006/'//                  KK
     #             'a06.dpdfs_'//pdford(kord)//'_'//pdfschem(kschem))    KK

         DO n=1,nxb-1
           DO m=1,nq
             DO i=0,npdf
               READ(nport,*) (df(k,i,n,m), k=1,npar)
          endDO
        endDO
      endDO
         CLOSE(unit=nport)
         DO k=1,npar
           DO i=0,npdf
             DO m=1,nq
               IF (i.ne.0) THEN
                 df(k,i,nxb,m)=0
                           ELSE
                 df(k,i,nxb,m)=df(k,i,nxb-1,m)
            endIF
               DO n=1,nxb
                 fsp(n)=df(k,i,n,m)
            endDO
               CALL spline(nxb,xx,fsp,bs,cs,ds)
               DO n=1,nxb
                 bspd(k,n,m,i)=bs(n)
                 cspd(k,n,m,i)=cs(n)
                 dspd(k,n,m,i)=ds(n)
            endDO
          endDO
        endDO
      endDO
         OPEN(unit=nport,status='old',err=199,
     #        FILE=IniPP//'PDFLIB-like/ALEKHIN/2006/'//                  KK
     #             'a06.pdfs_'//pdford(kord)//'_'//pdfschem(kschem))     KK
         DO n=1,nxb-1
           DO m=1,nq
             READ(nport,*) (f(n,m,i),i=0,npdf)
        endDO
      endDO
         DO i=0,npdf
           DO m=1,nq
             IF (i.ne.0) THEN
               f(nxb,m,i)=0
                         ELSE
               f(nxb,m,i)=f(nxb-1,m,i)
          endIF
             DO n=1,nxb-1
               f(n,m,i)=f(n,m,i)/(1-xx(n))**nexp(i)
          endDO
             DO n=1,nxb
               fsp(n)=f(n,m,i)
          endDO
             CALL SPLINE(nxb,xx,fsp,bs,cs,ds)
             DO n=1,nxb
               bsp(n,m,i)=bs(n)
               csp(n,m,i)=cs(n)
               dsp(n,m,i)=ds(n)
          endDO
        endDO
      endDO
         CLOSE(unit=nport)
   10    CONTINUE
*        ------------------------------------------------------------- * KK
         DO i=0,npdf                                                     KK
           DO m=1,nq                                                     KK
             DO n=1,nxb                                                  KK
               IF (bsp(n,m,i).gt.1.0d+05) bsp(n,m,i)=0.0d+00             KK
               IF (csp(n,m,i).gt.1.0d+05) csp(n,m,i)=0.0d+00             KK
               IF (dsp(n,m,i).gt.1.0d+05) dsp(n,m,i)=0.0d+00             KK
               IF (bsp(n,m,i).lt.1.0d-05) bsp(n,m,i)=0.0d+00             KK
               IF (csp(n,m,i).lt.1.0d-05) csp(n,m,i)=0.0d+00             KK
               IF (dsp(n,m,i).lt.1.0d-05) dsp(n,m,i)=0.0d+00             KK
          endDO                                                          KK
        endDO                                                            KK
      endDO                                                              KK
*        ------------------------------------------------------------- * KK

c KK     IF ((q2.lt.qsqmin).or.(q2.gt.qsqmax)) print 99,q2,qsqmin,qsqmax
c KK     IF ((xb.lt.xmin  ).or.(xb.gt.xmax  )) print 98,xb,  xmin,  xmax

         x  =max( xb,  xmin+epsilon(1.0d+00))                            KK
         x  =min(  x,  xmax-epsilon(1.0d+00))                            KK
         qsq=max( q2,qsqmin+epsilon(1.0d+00))                            KK
         qsq=min(qsq,qsqmax-epsilon(1.0d+00))                            KK

         IF (x.gt.x1) THEN
           xd=(1-x1)**2-(1-x)**2
           n =int(xd/delx1)+nxbb
                      ELSE
           xd=log(x)-xlog1
           n =nxbb+int(xd/DELX)-1
      endIF
         aa=x-xx(n)
c KK     ss=dlog(dlog(qsq/0.04d+00))-dlog(dlog(qsqmin/0.04d+00))
         ss=log(log(qsq*25))        -log(log(qsqmin*25))                 KK
         m =int(ss/dels)+1
         b =ss/dels-dble(m)+1

         DO i=0,npdf
c KK       f0=f(n,m,i)+
c KK #        aa*bsp(n,m,i)+aa**2*csp(n,m,i)+aa**3*dsp(n,m,i)
c KK       fp=f(n,m+1,i)+
c KK #        aa*bsp(n,m+1,i)+aa**2*csp(n,m+1,i)+aa**3*dsp(n,m+1,i)

           f0=f(n, m,   i)+aa*(bsp(n,m,  i)+                             KK
     #        aa*(csp(n, m,   i)+aa*dsp(n, m,   i)))                     KK
           fp=f(n,(m+1),i)+aa*(bsp(n,(m+1),i)+                           KK
     #        aa*(csp(n,(m+1),i)+aa*dsp(n,(m+1),i)))                     KK
           IF (m.ge.2) THEN
c KK         fm     =f(n,m-1,i)+aa*bsp(n,m-1,i)+aa**2*csp(n,m-1,i)+
c KK #               aa**3*dsp(n,m-1,i)
c KK         pdfs(i)=fm*b*(b-1)/2+f0*(1-b**2)+fp*b*(b+1)/2

             fm     =f(n,(m-1),i)+aa*(bsp(n,(m-1),i)+                    KK
     #               aa*(csp(n,(m-1),i)+aa*dsp(n,(m-1),i)))              KK
             pdfs(i)=fm*b*(b-1)*0.5+f0*(1-b**2)+fp*b*(b+1)*0.5           KK
                       ELSE
             pdfs(i)=f0*(1-b)+fp*b
        endIF
           pdfs(i)=pdfs(i)*(1-x)**nexp(i)
           DO k=1,npar
c KK         df0=df(k,i,n,m)+aa*bspd(k,n,m,i)+
c KK #           aa**2*cspd(k,n,m,i)+aa**3*dspd(k,n,m,i)
c KK         dfp=df(k,i,n,m+1)+aa*bspd(k,n,m+1,i)+
c KK #           aa**2*cspd(k,n,m+1,i)+aa**3*dspd(k,n,m+1,i)

             df0=df(k,i,n, m   )+aa*(bspd(k,n, m,   i)+                  KK
     #           aa*(cspd(k,n, m,   i)+aa*dspd(k,n, m,   i)))            KK
             dfp=df(k,i,n,(m+1))+aa*(bspd(k,n,(m+1),i)+                  KK
     #           aa*(cspd(k,n,(m+1),i)+aa*dspd(k,n,(m+1),i)))            KK

             IF (m.ge.2) THEN 
c KK           dfm       =df(k,i,n,m-1)+aa*bspd(k,n,m-1,i)+
c KK #                    aa**2*cspd(k,n,m-1,i)+aa**3*dspd(k,n,m-1,i)

               dfm       =df(k,i,n,(m-1))+aa*(bspd(k,n,(m-1),i)+         KK
     #                    aa*(cspd(k,n,(m-1),i)+aa*dspd(k,n,(m-1),i)))   KK

               dpdfs(i,k)=dfm*b*(b-1)*0.5+df0*(1-b**2)+dfp*b*(b+1)*0.5
                         ELSE
               dpdfs(i,k)=df0*(1-b)+dfp*b
          endIF
        endDO
      endDO

         RETURN
  199 PRINT *,'The PDF set is inavailable (FILE:',
     #        'a06.pdfs_'//pdford(kord)//'_'//pdfschem(kschem),')'
         RETURN
   99 FORMAT('A06 WARNING: Q^2 VALUE IS OUT OF RANGE',3g12.3)
   98 FORMAT('A06 WARNING: X   VALUE IS OUT OF RANGE',3g12.3)
      END SUBROUTINE ALEKHIN06

************************************************************************
      SUBROUTINE SPLINE(N,X,Y,B,C,D)
************************************************************************
*                                                                      *
*     Calculate the coefficients B,C,D in a cubic spline interpola-    *
*     tion. Interpolation subroutines are taken from G.E. Forsythe,    *
*     M.A. Malcolm and C.B. Moler, "Computer methods for mathemati-    *
*     cal computations," Prentice-Hall, 1977.                          *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         DIMENSION X(N),Y(N),B(N),C(N),D(N)

         NM1=N-1
         IF (N.LT.2) RETURN
         IF (N.LT.3) GOTO 250
         D(1)=X(2)-X(1)
         C(2)=(Y(2)-Y(1))/D(1)
         DO 210 K=2,NM1
           D(K)  =X(K+1)-X(K)
           B(K)  =2*(D(K-1)+D(K))
           C(K+1)=(Y(K+1)-Y(K))/D(K)
           C(K)  =C(K+1)-C(K)
  210    CONTINUE
         B(1)=-D(1)
         B(N)=-D(N-1)
         C(1)= 0.0d+00
         C(N)= 0.0d+00
         IF (N.EQ.3) GOTO 215
         C(1)= C(3)  /(X(4)-X(2))  -C(2)  /(X(3)  -X(1)  )
         C(N)= C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
         C(1)= C(1)*D(1)**2/(X(4)-X(1))
         C(N)=-C(N)*D(N-1)**2/(X(N)-X(N-3))
  215    CONTINUE
         DO 220 K=2,N
           T   =D(K-1)/B(K-1)
           B(K)=B(K)-T*D(K-1)
           C(K)=C(K)-T*C(K-1)
  220    CONTINUE
         C(N)=C(N)/B(N)
         DO 230 IB=1,NM1
           K   =N-IB
           C(K)=(C(K)-D(K)*C(K+1))/B(K)
  230    CONTINUE
         B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2*C(N))
         DO 240 K=1,NM1
           B(K)=(Y(K+1)-Y(K))/D(K)-D(K)*(C(K+1)+2*C(K))
           D(K)=(C(K+1)-C(K))/D(K)
           C(K)=3*C(K)
  240    CONTINUE
         C(N)=3*C(N)
         D(N)=D(N-1)
         RETURN
  250    CONTINUE
         B(1)=(Y(2)-Y(1))/(X(2)-X(1))
         C(1)=0.0d+00
         D(1)=0.0d+00
         B(2)=B(1)
         C(2)=0.0d+00
         D(2)=0.0d+00

        RETURN
      END SUBROUTINE SPLINE

************************************************************************
      SUBROUTINE ALEKHIN_TEST
************************************************************************
*                                                                      *
*     Template code to print the PDFs with their errors.               *
*                                                                      *
************************************************************************

         USE InpOutUnits

         IMPLICIT REAL (A-H,O-Z)

         REAL pdfs(0:9),dpdfs(0:9,23),delpdf(0:9)

         nbin =100
         q2   =1.0d+01
         xbmin=1.0d-08
         xbmax=1.0d+00
         bin  =(log(xbmax)-log(xbmin))/(nbin-1)

         OPEN(30,FILE=IniPP//'ALEKHIN06_TEST.dat')                       KK

         DO n=1,nbin
           xb=xbmin*exp(bin*(n-1))
*          select the NLO VFN PDFs (nominal set)
c KK       CALL ALEKHIN06(xb,q2,pdfs,dpdfs,npdf,npar)
           CALL ALEKHIN06(xb,q2,pdfs,npdf,npar)
*          PDFs are returned in PDFS(1:npdf),
*          value of alpha_s is PDFS(0)
           DO i=0,npdf
             delpdf(i)=0.
             DO k=1,npar
               delpdf(i)=delpdf(i)+dpdfs(i,k)**2
          endDO
        endDO
           WRITE(30,10)  xb, (pdfs  (i), i=0,npdf)                       KK
c KK       WRITE(30,11) (sqrt(delpdf(i)),i=0,npdf)
      endDO
         CLOSE(30)                                                       KK

         RETURN
c  10 FORMAT(f10.6,12e14.5) 
   11 FORMAT('          ',12e14.5) 
   10 FORMAT(1PE11.5,7(1PE12.5))                                         KK
      END SUBROUTINE ALEKHIN_TEST
