************************************************************************
      FUNCTION CTQ6PDF(Iparton,X,Q)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         LOGICAL Warn

         COMMON    /CtqPar2/Nx,Nt,NfMx
         COMMON   /QCDtable/Alambda,Nfl,Iorder

         DATA Warn/.TRUE./

         SAVE Warn

         IF (X.lt.0.0d+00 .or. X.gt.1.0d+00) THEN
           PRINT *, 'X out of range in CTQ6PDF: ', X
           STOP
      endIF
c        IF (Q.lt.Alambda                  ) THEN                        KK
c          PRINT *, 'Q out of range in CTQ6PDF: ', Q, Alambda            KK
c          STOP                                                          KK
c     endIF                                                              KK
         Q_CTEQ6=max(Q,Alambda+epsilon(Alambda))                         KK
c        IF (Q.lt.Alambda                  ) THEN                        KK
c          PRINT *, Q, Alambda, Q_CTEQ6                                  KK
c     endIF                                                              KK
         IF ((Iparton.lt.-NfMx .or. Iparton.gt.NfMx)) THEN
           IF (Warn) THEN
*            --------------------------------------------------------- *
*            PUT A WARNING FOR CALLING EXTRA FLAVOR.                   *
*            --------------------------------------------------------- *
             Warn=.FALSE.
             PRINT *, 'Warning: Iparton out of range in CTQ6PDF! '
             PRINT *, 'Iparton, MxFlvN0: ',Iparton,NfMx
        endIF
           CTQ6PDF=0.0d+00
           RETURN
      endIF
         CTQ6PDF=PartonX6(Iparton,X,Q_CTEQ6)                             KK
c        IF (CTQ6PDF.lt.0.0d+00) CTQ6PDF=0.0d+00                         KK

         RETURN
      END FUNCTION CTQ6PDF

************************************************************************
      SUBROUTINE ReadTbl(Nu)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         CHARACTER Line*80

         PARAMETER (MXX   =105,
     #              MXQ   = 25,
     #              MXF   =  6,
     #              MaxVal=  3,
     #              MXPQX =(MXF+1+MaxVal)*MXQ*MXX)

         COMMON    /CtqPar1/Al,XV(0:MXX),TV(0:MXQ),UPD(MXPQX)
         COMMON    /CtqPar2/Nx,Nt,NfMx
         COMMON    /XQrange/Qini,Qmax,Xmin
         COMMON   /QCDtable/Alambda,Nfl,Iorder
         COMMON    /Masstbl/Amass(6)
         COMMON    /Valence/MxVal

         READ(Nu,'(A)') Line
         READ(Nu,'(A)') Line
         READ(Nu,    *) Dr,Fl,Al,(Amass(I),I=1,6)
         Iorder =Nint(Dr)
         Nfl    =Nint(Fl)
         Alambda=Al
         READ(Nu,'(A)') Line

         IF (MxVal.eq.3) THEN ! This is the .pds (WKT) format
           READ(Nu,    *) N0,N0,N0,NfMx,N0,N0
           READ(Nu,'(A)') Line
           READ(Nu,    *) NX,NT,N0,N0,N0
           READ(Nu,'(A)') (Line,I=1,4)
           READ(Nu,    *) QINI,QMAX,(aa,TV(I),I=0,NT)
           READ(Nu,'(A)') Line
           READ(Nu,    *) XMIN,aa,(XV(I),I=1,NX)
           XV(0)=0.0d+00
                         ELSE ! This is the old .tbl (HLL) format
           READ(Nu,    *) NX,NT,NfMx
           READ(Nu,'(A)') Line
           READ(Nu,    *) QINI,QMAX,(TV(I),I=0,NT)
           READ(Nu,'(A)') Line
           READ(Nu,    *) XMIN,(XV(I),I=0,NX)
           DO 11 Iq=0,NT
             TV(Iq)=log(log(TV(Iq)/Al))
   11      CONTINUE
      endIF
         Nblk=(NX+1)*(NT+1)
         Npts=Nblk*(NfMx+1+Mxval)
         READ(Nu,        '(A)') Line
         READ(Nu,*,IOSTAT=IRET) (UPD(I),I=1,Npts)

         RETURN
      END SUBROUTINE ReadTbl

************************************************************************
      FUNCTION NextUn()
************************************************************************
*                                                                      *
*     This FUNCTION returns an unallocated FORTRAN i/o unit.           *
*                                                                      *
************************************************************************

         LOGICAL EX

         DO 10 N=10,300
           INQUIRE(UNIT=N, OPENED=EX)
           IF (.NOT.EX) THEN
             NextUn=N
             RETURN
        endIF
   10    CONTINUE

         STOP 'There is no available I/O unit.'

      END FUNCTION NextUn

************************************************************************
      FUNCTION PartonX6(IPRTN,XX,QQ)
************************************************************************
*                                                                      *
*     Given  the  parton distribution  function in  the array U in     *
*     COMMON/PEVLDT/, this routine interpolates to find the parton     *
*     distribution at an arbitray point in x and Q.                    *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         PARAMETER(MXX   =105,
     #             MXQ   = 25,
     #             MXF   =  6,
     #             MaxVal=  3,
     #             MXPQX =(MXF+1+MaxVal)*MXQ*MXX)
 
         COMMON   /CtqPar1/Al,XV(0:MXX),TV(0:MXQ),UPD(MXPQX)
         COMMON   /CtqPar2/Nx,Nt,NfMx
         COMMON   /XQrange/Qini,Qmax,Xmin
         COMMON   /Valence/MxVal
         COMMON /Setchange/Isetch

         DIMENSION fvec(4),fij(4),xvpow(0:mxx)

         DATA      OneP / 1.00001/
         DATA      xpow / 0.3d+00/                                       Choice of interpolation variable
         DATA     nqvec / 4/
c        DATA    ientry / 0/
         DATA X,Q,JX,JQ /-1.0d+00,-1.0d+00,0,0/

         SAVE xvpow,X,Q,JX,JQ,JLX,JLQ,ss,const1,const2,const3,const4,
     #        const5,const6,sy2,sy3,s23,tt,t12,t13,t23,t24,t34,ty2,ty3,
     #        tmp1,tmp2,tdet

          IF ((XX.eq.X) .and. (QQ.eq.Q)) GOTO 99
*         ------------------------------------------------------------ *
*         STORE  THE  POWERS USED  FOR INTERPOLATION ON FIRST CALL     *
*         ------------------------------------------------------------ *
          IF (Isetch.eq.1) THEN
            Isetch=0
            xvpow(0)=0.0d+00
            DO i=1,nx
              xvpow(i)=xv(i)**xpow
         endDO
       endIF
          X =XX
          Q =QQ
          tt=log(log(Q/Al))
*         ------------------------------------------------------------ *
*         Find lower end  of interval  containing x, i.e., get jx such *
*         that xv(jx) .le. x .le. xv(jx+1)...                          *
*         ------------------------------------------------------------ *
          JLx=  -1
          JU =Nx+1
   11     IF (JU-JLx.gt.1) THEN
            JM=0.5*(JU+JLx)
            IF (X.ge.XV(JM)) THEN
              JLx=JM
                             ELSE
              JU =JM
         endIF
            GOTO 11
       endIF
*         ------------------------------------------------------------ *
*         Ix    0   1   2      Jx  JLx         Nx-2     Nx             *
*               |---|---|---|...|---|-x-|---|...|---|---|              *
*         x     0  Xmin               x                 1              *
*         ------------------------------------------------------------ *
          IF (JLx.le.-1                  ) THEN
            PRINT '(A,1pE12.4)','Severe error: x<=0 in PartonX6! x=', x
            STOP
      ELSEIF (JLx.eq. 0                  ) THEN
            Jx=0
      ELSEIF (JLx.le. Nx-2               ) THEN
*           ---------------------------------------------------------- *
*           For interrior points, keep x in the middle, as shown above *
*           ---------------------------------------------------------- *
            Jx=JLx-1
      ELSEIF (JLx.eq. Nx-1 .or. x.lt.OneP) THEN
*           ---------------------------------------------------------- *
*           We tolerate a slight over-shoot of one (OneP=1.00001),     *
*           perhaps  due to roundoff  or whatever,  but not  more than *
*           that. Keep at least 4 points >= Jx.                        *
*           ---------------------------------------------------------- *
            Jx=JLx-2
                                           ELSE
            PRINT '(A,1pE12.4)','Severe error: x > 1 in PartonX6! x=',x
            STOP
       endIF
*         ------------------------------------------------------------ *
*         NOTE: JLx uniquely identIFies the x-bin; Jx does not
*         ------------------------------------------------------------ *

*         ------------------------------------------------------------ *
*         This is the variable to be interpolated in
*         ------------------------------------------------------------ *
          ss=x**xpow
          IF (JLx.ge.2 .and. JLx.le.Nx-2) THEN
*           ---------------------------------------------------------- *
*           Initiation work for "interior bins": store the lattice po- *
*           ints in s                                                  *
*           ---------------------------------------------------------- *
            svec1 =xvpow(jx)
            svec2 =xvpow(jx+1)
            svec3 =xvpow(jx+2)
            svec4 =xvpow(jx+3)
            s12   =svec1-svec2
            s13   =svec1-svec3
            s23   =svec2-svec3
            s24   =svec2-svec4
            s34   =svec3-svec4
            sy2   =ss   -svec2
            sy3   =ss   -svec3
*           ---------------------------------------------------------- *
*           Constants needed for interpolating in s at fixed t lattice *
*           points                                                     *
*           ---------------------------------------------------------- *
            const1=s13/s23
            const2=s12/s23
            const3=s34/s23
            const4=s24/s23
            s1213 =s12+s13
            s2434 =s24+s34
            sdet  =s12*s34-s1213*s2434
            tmp   =sy2*sy3/sdet
            const5=(s34  *sy2-s2434*sy3)*tmp/s12
            const6=(s1213*sy2-s12  *sy3)*tmp/s34
       endIF
*         ------------------------------------------------------------ *
*         Now find lower end of  interval containing Q,  i.e.,  get jq *
*         such that qv(jq) .le. q. le. qv(jq+1)                        *
*         ------------------------------------------------------------ *
          JLq=  -1
          JU =NT+1
   12     IF (JU-JLq.gt.1) THEN
            JM=0.5*(JU+JLq)
            IF (tt.ge.TV(JM)) THEN
              JLq=JM
                              ELSE
              JU =JM
         endIF
            GOTO 12
       endIF
          IF (JLq.le.   0) THEN
            Jq=    0
      ELSEIF (JLq.le.Nt-2) THEN ! Keep q in the middle, as shown above
            Jq=JLq-1
                           ELSE ! JLq.ge.Nt-1 case: Keep at least 4 points >= Jq
            Jq=Nt -3
       endIF
*         ------------------------------------------------------------ *
*         This is the interpolation variable in Q                      *
*         ------------------------------------------------------------ *
          IF (JLq.ge.1 .and. JLq.le.Nt-2) THEN
*           ---------------------------------------------------------- *
*           Store the lattice points in t                              *
*           ---------------------------------------------------------- *
            tvec1=Tv(jq)
            tvec2=Tv(jq+1)
            tvec3=Tv(jq+2)
            tvec4=Tv(jq+3)
            t12  =tvec1  -tvec2
            t13  =tvec1  -tvec3
            t23  =tvec2  -tvec3
            t24  =tvec2  -tvec4
            t34  =tvec3  -tvec4
            ty2  =tt     -tvec2
            ty3  =tt     -tvec3
            tmp1 =t12    +t13
            tmp2 =t24    +t34
            tdet =t12*t34-tmp1*tmp2
       endIF
*         ------------------------------------------------------------ *
*         Get the pdf function values at the lattice points            *
*         ------------------------------------------------------------ *
   99     IF (Iprtn.gt.MxVal) THEN
            Ip=-Iprtn
                              ELSE
            Ip= Iprtn
       endIF
          jtmp=((Ip+NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1
          DO it=1,nqvec
            J1=jtmp+it*(NX+1)
            IF (Jx.eq.    0) THEN
*             -------------------------------------------------------- *
*             For  the first 4  x points,  interpolate  x^2*f(x,Q).    *
*             This applies to the  two lowest bins JLx=0,1.  We can    *
*             not put the JLx.eq.1 bin into the "interrior" section    *
*             (as we do for q), since Upd(J1) is undefined.            *
*             -------------------------------------------------------- *
              fij(1)=0
              fij(2)=Upd(J1+1)*XV(1)**2
              fij(3)=Upd(J1+2)*XV(2)**2
              fij(4)=Upd(J1+3)*XV(3)**2
*             -------------------------------------------------------- *
*             Use Polint which allows x  to be anywhere w.r.t. grid    *
*             -------------------------------------------------------- *
              CALL Polint4 (XVpow(0),Fij(1),ss,Fx)
              IF (x.gt.0.0d+00) Fvec(it)=Fx/x**2
*             -------------------------------------------------------- *
*             Pdf is undefined for x.eq.0                              *
*             -------------------------------------------------------- *
        ELSEIF (JLx.eq.Nx-1) THEN
*             -------------------------------------------------------- *
*             This is the highest x bin:                               *
*             -------------------------------------------------------- *
              CALL Polint4 (XVpow(Nx-3),Upd(J1),ss,Fx)
              Fvec(it)=Fx
                             ELSE
*             -------------------------------------------------------- *
*             for all interior  points, use  Jon's in-line function    *
*             This applied to (JLx.ge.2 .and. JLx .le. Nx-2)           *
*             -------------------------------------------------------- *
              sf2=Upd(J1+1)
              sf3=Upd(J1+2)
              g1 = sf2*const1-sf3*const2
              g4 =-sf2*const3+sf3*const4
              Fvec(it)=(const5*(Upd(J1)-g1)+const6*(Upd(J1+3)-g4)
     #                 +sf2*sy3-sf3*sy2)/s23
         endIF
       endDO
*         ------------------------------------------------------------ *
*         We now have the four values Fvec(1:4) interpolate in t       *
*         ------------------------------------------------------------ *
          IF (JLq.le.   0) THEN
*           ---------------------------------------------------------- *
*           1st Q-bin, as well as extrapolation to lower Q             *
*           ---------------------------------------------------------- *
            CALL Polint4(TV(0),   Fvec(1),tt,ff)
      ELSEIF (JLq.ge.Nt-1) THEN
*           ---------------------------------------------------------- *
*           Last Q-bin, as well as extrapolation to higher Q           *
*           ---------------------------------------------------------- *
            CALL Polint4(TV(Nt-3),Fvec(1),tt,ff)
                           ELSE
*           ---------------------------------------------------------- *
*           Interrior bins: (JLq.ge.1 .and. JLq.le.Nt-2) which include *
*           JLq.eq.1 and JLq.eq.Nt-2,  since  Upd  is defined  for the *
*           full range QV(0:Nt) (in contrast to XV)                    *
*           ---------------------------------------------------------- *
            tf2=fvec(2)
            tf3=fvec(3)
            g1 =( tf2*t13-tf3*t12)/t23
            g4 =(-tf2*t34+tf3*t24)/t23
            h00=((t34 *ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     #          +(tmp1*ty2-t12 *ty3)*(fvec(4)-g4)/t34)
            ff =(h00*ty2*ty3/tdet+tf2*ty3-tf3*ty2)/t23
       endIF
          PartonX6=ff

         RETURN
      END FUNCTION PartonX6

************************************************************************
      SUBROUTINE POLINT4(XA,YA,X,Y)
************************************************************************
*                                                                      *
*     This SUBROUTINE is based on the POLINT routine from "Numeri-     *
*     cal recipes", but assuming N=4 and ignoring  the error esti-     *
*     mation suggested by Z. Sullivan.                                 *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         DIMENSION XA(*),YA(*)

          H1 =XA(1)-X
          H2 =XA(2)-X
          H3 =XA(3)-X
          H4 =XA(4)-X
          W  =YA(2)-YA(1)
          DEN=W/(H1-H2)
          D1 =H2*DEN
          C1 =H1*DEN
          W  =YA(3)-YA(2)
          DEN=W/(H2-H3)
          D2 =H3*DEN
          C2 =H2*DEN
          W  =YA(4)-YA(3)
          DEN=W/(H3-H4)
          D3 =H4*DEN
          C3 =H3*DEN
          W  =C2-D1
          DEN=W/(H1-H3)
          CD1=H3*DEN
          CC1=H1*DEN
          W  =C3-D2
          DEN=W/(H2-H4)
          CD2=H4*DEN
          CC2=H2*DEN
          W  =CC2-CD1
          DEN=W/(H1-H4)
          DD1=H4*DEN
          DC1=H1*DEN

          IF ((H3+H4).lt.0.0d+00) THEN
            Y=YA(4)+D3+CD2+DD1
      ELSEIF ((H2+H3).lt.0.0d+00) THEN
            Y=YA(3)+D2+CD1+DC1
      ELSEIF ((H1+H2).lt.0.0d+00) THEN
            Y=YA(2)+C2+CD1+DC1
                                  ELSE
            Y=YA(1)+C1+CC1+DC1
       endIF

         RETURN
      END SUBROUTINE POLINT4