************************************************************************
      PROGRAM spec2lm
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/04  *
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
              INTEGER,PARAMETER::
     #                Nfen    = 120,
     #                Nfea    = 121,
     #                Nfmn    = 122,
     #                Nfma    = 123,
     #                Nftn    = 124,
     #                Nfta    = 125,
     #                NE      = 100,     
     #                MinCal  = 100
                 REAL,PARAMETER::
     #                Xlow    = zero,
     #                Xupp    = one,
     #                E_nu_min= 1.0d-01,
     #                E_nu_max= 1.0d+04,
     #                RelErr  = 1.0d-08
         CHARACTER(*),PARAMETER::
     #                ext='.dat',nl='_x'
         CHARACTER*59
     #                fnen,fnea,fnmn,fnma,fntn,fnta
         CHARACTER*3
     #                Sp,
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                SA,
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/
                 REAL
     #                E(NE),Resen(NE),Resea(NE),
     #                Resmn(NE),Resma(NE),Restn(NE),Resta(NE),
     #                start_time,finish_time,time

         COMMON    /N/N                                                  Atmospheric neutino spectrum
         COMMON/C_lim/C_min,deltaC
         COMMON /E_nu/E_nu
         
         EXTERNAL underint,uifen,uifea,uifmn,uifma,uiftn,uifta
         
         CALL cpu_time(start_time)
         
         C_min =-1.0d+00
         C_max = 1.0d+00
         deltaC= C_max-C_min
         
         lgE_nu_min = log10(E_nu_min)
         lgE_nu_max = log10(E_nu_max)
         steplgE_nu = (lgE_nu_max-lgE_nu_min)/(NE-1)
         DO n_NE=1,NE
           E(n_NE)=10**(lgE_nu_min+(n_NE-1)*steplgE_nu)
      endDO

         setpar=underint(one)
         
         CALL GeMSet(underint,one,Xlow,Xupp,RelErr,MinCal,*99)
         
         N=1
         set=dFACN_dE(one,one,N)

         SELECTCASE(N)
               CASE(1)
                     WRITE(*,*) ' AN_Honda11 + AN_ISU_HE + AN_SHE '
                     Sp='H11'
                     WRITE(*,*) ' Maximal solar activity '
                     SA='x'
               CASE(2)
                     WRITE(*,*) ' CORTout '
                     Sp='CRT'
                     WRITE(*,*) ' Minimal solar activity '
                     SA='n'
      endSELECT

         n_DM=2
         n_hi=1

         n_NT= 1
         CALL OscMat3
         
		 DO n_NE=1,NE
           E_nu=E(n_NE)
           CALL GeMInt(uifen,Resen(n_NE),Xlow,Xupp,*100)
           CALL GeMInt(uifmn,Resmn(n_NE),Xlow,Xupp,*100)
           CALL GeMInt(uiftn,Restn(n_NE),Xlow,Xupp,*100)
           CALL cpu_time(finish_time)
           time=finish_time-start_time
           WRITE(*,*) n_NE,time
      endDO

         n_NT=-1
         CALL OscMat3
         
		 DO n_NE=1,NE
           E_nu=E(n_NE)
           CALL GeMInt(uifea,Resea(n_NE),Xlow,Xupp,*100)
           CALL GeMInt(uifma,Resma(n_NE),Xlow,Xupp,*100)
           CALL GeMInt(uifta,Resta(n_NE),Xlow,Xupp,*100)
           CALL cpu_time(finish_time)
           time=finish_time-start_time
           WRITE(*,*) n_NE,time
      endDO

         fnen=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(1)//
     #                                                   NTn(1)//nl//ext
         fnea=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(1)//
     #                                                   NTn(2)//nl//ext
         fnmn=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(2)//
     #                                                   NTn(1)//nl//ext
         fnma=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(2)//
     #                                                   NTn(2)//nl//ext
         fntn=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(3)//
     #                                                   NTn(1)//nl//ext
         fnta=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//fln(3)//
     #                                                   NTn(2)//nl//ext
         OPEN(Nfen,FILE=fnen)
         OPEN(Nfea,FILE=fnea)
         OPEN(Nfmn,FILE=fnmn)
         OPEN(Nfma,FILE=fnma)
         OPEN(Nftn,FILE=fntn)
         OPEN(Nfta,FILE=fnta)

		 DO n_NE=1,NE
           WRITE(Nfen,101) E(n_NE),Resen(n_NE)                          !/deltaC
           WRITE(Nfea,101) E(n_NE),Resea(n_NE)                          !/deltaC
           WRITE(Nfmn,101) E(n_NE),Resmn(n_NE)                          !/deltaC
           WRITE(Nfma,101) E(n_NE),Resma(n_NE)                          !/deltaC
           WRITE(Nftn,101) E(n_NE),Restn(n_NE)                          !/deltaC
           WRITE(Nfta,101) E(n_NE),Resta(n_NE)                          !/deltaC
      endDO
         CLOSE(Nfen)
         CLOSE(Nfea)
         CLOSE(Nfmn)
         CLOSE(Nfma)
         CLOSE(Nftn)
         CLOSE(Nfta)

!         CALL GeMInf
         
         CALL cpu_time(finish_time)
         time=finish_time-start_time
         WRITE(*,*) time

         STOP 'THE END OF PROGRAM spec2lm'

   99    STOP 'ERROR WITH GeMSet in spec2lm'
  100    STOP 'ERROR WITH GeMInt in spec2lm'

  101 FORMAT(2(1PE16.8))

      END PROGRAM spec2lm