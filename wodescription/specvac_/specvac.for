************************************************************************
      PROGRAM specvac
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/05  *
************************************************************************

         USE PhysMathConstants
         USE OscMatNeeds, ONLY: Nsa,Nhi,Nt,Nfl
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
              INTEGER,PARAMETER::
     #                Nfilof  = 110,
     #                NE      = 100,     
     #                MinCal  = 100
                 REAL,PARAMETER::
     #                Xlow    = zero,
     #                Xupp    = one,
     #                Enu_min = 1.0d-01,
     #                Enu_max = 1.0d+04,
     #                RelErr  = 1.0d-06
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*58
     #                filonf
         CHARACTER*3
     #                Sp,
     #                DMn(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                SA,
     #                san(Nsa)/'F','V'/,
     #                hin(Nhi)/'n','i'/,
     #                NTn(Nt)/'n','a'/,
     #                fln(Nfl)/'e','m','t'/
              INTEGER
     #                n_Enu,
     #                n_DM,n_hi,n_sa,n_NT,n_fl,nNTn
                 REAL
     #                E(NE),Res(NE),
     #                start_time,finish_time,time

         COMMON    /N/N                                                  Atmospheric neutino spectrum
         COMMON /n_sa/n_sa
         COMMON /n_hi/n_hi
         COMMON /n_fl/n_fl
         COMMON /n_NT/n_NT
         COMMON/C_lim/C_min,deltaC
         COMMON /E_nu/E_nu
         
         EXTERNAL underint,uif
         
         C_min =-1.0d+00
         C_max = 1.0d+00
         deltaC= C_max-C_min
         
         lgE_min=log10(Enu_min)
         lgE_max=log10(Enu_max)
         steplgE=(lgE_max-lgE_min)/(NE-1)
         DO n_NE=1,NE
           E(n_NE)=ten**(lgE_min+(n_NE-1)*steplgE)
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

         n_DM=1

         DO n_sa=1,Nsa
!         n_sa= 1
         DO n_hi=1,Nhi
         !n_hi= 1
         DO n_fl=1,Nfl
         !n_fl= 1
         DO nNTn=1,Nt
         !nNTn=1
           IF (nNTn.EQ.1) THEN
             n_NT= 1
                          ELSE
             n_NT=-1
        endIF

           Palphabeta=oscvac3(1,1,one,one)

		   DO n_NE=1,NE
             E_nu=E(n_NE)
             CALL GeMInt(uif,Res(n_NE),Xlow,Xupp,*100)
             CALL cpu_time(finish_time)
             time=finish_time-start_time
             WRITE(*,*) n_NE,time
        endDO

         filonf=Out//'spectrt/'//DMn(n_DM)//Sp//SA//san(n_sa)//
     #                              hin(n_hi)//fln(n_fl)//NTn(nNTn)//ext
         OPEN(Nfilof,FILE=filonf)
		 DO n_NE=1,NE
           WRITE(Nfilof,101) E(n_NE),Res(n_NE)                          !/deltaC
      endDO
         CLOSE(Nfilof)

         CALL GeMInf

         CALL cpu_time(finish_time)
         time=finish_time-start_time
         WRITE(*,*) time
      endDO
      endDO
      endDO
      endDO
         
         STOP 'THE END OF PROGRAM specvac'

   99    STOP 'ERROR WITH GeMSet in specvac'
  100    STOP 'ERROR WITH GeMInt in specvac'

  101 FORMAT(2(1PE16.8))

      END PROGRAM specvac
