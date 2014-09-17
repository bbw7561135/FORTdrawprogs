************************************************************************
      PROGRAM vacuum3
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/08/21  *
************************************************************************

         USE OscMatNeeds

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfilof = 110,
     #                NEnu   = 50,
     #                Nzth   = 50
                 REAL,PARAMETER::
     #                re     = R_Earth*1.0d-03,                          Earth radius [km]!*1.0d+05
     #                ra     = zero,                                     Atmosphere depth [km]!*1.0d+05
     #                Enu_fix= 1.0d+00,                                  Fixed energy [GeV]
     #                zth_fix= pi,                                       Fixed angle
     #                Enu_min= E_min*1.0d-09,                            Energy limits [GeV]
     #                Enu_max= E_max*1.0d-09,
     #                zth_min= ThetaMin,                                 Angle limits
     #                zth_max= ThetaMax
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*60
     #                filon
         CHARACTER*57
     #                filonf
         CHARACTER*3
     #                DMn(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                san(Nsa)/'F','V'/,
     #                hin(Nhi)/'n','i'/,
     #                NTn(Nt)/'n','a'/,
     #                fln(Nfl)/'e','m','t'/
              INTEGER
     #                n_Enu,n_zth,
     #                n_DM,n_hi,n_sa,n_NT,n_f1,n_f2,nNTn,
     #                varsw
                 REAL
     #                Enu(NEnu),lgEnu(NEnu),zth(Nzth),
     #                start_time,finish_time,time

         COMMON /n_sa/n_sa
         COMMON /n_hi/n_hi
         COMMON /n_NT/n_NT

         CALL cpu_time(start_time)

         n_DM= 1

         lgE_min=log10(Enu_min)
         lgE_max=log10(Enu_max)
         steplgE=(lgE_max-lgE_min)/(NEnu-1)
         stepzth=(zth_max-zth_min)/(Nzth-1)

         DO n_Enu=1,NEnu
           lgEnu(n_Enu)=lgE_min+(n_Enu-1)*steplgE
           Enu(n_Enu)=ten**lgEnu(n_Enu)
      endDO
         DO n_zth=1,Nzth
           zth(n_zth)=zth_min+(n_zth-1)*stepzth
      endDO

         varsw= 3                                                        1 for fixed E_nu, 2 for fixed cos, 3 for 3D

         SELECTCASE(varsw)
               CASE(1)
                     WRITE(*,*) ' Neutrino energy is fixed '
               CASE(2)
                     WRITE(*,*) ' Zenith angle is fixed '
               CASE(3)
                     WRITE(*,*)
     #                    ' Neutrino energy and zenith angle both vary '
       endSELECT

         DO n_sa=1,Nsa
         !n_sa= 1
         DO n_hi=1,Nhi
         !n_hi= 1
         DO nNTn=1,Nt
         !nNTn=1
           IF (nNTn.EQ.1) THEN
             n_NT= 1
                          ELSE
             n_NT=-1
        endIF
           WRITE(*,*) ' Analysis: ',san(n_sa),', hierarchy: ',hin(n_hi),
     #                ', type: ',NTn(nNTn)

           Palphabeta=oscvac3(one,one)

          SELECTCASE(varsw)
                CASE(1)
                      filonf=OutDir//'probanu/vc/'//DMn(n_DM)//NTn(nNTn)
     #                                       //san(n_sa)//hin(n_hi)//ext
                      OPEN(Nfilof,FILE=filonf)
                      DO n_zth=1,Nzth
                        root=sqrt((ra+re)**2-re**2*(sin(zth(n_zth)))**2)
                        bite=re*cos(zth(n_zth))
                        eL=root-bite
                        WRITE(Nfilof,300) cos(zth(n_zth)),
     #                         Pab(1,1,Enu_fix,eL),Pab(2,1,Enu_fix,eL),
     #                         Pab(3,1,Enu_fix,eL),Pab(1,2,Enu_fix,eL),
     #                         Pab(2,2,Enu_fix,eL),Pab(3,2,Enu_fix,eL),
     #                         Pab(1,3,Enu_fix,eL),Pab(2,3,Enu_fix,eL),
     #                         Pab(3,3,Enu_fix,eL)
                   endDO
                      CLOSE(Nfilof)
                CASE(2)
                      filonf=OutDir//'probanu/vE/'//DMn(n_DM)//NTn(nNTn)
     #                                       //san(n_sa)//hin(n_hi)//ext
                      OPEN(Nfilof,FILE=filonf)
                      root=sqrt((ra+re)**2-re**2*(sin(zth_fix))**2)
                      bite=re*cos(zth_fix)
                      eL=root-bite
                      DO n_Enu=1,NEnu
                        WRITE(Nfilof,300) Enu(n_Enu),
     #                  Pab(1,1,Enu(n_Enu),eL),Pab(2,1,Enu(n_Enu),eL),
     #                  Pab(3,1,Enu(n_Enu),eL),Pab(1,2,Enu(n_Enu),eL),
     #                  Pab(2,2,Enu(n_Enu),eL),Pab(3,2,Enu(n_Enu),eL),
     #                  Pab(1,3,Enu(n_Enu),eL),Pab(2,3,Enu(n_Enu),eL),
     #                  Pab(3,3,Enu(n_Enu),eL)
                   endDO
                      CLOSE(Nfilof)
                CASE(3)
                      DO n_f1=1,Nfl
                      !n_f1= 1
                        DO n_f2=1,Nfl
                        !n_f2= 1
                          filon=OutDir//'probanu/3D/'//DMn(n_DM)//
     #                    NTn(nNTn)//san(n_sa)//hin(n_hi)//'_'//
     #                    fln(n_f1)//fln(n_f2)//ext
                          OPEN(Nfilof,FILE=filon)
                          WRITE(Nfilof,301) zero,cos(zth)
                          DO n_Enu=1,NEnu
                            WRITE(Nfilof,302) lgEnu(n_Enu)
                            DO n_zth=1,Nzth
                              root=sqrt((ra+re)**2-re**2*
     #                                             (sin(zth(n_zth)))**2)
                              bite=re*cos(zth(n_zth))
                              eL=root-bite
                              WRITE(Nfilof,302)
     #                                      Pab(n_f1,n_f2,Enu(n_Enu),eL)
                         endDO
                            WRITE(Nfilof,*)
                       endDO
                          CLOSE(Nfilof)
                     endDO
                   endDO
       endSELECT
         
           CALL cpu_time(finish_time)
           time=finish_time-start_time
           WRITE(*,*) time
      endDO
      endDO
      endDO

  300 FORMAT(10(1PE16.8))
  301 FORMAT(501(1PE16.8))
  302 FORMAT(1PE16.8$)

         STOP 'THE END OF PROGRAM vacuum3'
      END PROGRAM vacuum3
