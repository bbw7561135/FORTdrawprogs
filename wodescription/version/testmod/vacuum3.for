************************************************************************
      PROGRAM vacuum3
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/08/21  *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one,ten,pi

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER(*),PARAMETER::
     #                Out    = '/home/redponick/phd/fortran/Output/'
              INTEGER,PARAMETER::
     #                Nfilo  = 100,
     #                NEnu   = 50,
     #                Nzth   = 50
                 REAL,PARAMETER::
     #                re     = 6.371d+03,                                Earth radius!*1.0d+05
     #                ra     = zero,                                     Atmosphere depth!*1.0d+05
     #                Enu_min= 1.00d-01,
     #                Enu_max= 1.00d+01,
     #                zth_min= 0.5d+00*pi,
     #                zth_max= pi
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*55
     #                filon
         CHARACTER*3
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                san(2)/'F','V'/,
     #                hin(2)/'n','i'/,
     #                NTn(2)/'n','a'/
              INTEGER
     #                n_DM,n_hi,n_sa,n_NT,nNTn
                 REAL
     #                Enu(NEnu),zth(Nzth),
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

         DO n_NEnu=1,NEnu
           Enu(n_NEnu)=ten**(lgE_min+(n_NEnu-1)*steplgE)
      endDO
         DO n_Nzth=1,Nzth
           zth(n_Nzth)=zth_min+(n_Nzth-1)*stepzth
      endDO

!         DO n_sa=1,2
         n_sa= 1
!         DO n_hi=1,2
         n_hi= 1
         DO nNTn=1,2
           IF (nNTn.EQ.1) THEN
             n_NT= 1
                          ELSE
             n_NT=-1
        endIF

           Palphabeta=oscvac3(1,1,one,one)
         
           filon=Out//'probanu/'//DM(n_DM)//'3'//NTn(nNTn)//san(n_sa)
     #                                                 //hin(n_hi)//ext
           OPEN(Nfilo,FILE=filon)

           DO n_NEnu=1,NEnu
             DO n_Nzth=1,Nzth
               root=sqrt((ra+re)**2-re**2*(sin(zth(n_Nzth)))**2)
               bite=re*cos(zth(n_Nzth))
               eL=root-bite
               WRITE(Nfilo,101) Enu(n_NEnu),cos(zth(n_Nzth)),
     #                 Pab(1,1,Enu(n_NEnu),eL),Pab(1,2,Enu(n_NEnu),eL),
     #                 Pab(1,3,Enu(n_NEnu),eL),Pab(2,1,Enu(n_NEnu),eL),
     #                 Pab(2,2,Enu(n_NEnu),eL),Pab(2,3,Enu(n_NEnu),eL),
     #                 Pab(3,1,Enu(n_NEnu),eL),Pab(3,2,Enu(n_NEnu),eL),
     #                 Pab(3,3,Enu(n_NEnu),eL)
          endDO
        endDO
           CLOSE(Nfilo)
!        endDO
!        endDO
         
           CALL cpu_time(finish_time)
           time=finish_time-start_time
           WRITE(*,*) time
      endDO
         
  101 FORMAT(11(1PE16.8))

         STOP 'THE END OF PROGRAM vacuum3'
      END PROGRAM vacuum3
