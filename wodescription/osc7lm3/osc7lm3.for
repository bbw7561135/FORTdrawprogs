************************************************************************
      PROGRAM osc7lm3
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/05  *
************************************************************************

         USE OscMatParameters

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfilo  = 100,
     #                Nfdt   = 116,
     #                NEnu   = 50,
     #                Nzth   = 50
                 REAL,PARAMETER::
     #                Enu_min=E_min,
     #                Enu_max=E_max,
     #                zth_min=ThetaMin,
     #                zth_max=ThetaMax
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*54
     #                filon
         CHARACTER*3
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/
              INTEGER
     #                n_DM,n_hi,n_NT,nNTn
                 REAL
     #                P(3,3),MaxD(3,3),
     #                Enu(NEnu),zth(Nzth),
     #                start_time,finish_time,time
     
               COMMON /MaxD/    MaxD
               COMMON /P/       P
               COMMON /n_NT/    n_NT                                     Switch for neutrino type
     
         CALL cpu_time(start_time)

         n_DM= 3

         n_hi= 1

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

         DO nNTn=1,2
           IF (nNTn.EQ.1) THEN
             n_NT= 1
                          ELSE
             n_NT=-1
        endIF

           CALL OscMat3

           filon=OutDir//'probanu/'//DM(n_DM)//'3'//NTn(nNTn)//hin(n_hi)
     #                                                            //ext
           OPEN(Nfilo,FILE=filon)
         
           IF (UTestSw) THEN 
             OPEN(Nfdt,FILE=OutDir//'2dt'//ext)
        endIF
           DO n_NEnu=1,NEnu
             DO n_Nzth=1,Nzth
               CALL Pmatrix(Enu(NEnu),zth(Nzth))
               WRITE(Nfilo,101) Enu(NEnu),cos(zth(Nzth)),P
               IF (UTestSw) THEN
                 WRITE(Nfdt,*) MaxD
            endIF
          endDO
        endDO
           CLOSE(Nfilo)
           IF (UTestSw) THEN
             CLOSE(Nfdt)
        endIF
         
           CALL cpu_time(finish_time)
           time=finish_time-start_time
           WRITE(*,*) time
      endDO
         
  101 FORMAT(11(1PE16.8))

         STOP 'THE END OF PROGRAM osc7lm3'
      END PROGRAM osc7lm3
