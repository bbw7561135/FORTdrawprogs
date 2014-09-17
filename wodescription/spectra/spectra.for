************************************************************************
      PROGRAM spectra
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/10/15  *
************************************************************************
 
         USE PhysMathConstants, ONLY: zero,one,ten
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfl  =   2,
     #                Nt   =   2,
     #                NC_sp=  13,                                       !2,6,10 or 13
     #                NE_sp= 100,
     #                NC_za= 100,
     #                NE_za=   6
                 REAL,PARAMETER::
     #                C_min=-1.0d+00,
     #                C_max= 1.0d+00,
     #                E_min= 1.0d-01,
     #                E_max= 1.0d+11
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*3
     #                Sp
         CHARACTER*1
     #                SA,
     #                fln(Nfl)/'e','m'/,
     #                NTn(Nt)/'n','a'/
                 REAL
     #                C2(2)  /zero,one/,
     #                C6(6)  /-1.00,-0.25,-0.05, 0.05, 0.25, 1.00/,
     #                C10(10)/ 0.05, 0.15, 0.25, 0.35, 0.45, 0.55,
     #                                          0.65, 0.75, 0.85, 0.95/,
     #                C13(13)/-1.00,-0.75,-0.50,-0.25,-0.15,-0.05,
     #                   0.00, 0.05, 0.15, 0.25, 0.50, 0.75, 1.00/,
     #                Cosine(NC_sp),zacos(NC_za),Energy(NE_za),
     #                PhiE3(NC_sp),Ratio(NC_za),dFv
     
         lgEmin =log10(E_min)
         lgEmax =log10(E_max)
         steplgE=(lgEmax-lgEmin)/(NE_sp-1)
         stepC  =(C_max-C_min)/(NC_za-1)

         DO n_NE=1,NE_za
           Energy(n_NE)=ten**(n_NE-2)
      endDO

         SELECTCASE(NC_sp)
               CASE(2)
                     DO n_NC=1,NC_sp
                       Cosine(n_NC)=C2(n_NC)
                  endDO
               CASE(6)
                     DO n_NC=1,NC_sp
                       Cosine(n_NC)=C6(n_NC)
                  endDO
               CASE(10)
                     DO n_NC=1,NC_sp
                       Cosine(n_NC)=C10(n_NC)!-
                  endDO
               CASE(13)
                     DO n_NC=1,NC_sp
                       Cosine(n_NC)=C13(n_NC)
                  endDO
               CASE DEFAULT
                     STOP 'ERROR IN PROGRAM spectra: WRONG NC_sp'
      endSELECT

         DO n_NC=1,NC_za
           zacos(n_NC)=C_min+(n_NC-1)*stepC
      endDO
         
         N  = 1                                                          (Honda11 spectrum)
         set=dFAN_dE(one,one,N)

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

         DO n_fl=1,Nfl
           DO nNTn=1,Nt
             OPEN(Ndata,FILE=Out//'specexp/'//Sp//SA//'sp'//
     #                                        fln(n_fl)//NTn(nNTn)//ext)
             WRITE(Ndata,101) zero,Cosine
             DO n_NE=1,NE_sp
               E =ten**(lgEmin+(n_NE-1)*steplgE)
               E3=E**3
               DO n_NC=1,NC_sp
                 PhiE3(n_NC)=dFnu_dE(n_fl,nNTn,E,Cosine(n_NC),N)*E3
            endDO
               WRITE(Ndata,101) E,PhiE3
          endDO
             CLOSE(Ndata)
        endDO
      endDO

         DO n_fl=1,Nfl
           DO nNTn=1,Nt
             OPEN(Ndata,FILE=Out//'specexp/'//Sp//SA//'za'//
     #                                        fln(n_fl)//NTn(nNTn)//ext)
             WRITE(Ndata,102) zero,zacos
             DO n_NE_za=1,NE_za
               E  =Energy(n_NE_za)
               dFv=dFnu_dE(n_fl,nNTn,E,one,N)
               DO n_NC=1,NC_za
                 Ratio(n_NC)=dFnu_dE(n_fl,nNTn,E,zacos(n_NC),N)/dFv
            endDO
               WRITE(Ndata,102) E,Ratio
          endDO
             CLOSE(Ndata)
        endDO
      endDO

  101 FORMAT(14(1PE16.8))
  102 FORMAT(101(1PE16.8))

         STOP 'THE END OF PROGRAM spectra'
      END PROGRAM spectra