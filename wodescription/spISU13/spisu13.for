************************************************************************
      PROGRAM spisu13
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/09/02  *
************************************************************************
 
         USE PhysMathConstants, ONLY: zero,one,ten
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

              INTEGER,PARAMETER::
     #                Nfspen  = 140,
     #                Nfspea  = 141,
     #                Nfzaen  = 142,
     #                Nfzaea  = 143,
     #                Nfspmn  = 144,
     #                Nfspma  = 145,
     #                Nfzamn  = 146,
     #                Nfzama  = 147,
     #                NC_ex   =   2,
!     #                NC_ex   =  10,
!     #                NC_ex   =  13,
!     #                NC_ex   =   6,
     #                NE_ex   =   6,
     #                NC      = 100,
     #                NE      = 100
                 REAL,PARAMETER::
     #                C_min=-1.0d+00,
     #                C_max= one,
     #                E_nu_min= 1.0d-01,
     #                E_nu_max= 1.0d+04
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*57
     #                fnen,fnea,fnmn,fnma
         CHARACTER*3
     #                Sp
         CHARACTER*1
     #                SA,
     #                fln(2)/'e','m'/,
     #                NTn(2)/'n','a'/
                 REAL
     #                Carr(NC_ex)/zero,one/,
!     #                Carr(NC_ex)/-0.05,-0.15,-0.25,-0.35,-0.45,-0.55,
!     #                            -0.65,-0.75,-0.85,-0.95/,
!     #                Carr(NC_ex)/0.05,0.15,0.25,0.35,0.45,0.55,
!     #                            0.65,0.75,0.85,0.95/,
!     #                Carr(NC_ex)/-1.00,-0.75,-0.5,-0.25,-0.15,-0.05,
!     #                            zero,0.05,0.15,0.25,0.5,0.75,one/,
!     #                Carr(NC_ex)/-1.00,-0.25,-0.05,0.05,0.25,one/,
     #                CarrE(NC)/NC*zero/,
     #                Earr(NE_ex)/1.0d-01,1.0d+00,1.0d+01,1.0d+02,
     #                         1.0d+03,1.0d+04/,
     #                PhienE3(NC_ex)/NC_ex*zero/, 
     #                PhieaE3(NC_ex)/NC_ex*zero/,
     #                PhimnE3(NC_ex)/NC_ex*zero/, 
     #                PhimaE3(NC_ex)/NC_ex*zero/,
     #                Ratiomn(NC)/NC*zero/,
     #                Ratioma(NC)/NC*zero/,
     #                Ratioen(NC)/NC*zero/,
     #                Ratioea(NC)/NC*zero/
     
         lnE_nu_min=log10(E_nu_min)
         lnE_nu_max=log10(E_nu_max)
         steplnE_nu=(lnE_nu_max-lnE_nu_min)/(NE-1)
         stepC     =(C_max-C_min)/(NC-1)

         DO n_NC=1,NC
           CarrE(n_NC)=C_min+(n_NC-1)*stepC
      endDO
         
         N  = 1                                                          (Honda11 spectrum)
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

         fnen=Out//'spISU13/'//Sp//SA//'sp'//fln(1)//NTn(1)//ext
         fnea=Out//'spISU13/'//Sp//SA//'sp'//fln(1)//NTn(2)//ext
         fnmn=Out//'spISU13/'//Sp//SA//'sp'//fln(2)//NTn(1)//ext
         fnma=Out//'spISU13/'//Sp//SA//'sp'//fln(2)//NTn(2)//ext
         OPEN(Nfspen,FILE=fnen)
         OPEN(Nfspea,FILE=fnea)
         OPEN(Nfspmn,FILE=fnmn)
         OPEN(Nfspma,FILE=fnma)

         WRITE(Nfspen,101) zero,Carr
         WRITE(Nfspea,101) zero,Carr
         WRITE(Nfspmn,101) zero,Carr
         WRITE(Nfspma,101) zero,Carr
         DO n_NE=1,NE
           E_nu=ten**(lnE_nu_min+(n_NE-1)*steplnE_nu)
           E3  =E_nu**3
           DO n_NC=1,NC_ex
             PhienE3(n_NC)=dFAen_dE(E_nu,Carr(n_NC),N)*E3
             PhieaE3(n_NC)=dFAea_dE(E_nu,Carr(n_NC),N)*E3
             PhimnE3(n_NC)=dFAmn_dE(E_nu,Carr(n_NC),N)*E3
             PhimaE3(n_NC)=dFAma_dE(E_nu,Carr(n_NC),N)*E3
        endDO
           WRITE(Nfspen,101) E_nu,PhienE3
           WRITE(Nfspea,101) E_nu,PhieaE3
           WRITE(Nfspmn,101) E_nu,PhimnE3
           WRITE(Nfspma,101) E_nu,PhimaE3
      endDO
         CLOSE(Nfspen)
         CLOSE(Nfspea)
         CLOSE(Nfspmn)
         CLOSE(Nfspma)

         fnen=Out//'spISU13/'//Sp//SA//'za'//fln(1)//NTn(1)//ext
         fnea=Out//'spISU13/'//Sp//SA//'za'//fln(1)//NTn(2)//ext
         fnmn=Out//'spISU13/'//Sp//SA//'za'//fln(2)//NTn(1)//ext
         fnma=Out//'spISU13/'//Sp//SA//'za'//fln(2)//NTn(2)//ext
         OPEN(Nfzaen,FILE=fnen)
         OPEN(Nfzaea,FILE=fnea)
         OPEN(Nfzamn,FILE=fnmn)
         OPEN(Nfzama,FILE=fnma)

         WRITE(Nfzaen,102) zero,CarrE
         WRITE(Nfzaea,102) zero,CarrE
         WRITE(Nfzamn,102) zero,CarrE
         WRITE(Nfzama,102) zero,CarrE
         DO n_NE_ex=1,NE_ex
           E_nu=Earr(n_NE_ex)
           denoren=dFAen_dE(E_nu,one,N)
           denorea=dFAea_dE(E_nu,one,N)
           denormn=dFAmn_dE(E_nu,one,N)
           denorma=dFAma_dE(E_nu,one,N)
           DO n_NC=1,NC
             Ratioen(n_NC)=dFAen_dE(E_nu,CarrE(n_NC),N)/denoren
             Ratioea(n_NC)=dFAea_dE(E_nu,CarrE(n_NC),N)/denorea
             Ratiomn(n_NC)=dFAmn_dE(E_nu,CarrE(n_NC),N)/denormn
             Ratioma(n_NC)=dFAma_dE(E_nu,CarrE(n_NC),N)/denorma
        endDO
           WRITE(Nfzaen,102) E_nu,Ratioen
           WRITE(Nfzaea,102) E_nu,Ratioea
           WRITE(Nfzamn,102) E_nu,Ratiomn
           WRITE(Nfzama,102) E_nu,Ratioma
      endDO
         CLOSE(Nfzaen)
         CLOSE(Nfzaea)
         CLOSE(Nfzamn)
         CLOSE(Nfzama)

!  101 FORMAT(3(1PE16.8))
  101 FORMAT(14(1PE16.8))
  102 FORMAT(101(1PE16.8))

         STOP 'THE END OF PROGRAM spISU13'
      END PROGRAM spisu13
