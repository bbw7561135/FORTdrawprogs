************************************************************************
      PROGRAM spCORTo
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2012/11/06  *
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
!     #                NC_ex   =   2,
     #                NC_ex   =  13,
     #                NE_ex   =   6,
     #                NC      = 100,
     #                NE      = 100
                 REAL,PARAMETER::
     #                C_min=-1.0d+00,
     #                C_max= one,
     #                E_nu_min= 1.0d-01,
     #                E_nu_max= 1.0d+04
                 REAL
!     #                Carr(NC_ex)/zero,one/,
     #                Carr(NC_ex)/-1.00,-0.75,-0.5,-0.25,-0.15,-0.05,
     #                            zero,0.05,0.15,0.25,0.5,0.75,1.00/,
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
         
         N  = 2                                                          (CORT spectrum)
         
         set=dFACN_dE(one,one,N)
         
         OPEN(Nfspen,FILE=Out//'CRTnspen.dat')
         OPEN(Nfspea,FILE=Out//'CRTnspea.dat')
         OPEN(Nfspmn,FILE=Out//'CRTnspmn.dat')
         OPEN(Nfspma,FILE=Out//'CRTnspma.dat')
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

         DO n_NC=1,NC
           CarrE(n_NC)=C_min+(n_NC-1)*stepC
      endDO
         OPEN(Nfzaen,FILE=Out//'CRTnzaen.dat')
         OPEN(Nfzaea,FILE=Out//'CRTnzaea.dat')
         OPEN(Nfzamn,FILE=Out//'CRTnzamn.dat')
         OPEN(Nfzama,FILE=Out//'CRTnzama.dat')
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

         STOP 'THE END OF PROGRAM spCORTo'
      END PROGRAM spCORTo
