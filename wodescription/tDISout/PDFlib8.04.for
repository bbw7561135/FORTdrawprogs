************************************************************************
      SUBROUTINE PDFlib8_04set
************************************************************************
*                                                                      *
*     ------------------------------------------------------------     *
*     "PDFLIB  Proton,  Pion and  Photon Parton Density Functions,     *
*     Parton Density Functions of the Nucleus, and $\alpha_s$ Cal-     *
*     culations". Users's Manual Version 8.04.                         *
*     ------------------------------------------------------------     *
*                                                                      *
*     ATTENTION !   These parameters should be set  with a call to     *
*                   subroutine PDFSET at the initialisation phase.     *
*     NPTYPE  is  number of particle type: 1 - protons, 2 - pions,     *
*     3 - photons.  NGROUP is number of author group ranging  from     *
*     1 to 9. NSET is number of a selected  structure function set     *
*     within the author group  ranging  from 1 to 58.  TMAS is the     *
*     optional user defined value of the top­quark mass in GeV/c^2     *
*     XMIN, XMAX  are  minimal,  maximal  allowed x values. Q2MIN,     *
*     Q2MAX are  minimal, maximal allowed Q^2 values in (GeV/c)^2.     *
*                                                                      *
*     ATTENTION !   NPTYPE = number of  particle type ranging from     *
*                   1 to 3 (Protons: NPTYPE=1); NGROUP = number of     *
*     author group ranging from 1 to 9;  NSET = number of a selec-     *
*     ted structure  function set  within the  author  group  ran-     *
*     ging from 1 to 58;  NFL = desired number of flavours  in the     *
*     $\alpha_s$ calculation ranging from 3 to 6 (Default: NFL=5);     *
*     LO = order of $\alpha_s$ s calculation, if LO = 1 $\alpha_s$     *
*     is calculated to first order only (Default: LO=2) TMAS = the     *
*     user  defined  value  of  the  top­quark  mass  in   GeV/c^2     *
*     (optional)  (Default:  TMAS=180.0d+00);  QCDL4 = QCD  scale;     *
*     \Lambda^[4]_QCD, in GeV for four flavours QCDL5 = QCD scale,     *
*     \Lambda^[5]_QCD, in GeV for five flavours  corresponding to      *
*     QCDL4 and XMIN = minimal allowed x value; XMAX = maximal al-     *
*     lowed x value; Q2MIN = minimal allowed Q^2 value (in (GeV/c)     *
*     ^2); Q2MAX = maximal allowed Q^2 value (in (GeV/c) 2).           *
*                                                                      *
*     ATTENTION !   The user has to provide the  following INPUTs:     *
*                   X = x value of parton,  SCALE = QCD  scale  in     *
*     GeV.  The  subroutine  STRUCTM returns the following OUTPUT:     *
*     UPV = up valence quark, DNV = down valence,  USEA =  sea up,     *
*     DSEA = sea down, STR = strange,  CHM = charm,  BOT = bottom,     *
*     TOP = top, GL = gluon.  In  case  up is not given separately     *
*     from down it is set USEA = DSEA.                                 *
*                                                                      *
*     ATTENTION !   PLEASE NOTE THAT IN ANY OF THE CALLING SEQUEN-     *
*                   CES FOR THE PROTON, THE PION AND THE PHOTON IT     *
*                   IS ALWAYS RETURNED x*PDF !                         *
*                                                                      *
*     ATTENTION !   TO  RUN  PDFLIB  A  LINK  TO THE  CERN LIBRARY     *
*                   (PACKLIB,  MATHLIB AND KERNLIB)  IS  REQUIRED.     *
*                                                                      *
*           Some sets of parton distributions for NPTYPE = 1.          *
*     ------------------------------------------------------------     *
*      NGROUP | NSET | XMIN | XMAX | Q2MIN | Q2MAX | NAME OF SET       *
*     ------------------------------------------------------------     *
*         3      99    10^-5    1    1.25             MRST-cdn         *
*         4      46    10^-5    1    1.00             CTEQ5-LO         *
*         5      14    10^-9    1    0.80             GRV98-HO-DIS     *
*     ------------------------------------------------------------     *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: m_t

         CHARACTER(20) PARAM(20)
               REAL VALUE(20)
            INTEGER NGROUP,NSET

         COMMON     /PDFLIB/NGROUP,NSET                                  Parameters for PDFLIB setup

         PARAM( 1)='NPTYPE'; VALUE( 1)= 1
         PARAM( 2)='NGROUP'; VALUE( 2)= NGROUP
         PARAM( 3)='NSET'  ; VALUE( 3)= NSET
         PARAM( 4)='QCDL4' ; VALUE( 4)= 0.20000d+00
         PARAM( 5)='QCDL5' ; VALUE( 5)= 0.17080d+00
         PARAM( 6)='XMIN'  ; VALUE( 6)= 0.00000d+00
         PARAM( 7)='XMAX'  ; VALUE( 7)= 0.99999d+00
         PARAM( 8)='Q2MIN' ; VALUE( 8)= 0.00000d+00
         PARAM( 9)='Q2MAX' ; VALUE( 9)= 1.00000d+07
         PARAM(10)='TMAS'  ; VALUE(10)= m_t
         PARAM(11)='NFL'   ; VALUE(11)= 4

         CALL PDFSET(PARAM,VALUE)
         CALL PDFSTA

         RETURN
      END SUBROUTINE PDFlib8_04set