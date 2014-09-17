************************************************************************
*                                                                      *
*            CTEQ PARTON DISTRIBUTION FUNCTIONS: VERSION 6             *
*                                                                      *
*                       April    10, 2002, v6.01                       *
*                       February 23, 2003, v6.10                       *
*                       August   06, 2003, v6.11                       *
*                       December 12, 2004, v6.12                       *
*                       December 04, 2006, v6.50                       *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] J. Pumplin,  D.R. Stump,  J. Huston, H.L. Lai, P. Nadol-     *
*         sky, and W.K. Tung,  "New generation of parton distribu-     *
*         tions  with  uncertainties  from global  QCD  analysis,"     *
*         JHEP 207 (2002) 012 [arXiv:hep-ph/0201195].                  *
*     [2] D. Stump, J. Huston,  J. Pumplin,  W.K. Tung,  H.L. Lai,     *
*         S. Kuhlmann,  and J. Owens,  "Inclusive  jet production,     *
*         parton distributions,  and the  search for new physics,"     *
*         JHEP 310 (2003) 046 [arXiv:hep-ph/0303013].                  *
*     [3] F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky,     *
*         H.L. Lai, S. Kretzer, J.F. Owens, and W.K. Tung,  "Neut-     *
*         rinodimuon  production and  strangeness asymmetry of the     *
*         nucleon," Eur. Phys. J. C 40 (2005) 145-156  [arXiv:hep-     *
*         ph/0312323].                                                 *
*     [4] S. Kretzer, H.L. Lai, F. Olness, and  W.K. Tung,  "CTEQ6     *
*         parton distributions  with  heavy quark  mass  effects,"     *
*         Phys. Rev. D 69 (2004) 114005 [arXiv:hep-ph/0307022].        *
*     [5] W.K. Tung, H.L. Lai,  A. Belyaev,  J. Pumplin, D. Stump,     *
*         and C.-P. Yuan,  "Heavy quark mass effects in deep ine-      *
*         lastic  scattering and  global QCD  analysis,"  J. High      *
*         Ener.Phys. 02 (2007) 053 [arXiv: hep-ph/0611254].            *
*                                                                      *
*     This package contains                                            *
*     (1) 4 standard sets of CTEQ6 PDF's  (CTEQ6M, CTEQ6D, CTEQ6L,     *
*         CTEQ6L1);                                                    *
*     (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty     *
*         studies from Ref.[1];                                        *
*     (3) updated version  of the above:  CTEQ6.1M  and its 40 up/     *
*         down eigenvector sets from Ref.[2];                          *
*     (4) 5 special sets for strangeness study from Ref.[3];           *
*     (5) 1 special set for heavy quark study from Ref.[4];            *
*     (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref.[5].   *
*                                                                      *
*     Details about calling convention are:                            *
*     ============================================================     *
*     Iset PDF-set Description    a_s(Mz)** Lam4  Lam5  Table File     *
*     ============================================================     *
*     Standard, "best-fit", sets, Ref.[1]                              *
*     ------------------------------------------------------------     *
*      1   CTEQ6M  Standard MSbar 0.118     326   226   cteq6m.tbl     *
*      2   CTEQ6D  Standard DIS   0.118     326   226   cteq6d.tbl     *
*      3   CTEQ6L  Leading Order  0.118**   326** 226   cteq6l.tbl     *
*      4   CTEQ6L1 Leading Order  0.130**   215** 165  cteq6l1.tbl     *
*     ------------------------------------------------------------     *
*     Special sets for strangeness study, Ref.[3]                      *
*     ------------------------------------------------------------     *
*     11   CTEQ6A  Class A        0.118     326   226  cteq6sa.pds     *
*     12   CTEQ6B  Class B        0.118     326   226  cteq6sb.pds     *
*     13   CTEQ6C  Class C        0.118     326   226  cteq6sc.pds     *
*     14   CTEQ6B+ Large [S-]     0.118     326   226 cteq6sb+.pds     *
*     15   CTEQ6B- Negative [S-]  0.118     326   226 cteq6sb-.pds     *
*     ------------------------------------------------------------     *
*     Special set for Heavy Quark study, Ref.[4]                       *
*     ------------------------------------------------------------     *
*     21   CTEQ6HQ                0.118     326   226  cteq6hq.pds     *
*     ------------------------------------------------------------     *
*     For uncertainty calculations using eigenvectors of the Hess-     *
*     ian: central +40up/down sets along 20 eigenvector directions     *
*     ------------------------------------------------------------     *
*     Original version,  Ref.[1]: central fit: CTEQ6M (=CTEQ6M.00)     *
*     ------------------------------------------------------------     *
*     1xx  CTEQ6M.xx  +/- sets    0.118     326   226 cteq6m1xx.tbl    *
*     ------------------------------------------------------------     *
*     where xx=01-40: 01/02 corresponds to +/-  for the 1st eigen-     *
*     vector, ... etc. e.g. 100 is CTEQ6M.00(=CTEQ6M), 101/102 are     *
*     CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.              *
*     ------------------------------------------------------------     *
*     Updated version, Ref.[2]:  central fit: CTEQ6.1M(=CTEQ61.00)     *
*     ------------------------------------------------------------     *
*     2xx CTEQ61.xx  +/- sets     0.118     326   226 ctq61.xx.tbl     *
*     ------------------------------------------------------------     *
*     where xx=01-40: 01/02 corresponds to +/- for  the 1st eigen-     *
*     vector, ... etc. e.g. 200 is  CTEQ61.00(=CTEQ6.1M),  201/202     *
*     are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.          *
*     ------------------------------------------------------------     *
*     Version with mass effects, Ref.[5]: central fit: CTEQ6.5M        *
*     (=CTEQ65.00)                                                     *
*     ------------------------------------------------------------     *
*     3xx CTEQ65.xx  +/- sets     0.118     326   226 ctq65.xx.pds     *
*     ------------------------------------------------------------     *
*     where xx=01-40: 01/02 corresponds to +/- for the 1st eigen-      *
*     vector, ... etc. e.g. 300 is CTEQ65.00(=CTEQ6.5M),  301/302      *
*     are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.          *
*     ============================================================     *
*     ----------                                                       *
*     ** ALL fits are obtained by using the same coupling strength     *
*        \alpha_s(Mz)=0.118 and the NLO running \alpha_s  formula,     *
*     except CTEQ6L1  which uses the LO running \alpha_s  and  its     *
*     value determined from the fit. For the LO fits, the evoluti-     *
*     on of the PDF and the hard cross sections are  calculated at     *
*     LO. More detailed  discussions are given  in the references.     *
*                                                                      *
*     The  table  grids  are  generated  for   10^-6 < x < 1   and     *
*     1.3 < Q < 10^4 GeV.  PDF values  outside of  the above range     *
*     are  returned  using  extrapolation. Lam5 (Lam4)  represents     *
*     Lambda value (in MeV) for 5(4) flavors. The matching alpha_s     *
*     between 4 and 5 flavors  takes place at Q=4.5 GeV,  which is     *
*     defined as the bottom quark mass, whenever it can be applied.    *
*                                                                      *
*     The Table_Files  are assumed to be in the working directory.     *
*     Before using the PDF, it is necessary  to do the initializa-     *
*     tion  by CALL SETCTQ6(Iset)  where Iset is  the desired  PDF     *
*     specified in the above table.                                    *
*     The function CTQ6PDF(Iparton,X,Q) returns the parton distri-     *
*     bution inside  the proton for parton [Iparton]  at [X] Bjor-     *
*     ken_x and scale [Q] (GeV) in PDF set [Iset]. Iparton  is the     *
*     parton label:                                                    *
*     (5, 4, 3, 2, 1,  0,  -1,    -2,....-3,    -4,    -5   )  for     *
*     (b, c, s, d, u,  g,  u_bar, d_bar, s_bar, c_bar, b_bar).         *
*     For detailed information on the parameters used,  e.q. quark     *
*     masses, QCD Lambda, etc., see info lines at the beginning of     *
*     the Table_Files. These programs, as provided,  are in double     *
*     precision. By removing the "IMPLICIT DOUBLEPRECISION" lines,     *
*     they can also be run in single precision.                        *
*                                                                      *
*     If you have detailed questions concerning these  CTEQ6 dist-     *
*     ributions, or if you find  problems/bugs using this package,     *
*     direct  inquires  to  Pumplin@pa.msu.edu or Tung@pa.msu.edu.     *
*                                                                      *
************************************************************************
      SUBROUTINE SETCTQ6(Iset)
************************************************************************

         USE InpOutUnits

         IMPLICIT REAL (A-H,O-Z)

         CHARACTER(*),PARAMETER::
     #   Tab=IniPP//'CTEQ6_04.12.2006/'                                  KK

         PARAMETER(Isetmax0=7)

         CHARACTER Flnm(Isetmax0)*6, nn*3, Tablefile*116

         DATA (Flnm(I),I=1,Isetmax0)
     #   /'cteq6m','cteq6d','cteq6l','cteq6l','ctq61.','cteq6s',
     #    'ctq65.'/
         DATA Isetold, Isetmin0,Isetmin1,Isetmax1 /-987,1,100,140/
         DATA Isetmin2,Isetmax2 /200,240/
         DATA Isetmin3,Isetmax3 /300,340/
         DATA IsetminS,IsetmaxS /11,15/
         DATA IsetHQ            /21/

         COMMON    /Valence/MxVal
         COMMON  /Setchange/Isetch

         SAVE

         IF (Iset.ne.Isetold) THEN
*          ----------------------------------------------------------- *
*          IF DATA FILE NOT INITIALIZED, DO SO.                        *
*          ----------------------------------------------------------- *
           MxVal=2
           IU   =NextUn()
           IF (Iset.ge.Isetmin0 .and. Iset.le.3       ) THEN             Iset=1,2,3 for 6m, 6d, 6l
             Tablefile=Tab//Flnm(Iset)//'.tbl'
       ELSEIF (Iset.eq.4                              ) THEN             4 (2nd LO fit)
             Tablefile=Tab//Flnm(Iset)//'1.tbl'
       ELSEIF (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) THEN             101-140
             WRITE(nn,'(I3)') Iset
             Tablefile=Tab//Flnm(1)//nn//'.tbl'
       ELSEIF (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) THEN             200-240
             WRITE(nn,'(I3)') Iset
             Tablefile=Tab//Flnm(5)//nn(2:3)//'.tbl'
       ELSEIF (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) THEN              11- 15
             MxVal=3
             IF (Iset.eq.11) THEN; Tablefile=Flnm(6)//'a.pds'
         ELSEIF (Iset.eq.12) THEN; Tablefile=Flnm(6)//'b.pds'
         ELSEIF (Iset.eq.13) THEN; Tablefile=Flnm(6)//'c.pds'
         ELSEIF (Iset.eq.14) THEN; Tablefile=Flnm(6)//'b+.pds'
         ELSEIF (Iset.eq.15) THEN; Tablefile=Flnm(6)//'b-.pds'
          endIF
       ELSEIF (Iset.eq.IsetHQ                         ) THEN              21
             MxVal=3
             TableFile='cteq6hq.pds'
       ELSEIF (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) THEN             (Cteq6.5)  300-340
             MxVal=3
             WRITE(nn,'(I3)') Iset
             Tablefile=Tab//Flnm(7)//nn(2:3)//'.pds'
c             print *, Tablefile
                                                        ELSE
             PRINT *, 'Invalid Iset number:', Iset
             STOP
        endIF
           OPEN(IU,File=Tablefile,Status='OLD',Err=100)
   21      CALL ReadTbl(IU)
           CLOSE(IU)
           Isetold=Iset
           Isetch =1
      endIF

         RETURN
  100    PRINT *,'DATA file',Tablefile,'cannot be opened!'
         STOP
      END SUBROUTINE SETCTQ6
