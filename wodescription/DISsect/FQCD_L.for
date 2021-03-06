************************************************************************
      FUNCTION FQCD_L_set(x,Q2)
************************************************************************
*                                                                      *
*                                                                      *
*                                      ITEP, Moscow, Russia 2005/06/16 *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants

         IMPLICIT REAL (A-H,K-M,O-Z), INTEGER (I,J,N)

         SAVE

         LOGICAL(2),PARAMETER::
     #              RESET = .TRUE.,
     #              TEST  = .FALSE.,                                     Test for spline interpolation
     #              Quiz  = .TRUE.                                       Coeff2 settings
            INTEGER,PARAMETER::
     #              Nx    = 100,                                         Number of points for x values
     #              NQ2   = 100,                                         Number of points for Q^2 values for Spline
     #              Mx    = 200,
     #              MQ2   = 100,                                         Number of points for Q^2 for test
     #              NxQ   = (Nx+2)*(NQ2+2),                              Number of points for Coeff2
     #              Mode  =   1,                                         Interpolation mode
     #              L     =   1,                                         Coeff2 settings
     #              MinCal= 105
               REAL,PARAMETER::
     #              x_ini = 1.0d-07,
     #              x_fin = one-10*Precision,
     #              Q2_ini= 1.0d-01,
     #              Q2_fin= 1.0d+05,
     #              RelErr= 1.0d-05

         COMMON      /n_NT/n_NT                                          Switch for neutrino type
         COMMON      /n_TT/n_TT                                          Switch for nuclear target type
         COMMON  /n_AG_DIS/n_AG_DIS                                      Switch for model of DIS structure functions
         COMMON  /n_FL_DIS/n_FL_DIS                                      Switch for model of function F_L
         COMMON  /n_RT_DIS/n_RT_DIS                                      Switch for model of function R
         COMMON  /n_Rc_DIS/n_Rc_DIS                                      Switch for model of modification of function R
         COMMON    /PDFLIB/NGROUP,NSET                                   Parameters for PDFLIB setup
         COMMON     /m_ini/m_ini,mm_ini                                  Mass of target nucleon
         COMMON         /x/z                                             Bjorken scaling variable x
         COMMON        /Q2/QQ                                            Square of mometum transfer (Q^2=-q^2)
         COMMON      /E_nu/E_nu                                          Neutrino energy

         EXTERNAL FunGeM,GeM_FQCD_L

         DIMENSION
     #          Fmn_p(Nx,NQ2),Fmn_n(Nx,NQ2),Fma_p(Nx,NQ2),Fma_n(Nx,NQ2),
     #          Cmn_p(NxQ),   Cmn_n(NxQ),   Cma_p(NxQ),   Cma_n(NxQ),
     #          Rmn_p(NQ2),   Rmn_n(NQ2),   Rma_p(NQ2),   Rma_n(NQ2),
     #          Amn_p(MQ2),   Amn_n(MQ2),   Ama_p(MQ2),   Ama_n(MQ2)

         lgx_ini = log10(x_ini)
         lgx_fin = log10(x_fin)
         lgQ2_ini= log10(Q2_ini)
         lgQ2_fin= log10(Q2_fin)
         IF (RESET) THEN
           n_AG_DIS= 1
           n_FL_DIS= 0
           n_RT_DIS= 1
           n_Rc_DIS= 1
           NGROUP  = 6
           NSET    = 2

           CALL GeMSet(FunGeM,Res,zero,one,RelErr,MinCal,*100)
           set= R_set(one,one,one)
           CALL NuclStructFuns(one,one,one,one,one,one,one,one,one)

           OPEN(Ndat01,FILE=OutSF//'F_L/FQCD_L_mn_p_R.dat')
           OPEN(Ndat02,FILE=OutSF//'F_L/FQCD_L_mn_n_R.dat')
           OPEN(Ndat03,FILE=OutSF//'F_L/FQCD_L_ma_p_R.dat')
           OPEN(Ndat04,FILE=OutSF//'F_L/FQCD_L_ma_n_R.dat')
           OPEN(Ndat11,FILE=OutSF//'F_L/FQCD_L_mn_p.dat')
           OPEN(Ndat12,FILE=OutSF//'F_L/FQCD_L_mn_n.dat')
           OPEN(Ndat13,FILE=OutSF//'F_L/FQCD_L_ma_p.dat')
           OPEN(Ndat14,FILE=OutSF//'F_L/FQCD_L_ma_n.dat')

           E_nu= 1.00d+01                                                TEST VALUE

           steplgx = (lgx_fin -lgx_ini )/(Nx -1)
           steplgQQ= (lgQ2_fin-lgQ2_ini)/(NQ2-1)
           DO i=1,Nx
             z =ten**(lgx_ini+(i-1)*steplgx)
             m= z**2/(2*pi)
             DO j=1,NQ2
               QQ  = ten**(lgQ2_ini+(j-1)*steplgQQ)

               n_TT= 1; m_ini= m_p; mm_ini= mm_p
               n_NT=+1; CALL GeMInt(GeM_FQCD_L,F_mn_p,z,one,*100)
               n_NT=-1; CALL GeMInt(GeM_FQCD_L,F_ma_p,z,one,*100)

               n_TT= 2; m_ini= m_n; mm_ini= mm_n
               n_NT=+1; CALL GeMInt(GeM_FQCD_L,F_mn_n,z,one,*100)
               n_NT=-1; CALL GeMInt(GeM_FQCD_L,F_ma_n,z,one,*100)

               Fmn_p(i,j)= m*F_mn_p; Rmn_p(j)= m*F_mn_p
               Fmn_n(i,j)= m*F_mn_n; Rmn_n(j)= m*F_mn_n
               Fma_p(i,j)= m*F_ma_p; Rma_p(j)= m*F_ma_p
               Fma_n(i,j)= m*F_ma_n; Rma_n(j)= m*F_ma_n
          endDO
             WRITE(Ndat01,101) z, (Rmn_p(j),   j=1,NQ2)
             WRITE(Ndat02,101) z, (Rmn_n(j),   j=1,NQ2)
             WRITE(Ndat03,101) z, (Rma_p(j),   j=1,NQ2)
             WRITE(Ndat04,101) z, (Rma_n(j),   j=1,NQ2)
             WRITE(Ndat11,101)    (Fmn_p(i,j), j=1,NQ2)
             WRITE(Ndat12,101)    (Fmn_n(i,j), j=1,NQ2)
             WRITE(Ndat13,101)    (Fma_p(i,j), j=1,NQ2)
             WRITE(Ndat14,101)    (Fma_n(i,j), j=1,NQ2)
        endDO
           CALL GeMInf
                    ELSE
           OPEN(Ndat11,FILE=OutSF//'F_L/FQCD_L_mn_p.dat')
           OPEN(Ndat12,FILE=OutSF//'F_L/FQCD_L_mn_n.dat')
           OPEN(Ndat13,FILE=OutSF//'F_L/FQCD_L_ma_p.dat')
           OPEN(Ndat14,FILE=OutSF//'F_L/FQCD_L_ma_n.dat')
           DO i=1,Nx
             READ(Ndat11,101) (Fmn_p(i,j), j=1,NQ2)
             READ(Ndat12,101) (Fmn_n(i,j), j=1,NQ2)
             READ(Ndat13,101) (Fma_p(i,j), j=1,NQ2)
             READ(Ndat14,101) (Fma_n(i,j), j=1,NQ2)
        endDO
      endIF
         CLOSE(Ndat01); CLOSE(Ndat02); CLOSE(Ndat03); CLOSE(Ndat04)
         CLOSE(Ndat11); CLOSE(Ndat12); CLOSE(Ndat13); CLOSE(Ndat14)

         CALL Coeff2(Mode,11,Nx,NQ2,lgx_ini,lgQ2_ini,lgx_fin,lgQ2_fin,
     #        Fmn_p,Cmn_p,Quiz,L)
         CALL Coeff2(Mode,12,Nx,NQ2,lgx_ini,lgQ2_ini,lgx_fin,lgQ2_fin,
     #        Fmn_n,Cmn_n,Quiz,L)
         CALL Coeff2(Mode,13,Nx,NQ2,lgx_ini,lgQ2_ini,lgx_fin,lgQ2_fin,
     #        Fma_p,Cma_p,Quiz,L)
         CALL Coeff2(Mode,14,Nx,NQ2,lgx_ini,lgQ2_ini,lgx_fin,lgQ2_fin,
     #        Fma_n,Cma_n,Quiz,L)

         IF (TEST) THEN
           OPEN(Ndat01,FILE=OutSF//'F_L/FQCD_L_mn_p_A.dat')
           OPEN(Ndat02,FILE=OutSF//'F_L/FQCD_L_mn_n_A.dat')
           OPEN(Ndat03,FILE=OutSF//'F_L/FQCD_L_ma_p_A.dat')
           OPEN(Ndat04,FILE=OutSF//'F_L/FQCD_L_ma_n_A.dat')
           steplgx = (lgx_fin -lgx_ini )/(Mx -1)
           steplgQQ= (lgQ2_fin-lgQ2_ini)/(MQ2-1)
           DO i=1,Mx
             lgz= lgx_ini+(i-1)*steplgx
             z  = ten**lgz
             DO j=1,MQ2
               lgQQ= lgQ2_ini+(j-1)*steplgQQ
               QQ  = ten**lgQQ
               Amn_p(j)=Sp2(11,Cmn_p,lgz,lgQQ)
               Amn_n(j)=Sp2(12,Cmn_n,lgz,lgQQ)
               Ama_p(j)=Sp2(13,Cma_p,lgz,lgQQ)
               Ama_n(j)=Sp2(14,Cma_n,lgz,lgQQ)
          endDO
             WRITE(Ndat01,102) z, (Amn_p(j), j=1,MQ2)
             WRITE(Ndat02,102) z, (Amn_n(j), j=1,MQ2)
             WRITE(Ndat03,102) z, (Ama_p(j), j=1,MQ2)
             WRITE(Ndat04,102) z, (Ama_n(j), j=1,MQ2)
        endDO
           CLOSE(Ndat01); CLOSE(Ndat02); CLOSE(Ndat03); CLOSE(Ndat04)
      endIF
         FQCD_L_set= one

         RETURN
  100    STOP 'Error with GeM in FUNCTION FQCD_L_set'

  101 FORMAT(1PE9.3,101(1PE10.3))
  102 FORMAT(1PE9.3,101(1PE11.3))

c 101 FORMAT(1PE9.3,300(1PE10.3))                                        GNU FORTRAN Compiler
c 102 FORMAT(1PE9.3,300(1PE11.3))                                        GNU FORTRAN Compiler

*     ==================================================================
      ENTRY FQCD_L(x,Q2)
*     ==================================================================
*        IF (x  < x_ini ) PRINT *, 'x  =', x,  'x_ini  =', x_ini
*        IF (x  > x_fin ) PRINT *, 'x  =', x,  'x_fin  =', x_fin
*        IF (Q2 < Q2_ini) PRINT *, 'Q2 =', Q2, 'Q2_ini =', Q2_ini
*        IF (Q2 > Q2_fin) PRINT *, 'Q2 =', Q2, 'Q2_fin =', Q2_fin
         SELECTCASE(n_TT)
               CASE(   1)                                                PROTON
               SELECTCASE(n_NT)
                     CASE(   1); FQCD_L=Sp2(11,Cmn_p,log10(x),log10(Q2)) NEUTRINO
                     CASE(  -1); FQCD_L=Sp2(13,Cma_p,log10(x),log10(Q2)) ANTINEUTRINO
            endSELECT
               CASE(   2)                                                NEUTRON
               SELECTCASE(n_NT)
                     CASE(   1); FQCD_L=Sp2(12,Cmn_n,log10(x),log10(Q2)) NEUTRINO
                     CASE(  -1); FQCD_L=Sp2(14,Cma_n,log10(x),log10(Q2)) ANTINEUTRINO
            endSELECT
      endSELECT
         RETURN
*     ==================================================================

      END FUNCTION FQCD_L_set