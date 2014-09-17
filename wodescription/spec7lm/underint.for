************************************************************************
      FUNCTION underint(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE OscMatParameters

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
         
                 REAL
     #                cosin,Theta,EneV,
     #                P(3,3)

         COMMON      /N/N
         COMMON  /C_lim/C_min,deltaC
         COMMON   /E_nu/E_nu
         COMMON      /P/P
         COMMON   /n_NT/n_NT

         underint=one
         
         RETURN
         
*     ==================================================================
      ENTRY uifen(var)
*     ==================================================================
         n_NT= 1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uifen=P(1,1)*dFAen_dE(E_nu,cosin,N)+
     #           P(2,1)*dFAmn_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uifen=dFAen_dE(E_nu,cosin,N)                                 !*deltaC
        endIF
         RETURN

*     ==================================================================
      ENTRY uifea(var)
*     ==================================================================
         n_NT=-1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uifea=P(1,1)*dFAea_dE(E_nu,cosin,N)+
     #           P(2,1)*dFAma_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uifea=dFAea_dE(E_nu,cosin,N)                                 !*deltaC
        endIF
        
         RETURN

*     ==================================================================
      ENTRY uifmn(var)
*     ==================================================================
         n_NT= 1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uifmn=P(2,2)*dFAmn_dE(E_nu,cosin,N)+
     #           P(1,2)*dFAen_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uifmn=dFAmn_dE(E_nu,cosin,N)                                 !*deltaC
        endIF
         
         RETURN

*     ==================================================================
      ENTRY uifma(var)
*     ==================================================================
         n_NT=-1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uifma=P(2,2)*dFAma_dE(E_nu,cosin,N)+
     #           P(1,2)*dFAea_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uifma=dFAma_dE(E_nu,cosin,N)                                 !*deltaC
        endIF
         
         RETURN

*     ==================================================================
      ENTRY uiftn(var)
*     ==================================================================
         n_NT= 1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uiftn=P(1,3)*dFAen_dE(E_nu,cosin,N)+
     #           P(2,3)*dFAmn_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uiftn=zero                                                   !*deltaC
        endIF
         
         RETURN

*     ==================================================================
      ENTRY uifta(var)
*     ==================================================================
         n_NT=-1
         cosin=deltaC*var+C_min
         Theta=acos(cosin)
         EneV =1.0d+09*E_nu
         IF (cosin.LE.zero) THEN
           CALL Pmatrix(EneV,Theta)
           uifta=P(1,3)*dFAea_dE(E_nu,cosin,N)+
     #           P(2,3)*dFAma_dE(E_nu,cosin,N)                          !*deltaC
                            ELSE
           uifta=zero                                                   !*deltaC
        endIF
         
         RETURN

      END FUNCTION underint
