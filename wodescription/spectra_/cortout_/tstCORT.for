************************************************************************
      PROGRAM tstCORT
************************************************************************

      INCLUDE 'Include.f'

cc    -----------------------------------------------------------------
cc    ATTENTION! Before running the program, user MUST specify the
cc               following 4 variables:
cc    -----------------------------------------------------------------

      INTEGER*4    ExpV/5/       ! Neutrino experiment number
      CHARACTER*3  SA/'Min'/    ! Level of solar activity
      INTEGER*4    PNmod/0/     ! Prompt Neutrino production Model
      CHARACTER*43 Output/'/home/redponick/phd/fortran/Output/CORTout/'/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------------------- c
c     Variable 'Output'                                                c
c     ---------------------------------------------------------------- c
c                                                                      c
c     Variable 'Output' defines the path the output table 'Test.dat'   c
c     If user does not want to see this table the first STOP operator  c
c     in this demonstration program must be uncommented.               c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc    -----------------------------------------------------------------
cc    ATTENTION! Before running the program, user also MUST define 5
cc               variables ('Tab', 'Out', 'Test', 'Res', and 'B')
cc               described in the subroutine 'CORTout'. 
cc               In fact, only variable 'Tab' is important (it defines
cc               the path to the AN flux tables), while 4 others can
cc               be lived with their default values.
cc    -----------------------------------------------------------------

      COMMON /InOut/Enu,cosZA,Fnu_e1,Fnu_e2,Fnu_m1,Fnu_m2

      CALL CORTout(ExpV,SA,PNmod)    ! Initializing the subroutine

c     -----------------------------------------------------------------
c
c     After one call of the subroutine 'CORTout' for given values
c     of ExpV, SA, and PNmod, one can repeatedly call the subroutine
c     (entry) 'Flux' (by DO loops, for example), changing values of
c     (anti)neutrino energy and cosine of zenith angle.
c
c     NOTE: cosZA = -1, 0, and +1 are for upward-going, horizontal,
c           and downward-going neutrinos, respectively.
c
c     -----------------------------------------------------------------

      Enu   =  2.500D10              ! Neutrino energy [GeV]
      cosZA =  0.000D0              ! Cosine of zenith angle

      CALL Flux(PNmod)              ! Getting the result

c     -----------------------------------------------------------------
c
c     In the example below, the output is printed into the screen.
c     Clearly, user can (re)organize the output for his/her own needs.
c
c     -----------------------------------------------------------------

      PRINT *, ' '
      PRINT *, ' Prompt Neutrino model (PNmod)',PNmod
      PRINT *, ' Neutrino energy [GeV]:       ',Enu
      PRINT *, ' Cosine of zenith angle:      ',cosZA
      PRINT *, ' '

c     Output of the differential energy spectra

      PRINT *, ' -------------------------'
      PRINT *, ' dN/dE [1/(m**2 s sr GeV)]'
      PRINT *, ' -------------------------'
      PRINT *, ' Electron neutrinos:      ',Fnu_e1
      PRINT *, ' Electron antineutrinos:  ',Fnu_e2
      PRINT *, ' Muon neutrinos:          ',Fnu_m1
      PRINT *, ' Muon antineutrinos:      ',Fnu_m2
      PRINT *, '                          '

c     -----------------------------------------------------------------
c
c     The next example is for an output of the AN differential energy
c     spectra as a function of energy at a given zenith angle into a
c     file [neutrino energy range is 50 MeV to 25 EeV with equidistant
c     values of LOG10(Energy)].
c
c     STOP                   ! uncomment if you don't want this output
c
c     -----------------------------------------------------------------

      cosZA=0.0D0            ! Cosine of zenith angle
      NE=601                 ! Number of neutrino energies for output
      LgEmin=DLOG10(0.05D00) ! 50 MeV is the minimum available energy
      LgEmax=DLOG10(2.50D10) ! 25 EeV is an energy near the GZK cutoff

      OPEN (1,FILE=Output//'CORTest.dat')

      StLgE=(LgEmax-LgEmin)/(NE-1)

         DO i=1,NE
         Log10E=LgEmin+StLgE*(i-1)
         Enu=Ten**Log10E                  ! Neutrino energy [GeV] 
         CALL Flux(PNmod)                 ! Getting the result
         WRITE(1,10) Enu,Fnu_e1,Fnu_e2,Fnu_m1,Fnu_m2
      endDO

      CLOSE(1)

      STOP

   10 FORMAT(5(1pE11.4))

      END
