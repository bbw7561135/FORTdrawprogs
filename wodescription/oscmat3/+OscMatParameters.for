************************************************************************
      MODULE OscMatParameters
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2013/02/25  *
************************************************************************
*     ---------------------------------------------------------------- *
*     PATHS TO INPUT/OUTPUT FOLDERS                                    *
*     ---------------------------------------------------------------- *
      CHARACTER(*),PARAMETER::	
     #             InpDir='/home/redponick/phd/fortran/Input/',
     #             OutDir='/home/redponick/phd/fortran/Output/'
*     ---------------------------------------------------------------- *
              REAL,PARAMETER::     
*     ---------------------------------------------------------------- *
*       NUMERICAL CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             zero    = 0.00d+00,
     #             one     = 1.00d+00,
     #             two     = 2.00d+00,
     #             three   = 3.00d+00,
     #             four    = 4.00d+00,
     #             five    = 5.00d+00,
     #             six     = 6.00d+00,
     #             ten     = 1.00d+01,
     #             half    = 5.00d-01,                                   1/2
     #             third   = one/three,                                  1/3
     #             sixth   = one/six,                                    1/6
     #             quarter = 2.50d-01,                                   1/4
     #             sqrt2   = 1.414213562373095048801688724209698d+00,    sqrt(2)
     #             sqrt3   = 1.732050807568877293527446341505872d+00,    sqrt(3)
     #             sqrt6   = 2.449489742783178098197284074705892d+00,    sqrt(6)
     #             pi      = 3.141592653589793238462643383279503d+00,    pi constant
     #             twopi   = two*pi,                                     2*pi
     #             halfpi  = half*pi,                                    pi/2
     #             Precision= epsilon(one)
           COMPLEX,PARAMETER::
     #             i       = (zero,one)
              REAL,PARAMETER:: 
*     ---------------------------------------------------------------- *
*        PHYSICAL CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             hbarc   = 1.9732696817182680069447d-07,               Reduced Plank constant [GeV cm]
*     ---------------------------------------------------------------- *
*           EARTH CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             R_Earth = 6.371d+06,                                  Earth radius [m]
     #             Coeff   = 2.543333333333333333333333d-14,             !WTF?!
*     ---------------------------------------------------------------- *
*     OSCILLATION PARAMETERS                                           *
*     ---------------------------------------------------------------- *
!      INCLUDE 'OscParValle.f'
      INCLUDE 'OscParFogli.f'
     #             dm32    =-(dm21+dm13),
*     ---------------------------------------------------------------- *
*        SPECIFIC PARAMETERS                                           *
*     ---------------------------------------------------------------- *
     #             DetDep  = zero,                                       Detector depression
     #             RadDet  = R_Earth-DetDep,
     #             ThetaMin= halfpi,                                     Lower Angle 
     #             ThetaMax= pi,                                         Upper Angle
     #             E_min   = 1.00d+08,                                   Lower Energy [eV]
     #             E_max   = 1.00d+09,                                   Upper Energy [eV]
     #             VarDens = one                                         Coefficient for the Earth density variation
           INTEGER,PARAMETER::
     #             N1      = 1,
     #             N2      = 1,
     #             N3      = 1,
     #             N4      = 1,
     #             N5      = 1,
     #             N6      = 1,
     #             N7      = 1,
     #             N8      = 1,
     #             N9      = 1,
     #             N10     = 1,
     #             KN      = N1+N2+N3+N4+N5+N6+N7+N8+N9+N10              Total number of Earth sublayers
*     ---------------------------------------------------------------- *
*                 SWITCHES                                             *
*     ---------------------------------------------------------------- *
        LOGICAL(2),PARAMETER::
     #             UTestSw =.FALSE.                                      Switch for unitarity test
      END MODULE OscMatParameters
************************************************************************
