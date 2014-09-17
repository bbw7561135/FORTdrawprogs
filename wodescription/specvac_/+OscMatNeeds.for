************************************************************************
      MODULE OscMatNeeds
************************************************************************
*                                                                      *
*                                            Tula, Russia, 2012/11/15  *
*                                BLTP JINR, Dubna, Russia, 2013/08/21  *
************************************************************************
*     ---------------------------------------------------------------- *
*     PATHS TO INPUT/OUTPUT FOLDERS                                    *
*     ---------------------------------------------------------------- *
      CHARACTER(*),PARAMETER::	
     #             InpDir   ='/home/redponick/phd/fortran/Input/',
     #             OutDir   ='/home/redponick/phd/fortran/Output/'
*     ---------------------------------------------------------------- *
              REAL,PARAMETER::     
*     ---------------------------------------------------------------- *
*       NUMERICAL CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             zero     = 0.00d+00,
     #             one      = 1.00d+00,
     #             two      = 2.00d+00,
     #             three    = 3.00d+00,
     #             four     = 4.00d+00,
     #             five     = 5.00d+00,
     #             six      = 6.00d+00,
     #             ten      = 1.00d+01,
     #             half     = 5.00d-01,                                  1/2
     #             third    = one/three,                                 1/3
     #             sixth    = one/six,                                   1/6
     #             quarter  = 2.50d-01,                                  1/4
     #             sqrt2    = 1.414213562373095048801688724209698d+00,   sqrt(2)
     #             sqrt3    = 1.732050807568877293527446341505872d+00,   sqrt(3)
     #             sqrt6    = 2.449489742783178098197284074705892d+00,   sqrt(6)
     #             pi       = 3.141592653589793238462643383279503d+00,   pi constant
     #             twopi    = two*pi,                                    2*pi
     #             halfpi   = half*pi,                                   pi/2
     #             Precision= epsilon(one)
           COMPLEX,PARAMETER::
     #             i        = (zero,one)
              REAL,PARAMETER:: 
*     ---------------------------------------------------------------- *
*        PHYSICAL CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             hbarc    = 1.973269681718d-14,                        Reduced Plank constant [eV m] (PDG12)
*     ---------------------------------------------------------------- *
*           EARTH CONSTANTS                                            *
*     ---------------------------------------------------------------- *
     #             R_Earth  = 6.371d+06,                                 Earth radius [m]
*     ---------------------------------------------------------------- *
*        SPECIFIC PARAMETERS                                           *
*     ---------------------------------------------------------------- *
     #             DetDep   = zero,                                      Detector depression
     #             RadDet   = R_Earth-DetDep,
     #             ThetaMin = halfpi,                                    Lower Angle 
     #             ThetaMax = pi,                                        Upper Angle
     #             E_min    = 1.00d+07,                                  Lower Energy [eV]
     #             E_max    = 1.00d+11                                   Upper Energy [eV]
           INTEGER,PARAMETER::
     #             Nfl      = 3,
     #             Nt       = 2,
     #             Nhi      = 2,
     #             Nsa      = 2
      END MODULE OscMatNeeds
************************************************************************
