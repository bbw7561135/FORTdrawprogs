************************************************************************
      MODULE OscParameters
************************************************************************
*                                                                      *
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/08/21  *
************************************************************************

         USE PhysMathConstants, ONLY: zero,half,pi

              REAL
     #             s2th12,s2th13,s2th23,                                 sin^2(mixing angle)
     #             dm21,dm13,dm32,delm,                                  Neutrino mass squared splittings [eV^2]!*1.0d-18
     #             CPVdel,                                               CP-violation phase
     #             dm(3,3)/9*zero/,
     #             o12(3,3)/9*zero/,
     #             o23(3,3)/9*zero/,
     #             o31(3,3)/9*zero/
              INTEGER
     #             n_hi,n_sa

               COMMON /n_sa/n_sa
               COMMON /n_hi/n_hi

      CONTAINS
        SUBROUTINE OscPar

          SELECTCASE(n_sa)
                CASE(1)                                                  Fogli et al.
                      s2th12= 3.07d-01
                      dm21  = 7.54d-05
                      IF (n_hi.EQ.1) THEN                                Normal neutrino mass hierarchy
                        s2th13  = 2.41d-02
                        s2th23  = 3.86d-01
                        delm    =-2.43d-03
                        CPVdel  = 1.08d+00*pi
                                     ELSE                                Inverse neutrino mass hierarchy
                        s2th13  = 2.44d-02
                        s2th23  = 3.92d-01
                        delm    = 2.42d-03
                        CPVdel  = 1.09d+00*pi
                   endIF
                      dm13=-delm-half*dm21
                CASE(2)                                                  Valle et al.
                      s2th12= 3.20d-01
                      dm21  = 7.62d-05
                      IF (n_hi.EQ.1) THEN                                Normal neutrino mass hierarchy
                        s2th13  = 2.46d-02
                        s2th23  = 6.13d-01
                        dm13    =-2.55d-03
                        CPVdel  = 8.00d-01*pi
                                     ELSE                                Inverse neutrino mass hierarchy
                        s2th13  = 2.50d-02
                        s2th23  = 6.00d-01
                        dm13    = 2.43d-03
                        CPVdel  =-3.00d-02*pi
                   endIF
       endSELECT

          dm32   = -(dm21+dm13)
          dm(1,2)= - dm21
          dm(1,3)=   dm13
          dm(2,1)=   dm21
          dm(2,3)= - dm32
          dm(3,1)= - dm13
          dm(3,2)=   dm32
         
          th12=asin(sqrt(s2th12))
          th23=asin(sqrt(s2th23))
          th13=asin(sqrt(s2th13))
         
          o12(1,1)= cos(th12); o12(1,2)= sin(th12)
          o12(2,1)=-sin(th12); o12(2,2)= cos(th12)                       Solar matrix
          o12(3,3)= 1

          o23(1,1)= 1
          o23(2,2)= cos(th23); o23(2,3)= sin(th23)                       Atmospheric matrix
          o23(3,2)=-sin(th23); o23(3,3)= cos(th23)

          o31(1,1)= cos(th13); o31(1,3)= sin(th13)
          o31(2,2)= 1                                                    Reactor matrix
          o31(3,1)=-sin(th13); o31(3,3)= cos(th13)

        END SUBROUTINE OscPar
      END MODULE OscParameters
************************************************************************
