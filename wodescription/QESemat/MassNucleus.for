************************************************************************
      FUNCTION MassNucleus(Z,A)
************************************************************************
*                                                                      *
*     This FUNCTION  returns the mass of the  nucleus with  atomic     *
*     number A and charge Z.                                           *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] Physics Values, Moscow, 1991, 1232 pp.,  edited  by I.S.     *
*         Grigorev and E.Z. Melikhov,  chapter 39 "Nuclear reacti-     *
*         ons", Table 39.1, pp. 1069-1085.                             *
*     [2] A.H. Wapstra and K. Bos, "The 1977 atomic mass evaluati-     *
*         on: in four parts part I. Atomic mass table," Atom. Dat.     *
*         Nucl. Dat. Tab. 19 (1977) 177-214.                           *
*     [3] G. Audi,  A.H. Wapstra,  and  C. Thibault,  "The AME2003     *
*         atomic mass evaluation (II). Tables, graphs and referen-     *
*         ces," Nucl. Phys. A 729 (2003) 337-676.                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,UAMU

         IMPLICIT REAL (A-Z)

          D=zero
          IF (A.eq.  6.0) THEN
            IF (Z.eq.  1.0) THEN; D=+41860.00000                         H  (HYDROGEN)
        ELSEIF (Z.eq.  2.0) THEN; D=+17595.10000                         He (HELIUM)
        ELSEIF (Z.eq.  3.0) THEN; D=+14086.79300                         Li (LITHIUM)
        ELSEIF (Z.eq.  4.0) THEN; D=+18375.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+43600.00000                         B  (BORON)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq.  7.0) THEN
            IF (Z.eq.  1.0) THEN; D=+49140.00000                         H  (HYDROGEN)
        ELSEIF (Z.eq.  2.0) THEN; D=+26101.00000                         He (HELIUM)
        ELSEIF (Z.eq.  3.0) THEN; D=+14908.14000                         Li (LITHIUM)
        ELSEIF (Z.eq.  4.0) THEN; D=+15770.03000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+27870.00000                         B  (BORON)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 11.0) THEN
            IF (Z.eq.  3.0) THEN; D=+40797.00000                         Li (LITHIUM)
        ELSEIF (Z.eq.  4.0) THEN; D=+20174.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+ 8667.90000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+10650.30000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+24300.00000                         N  (NITROGEN)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 12.0) THEN
            IF (Z.eq.  3.0) THEN; D=+50100.00000                         Li (LITHIUM)
        ELSEIF (Z.eq.  4.0) THEN; D=+25077.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+13368.90000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+00000.00000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+17338.10000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+32048.00000                         O  (OXIGEN)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 13.0) THEN
            IF (Z.eq.  4.0) THEN; D=+33250.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+16562.20000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+ 3125.01130                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+ 5345.48000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+23112.00000                         O  (OXIGEN)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 14.0) THEN
            IF (Z.eq.  4.0) THEN; D=+39950.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+23664.00000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+ 3019.89300                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+ 2863.41700                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+ 8007.36000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+32660.00000                         F  (FLUORINE)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 15.0) THEN
            IF (Z.eq.  4.0) THEN; D=+49800.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+28972.00000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+ 9873.10000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+  101.43800                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+ 2855.60000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+16780.00000                         F  (FLUORINE)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 16.0) THEN
            IF (Z.eq.  4.0) THEN; D=+57680.00000                         Be (BERYLLIUM)
        ELSEIF (Z.eq.  5.0) THEN; D=+37080.00000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+13694.00000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+ 5683.70000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=- 4737.00141                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+10680.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+23996.00000                         Ne (NEON)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 18.0) THEN
            IF (Z.eq.  5.0) THEN; D=+52320.00000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+24930.00000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+13114.00000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=-  781.50000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+  873.70000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+ 5317.17000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=+24190.00000                         Na (SODIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 19.0) THEN
            IF (Z.eq.  5.0) THEN; D=+59360.00000                         B  (BORON)
        ELSEIF (Z.eq.  6.0) THEN; D=+32420.00000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+15862.00000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+ 3334.90000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=- 1487.39000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+ 1751.44000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=+12927.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=+33040.00000                         Mg (MAGNESIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 20.0) THEN
            IF (Z.eq.  6.0) THEN; D=+37560.00000                         C  (CARBON)
        ELSEIF (Z.eq.  7.0) THEN; D=+21770.00000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+ 3797.50000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=-   17.40000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=- 7041.93130                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=+ 6848.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=+17570.00000                         Mg (MAGNESIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 23.0) THEN
            IF (Z.eq.  7.0) THEN; D=+38400.00000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+14610.00000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+ 3330.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=- 5154.05000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=- 9529.85360                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=- 5473.80000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=+ 6770.00000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=+23770.00000                         Si (SILICON)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 24.0) THEN
            IF (Z.eq.  7.0) THEN; D=+47540.00000                         N  (NITROGEN)
        ELSEIF (Z.eq.  8.0) THEN; D=+19070.00000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+ 7560.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=- 5951.50000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=- 8418.11000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=-13933.56700                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=-   56.90000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=+10755.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=+32000.00000                         P  (PHOSPHORUS)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 26.0) THEN
            IF (Z.eq.  8.0) THEN; D=+35710.00000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+18270.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+  430.00000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=- 6862.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=-16214.58200                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=-12210.31000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=- 7145.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=+10970.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=+25970.00000                         S  (SULFUR)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 27.0) THEN
            IF (Z.eq.  8.0) THEN; D=+44950.00000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+24930.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+ 7070.00000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=- 5517.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=-14586.65000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=-17196.66000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=-12384.30000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=-  717.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=+17540.00000                         S  (SULFUR)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 28.0) THEN
            IF (Z.eq.  8.0) THEN; D=+53850.00000                         O  (OXIGEN)
        ELSEIF (Z.eq.  9.0) THEN; D=+33230.00000                         F  (FLUORINE)
        ELSEIF (Z.eq. 10.0) THEN; D=+11240.00000                         Ne (NEON)
        ELSEIF (Z.eq. 11.0) THEN; D=-  989.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=-15018.60000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=-16850.44000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=-21492.79680                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=- 7159.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=+ 4070.00000                         S  (SULFUR)
        ELSEIF (Z.eq. 17.0) THEN; D=+26560.00000                         Cl (CHLORINE)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 35.0) THEN
            IF (Z.eq. 11.0) THEN; D=+39580.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=+16150.00000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=-  130.00000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=-14360.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=-24857.70000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=-28846.36000                         S  (SULFUR)
        ELSEIF (Z.eq. 17.0) THEN; D=-29013.54000                         Cl (CHLORINE)
        ELSEIF (Z.eq. 18.0) THEN; D=-23047.40000                         Ar (ARGON)
        ELSEIF (Z.eq. 19.0) THEN; D=-11169.00000                         K  (POTASSIUM)
        ELSEIF (Z.eq. 20.0) THEN; D=+ 4600.00000                         Ca (CALCIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 36.0) THEN
            IF (Z.eq. 11.0) THEN; D=+47950.00000                         Na (SODIUM)
        ELSEIF (Z.eq. 12.0) THEN; D=+21420.00000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=+ 5780.00000                         Al (ALUMINIUM)
        ELSEIF (Z.eq. 14.0) THEN; D=-12480.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=-20251.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=-30664.07000                         S  (SULFUR)
        ELSEIF (Z.eq. 17.0) THEN; D=-29521.86000                         Cl (CHLORINE)
        ELSEIF (Z.eq. 18.0) THEN; D=-30231.54000                         Ar (ARGON)
        ELSEIF (Z.eq. 19.0) THEN; D=-17426.00000                         K  (POTASSIUM)
        ELSEIF (Z.eq. 20.0) THEN; D=- 6440.00000                         Ca (CALCIUM)
        ELSEIF (Z.eq. 21.0) THEN; D=+13900.00000                         Sc (SCANDIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 39.0) THEN
            IF (Z.eq. 12.0) THEN; D=+43570.00000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=+21400.00000                         Al (ALUMINIUM)
        ELSEIF (Z.eq. 14.0) THEN; D=+ 1930.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=-12870.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=-23160.00000                         S  (SULFUR)
        ELSEIF (Z.eq. 17.0) THEN; D=-29800.20000                         Cl (CHLORINE)
        ELSEIF (Z.eq. 18.0) THEN; D=-33242.00000                         Ar (ARGON)
        ELSEIF (Z.eq. 19.0) THEN; D=-33807.01000                         K  (POTASSIUM)
        ELSEIF (Z.eq. 20.0) THEN; D=-27274.40000                         Ca (CALCIUM)
        ELSEIF (Z.eq. 21.0) THEN; D=-14168.00000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=+ 1500.00000                         Ti (TITANIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 40.0) THEN
            IF (Z.eq. 12.0) THEN; D=+50240.00000                         Mg (MAGNESIUM)
        ELSEIF (Z.eq. 13.0) THEN; D=+29300.00000                         Al (ALUMINUM)
        ELSEIF (Z.eq. 14.0) THEN; D=+ 5470.00000                         Si (SILICON)
        ELSEIF (Z.eq. 15.0) THEN; D=- 8110.00000                         P  (PHOSPHORUS)
        ELSEIF (Z.eq. 16.0) THEN; D=-22870.00000                         S  (SULFUR)
        ELSEIF (Z.eq. 17.0) THEN; D=-27560.00000                         Cl (CHLORINE)
        ELSEIF (Z.eq. 18.0) THEN; D=-35039.89600                         Ar (ARGON)
        ELSEIF (Z.eq. 19.0) THEN; D=-33535.20000                         K  (POTASSIUM)
        ELSEIF (Z.eq. 20.0) THEN; D=-34846.27000                         Ca (CALCIUM)
        ELSEIF (Z.eq. 21.0) THEN; D=-20523.20000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=- 8850.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=+10330.00000                         V  (VANADIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 55.0) THEN
            IF (Z.eq. 19.0) THEN; D=-  270.00000                         K  (POTASSIUM)
        ELSEIF (Z.eq. 20.0) THEN; D=-18120.00000                         Ca (CALCIUM)
        ELSEIF (Z.eq. 21.0) THEN; D=-29580.00000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=-41670.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=-49150.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-55107.50000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-57710.60000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-57479.40000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-54027.60000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-45336.00000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-31620.00000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-14920.00000                         Zn (ZINC)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 56.0) THEN
            IF (Z.eq. 20.0) THEN; D=-13440.00000                         Ca (CALCIUM)
        ELSEIF (Z.eq. 21.0) THEN; D=-25270.00000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=-38940.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=-46080.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-55281.20000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-56909.70000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-60605.40000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-56039.40000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-53904.00000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-38600.00000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-25730.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=- 4740.00000                         Ga (GALLIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 58.0) THEN
            IF (Z.eq. 21.0) THEN; D=-15170.00000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=-30770.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=-40210.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-51830.00000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-55910.00000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-62153.40000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-59845.90000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-60227.70000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-51662.10000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-42300.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-23990.00000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=- 8370.00000                         Ge (GERMANIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 59.0) THEN
            IF (Z.eq. 21.0) THEN; D=-10040.00000                         Sc (SCANDIUM)
        ELSEIF (Z.eq. 22.0) THEN; D=-25220.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=-37070.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-47890.00000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-55480.00000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-60663.10000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-62228.40000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-61155.70000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-56357.20000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-47260.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-34120.00000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=-17000.00000                         Ge (GERMANIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 63.0) THEN
            IF (Z.eq. 22.0) THEN; D=- 5200.00000                         Ti (TITANIUM)
        ELSEIF (Z.eq. 23.0) THEN; D=-20910.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-35530.00000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-46350.00000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-55550.00000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-61840.00000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-65512.60000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-65579.50000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-62213.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-56547.10000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=-46910.00000                         Ge (GERMANIUM)
        ELSEIF (Z.eq. 33.0) THEN; D=-33820.00000                         As (ARSENIC)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 64.0) THEN
            IF (Z.eq. 23.0) THEN; D=-15400.00000                         V  (VANADIUM)
        ELSEIF (Z.eq. 24.0) THEN; D=-33150.00000                         Cr (CHROMIUM)
        ELSEIF (Z.eq. 25.0) THEN; D=-42620.00000                         Mn (MANGANESE)
        ELSEIF (Z.eq. 26.0) THEN; D=-54770.00000                         Fe (IRON)
        ELSEIF (Z.eq. 27.0) THEN; D=-59793.00000                         Co (COBALT)
        ELSEIF (Z.eq. 28.0) THEN; D=-67099.30000                         Ni (NICKEL)
        ELSEIF (Z.eq. 29.0) THEN; D=-65424.20000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-66003.60000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-58834.30000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=-54350.00000                         Ge (GERMANIUM)
        ELSEIF (Z.eq. 33.0) THEN; D=-39520.00000                         As (ARSENIC)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq. 79.0) THEN
            IF (Z.eq. 29.0) THEN; D=-42330.00000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-53420.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-62510.00000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=-69490.00000                         Ge (GERMANIUM)
        ELSEIF (Z.eq. 33.0) THEN; D=-73637.00000                         As (ARSENIC)
        ELSEIF (Z.eq. 34.0) THEN; D=-75917.60000                         Se (SELENIUM)
        ELSEIF (Z.eq. 35.0) THEN; D=-76068.50000                         Br (BROMINE)
        ELSEIF (Z.eq. 36.0) THEN; D=-74443.00000                         Kr (KRYPTON)
        ELSEIF (Z.eq. 37.0) THEN; D=-70803.00000                         Rb (RUBIDIUM)
        ELSEIF (Z.eq. 38.0) THEN; D=-65477.00000                         Sr (STRONTIUM)
        ELSEIF (Z.eq. 39.0) THEN; D=-58360.00000                         Y  (YTTRIUM)
        ELSEIF (Z.eq. 40.0) THEN; D=-47360.00000                         Zr (ZIRCONIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq. 80.0) THEN
            IF (Z.eq. 29.0) THEN; D=-36450.00000                         Cu (COPPER)
        ELSEIF (Z.eq. 30.0) THEN; D=-51840.00000                         Zn (ZINC)
        ELSEIF (Z.eq. 31.0) THEN; D=-59140.00000                         Ga (GALLIUM)
        ELSEIF (Z.eq. 32.0) THEN; D=-69515.00000                         Ge (GERMANIUM)
        ELSEIF (Z.eq. 33.0) THEN; D=-72159.00000                         As (ARSENIC)
        ELSEIF (Z.eq. 34.0) THEN; D=-77759.90000                         Se (SELENIUM)
        ELSEIF (Z.eq. 35.0) THEN; D=-75889.50000                         Br (BROMINE)
        ELSEIF (Z.eq. 36.0) THEN; D=-77892.50000                         Kr (KRYPTON)
        ELSEIF (Z.eq. 37.0) THEN; D=-72173.00000                         Rb (RUBIDIUM)
        ELSEIF (Z.eq. 38.0) THEN; D=-70308.00000                         Sr (STRONTIUM)
        ELSEIF (Z.eq. 39.0) THEN; D=-61220.00000                         Y  (YTTRIUM)
        ELSEIF (Z.eq. 40.0) THEN; D=-55520.00000                         Zr (ZIRCONIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq.118.0) THEN
            IF (Z.eq. 43.0) THEN; D=-45200.000                           Tc ()
        ELSEIF (Z.eq. 44.0) THEN; D=-57920.000                           Ru (RUTHENIUM)
        ELSEIF (Z.eq. 45.0) THEN; D=-65140.000                           Rh (RHODIUM)
        ELSEIF (Z.eq. 46.0) THEN; D=-75470.000                           Pd (PALLADIUM)
        ELSEIF (Z.eq. 47.0) THEN; D=-79570.000                           Ag (SILVER)
        ELSEIF (Z.eq. 48.0) THEN; D=-86709.000                           Cd (CADMIUM)
        ELSEIF (Z.eq. 49.0) THEN; D=-87230.000                           In (INDIUM)
        ELSEIF (Z.eq. 50.0) THEN; D=-91656.100                           Sn (TIN)
        ELSEIF (Z.eq. 51.0) THEN; D=-87999.000                           Sb (ANTIMONY)
        ELSEIF (Z.eq. 52.0) THEN; D=-87721.000                           Te (TELLURIUM)
        ELSEIF (Z.eq. 53.0) THEN; D=-80971.000                           I  (IODINE)
        ELSEIF (Z.eq. 54.0) THEN; D=-78079.000                           Xe (XENON)
        ELSEIF (Z.eq. 55.0) THEN; D=-68409.000                           Sc (CESIUM)
        ELSEIF (Z.eq. 56.0) THEN; D=-62370.000                           Ba (BARIUM)
        ELSEIF (Z.eq. 57.0) THEN; D=-49620.000                           La (LANTHANUM)
                            ELSE; GOTO 1111
         endIF

      ELSEIF (A.eq.119.0) THEN
            IF (Z.eq. 44.0) THEN; D=-53240.000                           Ru (RUTHENIUM)
        ELSEIF (Z.eq. 45.0) THEN; D=-63240.000                           Rh (RHODIUM)
        ELSEIF (Z.eq. 46.0) THEN; D=-71620.000                           Pd (PALLADIUM)
        ELSEIF (Z.eq. 47.0) THEN; D=-78560.000                           Ag (SILVER)
        ELSEIF (Z.eq. 48.0) THEN; D=-83910.000                           Cd (CADMIUM)
        ELSEIF (Z.eq. 49.0) THEN; D=-87704.000                           In (INDIUM)
        ELSEIF (Z.eq. 50.0) THEN; D=-90068.400                           Sn (TIN)
        ELSEIF (Z.eq. 51.0) THEN; D=-89477.000                           Sb (ANTIMONY)
        ELSEIF (Z.eq. 52.0) THEN; D=-87184.000                           Te (TELLURIUM)
        ELSEIF (Z.eq. 53.0) THEN; D=-83766.000                           I  (IODINE)
        ELSEIF (Z.eq. 54.0) THEN; D=-78794.000                           Xe (XENON)
        ELSEIF (Z.eq. 55.0) THEN; D=-72305.000                           Sc (CESIUM)
        ELSEIF (Z.eq. 56.0) THEN; D=-64590.000                           Ba (BARIUM)
        ELSEIF (Z.eq. 57.0) THEN; D=-54970.000                           La (LANTHANUM)
        ELSEIF (Z.eq. 58.0) THEN; D=-44000.000                           Ce (CERIUM)
                            ELSE; GOTO 1111
         endIF
*
*
*
      ELSEIF (A.eq.206.0) THEN
            IF (Z.eq. 80.0) THEN; D=-20946.000                           Hg (MERCURY)
        ELSEIF (Z.eq. 81.0) THEN; D=-22253.100                           Tl (THALLIUM)
        ELSEIF (Z.eq. 82.0) THEN; D=-23785.400                           Pb (LEAD)
        ELSEIF (Z.eq. 83.0) THEN; D=-20028.000                           Bi (BISMUTH)
        ELSEIF (Z.eq. 84.0) THEN; D=-18182.000                           Po (POLONIUM)
        ELSEIF (Z.eq. 85.0) THEN; D=-12420.000                           At (ASTATINE)
        ELSEIF (Z.eq. 86.0) THEN; D=- 9116.000                           Rn (RADON)
        ELSEIF (Z.eq. 87.0) THEN; D=- 1243.000                           Fr (FRANCIUM)
        ELSEIF (Z.eq. 88.0) THEN; D=+ 3565.000                           Ra (RADIUM)
        ELSEIF (Z.eq. 89.0) THEN; D=+13510.000                           Ac (ACTINIUM)
                            ELSE; GOTO 1111
         endIF
      ELSEIF (A.eq.207.0) THEN
            IF (Z.eq. 80.0) THEN; D=-16220.000                           Hg (MERCURY)
        ELSEIF (Z.eq. 81.0) THEN; D=-21034.000                           Tl (THALLIUM)
        ELSEIF (Z.eq. 82.0) THEN; D=-22451.900                           Pb (LEAD)
        ELSEIF (Z.eq. 83.0) THEN; D=-20054.400                           Bi (BISMUTH)
        ELSEIF (Z.eq. 84.0) THEN; D=-17146.000                           Po (POLONIUM)
        ELSEIF (Z.eq. 85.0) THEN; D=-13243.000                           At (ASTATINE)
        ELSEIF (Z.eq. 86.0) THEN; D=- 8631.000                           Rn (RADON)
        ELSEIF (Z.eq. 87.0) THEN; D=- 2840.000                           Fr (FRANCIUM)
        ELSEIF (Z.eq. 88.0) THEN; D=+ 3540.000                           Ra (RADIUM)
        ELSEIF (Z.eq. 89.0) THEN; D=+11130.000                           Ac (ACTINIUM)
                            ELSE; GOTO 1111
         endIF
                          ELSE; GOTO 1111
       endIF

         MassNucleus=(D+A*UAMU)/1.0d+06                                  [keV]/10^6=[GeV]

         RETURN
 1111    PRINT *, A, Z; STOP 'UNKNOWN NUCLEUS'
      END FUNCTION MassNucleus