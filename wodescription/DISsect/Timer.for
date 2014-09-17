************************************************************************
*                                                                      *
*                       C H R O N O M E T R Y                          *
*                                                                      *
*                          by V. A. Naumov                             *
*           Version of 01/06/2000 for DEC Visual Fortran 5.0           *
*                                                                      *
*     ------------------------------------------------------------     *
*     SUBROUTINE Chronometer                                           *
*          ENTRY TimeRef,TimeRef1,RunTime,RunTime1                     *
*     SUBROUTINE SystemUserTime                                        *
*          ENTRY SystemUserTime_D,SystemUserTime_E                     *
*     SUBROUTINE ShowDate                                              *
*     SUBROUTINE ShowTime                                              *
*     ------------------------------------------------------------     *
*                                                                      *
************************************************************************

************************************************************************
      SUBROUTINE Chronometer(Put,DateTime,Nout)                          Uses DFPORT
************************************************************************

        IMPLICIT INTEGER(2)(C,N)
            SAVE N,s0,s1,s

      LOGICAL(2) Put              ! to manage the output into a logfile
      LOGICAL(2) DateTime         ! to manage the current date and time

         REAL*4 s0/0.0/
         REAL s1
         REAL s,TIMEF

*     ================================================================ *
      ENTRY TimeRef(Put,DateTime,Nout)  ! zero-time reference (RunTime)
*     ================================================================ *
      s=TIMEF()                   ! setting up the zero-time reference
      N=Nout
         IF (DateTime) THEN       ! output the current date and time
           CALL GetDat(Cyear,Cmonth,Cday)
           CALL GetTim(Chour,Cminute,Csecond,Ccent)
           WRITE(*,1) Cday,Cmonth,Cyear,Chour,Cminute,Csecond,Ccent
           IF (Put) WRITE(N,1) Cday,Cmonth,Cyear,
     #                         Chour,Cminute,Csecond,Ccent
      endIF
      RETURN
*     ================================================================ *
      ENTRY TimeRef1(Put,DateTime,Nout) ! Zero-time reference (RunTime1)
*     ================================================================ *
      s0=SECNDS(s0)               ! setting up the zero-time reference
      N=Nout
         IF (DateTime) THEN       ! output the current date and time    
           CALL GetDat(Cyear,Cmonth,Cday)
           CALL GetTim(Chour,Cminute,Csecond,Ccent)
           WRITE(*,1) Cday,Cmonth,Cyear,Chour,Cminute,Csecond,Ccent
           IF (Put) WRITE(N,1) Cday,Cmonth,Cyear,
     #                         Chour,Cminute,Csecond,Ccent
      endIF
      RETURN
*     ================================================================ *
      ENTRY RunTime(Put,DateTime)                  ! Runtime evaluating
*     ================================================================ *
                   ! A PortLib function. Returns the number of seconds
      s=TIMEF()    ! since the first time it is called (precise a hun-
                   ! dredth of a second), or zero.

         Nday=INT(s/86400)         ! number of days
            s=s-Nday*86400
        Nhour=INT(s/3600)          ! number of hours
            s=s-Nhour*3600
      Nminute=INT(s/60)            ! number of minutes
            s=s-Nminute*60
      Nsecond=INT(s)               ! number of seconds
            s=s-Nsecond
        Ncent=s*100                ! number of centesimals of a second

         IF (DateTime) THEN        ! output current date, time, runtime
           CALL GetDat(Cyear,Cmonth,Cday)
           CALL GetTim(Chour,Cminute,Csecond,Ccent)
           WRITE(*,2) Cday,Cmonth,Cyear,Chour,Cminute,Csecond,Ccent,
     #                            Nday,Nhour,Nminute,Nsecond,Ncent
           IF (Put) WRITE(N,2) Cday,Cmonth,Cyear,
     #                         Chour,Cminute,Csecond,Ccent,
     #                         Nday,Nhour,Nminute,Nsecond,Ncent
                       ELSE        ! output runtime only
           WRITE(*,3) Nday,Nhour,Nminute,Nsecond,Ncent
           IF (Put) WRITE(N,3) Nday,Nhour,Nminute,Nsecond,Ncent
      endIF
      RETURN
*     ================================================================ *
      ENTRY RunTime1(Put,DateTime)                 ! Runtime evaluating
*     ================================================================ *
                    ! A PortLib function. Returns the number of seconds
      s1=SECNDS(s0) ! that have elapsed since midnight, minus s0  (pre-
                    ! cise a hundredth of a second).

         Nday=INT(s1/86400)        ! number of days
           s1=s1-Nday*86400
        Nhour=INT(s1/3600)         ! number of hours
           s1=s1-Nhour*3600
      Nminute=INT(s1/60)           ! number of minutes
           s1=s1-Nminute*60
      Nsecond=INT(s1)              ! number of seconds
           s1=s1-Nsecond
        Ncent=s1*100               ! number of centesimals of a second

         IF (DateTime) THEN        ! output current date, time, runtime
           CALL GetDat(Cyear,Cmonth,Cday)
           CALL GetTim(Chour,Cminute,Csecond,Ccent)
           WRITE(*,2) Cday,Cmonth,Cyear,Chour,Cminute,Csecond,Ccent,
     #                            Nday,Nhour,Nminute,Nsecond,Ncent
           IF (Put) WRITE(N,2) Cday,Cmonth,Cyear,Chour,Cminute,Csecond,
     #                         Ccent,Nday,Nhour,Nminute,Nsecond,Ncent
                       ELSE        ! output runtime only
           WRITE(*,3) Nday,Nhour,Nminute,Nsecond,Ncent
           IF (Put) WRITE(N,3) Nday,Nhour,Nminute,Nsecond,Ncent
      endIF

         RETURN

    1 FORMAT( '  Date: ',I3,'/',I2.2,'/',I4.4/
     #        '  Time: ',I2,':',I2.2,':',I2.2,':',I2.2/)
    2 FORMAT(/'  Date: ',I3,'/',I2.2,'/',I4.4,
     #        '  Time: ',I2,':',I2.2,':',I2.2,':',I2.2,
     #        '  (RunTime: ',I2.2,'d ',I2.2,'h ',I2.2,'m ',
     #                                 I2.2,'.',I2.2,'s)'/)
    3 FORMAT(60x,'(',I2.2,'d ',I2.2,'h ',I2.2,'m ',I2.2,'.',I2.2,'s)')

      END SUBROUTINE Chronometer

************************************************************************
      SUBROUTINE SystemUserTime(Nout)
************************************************************************

      REAL(4) S,TA(2)

*     ================================================================ *
      ENTRY SystemUserTime_D(Nout)                                       Uses DFPORT
*     ================================================================ *
         S=DTIME(TA)                                                     Portability function
         IF (Nout>0) WRITE(Nout,1) S,TA(1),TA(2)
                     WRITE(*   ,1) S,TA(1),TA(2)

         RETURN

*     ================================================================ *
      ENTRY SystemUserTime_E(Nout)                                       Uses DFPORT
*     ================================================================ *
         S=ETIME(TA)                                                     Portability function
         IF (Nout>0) WRITE(Nout,2) S,TA(1),TA(2)
                     WRITE(*   ,2) S,TA(1),TA(2)

         RETURN

    1 FORMAT(1X,60('_')/62X,/,1X,'Program has been running for',1pE9.2,
     ,                ' seconds. This includes'/1X,E10.3,' sec of user',
     ,      ' time and ',E9.2,' sec of system time'/1X,60('_')/)
    2 FORMAT(1X,60('_')/62X,/,1X,'Program has used',1pE9.2,
     ,    ' seconds of CPU time. This includes'/1X,E10.3,' sec of user',
     ,    ' time and ',E9.2,' sec of system time'/1X,60('_')/)
      END SUBROUTINE SystemUserTime

************************************************************************
      SUBROUTINE ShowDate(year,month,day)                                Uses DFLIB
************************************************************************

         INTEGER(2)  year,month,day

         CALL GetDat(year,month,day)                                     Runtime routine
         WRITE(*,'("  DATE: ",I3,"/",I2.2,"/",I4.4/)') day,month,year

         RETURN
      END SUBROUTINE ShowDate

************************************************************************
      SUBROUTINE ShowTime(hour,minute,second,cent)                       Uses DFLIB
************************************************************************

         INTEGER(2)  hour,minute,second,cent

         CALL GetTim(hour,minute,second,cent)                            Runtime routine
         WRITE(*,'("  TIME: ",I2,":",I2.2,":",I2.2,":",I2.2/)')
     #            hour,minute,second,cent

         RETURN
      END SUBROUTINE ShowTime