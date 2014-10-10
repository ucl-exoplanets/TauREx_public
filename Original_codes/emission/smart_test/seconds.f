         SUBROUTINE SECONDS(T)
C --
C -- May 5 2006, Javier M-T
C Description: DATE_AND_TIME(DATE, TIME, ZONE, VALUES) gets the corresponding 
C date and time information from the real-time system clock. DATE is INTENT(OUT) 
C and has form ccyymmdd. TIME is INTENT(OUT) and has form hhmmss.sss. 
C ZONE is INTENT(OUT) and has form (+-)hhmm, representing the difference with 
C respect to Coordinated Universal Time (UTC). Unavailable time and date parameters 
C return blanks. 
C VALUES is INTENT(OUT) and provides the following: 
C
C --  VALUE(1):  The year 
C --  VALUE(2):  The month 
C --  VALUE(3):  The day of the month 
C --  VALUE(4):  Time difference with UTC in minutes 
C --  VALUE(5):  The hour of the day 
C --  VALUE(6):  The minutes of the hour 
C --  VALUE(7):  The seconds of the minute 
C --  VALUE(8):  The milliseconds of the second 
  
              character(8)  :: date
              character(10) :: time
              character(5)  :: zone
              integer,dimension(8) :: values
              integer secnds_
              ! using keyword arguments
              call date_and_time(date,time,zone,values)
              call date_and_time(DATE=date,ZONE=zone)
              call date_and_time(TIME=time)
              call date_and_time(VALUES=values)
              print '(a,2x,a,2x,a)', date, time, zone
              print '(8i5))', values(1),values(2)
              secnds_ = values(7) +values(6)*60+values(5)*3600- T
    
              return
              end
              
