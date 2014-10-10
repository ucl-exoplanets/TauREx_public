
!	VAX FORTRAN secnds function
!
!	usage:  REAL secnds, secMidnight
!		secMidnight = secnds(.0)
!
!	Sets secMidnight to the number of seconds since midnight.
!	The argument is subtracted from the seconds.

	INCLUDE "time.inc"

	REAL FUNCTION secnds_(t)
	IMPLICIT NONE
	REAL t
	INTEGER localtime, now
	RECORD /tm/clock
	POINTER (p_clock, clock)

	CALL time(now)
	p_clock=localtime(now)

	secnds_ = clock.tm_sec +
     +           clock.tm_min * 60 +
     +		 clock.tm_hour * 3600 -t
	END

$IF DEFINED(UCASE_NAMES)
	REAL FUNCTION SECNDS(t)
	IMPLICIT NONE
	REAL secnds_,t
	SECNDS = secnds_(t)
	END
$ENDIF

$IF DEFINED(UCASE_USCORE_NAMES)
	REAL FUNCTION SECNDS_(t)
	IMPLICIT NONE
	REAL secnds_,t
	SECNDS_ = secnds_(t)
	END
$ENDIF
