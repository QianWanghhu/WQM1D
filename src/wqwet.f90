SUBROUTINE WQWET

  ! *** RWQATM
  ! ** COMPUTES WET ATMOSPHERIC DEPOSITION USING CONSTANT CONCENTRATIONS
  ! ** FOR THE 22 STATE VARIABLES MULTIPLIED BY THE RAINFALL FLOW RATE      !VB CHANGED 21 TO 22
  ! ** ENTERING EACH GRID CELL.  COMPUTED LOADS ARE IN G/DAY.

  ! CHANGE RECORD

  USE GLOBAL
  IMPLICIT NONE

  INTEGER :: L,NW
  REAL    :: CV2,TIME

  !  CV2 = CONVERSION TO GET UNITS OF G/DAY
  !  WQATM(NW) HAS UNITS OF MG/L
  !  RAINT(L) HAS UNITS OF M/SEC
  !  DXYP(L) HAS UNITS OF M2
  !  WQATML(L,KC,NW) HAS UNITS OF G/DAY

  CV2=86400.0
  DO NW=1,NWQV
    DO L=2,LA
      WQATML(L,KC,NW)=WQATM(NW)*RAINT(L)*DXYP(L)*CV2
    ENDDO
  ENDDO
  IF( ITNWQ == 0 .AND. DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'WQATM.DIA',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'WQATM.DIA',STATUS='UNKNOWN')
    IF( ISDYNSTP == 0 )THEN
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
    ELSE
      TIME=TIMESEC/86400.
    ENDIF
    WRITE(1,112) N,TIME
    DO L=2,LA
      WRITE(1,110) IL(L),JL(L),(WQATML(L,KC,NW),NW=1,NWQV)
    ENDDO
    CLOSE(1)
  ENDIF
    110 FORMAT(1X,2I4,2X,1P,7E11.3,/,15X,7E11.3,/,15X,7E11.3)
    112 FORMAT('# WET ATMOSPHERIC DEPOSITION DIAGNOSTIC FILE',/, &
      ' N, TIME = ', I10, F12.5/)
  RETURN
END

