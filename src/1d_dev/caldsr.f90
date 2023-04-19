
SUBROUTINE CALDSR(DAYNEXT)
  ! Calculate the daily average solar radiation from the aser data
  ! Define variables
  ! TODO: DECIDE what variables are to return 

  IMPLICIT NONE

  REAL(RKD), STATIC :: DAYNEXT
  REAL(RKD), STATIC :: DAYNEXT
  REAL(RKD), STATIC :: SUNDAY1, SUNDAY2
  REAL  , STATIC :: SUNSOL1, SUNSOL2
  REAL  , STATIC :: SUNFRC1, SUNFRC2, WQTSDTINC
  INTEGER, STATIC :: M

   IF( ITNWQ == 0 .AND. IWQSUN > 1 .AND. NASER > 0 )THEN
        ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
        SUNDAY1 = DAYNEXT-1.
        SUNDAY2 = DAYNEXT

        ! *** FIND 1ST POINT
        ! TSATM is a global variable
        M = 1
        DO WHILE (TSATM(1).TIM(M) < SUNDAY1)
          M = M+1
        END DO
        
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL1 = 0.0
        DO WHILE (TSATM(1).TIM(M) < SUNDAY1)
          M1 = M1+1
          IF( TSATM(1).VAL(M,6) > 0. )THEN
            M2 = M2+1
            SUNSOL1=SUNSOL1+TSATM(1).VAL(M,6)
          ENDIF
          M = M+1
        END DO
        IF( M1 > 0 )THEN
          SUNFRC1=FLOAT(M2)/FLOAT(M1)
          SUNSOL1=SUNSOL1/FLOAT(M1)
        ELSE
         SUNFRC1=1.0
        ENDIF
        
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL2 = 0.
        DO WHILE (TSATM(1).TIM(M) < SUNDAY2)
          M1 = M1+1
          IF( TSATM(1).VAL(M,6) > 0. )THEN
            M2 = M2+1
            SUNSOL2=SUNSOL2+TSATM(1).VAL(M,6)
          ENDIF
          M = M+1
        END DO
        IF( M1 > 0 )THEN
          SUNFRC2=FLOAT(M2)/FLOAT(M1)
          SUNSOL2=SUNSOL2/FLOAT(M1)
        ELSE
          SUNFRC2=1.
        ENDIF
    ENDIF

    RETURN

END