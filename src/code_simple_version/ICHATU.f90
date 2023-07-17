INTEGER FUNCTION ICHATU(NSEC,NTUN0)
    ! c*******************************************************************
    ! c        COMMON /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
    ! c     &,                 NTUNTU(33)
    ! c**********************************************************************
    use global
    DO 500 J=1,NTUN0
        K=NTUNSS(J)
        KK=NTUNSS(J+1)
    301     K=K+1
        IF(NSEC .EQ. NTUNSE(K)) THEN
            ICHATU=J
        RETURN
        END IF
        IF (K .LT. KK) GOTO 301
    500   CONTINUE
END