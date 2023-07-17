SUBROUTINE ORINSM(NSH,NSE,KDI,DT,EX)
    ! c        REAL KI
    ! c*******************************************************************
    ! c*******************************************************************
    ! c        COMMON /INIT/   RHH(49),RQQ(49),C(49)
    ! c     &,        /CHAP/   SECFOR(49,3),DX(49)
    ! c     &,        /SOUR/   SC(49),SQ(49)
    ! c     &,        /WWW/    WORKM(49,3)
    ! c**********************************************************************
    ! c        COMMON /KI/KI(49)
    use global
    IF (KDI .GT. 0) THEN ! Return the water section area and width of the section after the head section.
        X1=DX(NSH+1)
        ! SERSEC: Calculate the water section area and width of the section after the head section.
        CALL SERSEC(NSH,A1,B1)
    ELSE
        X1=DX(NSE)
        CALL SERSEC(NSE,A1,B1)
    END IF
    DO 100 I=1,NSE-NSH !Loop over sections of a DX
        IF (KDI .GT. 0) THEN
            J=NSH+I
            JM1=J-1
            JS1=J+1
            CALL SERSEC(J,A2,B2)
            IF (JS1 .GT. NSE) THEN
                X2=3.*B2
C                 X2=X1
            ELSE
                X2=DX(JS1)
            END IF
        ELSE
            J=NSE-I
            JM1=J+1
            JS1=J-1
            CALL SERSEC(J,A2,B2)
            IF (JS1 .LT. NSH) THEN
                 X2=3.*B2
C                 X2=X1
            ELSE
                X2=DX(J)
            END IF
        END IF
        Q1=ABS(RQQ(JM1))
        Q2=ABS(RQQ(J))
        U1=Q1/A1
        U2=Q2/A2
        U2=.5*(U1+U2)
        A=DT*(U2+2.*EX/(X1+X2))/X1
        IF (KI(NSE).GT.03) THEN
		    PRINT *, KI(NSE)
            PAUSE
	    ENDIF		
        B=1.+DT*(U2/X1+KI(NSE)/86400.+2.*EX/X2/X1)
        G=2.*DT*EX/X2/(X1+X2)
        P=2.*DT/X1/(A1+A2)*SQ(J)*SC(J)+C(J) ! Concentrations including point source inputs.
C          A=DT*(U1+EX/X2)/X1
C          B=1.+DT*(U2/X1+KI(J)/86400.+2.*EX/X2/X1)
C          G=DT*EX/X2/X1
C          P=2.*DT/X1/(A1+A2)*SQ(J)*SC(J)+C(J)
        WORKM(I,1)=P/B
        WORKM(I,2)=A/B
        WORKM(I,3)=G/B
        X1=X2
        A1=A2
        B1=B2
100     CONTINUE
END