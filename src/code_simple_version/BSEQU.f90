SUBROUTINE BSEQU(NBC0,DT,EX)
    ! c        REAL ISECM
    ! c        INTEGER BCHASE,BCHADI
    ! c*******************************************************************
    ! c        COMMON /BC/     RBC(13,4296),DBC(13)
    ! c     &,        /INIT/   RHH(49),RQQ(49),C(49)
    ! c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
    ! c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
    ! c     &,        /WWW/    WORKM(49,3)
    ! c**********************************************************************
    use global
    DO 300 I=1,NBC0
        J=2*I
        J1=BCHASE(J-1) !The initial section of the outer channel.
        J2=BCHASE(J)  !The end section of the outer channel.
        IF (RQQ(J1) .GT. 0.) THEN
            BCHADI(I)=1
            IF (RQQ(J2) .LE. 0.) THEN
                PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
c                 WRITE(*,*) CHAR(7),'DIRECTION WRONG IN CHANLE',I,'B'
C                 STOP 'B1'
                DO 10 IJ=J1,J2
                    IF (RQQ(IJ) .LE. 0.) RQQ(IJ)=1.E-5
10              CONTINUE
            END IF
            CALL ORINSM(J1,J2,1,DT,EX)
            BSECM(J1-I+1,1)=WORKM(1,1)+WORKM(1,2)*DBC(I)
            BSECM(J1-I+1,2)=WORKM(1,3)
            DO 100 J=J1-I+2,J2-I
                JJ=J+I-J1
                EPF=1./(1.-BSECM(J-1,2)*WORKM(JJ,2))
                BSECM(J,1)=(WORKM(JJ,1)+BSECM(J-1,1)*WORKM(JJ,2))*EPF
                BSECM(J,2)=WORKM(JJ,3)*EPF
100         CONTINUE
        ELSE IF (RQQ(J1) .LT. 0.) THEN
            BCHADI(I)=-1
            IF (RQQ(J2) .GE. 0.) THEN
                PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
c                 WRITE(*,*) CHAR(7),'DIRECTION WRONG IN CHANLE',I,'B'
C                 STOP 'B2'
                DO 11 IJ=J1,J2
                    IF (RQQ(IJ) .GE. 0.) RQQ(IJ)=-1.E-5
11              CONTINUE
            END IF
            CALL ORINSM(J1,J2,-1,DT,EX)
            JJ=J2-J1
            BSECM(J1-I+1,1)=WORKM(JJ,1)+WORKM(JJ,3)*DBC(I)
            BSECM(J1-I+1,2)=WORKM(1,2)
C             WRITE(*,'(2I5,5E12.4)')J1-I+1,JJ,WORKM(1,1),WORKM(1,2)
C    &                              ,WORKM(1,3),BSECM(3,1),BSECM(3,2)
            DO 200 J=J1-I+2,J2-I
                JJ=JJ-1
                EPF=1./(1.-BSECM(J-1,2)*WORKM(JJ,3))
                BSECM(J,1)=(WORKM(JJ,1)+BSECM(J-1,1)*WORKM(JJ,3))*EPF
                BSECM(J,2)=WORKM(JJ,2)*EPF
200         CONTINUE
        ELSE
            PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
c              WRITE(*,*) CHAR(7),'Q=0 WRONG IN CHANLE',I,'B'
            STOP 'B3'
        END IF
300     CONTINUE
	  END