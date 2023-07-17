SUBROUTINE ISEQU(NCHAI0,DT,MBSEC0,EX)
   ! c        REAL ISECM
   ! c        INTEGER BCHASE,BCHADI
   ! c*******************************************************************
   ! c*******************************************************************
   use global
   ! c        COMMON /INIT/   RHH(49),RQQ(49),C(49)
   ! c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
   ! c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
   ! c     &,        /WWW/    WORKM(49,3)
   ! c**********************************************************************
   DO 400 I=1,NCHAI0
      J=2*I
      J1=ICHASE(J-1)
      J2=ICHASE(J)
      IF (RQQ(J1) .GT. 0.) THEN
         ICHADI(I)=1
         IF (RQQ(J2) .LT. 0.) THEN
            PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
   c                 WRITE(*,*) CHAR(7),'DIRECTION WRONG IN CHANLE',I,'I'
   C                 STOP 'I1'
            DO 10 IJ=J1,J2
               IF (RQQ(IJ) .LT. 0.) RQQ(IJ)=1.E-5
10          CONTINUE
         END IF
      ELSE IF (RQQ(J1) .LT. 0.) THEN
         ICHADI(I)=-1
         IF (RQQ(J2) .GT. 0.) THEN
            PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
            DO 11 IJ=J1,J2
               IF (RQQ(IJ) .GT. 0.) RQQ(IJ)=-1.E-5
11          CONTINUE
         END IF
      ELSE
         PRINT *,'J1= ',J1,'  J2=',J2,'  I=',I
         PAUSE 'I3'
      END IF
      CALL ORINSM(J1,J2,ICHADI(I),DT,EX)
      J=J1-MBSEC0-I+1
      ISECM(J,1)=WORKM(1,1)
      ISECM(J,2)=WORKM(1,2)
      ISECM(J,3)=WORKM(1,3)
      DO 300 J=J1-MBSEC0-I+2,J2-MBSEC0-I
         JJ=J+MBSEC0+I-J1
         EPF=1./(1.-ISECM(J-1,3)*WORKM(JJ,2))	
         ISECM(J,1)=(WORKM(JJ,1)+WORKM(JJ,2)*ISECM(J-1,1))*EPF
         ISECM(J,2)=WORKM(JJ,2)*ISECM(J-1,2)*EPF
         ISECM(J,3)=WORKM(JJ,3)*EPF
300   CONTINUE
400   CONTINUE
END