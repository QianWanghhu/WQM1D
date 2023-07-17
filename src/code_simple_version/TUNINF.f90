SUBROUTINE TUNINF(NBC0,NCHAI0,NTUN0) ! Specify the channels and nodes for each section.
   use global
   ! c        INTEGER BCHASE,BCHADI
   ! c*******************************************************************
   ! c        COMMON /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
   ! c     &,        /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
   ! c     &,                 NTUNTU(33)
   ! c**********************************************************************
   DO 200 I=1,NTUN0
   DO 100 KS=NTUNSS(I)+1,NTUNSS(I+1)
      NS=NTUNSE(KS) ! Get the index of the KS-th section.
      NC=0
   10	   NC=NC+1
      IF (NC.LE.NBC0 .AND. BCHASE(2*NC).EQ.NS) THEN !LE <=; if NS equals the end section index of the NC-th outer channel
         NTUNCH(KS)=NC ! The KS-th section belongs to the NC-th outer channel.
         NTUNTU(KS)=0
         MSIGN(KS)=1
      ELSE IF (NC.LE.NCHAI0 .AND. ICHASE(2*NC).EQ.NS) THEN ! if NS equals the end section index of the NC-th inner channel
         NTUNCH(KS)=NC ! The KS-th section belongs to the NC-th inner channel.
         NSH=ICHASE(2*NC-1) ! Head section index of the NC-th inner channel.
         NTT=ICHATU(NSH,NTUN0) ! Node of the head section index NSH.
         NTUNTU(KS)=NTT 
         MSIGN(KS)=1
      ELSE IF (NC.LE.NCHAI0 .AND. ICHASE(2*NC-1).EQ.NS) THEN ! if NS equals the head section index of the NC-th inner channel
         NTUNCH(KS)=NC
         NSE=ICHASE(2*NC) ! End section index of the NC-th inner channel.
         NTT=ICHATU(NSE,NTUN0) ! Node of the End section index NSH.
         NTUNTU(KS)=NTT
         MSIGN(KS)=-1
   ! C                ELSE IF (NC .GT. NCHAI0) THEN
   ! C                   STOP 100
      ELSE
         GOTO 10
      END IF
   100	CONTINUE
   200	CONTINUE
END