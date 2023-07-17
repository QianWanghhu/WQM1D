C================================
C         hdsz.FOR              $
C================================

module global

    PARAMETER (MKK=9)		   !No. of sources 
    PARAMETER (NBC=8		   !No. of boundaries.	
    &,            NSTUN=4	   
    &,            NCHAI=32	   ! No. of inner channel.
    &,            MBSEC=16)	   ! End section index of the boundary.
C
    PARAMETER (MTSEC=80	    ! The total number of sections.
    &,            NTUN=23		! No. of nodes.
    &,            NDT=24*10	
    &,            NVC=80			! End section index of the inner channels.
    &,            NB22=2*NBC)
C
    PARAMETER (NI22=2*NCHAI
    &,            NCHATO=NBC+NCHAI
    &,            NB21=MBSEC-NBC
    &,            NI21=MTSEC-MBSEC-NCHAI)
C
    PARAMETER (NDY=2*NCHATO-NBC
    &,            N11=NTUN+1
    &, 	         NTUS=2*NCHATO-NBC
    &,            NM=NI22+NTUN)
C
    CHARACTER*20 BNAME(NBC)
    CHARACTER*30 CCCC
C
    REAL RHH(MTSEC)
    &,      RQQ(MTSEC)
    &,      C(MTSEC)
    &,      RQ1(MTSEC)
    &,      RBC(NBC,NDT)
    &,      DBC(NBC)
    &,      SECFOR(MTSEC,3)
    &,      DX(MTSEC)
    &,      KI(MTSEC) 
    &,      PAS(NSTUN),PAH(NSTUN)
C
    REAL BSECM(NB21,2)
    &,      ISECM(NI21,3)
    &,      TUNMAT(NM)
    &,      RTUNM(NTUN)
    &,      VERIC(NVC)
    &,      SC(MTSEC)
    &,      SQ(MTSEC)
C
    INTEGER BCHASE(NB22)
    &,         BCHADI(NBC)
    &,         ICHADI(NCHAI)
    &,         ICHASE(NI22)
    &,         VCSEC(NVC)
    &,         NNNS(NSTUN)
C
    INTEGER NTUNSS(N11)
    &,         NTUNSE(NTUS)
    &,         MSIGN(NTUS)
    &, 	      NTUNCH(NTUS)
    &, 	      NTUNTU(NTUS)
    &,         ITUNM1(NM)
    &,         ITUNM2(NTUN)
C
    INTEGER SEC,MIN,HOU,DAY,MON,YEA,HOUR,DAYO
C
    REAL WORKM(MTSEC,6)

    COMMON 
c	  /BC/
    $  RBC,DBC
c       /INIT/ 
    &, RHH,RQQ,C
c       /CHAP/
    &, SECFOR,DX
c       /CHA/ 
    &, BCHASE,ICHASE,BCHADI,ICHADI
c       /SSSSSS/
    &, BSECM,ISECM
c       /TUN/
    &, NTUNSS,NTUNSE,MSIGN,NTUNCH
    &, NTUNTU
c       /TUNM/  
    &, TUNMAT,RTUNM,ITUNM1,ITUNM2
c       /SOUR/ 
    &, SC,SQ
c       /PAH/ 
    &, PAH,PAS,NNNS
c       /WWW/ 
    &, WORKM
    c       /KI/
    &, KI


c*******************************************************************
c        COMMON /BC/     RBC(13,4296),DBC(13)
c     &,        /INIT/   RHH(49),RQQ(49),C(49)
c     &,        /CHAP/   SECFOR(49,3),DX(49)
c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
c     &,        /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
c     &,                 NTUNTU(33)
c     &,        /TUNM/   TUNMAT(31),RTUNM(11),ITUNM1(31),ITUNM2(11)
c     &,        /SOUR/   SC(49),SQ(49)
c     &,        /PAH/    PAH(4),PAS(4),NNNS(4)
c     &,        /WWW/    WORKM(49,3)
c**********************************************************************
    INTEGER SCQ(MKK)
    REAL RBCBOU(NBC)
end module
c**********************************************************************
C
C        OPEN(6,FILE='TU.DAT')
c        PRINT *,'MKK=',MKK,'NBC=',NBC, 'NSTUN=',NSTUN
c        PRINT *,'NCHAI=',NCHAI,'MBSEC=',MBSEC,'MTSEC=',MTSEC
c        PRINT *,'NTUN=',NTUN,'NDT=',NDT,'NVC=',NVC,'NB22=2*NBC',NB22
c	  PRINT *,'NI22=2*NCHAI',NI22
c        PRINT *,'NCHATO=NBC+NCHAI',NCHATO
c        PRINT *,'NB21=MBSEC-NBC',NB21
c        PRINT *,'NI21=MTSEC-MBSEC-NCHAI',NI21
c	  PRINT *,'NDY=2*NCHATO-NBC',NDY
c        PRINT *,'N11=NTUN+1',N11
c        PRINT *,'NTUS=2*NCHATO-NBC',NTUS
c        PRINT *,'NM=NI22+NTUN',NM
C        GOTO 5432
PROGRAM wq1d
    use global
    DATA KT /0/ KM /1/

    DO 315 I=1,NVC
315   VCSEC(I)=I
C
    OPEN (10,FILE='DATA\DX.DAT')
	READ (10,*)
	READ (10,*) (DX(I),I=1,MTSEC)
	CLOSE(10)

C Read section shape elements
    OPEN (10,FILE='DATA\SECFOR.DAT')
    READ (10,*) ((SECFOR(I,J),J=1,3),I=1,MTSEC) 
	CLOSE(10)
C Read inital concentrations of each sections.
    OPEN (10,FILE='DATA\INITIAL.DAT')
	READ (10,*) SEC,MIN,HOU,DAY,MON,YEA
    READ (10,*) (C(I),I=1,MTSEC) ! �������ʼ���?��Ũ��
	CLOSE(10)
!      WRITE(6,555)DAY,HOU,(C(I),I=1,MTSEC) 
! Initialize SC and SQ being zeros.
    CALL ZERO(SC,MTSEC) 
    CALL ZERO(SQ,MTSEC)
! Read information of point source: SCQ: section index, SQ: discharge, SC: concentrations
    OPEN (10,FILE='DATA\SCQcod.DAT') 
    READ(10,*)
    READ(10,*) (SCQ(I),I=1,MKK) 
    READ(10,*)
    READ (10,*) (SQ(SCQ(I)),I=1,MKK) 
    READ(10,*)
    READ (10,*) (SC(SCQ(I)),I=1,MKK) 
    CLOSE(10)

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c        OPEN (10,FILE='DATA\BOUND.NAM')
c	READ (10,*) RBCDT
c	READ (10,*) (RBCBOU(I),I=1,NBC)
c       CLOSE(10)							  ԭ��
c
c         DO 2 I=1,NBC
c         DO 2 J=1,ndt
c         RBC(I,J)=RBCBOU(I)
c        CONTINUE	  
c       write(*,*)( rbc(2,mm),mm=1,24)    	,1
c 	 write(*,*) ( rbc(kk,1),kk=1,12)	    ,2
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	 
!TODO: TO sort out reading processes of .NAM files
	OPEN (10,FILE='DATA\BOUND.NAM')
    READ (10,*) RBCDT
    READ (10,*) (RBCBOU(I),I=1,NBC)
    CLOSE(10)

    DO 2 I=1,NBC
        DO 2 J=1,24
            RBC(I,J)=RBCBOU(I)
2       CONTINUE	  	  
        write(*,*) 1				      	     	     	    
        write(*,*) ( rbc(kk,1),kk=1,12)	    ,3			     
    
    OPEN (10,FILE='DATA\BOUND1.NAM')
    READ (10,*) RBCDT
    READ (10,*) (RBCBOU(I),I=1,NBC)
    CLOSE(10)

    DO 3 I=1,NBC
        DO 3 J=25,48
            RBC(I,J)=RBCBOU(I)
3       CONTINUE	  	  
        write(*,*) 2
    write(*,*) ( rbc(kk,25),kk=1,12)	    ,3

	OPEN (10,FILE='DATA\BOUND2.NAM')
    READ (10,*) RBCDT
    READ (10,*) (RBCBOU(I),I=1,NBC)
    CLOSE(10)

    DO 4 I=1,NBC
        DO 4 J=49,72
            RBC(I,J)=RBCBOU(I)
4       CONTINUE	
        write(*,*)   3

    write(*,*) ( rbc(kk,49),kk=1,12)	    ,3

    OPEN (10,FILE='DATA\BOUND3.NAM')
    READ (10,*) RBCDT
    READ (10,*) (RBCBOU(I),I=1,NBC)
    CLOSE(10)

    DO 5 I=1,NBC
        DO 5 J=73,96
            RBC(I,J)=RBCBOU(I)
5       CONTINUE	
        write(*,*) 4
	write(*,*) ( rbc(kk,73),kk=1,12)	    ,3

	OPEN (10,FILE='DATA\BOUND4.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 6 I=1,NBC
         DO 6 J=97,120
         RBC(I,J)=RBCBOU(I)
6        CONTINUE	
	               write(*,*) 5	   
	write(*,*) ( rbc(kk,97),kk=1,12)	    ,3


	OPEN (10,FILE='DATA\BOUND5.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 7 I=1,NBC
         DO 7 J=121,144
         RBC(I,J)=RBCBOU(I)
7        CONTINUE	
	              write(*,*) 6
	write(*,*) ( rbc(kk,121),kk=1,12)	    ,3


	OPEN (10,FILE='DATA\BOUND6.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 8 I=1,NBC
         DO 8 J=145,168
         RBC(I,J)=RBCBOU(I)
8        CONTINUE	
	                  write(*,*) 7
	write(*,*) ( rbc(kk,145),kk=1,12)	    ,3



	OPEN (10,FILE='DATA\BOUND7.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 9 I=1,NBC
         DO 9 J=169,192
         RBC(I,J)=RBCBOU(I)
9        CONTINUE	
	             write(*,*)   8

	write(*,*) ( rbc(kk,169),kk=1,12)	    ,3



	OPEN (10,FILE='DATA\BOUND8.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 11 I=1,NBC
         DO 11 J=193,216
         RBC(I,J)=RBCBOU(I)
11        CONTINUE	
	           write(*,*) 	9
	write(*,*) ( rbc(kk,193),kk=1,12)	    ,3



	OPEN (10,FILE='DATA\BOUND9.NAM')
	  READ (10,*) RBCDT
	  READ (10,*) (RBCBOU(I),I=1,NBC)
        CLOSE(10)

	   DO 12 I=1,NBC
         DO 12 J=217,240
         RBC(I,J)=RBCBOU(I)
12        CONTINUE
	           write(*,*)  10
	write(*,*) ( rbc(kk,217),kk=1,12)	    ,3 
c	write(*,*)  rbc

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Read the NETCODE.DAT file
    OPEN (10,FILE='DATA\NETCODE.DAT')
    READ (10,*) 
    READ (10,*) (NNNS(I),NNN,I=1,NSTUN)
    READ (10,*)
    READ (10,*) (BCHASE(I),I=2,NBC+1)
	READ (10,*)
	READ (10,*) (ICHASE(I),I=2,NCHAI+1)
	READ (10,*)
	READ (10,*) (NTUNSS(I),I=2,NTUN+1)
	READ (10,*)
	READ (10,*) (NTUNSE(I),I=1,NTUS)
	CLOSE(10)
    ! Input data from hydrodynamic model�� 1�� water depth, 2: discharge
    OPEN(1,FILE='DATA\H.DAT')
    !	write(*,*) 1  
	OPEN(2,FILE='DATA\Q.DAT')
    !	write(*,*) 1

	OPEN (20,FILE='k2O.dat')
    OPEN (21,FILE='k21.dat') 
    OPEN (22,FILE='22.dat')
c      OPEN (22,FILE='D:\nsbdouhe\001\szmx\data\bound.nam')


!   READ(*,*) DT,EX
    DT=3600.0
    EX=2.5
c      PRINT *, DT,EX
C	NEND=3600/INT(DT)*NEND
	HOUR=HOU
	DAYO=DAY
    MONO=MON
    KCOUNT=0
    KCO1=0
C
    NTUNSS(1)=0
	DO 10 I=2,NTUN+1 ! Calcualte the cumulative number of sections belonging to the 1 to the ith nodes.
10	   NTUNSS(I)=NTUNSS(I)+NTUNSS(I-1) 
    CALL CHAINF(NBC,NCHAI) ! CHAINF: Function for indexing channel sections.
    CALL TUNINF(NBC,NCHAI,NTUN) !   TUNINF: Function for specifying the nodes for each section.
C
C      WRITE(20,'(/I4,1H.,I2)') YEA,MON
C      WRITE(21,'(/I4,1H.,I2,1H.,I2)') YEA,MON,DAY
c      WRITE(*,'(I5,1H.,I2,1H.,I2,12H ----- Begin)') YEA,MON,DAY
    OPEN (10,FILE='DATA\KKIc.DAT') ! Read degradation coefficient for all sections
    READ (10,*)
    READ (10,*) (KI(I),I=1,MTSEC)
    CLOSE(10)

    KPP=-2
    KG=0
    SHUIZ=6.0
c	  write(*,*)  "��һ�� ˮ��ģ��------���� �� ����	  "
    WRITE(21,'(i5,20i6)') 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19
    WRITE(22,'(i5,20i6)') 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19
c	    WRITE(22,'(i5,13i6)') 0,1,2,3,4,5,6,7,8,9,10,11,12,13

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC
CCC    CALCULATION BEGIN
CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

10000   KT=KT+1
    if (mod(kt,100).eq.0) then 
        WRITE(*,*) Kt
    endif
	    
    ! Read water depth and discharge for each int. hour.
    IF(KT-KT/24*24.EQ.1) THEN  
        READ (1,*) CCCC
	    READ (1,'(10F14.3)') (RHH(I),I=1,MTSEC)
	    READ (2,*) CCCC
        READ (2,'(10F14.2)') (RQQ(I),I=1,MTSEC)
	    DO 789 I=1,MTSEC
            IF(RQQ(I).EQ.0) RQQ(I)=1E-5
789     CONTINUE
    ENDIF

    CALL SERBOU(REAL(KT)*DT,RBCDT,NDT,NBC)  !  SERBOU: Linear interpolation for the initial conentrations.
    CALL BSEQU(NBC,DT,EX) ! BSEQU: Calculate the coefficient matrix for outer channels.
    CALL ISEQU(NCHAI,DT,MBSEC,EX) !  BSEQU: Calculate the coefficient matrix for inner channels.
    CALL TUNEQU(NTUN,MBSEC,DT) ! Calculate the coefficient matrix for nodes.
    CALL SOLVER(TUNMAT,ITUNM1,ITUNM2,RTUNM,NM,NTUN,1.E-12,IE) ! Solving the coefficient matrix for nodes.

    CX=0.001
    DO 1200 I=1,NTUN
        if(rtunm(i).ge.100.0.or.rtunm(i).le.0.) then
            TQ=0.0
            TQC=0.0
            NP=NTUNSS(I)+1
            NPP=NTUNSS(I+1)
            DO 1100 J=NP,NPP
                NTT=NTUNTU(J)
                NTC=NTUNCH(J)
                IGN=MSIGN(J)
                IF (NTT .EQ. 0) THEN
                    KGN=BCHADI(NTC)
                    NSH=BCHASE(2*NTC-1)
                    NSE=BCHASE(2*NTC)
                    IF(KGN.GT.0) THEN
                        TQ=TQ+ABS(RQQ(NSE))
                        TQC=TQC+ABS(RQQ(NSE)*C(NSE))
                    ENDIF
                ELSE
                    KGN=ICHADI(NTC)
                    NSH=ICHASE(2*NTC-1)
                    NSE=ICHASE(2*NTC)
                    IF(KGN.GT.0.AND.IGN.GT.0) THEN
                        TQ=TQ+ABS(RQQ(NSE))
                        TQC=TQC+ABS(RQQ(NSE)*C(NSE))
                    ENDIF
                    IF(KGN.LT.0.AND.IGN.LT.0) THEN
                        TQ=TQ+ABS(RQQ(NSH))
                        TQC=TQC+ABS(RQQ(NSH)*abs(C(NSH)))
                    ENDIF
                END IF

1100            CONTINUE
            if(tq.le.0.05) then
                rtunm(i)=10.0
                goto 1200
            endif
            rtunm(i)=TQC/TQ
        endif
1200    CONTINUE 

    CALL RBTUNL(NTUN,MBSEC)
    CALL RBBSEC(NBC)
    CALL RBISEC(NCHAI,MBSEC)
    CALL CLOCK(DT,SEC,MIN,HOU,DAY,MON,YEA)

    DO 99 I=1,MTSEC
        IF (C(I) .LT. 0.) C(I)=0.
        IF (C(I) .GT. 100.) C(I)=100.0
99  CONTINUE

    DO 101 I=1,MTSEC
        veric(i)=veric(I)+c(I)*dt/86400.
101 CONTINUE
    SHUIZ=(SHUIZ+c(11))/2. 
                 
    IF (DAY.NE.DAYO) THEN
        KG=KG+1
        if (mod(kg,50).eq.0) then 
            WRITE(*,*) KG
        endif
        WRITE(20,*) kg
        WRITE(20,'(5F10.3)') veric
        WRITE(21,'(i5,19F6.1)') kg,veric(80),
    *							     veric(31),
	*							     veric(39),
    *		                         veric(67),
    *                                 veric(27),
    *                                 veric(33),
    *                                 veric(30),
    *                                 veric(43),
    *                                 veric(53),
    *                                 veric(47),
    *                                 veric(24),
    *                                 veric(17),
    *                                 veric(25),
    *                                 veric(51),
    *                                 veric(22),
    *	                             (veric(63)+veric(64))/2,
    *	                             veric(79),
    *                                 veric(91),
    *                                 veric(112)
                                     

        WRITE(22,'(i5,19F6.1)') kg, (veric(79)+veric(80))/2,	                                    
    *                                (veric(31)+veric(32))/2,
    *                                (veric(39)+veric(40))/2,
    *                                (veric(67)+veric(68))/2,
    *                                (veric(27)+veric(28))/2,   	                                     
    *                                (veric(33)+veric(34))/2,
    *                                (veric(29)+veric(30))/2,
    *                                (veric(43)+veric(44))/2,
    *                                (veric(53)+veric(54))/2,   
    *                                (veric(47)+veric(48))/2,
    *	 							(veric(23)+veric(24))/2,
    *	 							(veric(17)+veric(18))/2,
    *	                  		    (veric(25)+veric(26))/2,
    *	                            (veric(51)+veric(52))/2,
    *	                            (veric(21)+veric(22))/2,
    *	                            veric(63),
    *	                            (veric(79)+veric(80))/2,
    *	                            (veric(91)+veric(92))/2,
    *	                            (veric(111)+veric(112))/2


		 
		 
        DO 108 I=1,MTSEC
            veric(I)=0.
108     CONTINUE
        DAYO=DAY
    ENDIF

    IF(KT.LT.NDT) GOTO 10000
    close(20)
    close(21)
    WRITE(22,'(i5)') 3600
    WRITE(22,'(5F10.3)') 6.0,shuiz,6.0,6.0,6.0
    close(22)

    STOP 'successful end'
END
