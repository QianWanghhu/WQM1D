
SUBROUTINE RBTUNL(NTUN0,MBSEC0)
    ! c        REAL ISECM
    ! c        INTEGER BCHASE,BCHADI
    ! c*******************************************************************
    ! c*******************************************************************
    ! c        COMMON /INIT/   RHH(49),RQQ(49),C(49)
    ! c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
    ! c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
    ! c     &,        /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
    ! c     &,                 NTUNTU(33)
    ! c     &,        /TUNM/   TUNMAT(31),RTUNM(11),ITUNM1(31),ITUNM2(11)
    ! c**********************************************************************
    use global

    DO 200 I=1,NTUN0
        NP=NTUNSS(I)+1
        NPP=NTUNSS(I+1)
        DO 100 J=NP,NPP
        NTT=NTUNTU(J)
        NTC=NTUNCH(J)
        IGN=MSIGN(J)
        IF (NTT .EQ. 0) THEN
            KGN=BCHADI(NTC)
            NSH=BCHASE(2*NTC-1)
            NSE=BCHASE(2*NTC)
            IF (KGN .GT. 0) THEN
                NSM=NSE-NTC
                C(NSE)=BSECM(NSM,1)+BSECM(NSM,2)*RTUNM(I)
            ELSE
                C(NSE)=RTUNM(I)
            END IF
        ELSE
            KGN=ICHADI(NTC)
            NSH=ICHASE(2*NTC-1)
            NSE=ICHASE(2*NTC)
            IF (IGN .GT. 0) THEN
                IF (KGN .GT. 0) THEN
                    NSM=NSE-NTC-MBSEC0
                    C(NSE)=ISECM(NSM,1)+ISECM(NSM,2)*RTUNM(NTT)
     &                           +ISECM(NSM,3)*RTUNM(I)
                ELSE
                    C(NSE)=RTUNM(I)
                END IF
            ELSE
                IF (KGN .LT. 0) THEN
                    NSM=NSE-NTC-MBSEC0
                    C(NSH)=ISECM(NSM,1)+ISECM(NSM,2)*RTUNM(NTT)
     &                           +ISECM(NSM,3)*RTUNM(I)
                ELSE
                    C(NSH)=RTUNM(I)
                END IF
            END IF
        END IF
100	  CONTINUE
200   CONTINUE
END