
SUBROUTINE TUNEQU(NTUN0,MBSEC0,DT)
    ! c 	  REAL ISECM
    ! c        INTEGER BCHASE,BCHADI
    ! c*******************************************************************
    ! c*******************************************************************
    ! c        COMMON /BC/     RBC(13,4296),DBC(13)
    ! c     &,        /INIT/   RHH(49),RQQ(49),C(49)
    ! c     &,        /CHAP/   SECFOR(49,3),DX(49)
    ! c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
    ! c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
    ! c     &,        /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
    ! c     &,                 NTUNTU(33)
    ! c     &,        /TUNM/   TUNMAT(31),RTUNM(11),ITUNM1(31),ITUNM2(11)
    ! c     &,        /SOUR/   SC(49),SQ(49)
    ! c     &,        /PAH/    PAH(4),PAS(4),NNNS(4)
    ! c**********************************************************************
    use global
    KS=1
    NM0=1
    I=0
100 I=I+1
    IF (I .GT. NTUN0) RETURN
        COM=0.
        NCOM=0
        R=0.
        NP=NTUNSS(I)+1
        NPP=NTUNSS(I+1)
        ITUNM2(I)=NPP-NP+2
        IF(I.EQ.NNNS(KS)) THEN
            NS=NTUNSE(NP)
            DH=RHH(NS)-PAH(KS)
            COM=COM-DH*PAS(KS)/DT
            PAH(KS)=RHH(NS)
            KS=KS+1
        END IF
        DO 200 J=NP,NPP
            NTT=NTUNTU(J)
            NTC=NTUNCH(J) ! Channel index
            IF (NTT .EQ. 0) THEN  ! Calculate mass balance at the end sections of outer channels.
	            ITUNM2(I)=ITUNM2(I)-1
                NSE=BCHASE(2*NTC)
                NSM=NSE-NTC ! River reach index
                IF (BCHADI(NTC) .LT. 0) THEN
                    COM=COM+RQQ(NSE)
                    R=R-SC(NSE)*SQ(NSE)
                ELSE
                    COM=COM+RQQ(NSE)*BSECM(NSM,2)
                    R=R-RQQ(NSE)*BSECM(NSM,1)
                END IF
            ELSE
                NSE=ICHASE(2*NTC)
                NSH=ICHASE(2*NTC-1)
                NSM=NSE-NTC-MBSEC0
                IF (NTT.GT.I .AND. NCOM.EQ.0) THEN
                    NCOM=NM0
                    NM0=NM0+1
                END IF
	            IGN=MSIGN(J)
                KGN=ICHADI(NTC)
                IF (IGN .GT. 0) THEN
                    Q=RQQ(NSE)
                    IF (KGN .GT. 0) THEN
                        COM=COM+Q*ISECM(NSM,3)
                        R=R-Q*ISECM(NSM,1)
                        TUNMAT(NM0)=Q*ISECM(NSM,2)
                    ELSE
                        COM=COM+Q
                        TUNMAT(NM0)=0.
                        R=R-SC(NSE)*SQ(NSE)
                    END IF
                ELSE
                    Q=RQQ(NSH)
                    IF (KGN .GT. 0) THEN
                        COM=COM-Q
                        TUNMAT(NM0)=0.
                        R=R-SC(NSH)*SQ(NSH)
                    ELSE
                        COM=COM-Q*ISECM(NSM,3)
                        R=R+Q*ISECM(NSM,1)
                        TUNMAT(NM0)=-Q*ISECM(NSM,2)
                    END IF
                END IF
                ITUNM1(NM0)=NTT
                NM0=NM0+1
            END IF
200     CONTINUE
        IF (NCOM .EQ. 0) THEN
            NCOM=NM0
	        NM0=NM0+1
    END IF
    TUNMAT(NCOM)=COM
    ITUNM1(NCOM)=I
    RTUNM(I)=R
    GOTO 100
END