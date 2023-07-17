SUBROUTINE RBISEC(NCHAI0,MBSEC0)
    ! c     	  REAL ISECM
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
    DO 300 I=1,NCHAI0
        J=2*I
        NSH=ICHASE(J-1)
        NSE=ICHASE(J)
        IF (ICHADI(I) .GT. 0) THEN
            C1=C(NSH)
            NSM=NSE-I-MBSEC0-1
            DO 100 J=NSE-1,NSH+1,-1
                C(J)=ISECM(NSM,1)+ISECM(NSM,2)*C1+ISECM(NSM,3)*C(J+1)
                NSM=NSM-1
100         CONTINUE
        ELSE
            C1=C(NSE)
            NSM=NSE-I-MBSEC0-1
            DO 200 J=NSH+1,NSE-1
                C(J)=ISECM(NSM,1)+ISECM(NSM,2)*C1+ISECM(NSM,3)*C(J-1)
                NSM=NSM-1
200         CONTINUE
        END IF
300 CONTINUE
END