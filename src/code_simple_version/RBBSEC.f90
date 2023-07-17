SUBROUTINE RBBSEC(NBC0)
    ! c  	  REAL ISECM
    ! c        INTEGER BCHASE,BCHADI
    ! c*******************************************************************
    ! c*******************************************************************
    ! c        COMMON /BC/     RBC(13,4296),DBC(13)
    ! c     &,        /INIT/   RHH(49),RQQ(49),C(49)
    ! c     &,        /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
    ! c     &,        /SSSSSS/ BSECM(16,2),ISECM(10,3)
    ! c     &,        /TUN/    NTUNSS(12),NTUNSE(33),MSIGN(33),NTUNCH(33)
    ! c     &,                 NTUNTU(33)
    ! c     &,        /TUNM/   TUNMAT(31),RTUNM(11),ITUNM1(31),ITUNM2(11)
    ! c**********************************************************************
    use global
    DO 300 I=1,NBC0
        NSE=BCHASE(2*I) 
        NSH=BCHASE(2*I-1)
        IF (BCHADI(I) .GT. 0) THEN
            NSM=NSE-I-1
            C(NSH)=DBC(I)
            DO 100 J=NSE-1,NSH+1,-1
                C(J)=BSECM(NSM,1)+BSECM(NSM,2)*C(J+1)
                NSM=NSM-1
100         CONTINUE
        ELSE
            NSM=NSE-I
            DO 200 J=NSE-1,NSH,-1
                C(J)=BSECM(NSM,1)+BSECM(NSM,2)*C(J+1)
                NSM=NSM-1
200          CONTINUE
        END IF
300 CONTINUE
END