
SUBROUTINE SERSEC(NS,A,B)
   ! c*******************************************************************
   ! c        COMMON  /INIT/   RHH(49),RQQ(49),C(49)
   ! c     &,        /CHAP/   SECFOR(49,3),DX(49)
   ! c**********************************************************************
   use global
   DH=RHH(NS)-SECFOR(NS,1) ! Calculate the water depth.
   B=2.*DH*SECFOR(NS,3)+SECFOR(NS,2) ! Calculate the discharge section width.
   A=0.5*(B+SECFOR(NS,2))*DH ! Calculate the discharge section area.
END