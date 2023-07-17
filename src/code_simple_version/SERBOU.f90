
SUBROUTINE SERBOU(X,DX0,NDX,N)
   ! c*******************************************************************
   ! c        COMMON /BC/     RBC(13,4296),DBC(13)
   ! c**********************************************************************
   use global
   DO 100 I=1,N !Loop over boundaries.
      D=DX0 !3600s
      NX=INT(X/D+1.E-4)
      DDX=X-REAL(NX)*D !
      DBC(I)=RBC(I,NX+1)+(RBC(I,NX+2)-RBC(I,NX+1))*DDX/D  !Linear interpolation of initial concentrations.
   100	CONTINUE
END