

SUBROUTINE SOLVER(A,IA,MA,B,M,N,EPS,IE)
	use global
	DIMENSION A(M),IA(M),MA(N),B(N), Z(83),MV(83),V(5000)
	IK=1
	DO 10 I=1,N
		DO 21 J=1,N
21	   		Z(J)=0.
	   	C=0.
	   	MI=MA(I)
	   	DO 20 L=1,MI
	      	J=IA(IK)
	      	Z(J)=A(IK)
20	 		IK=IK+1
	   	IL=1
	   	IF(I.EQ.1)GOTO 15
	   	KK=I-1
	   	DO 30 L=1,KK
	      	MI=MA(L)
	      	IF(Z(MI).EQ.0)GOTO 40
	      	ML=MV(L)
	      	IF(ML.EQ.0)GOTO 45
	      	DO 50 K=1,ML
				J=INT(V(IL+1))
				Z(J)=Z(J)-Z(MI)*V(IL)
50	      		IL=IL+2
45	      	B(I)=B(I)-Z(MI)*B(L)
	      	Z(MI)=0.
	      	GOTO 30
40	      	IL=IL+2*MV(L)
30	   	CONTINUE

15	   	L=0
	   	DO 11 J=1,N
	      	IF(Z(J).EQ.0)GOTO 11
	      	IF(ABS(Z(J)).LE.ABS(C))GOTO 11
	      	C=Z(J)
	      	J0=J
11	   	CONTINUE
		IF (ABS(C) .LE. EPS) THEN
            WRITE(*,*) CHAR(7)
            WRITE(*,*) CHAR(7),'There is a error in SOLVER : IE=1'
            PAUSE ' '
            RETURN
		END IF
 	   	MA(I)=J0
	   	DO 60 J=1,N
	      	IF(Z(J).EQ.0)GOTO 60
	      	IF(J.EQ.J0)GOTO 60
	      	V(IL)=Z(J)/C
	      	V(IL+1)=FLOAT(J)
	      	IL=IL+2
	      	L=L+1
60	   	CONTINUE
		MV(I)=L
10		B(I)=B(I)/C
	DO 80 I=1,N
		IJ=N-I
		L=MA(IJ+1)
		Z(L)=B(IJ+1)
		IF(I.EQ.N)GOTO 80
		MI=MV(IJ)
		IF(MI.EQ.0)GOTO 80
	   	DO 70 L=1,MI
	      	IL=IL-2
	      	J=INT(V(IL+1))
70	   		B(IJ)=B(IJ)-Z(J)*V(IL)
80	CONTINUE
	DO 110 I=1,N
110		B(I)=Z(I)
	IE=0
END