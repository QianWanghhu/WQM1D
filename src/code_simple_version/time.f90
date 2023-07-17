SUBROUTINE CLOCK(DT,SE,MI,HR,DA,MN,YR)
    INTEGER SE,MI,HR,DA,MN,YR
    NT=INT(DT+0.01)
    SE=SE+NT
    IF (SE .GE. 60) THEN
        NT=SE/60
        SE=SE-NT*60
        MI=MI+NT
        IF (MI .GE. 60) THEN
        NT=MI/60
        MI=MI-NT*60
        HR=HR+NT
            IF (HR .GE. 24) THEN
                NT=HR/24
                HR=HR-NT*24
                DA=DA+NT
                NT=JMD(YR,MN)
                IF (DA .GT. NT) THEN
                    DA=DA-NT
                    MN=MN+1
                    IF (MN .EQ. 13) THEN
                        MN=1
                        YR=YR+1
                    END IF
                END IF
            END IF
        END IF
    END IF
END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          JMD                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
INTEGER FUNCTION JMD(YR,MN)
	INTEGER A,B,C,D,MN,YR
	A=INT(MN/2.)
	B=INT(MN/8.)
	C=INT(YR/4.)
	C=INT((C*4+4-YR)/4.)
	D=INT((12-MN)/10.)
	JMD=30-D*A*(2-C)+(MN-2*A)*(1-2*B)+B
END