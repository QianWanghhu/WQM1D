SUBROUTINE CHAINF(NBC0,NCHAI0)
    use global
    ! c        INTEGER BCHASE,BCHADI
    ! c*******************************************************************
    ! c        COMMON /CHA/    BCHASE(26),ICHASE(20),BCHADI(13),ICHADI(10)
    ! c**********************************************************************
    BCHASE(1)=0 ! The initial index of outer chain sections.
    ICHASE(1)=BCHASE(NBC0+1) !The initial index of inner chain sections, which is the index of the last outer section.
    DO 10 I=NBC0,1,-1
        BCHASE(2*I)=BCHASE(I+1) !The continuous section index of outer channels. This line assigns the index of end sections.
10      BCHASE(2*I-1)=BCHASE(I)+1 ! This line assigns the index of start sections.

    DO 30 I=NCHAI0,1,-1
	    ICHASE(2*I)=ICHASE(I+1) ! This does the same index assignments for inner channel sections.
30	    ICHASE(2*I-1)=ICHASE(I)+1

END