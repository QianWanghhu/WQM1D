SUBROUTINE CALWQC(ISTL_,IS2TL_)  

  ! CHANGE RECORD  
  ! **  SUBROUTINE CALWQC CALCULATES THE CONCENTRATION OF DISSOLVED AND  
  ! **  SUSPENDED WATER QUALITY CONSTITUTENTS AT TIME LEVEL (N+1).  
  ! **  CALLED ONLY ON ODD THREE TIME LEVEL STEPS  
  ! **  DS-INTL: 2010-12  ADDED OMP, CLEANED UP CODE AND CONVERTED TO F90
  ! 2014-08           D H CHUNG        SET EXPLICIT PRECISIONS OF INTEGER & REAL

  USE GLOBAL  
  USE OMP_LIB
  
  IMPLICIT NONE
  
  INTEGER :: ISTL_,IS2TL_
  INTEGER :: L,K,IP,ITHD,NW,ND,LF,LL,LP
  
  INTEGER, SAVE             :: NACTIVEWQ
  INTEGER, SAVE,ALLOCATABLE :: IACTIVEWQ(:)
  
  REAL(RKD)   :: DSTIME
  REAL(RKD)   :: RCDZKK, RCDZKMK, CCLBTMP, CCUBTMP, CCMBTMP, EEB
  REAL(RKD)   :: TTDS         ! MODEL TIMING TEMPORARY VARIABLE
  
  IF( .NOT. ALLOCATED(IACTIVEWQ) )THEN
    ALLOCATE( IACTIVEWQ(NWQV) )
    IACTIVEWQ = 0
    NACTIVEWQ = 0
    
    DO NW = 1,NWQV
      IF( ISTRWQ(NW) == 1 )THEN
        NACTIVEWQ = NACTIVEWQ+1
        IACTIVEWQ(NACTIVEWQ) = NW
      ENDIF
    ENDDO
  ENDIF
  
  DELT=DT2  
  IF( IS2TIM >= 1 )THEN  
    IF( ISDYNSTP == 0 )THEN  
      DELT=DT  
    ELSE  
      DELT=DTDYN  
    END IF  
  ENDIF  
  ITHD=1
  
  ! **  UPDATED TIME SERIES CONCENTRATION BOUNDARY CONDITIONS  
  ! **  3D ADVECTI0N TRANSPORT CALCULATION  
  TTDS=DSTIME(0)
  DO ND=1,NACTIVEWQ  
    NW=IACTIVEWQ(ND)
    CALL CALTRAN(ISTL_,IS2TL_,8,NW,WQV(1,1,NW),WQV(1,1,NW),ITHD)
    
    IF( ISICE > 2 .AND. NW == 19 )THEN !Qian : ISICE determines whether to use ice sub-module
      ! *** ZERO SURFACE MELT FLUX
      DO L=1,LC  
        FQC(L,KC,ITHD)=0.  
      ENDDO  
    ENDIF
  ENDDO  
  
  TWQADV=TWQADV+(DSTIME(0)-TTDS)  

  ! **  CALLS TO SOURCE-SINK CALCULATIONS  
  ! **  BYPASS OR INITIALIZE VERTICAL DIFFUSION CALCULATION  
  IF( KC == 1 ) GOTO 2000  

  ! **  VERTICAL DIFFUSION CALCULATION LEVEL 1  
  IF( ISWQLVL >= 1 .AND. ISWQLVL <= 3 )THEN
    TTDS=DSTIME(0)  

    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      ! *** BOTTOM LAYER
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        RCDZKK=-DELT*CDZKK(L,KSZ(L))  
        CCUBTMP=RCDZKK*HPI(L)*AB(L,KSZ(L))  
        CCMBTMP=1._8-CCUBTMP  
        EEB=1._8/CCMBTMP  
        CU1(L,KSZ(L))=CCUBTMP*EEB  
        DO IP=1,NWQV
          WQV(L,KSZ(L),IP) = WQV(L,KSZ(L),IP)*EEB  
        ENDDO
      ENDDO  
    
      ! *** MIDDLE LAYERS
      DO K=2,KS  
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND) 
          RCDZKMK=-DELT*CDZKMK(L,K)  
          RCDZKK=-DELT*CDZKK(L,K)  
          CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)  
          CCUBTMP=RCDZKK*HPI(L)*AB(L,K)  
          CCMBTMP=1._8-CCLBTMP-CCUBTMP  
          EEB=1._8/(CCMBTMP-CCLBTMP*CU1(L,K-1))  
          CU1(L,K)=CCUBTMP*EEB  
          DO IP=1,NWQV
            WQV(L,K,IP) = (WQV(L,K,IP) - CCLBTMP*WQV(L,K-1,IP))*EEB  
          ENDDO
        ENDDO  
      ENDDO  
    
      ! *** TOP LAYER !Qian: only need top layer.
      K=KC  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        RCDZKMK=-DELT*CDZKMK(L,K)  
        CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)  
        CCMBTMP=1._8-CCLBTMP  
        EEB=1._8/(CCMBTMP-CCLBTMP*CU1(L,K-1))  
        DO IP=1,NWQV
          WQV(L,K,IP) = (WQV(L,K,IP) - CCLBTMP*WQV(L,K-1,IP))*EEB  
        ENDDO
      ENDDO  

      ! *** FINAL PASS
      DO IP=1,NWQV
        DO K=KS,1,-1  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQV(L,K,IP) = WQV(L,K,IP) - CU1(L,K)*WQV(L,K+1,IP)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    TWQDIF=TWQDIF+(DSTIME(0)-TTDS)  
    
  ENDIF 

2000 CONTINUE  

  RETURN  
  
  END  

