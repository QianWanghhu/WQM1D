SUBROUTINE CALTRAN (ISTL_,IS2TL_,MVAR,MO,CON,CON1,IT)  

  ! **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO  
  ! **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! **  THE NUMBER OF TIME LEVELS IN THE STEP  

  USE GLOBAL  
  !  USE OMP_LIB
  
  IMPLICIT NONE  
    
  INTEGER, INTENT(IN) :: ISTL_,IS2TL_,MVAR,MO,IT  
  REAL    :: CON(LCM,KCM),CON1(LCM,KCM)  
  
  REAL    :: BSMALL, DDELT
  REAL    :: CTMP, CBT, AUHU, AVHV, UTERM, VTERM, WTERM  
  REAL    :: CBSTMP, CBWTMP, CBETMP, CBNTMP, UHU, VHV, AWW, WW  
  REAL    :: CWMAX, CEMAX, CSMAX, CNMAX, CMAXT  
  REAL    :: CWMIN, CEMIN, CSMIN, CNMIN, CMINT  
  
  INTEGER :: M
  INTEGER :: ISUD, K, NSID, IOBC, NMNLOD
  INTEGER :: LP,L, LN, LS, LE, LW, LSE, LNW, LL, ITRANFLOC  

  ! *** SET UP FLOC TRANSPORT
  ITRANFLOC=0
  
  BSMALL=1.0E-6  
  ISUD=1  ! Qian: Flag for determining whether conc is updated
  DDELT=DT 
    
  M=MO  
  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      WQBCCON(IOBC,K,IT)  = CON(L,K)  
      WQBCCON1(IOBC,K,IT) = CON1(L,K)  
    ENDDO  
  ENDDO  
  
  ! **  CALCULATED EXTERNAL SOURCES AND SINKS  
  CALL CALFQC (ISTL_,IS2TL_,MVAR,M,CON,CON1,IT)  
  
  ! **  SELECT TRANSPORT OPTION, ISPLIT=1 FOR HORIZONTAL-VERTICAL  
  ! **  OPERATOR SPLITTING  
  ! **  BEGIN COMBINED ADVECTION SCHEME  
  ! **  ADVECTIVE FLUX CALCULATION  

  !Qian: only one option for 1D
  IF( ISTL_ == 2 ) GOTO 300  
  ! *** UHDY2 AND VHDX2 ARE LAYER FLOWS (SIGMA-Z VERSION)

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED  
  ! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY  

  300 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      FUHUD(L,K,IT) = UHDY2(L,K)*CON1(LUPU(L,K),K)  
      FVHUD(L,K,IT) = VHDX2(L,K)*CON1(LUPV(L,K),K)  
    ENDDO  
  ENDDO
  ! KC == 1 for 1D and thus code for KC > 1 is not needed blow 
  GOTO 500  

  ! **  STANDARD ADVECTION CALCULATION  
  500 CONTINUE  

  ! *** CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX (PMC MOVED)  
  ! TODO: work on 1D diffusive
  IF( ISHDMF == 2 ) CALL CALDIFF (CON1,IT)
    
  ! *** BEGIN IF ON TRANSPORT OPTION CHOICE  
  ! *** IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE  
  !  ISCDCA:  0 FOR STANDARD DONOR CELL UPWIND DIFFERENCE ADVECTION (3TL ONLY)
  !           1 FOR CENTRAL DIFFERENCE ADVECTION FOR THREE TIME LEVEL STEPS (3TL ONLY)
  !           2 FOR EXPERIMENTAL UPWIND DIFFERENCE ADVECTION (FOR RESEARCH) (3TL ONLY)
  IF( ISCDCA(MVAR) == 0 )THEN
    DO K=1,KC  
    DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                    FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                                + (FWUU(L,K-1,IT)-FWUU(L,K,IT)) )  
    ENDDO  
    ENDDO  

    ! 1 TO ADD FLUX LIMITING TO ANTI-NUMERICAL DIFFUSION CORRECTION
    IF( ISFCT(MVAR) >= 1 .AND. ISADAC(MVAR) > 0 )THEN 
    DO K=1,KC  
        DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CON2(L,K,IT) = MAX(CON1(L,K),0.0)
        ENDDO  
    ENDDO  
    ENDIF  
    
  
    IF( IS2TL_ == 0  .AND.  ISUD == 1 )THEN  
      ! *** ADVANCE CON1 TO CON (3TL)
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CON1(L,K) = CON(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  
      
    ! *** UPDATE NEW CONCENTRATIONS  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)
        CON(L,K) = CD(L,K,IT)*HPKI(L,K)
      ENDDO
    ENDDO  
    
  !Qian: Deleted code for TRANSPORT OPTION CHOICE: ISCDCA(MVAR)/=0
  ! *** ELSE ON TRANSPORT OPTION CHOICE: ISCDCA(MVAR)/=0
    
  ENDIF ! *** END OF TRANSPORT OPTION CHOICE  

  ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S - 2TL & 3TL
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      CON1(L,K) = WQBCCON1(IOBC,K,IT)
    ENDDO  
  ENDDO  

  ! *** ALL OTHER WATER CONSTITUENTS  
  IF( MVAR == 8 )THEN  ! .AND. IWQPSL == 2 )THEN
    M=4+NTOX+NSED+NSND+MO  
  ENDIF
  
  ! ******************************************************************************************
  ! *** APPLY OPEN BOUNDARY CONDITIONS, BASED ON DIRECTION OF FLOW  

  ! *** SOUTH OPEN BC, WITHOUT FLOCS
  IF( NCBS > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBS  
        NSID=NCSERS(LL,M)  
        L=LCBS(LL) 
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LN=LNC(L)  
        IF( VHDX2(LN,K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K)+CON(L,K)) + 0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)  &  
                                           + DDELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          IF( M == 1 )THEN
            ! *** LIMIT CONCENTRATIONS TO MAXIMUM BC CONCENTRATIONS AT BOTTOM LAYER (SALINITY ONLY)
            CBSTMP = CBS(LL,1,M) + CSERT(1,NSID,M)  
            IF( CON(L,K) > CBSTMP )THEN
              CON(L,K) = CBSTMP  
            ENDIF
          ENDIF  
          CLOS(LL,K,M)=CON(L,K)  
          NLOS(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLOS(LL,K,M)  
          IF( NMNLOD >= NTSCRS(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBSTMP = CLOS(LL,K,M) + (CBT-CLOS(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRS(LL))  
            CON(L,K) = MAX(CBSTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** WEST OPEN BC, WITHOUT FLOCS  
  IF( NCBW > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBW  
        NSID=NCSERW(LL,M)  
        L=LCBW(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        IF( UHDY2(LEC(L),K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR)/=2) CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K)
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)         &  
                               +DDELT*(0.5*UHDY2(LEC(L),K)*(CON1(L,K)+CON(L,K))-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBWTMP=CBW(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBWTMP) CON(L,K) = CBWTMP  
          CLOW(LL,K,M)=CON(L,K)  
          NLOW(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLOW(LL,K,M)  
          IF( NMNLOD >= NTSCRW(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBWTMP=CLOW(LL,K,M)+(CBT-CLOW(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRW(LL))  
            CON(L,K) = MAX(CBWTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** EAST OPEN BC, WITHOUT FLOCS  
  IF( NCBE > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBE  
        NSID=NCSERE(LL,M)  
        L=LCBE(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) ) CYCLE 
        IF( UHDY2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2) CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)+DDELT*(FUHUD(L,K,IT) &  
                                      -0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBETMP = CBE(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBETMP) CON(L,K) = CBETMP  
          CLOE(LL,K,M)=CON(L,K)  
          NLOE(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBE(LL,1,M) + WTCI(K,2)*CBE(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER-NLOE(LL,K,M)  
          IF( NMNLOD >= NTSCRE(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBETMP = CLOE(LL,K,M) + (CBT-CLOE(LL,K,M)) * FLOAT(NMNLOD)/FLOAT(NTSCRE(LL))  
            CON(L,K) = MAX(CBETMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF
    
  ! *** NORTH OPEN BC, WITHOUT FLOCS  
  IF( NCBN > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBN  
        NSID=NCSERN(LL,M)  
        L=LCBN(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LS=LSC(L)  
        IF( VHDX2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR)/=2) CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K) 
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2PK(L,K)*HPKI(L,K) + DDELT*(FVHUD(L,K,IT) &  
                                      -0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBNTMP=CBN(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBNTMP) CON(L,K) = CBNTMP  
          CLON(LL,K,M)=CON(L,K)  
          NLON(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLON(LL,K,M)  
          IF( NMNLOD >= NTSCRN(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBNTMP=CLON(LL,K,M)+(CBT-CLON(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRN(LL))  
            CON(L,K) = MAX(CBNTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! ****************************************************************************************
  ! **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATIONS WITH FLUX CORRECTOR  
  ! Qian: GOTO 2000
  
  ! ----------------------------------------------------------------------------------------
  ! *** CALTRAN EXIT 
 2000 CONTINUE  
 
  ! *** ZERO HEAT FLUXES
  IF( MVAR == 2 )THEN        
    ! *** ZERO EVAP/RAINFALL
    DO L=1,LC  
      FQC(L,KC,IT) = 0.  
    ENDDO  
    IF( ISADAC(MVAR) >= 2 )THEN
      DO L=1,LC  
        FQCPAD(L,KC,IT) = 0.  
      ENDDO  
    ENDIF
    IF( ISADAC(MVAR) > 0 )THEN
      DO L=1,LC  
        QSUMPAD(L,KC,IT) = 0.  
      ENDDO  
    ENDIF
  ENDIF

  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  

  RETURN  
END  
