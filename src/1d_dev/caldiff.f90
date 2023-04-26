SUBROUTINE CALDIFF (CON1,IT)
                                                                                                                         
  ! **  SUBROUTINE CALDIFF CALCULATES THE HORIZONTAL DIFFUSIVE                                                            
  ! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO                                                     
  ! **  A REVISEDED VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL                                                          
  ! **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP

  ! *** VARIABLES   DESCRIPTION                                 UNITS
  ! *** FUHU,FVHU   GROSS MOMENTUM, U COMPONENTS                M4/S2
  ! *** FVHV,FUHV   GROSS MOMENTUM, V COMPONENTS                M4/S2
  ! *** FX, FY      INTERNAL MODE FORCING BY LAYER              M4/S2
  ! *** FBBX, FBBY  INTERNAL MODE BOUYANCY FORCING BY LAYER     M4/S2
  ! *** FCAX,FCAY   CORIOLIS FORCING BY LAYER                   M4/S2
  ! *** DU, DV      INTERNAL SHEARS BY LAYER                    M2/S2    
  ! AH: Horizontal Turbulent Viscosity, depth normalized (m/s)
  ! DYU: 0.5*(DYP+DYP(LS)) (m)
  ! DXIU: 1./DXU (1/m)
  ! DXIV: 1./DXV (1/m)
  ! DXU: 0.5*(DXP+DXP(L-1)) (m)
  ! DXV: 0.5*(DXP+DXP(LS)) (m)
  ! FUHUD(L,K,IT) = UHDY2(L,K)*CON1(LUPU(L,K),K)  
  ! FVHUD(L,K,IT) = VHDX2(L,K)*CON1(LUPV(L,K),K)  



  ! CHANGE RECORD      
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH

  USE GLOBAL
  IMPLICIT NONE   
  INTEGER,INTENT(IN) :: IT                                                                                                                                                                                                       
  REAL,   INTENT(IN) :: CON1(LCM,KCM)                                                                                                          
  INTEGER            :: ND,K,LP,L,LS,LW

  ! **  HORIZONTAL DIFFUSIVE FLUX CALCULATION                                                                             
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLHDMF(K,ND)
        L = LKHDMF(LP,K,ND)
        LW=LWC(L)
        FUHUD(L,K,IT)=FUHUD(L,K,IT)+0.5*SUB3D(L,K)*DYU(L)*HU(L)*DZC(L,K)*(AH(L,K)+AH(LW,K))*(CON1(LW,K)-CON1(L,K))*DXIU(L)
        LS=LSC(L)
        FVHUD(L,K,IT)=FVHUD(L,K,IT)+0.5*SVB3D(L,K)*DXV(L)*HV(L)*DZC(L,K)*(AH(L,K)+AH(LS,K))*(CON1(LS,K)-CON1(L,K))*DYIV(L)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE

