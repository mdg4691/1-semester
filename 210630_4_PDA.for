      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C     
      INCLUDE 'ABA_PARAM.INC'
C     
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1       DDSDDE(NTENS,NTENS),
     2       DDSDDT(NTENS),DRPLDE(NTENS),
     3       STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4       PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      DIMENSION STRANT(6)
      DIMENSION CFULL(6,6), CDFULL(6,6)
      DIMENSION DDFT(6),DDFC(6),DDDMT(6),DDDMC(6),
     1 DDOMT(6),DDOMC(6)
      DIMENSION DDFF(6),DDDMF(6),DDOMF(6),DDFMS(6)
      DIMENSION DCDFT(6),DCDFC(6),DCDDMT(6),DCDDMC(6),
     1 DCDOMT(6),DCDOMC(6)
      DIMENSION DCDFF(6,6),DCDDMF(6,6),DCDOMF(6,6),DCDFMS(6,6)
      DIMENSION ATEMPA(6), ATEMPB(6), ATEMPC(6), ATEMPD(6)
      DIMENSION OLD_STRESS(6)
      DIMENSION DOLD_STRESS(6),D_STRESS(6)
      PARAMETER (ZERO = 0.D0, ONE=1.D0, TWO=2.D0, HALF = 0.5D0)
C
C
C      GET THE MATERIAL PROPERTIES---USER MATERIAL
C	  
      EON = PROPS(1)        !E1
      ETW = PROPS(2)        !E2
      ETH = ETW             !E3
      GOT = PROPS(3)        !G12
      GOTH = GOT            !G13
      GTT = PROPS(4)        !G23
      NUOT = PROPS(5)       !V12
      NUTTH = PROPS(6)      !V23
      NUTO = NUOT / EON * ETW       !V21
C
C     GET THE FAILURE PROPERTIES
C
      XT = PROPS(7)         ! 1 DIRECTION TENSILE STRENGTH 
      XC = PROPS(8)         ! 1 DIRECTION COMPRESSIVE STRENGTH
      YT = PROPS(9)        ! 2 DIRECTION TENSILE STRENGTH 
      YC = PROPS(10)        ! 2 DIRECTION COMPRESSIVE STRENGTH
      ZT = YT               ! 3 DIRECTION TENSILE STRENGTH
      ZC = YC               ! 3 DIRECTION COMPRESSIVE STRENGTH
      SOT = PROPS(11)       ! SHEAR STRENGTH 12
      SOTH = SOT            ! SHEAR STRENGTH 13
      STTH = PROPS(12)      ! SHEAR STRENGTH 23
      GF = PROPS(13)        !FIBER FRACTURE ENERGY
      GM = PROPS(14)        !MATRIX FRACTURE ENERGY
      ETA = PROPS(15)       !VISCOSITY
C
C      CALCULATE THE STRAIN AT THE END OF THE INCREMENT
C
      DO I = 1, NTENS
	    STRANT(I) = STRAN(I) + DSTRAN(I)
      END DO
C
C      FILL THE 6X6 FULL STIFFNESS MATRIX
C
       DO I = 1, 6
          DO J = 1, 6
		   CFULL(I,J)=ZERO
          END DO
       END DO
      ATEMP = ONE - (TWO * NUOT * NUTO) - (NUTTH ** TWO)
     1     - (TWO * NUOT * NUTO * NUTTH)
      CFULL(1,1) = EON * (ONE - NUTTH ** TWO) / ATEMP
      CFULL(2,2) = ETW * (ONE - NUOT * NUTO) / ATEMP
      CFULL(3,3) = CFULL(2,2)
      CFULL(1,2) = ETW * (NUTO + NUOT * NUTTH) / ATEMP
      CFULL(1,3) = CFULL(1,2)
      CFULL(2,3) = ETW * (NUTTH + NUOT * NUTO) / ATEMP
      CFULL(4,4) = (CFULL(2,2)-CFULL(2,3))/TWO
      CFULL(5,5) = GOTH
      CFULL(6,6) = GTT
      DO I = 2,6
         DO J = 1,I-1
		   CFULL(I,J) = CFULL(J,I)
         END DO
      END DO
C
C  CALCULATE THE FAILURE STRAIN BY FAILURE STRESS
C
      XTFE = XT / CFULL(1,1)           !FAILURE STRAIN 1 DIRECTION IN TENSION
      XCFE = XC / CFULL(1,1)           !FAILURE STRAIN 1 DIRECTION IN COMPRESSION
      YTFE = YT / CFULL(2,2)           !FAILURE STRAIN 2 DIRECTION IN TENSION
      YCFE = YC / CFULL(2,2)           !FAILURE STRAIN 2 DIRECTION IN COMPRESSION
      ZTFE = ZT / CFULL(3,3)           !FAILURE STRAIN 3 DIRECTION IN TENSION
      ZCFE = ZC / CFULL(3,3)           !FAILURE STRAIN 3 DIRECTION IN COMPRESSION
      XCFMO = XC / CFULL(1,1)          !FAILURE STRAIN FIBER MATRIX SHEAR OUT
      SOTFE = SOT / CFULL(4,4)         !FAIURE SHEAR STRAIN 12
      SOTHFE = SOTH / CFULL(5,5)       !FAIURE SHEAR STRAIN 13
      STTFE = STTH / CFULL(6,6)        !FAIURE SHEAR STRAIN 23
C
C    CHECK THE FAILURE INITIATION CONDITION
C
      FTNOLD = STATEV(1)
      FCNOLD = STATEV(2)
      DMTNOLD = STATEV(3)
      DMCNOLD = STATEV(4)
      OMTNOLD = STATEV(5)
      OMCNOLD = STATEV(6)
      FNSNOLD = STATEV(7)
      DFTOLD = STATEV(8)
      DFCOLD = STATEV(9)
      DDMTOLD = STATEV(10)
      DDMCOLD = STATEV(11)
      DOMTOLD = STATEV(12)
      DOMCOLD = STATEV(13)
      DFMSOLD = STATEV(14)
      DFTVOLD = STATEV(15)
      DFCVOLD = STATEV(16)
      DDMTVOLD = STATEV(17)
      DDMCVOLD = STATEV(18)
      DOMTVOLD = STATEV(19)
      DOMCVOLD = STATEV(20)
      DFMSVOLD = STATEV(21)   
C	  
C
C
      CALL CheckFailureIni(XTFE,XCFE,YTFE,YCFE,ZTFE,ZCFE,XCFMO,
     1  SOTFE,SOTHFE,STTFE,STRANT,GF,GM, CELENT, CFULL,DFT, DFC, DFMS, 
     2 DDMT, DDMC, DOMT, DOMC, DDFT, DDFC, DDFMS, DDDMT, DDDMC, DDOMT,
     3 DDOMC, NTENS, DFTOLD, DFCOLD, DDMTOLD, DDMCOLD, DOMTOLD,
     4 DOMCOLD, DFMSOLD, DDFF, DDDMF, DDOMF,
     5 FTN,FCN,DMTN,DMCN,OMTN,OMCN,FMSN)
C
C
C       DFFV, DDMFV, DOMFV, DFMSV
C
        DFF = ZERO
        DDMF = ZERO
        DOMF = ZERO
      DO I = 1,6
          DDFF(I) = ZERO
          DDDMF(I) = ZERO
          DDOMF(I) = ZERO
      END DO
C
      IF (DFT .GE. ZERO) THEN
	  DFF = DFT
      DO I = 1,6
	  DDFF(I) = DDFT(I)
      END DO
      ELSE IF (DFC .GE. ZERO) THEN
	  DFF = DFC
      DO I = 1,6
	  DDFF(I) = DDFC(I)
      END DO
      END IF
C
      IF (DDMT .GE. ZERO) THEN
	  DDMF = DDMT
      DO I = 1,6
	  DDDMF(I) = DDDMT(I)
      END DO
      ELSE IF (DDMC .GE. ZERO) THEN
	  DDMF = DDMC
      DO I = 1,6
	  DDDMF(I) = DDDMC(I)
      END DO
      END IF
C
      IF (OMT .GE. ZERO) THEN
	  DOMF = DOMT
      DO I = 1,6
	  DDOMF(I) = DDOMT(I)
      END DO
      ELSE IF (OMC .GE. ZERO) THEN
	  DOMF = DOMC
      DO I = 1,6
	  DDOMF(I) = DDOMC(I)
      END DO
      END IF
      DFF = MAX (DFF, DFFOLD)
      DDMF = MAX (DDMF, DDMFOLD)
      DOMF = MAX(DOMF, DOMFOLD)
      DFMS = MAX(DFMS, DFMSLOD)
C
C     VISCOUS REGULARIZATION
C
      DFFV = ETA / (ETA + DTIME) * DFFVOLD + DTIME / (ETA + DTIME) * DFF
      DFMSV = ETA / (ETA + DTIME) * DFMSVOLD + DTIME / (ETA + DTIME)
     1 * DFMS
      DDMFV = ETA / (ETA + DTIME) * DDMFVOLD + DTIME / (ETA + DTIME)
     1 * DDMF
      DOMFV = ETA / (ETA + DTIME) * DOMFVOLD + DTIME / (ETA + DTIME)
     1 * DOMF
C
C     SAVE THE OLD STRESS TO OLD_STRESS
C
      DO I = 1, NTENS
	     OLD_STRESS(I) = STRESS(I)
      END DO
C
C     CALL ROUTINE TO CALCULATE THE STRESS
C
C     CALCULATE THE STRESS IF THERE'S NO VISCOUS REGULARIZATION
      CALL GetStress(CFULL,CDFULL,DFF,DDMF,DOMF,DFMS,D_STRESS,STRANT,
     1 NTENS)
C
C     CALCULATE THE STRESS IF THERE'S VISCOUS REGULARIZATION
      CALL GetStress(CFULL,CDFULL,DFFV,DDMFV,DOMFV,DFMSV,STRESS,STRANT,
     1 NTENS)
C     GET THE OLD STRESS IF THERE'S NO VISCOUS REGULARIZATION
      DO I=1,NTENS
	     DOLD_STRESS(I)=STATEV(I+27)
      END DO
C
C     SAVE THE CURRENT STRESS
      DO I=1,NTENS
	     STATEV(I+27)=D_STRESS(I)
      END DO
C
C     CALCULATE THE DERIVATIVE MATRIX DAMAGE MATRIX
C
      CALL ElasticDerivative(CFULL,DFFV,DDMFV,DOMFV,DFMSV, DCDFF,DCDDMF,
     1 DCDOMF,DCDFMS)
C
C     UPDATE THE JACOBIAN
C
      DO I = 1, NTENS
	    ATEMPA(I) = ZERO
		  DO J = 1, NTENS
		   ATEMPA(I) = ATEMPA(I) + DCDFF(I,J) * STRANT(J)
		  END DO
      END DO
C
      DO I = 1, NTENS
	    ATEMPB(I) = ZERO
		DO J = 1, NTENS
		   ATEMPB(I) = ATEMPB(I) + DCDDMF(I,J) * STRANT(J)
		END DO
	 END DO
C
      DO I = 1, NTENS
	    ATEMPC(I) = ZERO
		DO J = 1, NTENS
		   ATEMPC(I) = ATEMPC(I) + DCDOMF(I,J) * STRANT(J)
		END DO
	 END DO
C
      DO I = 1, NTENS
	    ATEMPD(I) = ZERO
		DO J = 1, NTENS
		   ATEMPD(I) = ATEMPD(I) + DCDFMS(I,J) * STRANT(J)
		END DO
	 END DO
C
      DO I = 1, NTENS
	    DO J = 1, NTENS
		   DDSDDE(I,J)=CDFULL(I,J) + (ATEMPA(I) * DDFF(J)
     1        + ATEMPB(I) * DDDMF(J) + ATEMPC(I) * DDOMF(J)
     2        + ATEMPD(I) * DDFMS(J)) * DTIME / (DTIME + ETA)
		END DO
      END DO
C
C     TO UPDATE THE STATE VARIABLE
C
      STATEV(1) = FTN
      STATEV(2) = FCN
      STATEV(3) = DMTN
      STATEV(4) = DMCN
      STATEV(5) = OMTN
      STATEV(6) = OMCN
      STATEV(7) = FMSN
      STATEV(8) = DFT
      STATEV(9) = DFC
      STATEV(10) = DDMT
      STATEV(11) = DDMC
      STATEV(12) = DOMT
      STATEV(13) = DOMC
      STATEV(14) = DFMS
      STATEV(15) = DFTV
      STATEV(16) = DFCV
      STATEV(17) = DDMTV
      STATEV(18) = DDMCV
      STATEV(19) = DOMTV
      STATEV(20) = DOMCV
      STATEV(21) = DFMSV
C
C      TO COMPUTE THE ENERGY
C
      DO I = 1, NDI
         SSE = SSE + HALF * (STRESS(I) + OLD_STRESS(I)) * DSTRAN(I)
      END DO
      DO I = NDI+1, NTENS
         SSE = SSE + (STRESS(I) + OLD_STRESS(I)) * DSTRAN(I)
      END DO
C
C     TO COMPUTE THE INTERNAL ENERGY WITHOUT VISCOUS REGULARIZATION
C
      DO I = 1, NDI
         SCD = SCD + HALF * (STRESS(I) + OLD_STRESS(I)
     1        -D_STRESS(I)-DOLD_STRESS(I)) * DSTRAN(I)
      END DO
      DO I = NDI+1, NTENS
         SCD = SCD + (STRESS(I) + OLD_STRESS(I)
     1        -D_STRESS(I)-DOLD_STRESS(I)) * DSTRAN(I)
      END DO
      RETURN
      END
C***************************************************************************
C        CALCULATE THE STRESS BASED ON THE DAMAGE VARAIBLES ****************
C***************************************************************************
      SUBROUTINE GetStress(CFULL,CDFULL,DFFV,DDMFV,DOMFV,DFMSV,STRESS,
     1 STRANT,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION CFULL(6,6),CDFULL(6,6),STRESS(NTENS),
     1     STRANT(6)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO = 2.D0)
C
C
      DO I = 1, 6
         DO J = 1, 6
            CDFULL(I,J)=CFULL(I,J)
         END DO
      END DO
      IF ((DFFV .NE. ZERO) .OR. (DDMFV .NE. ZERO) .OR. (DOMFV .NE. ZERO)
     1   .OR. (DFMSV .NE. ZERO)) THEN
      CDFULL(1,1) = (ONE - DFFV) ** TWO * CFULL(1,1)
      CDFULL(1,2) = (ONE - DFFV) * (ONE - DDMFV) * CFULL(1,2)
      CDFULL(1,3) = (ONE - DFFV) * (ONE - DOMFV) * CFULL(1,3)
      CDFULL(2,2) = (ONE - DDMFV) ** TWO * CFULL(2,2)
      CDFULL(2,3) = (ONE - DDMFV) * (ONE - DOMFV) * CFULL(2,3)
      CDFULL(3,3) = (ONE - DOMFV) ** TWO * CFULL(3,3)
      CDFULL(4,4) = (ONE - DFMSV) ** TWO * CFULL(4,4)
      CDFULL(5,5) = (ONE - DFFV) * (ONE - DOMFV) * CFULL(5,5)
      CDFULL(6,6) = (ONE - DDMFV) * (ONE - DOMFV) * CFULL(6,6)
      END IF
C
C     UPDATE THE STRESS STATE
C
      DO I = 1,NTENS
	    STRESS(I)=ZERO
	      DO J = 1, NTENS
		   STRESS(I)=STRESS(I)+CDFULL(I,J) * STRANT(J)
		  END DO
      END DO
      RETURN
      END
C*********************************************************
C     TO CHECK THE FAILURE INITIATION*********************
C*********************************************************
      SUBROUTINE CheckFailureIni(XTFE,XCFE,YTFE,YCFE,ZTFE,ZCFE,XCFMO,
     1  SOTFE,SOTHFE,STTFE,STRANT,GF,GM, CELENT, CFULL,DFT, DFC, DFMS, 
     2 DDMT, DDMC, DOMT, DOMC, DDFT, DDFC, DDFMS, DDDMT, DDDMC, DDOMT,
     3 DDOMC, NTENS, DFTOLD, DFCOLD, DDMTOLD, DDMCOLD, DOMTOLD,
     4 DOMCOLD, DFMSOLD, DDFF, DDDMF, DDOMF,
     5 FTN,FCN,DMTN,DMCN,OMTN,OMCN,FMSN) 
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION DDFT(6),DDFC(6),DDFMS(6),DDDMT(6),DDDMC(6),DDOMT(6),
     1 DDOMC(6),STRANT(6),CFULL(6,6)
      DIMENSION DDFF(6),DDDMF(6),DDOMF(6)
      DIMENSION DFTDE(6),DFCDE(6),DFMSDE(6),DDMTDE(6),DDMCDE(6),
     1 DOMTDE(6),DOMCDE(6)
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)
C
C       fiber tensile (FT) failure
C
      IF (STRANT(1) .GT. ZERO) THEN
      FT = (STRANT(1) / XTFE)**TWO + (STRANT(4) / SOTFE)**TWO
     1 + (STRANT(5) / SOTHFE)**TWO
      END IF
C
      IF (FT .GT. ZERO) THEN
	   FTN = SQRT(FT)
      ELSE
	   FTN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE   (DFMT)
C
      DFT =ZERO
      DDDFTN = ZERO
      DO I = 1,6
          DFTDE(I) = ZERO
          DDFT(I) = ZERO
      END DO
      IF (FTN .GT. ONE) THEN
      CALL DamageEvaluation(CFULL(1,1), FTN, GT, CELENT, STRANT(1), DFT,
     1 DDDFTN)
          IF (DFT .GT. DFTOLD) THEN
              DFTDE(1) = (HALF*FTN)*TWO*STRANT(1) / (XTFE**TWO)
              DFTDE(4) = (HALF*FTN)*TWO*STRANT(4) / (SOTFE**TWO)
              DFTDE(5) = (HALF*FTN)*TWO*STRANT(5) / (SOTHFE**TWO)
              DO I = 1, 6
               DDFT(I) = DFTDE(I) * DDDFTN
              END DO
          END IF
      END IF
      DFT = MAX(DFT, DFTOLD)
C
C      fiber compressive (FC) failure AND fiber-matrix shear-out (FMS)
C
      IF (STRANT(1) .LT. ZERO) THEN
      FC = (STRANT(1) / XCFE)**TWO
      FMS = (STRANT(1) / XCFMO)**TWO + (STRANT(4) / SOTFE)**TWO
     1 + (STRANT(5) / SOTHFE)**TWO
      END IF
C
      IF (FC .GT. ZERO) THEN
	   FCN = SQRT(FC)
      ELSE
	   FCN = ZERO
      END IF
C
      IF (FMS .GT. ZERO) THEN
	   FMSN = SQRT(FMS)
      ELSE
	   FMSN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE   (DFMC)
C
      DFC =ZERO
      DDDFCN = ZERO
      DO I = 1,6
          DFCDE(I) = ZERO
          DDFC(I) = ZERO
      END DO
      IF (FCN .GT. ONE) THEN
          CALL DamageEvaluation(CFULL(1,1), FCN, GT, CELENT, STRANT(1),
     1      DFC, DDDFCN)
	      IF (DFC .GT. DFCOLD) THEN
              DFCDE(1) = (HALF*FCN)*TWO*STRANT(1) / (XCFE**TWO)
              DO I = 1, 6
               DDFC(I) = DFCDE(I) * DDDFCN
              END DO
            END IF
      END IF
      DFC = MAX(DFC, DFCOLD)
C
C     INITIALIZE THE ARRAY AND VARIABLE  (DFMS)
C
      DFMS =ZERO
      DDDFMSN = ZERO
      DO I = 1,6
          DFMSDE(I) = ZERO
          DDFMS(I) = ZERO
      END DO
	   IF (FMSN .GT. ONE) THEN
	     CALL DamageEvaluation(CFULL(1,1), FMSN, GM, CELENT, STRANT(1),
     1      DFMS, DDDFMSN)
	      IF (DFMS .GT. DFMSOLD) THEN
              DFMSDE(1) = (HALF*FMSN)*TWO*STRANT(1) / (XCFMO**TWO)
              DFMSDE(4) = (HALF*FMSN)*TWO*STRANT(4) / (SOTFE**TWO)
              DFMSDE(5) = (HALF*FMSN)*TWO*STRANT(5) / (SOTHFE**TWO)
              DO I = 1, 6
              DDFMS(I) = DFMSDE(I) * DDDFMSN
              END DO
            END IF
      END IF
      DFMS = MAX(DFMS, DFMSOLD)
C
C      in-plane matrix cracking (DMT)
C
      IF (STRANT(2) .GT. ZERO) THEN
      DMT = (STRANT(2) / YTFE)**TWO + (STRANT(4) / SOTFE)**TWO
     1 + (STRANT(6) / STTFE)**TWO
      END IF
C
      IF (DMT .GT. ZERO) THEN
      DMTN = SQRT(DMT)
      ELSE
      DMTN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE  (DDMT)
C
      DDMT =ZERO
      DDDDMTN = ZERO
      DO I = 1,6
          DDMTDE(I) = ZERO
          DDDMT(I) = ZERO
      END DO
      IF (DMTN .GT. ONE) THEN
          CALL DamageEvaluation(CFULL(2,2), DMTN, GM, CELENT, STRANT(2),
     1     DDMT, DDDDMTN)
	      IF (DDMT .GT. DDMTOLD) THEN
              DDMTDE(2) = (HALF*DMTN)*TWO*STRANT(2) / (YTFE**TWO)
              DDMTDE(4) = (HALF*DMTN)*TWO*STRANT(4) / (SOTFE**TWO)
              DDMTDE(6) = (HALF*DMTN)*TWO*STRANT(6) / (STTFE**TWO)
              DO I = 1, 6
              DDDMT(I) = DDMTDE(I) * DDDDMTN
              END DO
           END IF
      END IF
      DDMT = MAX(DDMT, DDMTOLD)
C
C      in-plane matrix crushing (DMC)
C
      IF (STRANT(2) .LT. ZERO) THEN
      DMC = (STRANT(2) / YCFE)**TWO + (STRANT(4) / SOTFE)**TWO
     1 + (STRANT(6) / STTFE)**TWO
      END IF
C
      IF (DMC .GT. ZERO) THEN
	   DMCN = SQRT(DMC)
      ELSE
	   DMCN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE  (DDMC)
C
      DDMC =ZERO
      DDDDMCN = ZERO
      DO I = 1,6
          DDMCDE(I) = ZERO
          DDDMC(I) = ZERO
      END DO
      IF (DMCN .GT. ONE) THEN
          CALL DamageEvaluation(CFULL(2,2), DMCN, GM, CELENT, STRANT(2),
     1     DDMC, DDDDMCN)
	      IF (DDMC .GT. DDMCOLD) THEN
              DDMCDE(2) = (HALF*DMCN)*TWO*STRANT(2) / (YCFE**TWO)
              DDMCDE(4) = (HALF*DMCN)*TWO*STRANT(4) / (SOTFE**TWO)
              DDMCDE(6) = (HALF*DMCN)*TWO*STRANT(6) / (STTFE**TWO)
              DO I = 1, 6
              DDDMC(I) = DDMCDE(I) * DDDDMCN
              END DO
            END IF
      END IF
      DDMC = MAX(DDMC, DDMCOLD)
C
C      out-of-plane matrix cracking (OMT)
C
      IF (STRANT(3) .GT. ZERO) THEN
      OMT = (STRANT(3) / ZTFE)**TWO + (STRANT(5) / SOTHFE)**TWO
     1 + (STRANT(6) / STTFE)**TWO
      END IF
C
      IF (OMT .GT. ZERO) THEN
	   OMTN = SQRT(OMT)
      ELSE
	   OMTN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE  (DOMT)
C
      DOMT =ZERO
      DDDOMTN = ZERO
      DO I = 1,6
          DOMTDE(I) = ZERO
          DDOMT(I) = ZERO
      END DO
	   IF (OMTN .GT. ONE) THEN
	     CALL DamageEvaluation(CFULL(3,3), OMTN, GM, CELENT, STRANT(3),
     1      DOMT, DDDOMTN)
	      IF (DOMT .GT. DOMTOLD) THEN
              DOMTDE(3) = (HALF*OMTN)*TWO*STRANT(3) / (ZTFE**TWO)
              DOMTDE(5) = (HALF*OMTN)*TWO*STRANT(5) / (SOTHFE**TWO)
              DOMTDE(6) = (HALF*OMTN)*TWO*STRANT(6) / (STTFE**TWO)
              DO I = 1, 6
              DDOMT(I) = DOMTDE(I) * DDDOMTN
              END DO
            END IF
      END IF
      DOMT = MAX(DOMT, DOMTOLD)
C
C     out-of-plane matrix crushing (OMC)
C
      IF (STRANT(3) .LT. ZERO) THEN
      OMC = (STRANT(3) / ZCFE)**TWO + (STRANT(5) / SOTHFE)**TWO
     1 + (STRANT(6) / STTFE)**TWO
      END IF
C
      IF (OMC .GT. ZERO) THEN
	   OMCN = SQRT(OMC)
      ELSE
	   OMCN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE  (DOMC)
C
      DOMC =ZERO
      DDDOMCN = ZERO
      DO I = 1,6
          DOMCDE(I) = ZERO
          DDOMC(I) = ZERO
      END DO
      IF (OMCN .GT. ONE) THEN
	     CALL DamageEvaluation(CFULL(3,3), OMCN, GM, CELENT, STRANT(3),
     1      DOMC, DDDOMCN)
	      IF (DOMC .GT. DOMCOLD) THEN
              DOMCDE(3) = (HALF*OMCN)*TWO*STRANT(3) / (ZCFE**TWO)
              DOMCDE(5) = (HALF*OMCN)*TWO*STRANT(5) / (SOTHFE**TWO)
              DOMCDE(6) = (HALF*OMCN)*TWO*STRANT(6) / (STTFE**TWO)
              DO I = 1, 6
              DDOMC(I) = DOMCDE(I) * DDDOMCN
              END DO
            END IF
      END IF
      DOMC = MAX (DOMC, DOMCOLD)   
      RETURN
      END
C********************************************************************
C     SUBROUTINE TO EVALUATE THE DAMAGE AND THE DERIVATIVE***********
C********************************************************************
      SUBROUTINE DamageEvaluation(STIFF, FN, G, CELENT, STRAIN, D,
     1 DDDFN)
C     CALCULATE DAMAGE VARIABLE
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ONE = 1.D0,TWO = 2.D0, ZERO = 0.D0)
C
C       DAMAGE VARIABLES
C
       TERM1 = STIFF * (STRAIN ** TWO) * (ONE - FN) * CELENT / G
       D = ONE - EXP(TERM1) / FN
C
       DDDFN = (ONE / FN + TERM1) * (ONE - D)
      RETURN
      END
C**********************************************************************
C     SUBROUTINE TO GET THE DERIVATIVE MATRIX
C********************************************************************
      SUBROUTINE ElasticDerivative(CFULL,DFFV,DDMFV,DOMFV,DFMSV, DCDFF,
     1 DCDDMF,DCDOMF,DCDFMS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION CFULL(6,6), DCDFF(6,6),DCDDMF(6,6),
     1     DCDOMF(6,6), DCDFMS(6,6)
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)
C     initialize the data to zero
      DO I = 1, 6
         DO J = 1, 6
            DCDFF(I,J) = ZERO
            DCDFMS(I,J) = ZERO
            DCDDMF(I,J) = ZERO
            DCDOMF(I,J) = ZERO
         END DO
      END DO
C
C      DCDFF
C
      DCDFF(1,1) = -(ONE-DFFV)*TWO*CFULL(1,1)
      DCDFF(1,2) = -(ONE-DDMFV)*CFULL(1,3)
      DCDFF(1,3) = -(ONE-DOMFV)*CFULL(1,3)
      DCDFF(5,5) = -(ONE-DOMFV)*CFULL(5,5)
C
C      DCDDMF
C
      DCDDMF(1,2) = -(ONE-DFFV)*CFULL(1,2)
      DCDDMF(2,2) = -(ONE-DDMFV)*TWO*CFULL(2,2)
      DCDDMF(2,3) = -(ONE-DOMFV)*CFULL(2,3)
      DCDDMF(6,6) = -(ONE-DOMFV)*CFULL(6,6)
C
C      DCDOMF
C
      DCDOMF(1,3) = -(ONE-DFFV)*CFULL(1,3)
      DCDOMF(2,3) = -(ONE-DDMFV)*CFULL(2,3)
      DCDOMF(3,3) = -(ONE-DOMFV)*TWO*CFULL(3,3)
      DCDOMF(5,5) = -(ONE-DFFV)*CFULL(5,5)
      DCDOMF(6,6) = -(ONE-DDMFV)*CFULL(6,6)
C
C      DCDFMS
C
      DCDFMS(4,4) = -(ONE-DFMSV)*TWO*CFULL(4,4)
      RETURN
      END