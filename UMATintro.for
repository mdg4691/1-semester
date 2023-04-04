      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
	 

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
      1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
      2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
      3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
C     4 JSTEP(4),UPSTRAN(NTENS),dmg(6),e(6)
      PARAMETER (ZERO = 0.D0,ONE = 1.D0,TWO = 2.D0, HALF = 0.5D0)
      E11=PROPS(1)             ! Longitudinal elastic Young's modulus(Fiber: 11 axis)  
      E22=PROPS(2)             ! Transverse elastic Young's modulus   (Matrix: 22 axis)  
      E33=PROPS(3)             ! Interlaminar elastic Young's modulus (Interlaminar: 33 axis)
      G12=PROPS(4)             ! Elastic Shear Modulus  (1-2 plane)
      G13=PROPS(5)             ! Elastic Shear Modulus  (1-3 plane)   
      G23=PROPS(6)             ! Elastic Shear Modulus  (2-3 plane) 
      NU12=PROPS(7)            ! Poisson ratio   (1-2 plane)
      NU13=PROPS(8)            ! Poisson ratio   (1-3 plane)  
      NU23=PROPS(9)            ! Poisson ratio   (2-3 plane)                   
      1TS=PROPS(10)            ! Laminate tensile strength (11 fiber direction)   
      2TS=PROPS(11)            ! Laminate tensile strength  (22 matrix direction)
      3TS=PROPS(12)            ! Laminate interlaminar tensile strength (33 direction)
      1CS=PROPS(13)            ! Laminate compressive strength  (11 fiber direction)
      2CS=PROPS(14)            ! Laminate compressive strength (22 matrix direction)
      3CS=PROPS(15)            ! Laminate interlaminar compressive strength (33 direction)
      S12=PROPS(16)            ! Laminate Shear Strength  (1-2 plane) 
      S13=PROPS(17)            ! Laminate Shear Strength  (1-3 plane) 
      S23=PROPS(18)            ! Laminate Shear Strength  (2-3 plane)

C     SM IS STIFFNESS MATRIX
      DIMENSION SM(6,6)
      DO I = 1, 6
         DO J = 1, 6
            CFULL(I,J)=ZERO
         END DO
      END DO
      NU21 = NU12 / E11 * E22  ! Poisson ratio 2-1
      ATEMP = ONE - TWO * NU12 * NU21 - NU23 ** TWO
     1     - TWO * NU12 * NU21 * NU23
      SM(1,1) = E11 * (ONE - NU23 ** TWO) / ATEMP
      SM(2,2) = E22 * (ONE - NU12 * NU21) / ATEMP
      SM(3,3) = SM(2,2)
      SM(1,2) = E22 * (NU12 + NU12 * NU23) / ATEMP
      SM(1,3) = SM(1,2)
      SM(2,3) = TENT * (NU23 + NU12 * NU21) / ATEMP    
      SM(4,4) = G12
      SM(5,5) = G12
      SM(6,6) = G23
      DO I = 2, 6
         DO J = 1, I-1
            SM(I,J) = SM(J,I)
         END DO
      END DO
C     CALCULATE THE FAILURE STRAIN BY FAILUER STRESS
      1TSN = 1TS / SM(1,1)      ! FAILURE STRAIN 1 DIRECTION IN TENSION
      1CSN = 1CS / SM(1,1)      ! FAILURE STRAIN 1 DIRECTION IN COMPRESSION
      2TSN = 2TS / SM(2,2)      ! TENSILE FAILURE STRAIN 2 DIRECTION
      2CSN = 2CS / SM(2,2)      ! COMPRESSIVE FAILURE STRAIN 2 DIRECTION
      FSSN12 = S12 / G12        ! FAILURE SHEAR STRAIN 1-2     
C     +3D     
      3TSN = 3TS / SM(3,3)      ! FAILURE STRAIN 3 DIRECTION IN TENSION
      3CSN = 3SC / SM(3,3)      ! FAILURE STRAIN 3 DIRECTION IN COMPRESSION
      FSSN13 = S13 / G13        ! FAILURE SHEAR STRAIN 1-3 
      FSSN23 = S23 / G23        ! FAILURE SHEAR STRAIN 2-3 

      