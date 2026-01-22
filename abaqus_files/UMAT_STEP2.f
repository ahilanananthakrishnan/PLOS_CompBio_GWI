C***********************************************************************
C  UMAT_STEP2.f  (STEP 2: Elastic preload, hypertrophy frozen)
C
C***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,
     &                  NOEL,NPT,LAYER,KSPT)

      INCLUDE 'ABA_PARAM.INC'
      REAL*8 STATEV(NSTATV)

      RETURN
      END

C***********************************************************************

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     & RPL,DDSDDT,DRPLDE,DRPLDT,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      INTEGER I,J,K,L
      REAL*8 STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS)
      REAL*8 TIME(2),DTIME
      REAL*8 PROPS(NPROPS)
      REAL*8 DFGRD1(3,3)

      REAL*8 Iden(3,3),Fe(3,3),be(3,3),sigma(3,3)
      REAL*8 Cmod(3,3,3,3),JAC(3,3,3,3),ANISO_C(3,3,3,3)

      REAL*8 mu,kappa,kappa_f
      REAL*8 Je,logJe,I4e
      REAL*8 Tdir(3),tdir(3),norm

C --- passive parameters
      mu      = PROPS(25)
      kappa   = PROPS(26)
      kappa_f = PROPS(27)

C --- tangential direction in reference configuration
      Tdir(1) = PROPS(33)
      Tdir(2) = PROPS(34)
      Tdir(3) = PROPS(35)

C --- identity
      DO I=1,3
         DO J=1,3
            Iden(I,J)=0.D0
         END DO
         Iden(I,I)=1.D0
      END DO

C --- hypertrophy frozen: Fe = F
      DO I=1,3
         DO J=1,3
            Fe(I,J)=DFGRD1(I,J)
         END DO
      END DO

      CALL MDET(Fe,Je)
      IF (Je .LE. 0.D0) THEN
         WRITE(6,*) 'ERROR: Je <= 0. NOEL=',NOEL,' NPT=',NPT,
     &              ' TIME=',TIME(2),' KINC=',KINC,' Je=',Je
         STOP
      ENDIF
      logJe = DLOG(Je)

C --- be = Fe * Fe^T
      DO I=1,3
         DO J=1,3
            be(I,J)=0.D0
            DO K=1,3
               be(I,J)=be(I,J)+Fe(I,K)*Fe(J,K)
            END DO
         END DO
      END DO

C --- push-forward tangential direction: tdir = Fe*Tdir / ||Fe*Tdir||
      CALL MV3MULT(Fe,Tdir,tdir)
      norm = DSQRT(tdir(1)**2 + tdir(2)**2 + tdir(3)**2)

      IF (norm .GT. 1.D-12) THEN
         I4e = norm*norm
         DO I=1,3
            tdir(I)=tdir(I)/norm
         END DO
      ELSE
         I4e = 1.D0
         DO I=1,3
            tdir(I)=Tdir(I)
         END DO
      ENDIF

C --- passive Cauchy stress (iso + aniso)
      DO I=1,3
         DO J=1,3
            sigma(I,J)=
     &        ((kappa*logJe-mu)*Iden(I,J)+mu*be(I,J))/Je
     &        + (3.D0*kappa_f/Je)*(I4e-1.D0)**2
     &          * tdir(I)*tdir(J)
         END DO
      END DO

      STRESS(1)=sigma(1,1)
      STRESS(2)=sigma(2,2)
      STRESS(3)=sigma(3,3)
      STRESS(4)=sigma(1,2)
      STRESS(5)=sigma(1,3)
      STRESS(6)=sigma(2,3)

C --- tangent: iso + aniso + geometric
      DO I=1,3
         DO J=1,3
            DO K=1,3
               DO L=1,3

                  ANISO_C(I,J,K,L) =
     &              (12.D0*kappa_f/Je)*(I4e-1.D0)
     &              * tdir(I)*tdir(J)
     &              * tdir(K)*tdir(L)

                  Cmod(I,J,K,L)=
     &             (kappa/Je)*Iden(I,J)*Iden(K,L)
     &           - ((kappa*logJe-mu)/Je)
     &             * (Iden(I,K)*Iden(J,L) + Iden(I,L)*Iden(J,K))

                  JAC(I,J,K,L)=
     &             Cmod(I,J,K,L) + ANISO_C(I,J,K,L)
     &           + 0.5D0*(Iden(I,K)*sigma(J,L)
     &                  + Iden(I,L)*sigma(J,K)
     &                  + Iden(J,K)*sigma(I,L)
     &                  + Iden(J,L)*sigma(I,K))

               END DO
            END DO
         END DO
      END DO

      CALL JAC3D(JAC,DDSDDE)

      SSE=0.D0
      SPD=0.D0
      SCD=0.D0
      RPL=0.D0
      DRPLDT=0.D0

      RETURN
      END

C***********************************************************************
C  Helper routines
C***********************************************************************

      SUBROUTINE JAC3D(SPTAN,DDSDDE)
      REAL*8 SPTAN(3,3,3,3),DDSDDE(6,6)

      DDSDDE(1,1)=SPTAN(1,1,1,1)
      DDSDDE(1,2)=SPTAN(1,1,2,2)
      DDSDDE(1,3)=SPTAN(1,1,3,3)
      DDSDDE(1,4)=SPTAN(1,1,1,2)
      DDSDDE(1,5)=SPTAN(1,1,1,3)
      DDSDDE(1,6)=SPTAN(1,1,2,3)

      DDSDDE(2,1)=SPTAN(2,2,1,1)
      DDSDDE(2,2)=SPTAN(2,2,2,2)
      DDSDDE(2,3)=SPTAN(2,2,3,3)
      DDSDDE(2,4)=SPTAN(2,2,1,2)
      DDSDDE(2,5)=SPTAN(2,2,1,3)
      DDSDDE(2,6)=SPTAN(2,2,2,3)

      DDSDDE(3,1)=SPTAN(3,3,1,1)
      DDSDDE(3,2)=SPTAN(3,3,2,2)
      DDSDDE(3,3)=SPTAN(3,3,3,3)
      DDSDDE(3,4)=SPTAN(3,3,1,2)
      DDSDDE(3,5)=SPTAN(3,3,1,3)
      DDSDDE(3,6)=SPTAN(3,3,2,3)

      DDSDDE(4,1)=SPTAN(1,2,1,1)
      DDSDDE(4,2)=SPTAN(1,2,2,2)
      DDSDDE(4,3)=SPTAN(1,2,3,3)
      DDSDDE(4,4)=SPTAN(1,2,1,2)
      DDSDDE(4,5)=SPTAN(1,2,1,3)
      DDSDDE(4,6)=SPTAN(1,2,2,3)

      DDSDDE(5,1)=SPTAN(1,3,1,1)
      DDSDDE(5,2)=SPTAN(1,3,2,2)
      DDSDDE(5,3)=SPTAN(1,3,3,3)
      DDSDDE(5,4)=SPTAN(1,3,1,2)
      DDSDDE(5,5)=SPTAN(1,3,1,3)
      DDSDDE(5,6)=SPTAN(1,3,2,3)

      DDSDDE(6,1)=SPTAN(2,3,1,1)
      DDSDDE(6,2)=SPTAN(2,3,2,2)
      DDSDDE(6,3)=SPTAN(2,3,3,3)
      DDSDDE(6,4)=SPTAN(2,3,1,2)
      DDSDDE(6,5)=SPTAN(2,3,1,3)
      DDSDDE(6,6)=SPTAN(2,3,2,3)

      RETURN
      END

C***********************************************************************

      SUBROUTINE MDET(A,DET)
      REAL*8 A(3,3),DET
      DET = A(1,1)*A(2,2)*A(3,3)
     &    + A(1,2)*A(2,3)*A(3,1)
     &    + A(1,3)*A(2,1)*A(3,2)
     &    - A(1,3)*A(2,2)*A(3,1)
     &    - A(1,2)*A(2,1)*A(3,3)
     &    - A(1,1)*A(2,3)*A(3,2)
      RETURN
      END

C***********************************************************************

      SUBROUTINE MV3MULT(A,V,W)
      REAL*8 A(3,3),V(3),W(3)
      INTEGER I,K
      DO I=1,3
         W(I)=0.D0
         DO K=1,3
            W(I)=W(I)+A(I,K)*V(K)
         END DO
      END DO
      RETURN
      END
