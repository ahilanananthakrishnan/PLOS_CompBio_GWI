C***********************************************************************
C  UMAT_STEP1.f
C
C  STEP 1: Inflammation-driven hypertrophy (growth) 
C***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,
     &                  NOEL,NPT,LAYER,KSPT)

      INCLUDE 'ABA_PARAM.INC'

      REAL*8 STATEV(NSTATV)

      STATEV(1) = 1.D0
      STATEV(2) = 1.D0
      STATEV(3) = 0.D0

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

      REAL*8 Iden(3,3),F(3,3),Fe(3,3),Fh(3,3),FhInv(3,3)
      REAL*8 be(3,3),sigma(3,3)
      REAL*8 Cmod(3,3,3,3),JAC(3,3,3,3),ANISO_C(3,3,3,3)

      REAL*8 mu,kappa,kappa_f
      REAL*8 a,d1,d2
      REAL*8 theta_max0,b_h,c_h,K_h,n_h
      REAL*8 vartheta,Ival,intI
      REAL*8 Je,logJe,I4e
      REAL*8 Rdir(3),Tdir(3),tdir(3),norm
      REAL*8 Iold


      REAL*8 t_pb

C --- props
      t_pb       = PROPS(2)

      a          = PROPS(4)
      d1         = PROPS(5)
      d2         = PROPS(6)

      theta_max0 = PROPS(7)
      b_h        = PROPS(8)
      c_h        = PROPS(9)
      K_h        = PROPS(10)
      n_h        = PROPS(11)

      mu         = PROPS(25)
      kappa      = PROPS(26)
      kappa_f    = PROPS(27)

      Rdir(1)    = PROPS(30)
      Rdir(2)    = PROPS(31)
      Rdir(3)    = PROPS(32)

      Tdir(1)    = PROPS(33)
      Tdir(2)    = PROPS(34)
      Tdir(3)    = PROPS(35)

C --- identity
      DO I=1,3
         DO J=1,3
            Iden(I,J)=0.D0
         END DO
         Iden(I,I)=1.D0
      END DO

C --- fetch state
      vartheta = STATEV(1)
      Ival     = STATEV(2)
      intI     = STATEV(3)

C --- inflammation ODE + integral during exposure only
      IF (TIME(2) .LT. t_pb) THEN
         Iold = Ival
         Ival = Ival + DTIME*(a - d1*Iold)
         intI = intI + 0.5D0*(Iold + Ival)*DTIME
      ELSE
         Ival = Ival - DTIME*d2*Ival
      ENDIF

      STATEV(2)=Ival
      STATEV(3)=intI

C --- hypertrophy evolution
      vartheta = vartheta
     & + DTIME*(c_h*Ival**n_h/(K_h**n_h + Ival**n_h))
     & * (1.D0 - vartheta/(theta_max0 + b_h*intI))

      STATEV(1)=vartheta

C --- growth tensor Fh
      DO I=1,3
         DO J=1,3
            Fh(I,J)=Iden(I,J) + (vartheta-1.D0)*Rdir(I)*Rdir(J)
         END DO
      END DO

C --- deformation gradient F
      DO I=1,3
         DO J=1,3
            F(I,J)=DFGRD1(I,J)
         END DO
      END DO

C --- Fe = F * inv(Fh)
      CALL M3INV(Fh,FhInv)
      CALL M3MULT(F,FhInv,Fe)

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


C --- stress (passive iso + passive aniso)
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

      SUBROUTINE M3INV(A,AINV)
      REAL*8 A(3,3),AINV(3,3)
      REAL*8 DET,COF(3,3),ADJ(3,3)
      INTEGER I,J

      CALL MDET(A,DET)
      IF (DABS(DET) .LT. 1.D-20) THEN
         WRITE(6,*) 'ERROR: singular matrix in M3INV, det=',DET
         STOP
      ENDIF

      CALL MCOFAC(A,COF)
      CALL MTRANS(COF,ADJ)

      DO I=1,3
         DO J=1,3
            AINV(I,J)=ADJ(I,J)/DET
         END DO
      END DO
      RETURN
      END

C***********************************************************************

      SUBROUTINE MCOFAC(A,C)
      REAL*8 A(3,3),C(3,3)
      C(1,1)= A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C(1,2)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      C(1,3)= A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C(2,1)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      C(2,2)= A(1,1)*A(3,3)-A(1,3)*A(3,1)
      C(2,3)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      C(3,1)= A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C(3,2)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      C(3,3)= A(1,1)*A(2,2)-A(1,2)*A(2,1)
      RETURN
      END

C***********************************************************************

      SUBROUTINE MTRANS(A,AT)
      REAL*8 A(3,3),AT(3,3)
      INTEGER I,J
      DO I=1,3
         DO J=1,3
            AT(J,I)=A(I,J)
         END DO
      END DO
      RETURN
      END

C***********************************************************************

      SUBROUTINE M3MULT(A,B,C)
      REAL*8 A(3,3),B(3,3),C(3,3)
      INTEGER I,J,K
      DO I=1,3
         DO J=1,3
            C(I,J)=0.D0
            DO K=1,3
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
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
