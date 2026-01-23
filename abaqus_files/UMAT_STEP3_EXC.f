C***********************************************************************
C  UMAT_STEP3.f  
C
C***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,
     &                  NOEL,NPT,LAYER,KSPT)

      INCLUDE 'ABA_PARAM.INC'
      REAL*8 STATEV(NSTATV)

      STATEV(1) = 0.D0
      STATEV(2) = 0.D0
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

      REAL*8 Iden(3,3),Fe(3,3),be(3,3)
      REAL*8 sigma_pas(3,3),sigma_act(3,3),sigma_tot(3,3)
      REAL*8 Cmod(3,3,3,3),ANISO_C(3,3,3,3),JAC(3,3,3,3)

      REAL*8 mu,kappa,kappa_f
      REAL*8 day,t_pb,tc
      REAL*8 m_macro,y2_macro
      REAL*8 gamma_ch
      REAL*8 r_ch,k_ch,ch_ss
      REAL*8 k_ex,alpha_act_ex
      REAL*8 Tdir(3),t_dir(3),norm
      REAL*8 Je,logJe,I4e

      REAL*8 teval,step_t,ramp
      REAL*8 ch,chPB
      REAL*8 intM_norm,I1,I2,t_diff
      REAL*8 chi_ch,eta_ex,phi_ex,sigma_act_ex

C --- props
      day          = PROPS(1)
      t_pb         = PROPS(2)
      tc           = PROPS(3)

      m_macro      = PROPS(12)
      y2_macro     = PROPS(14)

      gamma_ch     = PROPS(15)

      r_ch         = PROPS(17)
      k_ch         = PROPS(18)
      ch_ss        = PROPS(19)

      k_ex         = PROPS(23)

      mu           = PROPS(25)
      kappa        = PROPS(26)
      kappa_f      = PROPS(27)

      alpha_act_ex = PROPS(28)

      Tdir(1)      = PROPS(33)
      Tdir(2)      = PROPS(34)
      Tdir(3)      = PROPS(35)

C --- evaluation time for biology: teval = day
      teval  = day

C --- TIME(2) is the current step time; DTIME advances to the end of increment
C --- Stept = elapsed time controlling LINEAR activation ramp
      step_t = TIME(2) + DTIME

C --- identity
      DO I=1,3
         DO J=1,3
            Iden(I,J)=0.D0
         END DO
         Iden(I,I)=1.D0
      END DO

C ----------------------------------------------------------------
C because Fh=I, F=Fe
C ----------------------------------------------------------------
      DO I=1,3
         DO J=1,3
            Fe(I,J)=DFGRD1(I,J)
         END DO
      END DO

      CALL MDET(Fe,Je)
      IF (Je .LE. 0.D0) THEN
         WRITE(6,*) 'ERROR: Je <= 0 in Step 3. NOEL=',NOEL,' NPT=',NPT,
     &              ' KINC=',KINC,' Je=',Je
         STOP
      ENDIF
      logJe = DLOG(Je)

C --- be = Fe Fe^T
      DO I=1,3
         DO J=1,3
            be(I,J)=0.D0
            DO K=1,3
               be(I,J)=be(I,J)+Fe(I,K)*Fe(J,K)
            END DO
         END DO
      END DO

C --- push-forward tangential direction: t_dir = Fe*Tdir / ||Fe*Tdir||
      CALL MV3MULT(Fe, Tdir, t_dir)
      norm = DSQRT(t_dir(1)**2 + t_dir(2)**2 + t_dir(3)**2)

      IF (norm .GT. 1.D-12) THEN
         I4e = norm*norm
         DO I = 1, 3
            t_dir(I) = t_dir(I)/norm
         END DO
      ELSE
         I4e = 1.D0
         DO I = 1, 3
            t_dir(I) = Tdir(I)
         END DO
      ENDIF

C ----------------------------------------------------------------
C passive Cauchy stress
C ----------------------------------------------------------------
      DO I=1,3
         DO J=1,3
            sigma_pas(I,J)=
     &        ((kappa*logJe-mu)*Iden(I,J)+mu*be(I,J))/Je
     &        + (3.D0*kappa_f/Je)*(I4e-1.D0)**2
     &          * t_dir(I)*t_dir(J)
         END DO
      END DO

C ----------------------------------------------------------------
C excitatory active stress
C ----------------------------------------------------------------

      IF (t_pb .LE. 1.D-12) THEN
         ch = 1.D0
      ELSEIF (teval .LT. t_pb) THEN
         ch = 1.D0 - r_ch*teval
      ELSE
         chPB = 1.D0 - r_ch*t_pb
         ch   = ch_ss - (ch_ss - chPB)*DEXP(-k_ch*(teval - t_pb))
      ENDIF
      IF (ch .LT. 0.D0) ch = 0.D0

      IF (teval .LE. 0.D0) THEN
         intM_norm = 0.D0
      ELSEIF (teval .LE. t_pb) THEN
         intM_norm = teval + 0.5D0*m_macro*teval*teval
      ELSE
         I1 = t_pb + 0.5D0*m_macro*t_pb*t_pb
         IF (y2_macro .GT. 1.D-12) THEN
            t_diff = teval - t_pb
            I2 = (1.D0 + m_macro*t_pb)/y2_macro
     &           * (1.D0 - DEXP(-y2_macro*t_diff))
         ELSE
            I2 = (1.D0 + m_macro*t_pb)*(teval - t_pb)
         ENDIF
         intM_norm = I1 + I2
      ENDIF

C chi_ch = exp( - (gamma_ch/tc) * intM_norm )
      IF (tc .GT. 1.D-12) THEN
         chi_ch = DEXP(-(gamma_ch/tc)*intM_norm)
      ELSE
         chi_ch = 1.D0
      ENDIF

C eta_ex = chi_ch * ch
      eta_ex = chi_ch*ch

C phi_ex = 1 / (1 + exp(k_ex * teval))
      phi_ex = 1.D0/(1.D0 + DEXP(k_ex*teval))

C scalar active stress (excitatory)
      sigma_act_ex = alpha_act_ex*eta_ex*phi_ex

C ramp in Step 3
      IF (step_t .EQ. 0.D0) THEN
         ramp = 0.D0
      ELSEIF (step_t .LT. 1.D0) THEN
         ramp = step_t
      ELSE
         ramp = 1.D0
      ENDIF

      DO I=1,3
         DO J=1,3
            sigma_act(I,J) = sigma_act_ex*ramp*t_dir(I)*t_dir(J)
            sigma_tot(I,J) = sigma_pas(I,J) + sigma_act(I,J)
         END DO
      END DO

C ----------------------------------------------------------------
C tangent: iso + aniso + geometric
C ----------------------------------------------------------------
      DO I=1,3
         DO J=1,3
            DO K=1,3
               DO L=1,3

                  ANISO_C(I,J,K,L) =
     &              (12.D0*kappa_f/Je)*(I4e-1.D0)
     &              * t_dir(I)*t_dir(J)*t_dir(K)*t_dir(L)

                  Cmod(I,J,K,L)=
     &             (kappa/Je)*Iden(I,J)*Iden(K,L)
     &           - ((kappa*logJe-mu)/Je)
     &             * (Iden(I,K)*Iden(J,L) + Iden(I,L)*Iden(J,K))

                  JAC(I,J,K,L)=
     &             Cmod(I,J,K,L) + ANISO_C(I,J,K,L)
     &           + 0.5D0*( Iden(I,K)*sigma_tot(J,L)
     &                    +Iden(I,L)*sigma_tot(J,K)
     &                    +Iden(J,K)*sigma_tot(I,L)
     &                    +Iden(J,L)*sigma_tot(I,K) )

               END DO
            END DO
         END DO
      END DO

      CALL JAC3D(JAC,DDSDDE)

C --- push stress to Voigt
      STRESS(1)=sigma_tot(1,1)
      STRESS(2)=sigma_tot(2,2)
      STRESS(3)=sigma_tot(3,3)
      STRESS(4)=sigma_tot(1,2)
      STRESS(5)=sigma_tot(1,3)
      STRESS(6)=sigma_tot(2,3)

      SSE=0.D0
      SPD=0.D0
      SCD=0.D0
      RPL=0.D0
      DRPLDT=0.D0

      RETURN
      END

C***********************************************************************
C Helper routines
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
C***********************************************************************
