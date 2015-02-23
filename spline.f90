!------------------------------------------------------------
      SUBROUTINE SPCOEF(N,X,F,B,C,D,FLAG)

      IMPLICIT NONE
      integer,intent(IN)::N
      integer,intent(OUT)::FLAG
      DOUBLE PRECISION,intent(IN):: X(N),F(N)
      double precision,intent(OUT):: B(N),C(N),D(N)
!
!  Input parameters:
!    N    = number of data points.
!    X    = vector of values of the independent variable ordered
!           SO THAT  x(i) < x(i+1)  for all I.
!    F    = vector of values of the dependent variable.
!  Output parameters:
!    B    = vector of S'(X(I)) values.
!    C    = vector of S"(X(I))/2 values.
!    D    = vector of S'''(X(I)+)/6 values (I .LT. N).
!    FLAG =  0  normal return;
!         = -1  input N .LE. 1;
!         = -2  X vector is incorrectly ordered.
!
!  The vectors X, F, B, C, D must be dimensioned at least N in the
!  calling program.
!
!  Local variables:
      INTEGER I,K
      DOUBLE PRECISION FP1,FPN,H,P,THREE,TWO,ZERO
      DATA THREE,TWO,ZERO/3.D0,2.D0,0.D0/
!
      IF (N .LE. 1) THEN
         FLAG = -1
         RETURN
      ENDIF
!
!     Calculate coefficients for the tri-diagonal system: store
!     sub-diagonal in B, diagonal in D, difference quotient in C.
!
      B(1) = X(2)-X(1)
      IF (B(1) .LE. ZERO) THEN
         FLAG = -2
         RETURN
      ENDIF
      C(1) = (F(2)-F(1))/B(1)
      IF (N .EQ. 2) THEN
         B(1) = C(1)
         C(1) = ZERO
         D(1) = ZERO
         B(2) = B(1)
         C(2) = ZERO
         FLAG = 0
         RETURN
      ENDIF
      D(1) = TWO*B(1)
      DO 20 I = 2,N-1
         B(I) = X(I+1)-X(I)
         IF (B(I) .LE. ZERO) THEN
            FLAG = -2
            RETURN
         ENDIF
         C(I) = (F(I+1)-F(I))/B(I)
         D(I) = TWO*(B(I)+B(I-1))
   20 CONTINUE
      D(N) = TWO*B(N-1)
!
!     Calculate estimates for the end slopes.  Use polynomials
!     interpolating data nearest the end.
!
      FP1 = C(1)-B(1)*(C(2)-C(1))/(B(1)+B(2))
      IF (N .GT. 3) FP1 = FP1+B(1)*((B(1)+B(2))*(C(3)-C(2)) &
       /(B(2)+B(3))-C(2)+C(1))/(X(4)-X(1))
      FPN = C(N-1)+B(N-1)*(C(N-1)-C(N-2))/(B(N-2)+B(N-1))
      IF (N .GT. 3) FPN = FPN+B(N-1)*(C(N-1)-C(N-2)-(B(N-2) &
       +B(N-1))*(C(N-2)-C(N-3))/(B(N-2)+B(N-3)))/(X(N)-X(N-3))
!
!     Calculate the right-hand-side and store it in C.
!
      C(N) = THREE*(FPN-C(N-1))
      DO 30 K = 2,N-1
         I = N-K+1
         C(I) = THREE*(C(I)-C(I-1))
   30 CONTINUE
      C(1) = THREE*(C(1)-FP1)
!
!     Solve the tridiagonal system.
!
      DO 40 K = 2,N
         P = B(K-1)/D(K-1)
         D(K) = D(K)-P*B(K-1)
         C(K) = C(K)-P*C(K-1)
   40 CONTINUE
      C(N) = C(N)/D(N)
      DO 50 K = N-1,1,-1
         C(K) = (C(K)-B(K)*C(K+1))/D(K)
   50 CONTINUE
!
!     Calculate the coefficients defining the spline.
!
      DO 60 I = 1,N-1
         H = X(I+1)-X(I)
         D(I) = (C(I+1)-C(I))/(THREE*H)
         B(I) = (F(I+1)-F(I))/H-H*(C(I)+H*D(I))
   60 CONTINUE
      B(N) = B(N-1)+H*(TWO*C(N-1)+H*THREE*D(N-1))
      FLAG = 0
      RETURN
      END
!----------------------------------------------------------------
      INTEGER FUNCTION LEFT(LXT,XT,X,MFLAG)
!   LEFT    FINDS INDEX LEFT OF AN ARRAY XT FOR WHICH XT(LEFT)
!           LIES IMMEDIATELY LEFT OF X

!   PURPOSE:
!           FINDS INDEX LEFT OF AN ARRAY XT FOR WHICH XT(LEFT)
!           LIES IMMEDIATELY LEFT OF X

!   INPUT ARGUMENTS:
!     LXT   : NUMBER OF ELEMENTS IN VECTOR XT
!     XT    : VECTOR OF LENGTH LXT STORING THE ABSCISSAE
!     X     : X-VALUE FOR WHICH THE INDEX LEFT IS TO BE FOUND

!   OUTPUT ARGUMENTS:
!     LEFT  : INDEX FOR WHICH XT(LEFT) LIES IMMEDIATELY LEFT OF X
!     MFLAG : FLAG SET IN THE FOLLOWING MANNER
!             LEFT  MFLAG
!              1     -1     IF               X .LT. XT(1)
!              I      0     IF  XT(I)   .LE. X .LT. XT(I+1)
!             LXT     1     IF  XT(LXT) .LE. X

!   METHOD:
!     THAT OF CARL DE BOOR AS DESCRIBED ON PAGE 91 FF. IN:
!     /1/ DE BOOR,C. (1978) A PRACTICAL GUIDE TO SPLINES.
!         APPLIED MATHEMATICAL SCIENCES, VOLUME 27.
!         NEW-YORK-HEIDELBERG-BERLIN: SPRINGER.

!   IMPLEMENTED BY:
!      KRAFT,D., DLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
!                D-8031 OBERPFAFFENHOFEN

!   STATUS: 15. JANUARY 1980

!   SUBROUTINES REQUIRED: NONE

      INTEGER LXT,MFLAG,IHI,ILO,ISTEP,MIDDLE
      DOUBLE PRECISION X,XT(LXT)
      SAVE ILO
      DATA ILO/1/

      IHI=ILO+1
      IF(IHI.LT.LXT)                   GOTO  10
      IF(X.GE.XT(LXT))                 GOTO 110
      IF(LXT.LE.1)                     GOTO  90
      ILO=LXT-1
      IHI=LXT
   10 IF(X.GE.XT(IHI))                 GOTO  40
      IF(X.GE.XT(ILO))                 GOTO 100
      ISTEP=1
   20 IHI=ILO
      ILO=IHI-ISTEP
      IF(ILO.LE.1)                     GOTO  30
      IF(X.GE.XT(ILO))                 GOTO  70
      ISTEP=ISTEP+ISTEP
                                       GOTO  20
   30 ILO=1
      IF(X.LT.XT(1))                   GOTO  90
                                       GOTO  70
   40 ISTEP=1
   50 ILO=IHI
      IHI=ILO+ISTEP
      IF(IHI.GE.LXT)                   GOTO  60
      IF(X.LT.XT(IHI))                 GOTO  70
      ISTEP=ISTEP+ISTEP
                                       GOTO  50
   60 IF(X.GE.XT(LXT))                 GOTO 110
      IHI=LXT
   70 MIDDLE=(ILO+IHI)/2
      IF(MIDDLE.EQ.ILO)                GOTO 100
      IF(X.LT.XT(MIDDLE))              GOTO  80
      ILO=MIDDLE
                                       GOTO  70
   80 IHI=MIDDLE
                                       GOTO  70
   90 MFLAG=-1
      LEFT=1
                                       GOTO 120
  100 MFLAG=0
      LEFT=ILO
                                       GOTO 120
  110 MFLAG=1
      LEFT=LXT
  120                                  RETURN

!   END OF LEFT

      END

