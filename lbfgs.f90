!     ----------------------------------------------------------------------
!     This file contains the LBFGS algorithm and supporting routines
!
!     ****************
!     LBFGS SUBROUTINE
!     ****************
!
  SUBROUTINE LBFGS(N,M,G,DIAG,W,ITER,POINT,NPT)
  include 'use.h' 
  IMPLICIT NONE
  include 'flibcalls.ph'  
!
  INTEGER N,M
  REAL(dp) G(N),DIAG(N),W(N*(2*M+1)+2*M),metric(N*N)
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
! 
!     This subroutine solves the unconstrained minimization problem
! 
!                      min F(x),    x= (x1,x2,...,xN),
!
!      using the limited memory BFGS method. The routine is especially
!      effective on problems involving a large number of variables. In
!      a typical iteration of this method an approximation Hk to the
!      inverse of the Hessian is obtained by applying M BFGS updates to
!      a diagonal matrix Hk0, using information from the previous M steps.
!      The user specifies the number M, which determines the amount of
!      storage required by the routine. The user may also provide the
!      diagonal matrices Hk0 if not satisfied with the default choice.
!      The algorithm is described in "On the limited memory BFGS method
!      for large scale optimization", by D. Liu and J. Nocedal,
!      Mathematical Programming B 45 (1989) 503-528.
! 
!      The user is required to calculate the function value F and its
!      gradient G. In order to allow the user complete control over
!      these computations, reverse  communication is used. The routine
!      must be called repeatedly under the control of the parameter
!      IFLAG. 
!
!      The steplength is determined at each iteration by means of the
!      line search routine MCVSRCH, which is a slight modification of
!      the routine CSRCH written by More' and Thuente.
! 
!      The calling statement is 
! 
!          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
! 
!      where
! 
!     N       is an INTEGER variable that must be set by the user to the
!             number of variables. It is not altered by the routine.
!             Restriction: N>0.
! 
!     M       is an INTEGER variable that must be set by the user to
!             the number of corrections used in the BFGS update. It
!             is not altered by the routine. Values of M less than 3 are
!             not recommended; large values of M will result in excessive
!             computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
!     X       is a DOUBLE PRECISION array of length N. On initial entry
!             it must be set by the user to the values of the initial
!             estimate of the solution vector. On exit with IFLAG=0, it
!             contains the values of the variables at the best point
!             found (usually a solution).
! 
! 
!     G       is a DOUBLE PRECISION array of length N. Before initial
!             entry and on a re-entry with IFLAG=1, it must be set by
!             the user to contain the components of the gradient G at
!             the point X.
! 
! 
!     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
!             then on initial entry or on re-entry with IFLAG=2, DIAG
!             it must be set by the user to contain the values of the 
!             diagonal matrix Hk0.  Restriction: all elements of DIAG
!             must be positive.
! 
!     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
!             workspace for LBFGS. This array must not be altered by the
!             user.
!
! 
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  REAL(dp) GTOL,ONE,ZERO,GNORM,STP1,FTOL,STPMIN, &
                      STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM,dotprod,YS2,YY2
  INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO, &
            BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
!
  SAVE
  DATA ONE,ZERO/1.0D+0,0.0D+0/
!
!     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!     ---------------------------------------
!     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!         OTHER TEMPORARY INFORMATION.
!     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!         IN THE FORMULA THAT COMPUTES H*G.
!     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!         STEPS.
!     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!         GRADIENT DIFFERENCES.
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.

  ISPT= N+2*M
  IYPT= ISPT+N*M  
   
  BOUND=ITER-1   
  IF (ITER .GT. M)BOUND=M

!  YS= MDDOT(9,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!  YY= MDDOT(9,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)

!  if (YY .ne. 0d0) then
!    DO  I=1,9
!      DIAG(I)= YS/YY
!    end do

!  else
!    DO  I=1,9
!      DIAG(I)= 0d0
!    end do
!  end if

!  YS2= MDDOT(N-9,W(IYPT+NPT+10),1,W(ISPT+NPT+10),1)
!  YY2= MDDOT(N-9,W(IYPT+NPT+10),1,W(IYPT+NPT+10),1)

  YS= MDDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
  YY= MDDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)

    DO  I=1,N
      DIAG(I)= YS/YY
    end do

!  if (YY2 .ne. 0d0) then
!    DO  I=10,N
!      DIAG(I)= YS2/YY2
!    end do

!    YS=YS2+YS
!  else
!    DO  I=10,N
!      DIAG(I)= 0d0
!    end do
!  end if
  !
  !     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
  !     "Updating quasi-Newton matrices with limited storage",
  !      Mathematics of Computation, Vol.24, No.151, pp. 773-782.
  !     --------------------------------------------------------- 
  !
  CP= POINT
  IF (POINT.EQ.0) CP=M
  W(N+CP)= ONE/YS

  DO I=1,N
    W(I)= -G(I)
  end do

  CP= POINT

  DO I= 1,BOUND
    CP=CP-1
    IF (CP.EQ. -1)CP=M-1

    SQ= MDDOT(N,W(ISPT+CP*N+1),1,W(1),1)
    INMC=N+M+CP+1
    IYCN=IYPT+CP*N
    W(INMC)= W(N+CP+1)*SQ
    CALL MDAXPY(N,-W(INMC),W(IYCN+1),1,W(1),1)
  end do

  DO I=1,N
    W(I)=DIAG(I)*W(I)
  end do

  DO I=1,BOUND

    YR= MDDOT(N,W(IYPT+CP*N+1),1,W(1),1)
    BETA= W(N+CP+1)*YR
    INMC=N+M+CP+1
    BETA= W(INMC)-BETA
    ISCN=ISPT+CP*N
    CALL MDAXPY(N,BETA,W(ISCN+1),1,W(1),1)
    CP=CP+1
    IF (CP.EQ.M)CP=0
  end do

  RETURN
  END
