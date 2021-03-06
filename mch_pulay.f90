  subroutine mch_pulay(vin,vout,iter,crys,gs,R0,w_in0,dw,dR,itmax)
  !
  !   2001 adapted by Dave Raczkowski from PEtot
  !
  !     DESCRIPTION:
  !     -----------
  !      
  !    Uses Pulay mixing to get a new Vin and Vout
  !
  use all_to_all_module
  include 'use.h'  
  IMPLICIT NONE
  include 'interface.h' 
  include 'all_to_all.h'  
  include 'flibcalls.ph'    
  !
  !     INPUT:  
  !     -----
  !
  integer, intent(in) :: iter,itmax            ! iteration number
  type(crystal), intent(in) :: crys            ! for vcell, nspin, ztot
  type(parallel_gspace), intent(in) ::  gs  
  complex(dp), intent(inout) :: &
    vin(gs%length,crys%nspin),& ! screening pot. the eigvectors were calc. with
    vout(gs%length,crys%nspin)  ! out screening pot. calc. from the eigvectors 
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp),  intent(inout) :: &
       dw(gs%length,crys%nspin,itmax),& !diff of past input pot from last in
       dR(gs%length,crys%nspin,itmax)   !diff of current vin-vout and past 
  !
  !     WORK:
  !     ------
  !
   complex(dp) :: R0(gs%length,crys%nspin)   
   complex(dp) :: w_in0(gs%length,crys%nspin) 
  !
  !     LOCAL:
  !     ------
  !
  integer, SAVE:: nreset
  real(dp), SAVE ::   AA(30,30)
  integer, SAVE :: nint
  integer is,m,m1,m2,ngrid,i,ierr,len
  real(dp) B(itmax),s,alpha2,w,s1
  real(dp) AA1(itmax,itmax)
  real(dp) :: dmax, dmin, cmax, cmin, warnmax,fac2,dfnorm  
  real(dp), parameter :: small = 1.0d-12 

  alpha2=done
!  alpha2 controls how many recent charge densities
!  to be used. If alpha2 > 1, then, the very old
!  charge density is not used effectivly. 
!  We find that alpha2=1 is O.K.
!******************************************************

  ngrid=gs%fftsize(1)*gs%fftsize(2)*gs%fftsize(3)
  len=gs%length
 
  if(iter.eq.1) then
    nreset=0
    nint=1
  end if 

  if (nreset .eq. 0) then
    nint=iter-nreset
    if(nint.gt.itmax+1) then
!      write(6,*) "restart pulay, iter,itmaxpul",iter,itmax
      nreset=iter-1
      nint=1
    endif
  else 
    nint=nint+1
    if(nint.gt.itmax) then
!      write(6,*) "restart pulay, iter,itmaxpul",iter,itmax
      nreset=iter-1
      nint=1
    endif
  end if

  if(nint.eq.1 .and. nreset.eq.0) then

    do is=1, crys%nspin  
     do i=1,len
       R0(i,is)=vout(i,is)-vin(i,is)
       w_in0(i,is)=vin(i,is)
     enddo
    enddo

  endif

!******************************************************
  if(nint.gt.1 .and. nreset.eq.0 ) then
  
    do is=1,crys%nspin  
     do i=1,len 
       dw(i,is,nint-1)=vin(i,is)-w_in0(i,is)
       dR(i,is,nint-1)=vout(i,is)-vin(i,is)- R0(i,is)
       R0(i,is)=vout(i,is)-vin(i,is)
       w_in0(i,is)=vin(i,is)
     enddo
    enddo


    do m=1,nint-1

     s=real(parallel_zdotc(len * crys%nspin, dR(1, 1,m), 1,&
                        R0(1,1), 1),dp)
      s=s*crys%vcell/ngrid
      B(m)=-s
    enddo 
!******************************************************
    do m=1,nint-1

      s1=real(parallel_zdotc(len * crys%nspin, dR(1, 1,m), 1,&
                        dR(1,1,nint-1), 1),dp)

      s1=s1* crys%vcell/ngrid
      AA(m,nint-1)=s1
      AA(nint-1,m)=s1
    enddo
!    pulay optimization
!**********************************************************
    do m1=1,nint-1
     do m2=1,nint-1
       AA1(m1,m2)=AA(m1,m2)
     enddo
    enddo

    w=done
    do m=nint-1,1,-1
      AA1(m,m)=AA1(m,m)*w
      w=w*alpha2
    enddo

!**********************************************************

    call gaussj(AA1,nint-1,itmax,B,1,1)

! vout=vout-vin
    call mzcopy(crys%nspin*len,R0(1,1),1,vout(1,1),1)   
      
    do m=1,nint-1

! vin=linear combination of input potential differences from initial vin
      call mzaxpy(crys%nspin*len,cmplx(B(m),dzero,dp),dw(1,1,m),1,vin(1,1),1)
! vout=linear combination of differences of differnces of (vin-vout)  
      call mzaxpy(crys%nspin*len,cmplx(B(m),dzero,dp),dR(1,1,m),1,vout(1,1),1)
    enddo

! vout=vout+vin : vout is now a guess at vout
    call mzaxpy(crys%nspin*len,zone,vin(1,1),1,vout(1,1),1)

!******************************************************
  else if(nreset.gt.0 ) then
  
    do is=1,crys%nspin  
     do i=1,len 
       dw(i,is,nint)=vin(i,is)-w_in0(i,is)
       dR(i,is,nint)=vout(i,is)-vin(i,is)- R0(i,is)
       R0(i,is)=vout(i,is)-vin(i,is)
       w_in0(i,is)=vin(i,is)
     enddo
    enddo

    do m=1,itmax
     s=real(parallel_zdotc(len * crys%nspin, dR(1, 1,m), 1,&
                        R0(1,1), 1),dp)
      s=s*crys%vcell/ngrid
      B(m)=-s
    enddo 
!******************************************************
    do m=1,itmax !nint-1

      s1=real(parallel_zdotc(len * crys%nspin, dR(1, 1,m), 1,&
                        dR(1,1,nint), 1),dp)


      s1=s1* crys%vcell/ngrid
      AA(m,nint)=s1
      AA(nint,m)=s1
    enddo
!    pulay optimization
!**********************************************************
    do m1=1,itmax !nint-1
     do m2=1,itmax !nint-1
       AA1(m1,m2)=AA(m1,m2)
     enddo
    enddo

!    w=done
!    do m=nint-1,1,-1
!      AA1(m,m)=AA1(m,m)*w
!      w=w*alpha2
!    enddo

!**********************************************************

    call gaussj(AA1,itmax,itmax,B,1,1)

! vout=vout-vin
    call mzcopy(crys%nspin*len,R0(1,1),1,vout(1,1),1)   
      
    do m=1,itmax !nint-1

! vin=linear combination of input potential differences from initial vin
      call mzaxpy(crys%nspin*len,cmplx(B(m),dzero,dp),dw(1,1,m),1,vin(1,1),1)
! vout=linear combination of differences of differnces of (vin-vout)  
      call mzaxpy(crys%nspin*len,cmplx(B(m),dzero,dp),dR(1,1,m),1,vout(1,1),1)
    enddo

! vout=vout+vin : vout is now a guess at vout
    call mzaxpy(crys%nspin*len,zone,vin(1,1),1,vout(1,1),1) 

  endif

  return
  end
!      
!
!---------------  MCH_KERK ----------------------------------
!
!
  subroutine mch_kerk(vin,vout,vion,iter,gs,alphamix)

  use all_to_all_module
  include 'use.h'  
  IMPLICIT NONE
  include 'interface.h' 
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) ::  gs 
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp), intent(inout) :: &
       vion(gs%length), &    ! the screened ionic potential
       vin(gs%length),&      ! screening pot. the eigvectors were calc. with   
       vout(gs%length)       ! out screening pot. calc. from the eigvectors  
  real(dp),intent(in) :: &
       alphamix ! mixing parameter for strength of preconditioning
  !
  !     LOCAL:
  !     ------
  !
  integer i,iter

!********************************************************
!*** Kerker mixing, this is like the precoditioning,
!*** and it adds the new component w_out-w_in into w_in
!********************************************************
  call mzaxpy(gs%length,zmone,vin(1),1,vout(1),1)

  do i = 1, gs%length  
    if (gs%ekin(i) .gt. dzero) then  
     vin(i)= vin(i) +  alphamix*vout(i)*gs%ekin(i)/(gs%ekin(i)+0.5)
     vion(i) = vion(i) + vin(i)
    end if

  end do

!**************************************************
  return
  end

       SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
! this is double precision gauss program
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=200)
      real*8 A(NP,NP),B(NP,MP)
      integer IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG= 0.d0
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                call mystop( 'Singular matrix 1' )
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW

        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.d0) then
          call mystop( 'Singular matrix. 2' )
        endif
        PIVINV=1.d0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.d0
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.d0
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END


      
      
      
      
  
     
