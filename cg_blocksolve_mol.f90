Subroutine cg_blocksolve_mol(imonitor,crys,ham,maxit,delta,n_uk,e_k,uk,&
     uktil,phi,nbandsfft)
  !
  !     1996 by Bernd Pfrommer and Francesco Mauri
  !     
  !     1999 Modified by Chris Pickard, molecular limit
  !     

  Use all_to_all_module
  Include 'use.h'
  Implicit None
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'
  !
  !     INPUT:
  !     -----

  Type (hamiltonian) :: ham  ! the hamiltonian at k
  Type (crystal)     :: crys ! need the bvec information for gradient

  Integer :: imonitor  ! monitor flag indicating additional printout
  Integer :: n_uk      ! number of u_k 
  Integer :: maxit     ! max number of iterations
  Integer :: nbandsfft ! max number of bands to FFT simultaneously

  Real(dp) :: e_k(n_uk) ! energy eigenvalues of u_k
  Real(dp) :: delta     ! convergence criterion
  Real(dp) :: time(6)   ! 

  Complex(dp) :: phi(ham%dim,n_uk) ! contains grad|u_k> Will be destroyed!
  Complex(dp) :: uk(ham%dim,n_uk)  ! the full occupied states at k
  !
  !     INPUT/OUTPUT:
  !     -------------

  Complex(dp) :: uktil(ham%dim, n_uk) ! the desired wave function uktilde

  !
  !     DESCRIPTION:
  !     -----------

  !
  !     solves the equation
  !     
  !     (H_k - e_k*I + shift* Q) |u~_k> = -(1-Q)grad|u_k>
  !     
  !     by means of a conjugate gradient algorithm
  !      
  !
  !     --------------------- local variables ----------------------
  !
  Integer :: i,j, k, l, iwork
  Integer :: nunconv   ! number of unconverged states
  Integer :: nconvp1   ! number of converged states plus one
  Integer :: nnconvp1  ! new number of unconverged states plus one
  Integer :: iter      ! the iteration number
  Integer :: len       ! length of full vectors
  Integer :: noccup    ! number of occupied states
  Integer :: ilumo     ! index of lumo
  Integer :: ihomo     ! index of homo

  Real(dp) :: dgamma   ! gamma 
  Real(dp) :: t0,gimmetime
  Real(dp) :: a,c

  Real(dp), Allocatable :: rhoold(:)   ! previous iterations rho
  Real(dp), Allocatable :: rho(:)      ! absolute of the gradient
  Real(dp), Allocatable :: eprecond(:)
  Real(dp), Allocatable :: shift(:) ! energy shift of the Hamiltonian[Ryd]

  Complex(dp) :: dcone, dczero, dcmone
  Parameter(dcone=(1.d0,0.d0), dczero=(0.d0,0.d0), dcmone=(-1.d0,0.d0))


  Complex(dp), Allocatable :: ukapp(:)     ! packed version of uk
  Complex(dp), Allocatable :: g(:,:)       ! the gradient
  Complex(dp), Allocatable :: gp(:,:)      ! the gradient
  Complex(dp), Allocatable :: h(:,:)       ! step direction
  Complex(dp), Allocatable :: hold(:,:)    ! prev iterations step dirn 
  Complex(dp), Allocatable :: t(:,:)       ! some work vector
  Complex(dp), Allocatable :: work(:)      ! work arrays
  Complex(dp), Allocatable :: mat(:)       ! intermediate matrix
  Complex(dp), Allocatable :: gradukapp(:) ! gradient of approximate uk
  Complex(dp), Allocatable :: me_k(:)      ! minus e_k, in complex(dp) form

  Complex(dp) :: ssum      ! debugging sum
  Complex(dp) :: tkinav    ! work variable for preconditioning
  Complex(dp) :: dcgamma   ! conjugation gamma
  Complex(dp) :: dclambda  ! step length lambda, in complex(dp) form

  Real(dp) :: average_ekin

  External gimmetime, average_ekin

  !
  !     ------------------------------------------------------
  !

  t0=gimmetime()

  len    = ham%gspace%length
  ilumo  = n_uk+1             ! index of lumo
  ihomo  = n_uk               ! index of homo
  noccup = n_uk


  !
  !     allocate a work array that is big enough for all purposes
  !
  iwork =Max(ham%vnloc%nvecs*noccup, &! for apply_ham_wfn
       noccup*noccup)                 ! for apply_proj_wfn

  Allocate(work(iwork))
  work=Cmplx(0.d0,0.d0)
  !
  !     -------- compute the rhs of the equation to solve. this is -------
  !
  !               psi = -(1-Q) * phi

  Call apply_proj_wavefn(len,-1.d0,1,uk(1,1),noccup,phi(1,1),noccup,&
       phi(1,1),work(1)) ! compute (1-Q) grad(u_k)

  Call mzdscal(len*noccup,-1.d0,phi(1,1),1)
  !
  !     Starting guess; uktil=0
  !

  uktil=Cmplx(0.d0,0.d0)

  !
  !     now uktil contains the full trial vectors
  !

  !
  !     ------------ here starts the conjugate gradient part --------
  !      

  !
  !     figure out how much the Hamiltonian must be shifted to make it 
  !     positive  definite.
  !
  Allocate(shift(noccup))

  shift(1:noccup) = (e_k(1:noccup)-e_k(1)) + 1.d0

  Allocate(g(len,noccup))          ! used as scratch array here
  Allocate(gp(len,noccup))
  Allocate(h(len,noccup))
  Allocate(hold(len,noccup))
  Allocate(t(len,noccup))
  Allocate(eprecond(noccup))
  Allocate(rho(noccup)) ; Allocate(rhoold(noccup))
  Allocate(me_k(noccup))
  me_k(1:noccup)=Cmplx(-e_k(1:noccup),0.d0,kind=8)

  Do i=1, noccup
     eprecond(i) = 1.35d0*average_ekin(uk(1,i),tkinav,ham%gspace%ekin(1),1,len)
  End Do

  nconvp1 = 1
  nunconv = noccup - nconvp1 +1

  Do iter=1, maxit
     !
     !        compute the gradient. can reuse information from previous step
     !     
     !

     If(iter.Eq.1) Then 
        !    g=H*uktil
        Call apply_ham_wavefn(time,imonitor,g(1,nconvp1),ham,uktil(1,nconvp1),&
             nunconv,work(1),nbandsfft)
        !    g=g + shift *Q   * uktil
        Call apply_proj_wavefn(len,shift(nconvp1),nunconv,uk(1,1), &
             noccup,uktil(1,nconvp1),nunconv,g(1,nconvp1),work(1)) 
        !    g=g -e_k * uktil
        Do i=nconvp1, noccup
           Call mzaxpy(len,me_k(i),uktil(1,i),1,g(1,i),1)
        End Do
        !    g = g - b
        Call mzaxpy(len*nunconv,dcmone,phi(1,nconvp1),1,g(1,nconvp1),1)
     Endif

     !
     !        compute residual
     !
     gp(1:len,nconvp1:noccup) = g(1:len,nconvp1:noccup)
     nnconvp1 = nconvp1
     Do i=nconvp1, noccup
        Call precondition(gp(1,i),eprecond(i),ham%ekinmod(1),&
             ham%gspace%ekin(1), 1,len)
        rho(i)=Dble(parallel_zdotc(len,gp(1,i),1,g(1,i),1))

        If((rho(i).Lt.delta).And.(i.Eq.nnconvp1)) Then
           nnconvp1=nnconvp1+1
        Endif
     End Do

     nconvp1=nnconvp1
     nunconv = noccup - nconvp1 +1

     If(nconvp1.Eq.noccup+1) Goto 100

     !
     !        compute the step direction h. Conjugate it to previous step
     !
     If(iter.Eq.1) Then
        h(1:len, nconvp1:noccup)=-gp(1:len, nconvp1:noccup)
     Else
        Do i=nconvp1,noccup
           dgamma = rho(i)/rhoold(i)
           h(1:len,i)=-gp(1:len,i) + dgamma*hold(1:len,i)
        End Do
     Endif

     rhoold(nconvp1:noccup)=rho(nconvp1:noccup)

     hold(1:len,nconvp1:noccup) = h(1:len,nconvp1:noccup)
     !
     !        compute t = A*h
     !
     Call apply_ham_wavefn(time,imonitor,t(1,nconvp1),ham,h(1,nconvp1),&
          nunconv,work(1),nbandsfft) ! t = H*h
     Call apply_proj_wavefn(len,shift(nconvp1),nunconv,uk(1,1),&
          noccup,h(1,nconvp1),nunconv,t(1,nconvp1),work(1)) ! t=t+shift*Q*h
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda

     Do i=nconvp1, noccup
        Call mzaxpy(len,me_k(i),h(1,i),1,t(1,i),1) !t=t-e_k*h

        a=Dble(parallel_zdotc(len,h(1,i),1,g(1,i),1))
        c=Dble(parallel_zdotc(len,h(1,i),1,t(1,i),1))
        dclambda=Cmplx(-a/c,0.d0,kind=8)
        !
        !           move to new position   
        !
        !           uktil = uktil+lambda*h
        Call mzaxpy(len,dclambda,h(1,i),1,uktil(1,i),1) 

        !           update to get the gradient
        Call mzaxpy(len,dclambda,t(1,i),1, g(1,i),1) ! g = g + lambda*t
     End Do

  End Do

  Write(9,*) '*** WARNING: cg_blocksolve not converged in',maxit,' iterations'
  Write(9,*) 'residual errors:', rho

100 Continue

  Deallocate(me_k)
  Deallocate(rho)   ;   Deallocate(rhoold)
  Deallocate(shift) ;   Deallocate(eprecond)
  Deallocate(g)     ;   Deallocate(gp)
  Deallocate(h)
  Deallocate(hold)
  Deallocate(t)
  Deallocate(work)
  If(Iand(imonitor,8).Eq.8) Then
     Write(9,930) iter, gimmetime()-t0
     Call myflush(9)
  Endif


900 Format('iter:',i3,' rho=', 1000g12.6)
910 Format('cg_solve converged at iteration:',i3,' with residual error rho=',&
         g12.6)
920 Format(' STATES ',i3,'-',i3,' TOOK ',i3,' ITERATIONS ')
930 Format(' CGSOLVE TOOK ',i3,' ITERATIONS AND ',f12.3, ' SECONDS')

End Subroutine cg_blocksolve_mol





