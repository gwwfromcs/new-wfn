!
!     find the lowest eigenvectors and eigenvalues for a given
!     Hamiltonian with a conjugate gradient minimization scheme
!
subroutine diagonalize_hamiltonian_gcg(time,iter,num_HPsi,pw_params, neigmax,maxit,&
     iconv, neigmin,  neig, delta, eval, ham, x)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(hamiltonian), intent(inout) :: ham    ! the hamiltonian to diagonalize
  type(pw_parameter), intent(in) :: pw_params            ! various parameters

  integer, intent(in) :: maxit, &                  ! max number of iterations
       neigmin, &               ! minimum number of eigenvalues to be computed
       neigmax                  ! maximum number of eigenvalues to be computed

  real(dp), intent(in) :: delta                 ! accuracy of diagonalization
  integer iter
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  integer, intent(inout) :: neig,&    ! max/real number of eigenvalues computed
            num_HPsi                  ! op count of H*Psi  

  complex(dp), intent(inout) :: &
       x(*)                      ! start/final wave functions on input/output
  real(dp), intent(inout) :: &
       eval(neig + 1)                 ! start/final eigenvalues, eval(neig+1)
  real(dp), intent(inout) ::time(num_time)
  !
  !     OUTPUT:
  !     -------
  !
  logical, intent(out) :: iconv               ! true if converged, false else
  !
  !
  !     --------------------- local variables ----------------------
  !
  real(dp) :: t0, t1, t2, t3, &                            ! timing variables
       oldrho, a, b, c, d, demax, de, lambda, eprecond, delrho,rho_init,theta,&
       rho, &                                                 ! current error
       nops, &                                         ! number of operations
       egap                ! gap between highest computed and next eigenvalue
  real(dp)   num,energy_old,energy1,denom
  integer :: myproc, ierrmall, mnew, mmax,nover_cnt, &
       diaginfo, &                  ! return value of diagonalization routine
       m, &                                       ! current number of vectors
       m_squared, &                                             ! square of m
       i, j, &                                                      ! dummies
       imax, &                                            ! largest gap index
       n, &                                                ! iteration number
       nreset, &                                 ! number of resets performed
       len                                                    ! gspace length
  complex(dp) :: dclambda, dclambdasq, dcgamma, phase
  !
  !     work arrays
  !
  complex(dp), allocatable :: &                               ! of size m x m
       eps(:), s(:)
  complex(dp), allocatable :: &                            ! of size nanl x m
       ews(:,:)  
  complex(dp), allocatable :: &                             ! of size len x m
       g(:,:), gold(:,:), gp(:,:), h(:,:), y(:,:)   
  integer, allocatable :: &                                       ! size of m
       iord(:)
  real(dp), allocatable :: rwork(:)  

  logical :: ishift, ireset, iconjug 
  !
  !          local routines
  real(dp), external :: trace, contract, average_ekin, gimmetime  


  t0 = gimmetime()  
  if (iand(pw_params%output(1), 8) == 8) t3 = t0
  m = neig                                     ! start with maximum allowed
  m_squared = m * m  
  len = ham%gspace%length  
  myproc = ham%gspace%myproc  
  !
  !     allocate all the work spaces
  !
  allocate(eps(m * m)) ; eps = zzero 
  allocate(s(m * m)) ; s = zzero  
  allocate(rwork(max(1, 3*m-2))) ; rwork = dzero  
  if (ham%vnloc%nvecs > 0 .and. .not. ham%pspot%NL_rspace(1)) then  
     allocate(ews(ham%vnloc%nvecs, m))  
     ews = zzero
  end if
  allocate(g(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'g', len, m, ierrmall  
     call mystop  
  end if
  g = zzero  
  allocate(gold(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gold', len, m, ierrmall  
     call mystop  
  end if
  gold = zzero
  allocate(h(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'h', len, m, ierrmall  
     call mystop  
  end if
  h = zzero
  allocate(y(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'y', len, m, ierrmall  
     call mystop  
  end if
  y = zzero  
  allocate(gp(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gp', len, m, ierrmall  
     call mystop  
  end if
  gp = zzero  

  allocate(iord(m))
  iord = 0  

  ireset = .true.                                ! first step is like a reset
  iconjug = .false.                      ! don't make conjugate in first step
  if (neig .gt. neigmin) then 
    ishift = .true.               ! want to shift if error has decreased enough
  else
    ishift = .false.
  end if

  rho = 1.0d4                                       ! give some initial error
  rho_init=dzero
  nops = 0  

  t1 = dzero
  egap = eval(neig + 1) - eval(neig)  

  if (egap < dzero) then                        ! no good gap as input, guess
     write(9, *) ' *** DIAGONALIZATION: INPUT EIGENVALUES BAD!'  
     ham%shift = dzero 
     call apply_ham_wavefn(time,pw_params%output(1), y(1, 1), ham, x(1), &
          m, ews(1, 1), pw_params%nbandsfft)                      ! y = H * x 
     num_HPsi=num_HPsi+1
     call overlap_matrix(eps(1), x(1), y(1, 1), m, len)         ! eps = x^T*y
     egap = dmfour * eps(neigmin + (neigmin - 1) * m)  
     write(9, *) ' *** GUESSED GAP:', egap  
  end if

  if (iter .eq. 1) then


    call cheap_overlap_matrix(s(1), x(1), m, len) 

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
    call occsp_sqrtinv_scal('U','n',myproc,s,m,ham%gspace%nproc)

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
    call mztrmm('R','U','N','N',len,m,zone,s(1),m,x(1),len)   

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
  end if
  !
  !     shift hamiltonian just enough to make it negative definite
  !
  ham%shift = dzero !eval(neig) + egap * dqtr + pw_params%shiftsafety  

  if (iand(pw_params%output(1), 536870912 ) == 536870912 ) then
     write(9,*) '               RHO               ENERGY'
     write(9,*) '            --------           ---------'    
  end if

  nreset = 0
  nover_cnt=0


  do n = 1, maxit                   ! ----------- start of big iteration loop
100  continue                                                       ! restart

     if (ireset) then  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        call apply_ham_wavefn(time,pw_params%output(1), y(1, 1), ham, x(1), &
             m, ews(1, 1), pw_params%nbandsfft)                     ! y=H * x  

        num_HPsi=num_Hpsi+1
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime() 
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
        call cheap_overlap_matrix(s(1), x(1), m, len)             ! s = x^T*x

        call overlap_matrix(eps(1), x(1), y(1, 1), m, len)      ! eps = x^T*y

        if (iand(pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + 1.2d1 * real(m * m * len, dp)
        end if
        ireset = .false.                              ! no resetting any more
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
   
        energy1=dzero
        do i=1,m
        energy1=energy1+eps(i+(i-1)*m)
        end do        

     else  

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)

        call cheap_overlap_matrix(s(1),x(1),m,len) ! s = x^T*x

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

        call occsp_sqrtinv_scal('U','n',myproc,s,m,ham%gspace%nproc)

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
        call mztrmm('R','U','N','N',len,m,zone,s(1),m,x(1),len)     

        call mztrmm('R','U','N','N',len,m,zone,s(1),m,y(1,1),len)

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

        !
        !  Trabsform search direction so that when ass to gradient is in same
        !  representation as wavefunctions (X) - not reaaly beneficial to do
        !
!        call mztrmm('R','U','N','N',len,m,zone,s(1),m,h(1,1),len)

        call overlap_matrix(eps(1),x(1),y(1,1),m,len) ! eps = x^T*y

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

        energy1=dzero
        do i=1,m
        energy1=energy1+eps(i+(i-1)*m)
        end do        

     end if
     !
     !        adjust shift, find gap. In principle, this is a hack,
     !        because the energy surface changes, and we should
     !        restart the algorithm. In practice, we find that it
     !        hardly has any impact on the performance.
     !
     if(ishift.and.(rho.lt.0.1).and. &            ! adjust shift of 
         iand(pw_params%optimize,256).eq.256) then       ! Hamiltonian
 
       call occsp_diag_scal(myproc,eps,m,ham%gspace%nproc,eval,s, &
                            m_squared,rwork)

       do i=1,m ; iord(i)=i; end do
       call sort(m,eval(1),iord(1)) ! sort with numerical recipes
       if(neigmin.gt.m) then
         write(9,*) 'diagonalize_hamiltonian: neigmin > m.!'
         write(9,*) neigmin, m
         call mystop
       elseif(neigmin.lt.m) then
         !
         !        find largest gap, and throw other vectors away
         !
         demax=0.d0
         imax=neigmin+1
         do i=neigmin+1,m
           de=(eval(i)-eval(i-1))/(dble(i)**2.5)
           if(demax.lt.de) then
             demax=de; imax=i-1
           endif
         end do
         mnew=imax        
         !
         !        superimpose x to form eigenvectors, and rearrange 
         if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
         call mzgemm('N','N',len,m,m,zone,x(1),len, &
                   eps(1),m,zzero,g(1,1),len)

         if(iand(pw_params%output(1),1).eq.1) then
           t1=t1+gimmetime()-t2
           nops = nops + 8*dble(m*m)*dble(len)
         endif
         call wavefn_rearrange(len,m,g(1,1),x(1),iord(1),mnew)

         egap = eval(mnew+1)-eval(mnew)
         ham%shift = dzero !ham%shift+eval(mnew) +egap/4 + pw_params%shiftsafety
         iconjug = .false. ! don't reuse old conjugate vector
         ireset  = .true. ! recompute H psi from scratch
         write(9,*) m,mnew,'  new subspace'
         m = mnew
         ishift = .false.

         goto 100 

       else                ! neigmin.eq.m, no rearrangement

         ham%shift = dzero !ham%shift+eval(m)+egap/4+ pw_params%shiftsafety
         ireset = .true.  ! recompute H psi from scratch
         ishift = .false.
         goto 100
       endif
     endif

     if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
     call mzgemm('N','N',len,m,m, zmone, x(1),len, eps(1), m,zzero,g(1,1),len) 

     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

     call mzaxpy(len*m, zone, y(1,1), 1, g(1,1), 1)

     !
     !        make conjugate
     !
     !prepare for precondi

     call mzcopy(len * m, g(1, 1), 1, gp(1, 1), 1)  

     !
     !        change the preconditioning energy during the iteration.
     !        this is not good style, since it effectively changes the
     !        matrix under consideration, and we lose the property of
     !        conjugation. we tried it out, and it doesn't hurt performance.
     !        the factor of 1.35 was found by numerical experiment to give
     !        best performance. Smaller factors are dangerous, as they kill
     !        good components, and that hurts a lot.
     !
          eprecond=average_ekin(x(1),s(1), &
             ham%gspace%ekin(1),m,len)*1.35d0

     call precondition(gp(1, 1), eprecond, ham%ekinmod(1), &
          ham%gspace%ekin(1), m, len)


     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     !
     !  Project out occupied space from gp (the conditioned gradient) 
     !  then inverse precondition to get the gradient equivalent
     !    Not so important for here, but in diag_ham_metal_gcg for metals it is
!     call overlap_matrix(s(1),x(1),gp(1,1),m,len) ! eps = x^T*y
!     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
!     call mzgemm('N','N',len,m,m, zmone, x(1),len,s(1), m, zone,gp(1,1),len) 
! 
!     call mzcopy(len*m,gp(1,1),1,g(1,1),1)
!     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
!     call precondition_inv(g(1,1),eprecond,ham%ekinmod(1), &
!          ham%gspace%ekin(1), m,len)

     rho = dble(parallel_zdotc(len*m,gp(1,1),1,g(1,1),1)) 
     if((n.gt.1).and. (.not.ireset)) then
!      step for Polak-Ribiere. Fletcher-Reeves would not require gold,
!      but is performing terribly in many cases. 
       delrho = dble(parallel_zdotc(len*m,gold(1,1),1,g(1,1),1))
       dcgamma =cmplx((rho-delrho)/oldrho,0.d0,kind=8)
!        save gradient into old gradient
       call mzcopy(len*m, gp(1,1), 1, gold(1,1), 1) 
!     dcgamma = cmplx(rho/oldrho) ! -------- do this for fletcher-Reeves
       call mzaxpy(len*m, dcgamma, h(1,1), 1, gp(1,1), 1)
     endif

     if (iand(pw_params%output(1), 536870912 ) == 536870912 ) &
        write(9,*) 'gcg iter #',n+nover_cnt, rho,energy1
     call myflush(9)
     call mzcopy(len*m, gp(1,1), 1, h(1,1), 1)
  
     oldrho = rho  
     if (rho < delta) then !! .or. rho .lt. .2*rho_init) then  
        if (n.eq.1) return
        iconv = .true.  
        goto 200  
     end if

     if (n .eq. 1) rho_init=rho 

     num = dble(parallel_zdotc(len*m,g(1,1),1,h(1,1),1))

     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     call apply_ham_wavefn(time,pw_params%output(1), g(1,1),ham, &
             h(1,1),m,ews(1,1),pw_params%nbandsfft)  ! g=H * h
            num_HPsi=num_Hpsi+1
     if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
              
     call mzgemm('N','N',len,m,m,zmone,h(1,1),len,eps(1),m,zzero,gp(1,1),len)  
        
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

     call mzaxpy(len*m, zone, g(1,1), 1, gp(1,1), 1)
          
     denom=dble(parallel_zdotc(len*m,h(1,1),1,gp(1,1),1))
     lambda=num/denom
     
     if (lambda .gt. done) lambda=done
     if (lambda .lt. 0d0) then
       write(9,*) lambda,num,denom,' reset'
       nreset = nreset+1
       ireset=.true.
       if(nreset+n.gt.maxit+pw_params%itdiag_add) goto 199 ! too many resets ... bye bye
       go to 100
     else
       dclambda=lambda      
       ireset=.false.
     end if
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     !
     !        move to new position
     !
     call mzaxpy(len*m, -dclambda, h(1,1), 1, x(1), 1)
     call mzaxpy(len*m,-dclambda, g(1,1), 1, y(1,1), 1)

     if (n .eq. maxit .and. (rho .gt. .3*rho_init .or. rho/m .gt. 1d-3) .and. nover_cnt .lt. pw_params%itdiag_add)then
       nover_cnt=nover_cnt+1
       go to 100
     end if
     !cjp write(9, 991) n, m, rho, ham%shift  
  end do                              ! ----------- end of big iteration loop

991  format('iter',I5,' m=',I5,' rho=',G14.8,' shift',F14.8)  
199  continue  


  iconv = .false.

200  continue  
  !  
  !     the vectors x we have now span the right subspace, but are not yet
  !     eigenvectors. to get those, we diagonalize the epsilon matrix.
  !     the eigenvalues are the eigenvalues we are looking for.
  !     the eigenvectors tell us how to superimpose the vectors x to get
  !     the true eigenvalues.
  !     There are parallel issues here: If the eps matrix is not
  !     exactly identical, the eigenvectors will acquire a different
  !     phase on different processors. Therefore, we diagonalize only
  !     on proc 0, and broadcast the result.
  !
  !
  !     in case the calculation did not converge, perform
  !     the orthonormalization of the  subspace: 
  !     
  !     a) Compute S^-1/2 = AT * D^-1/2 * A   
  !     b) Compute x_new= S^-1/2 * x
  !     
  !     The eigenvalue array is used as temporary storage
  !
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
  if (.not. iconv) &
    write(9,'(a,i3,a,g12.6)') ' after ',n+nover_cnt-1,' iter., residual error = ',rho

!  if (rho .gt. 1d0) then
    call cheap_overlap_matrix(s,x,m,len);
    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

    call occsp_sqrtinv_scal('U','n',myproc,s,m,ham%gspace%nproc)

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
    call mztrmm('R','U','N','N',len,m,zone,s(1),m,x(1),len)     
    call mztrmm('R','U','N','N',len,m,zone,s(1),m,y(1,1),len)
    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
!  end if

  call overlap_matrix(eps(1),x(1),y(1,1),m,len) ! eps = x^T*y

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

  call  occsp_diag_scal(myproc,eps,m,ham%gspace%nproc,eval,g,len*m,rwork)
  !
  !     combine the eigenvectors into g
  !
  if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
  call mzgemm('N','N',len,m,m, zone, x(1),len, eps(1),m, zzero, g(1,1), len) 
  if(iand(pw_params%output(1),1).eq.1) then
    t1=t1+gimmetime()-t2
    nops = nops + 8*dble(m*m)*dble(len)
  endif
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
  call mzcopy(len*m, g(1,1), 1, x(1), 1); !   x <- g 
  !
  !     save the eigenvalues
  !
  do i=1,m
    s(i)=cmplx(eval(i)+ham%shift,0.d0,kind=8)
  end do
  neig = m                  ! return number of eigenvectors calculated 
  !
  !     retrieve eigenvalues
  !
  eval(1:m)=dble(s(1:m))
  eval(m+1)=eval(m)+egap    ! and remember gap
  !
  !     performance printout
  !
  if (iand(pw_params%output(1), 1) == 1) write(9, 920) t1, nops / t1 * 1.0d-6
  !
  !     timing printout
  !
  if (iand (pw_params%output(1), 8) == 8) then  
     write(9, 910) n+nover_cnt-1, gimmetime() - t0  
     call myflush(9)  
  end if
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
  !
  !     deallocate the work spaces
  !
  deallocate(eps)
  deallocate(s) 
  deallocate(rwork)
  if (ham%vnloc%nvecs > 0 .and. .not. ham%pspot%NL_rspace(1) ) deallocate(ews)
  deallocate(g)
  deallocate(gold)
  deallocate(h)
  deallocate(gp)
  deallocate(y)
  deallocate(iord)

  return
  
910  format(' DIAGONALIZATION TOOK ',i3,' ITERATIONS AND ', &
          &     f12.3, ' SECONDS')
920  format(' TIME FOR MATRIX-MATRIX:',f12.6,', MFLOPS=',f12.6)  
930  format('*** DIAGONALIZE_HAMILTONIAN: ALLOCATE(',a,'(',i6, &
          &        ',',i4,') FAILED. ERROR CODE=',i4)

end subroutine diagonalize_hamiltonian_gcg

      subroutine precondition_inv(v, eprecond, ekinmod,ekin,m,len)
      use all_to_all_module
      implicit none
      include 'all_to_all.h'

!
!     INPUT:
!     -----

      integer &
          m,  &                 ! number of vectors to precondition
          len                  ! length of vectors

      real(dp)  &
          eprecond,      &      ! the preconditioning energy cutoff
          ekin(len),   &        ! kinetic energy of the gvectors |k+G|
          ekinmod(3)           ! kinetic energy modification parameters

!
!     INPUT/OUTPUT:
!     ------------

      complex(dp) &
          v(len,m)             ! vectors to be preconditioned

!
!     DESCRIPTION:
!     preconditions vectors v(len,m) with the preconditioner 
!     of Teter, Payne, and Allen
!

!     ---------- local variables -----------------
!
!
!
      integer i,j
      complex(dp) gt
      real(dp)  pcfn, precfn
      external precfn
      double precision myderf


      external myderf

      do i=1,m
         do j=1,len
            pcfn = precfn((ekin(j)+ekinmod(1) &
                *(1.0+myderf((ekin(j)-ekinmod(1)) &
                /ekinmod(2))))/eprecond)
            v(j,i)=v(j,i)*(1.0/pcfn)
         end do
      end do
      return

      end subroutine precondition_inv

