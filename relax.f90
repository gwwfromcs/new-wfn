! Originally relax.c
!
!  program to relax the atoms and the lattice vectors such that
!  the total energy acquires a minimum
!  
!  1996 Bernd Pfrommer
!
!  Modified Jan. 1999
!  No major changes -- more of a spring clean while porting
!  porting to a Hitach SR2201.
!
!  Greg McMullan
!  hpcf 
!  University of Cambridge
!  
!  Merged with existing code by Peter Haynes Jun 1999
!
!  Changed over to Fortran90 by Peter Haynes March 2000
!  
!  made into subroutin from main program. D. Raczkowski 2001
!
!------------------------------------------------------------------------
subroutine relax(f1,plen,H,H0,g,x,&
             t0, pressure, relax_myproc,pwcount,accuracy, &
            iupdateh, relaxepsfac, relaxcoordfac,decoup_rc,lambdalimit,&
            relaxstress, relaxforce, iwfnreuse_c, adjustpressure,&
           relax_method,stdout,&
          lambda_tolerance, lambda2_tolerance,itmax,fixed_vector,vector_plane,&
          pw_params,altkpoints,bands,altbands,syms,energs,crys,vinit)
  !
  !     ---------------------------------------------------------------
  !
  use constants
  use pw_parameter_module
  use kpoint_module
  use band_module
  use symmetry_module
  use energy_module
  use crystal_module
  use esdf
  use all_to_all_module
  use molecular_dynamics_module
  implicit none
  include 'flibcalls.ph' 
  include 'all_to_all.h'
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(pw_parameter) :: pw_params ! plane wave parameters
  type(kpoint) :: altkpoints  ! alternative BZ integration kpoints
  type(band) :: bands ! eigenvalues, occupation numbers etc
  type(band) :: altbands    ! another set of eigenvalues etc..
  type(symmetry) :: syms ! symmetry operations
  type(energy) :: energs    ! energies and other information
  type(crystal) :: crys     ! crystal structure
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) ::& 
        plen,relax_myproc,itmax,iupdateh,decoup_rc,stdout,relax_method

  real(dp), intent(in) :: &  
     relaxstress, relaxforce,lambda_tolerance, lambda2_tolerance,&
     fixed_vector(3),vector_plane(3),accuracy, relaxepsfac, relaxcoordfac,&
       lambdalimit,vinit
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  integer, intent(inout) :: iwfnreuse_c, adjustpressure,pwcount
  real(dp), intent(inout) :: pressure,t0,f1
  real(dp),intent(inout) :: H(plen*plen),H0(plen*plen),x(plen),g(plen)
  !
  !     -------------------------- local variables -----------------------
  !
  real(dp) funeval,f0,& ! function returns enthalpy , old enthalpy
           lambda,&  ! step size
           rho,&     ! sqrt(g*g) est. of error from relaxed configuration
           lvect_inv(3,3),& ! inverse of transpose of lvec 
           lvec_inv(3,3),&   ! inverse of lvec
           atam1(3,3),&    ! lvect_inv * lvec_inv
           dotprod,& ! dot product with given metric 
           matinvert,& ! inverts lvec matrix and returns volume
           ldelta,lg,g0,g1, gnew,lgnew,& ! inner products 
           lambdalin,&  ! lambda after 1st refinement
           cvar,avar,bvar,& ! used for 2nd refinement of step size
           discrim,L_p1,L_p2,& ! used in 2nd refinement of step size
           lambda_lbfgs,&  ! current lambda saved for lbfgs
           adot,bdot,adoti,bdoti,& ! used in update of H in BFGS
           temp 

  ! for extrapolation
  real(dp) in_pr11,in_pr12,in_pr13,in_pr22,in_pr23,in_pr33, inv_fac

  real(dp), allocatable :: &
           delta(:),&     ! search (update) direction 
           gold(:), &     ! old gtradient - stresses and forces
           xold(:), &     ! old positions - unit cell and atomic pos.
           pv(:),&        ! BFGS - work : forces in lat. coord. - for lbfgs
           q(:),&         ! BFGS - work : old forces in lat. units
           s(:),&         ! BFGS - work
           u(:),&         ! BFGS - work
           metric(:),&    ! metric for inner products
           del_r(:,:)     ! differnces in x - used for extrapoaltion

  integer max_lbfgs,&  !maximum # of past steps used for new search direction
   extrap,updateh,ISPT,IYPT,POINT,i,j,k,NPT,iter ! counters and flags
  integer, parameter :: logfile_c = 7
  integer, parameter :: convergefile = 8

  !
  !     -------------------------- BEGIN  -----------------------
  !
  allocate (delta(plen)); delta=zzero;
  allocate (metric(plen*plen)); metric=zzero;
  allocate (del_r(plen,4)); del_r=zzero;
  allocate (pv(plen)); pv=zzero;
  allocate (q(plen)); q=zzero;
  allocate (s(plen)); s=zzero;
  allocate (u(plen)); u=zzero;
  allocate (gold(plen)); gold=zzero;
  allocate (xold(plen)); xold=zzero;


!  write(34,*) int(real(plen**2-plen,dp)/real(2+2*plen)),' max_lbfgs'
  max_lbfgs = int(real(plen**2-plen,dp)/real(2+2*plen))

  call matvec(plen, H, g, delta)

  gold = g
  xold = x

  if (relax_method .eq. 2) then
    lvect_inv = transpose(crys%lvec)
    temp = matinvert(lvect_inv(1, 1))
    lvec_inv =  crys%lvec
    temp = matinvert(lvec_inv(1, 1))
    call matmul(atam1,lvect_inv,lvec_inv)

    q(1:9)=gold(1:9)
    call mdgemm('N','N',3,(plen-9)/3,3,done,atam1(1,1),3,g(10),&
               3,dzero,pv(10),3) 

  end if

  lambda = done

  call computemetric(plen, xold, xold, metric, 3, crys%lvec)
  rho = sqrt(dotprod(plen, gold, gold, metric))

  if (relax_myproc == 0) then
    write(logfile_c, '(A,F9.6)') 'rho = ', rho
    write(convergefile, '(A,I3,A,F9.6,A,F14.6)') 'iteration ', 0, &
          ' rho = ', rho, ' enthalpy = ', f1
  end if

  if (rho <= accuracy) then
    if (relax_myproc == 0) then
      write(logfile_c, '(A)') &
           'already done at iteration 0!'
      if (iand(pw_params%output(1), 134217728) == 134217728) &
           call print_xyz(0, crys%ntype, crys%mxdatm, crys%natom, &
           crys%nameat, 2, x(10), crys%lvec, f1)
    end if

   return
  end if

  if (itmax < 2 .and. relax_myproc == 0 .and. &
     iand(pw_params%output(1), 134217728) == 134217728) &
     call print_xyz(0, crys%ntype, crys%mxdatm, crys%natom, &
     crys%nameat, 2, x(10), crys%lvec, f1)

  ! ============================ begin main loop ===================

  ISPT= plen+2*max_lbfgs      !DBR
  IYPT= ISPT+plen*max_lbfgs
  POINT=0

  do iter = 1, itmax - 1

    pw_params%ilog = iter
    if (relax_myproc == 0) then
      write(logfile_c, '(//A,I3,A)') '--------- iteration ', iter, &
           ' -----------'
      if (iand(pw_params%output(1), 134217728) == 134217728) &
           call print_xyz(iter - 1, crys%ntype, crys%mxdatm, crys%natom, &
           crys%nameat, 2, x(10), crys%lvec, f1)
    end if

    x = xold

    call update(plen, x, delta, dmone)                  ! x = x - delta
    call fix_rotation(x, fixed_vector, vector_plane)

    f0 = f1
    call computemetric(plen, xold, xold, metric, 1, crys%lvec)   
    g0 = -dotprod(plen, delta, g, metric)
    call computemetric(plen, xold, xold, metric, 2, crys%lvec)
    ldelta = sqrt(dotprod(plen, delta, delta, metric))

    if (pw_params%extrapolation .ge. 0) then 

      extrap=mod( pw_params%extrapolation , 10)

      if (iter .eq. 1) then
        del_r(:,1)=x(:)-xold(:)
        pw_params%MD_alpha=dzero
        pw_params%MD_beta=dzero         
      else if (iter .eq. 2 .or. extrap .eq. 0) then
        del_r(:,2)= del_r(:,1)
        del_r(:,1)=x(:)-xold(:)
        pw_params%MD_alpha= done
        pw_params%MD_beta=dzero
      else 
        del_r(:,3)= del_r(:,2)
        del_r(:,2)= del_r(:,1)
        del_r(:,1)=x(:)-xold(:)

        in_pr11 = mddot(plen-9, del_r(10,1),1, del_r(10,1), 1)
        in_pr12 = mddot(plen-9, del_r(10,1),1, del_r(10,2), 1)
        in_pr13 = mddot(plen-9, del_r(10,1),1, del_r(10,3), 1)
        in_pr22 = mddot(plen-9, del_r(10,2),1, del_r(10,2), 1)
        in_pr23 = mddot(plen-9, del_r(10,2),1, del_r(10,3), 1)
        in_pr33 = mddot(plen-9, del_r(10,3),1, del_r(10,3), 1)

        inv_fac = 1/( in_pr22* in_pr33 - in_pr23* in_pr23)

        if (mod(iter,extrap) .eq. 0) then
          pw_params%MD_alpha= inv_fac* (in_pr33*in_pr12 - in_pr23*in_pr13)
          pw_params%MD_beta= inv_fac* (-in_pr23*in_pr12 + in_pr22*in_pr13)
        else
          pw_params%MD_alpha=done
          pw_params%MD_beta=dzero
        end if

      end if 

    else

      pw_params%MD_alpha=dzero
      pw_params%MD_beta=dzero 
    end if

    write(34,*) pw_params%MD_beta, pw_params%MD_alpha, ' alpha'


    f1 = funeval(t0, plen, x, pressure, g, relax_myproc, & 
           relaxstress, relaxforce, &
         iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
        bands,altbands,syms,energs,crys,vinit)
    pwcount = pwcount + 1

    call computemetric(plen, xold, x, metric, 1, crys%lvec)
    g1 = -dotprod(plen, delta, g, metric)

    call computemetric(plen, x, x, metric, 3, crys%lvec)
    rho = sqrt(dotprod(plen, g, g, metric))
    lg = rho

    if (rho <= accuracy) goto 999

    lambdalin = -g0 / (g1 - g0)
    lambda = lambdalin
    updateh = iupdateh

    if (lambda < 0) then
      if (relax_myproc == 0) write(logfile_c, '(A)') &
           '!!!! set lambda to 1 and initialized H!!!!!' 
      call initialize_hessian(plen, x, H, relaxepsfac, relaxcoordfac, -1, &
           decoup_rc, crys%lvec)
      H0 = H
      updateh = 0       ! don't update H
      lambda = done
    end if
    if (lambda > lambdalimit) then
      if (relax_myproc == 0) write(logfile_c, '(A,2(F9.6,A))') &
           '!!!! reduced lambda from ', lambda, ' to ', lambdalimit, ' !!!!!'
      lambda = lambdalimit
    end if

    if (relax_myproc == 0) then
      write(logfile_c, '(2(A,F12.6,A,G12.6))') &
           'f0=', f0, ' g0=', g0, ' f1=', f1, ' g1=', g1
      write(logfile_c, '(A,F9.6)') 'lambda from linear=', lambda

      write(logfile_c, '(A,G12.6)') &
           'trialstep: normalized product between vectors (g1,delta): ', &
           -g1 / ldelta / lg
    end if
    !
    ! if step size (lambda) is not within tolerance than do refinement
    ! Also if energy is lower do refinement
    !
    write(34,*) lambda,' lambda 0'
    if(abs(lambda - done) > lambda_tolerance .or. f0 .lt. f1 ) then
      ! take shorter step if necessary
      x = xold
      call update(plen, x, delta, -lambda)
      call fix_rotation(x, fixed_vector, vector_plane)
      if (relax_myproc == 0) write(logfile_c, '(A)') '####### move to new x'

      f1 = funeval(t0, plen, x, pressure, g, relax_myproc, & 
           relaxstress, relaxforce, &
         iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
        bands,altbands,syms,energs,crys,vinit)

      call computemetric(plen, xold, x, metric, 1, crys%lvec)
      gnew = -dotprod(plen, delta, g, metric)
      call computemetric(plen, x, x, metric, 3, crys%lvec)
      rho = sqrt(dotprod(plen, g, g, metric))
      lgnew = rho
      if (rho <= accuracy) then
         if (relax_myproc == 0) &
              write(logfile_c, '(A)') 'converged at new x!'
         goto 999
      end if

      if (relax_myproc == 0) write(logfile_c, '(A,G12.6)') &
           'linear: normalized product between vectors (gnew,delta): ', &
           -gnew / ldelta / lgnew
      pwcount = pwcount + 1

      cvar = g0
      avar = (g1 - g0) * (g1 - g0) * gnew / (g0 * g1)
      bvar = g1 - avar - g0

      if (relax_myproc == 0) write(logfile_c, '(3(A,F9.6))') &
           'a=', avar, ' b=', bvar, ' c=', cvar

      discrim = bvar * bvar - dfour * avar * cvar
      if (discrim < dzero) then
         if (relax_myproc == 0) write(logfile_c, '(A)') &
              'Error: complex roots. Taking linear lambda.'
         lambda = lambdalin
      else
         if (abs(avar) < abs(cvar * 1.0d-6)) then
            if (relax_myproc == 0) write(logfile_c, '(A)') &
                 'quadratic term very small, doing linear fit.'
            lambda = -cvar / bvar
         else
            L_p1 = (-bvar + sqrt(discrim)) / (dtwo * avar)
            L_p2 = (-bvar - sqrt(discrim)) / (dtwo * avar)
            if (abs(L_p1 - lambdalin) < abs(L_p2 - lambdalin)) then
               lambda = L_p1
            else
               lambda = L_p2
            end if
         end if
      end if
      if (relax_myproc == 0) write(logfile_c, '(A,F9.6)') &
           'lambda from parabola=', lambda
      !
      ! If refined step size is not within 2_tolerance than do refinement 
      !
      write(34,*) lambda,' lambda 1'
      if (abs(lambda - lambdalin) > lambda2_tolerance * lambda ) then
         if (relax_myproc == 0) write(logfile_c, '(A)') &
              '################# moving to parabolic x'
         x = xold
         call update(plen, x, delta, -lambda)
         call fix_rotation(x, fixed_vector, vector_plane)

         f1 = funeval(t0, plen, x, pressure, g, relax_myproc, & 
            relaxstress, relaxforce, &
            iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
            bands,altbands,syms,energs,crys,vinit)
         pwcount = pwcount + 1

         call computemetric(plen, xold, x, metric, 1, crys%lvec)
         gnew = -dotprod(plen, delta, g, metric)
         call computemetric(plen, x, x, metric, 3, crys%lvec)  
         rho = sqrt(dotprod(plen, g, g, metric))
         lgnew = rho

         if (relax_myproc == 0) write(logfile_c, '(A,G12.6)') &
              'parabolic: normalized product between vectors &
              &(gnew,delta): ', -gnew / ldelta / lgnew
         if (rho <= accuracy) then
            if (relax_myproc == 0) write(logfile_c, '(A)') &
                 'converged at new parabolic x!'
            goto 999
         end if
         lambda_lbfgs=lambda  !DBR
      else
         lambda = lambdalin
         lambda_lbfgs=lambda  !DBR
         if (relax_myproc == 0) write(logfile_c, '(A)') &
              'Did not move to lambda quadratic!'
      end if

    else

      lambda_lbfgs=done  !DBR
      if (relax_myproc == 0) then
        write(logfile_c, '(A)') 'Trial step is close enough to &
              &estimated minimum!'
        write(logfile_c, '(A,G12.6)') &
              'Normalized dot product between (g1,delta): ', -g1 / ldelta / lg
      end if
    end if

999 continue

    if (relax_myproc == 0) then
      write(convergefile, '(A,I3,A,F9.6,A,F14.6)') 'iteration ', iter, &
             ' rho = ', rho, ' enthalpy = ', f1
      write(logfile_c, '(A,F9.6)') 'rho of new minimum position = ', rho
      write(logfile_c, '(A,I3,A)') 'So far ', pwcount, &
             ' call(s) to the pw program!'
    end if
    !
    ! ------------ update H ---------------------
    !
    if (relax_method .eq. 1) then

      if (updateh == 1) then
         pv = -lambda * delta
         q = g - gold
         call matvec(plen, H, q, s)          ! s = H*q

         adot = dzero ; bdot = dzero
         do i = 1, plen
	   adot = adot + pv(i) * q(i)
	   bdot = bdot + q(i) * s(i)
         end do
         adoti = done / adot ; bdoti = done / bdot

         u = adoti * pv - bdoti * s

         k = 0 
         do i = 1, plen
           do j = 1, plen
              k = k + 1
	      H(k) = H(k) + pv(i) * pv(j) * adoti - s(i) * s(j) * bdoti + &
			u(i) * u(j) * bdot
	   end do
         end do
      end if

      if (relax_myproc == 0) call checkpoint(plen, x, g, H, H0, f1)
      if (rho <= accuracy) exit

      ! ------------ apply H ---------------------

      call matvec(plen, H, g, delta)

      if (relax_myproc == 0) then
        write(logfile_c, '(A)') 'negative of suggested step direction:'
        call printpacked(plen, delta, crys%lvec)
      end if

    else

      call mdscal(plen,-lambda_lbfgs,delta(1),1)

      NPT=POINT*plen
      call mdcopy(plen,delta(1),1,H(ISPT+NPT+1),1)

      lvect_inv = transpose(crys%lvec)
      temp = matinvert(lvect_inv(1, 1))
      lvec_inv =  crys%lvec
      temp = matinvert(lvec_inv(1, 1))
      call matmul(atam1,lvect_inv,lvec_inv)

      pv(1:9)=g(1:9)
      call mdgemm('N','N',3,(plen-9)/3,3,done,atam1(1,1),3,g(10),&
               3,dzero,pv(10),3)  

      do I=1,plen
        H(IYPT+NPT+I)= pv(I)-q(I)
      end do

      q=pv
      POINT=POINT+1
      if (POINT.EQ.max_lbfgs)POINT=0

      if (relax_myproc == 0) call checkpoint(plen, x, g, H, H0, f1)
      if (rho <= accuracy) exit

      call LBFGS(plen,max_lbfgs,pv(1),gold,H,ITER+1,POINT, NPT)
  
      call mdcopy(plen,H(1),1,delta(1),1)  
      call mdscal(plen,dmone,delta(1),1)

    end if

    ! ----------- save old variables ----------------

    gold = g
    xold = x

    if (relax_myproc == 0) then
        call myflush(stdout)
        call myflush(logfile_c)
        call myflush(convergefile)
    end if

  end do

  deallocate(delta,metric,del_r,pv,q,s,u,gold,xold)


  return
  end subroutine relax

! ============================ support routines ===================

subroutine pack(ntyp, mxdatm, plen, natom, statevec, coord, eps)

  use constants
  implicit none

  integer, intent(in) :: ntyp, mxdatm, plen, natom(ntyp)
  real(dp), intent(out) :: statevec(plen)
  real(dp), intent(in) :: coord(3 * mxdatm * ntyp), eps(9)
  integer :: i, j, k, l, m

  ! write epsilon into right position

  statevec(1:9) = eps(1:9)

  l = 9
  do i = 1, ntyp
     do j = 1, natom(i)
        do k = 1, 3
           l = l + 1
           statevec(l) = coord(3 * (i - 1) * mxdatm + 3 * (j - 1) + k)
        end do
     end do
  end do

end subroutine pack
! ---------------------------------------------------------------
subroutine unpack( plen, statevec, crys, eps)

  use constants
  use crystal_module
  implicit none
  type(crystal) :: crys     ! crystal structure

  integer, intent(in) :: plen
  real(dp), intent(in) :: statevec(plen)
  real(dp), intent(out) ::  eps(9)


  integer :: i, j, k, l, m

  eps(1:9) = statevec(1:9)

  l = 9
  do i = 1, crys%ntype
     do j = 1, crys%natom(i)
        do k = 1, 3
           l = l + 1
!           coord(3 * (i - 1) * mxdatm + 3 * (j - 1) + k) = statevec(l)
          crys%coord(k,j,i)=statevec(l)
        end do
     end do
  end do

  crys%rat=crys%coord

end subroutine unpack
!  ---------------------------------------------------------------
integer function packlen(ntyp, natom)

  implicit none

  integer, intent(in) :: ntyp, natom(ntyp)
  integer :: i

  packlen = 9

  do i = 1, ntyp
     packlen = packlen + 3 * natom(i)
  end do

end function packlen
! ---------------------------------------------------------------
subroutine printpacked(xl, statevec, lvec)

  use constants
  implicit none

  integer, parameter :: logfile_c = 7
  integer, intent(in) :: xl
  real(dp), intent(in) :: statevec(xl), lvec(3, 3)

  integer :: i, j, k
  real(dp) :: rc(3), a(3, 3)

  write(logfile_c, '(A)') 'strain part:'
  write(logfile_c, '(3(G12.6,X)/,3(G12.6,X)/,3(G10.4,X))') statevec(1:9)

  if (xl > 9) then
     write(logfile_c, '(A)') 'coordinate components:'
     do i = 10, xl, 3
        write(logfile_c, '(3(G12.6,X))') statevec(i:i + 2)
     end do

     write(logfile_c, '(A)') 'in cartesian coordinates:'
     call apply_eps_lat(statevec, lvec, a)
     do i = 9, xl-1, 3
        do j = 1, 3
           rc(j) = dzero
           do k = 1, 3
              rc(j) = rc(j) + a(k, j) * statevec(i + k)
           end do
        end do
        write(logfile_c, '(3(G12.6,X))') rc(1:3)
     end do
  end if

end subroutine printpacked
! ---------------------------------------------------------------
subroutine printunpacked(ntyp, mxdatm, natom, coord, eps)

  use constants
  implicit none

  integer, parameter :: logfile_c = 7
  integer, intent(in) :: ntyp, mxdatm, natom(ntyp)
  real(dp), intent(in) :: coord(3 * mxdatm * ntyp), eps(9)
  integer :: i, j

  write(logfile_c, '(A)') 'strain part:'
  write(logfile_c, '(3(G10.4,X)/,3(G10.4,X)/,3(G10.4,X))') eps(1:9)

  do i = 1, ntyp
     do j = 1, natom(i)
        write(logfile_c, '(3(F9.5,X))') &
             coord(3 * (i - 1) * mxdatm + 3 * (j - 1) + 1 : &
             3 * (i - 1) * mxdatm + 3 * (j - 1) + 3)
     end do
  end do

end subroutine printunpacked
! ---------------------------------------------------------------
subroutine apply_eps_lat(eps, a, c)

  use constants
  implicit none

  real(dp), intent(in) :: eps(9), a(3, 3)
  real(dp), intent(out) :: c(3, 3)

  !  epsilon is assumed to have the FORTRAN convention:
  !
  !      e_11, e_12, e_13, e_21, e_22, e_23, e_31, e_32, e_33
  !
  !       This was different before. Changed 11/16/96. BP

  integer :: i, j, k
  real(dp) :: hilf(3, 3)

  k = 0
  do i = 1, 3
     do j = 1, 3
        k = k + 1
        hilf(j, i) = eps(k)
     end do
     hilf(i, i) = hilf(i, i) + done
  end do

  call matmul(c, hilf, a)

end subroutine apply_eps_lat
! ---------------------------------------------------------------
function dotprod(xl, a, b, met)

  use constants
  implicit none

  real(dp) :: dotprod
  integer, intent(in) :: xl
  real(dp), intent(in) :: a(xl), b(xl), met(xl * xl)
  integer :: i, j, k

  dotprod = dzero
  k = 0
  do i = 1, xl
     do j = 1, xl
        k = k + 1
        dotprod = dotprod + a(i) * met(k) * b(j)
     end do
  end do

end function dotprod
! ---------------------------------------------------------------
subroutine matvec(xl, H, b, c)

  use constants
  implicit none

  integer, intent(in) :: xl
  real(dp), intent(in) :: H(xl * xl), b(xl)
  real(dp), intent(out) :: c(xl)
  integer :: i, j, k

  k = 0
  do i = 1, xl
     c(i) = dzero
     do j = 1, xl
        k = k + 1
        c(i) = c(i) + H(k) * b(j)
     end do
  end do

end subroutine matvec
! ---------------------------------------------------------------
subroutine sym2full(symmat, fullmat)

  use constants
  implicit none

  real(dp), intent(in) :: symmat(6)
  real(dp), intent(out) :: fullmat(3, 3)

  fullmat(1, 1) = symmat(1)
  fullmat(2, 2) = symmat(2)
  fullmat(3, 3) = symmat(3)
  fullmat(1, 2) = symmat(4)
  fullmat(2, 1) = symmat(4)
  fullmat(2, 3) = symmat(5)
  fullmat(3, 2) = symmat(5)
  fullmat(1, 3) = symmat(6)
  fullmat(3, 1) = symmat(6)

end subroutine sym2full
! ---------------------------------------------------------------
subroutine full2sym(symmat, fullmat)

  use constants
  implicit none

  real(dp), intent(out) :: symmat(6)
  real(dp), intent(in) :: fullmat(3, 3)

  symmat(1) = fullmat(1, 1)
  symmat(2) = fullmat(2, 2)
  symmat(3) = fullmat(3, 3)
  symmat(4) = dhalf * (fullmat(1, 2) + fullmat(2, 1))
  symmat(5) = dhalf * (fullmat(2, 3) + fullmat(3, 2))
  symmat(6) = dhalf * (fullmat(1, 3) + fullmat(3, 1))

end subroutine full2sym
! ---------------------------------------------------------------
subroutine matmul(c, a, b)

  use constants
  implicit none

  real(dp), intent(out) :: c(3, 3)
  real(dp), intent(in) :: a(3, 3), b(3, 3)
  integer :: i, j, k

  do i = 1, 3
     do j = 1, 3
        c(j, i) = dzero
        do k = 1, 3
           c(j, i) = c(j, i) + a(k, i) * b(j, k)
        end do
     end do
  end do

end subroutine matmul
! ---------------------------------------------------------------
subroutine mattrans(a,b)

  use constants
  implicit none

  real(dp), intent(in) :: a(3, 3)
  real(dp), intent(out) :: b(3, 3)
  integer :: i, j, k

  do i = 1, 3
     do j = 1, 3
        b(j, i) = a(i, j)
     end do
  end do

end subroutine mattrans
! ---------------------------------------------------------------
subroutine pstress(eps, p, pst, vinit)

  use constants
  implicit none

  real(dp), intent(in) :: eps(3, 3), vinit, p
  real(dp), intent(out) :: pst(3, 3)
  real(dp) :: v0  ! inverse of scaled volume
  integer :: i, j
  real(dp) :: matinvert

  do i = 1, 3
     do j = 1, 3
        pst(j, i) = eps(i, j) 
     end do
     pst(i, i) = pst(i, i) + done
  end do

  ! calculate 1/V * p * dV/d(eps)

  v0 =  matinvert(pst)
  v0 = v0 * vinit

  ! convert pressure from GPa to Ry/a.u.^3

  pst = pst * v0 * p / 14710.75d0

end subroutine pstress
! ---------------------------------------------------------------
subroutine checkpoint(plen, x, g, H, H0, f1)

  use constants
  implicit none

  integer, intent(in) :: plen
  real(dp), intent(in) :: x(plen), g(plen), H(plen * plen), H0(plen * plen), f1

  character(len=10), parameter :: checkpointfile = 'CHECKPOINT'
  integer, parameter :: ckpfile = 10
  integer :: i, ios

  open(unit = ckpfile, file = checkpointfile, iostat = ios, &
       form = 'formatted', status = 'replace')
  if (ios /= 0) call mystop( 'Cannot open checkpoint file!' )

  write(ckpfile, '(A)') 'x='
  do i = 1, plen
     write(ckpfile, '(G30.20)') x(i)
  end do

  write(ckpfile, '(A)') 'g='
  do i = 1, plen
     write(ckpfile, '(G30.20)') g(i)
  end do

  write(ckpfile, '(A)') 'H='
  do i = 1, plen * plen
     write(ckpfile, '(G30.20)') H(i)
  end do

  write(ckpfile, '(A)') 'H0='
  do i = 1, plen * plen
     write(ckpfile, '(G30.20)') H0(i)
  end do

  write(ckpfile, '(A10,G30.20)') 'enthalpy =', f1

  close(unit = ckpfile)

end subroutine checkpoint
! ---------------------------------------------------------------
subroutine update(plen, x, delta, lambda)

  use constants
  implicit none

  integer, intent(in) :: plen
  real(dp), intent(inout) :: x(plen)
  real(dp), intent(in) :: delta(plen), lambda
  real(dp) :: extraterm(9)

  x = x + lambda * delta

  call matmul(extraterm,delta(1:9),x(1:9))
  x(1:9) = x(1:9) + lambda*extraterm

end subroutine update
! ---------------------------------------------------------------
subroutine cross_product(c,a,b)
  
  use constants
  implicit none
  
  real(dp), intent(in) :: a(3), b(3)
  real(dp), intent(out) :: c(3)
  
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
  
end subroutine cross_product
! ---------------------------------------------------------------
function vector_length(vector)
  
  use constants
  implicit none
  
  real(dp) :: vector_length
  real(dp), intent(in) :: vector(3)
  vector_length = sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
end function vector_length
! ---------------------------------------------------------------
subroutine fix_rotation(eps, fixed_vector, vector_plane)
  
  use constants
  implicit none
  
  real(dp), intent(inout) :: eps(3,3)
  real(dp), intent(in) :: fixed_vector(3), vector_plane(3)
  
  real(dp) :: rot(3,3), nhat(3), other(3), &
       fixed_vector_prime(3), other_prime(3), epsrot(3,3), &
       vector_plane_prime(3), nhat_prime(3)
  integer :: i, j
  
  real(dp) :: vector_length
  
  integer, parameter :: logfile_c = 7
  
  if (vector_length(fixed_vector) .eq. dzero) then
     return
  end if
  if (vector_length(vector_plane) .eq. dzero) then
     return
  end if
  
  do i = 1, 3
     eps(i,i) = eps(i,i) + done
  end do
  
  call matvec(3,eps,fixed_vector, fixed_vector_prime)
  fixed_vector_prime = fixed_vector_prime/vector_length(fixed_vector_prime)
  
  call cross_product(nhat,fixed_vector, fixed_vector_prime)
  if (vector_length(nhat) > 1e-6) then
     nhat = nhat/vector_length(nhat)
     call cross_product(other,fixed_vector,nhat)
     call cross_product(other_prime,fixed_vector_prime,nhat)
     do i=1,3
        do j=1,3
           rot(j,i) = nhat(i)*nhat(j) + fixed_vector(i)*fixed_vector_prime(j) + &
                other(i)*other_prime(j)
        end do
     end do
     call matmul(epsrot,rot,eps)
  else
     epsrot = eps
  end if
  call matvec(3,epsrot,vector_plane, vector_plane_prime)
  call cross_product(nhat, fixed_vector, vector_plane)
  nhat = nhat/vector_length(nhat)
  call cross_product(nhat_prime, fixed_vector, vector_plane_prime)
  nhat_prime = nhat_prime/vector_length(nhat_prime)
  call cross_product(other,fixed_vector,nhat)
  call cross_product(other_prime,fixed_vector,nhat_prime)
  do i=1,3
     do j=1,3
        rot(j,i) = nhat(i)*nhat_prime(j) + fixed_vector(i)*fixed_vector(j) + &
             other(i)*other_prime(j)
     end do
  end do
  call matmul(eps,rot,epsrot)
  
  do i = 1, 3
     eps(i,i) = eps(i,i) - done
  end do
  
end subroutine fix_rotation
! ---------------------------------------------------------------
subroutine computemetric(plen, x, xp, met, imode, lvec)

  use constants
  implicit none

  integer, intent(in) :: plen, imode
  real(dp), intent(in) :: x(plen), xp(plen), lvec(3, 3)
  real(dp), intent(out) :: met(plen * plen)
  integer :: i, j, k
  real(dp) :: eps(3, 3), epsp(3, 3), a(3, 3), ap(3, 3), mmat(3, 3), &
       dum(3, 3),dummy
  real(dp) :: matinvert

  k = 0
  do i = 1, 3
     do j = 1, 3
        k = k + 1
        eps(j, i) = x(k) 
        epsp(j, i) = xp(k)
     end do
  end do

  call apply_eps_lat(eps, lvec, a)
  call apply_eps_lat(epsp, lvec, ap)

  !   mode 1: compute metric between covariant*contravariant =   
  !            (at)*(a't)^-1                                      
  !   mode 2: compute metric between contravariant*contravariant =   
  !            (at)*(a')                                      
  !   mode 3: compute metric between covariant*covariant =   
  !            ((at)*(a'))^-1

  select case(imode)
  case(1)
     dummy = matinvert(ap)
     call matmul(dum, ap, a)
     call mattrans(dum, mmat)
  case(2)
     call mattrans(a, dum)
     call matmul(mmat, dum, ap)
  case(3)
     call mattrans(a, dum)
     call matmul(mmat, dum, ap)
     dummy = matinvert(mmat)
  case default
     call mystop( 'Invalid metric mode in computemetric!' )
  end select

  met = dzero
  do i = 1, 9
     met((i - 1) * plen + i) = done
  end do

  do i = 10, plen, 3
     do j = 1, 3
        do k = 1, 3
           met((i + j - 2) * plen + i + k - 1) = mmat(k, j)
        end do
     end do
  end do

  !  FPRINTF(logfile_c,"Metric matrix for mode %d\n",imode);
  !  for(i=0;i<plen;i++) {
  !    for(j=0;j<plen;j++) {
  !      FPRINTF(logfile_c,"%5.2f ", met[i*plen+j]);
  !      FPRINTF(logfile_c,"\n");
  !    }
  !  }

end subroutine computemetric
! ---------------------------------------------------------------
subroutine initialize_hessian(plen, x, H, relaxepsfac, relaxcoordfac, mode, &
     decoup_rc, lvec)

  use constants
  implicit none

  integer, intent(in) :: plen, mode, decoup_rc
  real(dp), intent(in) :: x(plen), relaxepsfac, relaxcoordfac, lvec(3, 3)
  real(dp), intent(out) :: H(plen * plen)
  real(dp) :: a(3, 3), atam1(3, 3), at(3, 3), ainv(3, 3), rc(3), dummy
  integer :: i, j, k, natoms, ic, iat, l
  real(dp) :: matinvert

  call apply_eps_lat(x, lvec, a)

  ainv = a
  dummy = matinvert(ainv)

  call mattrans(a, at)

  call matmul(atam1, at, a)

  if (mode == -1) dummy = matinvert(atam1)
  if (mode == 0) atam1 = a

  H = dzero
  do i = 1, 9
     H((i - 1) * plen + i) = relaxepsfac
  end do

  do i = 10, plen, 3
     do j = 1, 3
        do k = 1, 3
           H((i + j - 2) * plen + i + k - 1) = relaxcoordfac * atam1(k, j)
        end do
     end do
  end do

  ! initialize the off-diagonal elements if decoupling is desired 
  ! 
  !    notice that  epsilon_ij =>    epsilon[3*j+i]
  !   (fortran convention)

  if (decoup_rc == 1) then
     natoms = (plen - 9) / 3
     do iat = 1, natoms  ! loop over atoms
        do l = 1, 3      ! compute cartesian coordinate
           rc(l) = dzero
           do k = 1, 3
              rc(l) = rc(l) + a(k, l) * x(9 + (iat - 1) * 3 + k)
           end do
        end do
        do ic = 1, 3     ! loop over the vector components
           do k = 1, 3
              do l = 1, 3
                 H(9 + (iat - 1) * 3 + ic + plen * (3 * (k - 1) + l - 1)) = &
                      -relaxepsfac * ainv(k, ic) * rc(l)
                 H((9 + (iat - 1) * 3 + ic - 1) * plen + 3 * (k - 1) + l) = &
                      -relaxepsfac * ainv(k, ic) * rc(l)
              end do
           end do
        end do
     end do
  end if


end subroutine initialize_hessian
! ---------------------------------------------------------------
subroutine print_xyz(iter, ntyp, mxdatm, natom, atm, alen, coord, &
     lvec, enthalpy)

  use constants
  implicit none

  integer, intent(in) :: iter, ntyp, mxdatm, natom(ntyp), alen
  real(dp), intent(in) :: coord(3, mxdatm, ntyp), lvec(3, 3), &
       enthalpy
  character(len=alen), intent(in) :: atm(ntyp)

  character(len=9), parameter :: xyz_file = 'STRUC_XYZ'

  integer :: ityp, natoms, iatom, i, j
  real(dp) :: cart(3)

  natoms = 0
  do ityp = 1, ntyp
     natoms = natoms + natom(ityp)
  end do

  if (iter == 0) then
     open(31, file = xyz_file, status = 'replace')
  else
     open(31, file = xyz_file, status = 'old', position = 'append')
  end if

  write(31, *) natoms
  write(31, 100) iter, enthalpy

  do ityp = 1, ntyp
     do iatom = 1, natom(ityp)
        cart = dzero
        do i = 1, 3
           do j = 1, 3
              cart(i) = cart(i) + 0.529177274d0 * lvec(j, i) * &
                   modulo(coord(j, iatom, ityp), done)
           end do
        end do
        write(31, 200) atm(ityp)(1:2), cart(1:3)
     end do
  end do

  close(31)

  return

100 format('Iteration: ',i3,';  enthalpy = ',f14.6)
200 format(a2,3(2x,f15.10))

end subroutine print_xyz



