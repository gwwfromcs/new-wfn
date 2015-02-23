!     @process extchk
!
subroutine mixer(iter,errscfconv, pw_params, crys, gs, mixgs, vin, vout, &
     energs, vion, ajac)
  !
  !     1996 by Bernd Pfrommer, while at UCB
  !
  !     Based on a code by JLM.
  !
  use all_to_all_module
  include 'use.h'  
  implicit none         ! implicit? Just say no!
  include 'interface.h' 
  include 'all_to_all.h' 
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: iter                  ! iteration number
  real(dp), intent(in):: errscfconv            ! largest error in the potential
  type(crystal), intent(in) :: crys            ! for vcell, nspin, ztot
  type(pw_parameter), intent(in) :: pw_params  ! for alphamix
  type(parallel_gspace), intent(in) :: &
       gs, &                                   ! the gspace for the potentials
       mixgs                                   ! the (smaller mixing gspace)
  complex(dp), intent(in) :: vout(gs%length * crys%nspin)  
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: energs        ! for vxc0
  complex(dp), intent(inout) :: &
       vion(gs%length * crys%nspin), &         ! the screened ionic potential
       ajac(crys%nspin * mixgs%length,&        ! inverse of jacobian matrix
               crys%nspin * mixgs%length), &           
       vin(gs%length * crys%nspin)             ! the screening potential with
                                               ! which present eigenvectors
                                               ! were calculated
  !     DESCRIPTION:
  !     -----------
  !
  !     computes the new screened potential given
  !
  !     * the input and the  output potential from the previous
  !       iteration.
  !     * the approximate inverse of the hessian matrix
  !
  !     A Broyden Quasi-Newton method is used to estimate the dielectric
  !     matrix for the mixgs%length smallest G-vectors.
  !
  !     A simple linear mixing can also be used by setting
  !     pw_params%alphamix(2)>0.
  !     G-vec components of the potential are intermixed with the previous
  !     step more strongly.
  !
  !     For the spin polarized case, the  assignment of the
  !     inverse Hessian is:
  !
  !
  !           |1 | 2 3 .... nj | nj+1 .... 2*nj-1
  !     ------+--+-------------+-------------------
  !       1   |  |             |
  !     ------+--+-------------+-------------------
  !       2   |  |             |
  !       .   |  |  up,up      |     up,down
  !       .   |  |             |
  !       .   |  |             |
  !      nj   |  |             |
  !     ------+--+-------------+-------------------
  !     nj+1  |  |             |
  !       .   |  |             |
  !       .   |  |  down,up    |     down,down
  !       .   |  |             |
  !     2*nj-1|  |             |
  !
  !
  !     ------------------ local variables ------------------------
  !
  real(dp) :: t0, &
       qtf2, &                             ! square of thomas-fermi vector
       amix, alpha, fac, qf, x, &
       x2, gx, xmix, xmix1, vnorm, ekmmax, &
       vxc0_tmp(2), &   ! temp for energs%vxc0 needed to sum and broadcast
       changemax, changemaxmax   ! max change in mixed potential

  integer :: init, &                  ! flag, what initialization is required
       i, j, k, j0, isub, gsvi(3), isubi, is, it, jt, ind, &
       nchanged, ichangemax   ! number changed components in potential
  complex(dp), allocatable :: av(:), va(:), vinm(:), voutm(:), wmix(:)
  complex(dp) :: vioni, vnew  
  integer, external :: findvec
  real(dp), external :: gimmetime  
  !
  !     ***************************************************************
  !
  t0 = gimmetime()  

  qtf2 = dfour * (dthree * crys%ztot / (pi * crys%vcell))**dthird

  init = 0

  nchanged = 0
  changemax = 0.d0
  vxc0_tmp=0d0  ! init to zero so sum on all proc. contains only new result 

  if (iter == 1) init = 2  
  !
  !     --------------------- use linear mixing --------------------------
  !
  !
  if (pw_params%alphamix(2) > dzero) then  
     do is = 1, crys%nspin  
        do i = 1, gs%length  
           amix = pw_params%alphamix(2) + pw_params%alphamix(3) * &
                exp(-pw_params%alphamix(4) * gs%ekin(i))
           ind = i + (is - 1) * gs%length
           if (gs%ekin(i) > dzero) then  
              vion(ind) = vion(ind) - vin(ind)              ! unscreen
              vin(ind) = amix * vin(ind) + (done - amix) * &
                   vout(ind)                                ! mix new and old
              vion(ind) = vion(ind) + vin(ind)    ! screen with new potential
           else  
              xmix = pw_params%alphamix(2)  
              if (init /= 2) then  
                 energs%vxc0(is) = xmix * energs%vxc0(is) + &
                      (done - xmix) * vin(ind)
              end if
              vin(ind) = energs%vxc0(is)  
           end if
        end do
     end do
     write(9, '('' Mixing new and old potentials with a1='''//',f7.4, &
          &'', a2='',f7.4, '', a3='',f7.4)') pw_params%alphamix(2), &
          pw_params%alphamix(3), pw_params%alphamix(4)

     return  
  end if
  !
  !     ----------------------- Broyden Mixing ---------------------------
  !
  !
  !     allocate work arrays
  !
  allocate(av(mixgs%length * crys%nspin))  
  allocate(va(mixgs%length * crys%nspin))  
  allocate(vinm(mixgs%length * crys%nspin))  
  allocate(voutm(mixgs%length * crys%nspin))  
  allocate(wmix(mixgs%nstar))  
  !
  !     figure out largest mixing energy and index of G=0 component
  !
  ekmmax = dzero
  do i = 1, mixgs%length  
     if (mixgs%ekin(i) == dzero) j0 = i  
     if (mixgs%ekin(i) > ekmmax) ekmmax = mixgs%ekin(i)  
  end do
  alpha = pw_params%alphamix(1)  
  fac = -dsix * alpha / (pi * pi * qtf2)  

  qf = pi * qtf2 * dqtr
  if (mixgs%length > 1) then  
     !
     !     initialize inverse jacobian matrix
     !
     if (init == 2) then  
        do i = 1, mixgs%length * crys%nspin  
           it = mod(i - 1, mixgs%length) + 1  
           do j = 1, mixgs%length * crys%nspin  
              jt = mod(j - 1, mixgs%length) + 1  
              ajac(i, j) = zzero
              if (i == j .and. it /= j0) then  
                 x = sqrt(mixgs%ekin(it)) / qf  
                 x2 = x * dhalf
                 gx = dhalf
                 if (x2 /= done) gx = (done - (done - x2 * x2) * &
                      log(abs((done - x2) / (done + x2))) / x) * dhalf
                 xmix = done / (done + (fac + qtf2 / mixgs%ekin(it)) * gx)
                 ajac(i, j) = cmplx(-xmix, dzero, dp)  
              end if
           end do
        end do
        !         write(9,*) 'initialized Broyden mixing matrix:'
        !         do i=1, mixgs%length*crys%nspin
        !            write(9,'(1000f8.4)') dble(ajac(i,:))
        !         end do
     end if
     !
     !     restore the symmetry V(G)=cc(V(-G))      explicitly
     !
     do is = 1, crys%nspin  
        call regspace(gs, vin(1 + (is - 1) * gs%length), mixgs, &
             vinm(1 + (is - 1) * mixgs%length))
        call regspace(gs, vout(1 + (is - 1) * gs%length), mixgs, &
             voutm(1 + (is - 1) * mixgs%length))
        call symmetrize_local(mixgs%length, &
             vinm(1 + (is - 1) * mixgs%length), &
             mixgs%nstar, mixgs%inds, mixgs%mstar, mixgs%phase, wmix(1), &
             mixgs%inverse)
        call symmetrize_local(mixgs%length, &
             voutm(1 + (is - 1) * mixgs%length), &
             mixgs%nstar, mixgs%inds, mixgs%mstar, mixgs%phase, wmix(1), &
             mixgs%inverse)
     end do
     !
     !     update inverse jacobian matrix a
     !
     if (init == 0) then  
        !         write(9,*) 'voutm     =',dble(voutm)
        !         write(9,*) 'vinm      =',dble(vinm)
        !         write(9,*) 'ajac(j0,j)=',dble(ajac(j0,:))
        !         write(9,*) 'ajac(j,j0)=',dble(ajac(:,j0))
        !         write(9,*) 'all       =',dble(voutm(:)-vinm(:)-
        !     $        ajac(j0,:)+ajac(:,j0))
        do i = 1, mixgs%length * crys%nspin  
           it = mod(i - 1, mixgs%length) + 1  
           if (it /= j0) then  
              av(i) = zzero
              va(i) = zzero
              do j = 1, mixgs%length * crys%nspin  
                 jt = mod(j - 1, mixgs%length) + 1  
                 if (jt /= j0) then  
                    !                              -1
                    !                    av = f = J  (F_m - F_(m-1))
                    !
                    av(i) = av(i) + ajac(i, j) * (voutm(j) - vinm(j) - &
                         ajac(j0, j) + ajac(j, j0))
                    !                          * -1
                    !                    va = g J
                    !
                    va(i) = va(i) + conjg(vinm(j) - ajac(j, j0)) * ajac(j, i)
                 end if
              end do
           end if
        end do
        vnorm = dzero
        do i = 1, mixgs%length * crys%nspin  
           it = mod(i - 1, mixgs%length) + 1  
           if (it /= j0) then  
              vnorm = vnorm + real(conjg(vinm(i) - ajac(i, j0)) * av(i), dp)
           end if
        end do
        !         write(9,*) 'va=',dble(va)
        !         write(9,*) 'av=',dble(av)
        !         write(9,*) 'vinm=',dble(vinm)
        !         write(9,*) 'vnorm=',dble(vnorm)
        !         write(9,*) 'ajac(:,j0)=',dble(ajac(:,j0))
        do i = 1, mixgs%length * crys%nspin  
           it = mod(i - 1, mixgs%length) + 1  
           if (it /= j0) then  
              do j = 1, mixgs%length * crys%nspin  
                 jt = mod(j - 1, mixgs%length) + 1  
                 if (jt /= j0) then  
                    ajac(i, j) = ajac(i, j) + (vinm(i) - ajac(i, j0) - &
                         av(i)) * va(j) / vnorm
                 end if
              end do
           end if
        end do
        !         write(9,*) 'updated Broyden mixing matrix to:'
        !         do i=1, mixgs%length*crys%nspin
        !            write(9,'(1000f8.4)') dble(ajac(i,:))
        !         end do
        ! end of (if init.eq.0)
     end if
     !
     !     save vin and vout for next iteration
     !
     do i = 1, mixgs%length * crys%nspin  
        it = mod(i - 1, mixgs%length) + 1  
        if (it /= j0) then  
           ajac(i, j0) = vinm(i)  
           ajac(j0, i) = voutm(i)  
        end if
     end do
     !
     !     compute correction to vin from inverse jacobian
     !
     do i = 1, mixgs%length * crys%nspin  
        it = mod(i - 1, mixgs%length) + 1  
        if (it /= j0) then  
           av(i) = zzero
           do j = 1, mixgs%length * crys%nspin  
              jt = mod(j - 1, mixgs%length) + 1  
              if (jt /= j0) then  
                 av(i) = av(i) + ajac(i, j) * (voutm(j) - vinm(j))  
              end if
           end do
        end if
     end do
  end if     ! of if(mixgs%length.gt.1)
  !
  !     loop over gspace
  !
  do is = 1, crys%nspin  
     do i = 1, gs%length  
        if (gs%ekin(i) .gt.0.d0) then  
           if (gs%ekin (i) .gt.ekmmax) then     ! do linear mixing for this
              x = sqrt (gs%ekin (i) ) / qf  
              x2 = x / 2.d0  
              gx = 1.d0 / 2.d0  
              if (x2.ne.1.d0) gx = (1.d0 - (1.d0 - x2 * x2) * log (abs ( ( &
                   1.d0 - x2) / (1.d0 + x2) ) ) / x) / 2.d0
              xmix = 1.d0 / (1.d0 + (fac + qtf2 / gs%ekin (i) ) * gx)  
              xmix1 = 1.d0 - xmix  
              vnew = xmix * vout (i + (is - 1) * gs%length) + xmix1 * vin &
                   (i + (is - 1) * gs%length)
           else                                 ! and Broyden mixing here
              isub = findvec(gs%gvec(1, i), mixgs)  
              if (isub <= 0) then  
                 write(9, *) 'mixer: subspace mismatch!', isub  
                 call mystop  
              end if
              vnew = vin(i + (is - 1) * gs%length) - &
                   av(isub + (is - 1) * mixgs%length)
              gsvi(1) = -gs%gvec(1, i)  
              gsvi(2) = -gs%gvec(2, i)  
              gsvi(3) = -gs%gvec(3, i)  
              isubi = findvec(gsvi, mixgs)  
              if (isubi <= 0) then  
                 write(9, *) 'mixer: subspace mismatch!', isubi  
                 call mystop  
              end if
              if (abs(aimag(vinm(isub + (is - 1) * mixgs%length) + &
                   vinm(isubi + (is - 1) * mixgs%length))) > 5.0d-17) &
                   write(9, *) 'mixer-error:', isub + (is - 1) * &
                   mixgs%length, isubi + (is - 1) * mixgs%length, &
                   gs%ekin(i), aimag(vinm(isub + (is - 1) * mixgs%length) + &
                   vinm(isubi + (is - 1) * mixgs%length)), &
                   aimag(voutm(isub + (is - 1) * mixgs%length) + &
                   voutm(isubi + (is - 1) * mixgs%length))
           end if
           !
           !           replace old potential with new potential and compute
           !           the dielectric function.
           !
           vioni = vion(i + (is - 1) * gs%length) - &
                vin(i + (is - 1) * gs%length)
           vion(i + (is - 1) * gs%length) = vioni + vnew 

           if (abs(vin(i+(is-1)*gs%length) - vnew).gt.changemax) then
             changemax = abs(vin(i+(is-1)*gs%length) - vnew)
             ichangemax = i;
           endif
           if (abs(dble(vin(i+(is-1)*gs%length) - vnew))  &
                .gt.pw_params%epscv) then
               nchanged = nchanged + 1
           endif
 
           vin(i + (is - 1) * gs%length) = vnew  
         
        else      ! linear mixing for the G=0 component
           xmix = 0.2d0  
           if (init /= 2) then  
              energs%vxc0(is) = xmix * energs%vxc0(is) + &
                   (done - xmix) * vin(i + (is - 1) * gs%length)
           end if
           vin(i + (is - 1) * gs%length) =energs%vxc0(is)  
           vxc0_tmp(is)= energs%vxc0(is)   
        end if
     end do
  end do

  call all_sum_all(vxc0_tmp,2) 
  energs%vxc0=vxc0_tmp   ! ensure all proc. have same energs%vxc0 DBR

  deallocate(av)
  deallocate(va)  
  deallocate(vinm)
  deallocate(voutm)  
  deallocate(wmix)  

  if (iand(pw_params%output(1), 8) == 8) write(9, 950) gimmetime() - t0

  call all_sum_all(nchanged)
  changemaxmax = changemax
  call all_max_all(changemaxmax)
  write(9,103) nchanged
  write(9,104) changemaxmax
!     ------- print out energy of biggest change -- David Roundy
  write(9,107) gs%ekin(ichangemax)


  return  


103 format(/'Broyden Mixing:'/  &
         '   NUMBER OF CHANGED POTENTIAL COMPONENTS:',i10)
104 format('   LARGEST POTENTIAL CHANGE: ',g12.6)
107 format('            AT AN ENERGY OF: ',g12.6)
100 FORMAT('MIXING MATRIX',//,9(9(1X,2F6.2),/))  

950 format(' TIME FOR MIXING:',f12.3)  

end subroutine mixer
