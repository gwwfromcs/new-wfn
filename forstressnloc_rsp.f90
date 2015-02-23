!     @process extchk
!
subroutine forstressnloc_rsp(ipr, crys, bands, kgs, ekinmod, pspot, &
     rk, aadot, eigenvec, fsum, ssum, pot_gspace,ffts,pw_params)
     
  !
  !     computes the kinetic energy and nonlocal potential contributions
  !     to the Hellman-Feynman forces and the stress.
  !     Adapted from a code by J.L. Martins 1988
  !
  !     parallel version 1995/1996 Bernd Pfrommer
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                   ! never do implicit !
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       kgs(*)                            ! gspaces for the different k-points
  type(parallel_gspace), intent(in) :: pot_gspace
  type(fft_struc), intent(in) :: ffts 
  type(pseudo_potential), intent(in) :: pspot  
  type(band), intent(in) :: bands  
  type(crystal), intent(in) :: crys  
  type(pw_parameter), intent(in) :: pw_params  
  type(complex_gspace_array), intent(in) :: &
       eigenvec                               ! the groundstate wavefunctions
  integer, intent(in) :: ipr                                     ! print flag
  real(dp), intent(in) :: &
       rk(3, bands%nrk), &          ! k-points in units of reciprocal lattice
       ekinmod(3), &             ! parameters for the modified kinetic energy
       aadot(6)                           ! metric in realspace:
                                          ! adot(11,22,33,12,23,31)/(2*pi)**2
  !
  !     INPUT/OUTPUT:  (the nonlocal+kinetic part is added)
  !     -------------
  !
  real(dp), intent(inout) :: &
       fsum(3, crys%mxdatm, crys%ntype), & ! unsymm forces [latt. coordinates]
                                                ! from local part on each atom
       ssum(6)                     ! unsymmetrized, metric stress tensor [Ryd]
  !
  !     ---------------- local variables ---------------------------------
  !     dynamic arrays:
  !     --------------
  complex(dp), allocatable :: anlga(:,:,:),dhd(:,:)   
  real(dp), allocatable :: xnlkb(:,:)
  complex(dp) cphase(pspot%NLPP_rsp_pts),anlra(pot_gspace%r_size, 9) 
  real(dp) anlra_mask(9,pspot%NLPP_rsp_nrefpts)   
  complex(dp) eig_rsp(pot_gspace%r_size)
  !  
  !     other simple variables:
  !     -----------------------
  integer :: i, j, k, l, n, kk, kl, lm, ind, lo, ist, is2, is3, irk, &
       neig, ndim, jmin, jmax, nanl, nanl10, nanlga, nanlga10, &
       ni, lmmin, lmmax, is, g_len, nr, iatom, ico, ico1, ico2, imap,&
       ii1, jj, ir, ico_start, ico1_start, index
  real(dp) :: &
       delta, rkpt(3), qki(3), qk(3), qi, qinv, dqi(6), qcar(3), kcar(3),&
       dqcar(6, 3), forefac, zsumr, fi, xni, &
       vq, dvq, xi, enl, fac_rsp, xt, yt, zt, r, rr, r2, f1, f2, rmask,&
       drmask, damr0, damr1, F_atom(3)
  real(dp) F_NLcart(3, crys%mxdatm, crys%ntype), s_NLcart(6), s_enl(6),temp(6)
 
  logical :: isp  
  complex(dp)  st, flm(16), dflm(6, 16), y
  real(dp), parameter :: eps = 1.0d-12 
  real(dp) amr(201),ri(201)

  amr(:)=pspot%amr(:)
  ri(:)=pspot%ri(:)
  F_NLcart=dzero
  s_NLcart=dzero
  s_enl=dzero

  isp = .false.  
  if (crys%nspin == 2) isp = .true.  

  forefac = dtwo / (sqrt(pi) * ekinmod(3)) * ekinmod(1)  
  !     ------- first allocate dynamic arrays ------------
  g_len=pot_gspace%length
  !
  !     count the max number of nonlocal components nanlga for array dimen
  !
  nanlga = 0  
  do k = 1, crys%ntype  
     ind = 0  
     do l = 1, 5  
        if (pspot%nkb(l, k) /= 0) then  
           lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
           lmmax = (pspot%lo(l, k) + 1) *(pspot%lo(l, k) + 1)  
           do lm = lmmin, lmmax  
              ind = ind + 1  
           end do
        end if
     end do
     nanlga = max(nanlga, ind)
  end do

  nanlga10 = 10 * nanlga
  allocate(anlga(g_len, 9, nanlga)) ; anlga = zzero
  write(9,*) pot_gspace%r_size,nanlga,pspot%NLPP_rsp_nrefpts
  anlra = zzero
  allocate(dhd(10, nanlga)) ; dhd = zzero

  if (pspot%nanl > 0) allocate(xnlkb(nanlga,crys%ntype))  
  do is = 1, crys%nspin  
     do irk = 1, bands%nrk                           ! loop over all k-points
        neig = bands%nband(irk, is)  
        ndim = kgs(irk)%length  
        rkpt(1) = rk(1, irk)  
        rkpt(2) = rk(2, irk)  
        rkpt(3) = rk(3, irk)  
        !        ----------- kinetic energy contribution to the stress ---------
        !        find min and max occupied band
        jmin = 0  
        jmax = 0  
        do i = 1, neig  
           if (bands%occup(i, irk, is) /= dzero .and. jmin == 0) jmin = i  
           if (bands%occup(i, irk, is) /= dzero .and. jmin /= 0) jmax = i  
        end do
        !
        !      start loop over g(i)
        !
        do i = 1, ndim  
           !
           !           compute |g+k|
           !
           qki(1) = rkpt(1) + real(kgs(irk)%gvec(1, i), dp)  
           qki(2) = rkpt(2) + real(kgs(irk)%gvec(2, i), dp)  
           qki(3) = rkpt(3) + real(kgs(irk)%gvec(3, i), dp)  
           !
           !           do eigenvector multiplication and band sum
           !
           zsumr = dzero  
           do k = jmin, jmax  
              zsumr = zsumr + bands%occup(k, irk, is) / &
                   real(crys%nspin, dp) * &
                   (real(eigenvec%data((k - 1) * ndim + i, irk, is)) * &
                   real(eigenvec%data((k - 1) * ndim + i, irk, is)) + &
                   aimag(eigenvec%data((k - 1) * ndim + i, irk, is)) * &
                   aimag(eigenvec%data((k - 1) * ndim + i, irk, is)))
           end do
           !
           !           add the kinetic energy functional modifications
           !
           delta = (kgs(irk)%ekin(i) - ekinmod(2)) / ekinmod(3)  

           zsumr = zsumr * (done + forefac * exp(-delta * delta))  

           zsumr = dtwo * zsumr  
           ssum(1) = ssum(1) + zsumr * qki(1) * qki(1)  
           ssum(2) = ssum(2) + zsumr * qki(2) * qki(2)  
           ssum(3) = ssum(3) + zsumr * qki(3) * qki(3)  
           ssum(4) = ssum(4) + zsumr * qki(1) * qki(2)  
           ssum(5) = ssum(5) + zsumr * qki(2) * qki(3)  
           ssum(6) = ssum(6) + zsumr * qki(3) * qki(1)  
        end do                                     ! end of loop over k-points
     end do                                           ! end of loop over spins
  end do
  !     ----------- nonlocal potential contribution to forces and stress -
  !
  !     Consolidate across processors by summing. THIS CAN ONLY BE DONE ON
  !
  call all_sum_all(fsum, 3 * crys%mxdatm * crys%ntype)  
  call all_sum_all(ssum, 6)  

  nr=ffts%fftsize(1)*ffts%fftsize(2)*ffts%fftsize(3)
  fac_rsp=crys%vcell/real(nr,dp)

  iatom=0
  ico = 0
  ico1=0
  do k = 1, crys%ntype  
    do kk = 1, crys%natom(k)  
      iatom=iatom+1
      nanl=pspot%numref(iatom)
      do i = 1, g_len  
        qk(1) = real(pot_gspace%gvec(1, i), dp)  
        qk(2) = real(pot_gspace%gvec(2, i), dp)  
        qk(3) = real(pot_gspace%gvec(3, i), dp)  

        qi=sqrt(pot_gspace%ekin(i))

        if(qi > eps) then

         qinv = done / qi  
         !
         !         derivatives of qi with respect to strain
         !         dq/de(a,b) = -q(a)*q(b)/q^2
         !
         !
         !         qvector in cartesian coordinates and normalized
         !
         kcar(1)=crys%bvec(1,1)*qk(1)+crys%bvec(1,2)*qk(2)+crys%bvec(1,3)*qk(3)
         kcar(2)=crys%bvec(2,1)*qk(1)+crys%bvec(2,2)*qk(2)+crys%bvec(2,3)*qk(3)
         kcar(3)=crys%bvec(3,1)*qk(1)+crys%bvec(3,2)*qk(2)+crys%bvec(3,3)*qk(3)
         qcar(1) = kcar(1) * qinv  
         qcar(2) = kcar(2) * qinv  
         qcar(3) = kcar(3) * qinv  


         dqi(1) = -kcar(1) * kcar(1) * qinv  
         dqi(2) = -kcar(2) * kcar(2) * qinv  
         dqi(3) = -kcar(3) * kcar(3) * qinv  
         dqi(4) = -kcar(1) * kcar(2) * qinv  
         dqi(5) = -kcar(2) * kcar(3) * qinv  
         dqi(6) = -kcar(3) * kcar(1) * qinv  
         !
         !         derivatives of qcar with respect to strain
         !         so d(kcar/k) / de  = kcar * (d(1/k) / dk) * dk/de
         !          + 1/k*d(kcar)/de
         !         =  -kcar/k^2 * kcar(a)*kcar(b)/k + kcar(b)/k 
         !
         dqcar(1, 1) = - kcar(1)   
         dqcar(2, 1) = dzero    
         dqcar(3, 1) = dzero   
         dqcar(4, 1) = - kcar(2)   
         dqcar(5, 1) = dzero    
         dqcar(6, 1) = - kcar(3)   
         dqcar(1, 2) = dzero    
         dqcar(2, 2) = - kcar(2)   
         dqcar(3, 2) = dzero   
         dqcar(4, 2) = - kcar(1)   
         dqcar(5, 2) = - kcar(3)   
         dqcar(6, 2) = dzero   
         dqcar(1, 3) = dzero   
         dqcar(2, 3) = dzero   
         dqcar(3, 3) = - kcar(3)   
         dqcar(4, 3) = dzero    
         dqcar(5, 3) = - kcar(2)   
         dqcar(6, 3) = - kcar(1)  
         !
         !
         do kl = 1, 6
            dqcar(kl, 1) = (dqcar(kl, 1) - dqi(kl) * qcar(1)) * qinv
            dqcar(kl, 2) = (dqcar(kl, 2) - dqi(kl) * qcar(2)) * qinv
            dqcar(kl, 3) = (dqcar(kl, 3) - dqi(kl) * qcar(3)) * qinv
         end do
         !
         !         angular functions
         !
         flm(1) = done  
         flm(2) = cmplx(dzero,qcar(1),dp)  
         flm(3) = cmplx(dzero,qcar(2),dp)  
         flm(4) = cmplx(dzero,qcar(3),dp)  
         flm(5) = drt3 * qcar(1) * qcar(2)  
         flm(6) = drt3 * qcar(2) * qcar(3)  
         flm(7) = drt3 * qcar(3) * qcar(1)  
         flm(8) = dthree * qcar(3) * qcar(3) * dhalf - dhalf
         flm(9) = drt3 * (qcar(1) * qcar(1) - qcar(2) * qcar(2)) * dhalf

                 flm(10) = 0.5d0*drt5*oort2* (3.d0*qcar(1)*qcar(1)*qcar(2) - qcar(2)* qcar(2)*qcar(2) )
                 flm(11) = drt15*qcar(1)*qcar(2)*qcar(3)
                 flm(12) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(2)
                 flm(13) = 0.5d0*(5.d0*qcar(3)*qcar(3)*qcar(3)-3.d0*qcar(3))
                 flm(14) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(1)
                 flm(15) = 0.5d0*drt15*(qcar(1)*qcar(1)-qcar(2)*qcar(2))*qcar(3)
                 flm(16) =-0.5d0*drt5*oort2* (3.d0*qcar(2)*qcar(2)*qcar(1) - qcar(1)* qcar(1)*qcar(1) )


         !
         !                    strain derivatives of angular functions
         !                    in lattice coordinates
         !
         do kl = 1, 6  
            dflm(kl, 1) = zzero  
            dflm(kl, 2) = dqcar(kl, 1)  
            dflm(kl, 3) = dqcar(kl, 2)  
            dflm(kl, 4) = dqcar(kl, 3)  
            dflm(kl, 5) = drt3*(dqcar(kl,1) * qcar(2) + qcar(1) * dqcar(kl,2))
            dflm(kl, 6) = drt3*(dqcar(kl,2) * qcar(3) + qcar(2) * dqcar(kl,3))
            dflm(kl, 7) = drt3*(dqcar(kl,3) * qcar(1) + qcar(3) * dqcar(kl,1))
            dflm(kl, 8) = dthree * dqcar(kl, 3) * qcar(3)  
            dflm(kl, 9) = drt3*(dqcar(kl,1) * qcar(1) - dqcar(kl,2) * qcar(2))

                 dflm(kl,10) = 0.5d0*drt5*oort2*  &
                              (6.d0*dqcar(kl,1)*qcar(1)*qcar(2)+3.d0*qcar(1)*qcar(1)*dqcar(kl,2) &
                               - 3.d0*qcar(2)* qcar(2)*dqcar(kl,2) )

                 dflm(kl,11) = drt15*(dqcar(kl,1)*qcar(2)*qcar(3)+ &
                                      dqcar(kl,2)*qcar(1)*qcar(3)+ &
                                      dqcar(kl,3)*qcar(1)*qcar(2))

                 dflm(kl,12) = 0.5d0*drt3*oort2* &
                               (10.d0*dqcar(kl,3)*qcar(3)*qcar(2)+ &
                                (5.d0*qcar(3)*qcar(3)-1.d0)*dqcar(kl,2))

                 dflm(kl,13) = 0.5d0*(15.d0*qcar(3)*qcar(3)-3.d0)*dqcar(kl,3)

                 dflm(kl,14) = 0.5d0*drt3*oort2* &
                              ( (5.d0*qcar(3)*qcar(3)-1.d0)*dqcar(kl,1)+ &
                                10.d0*qcar(3)*dqcar(kl,3)*qcar(1) )

                 dflm(kl,15) = 0.5d0*drt15* &
                               ( (2.d0*dqcar(kl,1)*qcar(1)-2.d0*dqcar(kl,2)*qcar(2))*qcar(3) + &
                                 (qcar(1)*qcar(1)-qcar(2)*qcar(2))*dqcar(kl,3) )

                 dflm(kl,16) =-0.5d0*drt5*oort2* &
                              ( 6.d0*dqcar(kl,2)*qcar(2)*qcar(1) + 3.d0*qcar(2)*qcar(2)*dqcar(kl,1) - &
                                3.d0*dqcar(kl,1)*qcar(1)*qcar(1) )

         end do
        else
        !
        !  treat G=0 component
        !
          qi = dzero; dqi(:)=dzero; dflm(:,:)=dzero; flm(:)=dzero;
          flm(1)=done
        end if 
        !
        !   structure factor
        !
        fi = qk(1) * crys%rat(1, kk, k) + qk(2) * crys%rat(2, kk, k) + &
             qk(3) * crys%rat(3, kk, k)
        st = exp(cmplx(dzero, -fi, dp))  
        !
        ind = 0  
        do l = 1, 5  
           if (pspot%nkb(l, k) /= 0) then  
              xni = qi / pspot%delqnl(k) + dtwo  
              ni = int(xni + dhalf)  
              if (ni <= 3) ni = 4  
              vq = dzero  
              dvq = dzero  
              if (ni < pspot%nqnl(k)) then  
                 xi = xni - real(ni, dp)  
                 vq = pspot%vkb_mask(ni, l, k) * (done + xi) * (done - xi) + &
                      dhalf * (pspot%vkb_mask(ni + 1, l, k) * (done + xi) - &
                      pspot%vkb_mask(ni - 1, l, k) * (done - xi)) * xi
                 dvq = dmtwo * xi * pspot%vkb_mask(ni, l, k) + &
                      (xi + dhalf) * pspot%vkb_mask(ni + 1, l, k) + &
                      (xi - dhalf) * pspot%vkb_mask(ni - 1, l, k)
                 dvq = dvq / pspot%delqnl(k)  
              end if
              lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  

              lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)
              !
              do lm = lmmin, lmmax  
                 ind = ind + 1  
                 !
                 !                          forces
                 !
                 anlga(i,1,ind) = -cmplx(dzero, dmone, dp) * &
                      st * flm(lm) * vq  * kcar(1)/crys%vcell
                 anlga(i,2,ind) = -cmplx(dzero, dmone, dp) * &
                      st * flm(lm) * vq * kcar(2)/crys%vcell
                 anlga(i,3,ind) = -cmplx(dzero, dmone, dp) * &
                      st * flm(lm) * vq * kcar(3)/crys%vcell
                 !
                 !                          stresses
                 !
                 do ist = 1, 6
                    anlga(i,3 + ist, ind) = st*  ( dflm(ist, lm) * vq   + &
                     flm(lm)   * dvq * dqi(ist)   )/crys%vcell
                 end do
              end do
           end if
           !
        end do                     ! loop over lo
        !
 
      end do                        ! loop over gspace
      ! 
      !  transform pjoectors and derivatives to real space
      !
      ico_start=ico
      ico1_start=ico1
      do ind=1,nanl
        ico1=ico1_start
        do ist=1,9
          call fourier_transform(-1,ffts,pot_gspace,anlga(1,ist,ind), &
                                                     anlra(1,ist),1 )
        end do
        ! 
        !  pack in to place : anlga -> anlra ;  and apply mask function
        ! 
        do imap=1,pspot%nmap(iatom)
          ico1=ico1+1
       
          xt=pspot%xyzmap(3*ico1-2)      
          yt=pspot%xyzmap(3*ico1-1) 
          zt=pspot%xyzmap(3*ico1)

          rr=xt**2+yt**2+zt**2
          r=dsqrt(rr)
          r2=r/pw_params%NLPP_rcut(k)

          ir=1+r2*200.d0
          f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
          f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))

          rmask=amr(ir)*f1+amr(ir+1)*f2

          if(ir.eq.1) ir=ir+1
          if(ir.ge.200) ir=200-1
          damr0=(amr(ir+1)-amr(ir-1))/(ri(ir+1)-ri(ir-1))
          damr1=(amr(ir+2)-amr(ir))/(ri(ir+2)-ri(ir))
          drmask=(damr0*f1+damr1*f2)/pw_params%NLPP_rcut(k)
          if(r.lt.1.E-9) then
            drmask=dzero
          endif
 
          ii1=pspot%indm(ico1)
          index=ico_start+(imap-1)*nanl+ind

          if(r.lt.1.E-9) then
            drmask=dzero
          else
            drmask=drmask/r
          end if

          anlra_mask(1,index)=&
            -pspot%kb_rsp(index)*drmask*xt + real(anlra(ii1,1),dp)*rmask
          anlra_mask(2,index)=&
            -pspot%kb_rsp(index)*drmask*yt + real(anlra(ii1,2),dp)*rmask
          anlra_mask(3,index)=&
            -pspot%kb_rsp(index)*drmask*zt + real(anlra(ii1,3),dp)*rmask
          anlra_mask(4,index)=pspot%kb_rsp(index)&
                 *drmask*xt*xt+ &
             real(anlra(ii1,4),dp)*rmask
          anlra_mask(5,index)=pspot%kb_rsp(index)&
                 *drmask*yt*yt + &
          real(anlra(ii1,5),dp)*rmask
          anlra_mask(6,index)=pspot%kb_rsp(index)&
                    *drmask*zt*zt + & 
              real(anlra(ii1,6),dp)*rmask
          anlra_mask(7,index)=pspot%kb_rsp(index)&
                  *drmask*xt*yt + & 
          real(anlra(ii1,7),dp)*rmask
          anlra_mask(8,index)=pspot%kb_rsp(index)&
                   *drmask*yt*zt + &
            real(anlra(ii1,8),dp)*rmask
          anlra_mask(9,index)=pspot%kb_rsp(index)&
                    *drmask*zt*xt + &
              real(anlra(ii1,9),dp)*rmask

        end do  !nmap
       end do ! nanl
      ico=ico+pspot%nmap(iatom)*nanl
    end do  ! # atoms of ntype
  end do    ! ntypes
  !
  !   calculates the matrix elements
  !
  do irk = 1, bands%nrk                               ! loop over all k-points

    ndim = kgs(irk)%length  
    rkpt(1) = rk(1, irk)  
    rkpt(2) = rk(2, irk)  
    rkpt(3) = rk(3, irk)  
  
    kcar(1) = crys%bvec(1, 1) * rkpt(1) + crys%bvec(1, 2) * rkpt(2) + &
              crys%bvec(1, 3) * rkpt(3)
    kcar(2) = crys%bvec(2, 1) * rkpt(1) + crys%bvec(2, 2) * rkpt(2) + &
             crys%bvec(2, 3) * rkpt(3)
    kcar(3) = crys%bvec(3, 1) * rkpt(1) + crys%bvec(3, 2) * rkpt(2) + &
             crys%bvec(3, 3) * rkpt(3)    
    !
    !  setup cphase
    !
    iatom=0
    ico1=0
    do k = 1, crys%ntype  
      do kk = 1, crys%natom(k)  
        iatom=iatom+1
        do i=1,pspot%nmap(iatom)
          ico1=ico1+1
          cphase(ico1)=exp(-zi*(pspot%xyzmap(ico1*3-2)*kcar(1)+ &
             pspot%xyzmap(ico1*3-1)*kcar(2)+pspot%xyzmap(ico1*3)*kcar(3))) 
        enddo
      end do
    end do
    !
    !  calculate the matrix vector product
    !
    do k = 1, crys%ntype  
      ind = 0 
      do l = 1, 5  
        if (pspot%nkb(l, k) /= 0) then  
          lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
          lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)  
          do lm = lmmin, lmmax  
            ind = ind + 1  
            xnlkb(ind,k) = real(pspot%nkb(l, k), dp)  
          end do
        end if
      end do
    end do

    do is = 1, crys%nspin  

      neig = bands%nband(irk, is)  

      do n=1,neig

        call fourier_transform(-1,ffts,kgs(irk),&
                            eigenvec%data((n-1)*ndim+1,irk,is), eig_rsp(1),1) 

        call mzdscal(pot_gspace%r_size,done/sqrt(crys%vcell),eig_rsp(1),1)
 
        ico1=0
        ico2=1
        iatom=1

        do k = 1, crys%ntype  
      
          nanl=pspot%numref(iatom)
         
          do kk = 1, crys%natom(k)        

            dhd=zzero
            do i=1,pspot%nmap(iatom)
              ico1=ico1+1
              y=eig_rsp(pspot%indm(ico1))*cphase(ico1)

              select case(nanl)
  
              case(8)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                dhd(1,3)=dhd(1,3) + pspot%kb_rsp(ico2+2) * y
                dhd(1,4)=dhd(1,4) + pspot%kb_rsp(ico2+3) * y
                dhd(1,5)=dhd(1,5) + pspot%kb_rsp(ico2+4) * y
                dhd(1,6)=dhd(1,6) + pspot%kb_rsp(ico2+5) * y
                dhd(1,7)=dhd(1,7) + pspot%kb_rsp(ico2+6) * y
                dhd(1,8)=dhd(1,8) + pspot%kb_rsp(ico2+7) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                  dhd(1+ist,3)=dhd(1+ist,3)+anlra_mask(ist,ico2+2)* y
                  dhd(1+ist,4)=dhd(1+ist,4)+anlra_mask(ist,ico2+3)* y
                  dhd(1+ist,5)=dhd(1+ist,5)+anlra_mask(ist,ico2+4)* y
                  dhd(1+ist,6)=dhd(1+ist,6)+anlra_mask(ist,ico2+5)* y
                  dhd(1+ist,7)=dhd(1+ist,7)+anlra_mask(ist,ico2+6)* y
                  dhd(1+ist,8)=dhd(1+ist,8)+anlra_mask(ist,ico2+7)* y
                enddo
                ico2=ico2+8
              case(6)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                dhd(1,3)=dhd(1,3) + pspot%kb_rsp(ico2+2) * y
                dhd(1,4)=dhd(1,4) + pspot%kb_rsp(ico2+3) * y
                dhd(1,5)=dhd(1,5) + pspot%kb_rsp(ico2+4) * y
                dhd(1,6)=dhd(1,6) + pspot%kb_rsp(ico2+5) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                  dhd(1+ist,3)=dhd(1+ist,3)+anlra_mask(ist,ico2+2)* y
                  dhd(1+ist,4)=dhd(1+ist,4)+anlra_mask(ist,ico2+3)* y
                  dhd(1+ist,5)=dhd(1+ist,5)+anlra_mask(ist,ico2+4)* y
                  dhd(1+ist,6)=dhd(1+ist,6)+anlra_mask(ist,ico2+5)* y
                enddo
                ico2=ico2+6
                ico2=ico2+8
              case(5)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                dhd(1,3)=dhd(1,3) + pspot%kb_rsp(ico2+2) * y
                dhd(1,4)=dhd(1,4) + pspot%kb_rsp(ico2+3) * y
                dhd(1,5)=dhd(1,5) + pspot%kb_rsp(ico2+4) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                  dhd(1+ist,3)=dhd(1+ist,3)+anlra_mask(ist,ico2+2)* y
                  dhd(1+ist,4)=dhd(1+ist,4)+anlra_mask(ist,ico2+3)* y
                  dhd(1+ist,5)=dhd(1+ist,5)+anlra_mask(ist,ico2+4)* y
                enddo
                ico2=ico2+5
              case(4)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                dhd(1,3)=dhd(1,3) + pspot%kb_rsp(ico2+2) * y
                dhd(1,4)=dhd(1,4) + pspot%kb_rsp(ico2+3) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                  dhd(1+ist,3)=dhd(1+ist,3)+anlra_mask(ist,ico2+2)* y
                  dhd(1+ist,4)=dhd(1+ist,4)+anlra_mask(ist,ico2+3)* y
                enddo
                ico2=ico2+4
              case(3)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                dhd(1,3)=dhd(1,3) + pspot%kb_rsp(ico2+2) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                  dhd(1+ist,3)=dhd(1+ist,3)+anlra_mask(ist,ico2+2)* y
                enddo
                ico2=ico2+3
              case(2)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                dhd(1,2)=dhd(1,2) + pspot%kb_rsp(ico2+1) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                  dhd(1+ist,2)=dhd(1+ist,2)+anlra_mask(ist,ico2+1)* y
                enddo
                ico2=ico2+2
              case(1)
                dhd(1,1)=dhd(1,1) + pspot%kb_rsp(ico2) * y
                do ist=1,9
                  dhd(1+ist,1)=dhd(1+ist,1)+anlra_mask(ist,ico2)* y
                enddo
                ico2=ico2+1
              case default
                write(0, *) 'ERROR # of projectors is wrong', nanl
                call mystop
              end select

            enddo

            call all_sum_all(dhd, 10*nanlga )  
            call mzdscal(10*nanlga,fac_rsp,dhd(1,1),1)

            F_atom=dzero       

       select case(nanl)

         case(8)
              do ist = 1, 3  
                F_atom(ist)=F_atom(ist)+&
                  xnlkb(1,k)*real(dhd(1,1)*conjg(dhd(1+ist,1)),dp)+&
                  xnlkb(2,k)*real(dhd(1,2)*conjg(dhd(1+ist,2)),dp)+&
                  xnlkb(3,k)*real(dhd(1,3)*conjg(dhd(1+ist,3)),dp)+&
                  xnlkb(4,k)*real(dhd(1,4)*conjg(dhd(1+ist,4)),dp)+&
                  xnlkb(5,k)*real(dhd(1,5)*conjg(dhd(1+ist,5)),dp)+&
                  xnlkb(6,k)*real(dhd(1,6)*conjg(dhd(1+ist,6)),dp)+&
                  xnlkb(7,k)*real(dhd(1,7)*conjg(dhd(1+ist,7)),dp)+&
                  xnlkb(8,k)*real(dhd(1,8)*conjg(dhd(1+ist,8)),dp)
              end do 

         case default
          write(0, *) 'ERROR # of projectors is wrong', nanl
          call mystop
         end select

      
            do ist = 1, 3  
              F_NLcart(ist, kk, k) = F_NLcart(ist, kk, k) + dtwo * &
                     bands%occup(n, irk, is) / real(crys%nspin, dp)*F_atom(ist)
   
            end do


       select case(nanl)

         case(8)

              enl = bands%occup(n, irk, is) / real(crys%nspin, dp) * ( &
                    xnlkb(1,k) * real(dhd(1,1)*conjg(dhd(1,1)),dp)+& 
                  xnlkb(2,k) * real(dhd(1,2)*conjg(dhd(1,2)),dp)+& 
                    xnlkb(3,k) * real(dhd(1,3)*conjg(dhd(1,3)),dp)+& 
                  xnlkb(4,k) * real(dhd(1,4)*conjg(dhd(1,4)),dp)+& 
                    xnlkb(5,k) * real(dhd(1,5)*conjg(dhd(1,5)),dp)+& 
                  xnlkb(6,k) * real(dhd(1,6)*conjg(dhd(1,6)),dp)+& 
                    xnlkb(7,k) * real(dhd(1,7)*conjg(dhd(1,7)),dp)+& 
                    xnlkb(8,k) * real(dhd(1,8)*conjg(dhd(1,8)),dp) )
 
              do ist = 1, 6  
                s_enl(ist) =  s_enl(ist) + enl * aadot(ist)
                s_NLcart(ist) = s_NLcart(ist) - dtwo * &
                            bands%occup(n, irk, is) / real(crys%nspin, dp) *( &
                           xnlkb(1,k)*real(dhd(4+ist,1)*conjg(dhd(1, 1)))+& 
                           xnlkb(2,k)*real(dhd(4+ist,2)*conjg(dhd(1, 2)))+& 
                           xnlkb(3,k)*real(dhd(4+ist,3)*conjg(dhd(1, 3)))+& 
                           xnlkb(4,k)*real(dhd(4+ist,4)*conjg(dhd(1, 4)))+& 
                           xnlkb(5,k)*real(dhd(4+ist,5)*conjg(dhd(1, 5)))+& 
                           xnlkb(6,k)*real(dhd(4+ist,6)*conjg(dhd(1, 6)))+& 
                           xnlkb(7,k)*real(dhd(4+ist,7)*conjg(dhd(1, 7)))+&
                           xnlkb(8,k)*real(dhd(4+ist,8)*conjg(dhd(1, 8))) )
              end do


         case default
          write(0, *) 'ERROR # of projectors is wrong', nanl
          call mystop
         end select

            iatom=iatom+1  

          end do   ! # natoms of type

        end do    ! # ntypes
      end do  ! # of bands
 
    end do               ! end of loop over spins
  end do                        ! end of loop over k-points

  !     free the dynamic arrays
  deallocate(anlga)  
  deallocate(dhd) 

  do k = 1, crys%ntype  
    do kk = 1, crys%natom(k)  
      do i = 1, 3  
           fsum(i,kk,k) = fsum(i,kk,k) +  (crys%avec(1,i)*F_NLcart(1, kk, k) + &
                crys%avec(2, i) * F_NLcart(2, kk, k) + &
                crys%avec(3, i) * F_NLcart(3, kk, k) )/(2*pi)
      end do
    end do
  end do

  temp=dzero

  do i = 1, 6    
     j = i  
     k = i  
     if (i > 3) j = i - 3  
     if (i > 3) k = j + 1  
     if (k > 3) k = 1  
     temp(i) = temp(i) + (crys%avec(1, j) * s_NLcart(1) * crys%avec(1, k) + &
          crys%avec(2, j) * s_NLcart(2) * crys%avec(2, k) + &
          crys%avec(3, j) * s_NLcart(3) * crys%avec(3, k) + &
          crys%avec(1, j) * s_NLcart(4) * crys%avec(2, k) + &
          crys%avec(2, j) * s_NLcart(5) * crys%avec(3, k) + &
          crys%avec(3, j) * s_NLcart(6) * crys%avec(1, k) + &
          crys%avec(1, k) * s_NLcart(4) * crys%avec(2, j) + &
          crys%avec(2, k) * s_NLcart(5) * crys%avec(3, j) + &
          crys%avec(3, k) * s_NLcart(6) * crys%avec(1, j) )/((2*pi)**2)
  end do

  ssum=ssum+temp

  ssum=ssum+s_enl
  
  if (pspot%nanl > 0) deallocate (xnlkb)  

  return  

end subroutine forstressnloc_rsp
