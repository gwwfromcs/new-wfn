  subroutine getwmask(xatom,nmap,indmtmp,wq,mxdlqp,norb, wmasktmp, &
       xyzmaptmp,mrb2,nref,  crys,ffts,gs,amr,ri,iloc,rcut,delql,d2vkbdq2,&
         nfn)
  include 'use.h'     
  implicit NONE
  include 'flibcalls.ph' 
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs  ! potenital gspace
  type(crystal), intent(in) :: crys      ! crystal structure
  type(fft_struc), intent(in) :: ffts    ! info for the fast fourier transform
  real(dp) xatom(3)              ! atomic coord in lattice vectors * 2pi

  integer mxdlqp, &            ! # of q points for 1d q-space potentials
          norb, &              ! # of orbital shells of projectors for atom
          mrb2, &              ! estimate of max of nmap
          iloc, &              ! the projector used as the local
          nkb, &
          nfn                  ! # bands for simultaneous FFT
  real(dp) wq(mxdlqp,norb),d2vkbdq2(mxdlqp,norb)
  real(dp) ri(201),&           ! radial grid of r/rcut
           amr(201),&          ! mask value on radial grid of r/rcut
           rcut, &             ! cutoff for real space NL potnetial
           delql               ! spacing in qspace for 1d q-space pot.
  !
  !     OUTPUT:
  !     -----
  !
  integer nref, &              ! # of total projecyors over all shells
          nmap, &              ! # grid points within rcut of atom
          indmtmp(mrb2)        ! # of array element for real rep. of 
                               ! WaveFunc correponding to a given wmask value

  real(dp) wmasktmp(9,mrb2), & ! # wmask (ie. pot.*mask) for atom
           xyzmaptmp(3,mrb2)   ! cart. coord. of grid points with rcut

  !     adapted for Paratec July, 2001  by Dave Raczkowski from
  !     PEtot ( written by LW Wang
  !
  !
  ! Does interpolation to 3d grid in Fourier space of the non-local projectors.
  ! Then a 3d fft to real space. The code then searches for grid points within
  ! rcut. The potential*mask function at these points are then stored in wmask.
  ! This algorithm will not work for small systems
  ! where 2*rcut > |any latic vector|. If this is necessary the code must be
  ! changed to do the 1d fft from q-> r and then do the interpolation of the
  ! projectors on a 3d grid in real space. Can look at Escan for the exact
  ! algorithm
  !
  !     --------------------- local variables ---------------------------- 
  !
  complex(dp), allocatable :: &
           workr_n(:),&           ! interpolated 3d potential
            fr(:)                 ! real space interpolated potenital
  complex(dp)            st,flm

  real(dp) x1,y1,z1,x11,y11,z11,r2,r,rr,xt,yt,zt,x10,z10,y10,f1,f2,y, &
           qk(3),qcar(3),qj,qinv,qinv1,qinv2,qinv3,fi,xn,vqj,xni,xa,xb

  integer ir,ierr,k,j,iz_start,iy_start,ncolx,icol_tot,ico,&
          ii1,imap,i,isp,lmin,lmax,n1,n2,n3,iref,n,ni,ni1
  !
  !      --------------------------------------------------------------
  !

  if (nfn .ge. 8) then
    allocate( workr_n(gs%r_size*9), fr(gs%length*9) )
  else
    allocate( workr_n(gs%r_size), fr(gs%length) )
  end if

  n1=ffts%fftsize(1) 
  n2=ffts%fftsize(2)
  n3=ffts%fftsize(3)
  !
  ! get atomic pos. in lat. vec. * grid pts. since xatom incl. 2pi must remove
  !
  x1=xatom(1)*n1/(2*pi)
  y1=xatom(2)*n2/(2*pi)
  z1=xatom(3)*n3/(2*pi)
  !
  !  need to invert the atom coordinates as the structure factor convention
  !  for the potentials invert the original axes
  !
  x1=-x1; y1=-y1; z1=-z1;
  !
  ! shift atomic postion to be in unit cell
  !
  if(x1.lt.0.d0) x1=x1+n1
  if(y1.lt.0.d0) y1=y1+n2
  if(z1.lt.0.d0) z1=z1+n3
  if(x1.gt.n1) x1=x1-n1
  if(y1.gt.n2) y1=y1-n2
  if(z1.gt.n3) z1=z1-n3
     
  x11=crys%avec(1,1)*x1/n1+crys%avec(1,2)*y1/n2+crys%avec(1,3)*z1/n3
  y11=crys%avec(2,1)*x1/n1+crys%avec(2,2)*y1/n2+crys%avec(2,3)*z1/n3
  z11=crys%avec(3,1)*x1/n1+crys%avec(3,2)*y1/n2+crys%avec(3,3)*z1/n3

  nref=0

  if (nfn .ge. 8) then  

    fr=zzero
    workr_n=zzero

    do 10 isp=1,norb

      if (isp .eq. 1) then
        lmin=1; lmax=1;
      else if (isp .eq. 2) then
        lmin=2; lmax=4;
      else if (isp .eq. 3) then
        lmin=5; lmax=9;
      end if

      do iref=lmin,lmax

        if(isp .eq. 1 .and. iloc .eq. 1) goto 10 
        if(isp .eq. 2 .and. iloc .eq. 2) goto 10
        if(isp .eq. 3 .and. iloc .eq. 3) goto 10
        nref=nref+1
        !
        !  INTERPOLATE 1d q-space projector to 3d q-space
        !
        do j=1,gs%length  

          qk(:) = real(gs%gvec(:,j),dp)  

          qcar(1) = crys%bvec(1, 1) * qk(1) + crys%bvec(1, 2) * qk(2) + &
             crys%bvec(1, 3) * qk(3)
          qcar(2) = crys%bvec(2, 1) * qk(1) + crys%bvec(2, 2) * qk(2) + &
             crys%bvec(2, 3) * qk(3)
          qcar(3) = crys%bvec(3, 1) * qk(1) + crys%bvec(3, 2) * qk(2) + &
             crys%bvec(3, 3) * qk(3)

          fi = qcar(1) * x11 + qcar(2) * y11 + qcar(3) * z11
          st = exp(cmplx(dzero, fi, dp))
          qj = sqrt(gs%ekin(j)) 

          xni = qj / delql + dtwo  
          vqj=dzero
          ni = xni  
          if (ni <= 2) ni = 3  
          ni1 = ni + 1  
          if (ni < mxdlqp ) then  
            xa = real(ni1, dp) - xni  
            xb = xni - real(ni, dp)  
            vqj = xa * wq(ni, isp) + xb * wq(ni1, isp) + ((xa**3 - xa) * &
                d2vkbdq2(ni - 2, isp) + (xb**3 - xb) * &
                d2vkbdq2(ni1 - 2, isp)) * dsixth
          end if

          if (iref .eq. 1 ) then
            flm  = zone 
          else
            if (qj .lt. 1d-8) then
              flm  = zzero
            else  
              qinv = done / qj  
              qinv2=qinv*qinv
              qinv3=qinv*qinv*qinv
              if(iref .eq. 2) then  
                flm =zi*qcar(1) * qinv  
              else if(iref .eq. 3) then
                flm =zi*qcar(2) * qinv  
              else if(iref .eq. 4) then
                flm =zi*qcar(3) * qinv  
              else if(iref .eq. 5) then
                flm = drt3 * qcar(1) * qcar(2) * qinv2  
              else if(iref .eq. 6) then
                flm = drt3 * qcar(2) * qcar(3) * qinv2  
              else if(iref .eq. 7) then
                flm = drt3 * qcar(3) * qcar(1) * qinv2  
              else if(iref .eq. 8) then
                flm = (dtrhalf * qcar(3) * qcar(3) * qinv2 - dhalf)
              else if(iref .eq. 9) then
                flm = (drt3 * (qcar(1) * qcar(1) - qcar(2) * qcar(2))*qinv2* &
                dhalf)

              else if(iref.eq.10) then
                flm = qinv3*0.5d0*drt5*oort2* (3.d0*qcar(1)*qcar(1)*qcar(2) - qcar(2)* qcar(2)*qcar(2) )
              else if(iref.eq.11) then
                flm = qinv3*drt15*qcar(1)*qcar(2)*qcar(3)
              else if(iref.eq.12) then
                flm = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)*qinv2-1.d0)*qcar(2)*qinv
              else if(iref.eq.13) then
                flm = 0.5d0*(5.d0*qcar(3)*qcar(3)*qcar(3)*qinv3-3.d0*qcar(3)*qinv)
              else if(iref.eq.14) then
                flm = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)*qinv2-1.d0)*qcar(1)*qinv
              else if(iref.eq.15) then
                flm = qinv3*0.5d0*drt15*(qcar(1)*qcar(1)-qcar(2)*qcar(2))*qcar(3)
              else if(iref.eq.16) then
                flm =-qinv3*0.5d0*drt5*oort2* (3.d0*qcar(2)*qcar(2)*qcar(1) - qcar(1)* qcar(1)*qcar(1) )


              end if
            end if
          end if
  
          vqj=vqj/crys%vcell
          fr(j+(nref-1)*gs%length) = vqj  * st   *flm
        end do   ! end of interpolation
        !
        !  transform to real space
        !

      end do
10  continue

    call fourier_transform(-1, ffts, gs, fr(1), workr_n(1), nref)
    !
    !**** x1,y1,z1 are real number grid indexes of the atom
    !
    ! each PE goes thru grid points that it holds. each has ncolx  x cols. 
    !
    imap=0
    ii1 = 0
    !
    ! determine how many cloumns are on lower number PEs
    !
    icol_tot=dzero
    do i=0,gs%myproc-1
      ncolx=(n2*n3)/ gs%nproc
      if (mod(n2*n3, gs%nproc) > i) ncolx = ncolx + 1
      icol_tot=icol_tot+ncolx
    end do
    !     
    !  determine how many columns are on this PE
    !
    ncolx=(n2*n3)/ gs%nproc
    if (mod(n2*n3, gs%nproc) > gs%myproc) ncolx = ncolx + 1
    !
    ! get starting y and z coordinate from total cols on lower PEs
    !
    iy_start=mod(icol_tot,n2)
    iz_start=icol_tot/n2

    do ico = 1,ncolx
    
      j  = iy_start 
      y10=j-y1
      if(dabs(y10-n2).lt.dabs(y10)) y10=y10-n2
      if(dabs(y10+n2).lt.dabs(y10)) y10=y10+n2

      k=iz_start
      z10=k-z1
      if(dabs(z10-n3).lt.dabs(z10)) z10=z10-n3
      if(dabs(z10+n3).lt.dabs(z10)) z10=z10+n3

      iy_start=iy_start+1
      if(iy_start .ge. n2) then
        iy_start=0
        iz_start=iz_start+1
      end if

      do i=0,n1-1         
        ii1 = ii1 + 1

        x10=i-x1
        if(dabs(x10-n1).lt.dabs(x10)) x10=x10-n1
        if(dabs(x10+n1).lt.dabs(x10)) x10=x10+n1
        !
        !  rr the distance from the atom to the grid point
        !
        xt=crys%avec(1,1)*x10/n1+crys%avec(1,2)*y10/n2+crys%avec(1,3)*z10/n3
        yt=crys%avec(2,1)*x10/n1+crys%avec(2,2)*y10/n2+crys%avec(2,3)*z10/n3
        zt=crys%avec(3,1)*x10/n1+crys%avec(3,2)*y10/n2+crys%avec(3,3)*z10/n3

        rr=xt**2+yt**2+zt**2
        r=dsqrt(rr)

        if(r.lt.rcut-1.D-6) then

          imap=imap+1
          indmtmp(imap)=ii1

          r2=r/rcut

          ir=1+r2*200.d0
          f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
          f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))

          y=amr(ir)*f1+amr(ir+1)*f2

          if(imap.gt.mrb2) then
            write(6,*) "imap > mrb2, stop", imap,mrb2
          endif

          xyzmaptmp(1,imap)=xt
          xyzmaptmp(2,imap)=yt
          xyzmaptmp(3,imap)=zt

          if(nref.gt.0) then     
            do iref=1,nref 
              wmasktmp(iref,imap)=real(workr_n(ii1+(iref-1)*gs%r_size))*y
            end do
          endif

        end if

      end do
    end do

    nmap=imap

    if(nmap.gt.mrb2) then
      write(6,*) "nmap > mrb2, stop", nmap,mrb2
    endif

  else

    do 20 isp=1,norb

      if (isp .eq. 1) then
        lmin=1; lmax=1;
      else if (isp .eq. 2) then
        lmin=2; lmax=4;
      else if (isp .eq. 3) then
        lmin=5; lmax=9;
      end if

      do iref=lmin,lmax

        if(isp .eq. 1 .and. iloc .eq. 1) goto 20 
        if(isp .eq. 2 .and. iloc .eq. 2) goto 20
        if(isp .eq. 3 .and. iloc .eq. 3) goto 20
        nref=nref+1
        !
        !  INTERPOLATE 1d q-space projector to 3d q-space
        !
        fr=dzero
        do j=1,gs%length

          qk(:) = real(gs%gvec(:,j),dp)  

          qcar(1) = crys%bvec(1, 1) * qk(1) + crys%bvec(1, 2) * qk(2) + &
             crys%bvec(1, 3) * qk(3)
          qcar(2) = crys%bvec(2, 1) * qk(1) + crys%bvec(2, 2) * qk(2) + &
             crys%bvec(2, 3) * qk(3)
          qcar(3) = crys%bvec(3, 1) * qk(1) + crys%bvec(3, 2) * qk(2) + &
             crys%bvec(3, 3) * qk(3)

          fi = qcar(1) * x11 + qcar(2) * y11 + qcar(3) * z11
          st = exp(cmplx(dzero, fi, dp))
          qj = sqrt(gs%ekin(j)) 

          xni = qj / delql + dtwo  
          vqj=dzero
          ni = xni  
          if (ni <= 2) ni = 3  
          ni1 = ni + 1  
          if (ni < mxdlqp ) then  
            xa = real(ni1, dp) - xni  
            xb = xni - real(ni, dp)  
            vqj = xa * wq(ni, isp) + xb * wq(ni1, isp) + ((xa**3 - xa) * &
                d2vkbdq2(ni - 2, isp) + (xb**3 - xb) * &
                d2vkbdq2(ni1 - 2, isp)) * dsixth
          end if

          if (iref .eq. 1 ) then
            flm  = zone 
          else
            if (qj .lt. 1d-8) then
              flm  = zzero
            else  
              qinv = done / qj  
              if(iref .eq. 2) then  
                flm =zi*qcar(1) * qinv  
              else if(iref .eq. 3) then
                flm =zi*qcar(2) * qinv  
              else if(iref .eq. 4) then
                flm =zi*qcar(3) * qinv  
              else if(iref .eq. 5) then
                qinv = qinv * qinv  
                flm = drt3 * qcar(1) * qcar(2) * qinv  
              else if(iref .eq. 6) then
                qinv = qinv * qinv  
                flm = drt3 * qcar(2) * qcar(3) * qinv  
              else if(iref .eq. 7) then
                qinv = qinv * qinv 
                flm = drt3 * qcar(3) * qcar(1) * qinv  
              else if(iref .eq. 8) then
                qinv = qinv * qinv 
                flm = (dtrhalf * qcar(3) * qcar(3) * qinv - dhalf)
              else if(iref .eq. 9) then
                qinv = qinv * qinv 
                flm = (drt3*(qcar(1) * qcar(1) - qcar(2) * qcar(2)) * qinv * &
                dhalf)

              else if(iref.eq.10) then
                qinv = qinv * qinv *qinv
                flm = qinv*0.5d0*drt5*oort2* (3.d0*qcar(1)*qcar(1)*qcar(2) - qcar(2)* qcar(2)*qcar(2) )
              else if(iref.eq.11) then
                qinv = qinv * qinv *qinv
                flm = qinv*drt15*qcar(1)*qcar(2)*qcar(3)
              else if(iref.eq.12) then
                qinv = qinv * qinv *qinv
                flm = qinv*0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(2)
              else if(iref.eq.13) then
                qinv = qinv * qinv *qinv
                flm = qinv*0.5d0*(5.d0*qcar(3)*qcar(3)*qcar(3)-3.d0*qcar(3))
              else if(iref.eq.14) then
                qinv = qinv * qinv *qinv
                flm = qinv*0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(1)
              else if(iref.eq.15) then
                qinv = qinv * qinv *qinv
                flm = qinv*0.5d0*drt15*(qcar(1)*qcar(1)-qcar(2)*qcar(2))*qcar(3)
              else if(iref.eq.16) then
                qinv = qinv * qinv *qinv
                flm =-qinv*0.5d0*drt5*oort2* (3.d0*qcar(2)*qcar(2)*qcar(1) - qcar(1)* qcar(1)*qcar(1) )


              end if
            end if
          end if
  
          vqj=vqj/crys%vcell
          fr(j) = vqj  * st   *flm

        end do   ! end of interpolation
        !
        !  transform to real space
        !
        call fourier_transform(-1, ffts, gs, fr(1), workr_n(1), 1)
        !
        !**** x1,y1,z1 are real number grid indexes of the atom
        !
        ! each PE goes thru grid points that it holds. each has ncolx  x cols. 
        !
        imap=0
        ii1 = 0
        !
        ! determine how many cloumns are on lower number PEs
        !
        icol_tot=dzero
        do i=0,gs%myproc-1
          ncolx=(n2*n3)/ gs%nproc
          if (mod(n2*n3, gs%nproc) > i) ncolx = ncolx + 1
          icol_tot=icol_tot+ncolx
        end do
        !     
        !  determine how many columns are on this PE
        !
        ncolx=(n2*n3)/ gs%nproc
        if (mod(n2*n3, gs%nproc) > gs%myproc) ncolx = ncolx + 1
        !
        ! get starting y and z coordinate from total cols on lower PEs
        !
        iy_start=mod(icol_tot,n2)
        iz_start=icol_tot/n2

        do ico = 1,ncolx
      
          j  = iy_start 
          y10=j-y1
          if(dabs(y10-n2).lt.dabs(y10)) y10=y10-n2
          if(dabs(y10+n2).lt.dabs(y10)) y10=y10+n2

          k=iz_start
          z10=k-z1
          if(dabs(z10-n3).lt.dabs(z10)) z10=z10-n3
          if(dabs(z10+n3).lt.dabs(z10)) z10=z10+n3

          iy_start=iy_start+1
          if(iy_start .ge. n2) then
            iy_start=0
            iz_start=iz_start+1
          end if
  
          do i=0,n1-1         
            ii1 = ii1 + 1

            x10=i-x1
            if(dabs(x10-n1).lt.dabs(x10)) x10=x10-n1
            if(dabs(x10+n1).lt.dabs(x10)) x10=x10+n1
            !
            !  rr the distance from the atom to the grid point
            !
            xt=crys%avec(1,1)*x10/n1+crys%avec(1,2)*y10/n2+&
                                     crys%avec(1,3)*z10/n3
            yt=crys%avec(2,1)*x10/n1+crys%avec(2,2)*y10/n2+&
                                     crys%avec(2,3)*z10/n3
            zt=crys%avec(3,1)*x10/n1+crys%avec(3,2)*y10/n2+&
                                     crys%avec(3,3)*z10/n3

            rr=xt**2+yt**2+zt**2
            r=dsqrt(rr)

            if(r.lt.rcut-1.D-6) then

              imap=imap+1
              indmtmp(imap)=ii1

              r2=r/rcut

              ir=1+r2*200.d0
              f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
              f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))

              y=amr(ir)*f1+amr(ir+1)*f2
  
              if(imap.gt.mrb2) then
                write(6,*) "imap > mrb2, stop", imap,mrb2
              endif

              xyzmaptmp(1,imap)=xt
              xyzmaptmp(2,imap)=yt
              xyzmaptmp(3,imap)=zt

              if(nref.gt.0) then     
                wmasktmp(nref,imap)=real(workr_n(ii1))*y
              endif

            end if

          end do
        end do
  
        nmap=imap
  
        if(nmap.gt.mrb2) then
          write(6,*) "nmap > mrb2, stop", nmap,mrb2
        endif
     
      end do

20  continue

  end if

  deallocate(fr)
  deallocate(workr_n)

  return
  end subroutine getwmask
	 

