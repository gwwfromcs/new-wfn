subroutine w_linerq(it,nt,zatom,vql0,vda,denc_atom,rhoq,vnl_atom,&
         norb,pspot,pw_params)
  !     ------------------------------------------------------------------
  use all_to_all_module  
  include 'use.h'  
  implicit none               ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'
  !
  !     INPUT:
  !     -----
  !     
  type(pw_parameter), intent(inout) :: pw_params 
  integer, intent(in) :: it,nt
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(pseudo_potential), intent(inout) :: pspot
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out)  ::  norb          ! # of orbital for pp_format=3
  real(dp), intent(out) ::  zatom,&       ! charge of ion
             vql0,&
      vda(pspot%mxdlqp),&                 ! local potential in q-space 
      rhoq(pspot%mxdlqp),&                 ! charge in q-space
      denc_atom(pspot%mxdlqp),&           ! core charge in q-space
      vnl_atom(pspot%mxdlqp)              ! nl avergae in q-space
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !    read in and convert real space pseudo to paratec fourier space
  !
  !     ------------ local variables ---------------------------------
  !     

! PZ
  integer, parameter :: max_nrr=3000
  real(dp) vloc(max_nrr),&    ! local potential on real space grid
           vlocT(max_nrr),&   ! local occupation-average of NL pot.
           rhoc(max_nrr),&    ! core charge
           rhoc1d(max_nrr),&  ! 1 derivative of core charge. for FHI format
           rhoc2d(max_nrr),&  ! 2 derivative of core charge. for FHI format
           vw(max_nrr),&      ! workspace for projector
           vw2(max_nrr)       ! workspace for projector*mask

  real(dp) ws(max_nrr),wp(max_nrr),wd(max_nrr),&  ! atomic reference WF
           vs(max_nrr),vp(max_nrr),vd(max_nrr)    ! atomic PPs

  real(dp) r(pspot%mxdlqp),&  ! radial mesh in real spaace
           occ(3),&           ! atomic occupations for the WF if  pp_format=3
           rcut,&             ! rcut if real space non-local used
           delql              ! spacing for 1-space 1d mesh


  integer iloc,&                 ! local pseudopotenital number if pp_format=3
           nrr                   ! # of radial points

  integer  iiatom,idum,isp,ic,i,iq,ierr,ir,nproj 
  real(dp) rhoi1,c11,b22 ,c1,b11 ,rhoi ,scale,fac_s,fac_p,fac_d ,aj2 ,&
            x ,aj0,aj1,ch,s ,z ,rst,x1,x2,r2,b2,b1,s1,g1, s2, r1,rst1,&
            s3,a,b,amrI,f1,f2,dum(4)
  !   
  !     ------------------------------------------------------------------
  !   
  rhoc=dzero
  rhoc1d=dzero
  rhoc2d=dzero

  delql=pspot%delqnl(nt)
  rcut=pw_params%NLPP_rcut(nt)
  !
  ! read in Psi/r and PPs for format for PEtot
  !
  if (pw_params%pp_format .eq. 2 ) then   

    read(it,*) nrr,ic,iiatom,zatom,iloc,occ(1), occ(2),occ(3)
    read(it,*)   pspot%nkb(1,nt),pspot%nkb(2,nt),pspot%nkb(3,nt)

    pw_params%pp_iloc(nt) = iloc

    if (pspot%nkb(3,nt) .eq. 0) then
      norb=2
    else if (pspot%nkb(2,nt) .eq. 0 .and. pspot%nkb(3,nt) .eq. 0) then
      norb=1
    else
      norb = 3
    end if

    if(nrr.gt.max_nrr) then
      write(6,*) "nrr > 2000, stop"
      call mystop
    endif   
    !
    ! vs,vp,vd are in Hartee
    !
    if(iloc.ne.0.and.ic.eq.0) then
      do i=1,nrr
        read(it,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i) 
      enddo
    endif

    if(iloc.eq.0.and.ic.eq.0) then
      do i=1,nrr
       read(it,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),vloc(i)
      enddo
    endif

    if(iloc.ne.0.and.ic.eq.1) then
      do i=1,nrr
       read(it,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),rhoc(i)
      enddo
    endif

    if(iloc.eq.0.and.ic.eq.1) then
      do i=1,nrr
       read(it,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),vloc(i),rhoc(i)
      enddo
    endif
  !
  !  reads in Psi. It is then converted to Psi/r for FHI format
  !
  else if ( pw_params%pp_format .eq. 3) then   

    write(9,*) ' USING an FHI pseudopotential'
    write(9,*)   

    read(it,*) zatom,nproj
    norb=nproj+0.1

    iloc= pw_params%pp_iloc(nt)
    occ(:)=pw_params%pp_pwfocc(:,nt)
    write(9,10) occ(1:norb)
10  format(3F6.2, '  atomic occupations')

    read(it,*) dum(1),dum(2),dum(3),dum(4)
    do i=1,9
      read(it,*) dum(1),dum(2),dum(3)
    end do
   
    read(it,*) nrr,dum(1)

    nrr=nrr+1
    do i=2,nrr
      read(it,*) idum,r(i),ws(i),vs(i)
      ws(i)= ws(i)/r(i)
    enddo   
    ws(1)=ws(2);   vs(1)=vs(2)

    if(nproj>1) then
       if(occ(2)==-1d0) then
          write(9,*) 'ERROR: Found projector for p but no occupation was given'
          call mystop
       end if
       read(it,*) nrr,dum(1)
       nrr=nrr+1
       do i=2,nrr
          read(it,*) idum,r(i),wp(i),vp(i)
          wp(i)= wp(i)/r(i)
       enddo
       wp(1)=dzero; vp(1)=vp(2)
    else
       occ(2)=0d0
    endif

    if(nproj>2) then
       if(occ(2)==-1d0) then
          write(9,*) 'ERROR: Found projector for d but no occupation was given'
          call mystop
       end if
       read(it,*) nrr,dum(1)
       nrr=nrr+1
       do i=2,nrr
          read(it,*) idum,r(i),wd(i),vd(i)
          wd(i)= wd(i)/r(i)
       enddo
       wd(1)=dzero; vd(1)=vd(2)
    else
       occ(3)=0d0
    endif

    r(1)=dzero

    do i=1,nrr
      rhoc(i)=dzero
    end do

    do i=2,nrr
      read(it,*,END=100)  r(i),rhoc(i),rhoc1d(i),rhoc2d(i)
       rhoc(i)=rhoc(i)/4.0/pi
    enddo
    r(1)=dzero;
    rhoc(1)=rhoc(2); rhoc1d(1)=rhoc1d(2);  rhoc2d(1)=rhoc2d(2); 

100 continue

  else
    call mystop( ' iformat not 2 or 3 in w_linerq' )
  end if
    
  close(it)
  !
  !  END READING IN DATA
  !
  ! --------------------------
  !
  ! take the local potential
  !
  do i=1,nrr
      if(iloc.eq.0) vloc(i)=vloc(i)
      if(iloc.eq.1) vloc(i)=vs(i)
      if(iloc.eq.2) vloc(i)=vp(i)
      if(iloc.eq.3) vloc(i)=vd(i)
      if(iloc.eq.12) vloc(i)=(vs(i)+vp(i))/2
      if(iloc.eq.13) vloc(i)=(vs(i)+vd(i))/2
      if(iloc.eq.23) vloc(i)=(vp(i)+vd(i))/2
  enddo
  !
  ! calculate vql0
  !
  s=dzero
  ch=zatom
  do i=2,nrr-1
    if(r(i).lt.15.d0) then
      s=s+(ch*r(i)+vloc(i)*r(i)**2)*(r(i+1)-r(i-1))/2
    endif
  enddo
  vql0=s*4*pi
  vql0=2*vql0

  write(9,*) "vql0=",vql0
  call myflush(9)
  !
  ! Calculate vlocT : is for the use of Thomas procedure
  !
  do i=1,nrr
    if(r(i).lt.4.d0) then

      a = vs(i)*occ(1)*ws(i)**2 + vp(i)*occ(2)*wp(i)**2 + vd(i)*occ(3)*wd(i)**2

      b = occ(1)*ws(i)**2 + occ(2)*wp(i)**2 + occ(3)*wd(i)**2

      if (a .eq. dzero .and. b .eq. dzero) then
       vlocT(i)=vp(i) 
      else
        vlocT(i)=a/b
      end if
    else
      call myflush(9)
      vlocT(i)=a/b
      vlocT(i)=vp(i)
    endif
  enddo
  !
  !  Do 1d fourier transform of 
  !  local vloc,vlocT,valence charge, and core charge
  !
  z=zatom
  rst=10.d0
  
  do iq=2,pspot%mxdlqp   ! loop starts at 2 since G=0 is shifted in read_pseudo

    g1=(iq-2)*delql 
    if(g1.lt.1.D-3) g1=1.D-3

    s=0.d0
    s1=0.d0
    s2=0.d0
    s3=dzero 

    do i=1,nrr-1

      if(r(i).gt.rst) then
        rst1=r(i)
        goto 97
      endif

      r1=r(i)
      r2=r(i+1)
      x1=g1*r1
      x2=g1*r2
      b1=vloc(i)*r1
      b2=(vloc(i+1)*r2-b1)/(x2-x1)

      c1=(cos(x1)-cos(x2))*b1+(x1-x2)*b2*cos(x2)+ &
          b2*(sin(x2)-sin(x1))
      s=s+c1

      rhoi=(occ(1)*ws(i)**2+occ(2)*wp(i)**2+ &
       occ(3)*wd(i)**2)/(4*pi)
      rhoi1=(occ(1)*ws(i+1)**2+occ(2)*wp(i+1)**2+ &
      occ(3)*wd(i+1)**2)/(4*pi)

      b11=rhoi*r1
      b22=(rhoi1*r2-b11)/(x2-x1)
      c11=(dcos(x1)-dcos(x2))*b11+(x1-x2)*b22*dcos(x2)+ &
        b22*(dsin(x2)-dsin(x1))
      s1=s1+c11

      b1=vlocT(i)*r1
      b2=(vlocT(i+1)*r2-b1)/(x2-x1)
      c1=(cos(x1)-cos(x2))*b1+(x1-x2)*b2*cos(x2)+ &
       b2*(sin(x2)-sin(x1))
      s2=s2+c1

      b1=rhoc(i)*r1 
      b2=(rhoc(i+1)*r2-b1)/( (x2-x1) )
      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+ &
       b2*(dsin(x2)-dsin(x1))
      s3=s3+c1

    enddo

97  continue

    s=s/g1**2
    s=s*4*pi
    s=s-z*4*pi*dcos(g1*rst1)/g1**2

    s1=s1/g1**2
    s1=s1*4*pi

    s2=s2/g1**2
    s2=s2*4*pi
    s2=s2-z*4*pi*dcos(g1*rst1)/g1**2

    s3=s3/g1**2
    s3=s3*4*pi

    vda(iq)=s
    rhoq(iq)=s1
    vnl_atom(iq)=s2-s
    denc_atom(iq)=s3

  end do
  !
  ! calc. the nonlocal potential Kleiman-Bylander ref. dv*psi
  !
  ! isp=1, s state; 2, p state; 3, d state.
  !
  if (pw_params%NL_rspace(1)) then
    rst=rcut/1.4    
  else
    rst=4.d0
  end if

  do isp=1,norb

    s=0.d0
    do i=1,nrr

      if (pw_params%NL_rspace(1)) then
        if(r(i).gt.rst) then
          amrI=0.d0
        else
          r2=r(i)/rcut
          ir=1+r2*200.d0
          f1=(pspot%ri(ir+1)-r2)/(pspot%ri(ir+1)-pspot%ri(ir))
          f2=(r2-pspot%ri(ir))/(pspot%ri(ir+1)-pspot%ri(ir))
          amrI=pspot%amr(ir)*f1+pspot%amr(ir+1)*f2
          amrI=done/amrI
        endif
      end if
  
      if(isp.eq.1) then
        vw(i)=(vs(i)-vloc(i))*ws(i)
        if(r(i).lt.rst.and.i.gt.1) then
          s=s+ws(i)**2*(vs(i)-vloc(i))*(r(i+1)-r(i-1)) /2*r(i)**2 
        endif
      endif
  
      if(isp.eq.2) then
        vw(i)=(vp(i)-vloc(i))*wp(i)
        if(r(i).lt.rst.and.i.gt.1) then
          s=s+wp(i)**2*(vp(i)-vloc(i))*(r(i+1)-r(i-1))  /2*r(i)**2 
        endif

      endif
  
      if(isp.eq.3) then
        vw(i)=(vd(i)-vloc(i))*wd(i)
        if(r(i).lt.rst.and.i.gt.1) then
          s=s+wd(i)**2*(vd(i)-vloc(i))*(r(i+1)-r(i-1)) /2*r(i)**2 
        endif
      endif
  
      if (pw_params%NL_rspace(1)) then
        vw2(i)=vw(i)*amrI
      end if 
 
    enddo

    s=4*pi*s

    if(isp.eq.iloc) then
      scale=0.d0
    else
      scale=1/dsqrt(dabs(s))
    endif

    if(s.ge.0.d0) then
      pspot%nkb(isp,nt)=1
    else
      pspot%nkb(isp,nt)=-1
    endif
    if(dabs(s).lt.1.D-10) pspot%nkb(isp,nt)=0

    if(isp.eq.1) then

      do iq=2,pspot%mxdlqp ! loop starts at 2 as G=0 is shifted in read_pseudo

        g1=(iq-2)*delql

        if(iq.eq.2) g1=1.D-6

        s=dzero
        s2=dzero
        do i=2,nrr-1
          if(r(i).gt.rst) goto 96
          x=r(i)*g1
          if(x.gt.0.01) then
            aj0=dsin(x)/x
          else
            aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
          endif
          s=s+aj0*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
          s2=s2+aj0*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2
        enddo
96      continue
        s=s*4*pi
        s=s*scale

        s2=s2*4*pi
        s2=s2*scale
  
        pspot%qi(iq,nt)=g1
        pspot%vkb(iq,1,nt)=s
        pspot%vkb_mask(iq,1,nt)=s2
      enddo
    endif

    if(isp.eq.2) then

      do iq=2,pspot%mxdlqp ! loop starts at 2 as G=0 is shifted in read_pseudo 
        g1=(iq-2)*delql
       
       if(g1.lt.1.D-3) g1=1.D-3

        s=dzero
        s2=dzero
        do i=2,nrr-1
          if(r(i).gt.rst) goto 99
          x=r(i)*g1

          if(x.gt.0.1) then
            aj1=dsin(x)/x**2-dcos(x)/x
          else
            aj1=x/3.d0-x**3/30.d0+x**5/840.d0-x**7/45360.d0
          endif
          s=s+aj1*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
          s2=s2+aj1*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2
        enddo
99      continue
        s=s*4*pi
        s=s*scale

        s2=s2*4*pi
        s2=s2*scale

        pspot%qi(iq,nt)=g1
        pspot%vkb(iq,2,nt)=s
        pspot%vkb_mask(iq,2,nt)=s2
      enddo
    endif

    if(isp.eq.3) then

      do iq=2,pspot%mxdlqp  ! loop starts at 2 as G=0 is shifted in read_pseudo

        g1=(iq-2)*delql
        if(iq.eq.2) g1=1.D-6

        s=dzero
        s2=dzero
        do i=2,nrr-1
          if(r(i).gt.rst) goto 98
          x=r(i)*g1

          if(x.gt.0.2) then
            aj2=(3/x**3-1/x)*dsin(x)-3*dcos(x)/x**2
          else
            aj2=x**2/15.d0-x**4/210.d0+x**6/7560.d0-x**8/498960.d0 
          endif
          s=s+aj2*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
          s2=s2+aj2*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2        
        enddo
98      continue
        s=s*4*pi
        s=s*scale

        s2=s2*4*pi
        s2=s2*scale

        pspot%qi(iq,nt)=g1
        pspot%vkb(iq,3,nt)=s
        pspot%vkb_mask(iq,3,nt)=s2
     enddo
    endif

  end do
  !
  !  END OF CALC. KB FORM
  !
  !  -------------------------------------------------
  !
  !  put in proper factors for Ylm and convert to Ryd.
  !
  fac_s=dsqrt(1/2.d0)
  fac_p=dsqrt((2*1+1)/2.d0)       ! 2*l+1 special for paratec conven.
  fac_d=dsqrt((2*2+1)/2.d0)

  do i=2,pspot%mxdlqp
     vda(i)=vda(i)*2      ! Ryd
     vnl_atom(i)=vnl_atom(i)*2
     pspot%vkb(i,1,nt)=pspot%vkb(i,1,nt)*fac_s*2
     pspot%vkb(i,2,nt)=pspot%vkb(i,2,nt)*fac_p*2 
     pspot%vkb(i,3,nt)=pspot%vkb(i,3,nt)*fac_d*2

     pspot%vkb_mask(i,1,nt)=pspot%vkb_mask(i,1,nt)*fac_s*2
     pspot%vkb_mask(i,2,nt)=pspot%vkb_mask(i,2,nt)*fac_p*2 
     pspot%vkb_mask(i,3,nt)=pspot%vkb_mask(i,3,nt)*fac_d*2 
  enddo
  !
  ! SET some temp cariables to their storage in pspot
  !
  pspot%nrr(nt)=nrr
  pspot%r(1:nrr,nt)=r(1:nrr)

  pspot%wr_pp(1:nrr,1,nt)=ws(1:nrr)
  pspot%wr_pp(1:nrr,2,nt)=wp(1:nrr)
  pspot%wr_pp(1:nrr,3,nt)=wd(1:nrr)

  pspot%vr_pp(1:nrr,1,nt)=vs(1:nrr)
  pspot%vr_pp(1:nrr,2,nt)=vp(1:nrr)
  pspot%vr_pp(1:nrr,3,nt)=vd(1:nrr)

  return
  end

