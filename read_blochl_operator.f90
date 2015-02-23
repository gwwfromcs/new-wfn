!-*-F90-*- 
!
! Subroutine to set up the non-local Blochl operator required
! in the calculation of all electron matrix elements from the
! psdeudo-wavefunction representation
!    
! (1999) Chris Pickard
!        Paris/Kiel

Subroutine read_blochl_operator(gs,blop,pspot,crys)

  Use all_to_all_module
  Use radin_mod
  Use bessfn_mod
  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  !
  ! INPUT:
  ! -----
  !   
  Type (parallel_gspace) :: gs
  Type (pseudo_potential):: pspot
  Type (crystal)         :: crys 
  !
  ! OUTPUT:
  ! ------
  !
  Type(blochl_operator)  :: blop
  !
  ! ==================================================
  !
  ! LOCAL:
  ! -----
  !
  Integer :: ntp,np,max

  Integer,  Dimension(:),     Allocatable :: nrpts
  Real(dp), Dimension(:,:,:), Allocatable :: aepswfn

  Logical :: exists

  ! ==================================================
  !

  Write (9,'(/" Initialising Blochl Operators:",/&
       &" ------------------------------")') 
  Call myflush(9)

  ! ** Initialise what we can now

  Call blop_init

  ! ** If BLOP.DAT exists then we don't need to calculate the
  !    Blochl operators from scratch

  Inquire(file='BLOP.DAT',exist=exists)

  If(exists) Then  
     
     ! ** Read in the Blochl operator 

     Call blop_read

  Else
     
     ! ** Loop over species and calculate the Blochl operator

     Do ntp=1,blop%ntype

        ! * Read in the wavefunctions

        Call read_atomic_wavefunctions

        ! * Generate the projectors

        Call construct_projectors

        ! * Calculate the current matrix elements
        Call recon_matrix_elements

        ! * Count the number of projectors

        Do np=1,blop%nlprjmx
           If(blop%nkb(np,ntp).Ne.0) blop%nanl = blop%nanl + &
                (2*blop%lo(np,ntp)+1)*crys%natom(ntp)
        End Do

     End Do
     Call blop_write
     Deallocate(aepswfn)

  End If

  Write(9,'(/" Number of projectors: ",i5)') blop%nanl
  Call myflush(9)

  Deallocate(nrpts)

  Return


Contains

  ! Subroutine to initialise the Blochl operators

  Subroutine blop_init

    ! ** Allocate what we can of the Blochl-op structure now

    blop%nlprjmx = 6            ! Set to the max number of projectors
    blop%ntype   = pspot%ntype  ! Number of atomic species
    blop%mxdlqp  = pspot%mxdlqp ! The q-space grid size

    ! Allocate the BLOP

    Allocate(blop%nqnl(blop%ntype))
    Allocate(blop%delqnl(blop%ntype))
    Allocate(blop%nkb(blop%nlprjmx,blop%ntype))
    Allocate(blop%lo(blop%nlprjmx,blop%ntype))
    Allocate(blop%r_rec(blop%nlprjmx,blop%ntype))
    Allocate(blop%numl(0:3,blop%ntype))
    Allocate(blop%ln2np(0:3,2,blop%ntype))
    Allocate(blop%vkb(blop%mxdlqp,blop%nlprjmx,blop%ntype))
    Allocate(blop%delta_dia(0:3,2,2,blop%ntype))
    Allocate(blop%delta_para(0:3,2,2,blop%ntype))

    ! Local array

    Allocate(nrpts(blop%nlprjmx))

    ! Set the radial reciprocal space grid to that of the pseudopotential

    blop%nqnl   = pspot%nqnl
    blop%delqnl = pspot%delqnl

    ! Initialise

    blop%nkb      = 0
    blop%nanl     = 0
    blop%numl     = 0
    blop%ln2np    = 0
    blop%lo       = 0
    blop%r_rec    = 0.d0
    blop%delta_dia = 0.d0
    blop%delta_para   = 0.d0

    ! Allocate d2vkbdq2

    max=Maxval(blop%nqnl(:))

    Allocate(blop%d2vkbdq2(max-2,blop%nlprjmx,blop%ntype))   

    blop%d2vkbdq2=0.d0

  End Subroutine blop_init


  ! Subroutine to read in the atomic all electron and pseudo
  ! wavefunctions

  Subroutine read_atomic_wavefunctions

    Integer, Parameter :: unit_atom=11
    Integer            :: nrmax,np,npj,nr
    Character(len=80)  :: cjunk,ctemp
    Logical            :: exists2

    ! * Open the file

     
    !Does the AEPS file exist?

    Inquire(file=Trim(Adjustl(crys%nameat(ntp)))//'_AEPS.DAT',exist=exists2)


    ! If it does exist then go ahead and read it

    If(exists2) Then  

    If(gs%myproc.Eq.0) Then
       Open(unit=unit_atom,file=Trim(Adjustl(crys%nameat(ntp)))//'_AEPS.DAT',&
            form='formatted',status='old',position='rewind',err=100)
    End If

    Write (9,'(/" Reading from file: ",a)') &
         Trim(Adjustl(crys%nameat(ntp)))//'_AEPS.DAT'
    Call myflush(9)

    ! * Start reading

    nrpts = 0  
    npj   = 0

    ! Read the header

    If(gs%myproc.Eq.0) Read(unit_atom,'(a)',End=101) cjunk

    ! First count the number of projectors (l count only - not m)
    !  ----  and the max number of realspace grid points 

    If(gs%myproc.Eq.0) Then

       Do
          Read(unit_atom,'(a)',End=101) cjunk
          ctemp=Trim(Adjustl(cjunk))

          ! Count projectors and gridpoints

          If(ctemp(1:3).Eq.'#l=') Then

             npj=npj+1

             ! Trap the the case of too many projectors. 
             ! Just increase blop%nlprjmx

             If(npj.Gt.blop%nlprjmx) Then
                Write (9,'(a,i3,a)') 'npj.gt.blop%nlprjmx=',blop%nlprjmx,&
                     ' STOPPING'
             Endif

          Else
             If(ctemp(1:1).Ne.'&') nrpts(npj)=nrpts(npj)+1
          End If

       End Do

    End If

101 Call my_broadcast(npj,0)
    Call my_broadcast(nrpts(1),blop%nlprjmx,0)



    nrmax = Maxval(nrpts(1:blop%nlprjmx))
    If(Allocated(aepswfn)) Deallocate(aepswfn)

    Allocate(aepswfn(nrmax,npj,3))

    ! * Rewind, and now actually read in the data

    If(gs%myproc.Eq.0) Rewind(unit_atom)

    ! Read header

    If((gs%myproc.Eq.0).And.(npj.Gt.0)) Read(unit_atom,'(a)') cjunk

    ! Loop over l projectors
    ! change to MPI stuff in following loop :jry

    Do np=1,npj
       If(gs%myproc.Eq.0) Read(unit_atom,'(a)') cjunk ! #l=?
       Call my_broadcast(cjunk,80,0)
       Read(cjunk(4:5),*) blop%lo(np,ntp)
       Read(cjunk(12:),*) blop%r_rec(np,ntp)
       Do nr=1,nrpts(np)
          If(gs%myproc.Eq.0) Read(unit_atom,*) aepswfn(nr,np,1:3)
       End Do
          Call my_broadcast(aepswfn(:,np,1),nrpts(np),0)
          Call my_broadcast(aepswfn(:,np,2),nrpts(np),0)
          Call my_broadcast(aepswfn(:,np,3),nrpts(np),0)
       If(gs%myproc.Eq.0) Read(unit_atom,'(a)') cjunk 
    End Do


    ! * Close the file

    Close(unit_atom)

    Return

    ! If there is no wavefunction file for this atomic type

100 STOP 'something is wrong with the AEPS file'
    

 else   !the AEPS file doesn't exist, we carry on regardless 
            
    If(Allocated(aepswfn)) Deallocate(aepswfn)
    write(9,*) 'No AEPS file for ',crys%nameat(ntp)
    Allocate(aepswfn(0,0,3)) ; npj = 0 ; nrpts = 0
    aepswfn=0.d0
    return
 end if



  End Subroutine read_atomic_wavefunctions

  ! Subroutine to take pseudowavefunctions (partial waves) and
  ! generate a corresponding set of projectors:
  ! 
  ! |o_i> = SUM_j f |phi~_j> alpha_ji
  ! 
  ! <o_i|phi~_j> = SUM_k alpha_ik <phi~_k|f|phi~_j> = delta_ij
  !
  !
  ! At the moment we follow van der Walle and Blochl and make the 
  ! simplification: f is a third order polynomial, zero beyond rc
  !                 single projector per l-channel

  Subroutine construct_projectors

    Integer :: nrmax,npj,np,l,n,nr,nq,nrc,i,j,np1,np2

    Real(dp) :: q,rc,rs,fpow

    Real(dp) :: norm,s(2,2),sinv(2,2)
    Real(dp), Allocatable, Dimension(:)   :: work
    Real(dp), Allocatable, Dimension(:,:) :: projr

    Write (9,'(/" Constructing projectors")') 
    Call myflush(9)

    npj = Size(aepswfn,2)

    nrmax = Maxval(nrpts(1:blop%nlprjmx))
    Allocate(work(nrmax))
    Allocate(projr(nrmax,blop%nlprjmx))

    Do np =1,npj

       ! Number of grid points for this projector

       nr = nrpts(np)

       ! Find the angular momentum for this projector

       l = blop%lo(np,ntp)

       blop%numl(l,ntp)=blop%numl(l,ntp)+1

       ! Flag a projector

       blop%nkb(np,ntp)=1

       If(blop%numl(l,ntp).Gt.2) Stop '** too many projectors per channel'

       ! Fill in ln2np - indexing array

       blop%ln2np(l,blop%numl(l,ntp),ntp) = np

    End Do

    ! Set the inner core radius here for all channels

    rs  = 0.d0 ! Keep as zero

    ! Initialise the projectors

    projr = 0.d0

    ! 
    ! *** Construct the projectors in real space
    ! 

    fpow = 1.0d0 ! fpow>=1.d0 for smoothing, fpow=0.d0 for step fn.

    Do l=0,3

       ! write (9,*) 'l:',l,'numl:',blop%numl(l,ntp)

       If(blop%numl(l,ntp).Eq.1) Then

          np  = blop%ln2np(l,1,ntp)
          nr  = nrpts(np)
          nrc = Count(aepswfn(1:nr,np,1).Le.blop%r_rec(np,ntp))
          rc  = aepswfn(nrc,np,1)

          work(1:nr) = aepswfn(1:nr,np,3)*aepswfn(1:nr,np,3)            
          work(1:nr) = apply_f(work(1:nr), aepswfn(1:nr,np,1),rs,rc,fpow)

          norm = radin(work(1:nr),aepswfn(1:nr,np,1),nrc)

          projr(1:nr,np) = apply_f(aepswfn(1:nr,np,3)/norm,&
               aepswfn(1:nr,np,1),rs,rc,fpow)

       Else If(blop%numl(l,ntp).Eq.2) Then

          ! ** Rescale the wavefunctions so that
          !         int_0^rc f|psi_ps|^2 = 1

          Do i=1,blop%numl(l,ntp)

             np  = blop%ln2np(l,i,ntp)
             nr  = nrpts(np)
             nrc = Count(aepswfn(1:nr,np,1).Le.blop%r_rec(np,ntp))
             rc  = aepswfn(nrc,np,1)

             work(1:nr) = apply_f(aepswfn(1:nr,np,3)**2,&
                  aepswfn(1:nr,np,1),rs,rc,fpow)

             norm = radin(work(1:nr),aepswfn(1:nr,np,1),nrc)

             aepswfn(1:nr,np,2:3) = aepswfn(1:nr,np,2:3)/Sqrt(norm)

          End Do

          ! ** Construct the overlap matrix

          ! Loop over 2x2 block

          Do i=1,2
             np1 = blop%ln2np(l,i,ntp)
             nr  = nrpts(np1) 
             nrc = Count(aepswfn(1:nr,np1,1).Le.blop%r_rec(np1,ntp))
             rc  = aepswfn(nrc,np1,1)
             Do j=1,2

                np2=blop%ln2np(l,j,ntp)

                If(nrpts(np1).Ne.nrpts(np2)) Stop '*** incompatible grids'

                If(aepswfn(nr,np1,1).Ne.aepswfn(nr,np2,1)) &
                     Stop '*** incompatible grids'

                work(1:nr) = aepswfn(1:nr,np1,3)*aepswfn(1:nr,np2,3)

                work(1:nr)=apply_f(work(1:nr), aepswfn(1:nr,np1,1),&
                     rs,rc,fpow)

                s(i,j) = radin(work(1:nr),aepswfn(1:nr,np1,1),nrc)

             End Do
          End Do

          ! ** Invert the overlap matrix

          Call invert_2x2(s,sinv)

          ! ** Construct the projectors in realspace

          !       Loop over 2x2 block

          Do i=1,2

             np1 = blop%ln2np(l,i,ntp)
             nr  = nrpts(np1)   
             nrc = Count(aepswfn(1:nr,np1,1).Le.blop%r_rec(np1,ntp))
             rc  = aepswfn(nrc,np1,1)

             projr(:,np1) = 0.d0

             Do j=1,2
                np2=blop%ln2np(l,j,ntp)

                projr(1:nr,np1) = projr(1:nr,np1)+sinv(i,j)*aepswfn(1:nr,np2,3)

             End Do

             projr(1:nr,np1) =  apply_f(projr(1:nr,np1),aepswfn(1:nr,np1,1),&
                  rs,rc,fpow)

          End Do


!!$          np1 = blop%ln2np(l,1,ntp)
!!$          np2 = blop%ln2np(l,2,ntp)
!!$
!!$          nr  = nrpts(np1)
!!$
!!$          write (9,*)
!!$          write (9,*) 'l=',l
!!$          write (9,*)
!!$
!!$          write (9,'(a,f15.10)') ' 11:',radin(projr(1:nr,np1)*&
!!$               aepswfn(1:nr,np1,3),aepswfn(1:nr,np1,1),nr)
!!$          write (9,'(a,f15.10)') ' 12:',radin(projr(1:nr,np1)*&
!!$               aepswfn(1:nr,np2,3),aepswfn(1:nr,np1,1),nr)
!!$          write (9,'(a,f15.10)') ' 22:',radin(projr(1:nr,np2)*&
!!$               aepswfn(1:nr,np2,3),aepswfn(1:nr,np1,1),nr)
!!$          do i=1,maxval(nrpts)
!!$             if(aepswfn(i,np1,1).le.rc) write (21,'(3f20.10)') &
!!$                  aepswfn(i,np1,1),projr(i,np1),projr(i,np2)
!!$          end do
!!$          write (21,*) '&'

       Else

          If(blop%numl(l,ntp).Ne.0) Stop ' *** too many projectors'

       End If

    End Do

    !
    ! *** Now transform into reciprocal space
    !

    Do np=1,npj

       nr  = nrpts(np)
       nrc = Count(aepswfn(1:nr,np,1).Le.blop%r_rec(np,ntp))
       rc  = aepswfn(nrc,np,1)

       Do nq = 2,blop%nqnl(ntp)

          q = Dble(nq-2)*blop%delqnl(ntp)

          Do n = 1,nrc
             work(n)=bessfn(aepswfn(n,np,1)*q,&
                  blop%lo(np,ntp))*aepswfn(n,np,1)*projr(n,np)
          End Do

          blop%vkb(nq,np,ntp)=radin(work(1:nr),aepswfn(1:nr,np,1),nrc)

       End Do

       !
       !   "Do some strange extrapolation to 0 i.e. extend vkb to small
       !   negative ql according to parity (-1)**l so that interpolation
       !   works around 0."
       !
       blop%vkb(1,np,ntp)=blop%vkb(3,np,ntp)*(-1)**blop%lo(np,ntp)

       ! Now find the derivatives

       Call reg_grid_spline(blop%vkb(3,np,ntp),blop%nqnl(ntp)-2,&
            blop%d2vkbdq2(1,np,ntp)) 

    End Do

    Deallocate(work,projr)

  End Subroutine construct_projectors

  ! Subroutine to generate the atomic matrix elements required 
  ! for the Blochl reconstruction 
  ! 

  Subroutine recon_matrix_elements

    Integer :: nrmax,nr,np1,np2,l,i,j,nrc

    Real(dp), Allocatable :: work(:)
    Real(dp) :: atomic

    Write (9,'(/" Calculating reconstruction matrix elements")') 
    Call myflush(9)

    nrmax = blop%nlprjmx
    Allocate(work(Maxval(nrpts(1:nrmax))))

    ! Loop over channels

    Do l=0,3

       Do i=1,blop%numl(l,ntp)
          np1 = blop%ln2np(l,i,ntp)
          nr  = nrpts(np1)   
          nrc = Count(aepswfn(1:nr,np1,1).Le.blop%r_rec(np1,ntp))

          Do j=1,blop%numl(l,ntp)
             np2=blop%ln2np(l,j,ntp)

             ! *** The charge augmentation term

             work(1)    = 0.d0
             work(2:nr) = (aepswfn(2:nr,np1,2)*aepswfn(2:nr,np2,2)-&
                  aepswfn(2:nr,np1,3)*aepswfn(2:nr,np2,3))/aepswfn(2:nr,np1,1)

             ! Factor of 1/c^2 *2 for spin *1/3 for Trace/3

             blop%delta_dia(l,i,j,ntp) = radin(work(1:nr),&
                  aepswfn(1:nr,np1,1),nrc)/(137.036d0**2.d0)*1.d6*2.d0/3.d0

             ! *** The current augmentation term

             work(1)    = 0.d0
             work(2:nr) = (aepswfn(2:nr,np1,2)*aepswfn(2:nr,np2,2)-&
                  aepswfn(2:nr,np1,3)*aepswfn(2:nr,np2,3))/&
                  aepswfn(2:nr,np1,1)**3

             ! Factor of 1/c^2 *2 for spin *2 for Ryd->Har  

             blop%delta_para(l,i,j,ntp) =  radin(work(1:nr),&
                  aepswfn(1:nr,np1,1),nrc)/(137.036d0**2.d0)*1.d6*4.d0

          End Do

       End Do

    End Do


    Deallocate(work)

  End Subroutine recon_matrix_elements

  ! Read in the Blochl operator for this species

  Subroutine blop_read

    Integer, Parameter :: unit_blop=12
    Integer  :: itemp

    If(gs%myproc.Eq.0) Then
       Open(unit=unit_blop,file='BLOP.DAT',&
            form='unformatted',status='old',position='rewind',err=111)  
       Read(unit_blop) itemp
       If(itemp.Ne.blop%nlprjmx) Goto 111
       Read(unit_blop) itemp
       If(itemp.Ne.blop%ntype) Goto 111
       Read(unit_blop) itemp 
       If(itemp.Ne.blop%mxdlqp) Goto 111
       Read(unit_blop) blop%nqnl
       Read(unit_blop) blop%delqnl
       Read(unit_blop) blop%nkb
       Read(unit_blop) blop%lo
       Read(unit_blop) blop%r_rec 
       Read(unit_blop) blop%numl
       Read(unit_blop) blop%ln2np
       Read(unit_blop) blop%vkb
       Read(unit_blop) blop%delta_dia
       Read(unit_blop) blop%delta_para
       Read(unit_blop) blop%d2vkbdq2
       Read(unit_blop) blop%nanl
       Close(unit_blop)
    Endif

    ! Now broadcast to all nodes

    Call my_broadcast(blop%nanl,0)
    Call my_broadcast(blop%nqnl(1),blop%ntype,0)
    Call my_broadcast(blop%delqnl(1),blop%ntype,0)
    Call my_broadcast(blop%nkb(1,1),blop%nlprjmx*blop%ntype,0)
    Call my_broadcast(blop%lo(1,1),blop%nlprjmx*blop%ntype,0)
    Call my_broadcast(blop%r_rec(1,1),blop%nlprjmx*blop%ntype,0)
    Call my_broadcast(blop%numl(0,1),4*blop%ntype,0)
    Call my_broadcast(blop%ln2np(0,1,1),8*blop%ntype,0)
    Call my_broadcast(blop%vkb(1,1,1),blop%mxdlqp*blop%nlprjmx*blop%ntype,0)
    Call my_broadcast(blop%delta_dia(0,1,1,1),16*blop%ntype,0)
    Call my_broadcast(blop%delta_para(0,1,1,1),16*blop%ntype,0)
    Call my_broadcast(blop%d2vkbdq2(1,1,1),&
         (Maxval(blop%nqnl(:))-2)*blop%nlprjmx*blop%ntype,0)

    Return

111 Stop 'Problems reading the Blochl operator'

  End Subroutine blop_read

  ! Write out the Blochl operator for this species  

  Subroutine blop_write

    Integer, Parameter :: unit_blop=12

    If(gs%myproc.Eq.0) Then
       Open(unit=unit_blop,file='BLOP.DAT',&
            form='unformatted',status='replace',position='rewind',err=222)  
       Write(unit_blop) blop%nlprjmx
       Write(unit_blop) blop%ntype
       Write(unit_blop) blop%mxdlqp
       Write(unit_blop) blop%nqnl
       Write(unit_blop) blop%delqnl
       Write(unit_blop) blop%nkb
       Write(unit_blop) blop%lo
       Write(unit_blop) blop%r_rec 
       Write(unit_blop) blop%numl
       Write(unit_blop) blop%ln2np
       Write(unit_blop) blop%vkb
       Write(unit_blop) blop%delta_dia
       Write(unit_blop) blop%delta_para
       Write(unit_blop) blop%d2vkbdq2
       Write(unit_blop) blop%nanl    
       Close(unit_blop)
    Endif

    Return

222 Stop 'Problems writing Blochl operator'

  End Subroutine blop_write

  ! Multiply by f=1-3(r/rc)^2+2(r/rc)^3 for r<rc

  Function apply_f(func,r,rs,rc,pow)

    Real(dp), Intent(in) :: func(:),r(:)
    Real(dp), Intent(in) :: rs,rc,pow

    Real(dp) :: apply_f(Size(r))

    ! Local

    Integer :: n,i,nrc,nrs
    Real(dp) :: rcp,rsp

    If(Size(r).Ne.Size(func)) Stop 'apply_f incorrectly called'

    n=Size(r)

    ! Put rc,rs to grid point

    nrc = Count(r(:n).Le.rc)
    nrs = Count(r(:n).Le.rs)
    rcp = r(nrc)
    rsp = r(nrs)

    Do i=1,n
       If(r(i).Le.rsp) Then
          apply_f(i) = func(i)
       Else
          If(r(i).Le.rcp) Then
             apply_f(i)=func(i)* (1.d0-3.d0*((r(i)-rsp)/(rcp-rsp))**2+ &
                  2.d0*((r(i)-rsp)/(rcp-rsp))**3)**pow
          Else
             apply_f(i)=0.d0
          Endif
       Endif

    End Do

  End Function apply_f

  ! Subroutine to find the inverse of a 2x2 matrix

  Subroutine invert_2x2(a,ainv)

    Real(dp), Intent(in)  :: a(2,2)
    Real(dp), Intent(out) :: ainv(2,2)

    ! Local

    Real(dp) :: det

    det = a(1,1)*a(2,2)-a(1,2)*a(2,1)

    ainv(1,1) = a(2,2)
    ainv(2,2) = a(1,1)
    ainv(1,2) = -a(2,1)
    ainv(2,1) = -a(1,2)

    ainv = ainv/det

    !  write (9,*)
    !  write (9,'(2f10.5)') a
    !  write (9,*) 
    !  write (9,'(2f10.5)') ainv
    !  write (9,*)

  End Subroutine invert_2x2

End Subroutine read_blochl_operator






