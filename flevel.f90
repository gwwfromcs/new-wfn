!
!     @process extchk
!
subroutine flevel(imode, ipr, pw_params, energs, crys, bands, kpoints)

  include 'use.h'  
  implicit none                       ! implicit? No!
  include 'interface.h'  
  !
  !     subroutine finds the fermi level and computes
  !     nband, ifmax, and the (possibly fractional)
  !     level occupancy for a metallic solid.
  !
  !     1996    Bernd Pfrommer, while at UC Berkeley
  !     2001    Dave Raczkowski, major revision based on code from
  !                              Nicola Marzari
  !
  !     for the spin-polarized case, the eigenvalues must be shifted
  !     by vin(G=0), which is different for spin up/spin down. vin(G=0)
  !     is v_xc(G=0), because G=0 for the Hartree- and ionic potential is
  !     set to 0.
  !     later in the charge subroutine, the fermi energy must be shifted back
  !
  !   FD, COLD SMEARING, SPLINE OF GAUSSIANS, AND A CORRECTION TO METH-PAX
  !   INTRODUCED BY N. MARZARI by way of D. Raczkowski 
  !
  !   FOR A DETAILED EXPLANATION OF THE SMEARING APPROACH, LOOK AT
  !   N. MARZARI PhD THESIS, CHAPTER 4, AND REFERENCES THEREIN:
  !   http://alfaromeo.princeton.edu/~marzari/preprints/#phd
  !
  !   Cold smearing reference: N. Marzari, D. Vanderbilt, A. De Vita
  !   and M. C. Payne, ``Thermal contraction and disordering of the
  !   Al(110) surface'', Phys. Rev. Lett. 82, 3296 (1999).
  !
  !   GIVEN A SET OF WEIGHTS AND THE EIGENVALUES ASSOCIATED TO THEIR
  !   ASSOCIATED K-POINTS FOR BZ SAMPLING, THIS SUBROUTINE PERFORMS
  !   the tasks of calculating smearings and the entropy:
  !
  !   THE SMEARING SCHEMES 
  !
  !   (1) GAUSSIAN:
  !   SEE: C-L FU AND K-M HO, PHYS. REV. B 28, 5480 (1983).
  !   THEIR IMPLEMENTATION WAS VARIATIONAL BUT *NOT* CORRECTED FOR
  !   SECOND ORDER DEVIATION IN SIGMA AS ALSO WAS THE SIMILAR SCHEME
  !   (WITH OPPOSITE SIGN DEVIATION) IN: R.J.NEEDS, R.M.MARTIN AND O.H.
  !   NIELSEN, PHYS. REV. B 33 , 3778 (1986).
  !   USING THE CORRECTION CALCULATED HEREAFTER EVERYTHING SHOULD BE OK.
  !   THE SMEARING FUNCTION IS A GAUSSIAN NORMALISED TO 2.
  !   THE OCCUPATION FUNCTION IS THE ASSOCIATED COMPLEMENTARY
  !   ERROR FUNCTION.
  !
  !   (2) FERMI-DIRAC:
  !   SEE: M.J.GILLAN J. PHYS. CONDENS. MATTER 1, 689 (1989), FOLLOWING
  !   THE SCHEME OUTLINED IN J.CALLAWAY AND N.H.MARCH, SOLID STATE PHYS. 38,
  !   136 (1984), AFTER D.N.MERMIN, PHYS. REV 137, A1441 (1965).
  !   THE OCCUPATION FUNCTION IS TWICE THE SINGLE ELECTRON
  !   FERMI-DIRAC DISTRIBUTION.
  !
  !   (3) HERMITE-DELTA_EXPANSION:
  !   SEE: METHFESSEL AND PAXTON, PHYS. REV.B 40, 3616 (1989).
  !   THE SMEARING FUNCTION IS A TRUNCATED EXPANSION OF DIRAC'S DELTA
  !   IN HERMITE POLINOMIALS.
  !   FOR THE SMEARING FUNCTION IMPLEMENTED HERE THE TRUNCATION IS
  !   AT THE FIRST NON-TRIVIAL EXPANSION TERM D1(X).
  !   THE OCCUPATION FUNCTION IS THE ASSOCIATED PRIMITIVE.
  !   (NOTE: THE OCCUPATION FUNCTION IS NEITHER MONOTONIC NOR LIMITED
  !   BETWEEN 0. AND 2)
  !
  !   THE ENTROPY CORRECTION HOLDS UP TO THE THIRD ORDER IN DELTA AT LEAST,
  !   AND IS NOT NECESSARY (PUT = 0.) FOR THE HERMITE_DELTA EXPANSION,
  !   SINCE THE LINEAR ENTROPY TERM IN SIGMA IS ZERO BY CONSTRUCTION
  !   IN THAT CASE. (well, we still need the correct free energy. hence
  !   esmear is set to its true value, nmar)
  !
  !   (4) GAUSSIAN SPLINES:
  !   similar to a Gaussian smearing, but does not require the
  !   function inversion to calculate the gradients on the occupancies.
  !   It is thus to be preferred in a scheme in which the occ are
  !   independent variables  
  !
  !   (5) POSITIVE HERMITE, or COLD SMEARING I
  !   similar to Methfessel-Paxton (that is killing the linear order
  !   in the entropy), but with positive-definite occupations (still
  !   greater than 1) ! 
  !
  !   (6) POSITIVE HERMITE, or COLD SMEARING II: the one to use.
  !   (5) and (6) are practically identical; (6) is more elegant.
  !
  !     ARGUMENTS:
  !     ---------
  !
  integer, intent(in) :: imode, &                         ! mode of operation
       ipr                                                       ! print flag
  type(pw_parameter), intent(in) :: pw_params  
  !     get ioccup from pw_params
  !            ioccup = 0 standard fermi level smearing
  !            ioccup = 1 fill lowest ztot bands
  !            ioccup = 2 fill bands antiferromagnetically
  !            ioccup = 3 excited state - HOMO and LUMO half-filled - not here
  !     also get the fermi level smearing dsmear from here
  !
  type(energy), intent(inout) :: energs  
  !            input: energs%vxc0
  !            ouput: energs%efermi
  type(crystal), intent(in) :: crys             !            input: crys%ztot
  type(kpoint), intent(in) :: kpoints                   ! band structure info
  type(band), intent(inout) :: bands  
  !
  !          input:  bands%energy
  !          input:  bands%nband
  !          output: bands%occup ! this is really 2*occup*weight(kpoint)
  !          output: bands%ifmax
  !     ---------------------- local variables ---------------------------
  !
  real(dp) :: emin, emins(2), sumt, sumw, xdum, &
       sumwprime, deriv_contrib, derivative,ef_min,ef_max, &
       inv_smear, &
       dsmear, &                     ! gaussian smearing parameter in Rydberg
       desm, fraction, weight_diff,weight_diff_old, &
       dnel, &                             ! sum for total number of electrons
       emax, e0, z1, x, FMID, F, RTBIS, DX, XMID, FI, sqrpi, zeta, eesh, sq2i,piesqq,a, ee, sqr2
  real(dp) FERMID, DELTHM, SPLINE, POSHM, POSHM2
  real(dp) :: elowerbracket, eupperbracket, eguess ! added by David Roundy
  real(dp), external :: myderfc 
  real(dp) nmerfc 
  logical :: iafmag                      ! flag if antiferromagnetic/magnetic
  integer :: nrk, &                                      ! number of k-points
       ioccup, &                          ! occupy levels according to ioccup
       ispn, &                                      ! number of spins: 1 or 2
       i, idum, j, nel, nelp, itcount, is, tnel, n, ismear, nmax, ab(2),&
       icnt,k, gw_band_max 
  integer, parameter :: itmax = 10000,&   ! max number of iterations for efermi
           jmax =200              ! THE MAX NUMBER OF BISECTIONS TO GET EF

  real(dp), parameter :: range = dfour, &  ! range beyond which erfc is undef
       delta = 1.0d-8,&
           xacc=1.0D-10             ! THE DESIRED ACCURACY ON EF
  integer, allocatable :: ind(:), irk(:), gw_num_band(:,:),ind_tmp(:)  
  real(dp), allocatable :: el(:), &   ! the eigenvalues as a contiguous array
       ws(:)                        ! the weights of the k-points, for ispn=2
  real(dp) en
  !
  !     ------------------------------------------------------------------

  ! **************************************************************************
  ! **            Insulator only -- fixed spin polarisation                 **
  ! **************************************************************************

  if(pw_params%occupy == 4) then

     ! ** Set the eigenvalue shift
   
     energs%evalshift = energs%vxc0  

     ab(1) = crys%nalpha
     ab(2) = crys%nbeta

     ! ** Zero the occupations

     bands%occup = 0.d0

     ! ** Set the non-zero occupations

     do is = 1,bands%nspin
        do i = 1, bands%nrk   

           do n = 1, bands%nband(i, is)

              if(bands%nspin == 1) then
                 if(n<=(ab(1)+ab(2)+1)/2) bands%occup(n,i,is) = dtwo * kpoints%w(i)  
              else
                 if(n<=ab(is)) bands%occup(n,i,is) = dtwo * kpoints%w(i)  
              end if
              
           end do

        end do
     end do

     ! ** Set ifmax, same for all k-points

     do is = 1,bands%nspin
        if(bands%nspin == 1) then
           bands%ifmax(:,is) = (ab(1)+ab(2)+1)/2
        else
           bands%ifmax(:,is) = ab(is) 
        end if
     end do

     ! ** Check that we have the right number of electrons for non-spin polarised case

     if(bands%nspin==1) then
        dnel = sum(bands%occup(1:bands%ifmax(1,1),1:bands%nrk,1)) - real(sum(ab),dp)
        do i=1,bands%nrk
           bands%occup(bands%ifmax(i,1),i,1) =  bands%occup(bands%ifmax(i,1),i,1) - dnel * kpoints%w(i)  
        end do
     end if

     ! ** Find the Fermi level

     if (imode == 0) then  
        energs%efermi(1) = energs%inputfermi
        energs%efermi(2) = energs%inputfermi
     else
        do is = 1,bands%nspin
           energs%efermi(is) = -huge(1.d0)
           do i = 1,bands%nrk
              do j=1,bands%ifmax(i,is)
                 en = bands%energy(j, i, is) + energs%vxc0(is)
                 if(en > energs%efermi(is)) energs%efermi(is) = en
              end do
           end do
        end do
     end if

     en = maxval(energs%efermi(:))
     energs%efermi(:) = en

     ! ** Print out a summary

  if (ipr >= 1) then  
        do is = 1, bands%nspin  
           if (bands%nspin == 2 .and. is == 1) write(9, '(1x,a)') 'Spin up'  
           if (bands%nspin == 2 .and. is == 2) write(9, '(1x,a)') 'Spin down'  
           do i = 1, bands%nrk   
              write(9, 500) i  
              do n = 1, bands%nband(i, is), 7  
                 write(9, 510) (bands%energy(j, i, is) + energs%vxc0(is), &
                      j = n, min(bands%nband(i, is), n + 6))
                 write(9, 515) (bands%occup(j, i, is) *dhalf / kpoints%w(i), &
                      j = n, min(bands%nband(i, is), n + 6))
              end do
           end do
        end do
     end if

     return
  end if  !F.Mauri

  if(pw_params%occupy == 5) then 

  if (ipr >= 1) then
        do is = 1, bands%nspin
           if (bands%nspin == 2 .and. is == 1) write(9, '(1x,a)') 'Spin up'
           if (bands%nspin == 2 .and. is == 2) write(9, '(1x,a)') 'Spin down'
           do i = 1, bands%nrk
              write(9, 500) i
              do n = 1, bands%nband(i, is), 7
                 write(9, 510) (bands%energy(j, i, is) + energs%vxc0(is), &
                      j = n, min(bands%nband(i, is), n + 6))
                 write(9, 515) (bands%occup(j, i, is) *dhalf / kpoints%w(i), &
                      j = n, min(bands%nband(i, is), n + 6))
              end do
           end do
        end do
     end if
   return
   end if
 



  
 
  ! **************************************************************************


  !
  !     transfer data from the structures to local variables
  !
  sqrpi = sqrt(pi)  
  ee=exp(1.0)
  eesh=sqrt(ee)*0.5
  sq2i=sqrt(2.0)*0.5
  sqr2=sqrt(2.0)
  piesqq=sqrt(ee*pi)*0.25

  ismear=pw_params%smearing_method
  ioccup = pw_params%occupy  
  if (ioccup == 2) then  
     iafmag = .true.  
     ioccup = 0  
  else  
     iafmag = .false.  
  end if
  if (ioccup == 3) ioccup = 0
  nrk = bands%nrk  
  ispn = bands%nspin  
  dsmear = pw_params%smearing / ryd  
  !     ------------------------------------------------------------------
  !                     allocate dynamic arrays
  ! total number of electrons possible
  !
  nelp = nrk * bands%max * ispn  

! P. Zhang
  if(iafmag) then
   bands%nband1=MINVAL(bands%nband(:,:))
   nelp=nrk*bands%nband1*ispn
   end if

  allocate(ind(nelp))  
  allocate(irk(nelp))  
  allocate(el(nelp))
  allocate(ws(nrk * ispn))  
  !
  energs%evalshift = energs%vxc0  
  ! initialize some arrays....
  nel = 0  
! P. Zhang
  if (iafmag) then

  do is = 1, ispn
     do i = 1, nrk
        do j = 1, bands%nband1
           irk(nel + j) = i
           el(nel + j) = bands%energy(j, i, is) + energs%vxc0(is)
        end do
        nel = nel + bands%nband1
        ws(i + nrk * (is - 1)) = kpoints%w(i)
     end do
  end do
!-------------------------------------------

  else

  do is = 1, ispn
     do i = 1, nrk
        do j = 1, bands%nband(i, is)
           irk(nel + j) = i
           el(nel + j) = bands%energy(j, i, is) + energs%vxc0(is)
        end do
        nel = nel + bands%nband(i, is)
        ! copy weights over into temp array
        ws(i + nrk * (is - 1)) = kpoints%w(i)
     end do
  end do
  end if

  !
  !     sort the energy levels
  !
  ! sort spin up/down separately
  if (iafmag) then  
     call sort(nel / 2, el, ind)  
     call sort(nel / 2, el(1 + nel / 2), ind(1 + nel / 2))  
  else  
     call sort(nel, el, ind)  
  end if
  !
  !     Naive section, straight counting.
  !
  !      sum up the weight factors until we reach the maximum number
  !      allowed
  !
!-----------------------------------------------------------------
!Peihong Zhang, allow charged system
!  sumt = crys%ztot * dhalf * real(ispn, dp)

  sumt = (crys%ztot+crys%net_charge) * dhalf * real(ispn, dp)
!-----------------------------------------------------------------
  if(iafmag) then

  sumt=sumt/2.d0
  nel=nel/2

  do is=1,2
  sumw = dzero
  do i = 1, nel
     sumw = sumw + ws(irk(ind(i)))  
     if (sumw > sumt + delta) exit             ! leave the do loop
  end do

  !
  !        check if enough bands have been computed
  !
  if (sumw > sumt - delta) then  
     if (i > nel) then  
        i = i - 1  
     else  
        sumw = sumw - ws(irk(ind(i)))          ! subtract overcounting
     end if
  else  
     if (i > nel) i = i - 1  
     write(9, 600) dtwo * (sumt - sumw)  
     sumt = sumw  
  end if

  emin = el(ind(i+(is-1)*nel)) - dsmear   
  emax=el(ind(i+(is-1)*nel)) + 0.5*dsmear      
  e0= el(ind(i+(is-1)*nel))
  inv_smear = done / max (dsmear, 1.0d-6)  
  !  
  ! ensure emin and emax are proper bounds. ie. the solution lies between
  !
  nmax=1
101   continue
  !
  ! FMID = FUNC(emax) in Numerical Recipes.
  !
  Z1=dzero
  do i=1,nel 
    X = (emax - el(i+(is-1)*nel))/dsmear

    if(ismear.eq.1) then
      Z1 = Z1 + ws(irk(i))*myderfc(-x) 
    elseif(ismear.eq.2) then
      Z1 = Z1 + ws(irk(i))*FERMID(-X)
    elseif(ismear.eq.3) then
      Z1 = Z1 + ws(irk(i))*DELTHM(X)
    elseif(ismear.eq.4) then
      Z1 = Z1 + ws(irk(i))*SPLINE(-X)
    elseif(ismear.eq.5) then
      Z1 = Z1 + ws(irk(i))*POSHM(X)
    elseif(ismear.eq.6) then
      Z1 = Z1 + ws(irk(i))*POSHM2(X)
    endif
  end do
  Z1=dhalf*Z1
  
  FMID= Z1-sumt  
  !
  ! F = FUNC(emin)
  !
  Z1=dzero
  do i=1,nel 
    X = (emin - el(i+(is-1)*nel))/dsmear

    if(ismear.eq.1) then
      Z1 = Z1 + ws(irk(i))*myderfc(-x) 
    elseif(ismear.eq.2) then
      Z1 = Z1 + ws(irk(i))*FERMID(-X)
    elseif(ismear.eq.3) then
      Z1 = Z1 + ws(irk(i))*DELTHM(X)
    elseif(ismear.eq.4) then
      Z1 = Z1 + ws(irk(i))*SPLINE(-X)
    elseif(ismear.eq.5) then
      Z1 = Z1 + ws(irk(i))*POSHM(X)
    elseif(ismear.eq.6) then
      Z1 = Z1 + ws(irk(i))*POSHM2(X)
    endif
  end do
  Z1=dhalf*Z1
  
  F= Z1-sumt

  if(F*FMID .gt. dzero) then
    if (nmax.ge.10000) then
      write(9,*) 'ERROR: NO FERMI ENERGY FOUND ',emin,emax,F,FMID,nmax
      write(9,*) ' '
      write(9,*) 'IS THE ELECTRONIC TEMPERATURE TOO SMALL ? ',dsmear
      write(*,*) ' '
      call mystop
    else
      nmax=nmax+1
      emin=e0-real(nmax,dp)*dsmear
      emax=e0+(real(nmax,dp)-0.5)*dsmear
      goto 101
    end if
  end if

  if(F .le. dzero) then
    RTBIS = emin
    DX = emax - emin
  else
    RTBIS = emax
    DX = emin - emax
  end if
  !
  ! now that we have proper bounds, do bisection 
  !
  do J = 1, jmax
    DX = DX * 0.5D0
    XMID = RTBIS + DX
    !
    ! FMID=FUNC(XMID)
    !  
    Z1=dzero
    do i=1,nel 
      X = (xmid - el(i+(is-1)*nel))/dsmear

      if(ismear.eq.1) then
        Z1 = Z1 + ws(irk(i))*myderfc(-x) 
      elseif(ismear.eq.2) then
        Z1 = Z1 + ws(irk(i))*FERMID(-X)
      elseif(ismear.eq.3) then
        Z1 = Z1 + ws(irk(i))*DELTHM(X)
      elseif(ismear.eq.4) then
        Z1 = Z1 + ws(irk(i))*SPLINE(-X)
      elseif(ismear.eq.5) then
        Z1 = Z1 + ws(irk(i))*POSHM(X)
      elseif(ismear.eq.6) then
        Z1 = Z1 + ws(irk(i))*POSHM2(X)
      endif
    end do
    Z1=dhalf*Z1

    FMID= Z1-sumt

    if(FMID .le. 0.D0) RTBIS=XMID

    if(abs(DX) .lt. xacc .or. abs(FMID) .lt. 1d-12) then
       goto 11
    end if
  
  end do

  write(9,*) 'CANNOT BISECT FOREVER, CAN I ?',FMID,DX,xacc
  CALL EXIT

11  energs%efermi(is) = RTBIS
    write(9,*) "fermi level ",is,RTBIS
  end do

  energs%esmear=dzero

  else  ! not AFM occupation,i.e., both spine have the same fermi level

  sumw = dzero
  do i = 1, nel  
     sumw = sumw + ws(irk(ind(i)))  
     if (sumw > sumt + delta) exit             ! leave the do loop
  end do

  !
  !        check if enough bands have been computed
  !
  if (sumw > sumt - delta) then  
     if (i > nel) then  
        i = i - 1  
     else  
        sumw = sumw - ws(irk(ind(i)))          ! subtract overcounting
     end if
  else  
     if (i > nel) i = i - 1  
     write(9, 600) dtwo * (sumt - sumw)  
     sumt = sumw  
  end if

  emin = el(ind(i)) - dsmear   
  emax=el(ind(i)) + 0.5*dsmear      
  e0= el(ind(i))
  inv_smear = done / max (dsmear, 1.0d-6)  
  !  
  ! ensure emin and emax are proper bounds. ie. the solution lies between
  !
  nmax=1
100   continue
  !
  ! FMID = FUNC(emax) in Numerical Recipes.
  !
  Z1=dzero
  do i=1,nel 
    X = (emax - el(i))/dsmear

    if(ismear.eq.1) then
      Z1 = Z1 + ws(irk(i))*myderfc(-x) 
    elseif(ismear.eq.2) then
      Z1 = Z1 + ws(irk(i))*FERMID(-X)
    elseif(ismear.eq.3) then
      Z1 = Z1 + ws(irk(i))*DELTHM(X)
    elseif(ismear.eq.4) then
      Z1 = Z1 + ws(irk(i))*SPLINE(-X)
    elseif(ismear.eq.5) then
      Z1 = Z1 + ws(irk(i))*POSHM(X)
    elseif(ismear.eq.6) then
      Z1 = Z1 + ws(irk(i))*POSHM2(X)
    endif
  end do
  Z1=dhalf*Z1
  
  FMID= Z1-sumt  
  !
  ! F = FUNC(emin)
  !
  Z1=dzero
  do i=1,nel 
    X = (emin - el(i))/dsmear

    if(ismear.eq.1) then
      Z1 = Z1 + ws(irk(i))*myderfc(-x) 
    elseif(ismear.eq.2) then
      Z1 = Z1 + ws(irk(i))*FERMID(-X)
    elseif(ismear.eq.3) then
      Z1 = Z1 + ws(irk(i))*DELTHM(X)
    elseif(ismear.eq.4) then
      Z1 = Z1 + ws(irk(i))*SPLINE(-X)
    elseif(ismear.eq.5) then
      Z1 = Z1 + ws(irk(i))*POSHM(X)
    elseif(ismear.eq.6) then
      Z1 = Z1 + ws(irk(i))*POSHM2(X)
    endif
  end do
  Z1=dhalf*Z1
  
  F= Z1-sumt

  if(F*FMID .gt. dzero) then
    if (nmax.ge.10000) then
      write(9,*) 'ERROR: NO FERMI ENERGY FOUND ',emin,emax,F,FMID,nmax
      write(9,*) ' '
      write(9,*) 'IS THE ELECTRONIC TEMPERATURE TOO SMALL ? ',dsmear
      write(*,*) ' '
      call mystop
    else
      nmax=nmax+1
      emin=e0-real(nmax,dp)*dsmear
      emax=e0+(real(nmax,dp)-0.5)*dsmear
      goto 100
    end if
  end if

  if(F .le. dzero) then
    RTBIS = emin
    DX = emax - emin
  else
    RTBIS = emax
    DX = emin - emax
  end if
  !
  ! now that we have proper bounds, do bisection 
  !
  do J = 1, jmax
    DX = DX * 0.5D0
    XMID = RTBIS + DX
    !
    ! FMID=FUNC(XMID)
    !  
    Z1=dzero
    do i=1,nel 
      X = (xmid - el(i))/dsmear

      if(ismear.eq.1) then
        Z1 = Z1 + ws(irk(i))*myderfc(-x) 
      elseif(ismear.eq.2) then
        Z1 = Z1 + ws(irk(i))*FERMID(-X)
      elseif(ismear.eq.3) then
        Z1 = Z1 + ws(irk(i))*DELTHM(X)
      elseif(ismear.eq.4) then
        Z1 = Z1 + ws(irk(i))*SPLINE(-X)
      elseif(ismear.eq.5) then
        Z1 = Z1 + ws(irk(i))*POSHM(X)
      elseif(ismear.eq.6) then
        Z1 = Z1 + ws(irk(i))*POSHM2(X)
      endif
    end do
    Z1=dhalf*Z1

    FMID= Z1-sumt

    if(FMID .le. 0.D0) RTBIS=XMID

    if(abs(DX) .lt. xacc .or. abs(FMID) .lt. 1d-12) then
       goto 10
    end if
  
  end do

  write(9,*) 'CANNOT BISECT FOREVER, CAN I ?',FMID,DX,xacc
  CALL EXIT

10  energs%efermi(1) = RTBIS

  energs%efermi(2) = energs%efermi(1)  
  energs%esmear=dzero
  end if  ! if iafmag

  !
  !     ---------- we have the fermi energy. now we determine the fracs --
  !
  if (imode == 0) then  
     energs%efermi(1) = energs%inputfermi
     energs%efermi(2) = energs%inputfermi
  end if
  !
  ! ORGANIZE OUTPUT FOR GW CODE
  !
  if (iand(pw_params%output(1), 1048576) == 1048576 .or.  &
       (iand(pw_params%output(1), 524288) == 524288)) then 
    if ( bands%gw_mid_energy .lt. 1d4) then 
      xmid=bands%gw_mid_energy
    else
      xmid=energs%efermi(1)
    end if



    nel = 0 
    
    if ( bands%gw_low_energy .lt. 1d4) then       

      allocate(gw_num_band(nrk,ispn))  
      allocate(ind_tmp(bands%max))  


      do is = 1, ispn  
        do i = 1, nrk  
          icnt=0
          do j = 1, bands%nband(i, is)  
           nel=nel+1      
            if (el(nel) .ge. bands%gw_low_energy .and. &
               el(nel) .le. bands%gw_high_energy) then
                icnt=icnt+1  
                bands%gwout_index(icnt,i,is)=j
            end if
          end do
          gw_num_band(i,is)=icnt
          if( gw_num_band(i,is) .eq. 0) then
            if (el(nel) .lt.  bands%gw_low_energy) then
               bands%gwout_index(1,i,is)= bands%nband(i, is)+1
            else
               bands%gwout_index(1,i,is)= 1
            end if
          end if
        end do
      end do  

   
      gw_band_max=0
      do is = 1, ispn  
        do i = 1, nrk
          if(gw_num_band(i,is) .gt. gw_band_max) gw_band_max=gw_num_band(i,is)
        end do
      end do

      if(mod(gw_band_max,2) .eq. 1) gw_band_max=gw_band_max+1
      bands%num_gwout=gw_band_max/2

      do is = 1, ispn  
        do i = 1, nrk  
          icnt=0
          if(gw_num_band(i,is) .lt. gw_band_max) then
            k=min(gw_band_max-gw_num_band(i,is), bands%gwout_index(1,i,is)-1)
            do j = 1,gw_num_band(i,is)
                ind_tmp(j+k) =  bands%gwout_index(j,i,is) 
            end do
            do j = 1, k 
                ind_tmp(j)=  j+bands%gwout_index(1,i,is)-k-1
            end do
            gw_num_band(i,is)= gw_num_band(i,is)+k
            do j = 1, gw_num_band(i,is)
              bands%gwout_index(j,i,is)= ind_tmp(j)
            end do            
          end if
          if(gw_num_band(i,is) .lt. gw_band_max) then
            k=gw_band_max-gw_num_band(i,is)                            
             do j = gw_num_band(i,is)+1,gw_band_max
              bands%gwout_index(j,i,is)=j
             end do
             gw_num_band(i,is)= gw_num_band(i,is)+k          
          end if     
        end do
        end do


      deallocate(gw_num_band)
       deallocate(ind_tmp) 

    else



      do is = 1, ispn  
       do i = 1, nrk  

         if (bands%num_gwout .eq. -1 .or. &
                      2*bands%num_gwout .ge. bands%nband(i, is) ) then
           do j = 1, bands%nband(i, is)  
             bands%gwout_index(j,i,is)=j          
           end do


         else 


           icnt=0
           do j = 1, bands%nband(i, is)  
             nel=nel+1

             if (el(nel) .lt. xmid) then
               if(j .le.bands%num_gwout) then  
                 bands%gwout_index(j,i,is)=j   
                 icnt=icnt+1
               else
                 do k=1,bands%num_gwout           
                   bands%gwout_index(k,i,is)= bands%gwout_index(k,i,is)+1
                 end do
               end if
             else if(el(nel) .ge. xmid .and. icnt .lt. 2*bands%num_gwout) then
               icnt=icnt+1
!                bands%gwout_index(bands%num_gwout+icnt,i,is)=j              
               bands%gwout_index(icnt,i,is)=j  
             end if

           end do

           if (icnt .lt. 2*bands%num_gwout) then
             do k=1,icnt
               bands%gwout_index(k,i,is)=bands%gwout_index(k,i,is)-&
                (2*bands%num_gwout-icnt)
             end do
             do k=icnt+1,2*bands%num_gwout 
               bands%gwout_index(k,i,is)=bands%nband(i, is)-&
                 (2*bands%num_gwout-k)
             end do
           end if

         end if  
  



       end do
      end do
    end if

  end if  ! end data organization for gwout


!    do is = 1, ispn  
!      do i = 1, nrk  
!        do j = 1, 2*bands%num_gwout
!!        do j = 1,bands%nband(i, is)  
!          write(9,*) bands%gwout_index(j,i,is),is,i,j
!        end do
!       end do
!    end do

  tnel = nel                    ! total nel
  nel = 0  
  dnel = dzero  
  do is = 1, ispn  
    do i = 1, nrk  
      !
      bands%ifmax(i, is) = 0  



      if(iafmag) then
  

      do j = 1, bands%nband1
        !
        nel = nel + 1  
        if (ioccup == 0) then  

          if (iafmag) then  
            desm = (el(nel) - energs%efermi(is)) * inv_smear  
          else  
            desm = (el(nel) - energs%efermi(1)) * inv_smear  
          end if

          if(ismear.eq.1) then
            bands%occup(j, i, is) = kpoints%w(i) * myderfc(desm) 
            energs%esmear=energs%esmear &
              -dsmear/sqrpi*dtwo*kpoints%w(i)*dhalf*exp(-desm*desm)
          elseif(ismear.eq.2) then
            bands%occup(j, i, is) =kpoints%w(i) * FERMID(desm)
            FI=FERMID(desm)/dtwo
            if(abs(FI) .gt. 1.E-06 .and. abs(FI-done) .gt. 1.E-06) then
              energs%esmear=energs%esmear+dtwo*dsmear* kpoints%w(i) * &
                                 (FI*log(FI)+(done-FI)*log(done-FI))
            end if
          elseif(ismear.eq.3) then
            bands%occup(j, i, is) = kpoints%w(i) * DELTHM(-desm)
            energs%esmear=energs%esmear+dsmear/dtwo*kpoints%w(i) &
              *(dtwo*desm*desm-1)*exp(-desm*desm)/sqrpi
          elseif(ismear.eq.4) then
            bands%occup(j, i, is) = kpoints%w(i) * SPLINE(desm)
            desm=abs(desm)
            zeta=eesh*desm*exp(-(desm+sq2i)**2) +piesqq * myderfc(desm+sq2i)
            energs%esmear=energs%esmear-dtwo*dsmear*kpoints%w(i) *zeta
          elseif(ismear.eq.5) then
            bands%occup(j, i, is) = kpoints%w(i) * POSHM(-desm)
            a=-0.5634
            !        a=-0.8165
            energs%esmear=energs%esmear-dsmear/2.0*kpoints%w(i)  &
            ! NOTE g's are all intended to be normalized to 1 
            ! this following line is -2*int_minf^x [t*g(t)]dt
            *(2.0*a*desm**3-2.0*desm*desm+1 )*exp(-desm*desm)/sqrpi
          elseif(ismear.eq.6) then
            bands%occup(j, i, is) = kpoints%w(i) * POSHM2(-desm)
            energs%esmear=energs%esmear-dtwo*dsmear/2.0* kpoints%w(i) &
            ! NOTE g's are all intended to be normalized to 1 !
            ! this following line is -2*int_minf^x [t*g(t)]dt
            *(1-sqr2*desm)*exp(-(desm-1/sqr2)**2)/sqrpi
          endif
          if (bands%occup(j, i, is) > 1.0d-12) bands%ifmax(i, is) = j   


        else  

! charged system. P. Zhang

          if (real(j, dp) <= (crys%ztot +crys%net_charge + 1.0d-6) / 2.d0) then  
            bands%ifmax(i, is) = j  
            bands%occup(j, i, is) = dtwo * kpoints%w(i)  
          else  
            bands%occup(j, i, is) = dzero 
          end if
        end if
        dnel = dnel + bands%occup(j, i, is) * dhalf

      end do



      else

      do j = 1, bands%nband(i, is)  
        !
        nel = nel + 1  
        if (ioccup == 0) then  

          if (iafmag) then  
            desm = (el(nel) - energs%efermi(is)) * inv_smear  
          else  
            desm = (el(nel) - energs%efermi(1)) * inv_smear  
          end if

          if(ismear.eq.1) then
            bands%occup(j, i, is) = kpoints%w(i) * myderfc(desm) 
            energs%esmear=energs%esmear &
              -dsmear/sqrpi*dtwo*kpoints%w(i)*dhalf*exp(-desm*desm)
          elseif(ismear.eq.2) then
            bands%occup(j, i, is) =kpoints%w(i) * FERMID(desm)
            FI=FERMID(desm)/dtwo
            if(abs(FI) .gt. 1.E-06 .and. abs(FI-done) .gt. 1.E-06) then
              energs%esmear=energs%esmear+dtwo*dsmear* kpoints%w(i) * &
                                 (FI*log(FI)+(done-FI)*log(done-FI))
            end if
          elseif(ismear.eq.3) then
            bands%occup(j, i, is) = kpoints%w(i) * DELTHM(-desm)
            energs%esmear=energs%esmear+dsmear/dtwo*kpoints%w(i) &
              *(dtwo*desm*desm-1)*exp(-desm*desm)/sqrpi
          elseif(ismear.eq.4) then
            bands%occup(j, i, is) = kpoints%w(i) * SPLINE(desm)
            desm=abs(desm)
            zeta=eesh*desm*exp(-(desm+sq2i)**2) +piesqq * myderfc(desm+sq2i)
            energs%esmear=energs%esmear-dtwo*dsmear*kpoints%w(i) *zeta
          elseif(ismear.eq.5) then
            bands%occup(j, i, is) = kpoints%w(i) * POSHM(-desm)
            a=-0.5634
            !        a=-0.8165
            energs%esmear=energs%esmear-dsmear/2.0*kpoints%w(i)  &
            ! NOTE g's are all intended to be normalized to 1 
            ! this following line is -2*int_minf^x [t*g(t)]dt
            *(2.0*a*desm**3-2.0*desm*desm+1 )*exp(-desm*desm)/sqrpi
          elseif(ismear.eq.6) then
            bands%occup(j, i, is) = kpoints%w(i) * POSHM2(-desm)
            energs%esmear=energs%esmear-dtwo*dsmear/2.0* kpoints%w(i) &
            ! NOTE g's are all intended to be normalized to 1 !
            ! this following line is -2*int_minf^x [t*g(t)]dt
            *(1-sqr2*desm)*exp(-(desm-1/sqr2)**2)/sqrpi
          endif
          if (bands%occup(j, i, is) > 1.0d-12) bands%ifmax(i, is) = j   

        else  

! charged system. P. Zhang

          if (real(j, dp) <= (crys%ztot +crys%net_charge + 1.0d-6) / 2.d0) then  
            bands%ifmax(i, is) = j  
            bands%occup(j, i, is) = dtwo * kpoints%w(i)  
          else  
            bands%occup(j, i, is) = dzero 
          end if
        end if
        dnel = dnel + bands%occup(j, i, is) * dhalf

      end do
      end if


    end do
  end do
  energs%esmear = energs%esmear / real(bands%nspin, dp)
  !
  !
  !
  if (pw_params%occupy == 3) then
     do is = 1, ispn
        do i = 1, nrk
           j = bands%ifmax(i, is)
           if (j < bands%nband(i, is)) then
              bands%occup(j, i, is) = dhalf * (bands%occup(j, i, is) + &
                   bands%occup(j + 1, i, is))
              bands%occup(j + 1, i, is) = bands%occup(j, i, is)
              bands%ifmax(i, is) = j + 1
           else
              write(9, *) 'flevel: error assigning excited occupancies!'
           end if
        end do
     end do
  end if

  !
  !
  !
  if (iafmag) then  
     if (ipr > 1) write(9, 9021) energs%efermi(1), energs%efermi(2)
9021 format   ( /' The Fermi level (includes vxc(0))', &
          &        ' for spin(up,down) is ' // &
          &        '(',f14.8,',',f14.8,') Ry')

!        do i = 1, nrk  
!           write(9, *) i 
!           do j = 1, bands%nband1
!              write(9, *) (bands%energy(j, i, is) + energs%vxc0(is),is=1,2)
!              write(9, *) (bands%occup(j, i, is) *dhalf / kpoints%w(i),is=1,2)
!           end do
!        end do

  else  
     if (ipr > 1) write(9, 9020) energs%efermi(1)  
  end if
     if (ipr > 1) then
       if (ispn .eq. 1) then
         write(9, 9022) energs%vxc0(1)  
       else
         write(9, 9024) energs%vxc0(1),energs%vxc0(2)   
       end if       
9020 format   ( /' The Fermi level (includes vxc(0)) is',f14.8,' Ry')  
9022 format   ( /' vxc0=',f14.8,' Ry')  
9024 format   ( /' vxc0(up)=',f14.8,' Ry','    vxc0(down)=',f14.8,' Ry' )  
  end if


  if (ioccup == 1) write(9, 530)  


  if (ipr >= 1) then  
     do is = 1, ispn  
        if (ispn == 2 .and. is == 1) write(9, '(1x,a)') 'Spin up'  
        if (ispn == 2 .and. is == 2) write(9, '(1x,a)') 'Spin down'  
        do i = 1, nrk  
           write(9, 500) i 
           do n = 1, bands%nband(i, is), 7  
              write(9, 510) (bands%energy(j, i, is) + energs%vxc0(is), &
                   j = n, min(bands%nband(i, is), n + 6))
              write(9, 515) (bands%occup(j, i, is) *dhalf / kpoints%w(i), &
                   j = n, min(bands%nband(i, is), n + 6))
           end do
        end do
     end do

  end if


  deallocate(ind)
  deallocate(irk)
  deallocate (el)
  deallocate (ws)

  return  
  !
500 format(' k-point:',i3)  
510 format(7f11.4)  
  ! 510  format(7f22.16)
515 format(3x,7(' (',f6.4,')  '))  
  ! 515  format(3x,7(' (',f22.16,')  '))

530 format('Fractional occupation numbers (fixed occupation!):')  
600 format(' *** WARNING: NOT ENOUGH BANDS ARE COMPUTED.',/ &
       &     ' *** THE CRYSTAL IS POSITIVELY CHARGED:',f12.6)

end subroutine flevel


!C
!C =========================================================================
!C
      FUNCTION NMERFC(XX)
!C
!C     COMPLEMENTARY ERROR FUNCTION
!C     FROM THE SANDIA MATHEMATICAL PROGRAM LIBRARY
!C
!C     XMAX IS THE VALUE BEYOND WHICH ERFC(X) = 0 .
!C     IT IS COMPUTED AS SQRT(LOG(RMIN)), WHERE RMIN IS THE
!C     SMALLEST REAL NUMBER REPRESENTABLE ON THE MACHINE.
!C     IBM VALUE: (THE INTRINSIC ERFC COULD ALSO BE USED)
!C     PARAMETER ( XMAX = 13.4 )
!C     VAX VALUE: (XMAX = 9.3)
!C -----------------------------------
!C     12-Mar-90  Obtained from B. Hammer
!C     12-MAR-90  Changed to single precision at the end XW
!C                also XX1
      include 'use.h' 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( XMAX = 9.3D0)

      real(dp) nmerfc
!C

      DIMENSION P1(4),Q1(4),P2(6),Q2(6),P3(4),Q3(4)
!C
      DATA P1 /242.6679552305318D0 , 21.97926161829415D0 , &
      6.996383488619136D0 , -3.560984370181539D-2/
      DATA Q1 /215.0588758698612D0 , 91.16490540451490D0, &
      15.08279763040779D0 , 1.0D0/
      DATA P2 /22.898992851659D0 , 26.094746956075D0 , &
       14.571898596926D0 , 4.2677201070898D0 , &
      0.56437160686381D0 , -6.0858151959688D-6/
      DATA Q2 /22.898985749891D0 , 51.933570687552D0 , &
      50.273202863803D0 , 26.288795758761D0 , &
      7.5688482293618D0 , 1.0D0/
      DATA P3 /-1.21308276389978D-2 , -0.1199039552681460D0 , &
      -0.243911029488626D0 , -3.24319519277746D-2/
      DATA Q3 /4.30026643452770D-2 , 0.489552441961437D0 , &
      1.43771227937118D0 , 1.0D0/
!C     1/SQRT(PI)
      DATA SQPI /0.564189583547756D0/
!C
!C----------------------------------------------------------------------
      IF (XX .GT.  XMAX)    GOTO 330
      IF (XX .LT. -XMAX)    GOTO 320
      X = ABS(XX)
      X2 = X*X
      IF (X .GT. 4.0D0)     GOTO 300
      IF (X .GT. 0.46875D0) GOTO 200
!C
!C     -46875 < X < 0.46875
      NMERFC = X*(P1(1) + X2*(P1(2) + X2*(P1(3) + X2*P1(4))))
      NMERFC = NMERFC/(Q1(1) + X2*(Q1(2) + X2*(Q1(3) + X2*Q1(4))))
      IF (XX .LT. 0.0) NMERFC = - NMERFC
      NMERFC = 1.0D0 - NMERFC
      GOTO 9999
!C
200   NMERFC = EXP( -X2)*(P2(1) + X*(P2(2) + X*(P2(3) + X*(P2(4) + &
      X*(P2(5) + X*P2(6))))))
      NMERFC = NMERFC/(Q2(1) + X*(Q2(2) + X*(Q2(3) + X*(Q2(4) + X*(Q2(5) + &
      X*Q2(6))))))
      IF (XX .LE. 0.0) NMERFC = 2.0D0 - NMERFC
      GOTO 9999
!C
300   XI2 = 1.0D0/X2
      NMERFC = XI2*(P3(1) + XI2*(P3(2) + XI2*(P3(3) + XI2*P3(4))))/ &
      (Q3(1) + XI2*(Q3(2) + XI2*(Q3(3) + XI2*Q3(4))))
      NMERFC = EXP( -X2)*(SQPI + NMERFC)/X
      IF (XX .LT. 0.0) NMERFC = 2.0D0 - NMERFC
      GOTO 9999
!C
320   NMERFC = 2.0D0
      GOTO 9999
330   NMERFC = 0.0D0
!C
9999  RETURN
      END


!================================================================
      FUNCTION FERMID(XX)
      double precision::XX, FERMID
!
      IF(XX .GT. 30.D0) THEN
        FERMID=0.D0
      ELSEIF(XX .LT. -30.D0) THEN
        FERMID=2.D0
      ELSE
        FERMID=2.D0/(1.D0+EXP(XX))
      ENDIF
!
      RETURN
      END
!================================================================
      FUNCTION DELTHM(XX)
!
      double precision::XX,DELTHM
      pi=3.14159265358979
      IF(XX .GT. 10.D0) THEN
        DELTHM=2.D0
      ELSEIF(XX .LT. -10.D0) THEN
        DELTHM=0.D0
      ELSE
        DELTHM=(2.D0-ERFC(XX))+XX*EXP(-XX*XX)/SQRT(PI)
      ENDIF
!
      RETURN
      END
!================================================================
      FUNCTION SPLINE(X)
      include 'use.h'
      implicit none 
      real(dp) X,spline

      real(dp) eesqh,sq2i,fx     
  
      eesqh=sqrt(exp(done))*0.5
      sq2i=sqrt(dtwo)*0.5
      if (x.ge.0.0) then
        fx=eesqh*exp(-(x+sq2i)**2)
      else
        fx=done-eesqh*exp(-(x-sq2i)**2)
      endif
      spline=dtwo*fx
!
      return
      end
!================================================================
      FUNCTION POSHM(X)
      double precision::X,POSHM
!
! NOTE g's are all intended to be normalized to 1 !
! function = 2 * int_minf^x [g(t)] dt
!
      pi=3.141592653589793238
      a=-0.5634
!     a=-0.8165
      IF(X .GT. 10.D0) THEN
        POSHM=2.D0
      ELSEIF(X .LT. -10.D0) THEN
        POSHM=0.D0
      ELSE
        POSHM=(2.D0-ERFC(X))+(-2.0*a*x*x+2*x+a)*EXP(-X*X)/SQRT(PI)/2.0
      ENDIF
!
      RETURN
      END
!================================================================
      FUNCTION POSHM2(X)
!
! NOTE g's are all intended to be normalized to 1 !
! function = 2 * int_minf^x [g(t)] dt
!
      double precision::X,POSHM2

      pi=3.141592653589793238
      IF(X .GT. 10.D0) THEN
        POSHM2=2.D0
      ELSEIF(X .LT. -10.D0) THEN
        POSHM2=0.D0
      ELSE
        POSHM2=(2.D0-ERFC(X-1./sqrt(2.)))+ &
       sqrt(2.)*exp(-x*x+sqrt(2.)*x-0.5)/sqrt(pi)
      ENDIF
!
      RETURN
      END
!===========================WELLWELLWEHAVEFINISHED.
