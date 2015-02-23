!-*-F90-*- 
subroutine setup_blochl_operator_q(ham,blop,crys)

  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'

  type (hamiltonian)     :: ham
  type (blochl_operator) :: blop
  type (crystal)         :: crys

  !    
  !     Set up the nonlocal part of the Hamiltonian for a given k-point.
  !     Kleinman-Bylander Pseudopotential.
  !
  !     all it does is setting up the arrays ham%vnloc and ham%xnorm
  !     at the kpoint specified by ham%gspace%rk
  !     

  !     --------------------- local variables ------------------------------

  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4),fftn(4),ffth(4),igv3,irod,iord,igs

  real(dp) :: gv(3)

  integer :: ind,i,k,kk,l,ni,ni1,lmmin,lmmax,lm 

  real(dp) :: qi,qk(3),qcar(3),qinv,fi,xni,sum,fpi
  real(dp) :: zero,one,two,three,half,eps,vq,xi,xdum,xa,xb


  parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0)
  parameter ( three = 3.0d0, half = 0.5d0)
  parameter ( eps = 1.0d-8)

  complex(dp) :: flm(9)
  complex(dp), allocatable :: st(:) ! work array for structure factors

  fpi = 16.d0 * atan(1.d0)
  !
  !     -----------------------------------------------------------------------

  allocate(st(crys%mxdatm))

  fftn(1:3) = ham%gspace%fftsize(1:3)
  fftn(4)   = ham%gspace%fftsize(3)
  ffth(:)   = fftn(:)/2
  !
  !      starts loop over g-vectors in small gspace 
  !

  igs = 0 
  do iord=1, ham%gspace%lorder    ! loop through x/y gspace 
     irod=ham%gspace%order(1,iord)
     igv(1)= irod/ ham%gspace%fftsize(2)
     igv(2)= mod(irod, ham%gspace%fftsize(2))
     igv(3)= ham%gspace%order(2,iord)
     igv(4)= ham%gspace%order(3,iord)
     igv(:)= mod(igv(:)+ffth(:),fftn(:))-ffth(:) 
     gv(1:2)=dble(igv(1:2))
     do igv3=igv(3),igv(4)     ! loop over z axis
        gv(3)=dble(igv3)
        igs = igs + 1
        qi    = sqrt(ham%gspace%ekin(igs)) ! get kinetic energy from array
        qk(:) = ham%gspace%rk(:) +gv(:)
        !
        !           cartesian of k+G
        !
        qcar(1) = crys%bvec(1,1)*qk(1)+crys%bvec(1,2)*qk(2)+&
             crys%bvec(1,3)*qk(3)
        qcar(2) = crys%bvec(2,1)*qk(1)+crys%bvec(2,2)*qk(2)+&
             crys%bvec(2,3)*qk(3)
        qcar(3) = crys%bvec(3,1)*qk(1)+crys%bvec(3,2)*qk(2)+&
             crys%bvec(3,3)*qk(3)
        !
        !           ----- compute angular functions
        !     
        !     ** NOTE The projectors are the c.c. of the standard defns.
        !     
        flm(1) = one
        if(qi .gt. eps) then
           qinv = one/qi
           flm(2) = +sqrt(1.5d0)*cmplx(+qcar(1),+qcar(2))*qinv
           flm(3) = +sqrt(3.0d0)*cmplx(+qcar(3),+0.000d0)*qinv
           flm(4) = -sqrt(1.5d0)*cmplx(+qcar(1),-qcar(2))*qinv
           qinv   = qinv*qinv
           flm(5) = +sqrt(15.d0/8.d0)*cmplx(+qcar(1),+qcar(2))**2*qinv
           flm(6) = +sqrt(15.d0/2.d0)*cmplx(qcar(1),+qcar(2))*qcar(3)*qinv
           flm(7) = +sqrt(5.d0/4.d0)*(3.d0*qcar(3)**2*qinv-1.d0)
           flm(8) = -sqrt(15.d0/2.d0)*cmplx(qcar(1),-qcar(2))*qcar(3)*qinv
           flm(9) = +sqrt(15.d0/8.d0)*cmplx(+qcar(1),-qcar(2))**2*qinv
        else
           do lm=2,9
              flm(lm) = zero
           end do
        endif
        !          
        ! starts loop over second index
        !

        ind = 0
        do k=1,crys%ntype

           do kk=1,crys%natom(k) ! compute complex phase factor
              fi = gv(1)*crys%rat(1,kk,k) + &
                   gv(2)*crys%rat(2,kk,k) + &
                   gv(3)*crys%rat(3,kk,k)
              st(kk) = exp(cmplx(0.d0, fi,dp))
           end do

           ! loop over lo is replaced by loop over l
           ! and l quantum number is given  explicitely by referencing blop%lo

           do l=1,blop%nlprjmx              ! loop over angular momenta l
              if(blop%nkb(l,k) .ne. 0) then ! found pot of that ang. mom
                 xni = qi/blop%delqnl(k) + two ! interpolate potential
                 !
                 ! cubic spline interpolation
                 !
                 vq=zero
                 ni = xni 
                 if(ni .le. 2) ni = 3
                 ni1 = ni+1

                 if(ni .lt. blop%nqnl(k)) then
                    xa = dble(ni1)-xni
                    xb = xni-dble(ni)
                    vq = xa*blop%vkb(ni,l,k) + xb*blop%vkb(ni1,l,k) + &
                         ((xa**3-xa)*blop%d2vkbdq2(ni-2,l,k) + &
                         (xb**3-xb)*blop%d2vkbdq2(ni1-2,l,k))/6.d0
                 endif


                 lmmin = blop%lo(l,k)*blop%lo(l,k) + 1
                 lmmax = (blop%lo(l,k)+1)*(blop%lo(l,k)+1)
                 !     
                 do lm=lmmin,lmmax ! loop over m quantum number
                    do kk=1,crys%natom(k)
                       ind = ind + 1
!                       blop%opnloc_q%data(igs,ind,1) = &
!                            conjg(st(kk)*flm(lm)*vq)*sqrt(fpi/crys%vcell)
                    end do
                 end do
              endif
           end do


        end do              ! end of loop over atomic types

     end do                 ! end of loop over 3rd dimension

  end do                    ! end of loop over small gspace

  deallocate(st)

  return

end subroutine setup_blochl_operator_q

