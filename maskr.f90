     subroutine maskr(ri,amr)
!ccccc this program generate the mask to be
!ccccc used in real space implementation of the
!ccccc nonlocal potential
!cccccccccccccccccccccccccccccccccccc
  include 'use.h'
      use constants
  use flibcalls_module
  use all_to_all_module
      implicit NONE
  include 'interface.h'
      include 'flibcalls.ph'
  include 'all_to_all.h'
      integer nq
      parameter (nq=200)
      real(dp) B(nq)
      real(dp) H(nq,nq)
      real(dp) HL(nq*(nq+1)/2)
      real(dp) Z(nq,nq)
      real(dp) aux(3*nq)
      real(dp) E(nq)
      real(dp) ri(nq),amr(nq)
      real(dp) eta,x1,x2,yy,Zm
      integer i,j,nc,info
!ccccccccccccccccccccccccccccccccccccccccc
!cccc somehow, the H is still ill conditioned.
!cccc Its eigen value is not very good, but 
!cccc the lowest one is O.K

!cccccc xm=2.5 for 20x20x20 small box grid
!cccccc xm=3.0 for 24x24x24 small box grid
!      write(6,*) "input xm"
!      read(6,*) xm

 eta=real(15,dp)

  do i=1,nq
    x1=i*done/nq
    do j=1,i
      x2=j*done/nq

      if (i.eq.j) then
        yy=(pi*nq-eta)+dsin((x1+x2)*eta)/(x1+x2)
      else
        yy=-dsin((x1-x2)*eta)/(x1-x2)+ dsin((x1+x2)*eta)/(x1+x2)
      end if

      H(i,j)=yy
      H(j,i)=yy
    enddo
  enddo

!  nc=0
!  do i=1,nq
!    do j=1,i
!      nc=nc+1
!      HL(nc)=H(i,j)
!    enddo
!  enddo

!  call MDSPEV(21,HL,E,Z,nq,nq,aux,3*nq)

  call mdsyev('V','L',nq,H(1,1),nq,E(1),Z(1,1),3*nq,info)

!  write(6,*) "E(1)=",E(1)

!ccccc use this as the real space  mask

  do i=1,nq
    ri(i)=i*done/nq
    if(i.eq.1) Zm=H(i,1)/ri(i)
    amr(i)=H(i,1)/ri(i)/Zm
  enddo

!ccccc this is the mask in q space, used only for 
!ccccc testing
!      open(11,file="maskq")
!      rewind(11)
!
!      do j=1,300
!      r=1.d0*j/100.d0
!      s=0.d0
!      do i=1,nq
!      x=i*(done/nq)*eta
!      x2=x*r
!      s=s+dsin(x2)/x2*x**2*Z(i,1)/x
!      enddo
!      if(j.eq.1) sm=s
!      write(11,200) r,r*eta,s/sm,abs(s/sm)
!      enddo
!
!200   format(3(E15.8,1x))
!
!      close(11)


      return 
      end subroutine


