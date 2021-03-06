! m4undef.m4
!
! resets various m4 commands to have a prefix m4_
! this can be achieved using option -P but only on some architectures
! any file to be preprocessed with m4 should use this set of macro definitions
!
! David Prendergast, June 6, 2006



! fft_macros.m4
!
!     fft_aux_space(aux,naux,aux2,naux2,dec1,dec2)
!     fft_local_init(N,A,LDA,LOT,DIREC)
!     fft_perm_init(N,AUX)
!
!     fft_multiple_backward(N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!     fft_multiple_forward (N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!
!     fft_backward(N,A,LOTSTRIDE)
!     fft_forward (N,A,LOTSTRIDE)
!
!     fft_convol(N,A,B,POS)
!
!     fft_multiple_scale(N,A,LOTSTRIDE,LOT,LOOPDUMMY,SCALE)
!
!     timeget (t0)
!     timediff(t1,t0)
!
!
!     fft_fcblock     defines blocking factor for fast convolute
!     fft_local_free  release memory
!

! Originally services.c
!
!   The following routines are called from the fortran code to do
!   various system jobs.
!
! Converted to Fortran90 by Peter Haynes March 2000

subroutine frepack(packsize, out, pack2, outsize, in, pack1, insize)

  use constants
  implicit none

  integer, intent(in) :: packsize, outsize, insize
  integer, intent(in) :: pack1(packsize), pack2(packsize)
  complex(dp), intent(out) :: out(outsize)
  complex(dp), intent(in) :: in(insize)
  integer :: i, j, k, k1, j1, p2






  do i = 1, packsize
     j = pack2(i)
     k = pack1(i)
     out(j) = in(k)
  end do




end subroutine frepack


subroutine fftlt2(dir, chunk, chunk_size, rspace, rspace_size, t1buf, &
     t1buf_size, nnx, nny, ldx, ldy)

  use constants
  implicit none
  include 'flibcalls.ph'

  integer, intent(in) :: dir, nnx, nny, ldx, ldy, rspace_size, t1buf_size, &
       chunk_size
  integer, intent(in) :: chunk(chunk_size)
  complex(dp), intent(out) :: rspace(rspace_size)
  complex(dp), intent(in) :: t1buf(t1buf_size)
  integer :: bailout
  integer :: ichunk, srod, irod, sy, slen, swidth, drod, dx, dlen, ix

  bailout = 1

  if (dir < 0) then                   ! repacking for the inverse transform

     ichunk = 1
     do
        if (bailout < 0) exit

        srod = chunk(ichunk)
        sy = chunk(ichunk + 1) + 1

        slen = chunk(ichunk + 2) - 1
        bailout = chunk(ichunk + 3)
        swidth = abs(bailout)

        drod = chunk(ichunk + 4)
        dx = chunk(ichunk + 5)
        !        dlen = chunk(ichunk + 6) - 1

        do irod = srod, srod + swidth - 1
           ix = mod(dx + (irod - srod), nnx) + 1
           call mzcopy(slen + 1, t1buf(irod * ldy + sy), 1, &
                rspace(drod * ldx + ix), ldx )
           !             rspace(drod*ldx+ix:(drod+dlen)*ldx+ix:ldx) = &
           !                  t1buf(irod*ldy+sy:irod*ldy+sy+slen)
        end do
        ichunk = ichunk + 8
     end do
  else 

     ichunk = 1
     do
        if (bailout < 0) exit

        srod = chunk(ichunk)
        sy = chunk(ichunk + 1) + 1

        slen = chunk(ichunk + 2) - 1
        bailout = chunk(ichunk + 3)
        swidth = abs(bailout)

        drod = chunk(ichunk + 4)
        dx = chunk(ichunk + 5)
        !        dlen = chunk(ichunk + 6) - 1
        do irod = srod, srod + swidth - 1
           ix = mod(dx + (irod - srod), nnx) + 1
           call mzcopy(slen + 1, rspace(drod * ldx + ix), ldx, &
                t1buf(irod * ldy + sy), 1)
           !              t1buf(irod*ldy+sy:irod*ldy+sy+slen) = &
           !                   rspace(drod*ldx+ix:(drod+dlen)*ldx+ix:ldx)
        end do
        ichunk = ichunk + 8
     end do
  end if

end subroutine fftlt2


subroutine fftunfoldc(folded, unfolded, ldz, gsorder, gslorder, gslength, nfn)

  use constants
  implicit none

  integer, intent(in) :: ldz, gslorder, gslength, nfn
  integer, intent(in) :: gsorder(3 * gslorder)
  complex(dp), intent(in) :: folded(gslength * nfn)
  complex(dp), intent(out) :: unfolded(ldz * gslorder * nfn)

  !     DESCRIPTION:
  !     -----------
  !
  !     Transfers the rods from the packed gspace arrangement onto
  !     a linear grid of length gs%fftsize(3), such that a subsequent
  !     1d FFT can be performed

  integer :: ii, iord, istart, iend, idiff, iid, nzdiff,iuffz, &
       ifn, iunfoldskip, ifoldskip

  unfolded = zzero

  do ifn = 0, nfn - 1
     iunfoldskip = ifn * gslorder * ldz
     ifoldskip = ifn * gslength

     ii = 0
     iuffz = 0

     do iord = 0, gslorder - 1
        istart = gsorder(3 * iord + 2)
        iend = gsorder(3 * iord + 3)
        idiff = iend - istart
        if (idiff >= 0) then  ! there is only one contiguous block to copy
           iid = ii + mod(idiff + ldz, ldz)
           unfolded(istart + iuffz + iunfoldskip + 1: &
                iend + iuffz + iunfoldskip + 1) = &
                folded(ii + ifoldskip + 1:ii + ifoldskip + idiff + 1)
        else  ! there are two blocks to copy
           nzdiff = ldz - istart
           iid = ii + nzdiff - 1
           unfolded(istart + iuffz + iunfoldskip + 1: &
                istart + iuffz + iunfoldskip + nzdiff) = &
                folded(ii + ifoldskip + 1:ii + ifoldskip + nzdiff)
           ii = iid + 1
           iid = ii + iend
           unfolded(iuffz + iunfoldskip + 1:iuffz + iunfoldskip + iend + 1) = &
                folded(ii + ifoldskip + 1:ii + ifoldskip + iend + 1)
        end if
        iuffz = iuffz + ldz
        ii = iid + 1
     end do
  end do

end subroutine fftunfoldc


! provide a fortran substitute for memory copy

subroutine fmemcpy(dest, source, num)

  implicit none

  integer, intent(in) :: num
  character, intent(out) :: dest(num)
  character, intent(in) :: source(num)

  dest = source

end subroutine fmemcpy


! clear double complex array

subroutine fcleardcarray(dest, num)

  use constants
  implicit none

  integer :: num
  complex(dp), intent(out) :: dest(num)

  dest = zzero

end subroutine fcleardcarray


!void mymalloc(int *size, void **poi)
!{
!  /* printf("mallocing %d bytes\n",*size); */
!  *poi=malloc(*size);
!  /* printf("malloced %d bytes\n",*size); */
!}


!void myfree(void **poi)
!{
!  free(*poi);
!}

!#ifdef RS6K
!
!#include <sys/time.h>
!
!void mygettimeofday(double *time) {
!struct timeval Tp;
!struct timezone Tzp;
!int idum;
!
!idum = gettimeofday(&Tp,&Tzp);
!*time = (double )Tp.tv_sec + 1e-6*(double) Tp.tv_usec;
!}
!
!#endif




!       this one is for debugging on the IBMs. it prints out statistics 
!	about current memory usage of the process. if you want to use 
!	this function, you have to use mymalloc and myfree to allocate 
!	and free memory.

!#ifdef RS6K
!void mymallinfo()
!{
!  FILE *memy;
!
!  struct mallinfo minfo;
!  minfo = mallinfo();
!
!  memy=fopen("mem_state","a");
!  fprintf(memy, "memory state:\n");
!  
!  fprintf(memy,"arena    = %d\n",minfo.arena);
!  fprintf(memy,"#ordblks = %d\n",minfo.ordblks);
!  fprintf(memy,"#smblks  = %d\n",minfo.smblks);
!  fprintf(memy,"#hblks   = %d\n",minfo.hblks);
!  fprintf(memy,"hblkhd   = %d\n",minfo.hblkhd);
!  fprintf(memy,"usmblks  = %d\n",minfo.usmblks);
!  fprintf(memy,"fsmblks  = %d\n",minfo.fsmblks);
!  fprintf(memy,"uordblks = %d\n",minfo.uordblks);
!  fprintf(memy,"fordblks = %d\n",minfo.fordblks);
!
!  fclose(memy);
!}
!
!#endif

!#ifdef CVX
!#include<math.h>
!double myderf(double *x)
!{
!return erf(*x);
!}
!double myderfc(double *x)
!{
!return erfc(*x);
!}
!#endif
