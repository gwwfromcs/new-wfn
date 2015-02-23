!-*-Fortran-*-
!     -------- scalapack global parameters -----------------
!     
!
!     parameters for the old version of scalapack
!
!      INTEGER   CSRC_, CTXT_, DLEN_, LLD_, MB_, M_, NB_, N_, RSRC_
!      PARAMETER ( DLEN_ = 8, 
!     $     M_    = 1, 
!     $     N_    = 2, 
!     $     MB_   = 3, 
!     $     NB_   = 4,
!     $     RSRC_ = 5, 
!     $     CSRC_ = 6, 
!     $     CTXT_ = 7, 
!     $     LLD_  = 8 )
!
!     for the new version of scalapack
!
  integer, parameter :: dlen_ = 9, dtype_ = 1, &
          ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, &
          rsrc_ = 7, csrc_ = 8, lld_ = 9
!
!     ---------- scalapack helper functions ------------------
!
  integer, external :: iceil, numroc
