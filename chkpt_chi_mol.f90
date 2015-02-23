Subroutine chkpt_chi_mol(irk,j_bare,rho,r_size,myproc,ioutput,ireadwrite)

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  !
  !     1997 Bernd Pfrommer while at UCB
  !

  !     INPUT:
  !     -----

  Integer ::  ireadwrite ! if 1 then write, else read
  Integer :: irk         ! at which kpoint the last chkpointing was
  Integer :: myproc      ! which proc I am
  Integer :: r_size      ! size of rhom array
  Integer :: ioutput     ! flag for printing output

  Complex(dp) :: j_bare(r_size) ! the data to be saved
  Complex(dp) :: rho(r_size/9)

  !
  !     DESCRIPTION:
  !     -----------
  !
  !     For each processor, a checkpoint file is written/read which
  !     saves/retreives the necessary state for the chi calculation

  !     ---------------- local variables ----------------------------------

  Real(dp) :: t0, gimmetime

  External gimmetime

  Character*(255) :: filename    ! name of file to be written

  t0=gimmetime()

  Write(filename,'(''CHKPTCHI.'',i3.3)') myproc

  Open(75,file=filename,form='unformatted')

  If(ireadwrite.Eq.1) Then
     Write(75) irk,r_size
     Write(75) j_bare
     Write(75) rho
     Close(75)
     If(Iand(ioutput,8).Eq.8) Write(9,100) gimmetime()-t0
  Else
     Read(75) irk,r_size
     Read(75) j_bare
     Read(75) rho
     Close(75)
     If(Iand(ioutput,8).Eq.8) Write(9,110) gimmetime()-t0
  Endif

  Call myflush(9)

  Return

100 Format(' TIME TO WRITE CHKPTCHI:',f12.3)
110 Format(' TIME TO READ CHKPTCHI:',f12.3)

End Subroutine chkpt_chi_mol



