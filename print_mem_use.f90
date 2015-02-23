subroutine print_mem_use(pot_gspace,ffts,crys,k_gspace,pspot,&
              pw_params,sub_gspace,ham,kpoints,bands)
  include 'use.h'  
  IMPLICIT NONE
  !
  !     INPUT:
  !     -----
  !

  type(crystal), intent(in) :: crys       ! crystal structure
  type(fft_struc), intent(in) :: ffts     ! information for the FFT
  type(pseudo_potential), intent(in) :: &
       pspot                              ! all info about the pseudopotentials
 type(parallel_gspace), intent(in) :: &
       pot_gspace, k_gspace(*)
  type(pw_parameter), intent(in) :: &
       pw_params  
  type(parallel_gspace) :: &
       sub_gspace                         ! for the starting guess
  type(hamiltonian), intent(in) :: ham    ! hamiltonian at a given k-point
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints
  type(band), intent(in) :: bands      ! eigenvalues, occupation numbers etc
  !
  !     --------------- local variables -------------------------
  !
  integer maxlength,ntype,nfn,n,nb,nrpoc,nprow,npcol,ngmax,nbc,nbr,ngc,ngr,&
          nanl,nproc,ngdiffmax,itmaxpul

  real(dp) tot_mem,start_mem,diag_mem,fac_diag,fac_mem,pp_mem

  itmaxpul=20
  maxlength=k_gspace(1)%length
  ntype=crys%ntype
  nfn = pw_params%nbandsfft
  nanl =  ham%vnloc%nvecs

  nproc = pot_gspace%nproc
  n = sub_gspace%length  

  call layout_scalapack(n, nb, nproc, nprow, npcol)

  nbc = n / (nb *npcol)
  nbc = nbc + 1
  nbr = n / (nb *nprow)
  nbr = nbr + 1

  ngmax = nb * max(nbc, nbr)

  ngc=nb * nbc
  ngr=nb * nbr

  ngdiffmax = ngc * ngr   

  if (iand(pw_params%optimize,8) == 8 .or. &
             iand(pw_params%optimize,64) == 64) then ! direct minimization
     fac_diag=real(kpoints%nrk*crys%nspin,dp)
  else
     fac_diag=done
  end if
  
  fac_mem= 16d0/1024d0**2   ! 16 bytes for complex divided by 1Mbyte

  tot_mem=  (ffts%r_size*nfn*6 &
             + (4 + 2*(ntype+crys%nspin))*pot_gspace%length &
             + 3*pot_gspace%lorder +3* pot_gspace%fftsize(1)  &
             +         (7+ntype)*pot_gspace%length +   pot_gspace%totlength & 
             + crys%nspin*pot_gspace%length &
             + ( 3*k_gspace(1)%lorder + 3* k_gspace(1)%fftsize(1)+ &
                                (4+ntype)*k_gspace(1)%length ) * kpoints%nrk & 
             + maxlength*bands%max*kpoints%nrk*crys%nspin &
             + pot_gspace%length*(crys%nspin*(itmaxpul+1)*2) ) * fac_mem &
             + 12d0  

  if (.not. pspot%NL_rspace(1)) then
    pp_mem = (k_gspace(1)%length*pspot%nanl + pspot%nanl*bands%max)*fac_mem 
  else
    pp_mem= (13*pspot%mrb2_matom_node + pspot%NLPP_rsp_pts) * fac_mem

  end if

  start_mem =(ngdiffmax + ngr*ngc + (ngc*nanl+ngr*nanl)) * fac_mem


  if(pw_params%diag_method == 'Grassmann' .or. &
            pw_params%diag_method == 'Grassmann_metal'  )then

    diag_mem = (maxlength*bands%max*5 +  bands%max**2*2) * fac_mem

    write(9,135)  12d0,&
        ffts%r_size*nfn * 6 * fac_mem,&
        (4 + 2*(ntype+crys%nspin))*pot_gspace%length * fac_mem, &  
        ( 3*pot_gspace%lorder +3* pot_gspace%fftsize(1) + &
              (7+ntype)*pot_gspace%length  +  pot_gspace%totlength) * fac_mem,&
        crys%nspin*pot_gspace%length * fac_mem,&
        (3*k_gspace(1)%lorder +3* k_gspace(1)%fftsize(1)+ &
                          (4+ntype)*k_gspace(1)%length)*kpoints%nrk*fac_mem, & 
        maxlength*bands%max*kpoints%nrk*crys%nspin * fac_mem, & 
        pot_gspace%length*(crys%nspin*(itmaxpul+1)*2) * fac_mem 

    if (.not. pspot%NL_rspace(1)) then
      write(9,136) &
          (k_gspace(1)%length*pspot%nanl + pspot%nanl*bands%max) * fac_mem  
    else
      write(9,137) & 
        (13*pspot%mrb2_matom_node + pspot%NLPP_rsp_pts) * fac_mem
    end if

    write(9,138) ngdiffmax * fac_mem, &
        ngr*ngc*fac_mem, &
        (ngc*nanl+ngr*nanl)*fac_mem, &
        fac_diag * maxlength * bands%max * 5 * fac_mem, &
        fac_diag * bands%max**2 * 2 * fac_mem

  else

    diag_mem = (maxlength*bands%max*5 + bands%max**2*6) * fac_mem

    write(9,135)  12d0,&
        ffts%r_size*nfn*6 * fac_mem,&
        (4 + 2*(ntype+crys%nspin))*pot_gspace%length * fac_mem, &  
        ( 3*pot_gspace%lorder +3* pot_gspace%fftsize(1) + &
            (7+ntype)*pot_gspace%length  +   pot_gspace%totlength) * fac_mem,&
        crys%nspin*pot_gspace%length* fac_mem,&
        (3*k_gspace(1)%lorder +3* k_gspace(1)%fftsize(1)+ &
                          (4+ntype)*k_gspace(1)%length)*kpoints%nrk*fac_mem, & 
        maxlength*bands%max*kpoints%nrk*crys%nspin * fac_mem, & 
        pot_gspace%length*(crys%nspin*(itmaxpul+1)*2) * fac_mem

    if (.not. pspot%NL_rspace(1)) then
      write(9,136) &
          (k_gspace(1)%length*pspot%nanl + pspot%nanl*bands%max) * fac_mem  
    else
      write(9,137) & 
        (13*pspot%mrb2_matom_node + pspot%NLPP_rsp_pts) * fac_mem
    end if

     write(9,138)  ngdiffmax * fac_mem, &
        ngr*ngc * fac_mem, &
        (ngc*nanl+ngr*nanl) * fac_mem, &
        fac_diag*maxlength* bands%max*5* fac_mem, &
        fac_diag*bands%max**2*6* fac_mem
                    
  end if

  write(9,139)  tot_mem + max(start_mem,diag_mem) + pp_mem 
  write(9,*) 
  call myflush(9)

135 format(/' MEMORY USAGE:',/1x,27('-'), &
       &           /20x,'STRUCTURE',15x,'Mb/proc'&
       &     /,11X,'                       =',          &
       &     /,11X,' est. executable size  =',2X,F12.6,&
       &     /,11X,' FFT structure         =',2X,F12.6, &
       &     /,11X,' pot decl read_pseudo  =',2X,F12.6, &
       &     /,11X,' pot_gspace            =',2X,F12.6, &
       &     /,11X,' vout                  =',2X,F12.6, &
       &     /,11X,' k_gspace              =',2X,F12.6, &
       &     /,11X,' wavefunction- all k   =',2X,F12.6, &
       &     /,11X,' mixing in scfloop     =',2X,F12.6)
136 format(    11X,' k-space NON-local     =',2X,F12.6)
137 format(    11X,' real space NON-local  =',2X,F12.6)
138 format(        ' +  start guess memory'  &
       &     /,11X,' glist-to map vloc     =',2X,F12.6, &
       &     /,11X,' dense ham             =',2X,F12.6, &
       &     /,11X,' vnlsub map Vnon-loc   =',2X,F12.6, &
       &     / ' or +  diagonalization memory'  &
       &     /,11X,' wf copy in CG         =',2X,F12.6, &
       &     /,11X,' occupied space in CG  =',2X,F12.6 )

139  format(/' Estimated total memory usage: ',F10.3,' MBytes / proc')

end subroutine print_mem_use
