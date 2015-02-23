!
! Module to hold keyword list. This must be updated as
! new keywords are brought into existence.
!
! The 'label' is the label as used in calling the esdf routines
! 'typ' defines the type, with the following syntax. It is 3 characters
! long. 
! The first indicates:
!  I - integer
!  S - single
!  D - double
!  P - physical
!  T - string (text)
!  E - defined (exists)
!  L - boolean (logical)
!  B - block
! The second is always a colon (:)
! The third indicates the "level" of the keyword
!  B - Basic
!  I - Intermediate
!  E - Expert
!  D - Dummy
!
! 'dscrpt' is a description of the variable. It should contain a (short) title
! enclosed between *! ... !*, and then a more detailed description of the 
! variable.
!
Module esdf_key

  Implicit None

  Integer, Parameter :: numkw=200
  
  Character(len=80)   :: kw_label(numkw)
  Character(len=3)    :: kw_typ(numkw)
  Character(len=3000) :: kw_dscrpt(numkw)

  ! Now define the keywords
  Data kw_label(1) /"coordinates"/
  Data kw_typ(1)   /"B:B"/
  Data kw_dscrpt(1) /"*! Atomic data and coordinates !*"/

  Data kw_label(2) /"latticevecs"/
  Data kw_typ(2)   /"B:B"/
  Data kw_dscrpt(2) /"*! Lattice vectors and cell volume !*"/

  Data kw_label(3) /"coordinates_absolute"/
  Data kw_typ(3)   /"E:B"/
  Data kw_dscrpt(3)/"*! Absolute coordinates !*"/

  Data kw_label(4) /"number_kpoints"/
  Data kw_typ(4)   /"I:E"/
  Data kw_dscrpt(4)/"*! Number of k-points !*"/

  Data kw_label(5) /"k_grid"/
  Data kw_typ(5)   /"T:B"/
  Data kw_dscrpt(5)/"*! k-point grid !*"/

  Data kw_label(6) /"k_grid_shift"/
  Data kw_typ(6)  /"T:B"/
  Data kw_dscrpt(6)/"*! k-point grid shift !*"/

  Data kw_label(7) /"gw_shift"/
  Data kw_typ(7)   /"T:E"/
  Data kw_dscrpt(7)/"*! GW shift of k-point grid !*"/

  Data kw_label(8) /"gaussian_smearing"/
  Data kw_typ(8)  /"P:I"/
  Data kw_dscrpt(8)/"*! Gaussian smearing width in eV !*"/

  Data kw_label(9) /"accuracy_diag"/
  Data kw_typ(9)   /"D:E"/
  Data kw_dscrpt(9)/"*! Diagonalization accuracy !*"/

  Data kw_label(10) /"diagsafety"/
  Data kw_typ(10)  /"P:E"/
  Data kw_dscrpt(10)/"*! Diagonalization safety in eV !*"/

  Data kw_label(11) /"randomize_diag_start_guess"/
  Data kw_typ(11)   /"D:E"/
  Data kw_dscrpt(11)/"*! Seed for random initialization !*"/

  Data kw_label(12) /"bandgap"/
  Data kw_typ(12) /"P:E"/
  Data kw_dscrpt(12)/"*! Estimated band gap in eV !*"/

  Data kw_label(13) /"energy_cutoff"/
  Data kw_typ(13) /"P:B"/
  Data kw_dscrpt(13)/"*! Plane-wave cutoff in Ry !*"/

  Data kw_label(14) /"submatrix_energy_cutoff"/
  Data kw_typ(14) /"P:B"/
  Data kw_dscrpt(14)/"*! Cutoff for submatrix diagonalization !*"/

  Data kw_label(15) /"mixing_energy_cutoff"/
  Data kw_typ(15) /"P:B"/
  Data kw_dscrpt(15)/"*! Energy cutoff for potential mixing !*"/

  Data kw_label(16) /"modify_kinetic_energy"/
  Data kw_typ(16) /"T:E"/
  Data kw_dscrpt(16)/"*! Modified kinetic energy functional parameters !*"/

  Data kw_label(17) /"number_radial_gridpoints"/
  Data kw_typ(17) /"I:E"/
  Data kw_dscrpt(17)/"*! Number of radial grid points !*"/

  Data kw_label(18) /"number_bands_fft"/
  Data kw_typ(18) /"I:E"/
  Data kw_dscrpt(18)/"*! Number of bands to FFT simultaneously !*"/

  Data kw_label(19) /"polished_energy_cutoff"/
  Data kw_typ(19) /"P:E"/
  Data kw_dscrpt(19)/"*! Cutoff for stress correction !*"/

  Data kw_label(20) /"number_of_spins"/
  Data kw_typ(20) /"I:B"/
  Data kw_dscrpt(20)/"*! Number of spins !*"/

  Data kw_label(21) /"max_iter_diag"/
  Data kw_typ(21) /"I:I"/
  Data kw_dscrpt(21)/"*! Now minimum number of diagonalization iterations !*"/

  Data kw_label(22) /"max_iter_scfloop"/
  Data kw_typ(22) /"I:I"/
  Data kw_dscrpt(22)/"*! Maximum number of SCF iterations !*"/

  Data kw_label(23) /"nproc_pvm"/
  Data kw_typ(23) /"I:D"/
  Data kw_dscrpt(23)/"*! Defunct !*"/

  Data kw_label(24) /"screening_type"/
  Data kw_typ(24) /"T:I"/
  Data kw_dscrpt(24)/"*! Screening type !*"/

  Data kw_label(25) /"output_flags"/
  Data kw_typ(25) /"T:E"/
  Data kw_dscrpt(25)/"*! Flags for additional output !*"/

  Data kw_label(26) /"no"/
  Data kw_typ(26) /"T:E"/
  Data kw_dscrpt(26)/"*! Switch off restoration of inversion symmetry !*"/

  Data kw_label(27) /"checkpoint"/
  Data kw_typ(27) /"T:E"/
  Data kw_dscrpt(27)/"*! NMR checkpoint !*"/

  Data kw_label(28) /"checkpoint_emin"/
  Data kw_typ(28) /"I:I"/
  Data kw_dscrpt(28)/"*! Checkpoint frequency !*"/

  Data kw_label(29) /"exchange_correlation"/
  Data kw_typ(29) /"T:B"/
  Data kw_dscrpt(29)/"*! Exchange-correlation functional !*"/

  Data kw_label(30) /"input_flags"/
  Data kw_typ(30) /"T:E"/
  Data kw_dscrpt(30)/"*! Additional input !*"/

  Data kw_label(31) /"optimize"/
  Data kw_typ(31) /"T:E"/
  Data kw_dscrpt(31)/"*! Optimization flags !*"/

  Data kw_label(32) /"occupy_levels"/
  Data kw_typ(32) /"T:I"/
  Data kw_dscrpt(32)/"*! Occupation scheme !*"/

  Data kw_label(33) /"fermi_level"/
  Data kw_typ(33) /"P:I"/
  Data kw_dscrpt(33)/"*! Fermi level in Ry !*"/

  Data kw_label(34) /"mixing_parameter"/
  Data kw_typ(34) /"D:I"/
  Data kw_dscrpt(34)/"*! Mixing parameter !*"/

  Data kw_label(35) /"linear_mixing"/
  Data kw_typ(35) /"T:I"/
  Data kw_dscrpt(35)/"*! Linear mixing parameters !*"/

  Data kw_label(36) /"potential_convergence_criterion"/
  Data kw_typ(36) /"D:E"/
  Data kw_dscrpt(36)/"*! Potential convergence criterion !*"/

  Data kw_label(37) /"nmr_q"/
  Data kw_typ(37) /"D:I"/
  Data kw_dscrpt(37)/"*! NMR wavevector !*"/

  Data kw_label(38) /"nmr_g0mask"/
  Data kw_typ(38) /"T:I"/
  Data kw_dscrpt(38)/"*! NMR parameters !*"/

  Data kw_label(39) /"visualize"/
  Data kw_typ(39) /"T:I"/
  Data kw_dscrpt(39)/"*! Visualization parameters !*"/

  Data kw_label(40) /"number_symmetry_ops"/
  Data kw_typ(40) /"I:E"/
  Data kw_dscrpt(40)/"*! Number of symmetry operations !*"/

  Data kw_label(41) /"number_bands"/
  Data kw_typ(41) /"I:B"/
  Data kw_dscrpt(41)/"*! Number of bands !*"/

  Data kw_label(42) /"magnetic_field_kvec"/
  Data kw_typ(42) /"T:E"/
  Data kw_dscrpt(42)/"*! Magnetic field wavevector !*"/

  Data kw_label(43) /"magnetic_field_strength"/
  Data kw_typ(43) /"D:E"/
  Data kw_dscrpt(43)/"*! Magnetic field strength !*"/

  Data kw_label(44) /"electric_potential_kvec"/
  Data kw_typ(44) /"T:E"/
  Data kw_dscrpt(44)/"*! Electric field wavevector!*"/

  Data kw_label(45) /"electric_potential_strength"/
  Data kw_typ(45) /"D:I"/
  Data kw_dscrpt(45)/"*! Electric field strength !*"/

  Data kw_label(46) /"plot_wave_function"/
  Data kw_typ(46) /"B:E"/
  Data kw_dscrpt(46)/"*! Wavefunction plot !*"/

  Data kw_label(47) /"tensor_bra"/
  Data kw_typ(47) /"T:E"/
  Data kw_dscrpt(47)/"*! Tensor bra !*"/

  Data kw_label(48) /"tensor_ket"/
  Data kw_typ(48) /"T:E"/
  Data kw_dscrpt(48)/"*! Tensor ket !*"/

  Data kw_label(49) /"line_plot"/
  Data kw_typ(49)/"B:E"/
  Data kw_dscrpt(49)/"*! Bandstructure lines !*"/

  Data kw_label(50) /"pw_jobs"/
  Data kw_typ(50)/"B:B"/
  Data kw_dscrpt(50)/"*! Plane-wave job list !*"/

  Data kw_label(51) /"energy_window"/
  Data kw_typ(51) /"B:I"/
  Data kw_dscrpt(51)/"*! Energy window list !*"/

  Data kw_label(52) /"stay_put"/
  Data kw_typ(52) /"B:B"/
  Data kw_dscrpt(52)/"*! Relaxation constraint list !*"/

  Data kw_label(53) /"starting_distortion"/
  Data kw_typ(53)/"T:I"/
  Data kw_dscrpt(53)/"*! Initial distortion !*"/

  Data kw_label(54) /"decouple_coordinates"/
  Data kw_typ(54)/"T:I"/
  Data kw_dscrpt(54)/"*! Decouple coordinates !*"/

  Data kw_label(55) /"estimated_bulk_modulus"/
  Data kw_typ(55) /"P:I"/
  Data kw_dscrpt(55)/"*! Estimated bulk modulus !*"/

  Data kw_label(56) /"relax_eps_factor"/
  Data kw_typ(56) /"D:I"/
  Data kw_dscrpt(56)/"*! Relaxation parameter !*"/

  Data kw_label(57) /"estimated_optical_phonon_frequency"/
  Data kw_typ(57)/"P:I"/
  Data kw_dscrpt(57)/"*! Estimated optical phonon frequency !*"/

  Data kw_label(58) /"relax_coord_factor"/
  Data kw_typ(58) /"D:I"/
  Data kw_dscrpt(58)/"*! Relaxation parameter !*"/

  Data kw_label(59) /"relax_accuracy"/
  Data kw_typ(59) /"D:B"/
  Data kw_dscrpt(59)/"*! Relaxation convergence criterion !*"/

  Data kw_label(60) /"relax_pressure"/
  Data kw_typ(60) /"T:I"/
  Data kw_dscrpt(60)/"*! Relaxation pressure !*"/

  Data kw_label(61) /"job"/
  Data kw_typ(61)/"T:B"/
  Data kw_dscrpt(61)/"*! Job !*"/

  Data kw_label(62) /"lambda_limit"/
  Data kw_typ(62) /"D:E"/
  Data kw_dscrpt(62)/"*! Relaxation parameter !*"/

  Data kw_label(63) /"relax_what"/
  Data kw_typ(63) /"T:B"/
  Data kw_dscrpt(63)/"*! Relaxation parameter !*"/

  Data kw_label(64) /"relax_how"/
  Data kw_typ(64) /"T:B"/
  Data kw_dscrpt(64)/"*! Relaxation method !*"/

  Data kw_label(65) /"relax_max_iter"/
  Data kw_typ(65) /"I:B"/
  Data kw_dscrpt(65)/"*! Maximum number of relaxation iterations !*"/

  Data kw_label(66) /"number_alt_bands"/
  Data kw_typ(66) /"I:B"/
  Data kw_dscrpt(66)/"*! Number of bands for band structure !*"/

  Data kw_label(67) /"number_bands_tddft"/
  Data kw_typ(67) /"I:B"/
  Data kw_dscrpt(67)/"*! Number of bands for TDDFT !*"/

  Data kw_label(68) /"nmrkpts"/
  Data kw_typ(68) /"I:I"/
  Data kw_dscrpt(68) /"*! Number of k-points to do in one go, for NMR !*"/

  Data kw_label(69) /"bandstructure"/
  Data kw_typ(69) /"B:I"/
  Data kw_dscrpt(69) /"*! Bandstructure information !*"/

  Data kw_label(70) /"diagonalization_method"/
  Data kw_typ(70) /"T:E"/
  Data kw_dscrpt(70) /"*! method for electronic minimization !*"/

  Data kw_label(71) /"mix_method"/
  Data kw_typ(71) /"T:E"/
  Data kw_dscrpt(71) /"*! for choosing different potential mixing methods !*"/

  Data kw_label(72) /"accuracy_diag_per_band"/
  Data kw_typ(72)   /"E:B"/
  Data kw_dscrpt(72)/"*! the accuracy needed per band !*"/
 
   Data kw_label(73) /"fixed_vector"/
   Data kw_typ(73)   /"T:E"/
   Data kw_dscrpt(73) /"*! the fixed vector !*"/
 
   Data kw_label(74) /"other_vector_in_fixed_plane"/
   Data kw_typ(74)   /"T:E"/
   Data kw_dscrpt(74) /"*! the other vector in the fixed plane !*"/
 
   Data kw_label(75) /"deformation"/
   Data kw_typ(75)   /"B:I"/
   Data kw_dscrpt(75) /"*! the initial deformation tensor !*"/

   Data kw_label(76) /"pp_format"/
   Data kw_typ(76)   /"I:B"/
   Data kw_dscrpt(76) /"*!format for pseupotenital file!*"/

   Data kw_label(77) /"pp_data"/
   Data kw_typ(77) /"T:B"/
   Data kw_dscrpt(77) /"*! local PP and atomic occ !*"/

   Data kw_label(78) /"pseudopotential"/
   Data kw_typ(78)   /"B:B"/
   Data kw_dscrpt(78) /"*! pseudopotential data !*"/

   Data kw_label(79) /"PTF_num_init_PK"/
   Data kw_typ(79)   /"I:E"/
   Data kw_dscrpt(79) /"*! number PK steps to do before PTF starts !*"/   

   Data kw_label(80) /"PTF_max_iter_cg"/
   Data kw_typ(80)   /"I:E"/
   Data kw_dscrpt(80) /"*! number cg steps to solve the TFW equations !*"/

  Data kw_label(81) /"energy_convergence_criterion"/
  Data kw_typ(81) /"D:E"/
  Data kw_dscrpt(81)/"*! energy convergence criterion !*"/

  Data kw_label(82) /"slice_plot"/
  Data kw_typ(82)/"B:E"/
  Data kw_dscrpt(82)/"*!  Slice plotting coordinates !*"/
                                                              
  Data kw_label(83) /"bond_length"/
  Data kw_typ(83)/"D:E"/
  Data kw_dscrpt(83)/"*! length of the bond for angle calculation !*"/

  Data kw_label(84) /"NLPP_rspace"/
  Data kw_typ(84) /"L:E"/
  Data kw_dscrpt(84)/"*! for real space multiplication of Vnl*Psi !*"/

  Data kw_label(85) /"NLPP_rcut"/
  Data kw_typ(85) /"T:E"/
  Data kw_dscrpt(85)/"*! for real space multiplication of Vnl*Psi !*"/

  Data kw_label(86) /"super_cell"/
  Data kw_typ(86) /"T:E"/
  Data kw_dscrpt(86)/"*! add layers of atoms along lattice vecs!*"/

  Data kw_label(87) /"super_cell_vac"/
  Data kw_typ(87) /"T:E"/
  Data kw_dscrpt(87)/"*! add layers of vacuum along lattice vecs  !*"/

  Data kw_label(88) /"NLPP_rspace_force"/
  Data kw_typ(88) /"L:E"/
  Data kw_dscrpt(88)/"*! add layers of vacuum along lattice vecs  !*"/

  Data kw_label(89) /"max_iter_diag_band"/
  Data kw_typ(89) /"I:I"/
  Data kw_dscrpt(89)/"*! Maximum number of diagonalization iterations !*"/

  Data kw_label(90) /"smearing_method"/
  Data kw_typ(90) /"I:I"/
  Data kw_dscrpt(90)/"*! method for determining occupations !*"/

  Data kw_label(91) /" io_scf "/
  Data kw_typ(91) /"L:E"/
  Data kw_dscrpt(91)/"*! restrict io in SCFLOOP !*"/

  Data kw_label(92) /"smearing_energy"/
  Data kw_typ(92)  /"P:I"/
  Data kw_dscrpt(92)/"*! occupation smearing width in eV !*"/

  Data kw_label(93) /"relax_method"/
  Data kw_typ(93)  /"I:E"/
  Data kw_dscrpt(93)/"*! method to relax atomic positions !*"/

  Data kw_label(94) /"MD_time_step"/
  Data kw_typ(94)  /"D:E"/
  Data kw_dscrpt(94)/"*! time step for MD in femtoseconds!*"/

  Data kw_label(95) /"MD_temp"/
  Data kw_typ(95)  /"D:E"/
  Data kw_dscrpt(95)/"*! temp for MD in K!*"/

  Data kw_label(96) /"MD_max_iter"/
  Data kw_typ(96)  /"I:B"/
  Data kw_dscrpt(96)/"*! number of molecular dynamics steps !*"/

  Data kw_label(97) /"ensemble_type"/
  Data kw_typ(97)  /"I:I"/
  Data kw_dscrpt(97)/"*! type of enemble for MD !*"/

  Data kw_label(98) /"extrapolation_method"/
  Data kw_typ(98)  /"I:E"/
  Data kw_dscrpt(98)/"*! extrapolation method!*"/

  Data kw_label(99) /"MD_Q_mass"/
  Data kw_typ(99)  /"D:I"/
  Data kw_dscrpt(99)/"*! mass of extended heat bath!*"/

  Data kw_label(100) /"number_of_alpha"/
  Data kw_typ(100) /"I:B"/
  Data kw_dscrpt(100)/"*! Number of alpha !*"/

  Data kw_label(101) /"number_of_beta"/
  Data kw_typ(101) /"I:B"/
  Data kw_dscrpt(101)/"*! Number of beta !*"/

  Data kw_label(102) /"number_of_bands_alpha"/
  Data kw_typ(102) /"I:B"/
  Data kw_dscrpt(102)/"*! Number of bands for alpha (spin up) !*"/

  Data kw_label(103) /"number_of_bands_beta"/
  Data kw_typ(103) /"I:B"/
  Data kw_dscrpt(103)/"*! Number of bands for beta (spin down)!*"/

  Data kw_label(104) /"num_gwout"/
  Data kw_typ(104) /"I:B"/
  Data kw_dscrpt(104)/"*! number of bands/2 for gwout !*"/

  Data kw_label(105) /"gwout_mid_energy"/
  Data kw_typ(105) /"D:I"/
  Data kw_dscrpt(105)/"*! energy level midpoint for gwout output!*"/

  Data kw_label(106) /"gwout_low_energy"/
  Data kw_typ(106) /"D:I"/
  Data kw_dscrpt(106)/"*! low energy level for gwout output!*"/

  Data kw_label(107) /"gwout_high_energy"/
  Data kw_typ(107) /"D:I"/
  Data kw_dscrpt(107)/"*! high energy level for gwout output!*"/

  Data kw_label(108) /"iter_diag_add"/
  Data kw_typ(108) /"I:E"/
  Data kw_dscrpt(108)/"*! additional iterations for diag*"/

  Data kw_label(109) /"iter_diag_add_metal"/
  Data kw_typ(109) /"I:E"/
  Data kw_dscrpt(109)/"*! additional iterations for Grassmann_metal diag*"/

  Data kw_label(110) /"min_iter_diag"/
  Data kw_typ(110) /"I:I"/
  Data kw_dscrpt(110)/"*! Now minimum number of diagonalization iterations !*"/

  Data kw_label(111) /"init_dir_econv"/
  Data kw_typ(111)   /"D:E"/
  Data kw_dscrpt(111)/"*! enery accuracy before H updated in direct method !*"/

  Data kw_label(112) /"lsda_u"/
  Data kw_typ(112) /"I:E"/
  Data kw_dscrpt(112)/"*! perform LSDA+U calculations*"/

  Data kw_label(113) /"wannier"/
  Data kw_typ(113) /"I:E"/
  Data kw_dscrpt(113)/"*! output for wannier program*"/

  Data kw_label(114) /"gw_band_reordering"/
  Data kw_typ(114) /"I:E"/
  Data kw_dscrpt(114)/"*! re-ordering band index in the wavefunction output*"/

  Data kw_label(115) /"starting_guess"/
  Data kw_typ(115) /"I:E"/
  Data kw_dscrpt(115)/"*! use atomic wavefunction as starting guess*"/

!Weiwei Gao 2015
  Data kw_label(116) /"gw_sum"/
  Data kw_typ(116) /"I:E"/
  Data kw_dscrpt(116)/"*! Ouput part of wavefunctions for GW energy integration !*"/

End Module esdf_key







