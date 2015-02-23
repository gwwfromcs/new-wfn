#
# Description:  Makefile for paratec executable
#
# Authors:	Bernd Pfrommer
#
# Version:	%Z%%M% %I% %E%
#
# 
#    This makefile looks intimidating at first, but really, it is
#    not that complicated. You should only have to  make changes 
#    until you hit the ==== line.
#
#    When you write a new routine, add it to the variable list
#    according to the classification, i.e. whether it is F90 or C,
#    m4-preprocessed or plain, optimized or unoptimized etc.
#    When you use special include files, you will have to make
#    sure that the dependencies are correct. Always test that 
#    touching the include file will cause a rebuild of the
#    file!

# ------------------------------------------------------------
#

#
# f90 sources that contain modules and need no preprocessing
#

F90SRCMOD   = constants.f90 iolib.f90 structures.f90 esdf_key_mod.f90        \
		esdf_mod.f90 lsdau_shared.f90

#
# f90 sources that contain modules and need preprocessing
#

F90PSRCMOD  = all_to_all.f90p flibcalls.f90p

#      
# f90 sources that don't require optimization or preprocessing
#

F90SRCUNOPT = write_eigvals.f90   cgemin.f90 \
		fermi_surface.f90 cgemin_gcg.f90 angular.f90 induced_current.f90  \
		randomize_startguess.f90 check_array.f90 findbound.f90      \
		switch.f90 string.f90 rest_inv.f90 group.f90 inpmat.f90     \
		symchk.f90 symm_ident.f90 symgen.f90 crstl.f90 dos.f90      \
		direct_emin.f90 lukman.f90 excorr.f90 pw91_excorr.f90       \
		gga91sub.f90 momentum_density.f90 chkpt_chi.f90 tpage.f90   \
		warn.f90 sort.f90 symmetrize_local.f90                      \
		magnetic_suscept_orig.f90                                   \
		cg_blocksolve.f90          \
		resolved_charge.f90 add_to_charge.f90 magnetic_kkterm.f90   \
		writegsdat.f90 writegsdat2.f90 add_to_j_bare.f90        \
		readgsdat.f90 vqmc.f90 gwreal.f90 gwout.f90  gwoutnewfmt.f90 \
		plot_bandstruc.f90 logfile.f90 ball_and_stick.f90           \
		reg_grid_spline.f90 neighbor.f90 force_stress.f90           \
		dxplot.f90 dxplot_wannier.f90 dxplot_field.f90  dx_latticevec.f90 \
		dx_reciprocal_lvec.f90 lineplot.f90  \
		symmetrize_tensor.f90 generate_potential_gspace.f90         \
		generate_k_gspaces.f90 print_gspace.f90 print_pwparam.f90   \
		 generate_kpoints.f90  scfloop.f90     \
		setup_local_potential.f90 layout_scalapack.f90 findroot.f90 \
		diagonalize_hamiltonian.f90  layout_occ_scalapack.f90     \
                diagonalize_hamiltonian_gcg.f90 diag_ham_metal_gcg.f90  \
                symmetrize_crystal.f90 generate_kpts.f90 tddft.f90          \
		magnetic_suscept_crys.f90 paramagnetic_correction_crys.f90 \
                diamagnetic_correction.f90 paramagnetic_correction_mol.f90  \
		magnetic_suscept_eqn3.f90 chkpt_chi_mol.f90 sliceplot.f90    \
		cg_blocksolve_mol.f90 magnetic_suscept_eqn8.f90              \
		chkpt_chi_mol2.f90 chkpt_j_crys.f90 chkpt_chi_crys.f90  \
		symmetrize_tensor_j.f90  direct_current_process.f90   \
		 print_mem_use.f90  lbfgs.f90 mol_dyn.f90 funeval.f90 relax.f90 
#
# f90 sources that require optimization, but no preprocessing
#

F90SRCOPT    = spec_tet.f90 ewald.f90 ewald_dipole.f90 put_gset.f90         \
		get_gset.f90 findvec.f90 nmr_shift.f90   \
		 calc_vnl_ave.f90 regspace.f90  compute_ekin_gspace.f90      \
		expand.f90  diagaux.f90 setup_nonlocal_potential.f90        \
		derivative_gspace.f90 wrapup_chi.f90 take_nonloc_deriv.f90  \
		 setup_nonlocal_derivative.f90 apply_ham_wavefn.f90  \
	        symmetrize_scalar_global.f90 update_hamiltonian.f90        \
                update_hamiltonian_metal.f90 update_hamiltonian_gcg.f90     \
		electric_field.f90 qualitykp.f90 radin_mod.f90              \
		bessfn_mod.f90 read_blochl_operator.f90 flevel.f90           \
		setup_blochl_operator.f90  init_rdm_wf.f90 velect.f90     \
		take_nonloc_deriv_kq_k.f90  forstressnloc_rsp.f90           \
		nmr_shift_new.f90 position_in_rspace.f90 apply_mom.f90      \
		apply_r2psi.f90 apply_r2rho.f90  w_linerq.f90 \
		charge.f90 charge2.f90  kinet.f90  read_pseudo.f90 maskr.f90  \
		cgemin_metal_gcg.f90 generate_gspace.f90  getwmask.f90 \
		etotal.f90 mixer.f90 mch_pulay.f90 thom_fer.f90 screen.f90  \
		forstressloc.f90 forstressnloc.f90 forstresssym.f90         \
		apply_rixdvnl.f90 j2chi_mol.f90 apply_dvnl.f90 lsdau.f90  \
		setzgpfa.f90  zgpfa2f.f90  zgpfa3f.f90  zgpfa5f.f90  zgpfa.f90 \
                project_pz.f90 project_dos.f90

#
# f90 sources that don't require optimization, but preprocessing
#

F90PSRCUNOPT = adjustfft.f90p fft_workspace.f90p angular_wavefn.f90p       \
		pwmain.f90p  data_gather.f90p            \
		isort_gather.f90p gspace_gather.f90p zdiag_scalapack.f90p  \
		rdiag_scalapack.f90p main.f90p inputread.f90p  \
		zdiag_occ_scalapack.f90p

#
# f90 sources that require optimization and preprocessing
#

F90PSRCOPT   = fourier_transform.f90p start_vectors.f90p                   \
		symmetrize_tensor_global.f90p setup_packinfo.f90p          \
		services.f90p  occsp_diag_scal.f90p occsp_sqrtinv_scal.f90p \
                lsdau_proj.f90p  lsdau_proj2.f90p wannier.f90p orthogonalization_PZ.f90p
              

#
# f90 include files that need to be generated by preprocessing the
#      corresponding .m4h file
#

F90INCFILES = flibcalls.ph

#
# C  sources, optimized and unoptimized separately
#

CSRCOPT   =
CSRCUNOPT = 
CINCFILES =

#
# Sources from the shared directory (they are shared with the benchmark)
#

#SHAREDSRC = fft_opcount.f90
# This is now included in the library $(SHAREDLIB) provided in $(DEPTH)/conf.mk

#
# ------------------------------------------------------------------
#            no need to touch stuff  below this line

DEPTH	= ../..
include $(DEPTH)/conf.mk
M4SRC = $(M4FFT)

#
#  the .f90=.o assignment converts all the filenames with .f90 into .o
#
#SHAREDOBJS	= $(F90SYSDEPSHAREDSRC:.f90=.o) $(SHAREDSRC:.f90=.o)

F90OBJSMOD	= $(F90SRCMOD:.f90=.o)
F90POBJSMOD	= $(F90PSRCMOD:.f90p=.o)
F90OBJSOPT	= $(F90SRCOPT:.f90=.o) $(FFTSRC:.f90=.o)
F90POBJSOPT	= $(F90PSRCOPT:.f90p=.o)

F90OBJSUNOPT	= $(F90SRCUNOPT:.f90=.o)
F90POBJSUNOPT	= $(F90PSRCUNOPT:.f90p=.o)

F90OBJS		= $(F90OBJSMOD) $(F90POBJSMOD) $(F90OBJSOPT) $(F90POBJSOPT) $(F90OBJSUNOPT) \
                  $(F90POBJSUNOPT) 

#COBJSOPT	= $(CSRCOPT:.c=.o)
#COBJSUNOPT	= $(CSRCUNOPT:.c=.o)
#COBJS 		= $(COBJSOPT) $(COBJSUNOPT)

#OBJS		= $(COBJS) $(F90OBJS) $(SHAREDOBJS)
#OBJS		= $(F90OBJS) $(SHAREDOBJS)
OBJS		= $(F90OBJS)
#$(COBJS) : $(CINCFILES)
$(F90OBJS) : $(F90INCFILES)


paramodules: $(F90INCFILES) $(F90INCFILES2)
	@$(MAKE) $(F90OBJSMOD) $(F90POBJSMOD) 

paratec: paraf parafopt parac paracopt
	$(F90) $(OBJS) $(LIBS) \
	-o $(EXECDIR)/$@$(FLAVOR)  $(FLINKSPECIAL) 


paraf:	paramodules 
	@$(MAKE) EXTFLAGS="$(F90DEBUGFLAGS)" ftarget

parafopt: paramodules 
	@$(MAKE) EXTFLAGS="$(F90OPTFLAGS)" ftargetopt

parac:  $(CINCFILES)
	@$(MAKE) EXTFLAGS="$(CDEBUGFLAGS)" ctarget

paracopt:  $(CINCFILES)
	@$(MAKE) EXTFLAGS="$(COPTFLAGS)" ctargetopt

ftarget:    $(F90OBJSUNOPT) $(F90POBJSUNOPT) 
ftargetopt: $(F90OBJSOPT) $(F90POBJSOPT)
ctarget:    $(COBJSUNOPT)
ctargetopt: $(COBJSOPT)

cleanparatec:
	@-rm -f *.o
	@-rm -f $(EXECDIR)/paratec$(FLAVOR)
	@-rm -f *.p.f90
	@-rm -f *.ph
	@-rm -f *.kmo
	@-rm -f *.mod
	@-rm -f *.M
	@-rm -f *.kif
	@-rm -f *.T
	@-rm -f *.ii
	@-rm -f core
	@-rm -f *~
	@-rm -f *.tar



