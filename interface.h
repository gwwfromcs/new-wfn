!-*-F90-*-
!
interface destroy
   module procedure destroy_parallel_gspace, &
        destroy_crystal, destroy_fftstruc, &
        destroy_pseudo_potential, &
        destroy_double_gspace_array, &
        destroy_complex_gspace_array, &
        destroy_complex_rspace_array, &
        destroy_force,destroy_kpoint, &
        destroy_band,destroy_hamiltonian, &
        destroy_dense_hamiltonian, &
        destroy_blochl_operator
end interface

interface create_gspace_array
   module procedure create_double_gspace_array, &
        create_complex_gspace_array
end interface

interface create_rspace_array
   module procedure create_complex_rspace_array
end interface

interface print_array
   module procedure print_complex_gspace_array
end interface

interface assignment(=)
   module procedure parallel_gspace_assignment, &
        energy_assignment, energy_assignment0
end interface

interface operator(-)
   module procedure energy_subtract
end interface




