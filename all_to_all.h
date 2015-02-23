!-*-Fortran-*-
!
interface all_sum_all

   module procedure &
        all_sum_all_complex, all_sum_all_complex1, &
        all_sum_all_complex2, all_sum_all_complex3, &
        all_sum_all_complex4, &
        all_sum_all_integer, all_sum_all_integer1, &
        all_sum_all_double, all_sum_all_double1, &
        all_sum_all_double2, all_sum_all_double3

end interface

interface all_sum_all_sub

   module procedure &
        all_sum_all_complex_sub, all_sum_all_complex1_sub

end interface

interface all_max_all

   module procedure &
        all_max_all_double1, all_max_all_int1,all_max_all_int2

end interface

interface all_min_all

   module procedure &
        all_min_all_double1, all_min_all_int1, all_min_all_int2

end interface

interface my_broadcast
  
   module procedure &
        my_broadcast_int1, my_broadcast_int, my_broadcast_dcmplx, &
        my_broadcast_dcmplx3, my_broadcast_dcmplx2, &
        my_broadcast_double1, &
        my_broadcast_double, my_broadcast_double3, my_broadcast_char

end interface

interface my_broadcast_small

   module procedure &
        my_broadcast_int1_smcomm, my_broadcast_double_smcomm,&
        my_broadcast_dcmplx_small

end interface
