
do n = 1, instances  !{  necessary but unused feature
include "bio.f90"

   do k = 1, nk  !{
   do j = jsc, jec  !{ 
    do i = isc, iec  !{

     t_prog(ind_dic)%th_tendency(i,j,k) = 0.
     t_prog(ind_alk)%th_tendency(i,j,k) = 0.    
     t_prog(ind_o2)%th_tendency(i,j,k)  = 0.
     t_prog(ind_po4)%th_tendency(i,j,k) = 0.
     t_prog(ind_no3)%th_tendency(i,j,k) = 0.

     enddo  !} i
    enddo  !} j
    enddo  !} k
enddo  !} n

