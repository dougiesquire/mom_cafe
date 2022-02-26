subroutine bio_v0a (Thickness)  
! rjm  null model

type(ocean_thickness_type), intent(in) :: Thickness


!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------
integer :: i
integer :: j
integer :: k
integer :: n
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second

!
! =====================================================================
!     begin executable code
! =====================================================================


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

return
end subroutine bio_v0a

