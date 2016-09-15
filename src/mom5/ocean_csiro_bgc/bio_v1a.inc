subroutine bio_v1a (Thickness)  
! rjm  Basic EP model
! csiro_bgc version 1 for mom4 

!    author:  r.j.matear  e-mail richard.matear@ml.csiro.au
!
!	t.m.roy updates from 22.8.03:
!* Employed a variable for the euphotic layer depth, n_eudepth 	
!* Global tracers averages are now calculated before rather than after  
! the calculation of the virtual fluxes
!* Added c14 to comply with OCMIP standards
!


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

integer :: nn, ntr_bgc

!chd a bunch ob biotic parameters
      real :: alpha    = 0.055/86400 !was 0.01
      real :: pmort    = 0.01/86400  !was 0.05 
      real :: presp    = 0.02/86400
      real :: pmax_0   = 0.8/86400
      real :: xnut_lim = 0.1 !was 1 
      real :: pfrac    = 0.1 
      real :: zmort    = 0.2/86400 !was 0.05 
      real :: graze    = 0.8/86400
      real :: zeff     = 0.8
      real :: zfrac    = 0.2
      real :: zhalf    = 0.75 !was 0.1
      real :: decay    = 0.1/86400 
      real :: pmax     = 0.8/86400 !pmax_0*1    **(temp-10.)*1e-1
      real :: w_detr   = 5/86400. 

!chd penetration of short wave radiation 
!chd jerlov water type I 

      real :: d1sol    = 0.35        
      real :: d2sol    = 23.         
      real :: rsol     = 0.58
        
      integer :: ts_npzd ! number of time steps within NPZD model
      integer :: tn      !number of timesteps within upstream det sink 
      real :: bion,biop,bioz,biod,bioo  !aux variable
      real :: dtsb       !timestep of npzd forward model
      real :: rdtts      !1/dtts

!chd auxiliary variables to prevent unnecessary computation
      real :: rup,rdown,gnut
      real :: xlight,glight,food 
      real :: gzoo,f12,f14,f23,f24,f34,f43
      real :: epsi      = 1e-12 

real    :: total_trc
!
! =====================================================================
!     begin executable code
! =====================================================================

!       get the model month

call get_date(time%model_time, year, month, day,                &
              hour, minute, second)

area_k(:,:)=0.
tmp(:,:)=0.
dummy1=0.
 
do n = 1, instances  !{  necessary but unused feature
 
include "bio.f90"

poc_tot(:,:)=0
pic_tot(:,:)=0
pic(:,:,:)=0
poc(:,:,:)=0

ratio_poc(1:5)=[1,106,-16,-170,16]
ratio_pic(1:5)=[0,1,2,0,0]
pp(:,:)=0
fmin_poc(:)=1
fmin_pic(:)=1

do k=n_eudepth+1,nk
 fmin_poc(k)= min((grid%zw(k)*1e-2)**(-0.9),1.)
 fmin_pic(k)= min (exp (-(grid%zw(k)-1e2) /3500.),1.) !no remin above 100m
enddo
fmin_poc(nk)=0
fmin_pic(nk)=0


!chd penetration of solar radiation
  rup   = 1;
  do k=1,nk-1
   rdown  =  rsol*exp(-grid%zw(k)/d1sol) +                       &
            (1-rsol)*exp(-grid%zw(k)/d2sol)
   ray(k) = rup - rdown
   rup    = rdown 
  enddo

! Compute New Production and CaCO3 shell production

 do j = jsc, jec  !{ 

!temporary store tracers at tau-1
        do k=1,nk
          do i=isc,iec
            biotr(i,k,1) = t_prog(ind_po4)%field(i,j,k,taum1) 
	     if (biotr(i,k,1).le. 0.02) biotr(i,k,1)=0.
            biotr(i,k,2) = t_prog(ind_o2)%field(i,j,k,taum1)
          enddo
        enddo
  
! compute the PP equation 
        do k=1,n_eudepth
          do i=isc,iec
   pmax_growth(i,k) = abio*bbio**  &
                   (cbio*t_prog(indtemp)%field(i,j,k,tau))
          pp(i,k) = stf_t(i,j)*pmax_growth(i,k)*fice_t(i,j)*   &
                   (1-exp(alpha_bio*ray(k)/pmax_growth(i,k)))

          enddo
        enddo




! Units are mmol/m3/s  m
! Calculate the export production over the whole euphotic zone 
! assumed to be the 4 surface layers (i.e. surface ~50m). 
           do i=isc,iec
           do k=1,min(n_eudepth,grid%kmt(i,j))
              poc(i,j,k) = (biotr(i,k,1) ) /  &
                   ( biotr(i,k,1) + p_k  )  &
                   *Thickness%dht(i,j,k,tau)*pp(i,k)*s_npp
              pic(i,j,k) = rain_ratio*poc(i,j,k)
           enddo  !} k
           enddo
! compute total export of PIC and POC
           do i=isc,iec
           do k=1,min(n_eudepth,grid%kmt(i,j))
              poc_tot(i,j)=poc_tot(i,j)+poc(i,j,k) 
              pic_tot(i,j)=pic_tot(i,j)+pic(i,j,k)
           enddo
           enddo  !} k

! Compute remineralization of POC and PIC   

     do i=isc,iec
       if (grid%kmt(i,j) .gt. n_eudepth) then
         do k=n_eudepth+1,grid%kmt(i,j)-1
           poc(i,j,k)=-poc_tot(i,j)*(fmin_poc(k-1)-fmin_poc(k) )
           pic(i,j,k)=-pic_tot(i,j)*(fmin_pic(k-1)-fmin_pic(k) )
         enddo
! redissolved the remainder in the bottom
           poc(i,j,grid%kmt(i,j))=-poc_tot(i,j) &
	      *(fmin_poc(grid%kmt(i,j)-1)-fmin_poc(nk) )
           pic(i,j,grid%kmt(i,j)) = -pic_tot(i,j) &
	      *(fmin_pic(grid%kmt(i,j)-1) - fmin_pic(nk) )
	elseif (grid%kmt(i,j).gt.0) then
! Section to remineralise POC and PIC in the bottom which
! is less than the euphotic zone i.e. conservation
           poc(i,j,grid%kmt(i,j)) = poc(i,j,grid%kmt(i,j)) -poc_tot(i,j)
           pic(i,j,grid%kmt(i,j)) = pic(i,j,grid%kmt(i,j)) -pic_tot(i,j)
	 else
	       poc(i,j,1:nk)=0.
		pic(i,j,1:nk)=0.
       endif


      enddo  !} i



!  -----
!
! Compute sources and sinks 
!  source(mmol/m^3/s)
dummy(:) = 0.
dummy1 = 0.

! rjm compute the bgc thickness weight time tendency
   do k = 1, nk  !{
    do i = isc, iec  !{
     t_prog(ind_dic)%th_tendency(i,j,k) =                    & 
     (-poc(i,j,k)*ratio_poc(id_dic) -pic(i,j,k)*ratio_pic(id_dic) )  

     t_prog(ind_alk)%th_tendency(i,j,k) =                    & 
     (-poc(i,j,k)*ratio_poc(id_alk) -pic(i,j,k)*ratio_pic(id_alk) )   

     t_prog(ind_o2)%th_tendency(i,j,k) =                    & 
     (-poc(i,j,k)*ratio_poc(id_o2) -pic(i,j,k)*ratio_pic(id_o2) )  

     t_prog(ind_po4)%th_tendency(i,j,k) =                    & 
     (-poc(i,j,k)*ratio_poc(id_po4) -pic(i,j,k)*ratio_pic(id_po4) )  

     t_prog(ind_no3)%th_tendency(i,j,k) =                    & 
     (-poc(i,j,k)*ratio_poc(id_no3) -pic(i,j,k)*ratio_pic(id_no3) )  

     enddo  !} i
    enddo  !} k

 enddo  !} j
enddo  !} n

return
end subroutine bio_v1a
