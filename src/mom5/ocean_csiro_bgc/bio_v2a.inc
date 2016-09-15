subroutine bio_v2a (Thickness)  
! rjm new npzd model

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
character(len=4) :: bgc_stf='stf1'
character(len=4) :: bgc_src='src1'


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

!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     n p z d model
!-----------------------------------------------------------------------

!chd   time step within NPZD model
!chd
! compute time-stepping for npzd model
  ts_npzd=max(1,nint(dtts / 900.))
  rdtts  = 1/dtts
      dtsb=dtts/float(ts_npzd)

!       Loop over multiple instances

do n = 1, instances  !{
 
include "bio.f90"

!chd penetration of solar radiation
  rup   = 1;
  do k=1,nk
   rdown  =  rsol*exp(-grid%zw(k)/d1sol) +                       &
            (1-rsol)*exp(-grid%zw(k)/d2sol)
   ray(k) = rup - rdown
   rup    = rdown 
  enddo

!chd calculate biotic source terms using Euler
!chd forward timestepping, loop order is j,time,k,i 

 do j = jsc, jec  !{ 

!chd    transfer tracers at tau-1 to temporary arrays
        do k=1,nk
          do i=isc,iec
            biotr(i,k,1) = t_prog(ind_no3)%field(i,j,k,taum1) 
            biotr(i,k,2) = t_prog(ind_phy)%field(i,j,k,taum1)
            biotr(i,k,3) = t_prog(ind_zoo)%field(i,j,k,taum1)
            biotr(i,k,4) = t_prog(ind_det)%field(i,j,k,taum1)
            biotr(i,k,5) = t_prog(ind_o2)%field(i,j,k,taum1) 
          enddo
        enddo

       do tn = 1,ts_npzd  !{
        do k=1,nk  !{
         do i=isc,iec  !{
          bion = max(0.0,biotr(i,k,1))
          biop = max(0.0,biotr(i,k,2))
          bioz = max(0.0,biotr(i,k,3))
          biod = max(0.0,biotr(i,k,4))
          bioo = max(0.0,biotr(i,k,5))

         pmax=pmax_0*1       !**(temp-10.)*1e-1
!chd is this really ment to be "pmax_0*1**(temp-10.)*1e-1" ? 
!chd if so, then:  
!chd         pmax=pmax_0 * 1                                       &
!chd              **(t_prog(indtemp)%field(i,j,k,tau)-10)*1e-1

!chd compute auxiliary variables to prevent unnecessary computation 
!chd using the same notation as in  the MOM3 version 

        gnut= bion/(bion+xnut_lim)
        xlight=ray(k)*stf_t(i,j) * fice_t(i,j) ! is already in W/m2
        glight=pmax*xlight/(pmax/alpha + xlight)

        food=biop+biod+  1e-15  ! ensure food is not zero
        gzoo=graze*bioz*food/(food+zhalf)
!chd put in ray(k) 29.3.2000
        f12= (gnut*glight-ray(k)*presp)*biop
        f14 = decay*biod

        f23 = gzoo*biop**2 / (food)
        f24=  pmort*biop

        f34 = zmort*bioz*bioz
        f43 = gzoo*biod**2 / (food)

!chd ensure positive tracer concentrations, whenever
!chd tracerconc. fall below epsi, outgoing biotic fluxes
!chd are stopped
        if (biotr(i,k,1) .lt. epsi) f12 = 0.
        if (biotr(i,k,2) .lt. epsi) f24 = 0.
        if (biotr(i,k,2) .lt. epsi) f23 = 0.
        if (biotr(i,k,3) .lt. epsi) f34 = 0.
        if (biotr(i,k,4) .lt. epsi) f14 = 0.
        if (biotr(i,k,4) .lt. epsi) f43 = 0.
 
!chd   -- nutrient equation
!chd
          biotr(i,k,1) = biotr(i,k,1) + dtsb *                  &
                (-f12 +f14 +f24*pfrac+f34*zfrac)

!chd   -- phytoplankton equation
!chd
          biotr(i,k,2) = biotr(i,k,2) + dtsb *                  &
           (f12 - f24 - f23) 

!chd   -- zooplankton equation
!chd
          biotr(i,k,3) = biotr(i,k,3) + dtsb *                  &
           (zeff*(f23+f43) -f34) 

!chd   -- detritus equation
!chd
          biotr(i,k,4) = biotr(i,k,4) + dtsb *                  &
           ((1-zeff)*(f23+f43) + f24*(1-pfrac)-f14              &
            +f34*(1-zfrac) -f43)

!chd   -- oxygen equation
!chd
          biotr(i,k,5) = biotr(i,k,5) - 10 * dtsb *             &
                (-f12 +f14 +f24*pfrac+f34*zfrac)
 
         enddo  !} i
        enddo  !} k
       enddo  !} tn 


!chd  add biotically induced tendency to biotracers,
        do k = 1, nk  !{
         do i = isc, iec  !{
        t_prog(ind_no3)%th_tendency(i,j,k) =                    & 
            rdtts *  Thickness%dht(i,j,k,tau) *                 &
         (biotr(i,k,1) - t_prog(ind_no3)%field(i,j,k,taum1))       
           
        t_prog(ind_phy)%th_tendency(i,j,k) =                    &
            rdtts *  Thickness%dht(i,j,k,tau) *                 &
           (biotr(i,k,2) - t_prog(ind_phy)%field(i,j,k,taum1))  
  
        t_prog(ind_zoo)%th_tendency(i,j,k) =                    &
            rdtts *  Thickness%dht(i,j,k,tau) *                 &
           (biotr(i,k,3) - t_prog(ind_zoo)%field(i,j,k,taum1))

        t_prog(ind_det)%th_tendency(i,j,k) =                    &
            rdtts *  Thickness%dht(i,j,k,tau) *                 &
           (biotr(i,k,4) - t_prog(ind_det)%field(i,j,k,taum1))
 
        t_prog(ind_o2)%th_tendency(i,j,k) =                     &
            rdtts *  Thickness%dht(i,j,k,tau) *                 &
            (biotr(i,k,5) - t_prog(ind_o2)%field(i,j,k,taum1))

!rjm biotic source terms stored in  t_prog(ind_no3)%th_tendency(i,j,k) 
   
         enddo  !} i
        enddo  !} k

!chd crude first order upstream advection mimicking sinking
!chd of detritus; detritus just falls out of the model
!chd at the bottom mimicking burial

     do k=1,nk
      do i=isc,iec
       biotr(i,k,4)= t_prog(ind_det)%field(i,j,k,taum1)
	biotr(i,k,5) = 0.
      enddo
     enddo 

! rjm only need the thickness weighted time tendency
     do k=2,nk  !{
      do i=isc,iec  !{
       biotr(i,k,5)  =   (-1.)* w_detr* (       &
               biotr(i,k,4) - biotr(i,k-1,4))                        
      enddo  !} i
     enddo  !} k

     do i=isc,iec  !{
       biotr(i,1,5) =  (-1.)*w_detr* biotr(i,2,4)
! rjm - fix for bottom to retain detritus in bottom layer
       biotr(i,grid%kmt(i,j),5 ) = w_detr* biotr(i,grid%kmt(i,j)-1,4)
     enddo  !} i

     do k =1,nk  !{
      do i =isc,iec  !{
      t_prog(ind_det)%th_tendency(i,j,k) =                   & 
            t_prog(ind_det)%th_tendency(i,j,k) +             &
         biotr(i,k,5) 
      enddo  !} i
     enddo  !} k

 enddo  !} j

enddo  !} n

return
end subroutine bio_v2a

