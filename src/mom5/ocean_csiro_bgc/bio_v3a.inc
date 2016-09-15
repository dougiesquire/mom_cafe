subroutine bio_v3a (Thickness)

! based on the pt_npzd.ao  (unoptimized version)

type(ocean_thickness_type), intent(in) :: Thickness

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

integer :: i
integer :: j
integer :: k
integer :: n
integer :: ind_no3
integer :: ind_phy
integer :: ind_zoo
integer :: ind_det
integer :: ind_o2
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second
        
!chd-------------- a bunch of biotic parameters------------------------
!chd Phytoplankton Coefficients
      real ::      alphabio = 0.025 /86400. ![s^-1(Wm^-2)^-1]ini. slope of P-I curve
      real ::      parbio   = 0.43          ![] photosynthe. active rad.
      real ::      kwbio    = 1/0.04        ![m] light att. due to water
      real ::      kcbio    = 1/0.03        ![m] light att. due to phyt
      real ::      abio     = 0.6/86400.    ![1/s]     max growth para.
      real ::      bbio     = 1.066         ![]        max growth para.
      real ::      cbio     = 1.0           ![degC^-1] max growth para.
      real ::      k1bio    = 0.5           ![mmolm^-3] half-sat. const. for N uptake
      real ::      muepbio  = 0.03 /86400.  ![1/s] spec. mort. rate
!chd Zooplankton Coefficients
      real ::      gam1bio  = 0.75          ![] assmimilation efficiency
      real ::      gbio     = 2.0 /86400.   ![1/s] max. grazing rate
      real ::      epsbio   = 1.0 /86400.   ![1/s(mmolm^-3)^-2] prey capture rate
      real ::      muezbio  = 0.2 /86400.   ![1/s(mmolm^-3)^-1] quad. mortality   
      real ::      gam2bio  = 0.03 /86400.  ![1/s] excretion
!chd Detrital Coefficients
      real ::      muedbio  = 0.05 /86400.  ![1/s] Remineralization rate
      real ::      wdetbio  = 5.0/86400.    ! [m/s] sinking velocity
!chd-------------------------------------------------------------------

      integer :: ts_npzd ! number of time steps within NPZD model
      integer :: kmeuph = 20 ! deepest level of  euphotic zone
      integer :: tn, trn

      real :: pi    = 3.14159265358979  !yes, this is pi
      real :: biotr(isc:iec,nk,5), bioma(isc:iec,nk,5)
      real :: bion,biop,bioz,biod,bioo
      real :: u_npz,g_npz
      real :: fx1,fx2,fx3,fx4,fu1,fu2,radbio
      real :: vpbio(isc:iec,nk)
      real :: daylen(isc:iec),dtsb,cobeta(isc:iec)
      real :: timeyear,declin
      real :: rayb(isc:iec,0:nk),avej(isc:iec,nk)
!chd auxiliary variables to prevent unnecessary computation
      real :: fbc
      real :: f21,f342,f42,f13,f43,f14
      real :: epsi      = 1e-15
      real :: rdtts      !1/dtts
      real :: adv_fb(isc:iec,1:nk+1)
!-------chd----------------------------------
!
! =====================================================================
!     begin executable code
! =====================================================================

!       get the model month

call get_date(time%model_time, year, month, day,                &
              hour, minute, second)
  timeyear = real(month)/12 + (  real(day) +                    &
                               ( real(hour)+                    &
                                (real(minute)+                  &
                                 real(second)/60)/60)/24)/365

!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------

  ts_npzd= max(1,nint(dtts*2. / 900.))
  rdtts  = 1/dtts

!  write (stdout(),*) ' AO-NPZD model will do ',ts_npzd,' time steps'
!  write (stdout(),*) ' time step in NPZD model will be ',      &
!                dtts/ts_npzd,'sec.'
       
!chd   time step within NPZD model
!chd
      dtsb=2*dtts/float(ts_npzd)

!
!       Loop over multiple instances
!



do n = 1, instances  !{

include "bio.f90" 



!chd calculate biotic source terms using Euler
!chd forward timestepping

 do j = jsc, jec  !{ 


!chd    transfer tracers at tau-1 to temporary arrays
        do k=1,nk
          do i=isc,iec
            biotr(i,k,1) = max(0.0,t_prog(ind_no3)%field(i,j,k,taum1)) 
            biotr(i,k,2) = max(0.0,t_prog(ind_phy)%field(i,j,k,taum1))
            biotr(i,k,3) = max(0.0,t_prog(ind_zoo)%field(i,j,k,taum1))
            biotr(i,k,4) = max(0.0,t_prog(ind_det)%field(i,j,k,taum1))
            biotr(i,k,5) = max(0.0,t_prog(ind_o2)%field(i,j,k,taum1)) 
          enddo
        enddo

!chd   create mask bioma, which is 0 for biotr < epsi and 1
!chd   otherwise
       do trn=1,5   
        do k=1,nk
          do i=isc,iec
           bioma(i,k,trn) = 0.
           if (biotr(i,k,trn) .gt. epsi) then
            bioma(i,k,trn) = 1.
           endif 
          enddo
        enddo
       enddo

!chd calculated length of day

      do i=isc,iec !{
       fx1    = 2.0*pi*timeyear
       declin = 0.006918                                       &
             -0.399912*cos(   fx1)+0.070257*sin(   fx1)        &
             -0.006758*cos(2.*fx1)+0.000907*sin(2.*fx1)        &
             -0.002697*cos(3.*fx1)+0.001480*sin(3.*fx1)
!chd
!chd   compute solar angle at noon (and assume that this is the equivalent
!chd   daily averaged incidence angle for direct+diffuse radiation)
!chd   1. cobeta is cos(incidence angle of solar radiation at noon)
!chd
       fx1    = pi/180. * grid%yt(i,j)
       cobeta(i) = max(0.,sin(fx1)*sin(declin)+cos(fx1)*cos(declin))
!chd
!chd   2. cobeta is cos(incidence angle in water)
!chd      using Snells law: sin(water)=sin(air)/1.33
!chd
       cobeta(i) = sqrt(1.-(1.-cobeta(i)**2.)/1.33**2.)
!chd
!chd   finally,length of day

       fx2    = max(-1., min( 1.,-tan(fx1)*tan(declin)) )
       daylen(i) = max(1.e-12,acos(fx2)/pi)
      enddo  !}i

!chd    compute daily averaged light-limited growth rate
!chd    analytical integration over layer thickness and day
!chd    after Evans and Parslow (1985)
!chd
       avej(:,:)= 0.0
       rayb(:,0)= 1.0

       do k=1,nk
        do i=isc,iec
         rayb(i,k)= exp(-grid%zw(k)/(cobeta(i)*kwbio))
        enddo
       enddo

       do k=1,kmeuph  
        do i=isc,iec
          fx1 = 0.5*daylen(i)*daylen(i)/alphabio  !chd dirty
          radbio = max(1.0,parbio*stf_t(i,j))
          vpbio(i,k) = abio*bbio**                                  &
                  (cbio*t_prog(indtemp)%field(i,j,k,tau))
          fx3 = fx1*vpbio(i,k)/(radbio*rayb(i,k-1))
          fx4 = fx1*vpbio(i,k)/(radbio*rayb(i,k))
          fu1 = sqrt(fx3**2.+daylen(i)**2.)
          fu1 = fu1-daylen(i)*log((daylen(i)+fu1)/fx3)
          fu2 = sqrt(fx4**2.+daylen(i)**2.)
          fu2 = fu2-daylen(i)*log((daylen(i)+fu2)/fx4)
       avej(i,k) = -1.*vpbio(i,k)/(grid%dzt(k)) * cobeta(i)*kwbio   &
                          *(fu1 - fu2 - fx3 + fx4)
        enddo
       enddo


!chd  begin time stepping over dtts in NPZD model with Euler forward
!chd     This is the NPZD model:
!chd      (P: phytoplankton, Z: Zooplankton, N: Nitrate and D: Detritus)
!chd
!chd      dP/dt = u(N,Temp.,Light) P - p_P P - g(P) P Z
!chd
!chd      dZ/dt = a g(P) P Z - d Z - p_Z Z^2
!chd
!chd      dN/dt = r D + d Z - u(N,Temp.,Light) P  [ + r_d DOC ]
!chd
!chd      dD/dt = (1-s)[ (1-a) g(P) P Z + p_P P + p_Z Z^2] -r D + w_D dD/dz

       do tn = 1,ts_npzd  !{
        do k=1,nk  !{
         do i=isc,iec  !{
          bion = max(0.0,biotr(i,k,1))
          biop = max(0.0,biotr(i,k,2))
          bioz = max(0.0,biotr(i,k,3))
          biod = max(0.0,biotr(i,k,4))
          bioo = max(0.0,biotr(i,k,5))

!chd   -- phytoplankton equation
!chdc
!chdc      use Liebigs Law of the Minimum (Liebig, 1845) for growth rate
!chdc      (minimum of light-limited and nutrient limited growth rates;
!chdc      although chlorophyll is not explicitly considered, this will
!chdc      later allow for a diagnostic determination of a Chl:N ratio
!chdc      depending on light- or nutrient-limited growth.
!chdc      --> Hurtt and Armstrong, 1996)
!chdc      saturation growth rate (infinite light, infinite nutrients)
!chd          vpbio = abio*bbio**(cbio*t_prog(indtemp)%field(i,j,k,tau))
!chd      growth rate
          u_npz = min(avej(i,k),vpbio(i,k)*bion/(k1bio+bion))
!chd          if (k .eq. 1) then
!chd           write(stdout(),*) 'ave'
!chd           write(stdout(),*) avej(i,k)
!chd          endif
!chd      grazing function
          g_npz = gbio*epsbio*biop*biop/(gbio+epsbio*biop*biop)
!chd      temperature dependance of growth rates
          fbc = bbio**(cbio*t_prog(indtemp)%field(i,j,k,tau))

!chd auxiliary variables to prevent unnecessary computations
!chd nutrient is compartment 1, phy is 2, zoo is 3, det is 4
!chd f21 := going into phy from nutrient

        f21  = u_npz*biop        * bioma(i,k,1) 
        f342 = g_npz*bioz        * bioma(i,k,2)
        f42  = muepbio*biop      * bioma(i,k,2)
        f13  = gam2bio*bioz      * bioma(i,k,3)
        f43  = muezbio*bioz*bioz * bioma(i,k,3)
        f14  = muedbio*biod      * bioma(i,k,4)
  



!chd   -- nutrient equation
          biotr(i,k,1)  = biotr(i,k,1)  + dtsb * (                     &
               f14 + f13 - f21)  

!chd   -- phyto plankton equation
          biotr(i,k,2)  = biotr(i,k,2)  + dtsb * (                     &
               f21 - f342 - f42) 

!chd   -- zooplankton equation
          biotr(i,k,3)  = biotr(i,k,3)  + dtsb * (                     &
               gam1bio*f342 - f13 - f43)

!chd   -- detritus equation
          biotr(i,k,4)  = biotr(i,k,4)  + dtsb * (                     &
               (1-gam1bio)*f342 + f42 + f43 - f14)    

!chd   -- oxygen equation
          biotr(i,k,5)  = biotr(i,k,5)  - bioma(i,k,5) * 10 * dtsb * ( &
               f14 + f13 - f21) 
 

         enddo  !} i
        enddo  !} k
       enddo  !} tn 


!chd  add biotically induced tendency to biotracers,
        do k = 1, nk  !{
         do i = isc, iec  !{

        t_prog(ind_no3)%th_tendency(i,j,k) = grid%tmask(i,j,k)* & 
            rdtts * 0.5 * Thickness%dht(i,j,k,tau) *            &
         (biotr(i,k,1) - max(0.0,t_prog(ind_no3)%field(i,j,k,taum1)))       
           
        t_prog(ind_phy)%th_tendency(i,j,k) = grid%tmask(i,j,k)* &
            rdtts *  0.5 * Thickness%dht(i,j,k,tau) *           &
           (biotr(i,k,2) - max(0.0,t_prog(ind_phy)%field(i,j,k,taum1)))  
  
        t_prog(ind_zoo)%th_tendency(i,j,k) = grid%tmask(i,j,k)* &
            rdtts *  0.5 * Thickness%dht(i,j,k,tau) *           &
           (biotr(i,k,3) - max(0.0,t_prog(ind_zoo)%field(i,j,k,taum1)))

        t_prog(ind_det)%th_tendency(i,j,k) = grid%tmask(i,j,k)* &
            rdtts *  0.5 * Thickness%dht(i,j,k,tau) *           &
           (biotr(i,k,4) - max(0.0,t_prog(ind_det)%field(i,j,k,taum1)))

        t_prog(ind_o2)%th_tendency(i,j,k) = grid%tmask(i,j,k)*  &
            rdtts *  0.5 * Thickness%dht(i,j,k,tau) *           &
            (biotr(i,k,5) - max(0.0,t_prog(ind_o2)%field(i,j,k,taum1)))

   
         enddo  !} i
        enddo  !} k

     do k=1,nk
      do i=isc,iec
       biotr(i,k,4) = t_prog(ind_det)%field(i,j,k,taum1) 
       biotr(i,k,5) = t_prog(ind_det)%field(i,j,k,taum1) 
      enddo
     enddo

!chd upstream sinking of detritus
     do k=2,nk
      do i=isc,iec
       adv_fb(i,k)=wdetbio*biotr(i,k-1,4)
      enddo
     enddo

!chd no flux boundary conditions
     do i=isc,iec
      adv_fb(i,1)               = 0.0
      adv_fb(i,nk+1)            = 0.0
      adv_fb(i,grid%kmt(i,j)+1) = 0.0
     enddo 


     do k =1,nk  !{
      do i =isc,iec  !{
       biotr(i,k,5) = grid%tmask(i,j,k) &
                      *(-adv_fb(i,k+1) + adv_fb(i,k))
      enddo  !} i
     enddo  !} k


     do k =1,nk  !{
      do i =isc,iec  !{
        t_prog(ind_det)%th_tendency(i,j,k) = t_prog(ind_det)%th_tendency(i,j,k) + &
            grid%tmask(i,j,k) *  biotr(i,k,5) 
      enddo  !} i
     enddo  !} k

 enddo  !} j

enddo  !} n

return
end subroutine bio_v3a
