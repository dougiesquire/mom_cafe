
! rjm modified ao parameter set

!chd-------------- a bunch of biotic parameters------------------------
!chd Phytoplankton Coefficients
      real ::      alphabio = 0.025 /86400. ![s^-1(Wm^-2)^-1]ini. slope of P-I curve
      real ::      parbio   = 0.43          ![] photosynthe. active rad.
      real ::      kwbio    = 1/0.04        ![m] light att. due to water
      real ::      kcbio    = 1/0.03        ![m] light att. due to phyt
      real ::      abio     = 0.6/86400.    ![1/s]     max growth para.
      real ::      bbio     = 1.066         ![]        max growth para.
      real ::      cbio     = 1.0           ![degC^-1] max growth para.
!     real ::      k1bio    = 0.5           ![mmolm^-3] half-sat. const. for N uptake
! increasing the nutrient half saturation, to allow the nutrient to increase without increasing phytoplankton/chlorophyll, mac jun10
      real ::      k1bio    = 1.0           ![mmolm^-3] half-sat. const. for N uptake
      real ::      muepbio  = 0.01 /86400.  ![1/s] spec. mort. rate
      real ::      muepsbio = 0.0 /86400.   ![1/s(mmolm^-3)^-1] quad. mort. rate
!chd Zooplankton Coefficients
      real ::      gam1bio  = 0.85          ![] assmimilation efficiency
      real ::      gbio     = 2.1 /86400.   ![1/s] max. grazing rate
      real ::      epsbio   = 1.1 /86400.   ![1/s(mmolm^-3)^-2] prey capture rate
      real ::      muezbio  = 0.06/86400.   ![1/s(mmolm^-3)^-1] quad. mortality   
      real ::      gam2bio  = 0.01 /86400.  ![1/s] excretion
!chd Detrital Coefficients
!      real ::      muedbio  = 0.05 /86400.  ![1/s] Remineralization rate
! change to remin to be twice the deep water remin rate, mac jun10.
      real ::      muedbio  = 0.02 /86400.  ![1/s] Remineralization rate
      real ::      wdetbio  = 5.0/86400.    ! [m/s] sinking velocity
! rjm CaCO3
      real ::      muecaco3  = 0.0035 /86400. ![1/s] Remineralization rate
      real ::      wcaco3  = 10./86400.     !35 [m/s] sinking velocity
! rjm Fe
      real ::      tscav_fe = 1./(86400*1.) ! Fe scavenging 1/s

!chd-------------------------------------------------------------------
