! -------------------------*---------*--------------------------
! These moudles contain the experimental device parameters for 
! a particular device corresponding to each module name.
! All values below are asummed to be fixed and any device-to-device
! or cycle-to-cycle variation is stored within the python dev class.
! -------------------------*---------*--------------------------

module SHE_MTJ_rng_params
    use iso_c_binding, only: sp=>c_float, dp=>c_double
    implicit none
    !these variables do not change per device
    real(dp),parameter :: pi    = 4.0*DATAN(real(1.0,dp))
    real(dp),parameter :: uB    = 9.274e-24
    real(dp),parameter :: h_bar = 1.054e-34          
    real(dp),parameter :: u0    = pi*4e-7               
    real(dp),parameter :: e     = 1.6e-19                
    real(dp),parameter :: kb    = 1.38e-23   
    real(dp),parameter :: gammall = 2.0*u0*uB/h_bar 
    real(dp),parameter :: gammab  = gammall/u0

    real(dp),parameter :: t_step = 5e-11
    real(dp) :: v_pulse = 0.0
    real(dp) :: Jshe = 5e11
    real(dp),parameter :: t_pulse = 10e-9
    real(dp),parameter :: t_relax = 15e-9
    real(dp),parameter :: a  = 50e-9
    real(dp),parameter :: b  = 50e-9
    real(dp),parameter :: tf = 1.1e-9
    real(dp) :: tox = 1.5e-9
    real(dp) :: P   = 0.6
    real(dp),parameter :: alpha = 0.03
    real(dp),parameter :: T     = 300.0
    real(dp),parameter :: Ms    = 1.2e6
    real(dp),parameter :: Bsat  = Ms*u0
    real(dp) :: ksi = 75e-15
    real(dp) :: gammap = gammall/(1+alpha*alpha)
    real(dp),parameter :: volume = tf*pi*b*a/4.0
    real(dp) :: Vh    = 0.5
    real(dp) :: delta = 40.0
    real(dp) :: A1    = a*b*pi/4
    real(dp) :: eta   = 0.3
    real(dp),parameter :: w = 100e-9
    real(dp),parameter :: l = 100e-9
    real(dp),parameter :: d = 3e-9
    real(dp) :: A2 = d*w
    real(dp),parameter :: rho = 200e-8
    real(dp) :: R2 = rho*l/(w*d)
    real(dp) :: eps_mgo = 4.0
    real(dp) :: cap_mgo = 8.854e-12
    real(dp) :: Hx = 0.0
    real(dp) :: Hy = 0.0
    real(dp) :: Hz = 0.0
    real(dp) :: Htherm = sqrt((2*u0*alpha*kb*T)/(Bsat*gammab*t_step*volume))/u0
    real(dp) :: F  = (gammall*h_bar)/(2*u0*e*tf*Ms)
    real(dp) :: Nx = 0.010613177892974
    real(dp) :: Ny = 0.010613177892974
    real(dp) :: Nz = 0.978773644214052
    integer,parameter :: pulse_steps = int(t_pulse/t_step)
    integer,parameter :: relax_steps = int(t_relax/t_step)
end module SHE_MTJ_rng_params

!FIXME: add module for VCMA driven mtj device
