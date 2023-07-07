! -------------------------*---------*--------------------------
! These moudles contain the experimental device parameters for 
! a particular device corresponding to each module name.
! All values below are asummed to be fixed and any device-to-device
! or cycle-to-cycle variation is stored within the python dev classes.
! -------------------------*---------*--------------------------

module SHE_MTJ_rng_params
    implicit none
    !these variables do not change per device
    real,parameter :: pi    = 4.0*ATAN(1.0)
    real,parameter :: uB    = 9.274e-24
    real,parameter :: h_bar = 1.054e-34          
    real,parameter :: u0    = pi*4e-7               
    real,parameter :: e     = 1.6e-19                
    real,parameter :: kb    = 1.38e-23   
    real,parameter :: gammall = 2*u0*uB/h_bar 
    real,parameter :: gammab  = gammall/u0

    real,parameter :: t_step = 5e-11
    real :: v_pulse = 0.0
    real :: Jshe = 5e11
    real,parameter :: t_pulse = 10e-9
    real,parameter :: t_relax = 15e-9
    real,parameter :: a  = 50e-9
    real,parameter :: b  = 50e-9
    real,parameter :: tf = 1.1e-9
    real :: tox = 1.5e-9
    real :: P   = 0.6
    real,parameter :: alpha = 0.03
    real,parameter :: T     = 300.0
    real,parameter :: Ms    = 1.2e6
    real,parameter :: Bsat  = Ms*u0
    real :: ksi = 75e-15
    real :: gammap = gammall/(1+alpha*alpha)
    real,parameter :: volume = tf*pi*b*a/4
    real :: Vh    = 0.5
    real :: delta = 40.0
    real :: A1    = a*b*pi/4
    real :: eta   = 0.3
    real,parameter :: w = 100e-9
    real,parameter :: l = 100e-9
    real,parameter :: d = 3e-9
    real :: A2 = d*w
    real,parameter :: rho = 200e-8
    real :: R2 = rho*l/(w*d)
    real :: eps_mgo = 4.0
    real :: cap_mgo = 8.854e-12
    real :: Hx = 0.0
    real :: Hy = 0.0
    real :: Hz = 0.0
    real :: Htherm = sqrt((2*u0*alpha*kb*T)/(Bsat*gammab*t_step*volume))/u0
    real :: F  = (gammall*h_bar)/(2*u0*e*tf*Ms)
    real :: Nx = 0.010613177892974
    real :: Ny = 0.010613177892974
    real :: Nz = 0.978773644214052
    integer,parameter :: pulse_steps = int(t_pulse/t_step)
    integer,parameter :: relax_steps = int(t_relax/t_step)
end module SHE_MTJ_rng_params

!FIXME: add module for VCMA driven mtj device
