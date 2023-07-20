module single_sample
    use iso_c_binding, only: sp=>c_float, dp=>c_double
    use she_mtj_rng_params
    use ziggurat
    implicit none
    contains
        ! --------------------------------*------------*-----------------------------------
        !   energy_usage,bit,theta_end,phi_end, should not be passed into this function
        !   they are declared as they are because subroutine arguments declared with intent(out)
        !   will return back to python as a tuple (not sure exactly why, see numpy f2py)
        !
        !
        !   mimics the python function call: out,energy = dev.single_sample(Jappl,Jshe_in, self.theta, self.phi, self.Ki....)
        !
        ! --------------------------------*------------*-----------------------------------
        subroutine pulse_then_relax(energy_usage,bit,theta_end,phi_end,&
                                    Jappl,Jshe_in,theta_init,phi_init,dev_Ki,dev_TMR,dev_Rp,dump_flag) 
        implicit none
        integer          :: i,t_iter
        real,intent(in)  :: Jappl, Jshe_in,theta_init,phi_init !input
        real,intent(in)  :: dev_Ki,dev_TMR, dev_Rp !device params
        logical, intent(in) :: dump_flag
        !return values
        real,intent(out) :: energy_usage,theta_end,phi_end
        integer,intent(out) :: bit
        !=================================================================
        !   time evolution for solve variables. uncomment if needed. array dump/print from this function is straightforward
        !
        !real,dimension(pulse_steps+relax_steps)    :: R,energy
        real,dimension(pulse_steps+relax_steps)    :: theta_over_time,phi_over_time  
        !==================================================================
        real,dimension(pulse_steps+relax_steps) :: cumulative_pow
        real(dp) :: V,J_SHE,J_STT,Hk,Ax,Ay,Az,dphi,dtheta,R1
        real(dp) :: phi_i,theta_i,power_i,seed!,energy_i

        !                            -*- some notes -*-                                                                       !
        !  static variables do not persist back in python so zigset (rng init function) is called each time this code runs                   !
        !  This code is a module for an existing object-oriented code so the return values from this function                 ! 
        !  should be made to update any existing objects if necessary                                                                    !

        call random_number(seed)
        call zigset(int(1+floor((1000001)*seed)))

        !solve init 
        t_iter  = 1 ! fortran has array indexing of 1, in math terms, t=0
        power_i = 0
        theta_i = theta_init
        phi_i   = phi_init
        if(dump_flag) then
            theta_over_time(t_iter) = theta_i
            phi_over_time(t_iter) = phi_i
        end if
        cumulative_pow(t_iter) = power_i
        !R(1) = dev_Rp
        !energy_i = 0
        !energy(1) = energy_i
        

        !====================== Pulse current and set device to be in-place ======================
        V     = v_pulse
        J_SHE = Jshe_in
        J_STT = Jappl
        Hk    = (2.0*dev_Ki)/(tf*Ms*u0)-(2.0*ksi*V)/(u0*tox*tf)
        do i = 1, pulse_steps
            !keep track of time steps for array navigation
            t_iter=i+1
            Ax = Hx-Nx*Ms*sin(theta_i)*cos(phi_i)     +rnor()*Htherm
            Ay = Hy-Ny*Ms*sin(theta_i)*sin(phi_i)     +rnor()*Htherm
            Az = Hz-Nz*Ms*cos(theta_i)+Hk*cos(theta_i)+rnor()*Htherm

            dphi   = gammap*(Ax*(-cos(theta_i)*cos(phi_i)-alpha*sin(phi_i))+Ay*(-cos(theta_i)*sin(phi_i)+alpha*cos(phi_i))+&
                Az*sin(theta_i))/(sin(theta_i))+J_SHE*F*eta*(sin(phi_i)-alpha*cos(phi_i)*cos(theta_i))/(sin(theta_i)*&
                (1+alpha*alpha))-((alpha*F*P*J_STT)/(1+alpha**2))
            dtheta = gammap*(Ax*(alpha*cos(theta_i)*cos(phi_i)-sin(phi_i))+Ay*(alpha*cos(theta_i)*sin(phi_i)+cos(phi_i))-Az*&
                alpha*sin(theta_i))-J_SHE*F*eta*(cos(phi_i)*cos(theta_i)+(alpha*sin(phi_i))/(1+alpha**2))+((F*P*J_STT)*&
                sin(theta_i)/(1+alpha**2))

            R1     = dev_Rp*(1+(V/Vh)**2+dev_TMR)/(1+(V/Vh)**2+dev_TMR*(1+(cos(theta_i)))/2)
            power_i= 0.5*cap_mgo*V**2+R2*(abs(J_SHE*A2))**2+R2*(abs(J_SHE*A2))**2+R1*(J_STT*A1)**2
            phi_i   = phi_i+t_step*dphi 
            theta_i = theta_i+t_step*dtheta
            cumulative_pow(t_iter) = power_i
            if(dump_flag) then
                theta_over_time(t_iter) = theta_i
                phi_over_time(t_iter) = phi_i
            end if

            !uneeded time evolution
            !R(t_iter)  = R1
            !energy_i= energy_i+t_step*power_i 
            !energy(t_iter) = energy_i
        end do

        !=================  Relax into a one of two low-energy states out-of-plane  ===================
        V=0
        J_SHE = 0
        J_STT = Jappl
        Hk = (2*dev_Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf)
        do i = 1, relax_steps
            t_iter=t_iter+1
            Ax = Hx-Nx*Ms*sin(theta_i)*cos(phi_i)     +rnor()*Htherm
            Ay = Hy-Ny*Ms*sin(theta_i)*sin(phi_i)     +rnor()*Htherm
            Az = Hz-Nz*Ms*cos(theta_i)+Hk*cos(theta_i)+rnor()*Htherm

            dphi = gammap*(Ax*(-cos(theta_i)*cos(phi_i)-alpha*sin(phi_i))+Ay*(-cos(theta_i)*sin(phi_i)+alpha*cos(phi_i))+&
                Az*sin(theta_i))/(sin(theta_i))+J_SHE*F*eta*(sin(phi_i)-alpha*cos(phi_i)*cos(theta_i))/(sin(theta_i)*&
                (1+alpha**2))-((alpha*F*P*J_STT)/(1+alpha**2))
            dtheta = gammap*(Ax*(alpha*cos(theta_i)*cos(phi_i)-sin(phi_i))+Ay*(alpha*cos(theta_i)*sin(phi_i)+cos(phi_i))-&
                Az*alpha*sin(theta_i))-J_SHE*F*eta*(cos(phi_i)*cos(theta_i)+(alpha*sin(phi_i))/(1+alpha**2))+((F*P*J_STT)*&
                sin(theta_i)/(1+alpha**2))

            R1 = dev_Rp*(1+(V/Vh)**2+dev_TMR)/(1+(V/Vh)**2+dev_TMR*(1+(cos(theta_i)))/2)
            power_i = 0.5*cap_mgo*V**2+R2*(abs(J_SHE*A2))**2+R2*(abs(J_SHE*A2))**2+R1*(J_STT*A1)**2
            phi_i   = phi_i+t_step*dphi
            theta_i = theta_i+t_step*dtheta
            cumulative_pow(t_iter) = power_i
            if(dump_flag) then
                theta_over_time(t_iter) = theta_i
                phi_over_time(t_iter) = phi_i
            end if
            !R(t_iter) = R1
            !energy_i= energy_i+t_step*power_i
            !energy(t_iter) = energy_i
        end do

        ! ===== array dump to file of theta/phi time evolution  ====
        if(dump_flag)then
            open(unit = 15, file = "time_evol_phi.txt", action = "write", status = "replace", &
                    form = 'formatted')
            open(unit = 20, file = "time_evol_theta.txt", action = "write", status = "replace", &
                    form = 'formatted')
            write(15,*),phi_over_time
            write(20,*),theta_over_time

            close(15)
            close(20)
        end if
        ! ========================================================== 

        ! ===== return final solve values: energy,bit,theta,phi ====
        theta_end = theta_i
        phi_end = phi_i
        if( cos(theta_end) > 0.0 ) then
            bit = 1
        else
            bit = 0
        end if
        energy_usage = sum(cumulative_pow)*t_step
        end subroutine pulse_then_relax
end module single_sample
