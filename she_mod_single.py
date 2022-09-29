import numpy as np

def she_mod(init,ph_init,Jappl,J_she):
    t_step = 2e-10
    v_pulse = 0
    J_she = J_she
    t_pulse = 20e-9
    t_relax = 20e-9
    bitstr = []

    uB = 9.274e-24
    h_bar = 1.054e-34          
    u0 = np.pi*4e-7               
    e = 1.6e-19                
    kb = 1.38e-23;              
    gamma = 2*u0*uB/h_bar           # Gyromagnetic ratio in m/As
    gamma_b = gamma/u0

    # MTJ Parameters- This is experimental values from real STT-SOT p-MTJ%
    a = 50e-9                       # Width of the MTJ in m
    b = 50e-9                       # Length of the MTJ in m
    tf = 1.1e-9                     # Thickness of the freelayer in m                           
    tox = 1.5e-9                    # Thickness of the MgO barrier in m
    P = 0.6                         # Spin polarization of the tunnel current
    alpha = 0.03                    # Gilbert damping damping factor
    T = 300                         # Temperature in K
    Ki = 1.0056364e-3               # The anisotropy energy in J/m2
    Ms = 1.2e6                      # Saturation magnetization in A/m
    Bsat = Ms*u0
    ksi = 75e-15                    # VCMA coefficient in J/vm
    gammap = gamma/(1+alpha*alpha)  # Reduced gyromagnetic ratio
    v = tf*np.pi*b*a/4              # Volume of free layer of the elliptical MTJ
    Vh = 0.5
    delta = 40                      # Thermal stability factor at V=0
    TMR = 1.2                       # TMR ratio at V=0,120%  
    Rp = 8e3                        # Magenetoresistance at parallel state, 8000 Ohm
    A = a*b*np.pi/4                 # MTJ area
    eta = 0.3                       # Spin hall angle
    w = 100e-9                    
    l = 100e-9
    d = 3e-9                        # Width,length and thichness of beta-W strip (heavy metal layer)
    A2 = d*w                        # Cross-sectional area of heavy metal layer 
    rho = 200e-8                    # Resistivity of beta-W
    R2 = rho*l/(w*d)                # Resistance of beta-W
    rho_ox = 1e15
    d_ox = 50e-9
    Rox = rho_ox*d_ox/(A)
    
    Hx = 0       
    Hy = 0                  
    Hz = 0 
    Htherm = np.sqrt((2*u0*alpha*kb*T)/(Bsat*gamma_b*t_step*v))/u0  
    
    # STT factor;
    F = (gamma*h_bar)/(2*u0*e*tf*Ms)     
    
    # Demagnetization field; This is experimental value
    Nx = 0.010613177892974
    Ny = 0.010613177892974
    Nz = 0.978773644214052

    #------------------------Initialization----------------------------#

    phi = []
    theta = []
    energy = []
    power = []
    R = []         

    if init:    
        if init == 0:
            print('Init theta cannot be 0, defaulted to pi/100')
            r_init = np.pi/100
            theta.append(r_init)  
        else:
            theta.append(init)
            phi.append(ph_init)
    else:
        print('No init theta provided, defaulted to pi/100')
        r_init = np.pi/100
        theta.append(r_init)  

    energy.append(0)
    power.append(0)
    R.append(Rp)                # TMR from Parallel state (Start from Rp=60)
    for i in range(int(t_pulse/t_step)):
        V = v_pulse
        J_SHE = J_she
        J_STT = Jappl
        Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);                                                                     # effective anisotropy field with VCMA effect
        Ax = Hx-Nx*Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*Htherm                                            # x axis component of the effective magnetic field
        Ay = Hy-Ny*Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*Htherm                                            # y axis component of the effective magnetic field
        Az = Hz-Nz*Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*Htherm   # z axis component of the effective magneitc field
        dphi = gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+alpha*np.cos(phi[-1]))+
                        Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*F*eta*(np.sin(phi[-1])-alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+alpha*alpha))-((alpha*F*P*J_STT)/(1+alpha*alpha))
        dtheta = gammap*(Ax*(alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                        Az*alpha*np.sin(theta[-1]))-J_SHE*F*eta*(np.cos(phi[-1])*np.cos(theta[-1])+(alpha*np.sin(phi[-1]))/(1+alpha*alpha))+((F*P*J_STT)*np.sin(theta[-1])/(1+alpha*alpha))   
        R1 = Rp*(1+(V/Vh)**2+TMR)/(1+(V/Vh)**2+TMR*(1+(np.cos(theta[-1]) ))/2)
        power.append(V**2/Rox+R2*(np.abs(J_SHE*A2))**2+R1*(J_STT*A)**2)
        phi.append(phi[-1]+t_step*dphi)                                     
        theta.append(theta[-1]+t_step*dtheta)
        energy.append(energy[-1]+t_step*power[-1])
        R.append(R1)     # MTJ Resistance 
    for i in range(int(t_relax/t_step)):
        V = 0
        J_SHE = 0
        J_STT = Jappl
        Hk = (2*Ki)/(tf*Ms*u0)-(2*ksi*V)/(u0*Ms*tox*tf);  # effective anisotropy field with VCMA effect
        Ax = Hx-Nx*Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*Htherm                                               # x axis component of the effective magnetic field
        Ay = Hy-Ny*Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*Htherm                                               # y axis component of the effective magnetic field
        Az = Hz-Nz*Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*Htherm   # z axis component of the effective magneitc field
        dphi = gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+alpha*np.cos(phi[-1]))+
                        Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*F*eta*(np.sin(phi[-1])-alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+alpha*alpha))-((alpha*F*P*J_STT)/(1+alpha*alpha))
        dtheta = gammap*(Ax*(alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                        Az*alpha*np.sin(theta[-1]))-J_SHE*F*eta*(np.cos(phi[-1])*np.cos(theta[-1])+(alpha*np.sin(phi[-1]))/(1+alpha*alpha))+((F*P*J_STT)*np.sin(theta[-1])/(1+alpha*alpha))      
        R1 = Rp*(1+(V/Vh)**2+TMR)/(1+(V/Vh)**2+TMR*(1+(np.cos(theta[-1]) ))/2)
        power.append(V**2/Rox+R2*(np.abs(J_SHE*A2))**2+R1*(J_STT*A)**2)
        phi.append(phi[-1]+t_step*dphi)                                     
        theta.append(theta[-1]+t_step*dtheta)
        energy.append(energy[-1]+t_step*power[-1])
        R.append(R1)     # MTJ Resistance 
    bitstr.append(1 if np.cos(theta[-1]) > 0 else 0)
    # R_angle = mod(R,2*pi)*180/pi;
    
    mx = np.sin(theta)*np.cos(phi)              #X axis component of magnetization vector m of free layer
    my = np.sin(theta)*np.sin(phi)              #Y axis component of magnetization vector m of free layer
    mz = np.cos(theta)                          #Z axis component of magnetization vector m of free layer
    G = 1/np.array(R)
    t = np.arange(0,len(mz)*t_step,t_step)
    return theta[-1],phi[-1],bitstr[0],np.sum(power)*t_step