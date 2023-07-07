import numpy as np

# constants
uB = 9.274e-24
h_bar = 1.054e-34          
u0 = np.pi*4e-7               
e = 1.6e-19                
kb = 1.38e-23   
gamma = 2*u0*uB/h_bar           # Gyromagnetic ratio in m/As
gamma_b = gamma/u0

def draw_norm(x,var,psig):
    if var:
        return x*np.random.normal(1,psig,1)
    else:
        return x

def draw_const(x,var,csig):
    if var:
        return x + np.random.normal(-csig,csig,1)
    else:
        return x
    
class SHE_MTJ_rng():

    def __init__(self,init_th,ph_init,var):
        if init_th == 0:
            print('Init theta cannot be 0, defaulted to pi/100')
            r_init = np.pi/100
            self.theta = r_init
        else:
            self.theta = init_th
        self.phi = ph_init
        self.t_step = 5e-11

        self.bitstr = []
        self.Ki = draw_norm(0.9056364e-3,var,0.05)               # The anisotropy energy in J/m2
        self.TMR = draw_norm(1.2,var,0.05)                       # TMR ratio at V=0,120%  
        self.Rp = draw_norm(8e3,var,0.05)                        # Magenetoresistance at parallel state, 8000 Ohm
        # Stimulation parameters
        self.v_pulse = 0
        self.J_she = 5e11
        self.t_pulse = 10e-9
        self.t_relax = 15e-9

        # MTJ Parameters- This is experimental values from real STT-SOT p-MTJ%
        self.a = 50e-9                       # Width of the MTJ in m
        self.b = 50e-9                       # Length of the MTJ in m
        self.tf = 1.1e-9                     # Thickness of the freelayer in m                           
        self.tox = 1.5e-9                    # Thickness of the MgO barrier in m
        self.P = 0.6                         # Spin polarization of the tunnel current
        self.alpha = 0.03                    # Gilbert damping damping factor
        self.T = 300                         # Temperature in K
        self.Ms = 1.2e6                      # Saturation magnetization in A/m
        self.Bsat = self.Ms*u0
        self.ksi = 75e-15                    # VCMA coefficient in J/vm
        self.gammap = gamma/(1+self.alpha*self.alpha)  # Reduced gyromagnetic ratio
        self.v = self.tf*np.pi*self.b*self.a/4              # Volume of free layer of the elliptical MTJ
        self.Vh = 0.5
        self.delta = 40                      # Thermal stability factor at V=0
        self.A = self.a*self.b*np.pi/4                 # MTJ area
        self.eta = 0.3                       # Spin hall angle
        self.w = 100e-9                    
        self.l = 100e-9
        self.d = 3e-9                        # Width,length and thichness of beta-W strip (heavy metal layer)
        self.A2 = self.d*self.w                        # Cross-sectional area of heavy metal layer 
        self.rho = 200e-8                    # Resistivity of beta-W
        self.R2 = self.rho*self.l/(self.w*self.d)                # Resistance of beta-W
        self.eps_mgo = 4.0                   # relative permittivity of MgO
        self.cap_mgo = 8.854e-12*self.eps_mgo*self.A/self.tox

        self.Hx = 0       
        self.Hy = 0                  
        self.Hz = 0 
        self.Htherm = np.sqrt((2*u0*self.alpha*kb*self.T)/(self.Bsat*gamma_b*self.t_step*self.v))/u0  

        # STT factor;
        self.F = (gamma*h_bar)/(2*u0*e*self.tf*self.Ms)     
    
        # Demagnetization field; This is experimental value
        self.Nx = 0.010613177892974
        self.Ny = 0.010613177892974
        self.Nz = 0.978773644214052
    
    def single_sample(self,Jappl,Jshe_in):
        #takes in two floats
        phi = []
        theta = []
        energy = []
        power = []
        R = []
        energy.append(0)
        #FIXME: idk why power is having this array appended
        #power.append(np.array([0]))
        #NOTE: changing the above to jsut append 0
        power.append(0)
        theta.append(self.theta)
        phi.append(self.phi)
        #print(f"theta init: {self.theta}")#\nphi init: {self.phi}")
        R.append(self.Rp)                # TMR from Parallel state (Start from Rp=60)
        for i in range(int(self.t_pulse/self.t_step)):
            V = self.v_pulse
            J_SHE = Jshe_in
            J_STT = Jappl
            Hk = (2*self.Ki)/(self.tf*self.Ms*u0)-(2*self.ksi*V)/(u0*self.Ms*self.tox*self.tf)                                                                     # effective anisotropy field with VCMA effect
            Ax = self.Hx-self.Nx*self.Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*self.Htherm                                            # x axis component of the effective magnetic field
            Ay = self.Hy-self.Ny*self.Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*self.Htherm                                            # y axis component of the effective magnetic field
            Az = self.Hz-self.Nz*self.Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*self.Htherm   # z axis component of the effective magneitc field
            dphi = self.gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-self.alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+self.alpha*np.cos(phi[-1]))+
                            Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*self.F*self.eta*(np.sin(phi[-1])-self.alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+self.alpha*self.alpha))-((self.alpha*self.F*self.P*J_STT)/(1+self.alpha*self.alpha))
            dtheta = self.gammap*(Ax*(self.alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(self.alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                            Az*self.alpha*np.sin(theta[-1]))-J_SHE*self.F*self.eta*(np.cos(phi[-1])*np.cos(theta[-1])+(self.alpha*np.sin(phi[-1]))/(1+self.alpha*self.alpha))+((self.F*self.P*J_STT)*np.sin(theta[-1])/(1+self.alpha*self.alpha))   
            R1 = self.Rp*(1+(V/self.Vh)**2+self.TMR)/(1+(V/self.Vh)**2+self.TMR*(1+(np.cos(theta[-1]) ))/2)
            power.append(0.5*self.cap_mgo*V**2+self.R2*(np.abs(J_SHE*self.A2))**2+self.R2*(np.abs(J_SHE*self.A2))**2+R1*(J_STT*self.A)**2)
            phi.append(phi[-1]+self.t_step*dphi)                                     
            theta.append(theta[-1]+self.t_step*dtheta)
            energy.append(energy[-1]+self.t_step*power[-1])
            R.append(R1)     # MTJ Resistance 
        for i in range(int(self.t_relax/self.t_step)):
            V = 0
            J_SHE = 0
            J_STT = Jappl
            Hk = (2*self.Ki)/(self.tf*self.Ms*u0)-(2*self.ksi*V)/(u0*self.Ms*self.tox*self.tf);  # effective anisotropy field with VCMA effect
            Ax = self.Hx-self.Nx*self.Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*self.Htherm                                               # x axis component of the effective magnetic field
            Ay = self.Hy-self.Ny*self.Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*self.Htherm                                               # y axis component of the effective magnetic field
            Az = self.Hz-self.Nz*self.Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*self.Htherm   # z axis component of the effective magneitc field
            dphi = self.gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-self.alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+self.alpha*np.cos(phi[-1]))+
                            Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*self.F*self.eta*(np.sin(phi[-1])-self.alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+self.alpha*self.alpha))-((self.alpha*self.F*self.P*J_STT)/(1+self.alpha*self.alpha))
            dtheta = self.gammap*(Ax*(self.alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(self.alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                            Az*self.alpha*np.sin(theta[-1]))-J_SHE*self.F*self.eta*(np.cos(phi[-1])*np.cos(theta[-1])+(self.alpha*np.sin(phi[-1]))/(1+self.alpha*self.alpha))+((self.F*self.P*J_STT)*np.sin(theta[-1])/(1+self.alpha*self.alpha))      
            R1 = self.Rp*(1+(V/self.Vh)**2+self.TMR)/(1+(V/self.Vh)**2+self.TMR*(1+(np.cos(theta[-1]) ))/2)
            power.append(0.5*self.cap_mgo*V**2+self.R2*(np.abs(J_SHE*self.A2))**2+self.R2*(np.abs(J_SHE*self.A2))**2+R1*(J_STT*self.A)**2)
            phi.append(phi[-1]+self.t_step*dphi)                                     
            theta.append(theta[-1]+self.t_step*dtheta)
            energy.append(energy[-1]+self.t_step*power[-1])
            R.append(R1)     # MTJ Resistance 
        bitstr = (1 if np.cos(theta[-1]) > 0 else 0)
        self.theta = theta[-1]
        self.phi = phi[-1]
        #print(f"theta end: {self.theta}")#\nphi end: {self.phi}")
        #FIXME: not sure why this was previously here, a concatenation of a single array
        #power = np.concatenate(power)
        return bitstr,np.sum(power)*self.t_step
    
    def cont_sample(self):
        pass

class VCMA_MTJ_rng():

    def __init__(self,init_th,ph_init,var):
        if init_th == 0:
            print('Init theta cannot be 0, defaulted to pi/100')
            r_init = np.pi/100
            self.theta = r_init
        else:
            self.theta = init_th
        self.phi = ph_init
        self.t_step = 1e-10

        # Stimulation parameters
        self.v_pulse = 1.5
        self.J_she = 0
        self.t_pulse = 50e-9
        self.t_relax = 15e-9
        self.bitstr = []

        # MTJ Parameters- This is experimental values from real STT-SOT p-MTJ%
        self.a = 50e-9                       # Width of the MTJ in m
        self.b = 50e-9                       # Length of the MTJ in m
        self.tf = 1.1e-9                     # Thickness of the freelayer in m                           
        self.tox = 1.5e-9                    # Thickness of the MgO barrier in m
        self.P = 0.6                         # Spin polarization of the tunnel current
        self.alpha = 0.03                    # Gilbert damping damping factor
        self.T = 300                         # Temperature in K
        self.Ki = draw_norm(0.9056364e-3,var,0.05)               # The anisotropy energy in J/m2
        self.Ms = 1.2e6                      # Saturation magnetization in A/m
        self.Bsat = self.Ms*u0
        self.ksi = 75e-15                    # VCMA coefficient in J/vm
        self.gammap = gamma/(1+self.alpha*self.alpha)  # Reduced gyromagnetic ratio
        self.v = self.tf*np.pi*self.b*self.a/4              # Volume of free layer of the elliptical MTJ
        self.Vh = 0.5
        self.delta = 40                      # Thermal stability factor at V=0
        self.TMR = draw_norm(1.2,var,0.05)                       # TMR ratio at V=0,120%  
        self.Rp = draw_norm(8e3,var,0.05)                        # Magenetoresistance at parallel state, 8000 Ohm
        self.A = self.a*self.b*np.pi/4                 # MTJ area
        self.eta = 0.3                       # Spin hall angle
        self.w = 100e-9                    
        self.l = 100e-9
        self.d = 3e-9                        # Width,length and thichness of beta-W strip (heavy metal layer)
        self.A2 = self.d*self.w                        # Cross-sectional area of heavy metal layer 
        self.rho = 200e-8                    # Resistivity of beta-W
        self.R2 = self.rho*self.l/(self.w*self.d)                # Resistance of beta-W
        self.eps_mgo = 4.0                   # relative permittivity of MgO
        self.cap_mgo = 8.854e-12*self.eps_mgo*self.A/self.tox

        self.Hx = 0       
        self.Hy = 0                  
        self.Hz = 0 
        self.Htherm = np.sqrt((2*u0*self.alpha*kb*self.T)/(self.Bsat*gamma_b*self.t_step*self.v))/u0  

        # STT factor;
        self.F = (gamma*h_bar)/(2*u0*e*self.tf*self.Ms)     
    
        # Demagnetization field; This is experimental value
        self.Nx = 0.010613177892974
        self.Ny = 0.010613177892974
        self.Nz = 0.978773644214052
    
    def single_sample(self,Jappl):
        phi = []
        theta = []
        energy = []
        power = []
        R = []
        energy.append(0)
        power.append(0)
        theta.append(self.theta)
        phi.append(self.phi)
        R.append(self.Rp)                # TMR from Parallel state (Start from Rp=60)
        for i in range(int(self.t_pulse/self.t_step)):
            V = self.v_pulse
            J_SHE = self.J_she
            J_STT = Jappl
            Hk = (2*self.Ki)/(self.tf*self.Ms*u0)-(2*self.ksi*V)/(u0*self.Ms*self.tox*self.tf)                                                                     # effective anisotropy field with VCMA effect
            Ax = self.Hx-self.Nx*self.Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*self.Htherm                                            # x axis component of the effective magnetic field
            Ay = self.Hy-self.Ny*self.Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*self.Htherm                                            # y axis component of the effective magnetic field
            Az = self.Hz-self.Nz*self.Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*self.Htherm   # z axis component of the effective magneitc field
            dphi = self.gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-self.alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+self.alpha*np.cos(phi[-1]))+
                            Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*self.F*self.eta*(np.sin(phi[-1])-self.alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+self.alpha*self.alpha))-((self.alpha*self.F*self.P*J_STT)/(1+self.alpha*self.alpha))
            dtheta = self.gammap*(Ax*(self.alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(self.alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                            Az*self.alpha*np.sin(theta[-1]))-J_SHE*self.F*self.eta*(np.cos(phi[-1])*np.cos(theta[-1])+(self.alpha*np.sin(phi[-1]))/(1+self.alpha*self.alpha))+((self.F*self.P*J_STT)*np.sin(theta[-1])/(1+self.alpha*self.alpha))   
            R1 = self.Rp*(1+(V/self.Vh)**2+self.TMR)/(1+(V/self.Vh)**2+self.TMR*(1+(np.cos(theta[-1]) ))/2)
            power.append(0.5*self.cap_mgo*V**2+self.R2*(np.abs(J_SHE*self.A2))**2+self.R2*(np.abs(J_SHE*self.A2))**2+R1*(J_STT*self.A)**2)
            phi.append(phi[-1]+self.t_step*dphi)                                     
            theta.append(theta[-1]+self.t_step*dtheta)
            energy.append(energy[-1]+self.t_step*power[-1])
            R.append(R1)     # MTJ Resistance 
        for i in range(int(self.t_relax/self.t_step)):
            V = 0
            J_SHE = 0
            J_STT = Jappl
            Hk = (2*self.Ki)/(self.tf*self.Ms*u0)-(2*self.ksi*V)/(u0*self.Ms*self.tox*self.tf);  # effective anisotropy field with VCMA effect
            Ax = self.Hx-self.Nx*self.Ms*np.sin(theta[-1])*np.cos(phi[-1])+np.random.normal()*self.Htherm                                               # x axis component of the effective magnetic field
            Ay = self.Hy-self.Ny*self.Ms*np.sin(theta[-1])*np.sin(phi[-1])+np.random.normal()*self.Htherm                                               # y axis component of the effective magnetic field
            Az = self.Hz-self.Nz*self.Ms*np.cos(theta[-1])+Hk*np.cos(theta[-1])+np.random.normal()*self.Htherm   # z axis component of the effective magneitc field
            dphi = self.gammap*(Ax*(-np.cos(theta[-1])*np.cos(phi[-1])-self.alpha*np.sin(phi[-1]))+Ay*(-np.cos(theta[-1])*np.sin(phi[-1])+self.alpha*np.cos(phi[-1]))+
                            Az*np.sin(theta[-1]))/(np.sin(theta[-1]))+J_SHE*self.F*self.eta*(np.sin(phi[-1])-self.alpha*np.cos(phi[-1])*np.cos(theta[-1]))/(np.sin(theta[-1])*(1+self.alpha*self.alpha))-((self.alpha*self.F*self.P*J_STT)/(1+self.alpha*self.alpha))
            dtheta = self.gammap*(Ax*(self.alpha*np.cos(theta[-1])*np.cos(phi[-1])-np.sin(phi[-1]))+Ay*(self.alpha*np.cos(theta[-1])*np.sin(phi[-1])+np.cos(phi[-1]))-
                            Az*self.alpha*np.sin(theta[-1]))-J_SHE*self.F*self.eta*(np.cos(phi[-1])*np.cos(theta[-1])+(self.alpha*np.sin(phi[-1]))/(1+self.alpha*self.alpha))+((self.F*self.P*J_STT)*np.sin(theta[-1])/(1+self.alpha*self.alpha))      
            R1 = self.Rp*(1+(V/self.Vh)**2+self.TMR)/(1+(V/self.Vh)**2+self.TMR*(1+(np.cos(theta[-1]) ))/2)
            power.append(0.5*self.cap_mgo*V**2+self.R2*(np.abs(J_SHE*self.A2))**2+self.R2*(np.abs(J_SHE*self.A2))**2+R1*(J_STT*self.A)**2)
            phi.append(phi[-1]+self.t_step*dphi)                                     
            theta.append(theta[-1]+self.t_step*dtheta)
            energy.append(energy[-1]+self.t_step*power[-1])
            R.append(R1)     # MTJ Resistance 
        bitstr = (1 if np.cos(theta[-1]) > 0 else 0)
        self.theta = theta[-1]
        self.phi = phi[-1]
        return bitstr,np.sum(power)*self.t_step
    
    def cont_sample(self):
        pass
