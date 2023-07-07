import single_sample as ss
import mtj_types_original as mtj
import numpy as np

th_init =  np.pi/2
ph_init = np.pi/2
num_iter = 2500
J_she = 1e-3
Jappl = 1e-3
device_o = mtj.SHE_MTJ_rng(th_init,ph_init,0)
device_f = mtj.SHE_MTJ_rng(th_init,ph_init,0)
ki_f  = device_f.Ki
tmr_f = device_f.TMR
rp_f  = device_f.Rp

fsum_t = 0
fsum_p = 0
osum_t = 0
osum_p = 0
for i in range(num_iter):
    #print(device_f.phi)
    #print(device_f.theta)
    e_f,_,theta_f,phi_f = ss.single_sample(Jappl,J_she,device_f.theta,device_f.phi,ki_f,tmr_f,rp_f)
    device_f.theta=theta_f
    device_f.phi=phi_f
    fsum_t += theta_f
    fsum_p += phi_f

for i in range(num_iter):
    #print(device_o.phi)
    #print(device_o.theta)
    _,e_o = device_o.single_sample(Jappl,J_she)
    osum_t += device_o.theta
    osum_p += device_o.phi

print(f"             fortran    ---*---    objects")
print(f"final theta: {fsum_t/num_iter:.4f}          {osum_t/num_iter:.4f}\n\
final phi:   {fsum_p/num_iter:.4f}          {osum_p/num_iter:.4f}")
