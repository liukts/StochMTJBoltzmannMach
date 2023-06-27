from tokenize import Double
import numpy as np
import matplotlib.pyplot as plt
from she_mod_single import she_mod
from mtj_types import SHE_MTJ_rng
import os
from tqdm import tqdm

# folder to save results
date = "09_30_22"
target_dir = ("RBM_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

# inject device variation function
def inject_add_dev_var(G_in,g_std):
    G_noise = np.random.normal(loc=0,scale=g_std,size=G_in.shape)
    G_out  = G_in + G_noise
    return G_out

# adding cycle-to-cycle noise is the same function 
# but separated for ease of comprehension
def inject_add_cyc_noise(G_in,g_std):
    G_noise = np.random.normal(loc=0,scale=g_std,size=G_in.shape)
    G_out  = G_in + G_noise
    return G_out

# boolean clauses: (X||Y||Z) & (X'||Y||Z) & (X'||Y'||Z) & (X||Y'||Z') & (X'||Y||Z')
# soln: 21/010101 28/011100
W = np.array([[-5, -1, -1, 10, -1, -1], 
              [-1, -7, -2, -2, 10, -1],
              [-1, -2, -7, -2, -1, 10], 
              [10, -2, -2, -7, -1, -1],
              [-1, 10, -1, -1, -5, -1], 
              [-1, -1, 10, -1, -1, -5]])

# synapse array scaling
rmin = 1000
rmax = 3000
gmax = 1/rmin
gmin = 1/rmax
g_dev_var = 100
g_cyc_noise = 50
G = W
G = (G/np.min(W))*(gmax-gmin)+gmin 
G = -G
G[W == 10] = 0
G_base = inject_add_dev_var(G,g_dev_var)

mag_dev_var = 0.01

# print(G)
# boolean clauses: (x||y||z) & (x'||y||z)
# soln: 21/010101 28/011100
# W = np.array([[-3, -1, -1, 10,  0,  0], 
#               [-1, -5, -2, -1, 10,  0],
#               [-1, -2, -5, -1,  0, 10], 
#               [10, -1, -1, -3,  0,  0],
#               [ 0, 10,  0,  0, -1,  0], 
#               [ 0,  0, 10,  0,  0, -1]])


V_start = 5e12
V_end = 2.5e11
steps = 21
V_arr = np.linspace(V_start,V_end,steps)
#print(V_arr)
Iter = 10 # number of Simulations to Run
sols = [] # empty array of solutions
allNeurs = [] # empty array to contain every travelled solution
allEnerg = [] # empty array to keep track on energies
scale = 1e13

# initialize neurons (The variables that make up our Boolean Clauses)
thetas = np.array([np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2])
phis = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))

devs = []
for z in range(len(thetas)):
    devs.append(SHE_MTJ_rng(thetas[z],phis[z],mag_dev_var))

neurs = np.array([[0,0,0,0,0,0]])
weighted = (neurs @ G) # stores the weighted neurons to determine activation probability
sysenergy = (neurs @ G @ neurs.T)

def convertToDec(args):
    sum = 0
    for k in range(0, len(args[0])):
        sum += (args[0][k] * (2**(len(args[0])-k-1)))
    return sum

for f in range(0, Iter):
    thetas = np.array([np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2])
    phis = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))
    neurs = np.array([[0,0,0,0,0,0]])
    G = inject_add_cyc_noise(G_base,g_cyc_noise)
    weighted = (neurs @ G) # stores the weighted neurons to determine activation probability
    sysenergy = (neurs @ G @ neurs.T)
    for h in range(0,6):
        out,energy = devs[h].single_sample(0,0)
        thetas[h] = devs[h].theta
        phis[h] = devs[h].phi
        neurs[0][h] = out
    
    SolArray = []
    energytemp = []

    for v in tqdm(V_arr,leave=False,ncols=80):
        for g in range(10): # iterations per temp
            for h in range(0,6): 
                # print(weighted*scale)
                out,energy = devs[h].single_sample(weighted[0,h]*scale,v)
                thetas[h] = devs[h].theta
                phis[h] = devs[h].phi
                neurs[0][h] = out
            SolArray.append(convertToDec(neurs))
            temp = (neurs @ G @ neurs.T)[0][0]
            energytemp.append(temp)
            G = inject_add_cyc_noise(G_base,g_cyc_noise)
            weighted = (neurs @ G)

    #Function to Convert Binary neurons to Decimal
    allEnerg.append(energytemp)
    allNeurs.append(SolArray)
    sum = convertToDec(neurs)
    sols.append(sum) #Save Solution
    print(sysenergy)
    print(f'Iteration {f+1}/{Iter}, {sum}, {bin(sum)}')



#Graphing of Histogram
# np.save('./sam_outputs/hist_' + date + '.npy',allNeurs)
# np.save('./sam_outputs/sols_' + date + '.npy',sols)
# np.save('./sam_outputs/energy_' + date + '.npy',allEnerg)
# col = ['purple']
# plt.xticks(range(0, 63, 3))
# plt.yticks(range(0, Iter, int(Iter/10)))
# plt.hist(sols,bins=64, color= col)
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
# plt.savefig('./sam_outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')