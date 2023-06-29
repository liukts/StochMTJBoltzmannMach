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
# optional: better way to translate problem to matrix?
W = np.array([[-5, -1, -1, 10, -1, -1], 
              [-1, -7, -2, -2, 10, -1],
              [-1, -2, -7, -2, -1, 10], 
              [10, -2, -2, -7, -1, -1],
              [-1, 10, -1, -1, -5, -1], 
              [-1, -1, 10, -1, -1, -5]])

# synapse array scaling
# important bit
rmin = 1000
rmax = 3000
gmax = 1/rmin
gmin = 1/rmax
g_dev_var = 0      # device to device variation
g_cyc_noise = 0    # cycle to cycle variation 
G = W
G = (G/np.min(W))*(gmax-gmin)+gmin 
G = -G
G[W == 10] = gmax
G_base = inject_add_dev_var(G,g_dev_var)

mag_dev_var = 0.0   # magnetic device variation

# current density start and end (J)
V_start = 5e12
V_end = 1e12
steps = 21 # you should change
V_arr = np.linspace(V_start,V_end,steps)
#print(V_arr)
Iter = 10 # number of Simulations to Run
sols = [] # empty array of solutions
allNeurs = [] # empty array to contain every travelled solution
allEnerg = [] # empty array to keep track on energies
# scaling factor for re-input into the system
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

    for h in range(0,6):
        out,energy = devs[h].single_sample(0,0)
        thetas[h] = devs[h].theta
        phis[h] = devs[h].phi
        neurs[0][h] = out
    
    SolArray = []
    energytemp = []

    for v in tqdm(V_arr,leave=False,ncols=80):
        for g in range(2): # iterations per temp  (you can play with)
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
    print(energytemp[-1])
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