from tokenize import Double
import numpy as np
import matplotlib.pyplot as plt
from she_mod_single import vcma_mod_single
import os
from tqdm import tqdm

# folder to save results
date = "09_30_22"
target_dir = ("RBM_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

# boolean clauses: (X||Y||Z) & (X'||Y||Z) & (X'||Y'||Z) & (X||Y'||Z') & (X'||Y||Z')
# soln: 21/010101 28/011100
W = np.array([[-5, -1, -1, 10, -1, -1], 
              [-1, -7, -2, -2, 10, -1],
              [-1, -2, -7, -2, -1, 10], 
              [10, -2, -2, -7, -1, -1],
              [-1, 10, -1, -1, -5, -1], 
              [-1, -1, 10, -1, -1, -5]])

# boolean clauses: (x||y||z) & (x'||y||z)
# soln: 21/010101 28/011100
# W = np.array([[-3, -1, -1, 10,  0,  0], 
#               [-1, -5, -2, -1, 10,  0],
#               [-1, -2, -5, -1,  0, 10], 
#               [10, -1, -1, -3,  0,  0],
#               [ 0, 10,  0,  0, -1,  0], 
#               [ 0,  0, 10,  0,  0, -1]])


V_start = 10
V_end = -2
step = -0.1
V_arr = np.arange(V_start,V_end-0.01,step)
#print(V_arr)
Iter = 1 # number of Simulations to Run
sols = [] # empty array of solutions
allNeurs = [] # empty array to contain every travelled solution
allEnerg = [] # empty array to keep track on energies
scale = 5e6

# initialize neurons (The variables that make up our Boolean Clauses)
thetas = np.array([np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2])
phis = np.array([0,0,0,0,0,0])
neurs = np.array([[0,0,0,0,0,0]])
weighted = (neurs @ W) # stores the weighted neurons to determine activation probability
sysenergy = (neurs @ W @ neurs.T)

def convertToDec(args):
    sum = 0
    for k in range(0, len(args[0])):
        sum += (args[0][k] * (2**(len(args[0])-k-1)))
    return sum

for f in range(0, Iter):
    for h in range(0,6):
        thetas[h],phis[h],out,energy = vcma_mod_single(thetas[h],phis[h],0,0)
        neurs[0][h] = out
    
    SolArray = []
    energytemp = []

    for v in tqdm(V_arr,leave=False,ncols=80):
        for g in range(0,1): # iterations per temp
            for h in range(0,5): 
                thetas[h],phis[h],out,energy = vcma_mod_single(thetas[h],phis[h],weighted[0,h]*scale,v)
                neurs[0][h] = out
            SolArray.append(convertToDec(neurs))
            temp = (neurs @ W @ neurs.T)[0][0]
            energytemp.append(temp)
            weighted = (neurs @ W)
    
    #Function to Convert Binary neurons to Decimal
    allEnerg.append(energytemp)
    allNeurs.append(SolArray)
    sum = convertToDec(neurs)
    sols.append(sum) #Save Solution
    print(f'Iteration {f+1}/{Iter}, {sum}, {bin(sum)}')



#Graphing of Histogram
np.save('./sam_outputs/hist_' + date + '.npy',allNeurs)
np.save('./sam_outputs/sols_' + date + '.npy',sols)
np.save('./sam_outputs/energy_' + date + '.npy',allEnerg)
col = ['purple']
plt.xticks(range(0, 63, 3))
plt.yticks(range(0, Iter, int(Iter/10)))
plt.hist(sols,bins=64, color= col)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./sam_outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')