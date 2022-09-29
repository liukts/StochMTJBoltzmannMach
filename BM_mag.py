import numpy as np
import matplotlib.pyplot as plt
from she_mod_single import she_mod
import os
from tqdm import tqdm

# folder to save results
date = "06_20_22"
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


###### MAXIMUM ACCEPTABLE START IS 5E12
###### MINIMUM SHOULD BE 2.3E11
V_start = 5e12
V_end = 2.3e11
steps = 21
V_arr = np.linspace(V_start,V_end,steps)
print(V_arr)
Iter = 10 # number of Simulations to Run
sols = [] # empty array of solutions
scale = 2e9 ##### SCALE CHANGED TO 5E9

# initialize neurons (The variables that make up our Boolean Clauses)
thetas = np.array([np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2,np.pi/2])
phis = np.array([0,0,0,0,0,0])
neurs = np.array([[0,0,0,0,0,0]])
weighted = (neurs @ W) # stores the weighted neurons to determine activation probability
sysenergy = (neurs @ W @ neurs.T)

for f in range(0, Iter):
    for h in range(0,6):
        thetas[h],phis[h],out,energy = she_mod(thetas[h],phis[h],0,0)
        neurs[0][h] = out

    for v in tqdm(V_arr,leave=False,ncols=80):
        for g in range(0,5): # iterations per neuron
            for h in range(0,5): 
                thetas[h],phis[h],out,energy = she_mod(thetas[h],phis[h],weighted[0,h]*scale,v)
                neurs[0][h] = out
        weighted = (neurs @ W)
    
    #Function to Convert Binary neurons to Decimal
    sum = 0
    for k in range(0, len(neurs[0])):
        sum += (neurs[0][k] * (2**(len(neurs[0])-k-1)))
    sols.append(sum) #Save Solution
    print(f'iteration {f+1}/{Iter}, {sum}, {bin(sum)}')

#Graphing of Histogram
np.save('./sam_outputs/hist.npy',sols)
plt.hist(sols,bins=10)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./sam_outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')