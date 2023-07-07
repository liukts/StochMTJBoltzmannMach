#from she_mod_single import she_mod
from mtj_types_original import SHE_MTJ_rng
import RRAM_types

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
from tokenize import Double
from tqdm import tqdm
import time
import os

# folder to save results
date = datetime.now().strftime("%m-%d-%y_%H:%M:%S")
target_dir = ("RBM_Sim_" + date)
def main():
    total_start = time.time()
    init_time = time.time()
    # if folder does not exist, create it
    if not os.path.isdir("./outputs/"):
        os.mkdir("./outputs/")

    # boolean clauses: (X||Y||Z) & (X'||Y||Z) & (X'||Y'||Z) & (X||Y'||Z') & (X'||Y||Z')
    # soln: 21/010101 28/011100
    # FIXME: optional: better way to translate problem to matrix?
    W = np.array([[-5, -1, -1, 10, -1, -1], 
                  [-1, -7, -2, -2, 10, -1],
                  [-1, -2, -7, -2, -1, 10], 
                  [10, -2, -2, -7, -1, -1],
                  [-1, 10, -1, -1, -5, -1], 
                  [-1, -1, 10, -1, -1, -5]])

    # synapse array scaling
    #FIXME look at resistance ranges for varying devices and the appropraite scale
    cb_array = RRAM_types.test
    gmax = 1/cb_array.LRS
    gmin = 1/cb_array.HRS
    #print(f"{gmax},{gmin}")
    #FIXME: off for now
    g_dev_var = 0      # device to device variation
    g_cyc_noise = 0    # cycle to cycle variation 

    G = -1 * ( (W/np.min(W))*(gmax-gmin)+gmin )
    #print(G)
    #G[W == 10] = gmax  #Make all elements in G where W is 10 correspond to a high energy
    G[W == 10] = 0  #Make all elements in G where W is 10 correspond to a high energy
    G_base = inject_add_dev_var(G,g_dev_var)
    #print(G_base)

    #FIXME: off for now
    mag_dev_var = 0.0   # magnetic device variation

    ################### SA init ##################
    Iter = 1 # number of Simulations to Run
    iter_per_temp = 5
    steps = 5 #NOTE: granularity of temperature, tweak this
    scale = cb_array.scale  # NOTE: scaling factor for re-input into the system,  refer to fig 3 a)
    # current density start and end (J), probably dont need to change
    J_start = 5e12
    J_end = 1e12
    J_arr = np.linspace(J_start,J_end,steps)
    sols     = [] # empty array of solutions
    allNeurs = [] # empty array to contain every travelled solution
    allEnerg = [] # empty array to keep track on energies

    # initialize neurons (The variables that make up our Boolean Clauses)
    thetas = np.full(6,np.pi/2)
    phis = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))
    devs = [ SHE_MTJ_rng(thetas[z], phis[z], mag_dev_var) for z in range(len(thetas))]
    neurs = np.full((6,6),0)
    weighted = (neurs @ G) # stores the weighted neurons to determine activation probability
    sysenergy = (neurs @ G @ neurs.T)
    init_end = time.time() - init_time

    total_SA_start = time.time()
    total_sample_time = 0
    total_conv_time = 0

    ################# execute SA ##################
    for f in range(Iter):
        for h in range(6):
            #f90 call file_name.module_name.subroutine(args)
            _ = time.time()
            out,energy = devs[h].single_sample(0,0)
            total_sample_time += time.time() - _
            thetas[h] = devs[h].theta
            phis[h] = devs[h].phi
            neurs[0][h] = out
        SolArray = []
        energytemp = []
        for J in tqdm(J_arr,leave=False,ncols=80):
            for g in range(iter_per_temp): # iterations per temp  (you can play with)
                for h in range(6): 
                    #print(amp(weighted[0,h]))
                    #weighted arr is result of VMM, scale that then use as input
                    _ = time.time()
                    out,energy = devs[h].single_sample(scale*weighted[0,h],J)
                    #print(energy)
                    total_sample_time += time.time() - _
                    thetas[h] = devs[h].theta
                    phis[h] = devs[h].phi
                    neurs[0][h] = out
                SolArray.append(convertToDec(neurs))
                print(f"sol {convertToDec(neurs)}")
                #return scalar
                temp = (neurs @ G @ neurs.T)[0][0]
                print(f"temp {temp}")
                input()
                energytemp.append(temp)
                G = inject_add_cyc_noise(G_base,g_cyc_noise)
                weighted = (neurs @ G)
        #Function to Convert Binary neurons to Decimal
        allEnerg.append(energytemp)
        allNeurs.append(SolArray)
        sum = convertToDec(neurs)
        sols.append(sum) #Save Solution
        #print(energytemp[-1])
        print(f"Iteration {f+1}/{Iter}, {sum}, {bin(sum)}")


    print("--- init time: %s seconds ---" % (init_end))
    print("--- total SA time: %s seconds ---" % (time.time() - total_SA_start))
    print("--- total sample time: %s seconds ---" % (total_sample_time))
    print("--- total program time: %s seconds ---" % (time.time() - total_start))

#device variation is gaussian for now, potentially change with something more experimental 
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

def convertToDec(args):
    sum = 0
    for k in range(0, len(args[0])):
        sum += (args[0][k] * (2**(len(args[0])-k-1)))
    return sum

def plot():
    print("FIXME")
    #FIXME: add save and plot functions
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
    
if __name__ == "__main__":
    main()
