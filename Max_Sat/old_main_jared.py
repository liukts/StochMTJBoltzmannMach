from   mtj_types import SHE_MTJ_rng
import RRAM_types
import funcs

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
from tokenize import Double
from tqdm import tqdm
from pathlib import Path
import time
import os

def main():
    # ========================== Problem definition =========================
    total_start = time.time()
    init_time = time.time()
    # FIXME: optional: better way to translate problem to matrix?
    # boolean clauses: (X||Y||Z) & (X'||Y||Z) & (X'||Y'||Z) & (X||Y'||Z') & (X'||Y||Z')
    # soln: 21/010101 28/011100
    W = np.array([[-5, -1, -1, 10, -1, -1], 
                  [-1, -7, -2, -2, 10, -1],
                  [-1, -2, -7, -2, -1, 10], 
                  [10, -2, -2, -7, -1, -1],
                  [-1, 10, -1, -1, -5, -1], 
                  [-1, -1, 10, -1, -1, -5]])

    #FIXME: off for now
    g_dev_var   = 0      # device to device variation
    g_cyc_noise = 0      # cycle to cycle variation 
    mag_dev_var = 0      # magnetic device variation

    # synapse array scaling
    cb_array = RRAM_types.HfO2
    scale = cb_array.scale
    gmax = 1/cb_array.LRS
    gmin = 1/cb_array.HRS

    G = -1 * (( (W-np.min(W)/(np.max(W)-np.min(W)))*(gmax-gmin))+gmin )
    # FIXME: not sure what to set as gmax, more negative or less? do all g's need to be pos or neg? 
    G[W == 10] = gmax  #Make all elements in G where W is 10 correspond to a high energy
    G_base = funcs.inject_add_dev_var(G,g_dev_var)
    # ===============================================================================
    # ////////////
    # ================ annealing schedule ====================
    Iter = 50 # number of Simulations to Run
    iter_per_temp = 1
    steps = 1000  #granularity of temperature

    # current density start and end (J), probably dont need to change
    J_start = 5e11 #NOTE: try 5e12->1e12?
    J_end   = 1e11
    J_arr   = np.linspace(J_start,J_end,steps)
    # ========================================================
    # ////////////
    # ================== initialize neurons (The variables that make up our Boolean Clauses)
    thetas    = np.full(6,np.pi/2)
    phis      = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))
    devs      = [ SHE_MTJ_rng(thetas[i], phis[i], mag_dev_var) for i in range(6)]
    neurs     = np.array([0,0,0,0,0,0])
    weighted  = (neurs @ G) # stores the weighted neurons to determine activation probability

    #FIXME: hmm?
    #sysenergy = (neurs @ W @ neurs.T)
    sysenergy = (neurs @ G @ np.array(neurs).T)
    init_end  = time.time() - init_time

    sols     = [] # empty array of solutions
    allNeurs = [] # empty array to contain every travelled solution
    allEnerg = [] # empty array to keep track on energies
    # ========================================================================================
    # ////////////
    # ////////////
    # ////////////
    # ================================== Exexcute SA =========================================
    total_SA_start = time.time()
    total_sample_time = 0
    for f in tqdm(range(Iter),leave=False,ncols=80):
        # =================================================
        #   Make an initial guess at random to start
        # =================================================
        neurs = funcs.sample_neurons(devs,0,scale,0)
        # =================================================
        SolArray   = []
        energytemp = []
        for J in J_arr: #   effective temperature
            for g in range(iter_per_temp): # iterations per temp  (you can play with)
                #============================
                #   weighted arr is the result of VMM --
                #   once scaled, it's used as input for the array of MTJs
                #
                #============================
                neurs = funcs.sample_neurons(devs,weighted,scale,J)
                #============================
                SolArray.append(funcs.convertToDec(neurs))
                #print(f"current sol arr sampled with J: {convertToDec(neurs)}")
                #return scalar
                #NOTE: temp or energy?
                #NOTE: not the same equation from the paper.
                temp = (neurs @ G @ np.array(neurs).T)
                energytemp.append(temp)
                #print(f"energytemp: {temp}")
                #input()

                #=================================================
                #   add cycle noise to new state before next loop
                #================================================
                G = funcs.inject_add_cyc_noise(G_base,g_cyc_noise)
                weighted = (neurs @ G)
        #///////////////////////////////////////////
        #///////////////////////////////////////////
        allEnerg.append(energytemp)
        allNeurs.append(SolArray)
        sum = funcs.convertToDec(neurs)
        sols.append(sum) #Save Solution
        
        #print(f"Iteration {f+1}/{Iter}, {sum}, {bin(sum)}")
    # ===================================================================================

    print("--- init time: %s seconds ---" % (init_end))
    print("--- total SA time: %s seconds ---" % (time.time() - total_SA_start))
    print("--- total sample time: %s seconds ---" % (total_sample_time))
    print("--- total program time: %s seconds ---" % (time.time() - total_start))
    funcs.my_hist("Max Sat",Iter,sols)


if __name__ == "__main__":
    main()
