import matplotlib.pyplot as plt
import random as rnd
import numpy as np

from   mtj_types import SHE_MTJ_rng
import RRAM_types
import funcs

from tqdm import tqdm
from re import S
import time
import os

sigmoid   =  lambda x: 1/(1+np.exp(x))
def main():
    total_start_time = time.time()
    # ========================== Problem definition =========================
    prob = "Max Sat" #  used to get correct scale, edge matrix, and plotting function

    if prob == "Max Sat":
        Edges = np.array([[-5, -1, -1, 10, -1, -1], 
                          [-1, -7, -2, -2, 10, -1],
                          [-1, -2, -7, -2, -1, 10], 
                          [10, -2, -2, -7, -1, -1],
                          [-1, 10, -1, -1, -5, -1], 
                          [-1, -1, 10, -1, -1, -5]])
    elif prob == "Max Cut":
        #Example Graph:
        #  0  1
        #  |/\|
        #  2  3     
        #  |/\|
        #  4  5
        #Solution is 110011/001100 = 51/12
        Edges = np.array([[10,10,-1,-1,10,10], 
                          [10,10,-1,-1,10,10], 
                          [-1,-1,10,10,-1,-1], 
                          [-1,-1,10,10,-1,-1], 
                          [10,10,-1,-1,10,10], 
                          [10,10,-1,-1,10,10]])
    # =======================================================================
    
   
   
    # ========================== device init  ===============================
    mag_dev_var = 1
    g_dev_var   = 0      # device to device variation
    g_cyc_noise = 0      # cycle to cycle variation 

    cb_array = RRAM_types.HfO2
    scale    = cb_array.scale.get(prob) #kinda jank, see RRAM_types.py
    gmin = 1.0/cb_array.HRS
    gmax = 1.0/cb_array.LRS
    
    #   map abstract weights to conductances
    Edges = (( (Edges-np.min(Edges)/(np.max(Edges)-np.min(Edges)))*(gmax-gmin))+gmin )
    Edges_base = funcs.inject_add_dev_var(Edges,g_dev_var)

    thetas    = np.full(6,np.pi/2)
    phis      = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))
    devs = [ SHE_MTJ_rng(thetas[i], phis[i], mag_dev_var) for i in range(6)]
    # ===============================================================================



    # ================ annealing schedule ====================
    Iter = 1000   #Number of Simulations to Run
    iter_per_temp = 1
    T_init = 10.00
    step = 0.01
    temp_to_J =  lambda t: ((t-1)/(T_init-1) * (5e11-1e11))+1e11
    # ========================================================



    # ////////////
    # ////////////
    # ////////////
    # ================================== Exexcute SA ==================================
    sample_time_sum = 0 #time MTJ sampling across sims
    sols    = [] 
    energys = []
    for f in tqdm(range(Iter),leave=False,ncols=80):
        # =================================================
        #   Random state to start
        # =================================================
        Vertices = funcs.sample_neurons(devs,0,scale,0)
        weighted = np.dot(Vertices, Edges) 

        Teff = T_init #Set/Reset Temperature    

        while(Teff >= 1): #J, effectively 
            for g in range(iter_per_temp):
                sample_start_time = time.time()
                Vertices = (funcs.sample_neurons(devs,weighted,scale,temp_to_J(Teff)))
                sample_time_sum += time.time() - sample_start_time

                energys.append(Vertices @ Edges @ np.array(Vertices).T)
                #============================
                #   weighted arr is the result of VMM --
                #   once scaled, it's used as input for the array of MTJs
                #   NOTE: i believe the weighted acts as a nearest neighbour function
                #
                #============================
                weighted = np.dot(Vertices, Edges)
                #=================================================
                #   add cycle noise to new state before next loop
                #   NOTE: verify that this shouldnt be between each iter
                #================================================
                Edges = funcs.inject_add_cyc_noise(Edges_base,g_cyc_noise)
            Teff -= step
        sols.append(funcs.convertToDec(Vertices))
    print("--- total sample time: %s seconds ---" % (sample_time_sum))
    print("--- total program time: %s seconds ---" % (time.time() - total_start_time))
    funcs.my_hist(prob,Iter,sols)

if __name__ == "__main__":
    main()
