import matplotlib.pyplot as plt
import random as rnd
import numpy as np

from   mtj_types import SHE_MTJ_rng
import RRAM_types
import funcs

import multiprocessing as mp
from tqdm import tqdm
from re import S
import time
import os

sigmoid   =  lambda x: 1/(1+np.exp(x))
def SA(sol_queue,get_run_data_flag):
    # ========================== Problem definition =========================
    prob = "Max Sat"

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
    else:
        # dummy Edges for call with get_run_data_flag
        Edges = np.zeros((6,6))
    # =======================================================================
    
   
   
    # ========================== device init  ===============================
    mag_dev_sig = 1
    g_dev_sig   = 0      # device to device variation
    g_cyc_sig = 0      # cycle to cycle variation 

    cb_array = RRAM_types.HfO2
    scale    = cb_array.scale.get(prob) #kinda jank, see RRAM_types.py
    gmin = 1.0/cb_array.HRS
    gmax = 1.0/cb_array.LRS
    
    #   map abstract weights to conductances
    Edges = (( (Edges-np.min(Edges)/(np.max(Edges)-np.min(Edges)))*(gmax-gmin))+gmin )
    Edges_base = funcs.inject_add_dev_var(Edges,g_dev_sig)

    thetas    = np.full(6,np.pi/2)
    phis      = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))
    devs = [ SHE_MTJ_rng(thetas[i], phis[i], mag_dev_sig) for i in range(6)]
    # ===============================================================================



    # ================ annealing schedule ====================
    total_iters = 1000   #Number of Simulations to Run
    iter_per_temp = 1
    T_init = 10.00
    T_step = 0.01
    temp_to_J =  lambda t: ((t-1)/(T_init-1) * (5e11-1e11))+1e11
    # ========================================================
    #   simple way to get run data for parallel process while keeping all variables contained within this function
    if get_run_data_flag == 0:
        pass
    else:
        return prob,total_iters,iter_per_temp, T_init, T_step, mag_dev_sig, g_dev_sig, g_cyc_sig



    # ////////////
    # ////////////
    # ////////////
    # ================================== Exexcute SA ==================================
    #   uncomment if analyzing evolution of algorithm across annealing schedule
    #sols_across_temp   = [] 
    #energy_across_temp = []
    # =================================================
    #   Random state to start
    # =================================================
    Vertices = funcs.sample_neurons(devs,0,scale,0)
    weighted = np.dot(Vertices, Edges) 

    Teff = T_init #Set/Reset Temperature    

    while(Teff >= 1): #J, effectively 
        for g in range(iter_per_temp):
            Vertices = (funcs.sample_neurons(devs,weighted,scale,temp_to_J(Teff)))

            #   uncomment if analyzing evolution of algorithm across annealing schedule
            #energy_across_temp.append(Vertices @ Edges @ np.array(Vertices).T)
            #sols_across_temp.apend(funcs.convertToDec(Vertices)) 
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
            Edges = funcs.inject_add_cyc_noise(Edges_base,g_cyc_sig)
        Teff -= T_step
    solution = funcs.convertToDec(Vertices)
    sol_queue.put(solution)

def print_simulation_setup(batch_size) -> None:
    prob,total_iters,iter_per_temp, T_init, T_step, mag_dev_sig, g_dev_sig, g_cyc_sig  = SA([],1)
    single_sample_time = 0.0002182
    print("====================================================")
    print(f"---------- starting parallel {prob} SA sim with:\n\
            parallel batch size of {batch_size}\n\
            {total_iters} total iterations\n\
            {iter_per_temp} iteration(s) per temp\n\
            {T_init} inital temp\n\
            MTJ device distribution std.dev = {mag_dev_sig}\n\
            G device distribution std.dev = {g_dev_sig}\n\
            G cycle distribution std.dev = {g_cyc_sig}\n\
          \r---------- estimated run time: {(total_iters * single_sample_time * iter_per_temp * T_init/T_step)/(60)}m")
    print("====================================================")
    return prob, total_iters

if __name__ == "__main__":
    total_start_time = time.time()
    # ================================================================
    #  parallelized wrapper for SA(), 
    #  run paramters are only accesible from the function body
    # ================================================================
    
    sols      = []
    #NOTE: fine if personal, change if lab puter
    batch_size = os.cpu_count() * 2
    #   function call to get run data for job creation and for user to see
    prob, total_iters = print_simulation_setup(batch_size)

    #   not to the best way to parallelize since
    #   batches are sequential, that is, even if an open
    #   core is available, it wont run till the slowest
    #   process finishes. good enough for now though.
    sims_to_run = total_iters
    pbar = tqdm(total=total_iters)
    while sims_to_run >= 1:        
        if sims_to_run < batch_size:
            batch_size = sims_to_run
        sol_queue = mp.Queue()  # parallel-safe queue
        processes = []
        #   create processes and start them
        for _ in range(batch_size):
            sim = mp.Process(target=SA, args=(sol_queue,0))
            processes.append(sim)
            sim.start()
        #   waits for solution to be available
        for sim in processes:
            sol = sol_queue.get()  #will block
            sols.append(sol)
        #   wait for all processes to wrap-up before continuing
        for sim in processes:
            sim.join()
        pbar.update(batch_size)
        sims_to_run -= batch_size 
    pbar.close()
    print("--- total program time: %s seconds ---" % (time.time() - total_start_time))
    #   plot
    funcs.my_hist(prob,total_iters,sols)
