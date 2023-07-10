import matplotlib.pyplot as plt
import random as rnd
import numpy as np

from   mtj_types import SHE_MTJ_rng
import RRAM_types
import funcs

from tqdm import tqdm
from re import S
import os

#Example Graph:
#  0  1
#  |/\|
#  2  3     
#  |/\|
#  4  5

#Solution is 110011/001100 = 51/12
#Graph G = (V, E) with edge set E and vertex set V

#map temperature to J from paper
sigmoid   =  lambda x: 1/(1+np.exp(x))
def main():
    Vertices = np.array([0,0,0,0,0,0])
    Edges = np.array([[10,10,-1,-1,10,10], #Connections of node 0
                      [10,10,-1,-1,10,10], #Connections of node 1
                      [-1,-1,10,10,-1,-1], #Connections of node 2
                      [-1,-1,10,10,-1,-1], #Connections of node 3
                      [10,10,-1,-1,10,10], #Connections of node 4
                      [10,10,-1,-1,10,10]]) #Connections of node 5
    thetas    = np.full(6,np.pi/2)
    phis      = np.ones_like(thetas)*np.random.uniform(0,2*np.pi,size=np.shape(thetas))

    cb_array = RRAM_types.HfO2
    mag_dev_var = 0
    g_dev_var   = 1      # device to device variation
    g_cyc_noise = 0      # cycle to cycle variation 
    gmin = 1.0/cb_array.HRS
    gmax = 1.0/cb_array.LRS
    
    Edges = (( (Edges-np.min(Edges)/(np.max(Edges)-np.min(Edges)))*(gmax-gmin))+gmin )
    Edges_base = funcs.inject_add_dev_var(Edges,g_dev_var)

    #Initialize Temperature & Step Size (Dictates the stochasticity of the model)
    T_init = 50.00
    temp_to_J =  lambda t: ((t-1)/(T_init-1) * (5e11-1e11))+1e11
    step = 0.001
    Iter = 250 #Number of Simulations to Run
    iter_per_temp = 1
    scale=cb_array.scale
    sols = [] #Empty array of solutions
    devs = [ SHE_MTJ_rng(thetas[i], phis[i], mag_dev_var) for i in range(6)]

    for f in tqdm(range(Iter),leave=False,ncols=80):
        #Generate Random Values to Start 
        for h in range(0,len(Vertices)):
            Vertices[h] = rnd.randint(0, 1)
        weighted = np.dot(Vertices, Edges) #Calculate Weights

        Teff = T_init #Reset Temperature    
        #Iterate until the system has cooled
        while(Teff >= 1):
            for g in range(iter_per_temp): #Iterations per temperature
                Vertices = (funcs.sample_neurons(devs,weighted,scale,temp_to_J(Teff)))
                weighted = np.dot(Vertices, Edges)
            Teff -= step
        sols.append(funcs.convertToDec(Vertices))
        Edges = funcs.inject_add_cyc_noise(Edges_base,g_cyc_noise)
    funcs.my_hist("Max Cut",Iter,sols )

if __name__ == "__main__":
    main()
