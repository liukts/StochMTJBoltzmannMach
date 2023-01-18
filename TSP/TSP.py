#Travelling Salesman Problem

import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "11_07_22"
target_dir = ("TSP_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

#Example Problem:
#        A-----20-----B
#        |\          /|
#        | \        / |
#        |  35    30  |
#        |   \    /   |
#        |    \  /    |
#        |     \/     |
#       42     /\     34
#        |    /  \    |
#        |   /    \   |
#        |  30    35  |
#        | /        \ |
#        |/          \|
#        C-----12-----D

#Solution is A->B->C->D->A
                      # A  B  C  D  
City_Dist = np.array([[-10,20,42,35], #A
                      [20,-10,30,34], #B
                      [42,30,-10,12], #C
                      [35,34,12,-10]])#D
City_Order = np.array([[0,0,0,0], #First City
                       [0,0,0,0], #Second City
                       [0,0,0,0], #Third City
                       [0,0,0,0]])#Fourth City
NumCit = 4

#Initialize Temperature & Step Size (Dictates the stochasticity of the model)
T_init = 10.00
step = 0.01
Iter = 10 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0, NumCit):
        for b in range(0, NumCit ):
            City_Order[h][b] = rnd.randint(0, 1)
    weighted = np.dot(City_Order, City_Dist) #Calculate Weights

    Teff = T_init #Reset Temperature    

    #Iterate until the system has cooled
    while(Teff >= 1):
        for g in range(0, 1): #Iterations per temperature
            for h in range(0,NumCit): #Do this to each Neuron
                for b in range(0, NumCit):
                    rand = rnd.uniform(0, 1) #rand num for determining set probability
                    if rand < sigmoid(weighted[h][b]/Teff):
                        City_Order[h][b] = 1
                    else:
                        City_Order[h][b] = 0
            weighted = np.dot(City_Order, City_Dist)
        Teff -= step
    
    #Function to Convert Binary neurons to city order
    print(City_Order)

    