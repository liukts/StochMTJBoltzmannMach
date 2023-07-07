import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "07_20_22"
target_dir = ("RBM_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

#Boolean Clauses: (X||Y||Z), (X'||Y||Z), (X'||Y'||Z), (X||Y'||Z'), (X'||Y||Z')
#Weight Matrix (Created from Boolean Clauses)
     # x   y   z  x'  y'  z'
W = ([-5, -1, -1, 10, -1, -1], 
     [-1, -7, -2, -2, 10, -1], 
     [-1, -2, -7, -2, -1, 10], 
     [10, -2, -2, -7, -1, -1], 
     [-1, 10, -1, -1, -5, -1], 
     [-1, -1, 10, -1, -1, -5])

#Initialize Temperature & Step Size (Dictates the stochasticity of the model)
T_init = 10.00
step = 0.01
Iter = 1000 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

#Initialize Neurons (The variables that make up our Boolean Clauses)
neurs = ([0,0,0,0,0,0])
weighted = np.dot(neurs, W) #Stores the weighted neurons to determine activation probability

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0,len(neurs)):
        neurs[h] = rnd.randint(0, 1)
    weighted = np.dot(neurs, W) #Calculate Weights

    Teff = T_init #Reset Temperature    

    #Iterate until the system has cooled
    while(Teff >= 0.01):
        for g in range(0, 1): #Iterations per temperature
            for h in range(0,len(neurs)): #Do this to each Neuron
                rand = rnd.uniform(0, 1) #rand num for determining set probability
                if rand < sigmoid(weighted[h]/Teff):
                    neurs[h] = 1
                else:
                    neurs[h] = 0
            weighted = np.dot(neurs, W)
        Teff -= step
    
    #Function to Convert Binary neurons to Decimal
    sum = 0
    for k in range(0, len(neurs)):
        sum += (neurs[k] * (2**(len(neurs)-k-1)))
    sols.append(sum) #Save Solution

#Graphing of Histogram
plt.xticks(range(0, 63, 3))
plt.yticks(range(0, 1000, 100))
plt.hist(sols, bins=64)
plt.xlabel('Value')
plt.ylabel('Frequencry')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')