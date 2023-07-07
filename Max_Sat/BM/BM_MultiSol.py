import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "06_29_22"
target_dir = ("RBM_MultiSol_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

#Boolean Clauses: (X||Y), (X'||Y')
#Weight Matrix (Created from Boolean Clauses)
     # x   y   x'  y'
W = ([-2, -1, 10,  -10], 
     [-1, -2,  -10, 10], 
     [10,  -10, -2, -1], 
     [ -10, 10, -1, -2])

#Initialize Temperature & Step Size (Dictates the _____ of the model)
T_init = 10.0
step = 0.1
Iter = 1000 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

#Initialize Neurons (The variables that make up our Boolean Clauses)
neurs = ([0,0,0,0])
weighted = np.dot(neurs, W) #Stores the weighted neurons to determine activation probability
#np.dot(np.transpose(neurs), np.dot(neurs, W))

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0,len(neurs)):
        neurs[h] = rnd.randint(0, 1)
    weighted = np.dot(neurs, W) #Stores the weighted neurons to determine activation probability

    Teff = T_init #Reset Temperature    

    #Iterate until the system has cooled
    while(Teff >= 0.1):
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
plt.hist(sols)
plt.xlabel('Value')
plt.ylabel('Frequencry')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')