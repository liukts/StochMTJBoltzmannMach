import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "06_20_22"
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
W = np.array([[-4, -1, -1, 10, -1, -1], [-1, -7, -2, -2, 10, -1], \
              [-1, -2, -6, -2, -1, 10], [10, -2, -2, -8, -1, -1], \
              [-1, 10, -1, -1, -5, -1], [-1, -1, 10, -1, -1, -5]])

#Initialize Temperature & Step Size (Dictates the _____ of the model)
Teff = 50
step = 5
Iter = 1000 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

#Initialize Neurons (The variables that make up our Boolean Clauses)
neurs = np.array([[1,0,0,1,1,1]])
weighted = (neurs @ W) #Stores the weighted neurons to determine activation probability

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0,6):
        neurs[0][h] = rnd.randint(0, 1)

    Teff = 50 #Reset Temperature    

    #Iterate until the system has cooled
    while(Teff >= 1):
        for g in range(0, 6): #Iterations per temperature
            for h in range(0,5): #Do this to each Neuron
                rand = rnd.uniform(0, 1) #rand num for determining set probability
                if rand < sigmoid(weighted[0][h]/Teff):
                    neurs[0][h] = 1
                else:
                    neurs[0][h] = 0
                weighted = (neurs @ W)
        Teff -= step
    
    #Function to Convert Binary neurons to Decimal
    sum = 0
    for k in range(0, len(neurs[0])):
        sum += (neurs[0][k] * (2**(len(neurs[0])-k-1)))
    sols.append(sum) #Save Solution

#Graphing of Histogram
plt.hist(sols)
plt.xlabel('Value')
plt.ylabel('Frequencry')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')