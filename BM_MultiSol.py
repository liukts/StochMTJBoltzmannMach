import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os
import copy

# folder to save results
date = "06_16_22"
target_dir = ("RBM_MultiSol_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

#Boolean Clauses: (X||Y), (X'||Y')
#Weight Matrix (Created from Boolean Clauses)
W = np.array([[-2, -1, 10, 0],  
              [-1, -2, 0, 10], 
              [10, 0, -2, -1],  
              [0, 10, -1, -2]])

#Initialize Temperature & Step Size (Dictates the _____ of the model)
Teff = 50
step = 5
Iter = 10 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

#Initialize Neurons (The variables that make up our Boolean Clauses)
neurs = np.array([[0,0,0,0]])

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0,len(neurs[0])):
        neurs[0][h] = rnd.randint(0, 1)

    Teff = 50 #Reset Temperature
    #Initialize Energy 
    temp = copy.deepcopy(neurs) #Temporary copy of neurons that holds the next possible solution
    curr = (neurs @ W @ np.transpose(neurs))[0][0]
    #Stores the calculated system energy of the current solution based on the weight matrix
    next = curr #Temporary variable that holds the calculated energy of the next possible solution

    #Iterate until the system has cooled
    while(Teff >= 5):
        for g in range(0,len(neurs[0])):
            temp[0][g] = (temp[0][g]+1) % 2 #Flips one of the neurons (0 to 1 or 1 to 0)
            next = (temp @ W @ np.transpose(temp))[0][0] #Calculates the energy of the next solution 
            rand = rnd.uniform(0.00,1.00) #Generate a random number from 0.00 - 1.00
            if next < curr: #If the next solution is better than the previous, save it
                curr = next
                neurs = copy.deepcopy(temp)
            elif rand < sigmoid(-((Teff*0.16)-4)): 
            #If the random number we generated is less than the number given from the sigmoid based on the temperature
            #Then save this solution even though it is worse
                curr = next
                neurs = copy.deepcopy(temp)
            else: #Reset temporary neurons
                temp = copy.deepcopy(neurs)
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