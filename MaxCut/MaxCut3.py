from re import S
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "11_01_22"
target_dir = ("MaxCut3_Sim_" + date)

# if folder does not exist, create it
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")
if not os.path.isdir("./outputs/"):
    os.mkdir("./outputs/")

#Example Graph:
#  0  1
#  |/\|
#  2  3     
#  |/\|
#  4  5

#Solution is 110011/001100 = 51/12
#Graph G = (V, E) with edge set E and vertex set V
Vertices = np.array([0,0,0,0,0,0])
Edges = np.array([[10,10,-1,-1,10,10], #Connections of node 0
                  [10,10,-1,-1,10,10], #Connections of node 1
                  [-1,-1,10,10,-1,-1], #Connections of node 2
                  [-1,-1,10,10,-1,-1], #Connections of node 3
                  [10,10,-1,-1,10,10], #Connections of node 4
                  [10,10,-1,-1,10,10]]) #Connections of node 5

#Initialize Temperature & Step Size (Dictates the stochasticity of the model)
T_init = 10.00
step = 0.01
Iter = 1000 #Number of Simulations to Run
sols = [] #Empty array of solutions

#Initialize Sigmoid (This is what we use to determine probability of activation)
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

for f in range(0, Iter):
    #Generate Random Values to Start 
    for h in range(0,len(Vertices)):
        Vertices[h] = rnd.randint(0, 1)
    weighted = np.dot(Vertices, Edges) #Calculate Weights

    Teff = T_init #Reset Temperature    

    #Iterate until the system has cooled
    while(Teff >= 1):
        for g in range(0, 1): #Iterations per temperature
            for h in range(0,len(Vertices)): #Do this to each Neuron
                rand = rnd.uniform(0, 1) #rand num for determining set probability
                if rand < sigmoid(weighted[h]/Teff):
                    Vertices[h] = 1
                else:
                    Vertices[h] = 0
            weighted = np.dot(Vertices, Edges)
        Teff -= step
    
    #Function to Convert Binary neurons to Decimal
    sum = 0
    for k in range(0, len(Vertices)):
        sum += (Vertices[k] * (2**(len(Vertices)-k-1)))
    sols.append(sum) #Save Solution

#Graphing of Histogram
plt.xticks(range(0, 63, 3))
plt.yticks(range(0, 1000, 100))
plt.hist(sols, bins=64)
plt.xlabel('Value')
plt.ylabel('Frequencry')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./outputs/' + target_dir + '_' + str(Iter) + '_Histogram.png')
