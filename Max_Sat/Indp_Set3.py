from re import S
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import os

# folder to save results
date = "1_24_23"
target_dir = ("Indp_Sim3_" + date)

# if folder does not exist, create it
if not os.path.isdir("./Indp_Set/"):
    os.mkdir("./Indp_Set/")
if not os.path.isdir("./Indp_Set/"):
    os.mkdir("./Indp_Set/")
if not os.path.isdir("./Indp_Set/"):
    os.mkdir("./Indp_Set/")

#Example Graph:
#           0---3
#              /|
#             1-2---4

#Solution is 11001 = 25
#Graph G = (V, E) with edge set E and vertex set V
                   # 0 1 2 3 4
Vertices = np.array([0,0,0,0,0])
Edges = np.array([[-10,-1,-1,10,-1], #Connections of node 0
                  [-1,-10,10,10,-1], #Connections of node 1
                  [-1,10,-10,10,10], #Connections of node 2
                  [10,10,10,-10,-1], #Connections of node 3
                  [-1,-1,10,-1,-10]]) #Connections of node 4        

#Initialize Temperature & Step Size (Dictates the stochasticity of the model)
T_init = 10.00
step = 0.01
Iter = 100 #Number of Simulations to Run
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
    while(Teff >= 0.01):
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
    print(Vertices)
    sols.append(sum) #Save Solution

#Graphing of Histogram
plt.xticks(range(0, 64, 1))
plt.yticks(range(0, 1000, 100))
plt.hist(sols, bins=64)
plt.xlabel('Value')
plt.ylabel('Frequencry')
plt.title('Solution Frequency Over ' + str(Iter) + ' Iterations')
plt.savefig('./Indp_Set/' + target_dir + '_' + str(Iter) + '_Histogram.png')
