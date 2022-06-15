import numpy as np
import random as rnd
import copy

#Boolean Clauses: (X||Y||Z), (X'||Y||Z), (X'||Y'||Z), (X||Y'||Z'), (X'||Y||Z')
#Weight Matrix 
W = np.array([[-4, -1, -1, 10, -1, -1], [-1, -7, -2, -2, 10, -1], \
     [-1, -2, -6, -2, -1, 10], [10, -2, -2, -8, -1, -1], \
     [-1, 10, -1, -1, -5, -1], [-1, -1, 10, -1, -1, -5]])


#Initialize Gate Bias Voltage
Vg = 20

#Initialize Temperature & Step Size
Teff = 50
step = 5

#Initialize Set Probability
#Prob = 1/(1+np.exp(-(Vte - Vte0)/Teff))

#Initialize Sigmoid
def sigmoid(x):
    sig = 1/(1+np.exp(x))
    return sig

#Initialize Neurons
neurs = np.array([0,0,0,0,0,0])
for h in range(0,6):
    neurs[h] = rnd.randint(0, 1)
print("Starting Solution:", neurs)

#Initialize Energy
temp = copy.deepcopy(neurs)
curr = np.sum(np.array(np.transpose(neurs) * W * neurs))
next = curr

#Iterate until Temperature is 0
while(Teff >= 0):
    for h in range(0, 6):
        temp[h] = (temp[h]+1) % 2
        next = np.sum(np.transpose(temp) * W * temp)
        rand = rnd.randint(0,10) / 10
        if next < curr:
            curr = next
            neurs = copy.deepcopy(temp)
        elif rand < sigmoid(-((Teff*0.16)-4)):
            curr = next
            neurs = copy.deepcopy(temp)
        else:
            temp = copy.deepcopy(neurs)
        Teff -= step
print("Ending Solution:  ", neurs)