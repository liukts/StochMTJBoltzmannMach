import numpy as np

#Boolean Clauses: (X||Y), (X'||Y')
#Weight Matrix (Created from Boolean Clauses)
W = np.array([[-2, -1, 10, 0],  
              [-1, -2, 0, 10], 
              [10, 0, -2, -1],  
              [0, 10, -1, -2]])

W2 = np.array([[-4, -1, -1, 10, -1, -1], [-1, -7, -2, -2, 10, -1], \
              [-1, -2, -6, -2, -1, 10], [10, -2, -2, -8, -1, -1], \
              [-1, 10, -1, -1, -5, -1], [-1, -1, 10, -1, -1, -5]])

#Initialize Neurons (The variables that make up our Boolean Clauses)
neurs = np.array([[0,0,0,0,0,0]])
for f in range(0, len(neurs[0])):
    neurs[0][f] = (neurs[0][f] + 1) %2
    for h in range(f+1, len(neurs[0])):
        neurs[0][h] = (neurs[0][h] + 1) %2
        curr = (neurs @ W2 @ np.transpose(neurs))[0][0]
        print(neurs)
        print(curr)

neurs = np.array([[0,1,1,1,0,0]])
curr = (neurs @ W2 @ np.transpose(neurs))[0][0]
print(neurs)
print(curr)
