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
City_Dist = np.array([[ 0,20,42,35], #A
                      [20, 0,30,34], #B
                      [42,30, 0,12], #C
                      [35,34,12, 0]])#D
City_Order = np.array([[0,0,0,0], #First City
                       [0,0,0,0], #Second City
                       [0,0,0,0], #Third City
                       [0,0,0,0]])#Fourth City