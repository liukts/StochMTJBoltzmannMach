from ctypes import sizeof
import numpy as np
import matplotlib.pyplot as plt
from torch import gradient

date = "09_30_22"
target_dir = ("RBM_Sim_" + date)
sols = np.load('./sam_outputs/hist_' + date + '.npy')
energy = np.load('./sam_outputs/energy_' + date + '.npy')
Iter = 1
#print(sols)

def Contour(sols):
    x = []
    y = []
    z = []
    map = []
    i = 0
    for array in sols:
        x.extend(np.linspace(0, 64, 64))
        y.extend(np.full(64, i))
        map = np.full(64, 0)
        for h in array:
            map[h] += 1
        z.extend(map)
        i += 1
    print(len(sols))
    print(len(x))
    print(len(y))
    print(len(z))
    c = np.linspace(1, 1000, 1000)
    ax = plt.axes(projection ='3d')
    ax.plot3D(x, y, z, c=c)
    ax.set_title('3D line plot geeks for geeks')
    plt.savefig('./sam_outputs/' + target_dir + '_' + str(Iter) + '_Contour.png')

def lineGraph(solSet):
    y = solSet
    x = np.linspace(1, 121, 121)
    plt.plot(x, y, color='blue', marker = 'o')
    plt.title('System Energy Vs Iteration', fontsize=14)
    plt.xlabel('Iteration', fontsize=14)
    plt.ylabel('System Energy', fontsize=14)
    plt.grid(True)
    plt.savefig('./sam_outputs/' + target_dir + '_' + str(Iter) + '_lineGraph.png')
    
lineGraph(energy[0])




