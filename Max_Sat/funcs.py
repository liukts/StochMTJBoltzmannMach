# =================================================================
# July 7,2023
#
# Various helper functions used in the simulated annealing scripts
# ==================================================================
import single_sample as ss

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
from pathlib import Path
import time
import os

def handle_w_path(prob):
    #create dir and write path
    date = datetime.now().strftime("%m-%d-%y_%H:%M:%S")
    out_dir = Path("./outputs")
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    w_file = ("RBM_Sim_" + prob + '_' + date + '_Histogram.png')
    #pathlib POSIX path creation
    w_path = Path( out_dir / w_file)
    return(w_path)

def my_hist(prob,num_iter,sols) -> None:
    # =================================
    # ===== Graphing of Histogram =====
    # =================================
    w_path = handle_w_path(prob)
    bar_col = 'cadetblue' # burnt orange lol
    sol_highlight = "red"
    #=====================================================================
    #   define x ticks and solution-tick highlighting unique for each problem
    if prob == "Max Cut":
        ticks = [0,5,11,12,13,20,25,30,35,40,45,50,51,52,60,68]
        labels = [0,5,' ',12,' ',20,25,30,35,40,45,' ',51,' ',60,68]
        plt.xticks(ticks=ticks,labels=labels)
        plt.gca().get_xticklabels()[3].set_color(sol_highlight)
        plt.gca().get_xticklabels()[-4].set_color(sol_highlight)
    elif prob == "Max Sat":
        ticks = [0,5,10,20,21,22,27,28,29,35,40,45,50,55,60,68]
        labels = [0,5,10,' ',21,' ',' ',28,' ',35,40,45,50,55,60,68]
        plt.xticks(ticks=ticks,labels=labels)
        plt.gca().get_xticklabels()[4].set_color(sol_highlight)
        plt.gca().get_xticklabels()[7].set_color(sol_highlight)
    #======================================================================
    plt.yticks(range(0, num_iter, int(num_iter/10)))
    plt.hist(sols,bins=64,facecolor=bar_col)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(prob + ' Solution Frequency Over ' + str(num_iter) + ' Iterations')
    plt.show()
    print("Save figure? y/n")
    user_input = input()
    if user_input == 'y' or user_input == 'Y':
        plt.savefig(w_path,format='png',dpi=1200)
    elif user_input == 'n' or user_input == 'N':
        pass
    else:
        print("invalid input... saving figure just in case")
        plt.savefig(w_path,format='png',dpi=1200)

def convertToDec(args) -> int:
    #can take optionally a list or a numpy 2D matrix
    sum = 0
    if type(args) == np.ndarray:
        for k in range(0, len(args[0])):
            sum += (args[0][k] * (2**(len(args[0])-k-1)))
    else:
        for k in range(0, len(args)):
            sum += (args[k] * (2**(len(args)-k-1)))
    return sum


#   NOTE: Fortran in here
def sample_neurons(devs,neurons_dot_w,scale,J_step) -> list:
    bits = []
    #   initial setting will be initialized here
    if type(neurons_dot_w) is int:
        neurons_dot_w = np.zeros(6)
    for h in range(6): 
        #   NOTE: python-fortran interface
        #   f90 call looks like import.module_name.function(args)
        #_, out, theta_end, phi_end = ss.single_sample.pulse_then_relax(scale*neurons_dot_w[0][h],J_step,\
        _, out, theta_end, phi_end = ss.single_sample.pulse_then_relax(scale*neurons_dot_w[h],J_step,\
                                                      devs[h].theta,devs[h].phi,                     \
                                                      devs[h].Ki,devs[h].TMR,devs[h].Rp)
        devs[h].theta = theta_end
        devs[h].phi = phi_end
        bits.append( out )
    return bits

#device variation is gaussian for now, potentially change with something more experimental 
# inject device variation function
def inject_add_dev_var(G_in,g_std) -> np.ndarray:
    G_noise = np.random.normal(loc=0,scale=g_std,size=G_in.shape)
    G_out  = G_in + G_noise
    return G_out

# adding cycle-to-cycle noise is the same function 
# but separated for ease of comprehension
def inject_add_cyc_noise(G_in,g_std) -> np.ndarray:
    G_noise = np.random.normal(loc=0,scale=g_std,size=G_in.shape)
    G_out  = G_in + G_noise
    return G_out
