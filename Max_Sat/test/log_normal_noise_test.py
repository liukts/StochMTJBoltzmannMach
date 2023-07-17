import numpy as np
import matplotlib.pyplot as plt
import itertools
import sys
# import functions to test directly
sys.path.append('../')
import funcs as f

# identically as in ../funcs.py
get_Log10R      = lambda G: np.log10(1.0/G)
sample_noise    = lambda logR,std_dev: np.random.normal(logR,std_dev)
to_G            = lambda noisy_logR: 1.0/pow(10,noisy_logR)
abs_arr         = lambda e: abs(e)
lmap = lambda func, *iterable: list(map(func, *iterable))

def main():
    # ======================
    y_select = "G"  #R,Rlog,G
    plot_select = "hist"    #hist, cdf
    # ======================
    G_LRS = 1.0/100000
    std_dev = 0.375
    G_LRS_mat = np.full((6,6),G_LRS)
    sample_num = 10

    y = []
    for i in range(sample_num):
        G_w_noise = f.inject_add_cyc_noise(G_LRS_mat, std_dev)
        R_w_noise = np.reshape([1.0/G for G in G_w_noise.flatten()],(6,6))
        R_log_w_noise = np.reshape([np.log10(R) for R in R_w_noise.flatten()],(6,6))
        if y_select == "R":
            for elem in R_w_noise.flatten():
                y.append(elem)
        elif y_select == "Rlog":
            for elem in R_log_w_noise.flatten():
                y.append(elem)
        elif y_select == "G":
            for elem in G_w_noise.flatten():
                y.append(elem)

    
    count, bins_count = np.histogram(y, bins=64)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)

    if plot_select == "cdf":
        plt.plot((bins_count[1:]), (cdf))
    elif plot_select == "hist":
        plt.hist(y,100)

    plt.savefig('./noise_test.png')

if __name__ == "__main__":
    main()

