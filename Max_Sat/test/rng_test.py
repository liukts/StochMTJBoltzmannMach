import pseudo_rngs as prngs
import time
import math as m
import numpy as np

num_iter = 1000000
mean_zig = 0
mean_bm  = 0

cum_sum_m=0
cum_sum_s=0
prngs.ziggurat.zigset(np.random.randint(1,100000))
for i in range(num_iter):
    x = prngs.ziggurat.rnor()
    cum_sum_m += x
    cum_sum_s += (x-mean_zig)**2

zig_actual_mean = cum_sum_m/num_iter
zig_std_dev = m.sqrt(cum_sum_s/num_iter)

cum_sum_m = 0
cum_sum_s = 0

for i in range(num_iter):
    x=prngs.box_muller.get_box_muller(0.0,1.0)
    cum_sum_m += x
    cum_sum_s += (x-mean_bm)**2

bm_actual_mean = cum_sum_m/num_iter
bm_std_dev = m.sqrt(cum_sum_s/num_iter)



print(f"---------- zig test ---------- ")
print(f"std. dev.: ", zig_std_dev)
print(f"expected mean: ", mean_zig)
print(f"actual mean: ", zig_actual_mean)
print(f"------------------------------ ")
print(f"---------- box_muller test ---------- ")
print(f"std. dev.: ", bm_std_dev)
print(f"expected mean: ", 0)
print(f"actual mean: ", bm_actual_mean)
print(f"------------------------------ ")
