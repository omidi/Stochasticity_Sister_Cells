import numpy as np

p_min = 0
p_max = 120
numb_proteins = (p_max - p_min + 1)
numb_states = numb_proteins*2
p_range =np.arange(p_min, p_max+1)
delta_t = 5

## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50**2
##
