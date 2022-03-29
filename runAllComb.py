
#/Users/eemanabbasi/anaconda/bin/python

import os
from microb_comm import * 
os.environ["OMP_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
os.environ["NUMEXPR_NUM_THREADS"]="1"
import numpy as np
from scipy.stats import halfnorm
from scipy import integrate
#import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
import sys
import warnings

def main():

    #Command line arguments 
    kappa = sys.argv[1]
    immune_microbial_load = sys.argv[2]
    immuneResponse = sys.argv[3]
    change_kappa = sys.argv[4]
    repeats = sys.argv[5]
    new_kappa = sys.argv[6]
    path_to_save = sys.argv[7]

    #repeats = 5
    #path_to_save = '/Users/eemanabbasi/Desktop/Comm_Dynamics_Results' 


    for j in range(0,repeats):
        print("Now running for repeat: ", j)

        # create lists all possible parameter values for Pc, Pm, C, h
        list_Pc = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0]
        list_Pm = [0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0]
    
        # create array of all possible combinations, subject to normalization
        numcombo = len(list_Pc) * len(list_Pm)
        combo = list(itertools.product(list_Pc, list_Pm)) # all possible combinations
        checknorm = []
        for i in range(0,numcombo): # need to check normalization
            checknorm.append(combo[i][0] + combo[i][1] <= 1)
        combo=np.array(combo)
        checknorm=np.array(checknorm)
        combo = combo[checknorm] # mask combo so that normalization is ensured

        fullcombo = list(itertools.product(combo))


        for i in fullcombo:
            params = i 
            Pc = params[0][0]
            Pm = params[0][1]
            print("Now running for comb: ", Pc, Pm)
            run_simulation(Pc,Pm,kappa,immune_microbial_load,immuneResponse,change_kappa,new_kappa,repeats,path_to_save) # 1. rmax, 2. Microbial Load

    return

main()
