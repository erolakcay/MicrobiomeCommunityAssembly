# MicrobiomeCommunityAssembly

The code is for simulation runs from "Host control and species interactions jointly determine microbiome community structure" by 
Eeman Abbasi and Erol Akçay, 2022.

The python script "microb_comb.py" simulates a single microbial community. This python script is imported as an instance of a class to the main python script
"runAllComb.py" which runs simualtions for all possible microbial communities each with a unique Pm, Pe, Pc ecological interaction proportions. 
The "runAllComb.py" script can be run on command line and takes in seven command line arguments:
1. Kappa (int) kappa = 0 no host control, kappa > 1 with host control 
2. Immune microbial load (int)  = 10000
3. Immune response (boolean) = True/False
4. Change kappa (boolean) = True/False
5. Repeats (int) = 5 
6. New kappa (int) = 0 no host control, kappa > 1 with host control 
7. Path to save (string) 

The output for all microbial communities will be saved as a numpy (.npz) file, and named based on each microbial community specific parameters. The data
from the simulations is used to generate the figures presented in the manuscript. 
