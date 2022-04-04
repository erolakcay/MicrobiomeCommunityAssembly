# MicrobiomeCommunityAssembly

The code is for simulation runs from "Host control and species interactions jointly determine microbiome community structure" by 
Eeman Abbasi and Erol Akçay, 2022.

The python script "microb_comb.py" simulates all possible microbial communities each with a unique combination of Pm, Pe, Pc ecological interaction types. 
The script can be run on command line and takes in seven command line arguments:

1. kappa (int) = 0 no host control, kappa > 0 with host control 
2. immune_microbial_load (int)  = 10000
3. immune_response (boolean) = True/False
4. change_kappa (boolean) = True/False
5. num_repeats (int) = 5 
6. updated_kappa (int) = 0 no host control, kappa > 0 with host control 
7. path_to_save (string) 

Create two folders in the path designated to save named:
a) "last_500_simm": saves an average of community ecological properties for the last 500 simulations. 
b) "final_community":  saves the community ecological properties of the last simulation 

The output for all microbial communities will be saved as a numpy (.npz) file to the designated path , and will be named based on each microbial community specific parameters. The data from the simulations is used to generate the figures presented in the manuscript. 
