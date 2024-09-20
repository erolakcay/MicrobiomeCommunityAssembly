Here's an updated README for the script, incorporating the details you provided:

---

# MicrobiomeCommunityAssembly

## Overview
This Python script simulates microbial communities based on the framework from the research paper:

**"Host control and species interactions jointly determine microbiome community structure"**  
by Eeman Abbasi and Erol Akçay (2022).  
It extends the code implemented by Qian and Akçay (2020).

The simulation explores the effects of host control and species interactions on microbiome community assembly, considering different ecological interaction types (Pm, Pe, Pc) and varying parameters.

## Requirements

Before running the script, ensure you have the following Python libraries installed:
- `numpy`
- `scipy`
- `matplotlib`
- `statsmodels`

You can install these dependencies using `pip`:
```bash
pip install numpy scipy matplotlib statsmodels
```

Additionally, ensure that you have set the following environment variables for efficient processing:
```bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

## Usage

The script `microb_comb.py` can be run from the command line and requires the following seven arguments:

1. **kappa** (int):  
   - `0` for no host control  
   - `> 0` for host control

2. **immune_microbial_load** (int):  
   - Default value: `10000`

3. **immune_response** (boolean):  
   - `True` or `False` to toggle immune response

4. **change_kappa** (boolean):  
   - `True` or `False` to indicate whether kappa should change during the simulation

5. **num_repeats** (int):  
   - Number of simulation repetitions (e.g., `5`)

6. **updated_kappa** (int):  
   - `0` for no host control  
   - `> 0` for host control in the updated simulations

7. **path_to_save** (string):  
   - Path to the directory where the output files will be saved

### Example Command

```bash
python3 microb_comb.py 1 10000 True False 5 1 /path/to/save/results
```

### Output
The output of the simulation is saved as `.npz` files in the designated directory. These files contain the data for each simulated microbial community, which can be used to generate figures and analyze community assembly properties.

The script will create two subdirectories in the `path_to_save` directory:

1. `last_500_simm`:  
   Stores an average of community ecological properties for the last 500 simulations.

2. `final_community`:  
   Stores the community ecological properties of the last simulation.

### Generated Files

The output files are saved with a naming convention that reflects the microbial community's specific parameters. These can be loaded using NumPy for further analysis or visualization.


## Figures
The simulation outputs can be used to generate the figures presented in the manuscript. Each microbial community's ecological properties can be visualized using standard plotting libraries like `matplotlib`.

## Contact
For questions, please contact the authors. 
---

