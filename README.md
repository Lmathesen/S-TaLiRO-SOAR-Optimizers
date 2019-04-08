# S-TaLiRO SOAR Optimizers

This repository contains additional stochastic optimizer packages to be used within the S-TaLiRO falsification toolbox for Matlab. 
The four optimizers in this repository are all instances of the Stochastic Optimization with Adaptive Restart (SOAR) optimization framework. 

This file outlines the installation of the repository within S-TaLiRO, and gives example usages of each of the optimizers.

For documentation on the SOAR algorithms implemented please refer to: *SOAR_S-TaLiRo_Algorithms.pdf* 

The S-TaLiRo toolbox is publicly available for download at (https://sites.google.com/a/asu.edu/s-taliro/s-taliro).

## Authors: 
 * **Logan Mathesen** lmathese@asu.edu - *Arizona State Univeristy*
 * **Giulia Pedrielli** giulia.pedrieli@asu.edu - *Arizona State Univeristy* 

## Contents:
* \[folder\]optimization 
  - \[folder\] auxiliary 
    - \[folder\] Global_Kriging 
  - SOAR_Taliro_LocalGPs.m
  - SOAR_Taliro_LocalGPs_parameters.m
  - SOAR_Taliro_2SPSA.m
  - SOAR_Taliro_2SPSA_parameters.m
  - SOAR_Taliro_SPSA.m
  - SOAR_Taliro_SPSA_parameters.m
  - SOAR_Taliro_FiniteDiff.m
  - SOAR_Taliro_FiniteDiff_parameters.m
* \[folder\]SOAR_examples 
  - \[folder\]models 
  - staliro_SOAR_navbench_example.m
  - staliro_SOAR_modulator_example.m
* setup_staliro.m
* SOAR_S-TaLiRo_Algorithms.pdf
* README.md
  
## Installation
_\*Note, the Matlab Global Optimization toolbox is needed. Please add this toolbox if has not been previously installed._

Quick Install: 
* Replace the `optimization` folder in your current S-TaLiRo `trunk` folder with the one from this repository. 
* Replace the `setup_staliro.m` file in your current S-TaLiRo `trunk` folder with the one from this repository.
* Add the `SOAR_examples` folder to your current S-TaLiRo `trunk` folder. 

Detailed Install: 
1. Ensure you have a working version of S-TaLiRo on your machine.
   - If you have never used S-TaLiRO, download and install from this [page](https://sites.google.com/a/asu.edu/s-taliro/s-taliro/download).
   - Execute a demo from the demos folder to verify correct installation; see this [quick guide](https://df1a2e36-a-0c9971f9-s-sites.googlegroups.com/a/asu.edu/s-taliro/s_taliro_v1_5_quick_guide_v1.pdf?attachauth=ANoY7crkl8lbPPFogHSzGu0vA7JenfGNm5_ZoIbnlN7dAcC1zy9ZAqQ_PPXzBy_vsR3Z3FsGqwMNTslvN2W4IQqNPH5JL0DNV-UTw5OKMxlFqjFD5vVFO-HdKfNP0kNHnAXWsx_MUm7T3Y6QgHHCGHauQjItbdOOZuTmemyOdo0mX5UuRI4yvzj2VfXT1PLrhgCozn-NxJV5IB13W37-z-XFZ_bIcB-tPT-F8UmthNZyN9RnlLGXRys%3D&attredirects=0) for more information.
   - *Note, for hybrid distance metric with distances to the location guards the Matlab package [MatlabBGL](https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl) is required.* 
2. Replace the entire `optimization` folder in the working version of S-TaLio with the folder from this repository.
   - The `optimization` folder from this repository includes all previous optimizers \(and their associated parameter files\) and includes an additional 8 files- an execution and parameter file for each of the new SOAR optimizers.
   - The `optimization` folder from this repository includes an updated `auxiliary` folder, that now includes a `sigmoid.m` function file and a `Global_Kriging` folder with all necessary Gaussian process modeling files. 
3. Replace the `setup_staliro.m` file in your working version of S-TaLiRo with the version from this repository (this file is in the `trunk` folder). 
   - This new setup file's execution is identical to the previous version but now additionally adds the `Global_Kriging` folder to the working path, such that the SOAR optimizers have access to the Gaussian process modeling files.    
4. Add the `SOAR_examples` folder to the `trunk` folder of S-TaLiRo. This folder should be at the same level as the `optimization` folder in your S-TaLiRo directory.
   - This folder contains example execution scripts that call S-TaLiRo with each of the SOAR optimizers. 
   - This folder also contains all of the Simulink models, in the `models` folder, that are necessary for the examples to be executed.
   
## Example Executions
Two example executions of the four SOAR optimizers are included in the `SOAR_examples` folder. First execute the `setup_staliro.m` script and then explore the examples below:
  
1. `staliro_SOAR_navbench_example.m` executes over the Navigation benchmarking model, which is a four dimensional problem consisting of a single input signal with three control points over time and an initial condition. 
This is a model with a hybrid distance metric, such that the `MatlabBGl` is required for successful execution. Further note that the SOAR optimizers currently only support scalar distances so the option `opt.map2line = 1;` in line 159 is critical.
In this case, the optimizer is set to run 10 macro replications with 100 calls to the simulator each. Moreover, SOAR parameters are all set to default values, making use of the crowded expected improvement scoring function for local search restart location decisions. See example 2 for changing SOAR algorithm parameters.
2.  `staliro_SOAR_modulator_example.m` executes over a 3rd order modulator benchmarking model, which is a three dimensional problem that attempts to keep the state of the system within a certain range. 
This is a model with a euclidean distance metric to define system trajectory robustness. As an example, all SOAR algorithms parameters are explicitly set by the user, this can be seen in lines 34-39. Note that line 34, switches SOAR to use the standard expected improvement for local search restart locations.


## Please cite the following papers if you use this code: 
```
Mathesen, L., Pedrielli, G. (2019). Stochastic Optimization With Adaptive Restart: A Framework for Integrated Local and Global Learning. Manuscript submitted for publication in the Journal of Global Optimization.
```
```
Mathesen, L., Yaghoubi, S., Pedrielli, G., & Fainekos, G. (2019). Falsification of Cyber-Physical Systems with Robustness Uncertainty Quantification Through Stochastic Optimization with Adaptive Restart. Manuscript submitted for publication in the 2019 IEEE Conference on Automation Science and Engineering. 
```
```
Mathesen, L., Pedrielli, G., & Ng, S. H. (2017, December). Trust region based stochastic optimization with adaptive restart: A family of global optimization algorithms. In Simulation Conference (WSC), 2017 Winter (pp. 2104-2115). IEEE.
```

##### Please contact Logan Mathesen (lmathese@asu.edu) with any question, comments, bugs, or issues found.