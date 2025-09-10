**Physics-Based-Dynamic Fuel Moisture content (FMC)**

This repository hosts the code related to the paper "A physics-based model of thermodynamically varying fuel moisture content for fire behavior prediction" authored by Ritambhara Raj Dubey and Neda Yaghoobian. The purpose of this repository is to provide access to the code implementation described in the paper. The repository is maintained by the authors of this paper. For more details please contact us by email at nyaghoobian@eng.famu.fsu.edu or rrd20z@fsu.edu.

**Usage- Dynamic FMC model**
Download all files from the folder Dynamic FMC model from Github.
Open fortran using Fortran IDE (such as CodeBlock) on your computer.
Run the file Dynamic_FMC_model.f90
The results will be saved in the same directory as the code.
Sample Input Data- Dynamic FMC model
The input parameters for air temperature and relative humidity are sourced from Matthews, S. (2006), specifically from their paper titled "A process-based model of fine fuel moisture," published in the International Journal of Wildland Fire (Vol. 15, No. 2, pp. 155-168). Additionally, wind speed and solar radiation data are sourced from the Himawari 2011-15: PSM v3 weather data file.

**Usage- Coupled FMC-Fire Dynamics Simulator Model**
To utilize the Coupled FMC-Fire Dynamics Simulator (FDS) model, follow these steps:

1. Download the FDS Code: Obtain the FDS code from the official repository at https://github.com/firemodels/fds.

2. Replace the 'Source' Folder: Within the FDS repository, locate and replace the existing 'Source' folder with the 'Source' folder from the Coupled FMC-FDS model repository.

3. Compile the Code: Compile the modified FDS code to ensure compatibility and integration with the Coupled FMC model.

4. Run the Cases: Execute the model by providing the necessary input files, specifying the simulation parameters and initial conditions for the desired cases.

5. Review Results: Upon completion, the simulation results will be saved in the same directory as the code, allowing for analysis and comparison with expected outcomes.

Sample Input Data- Coupled FMC-Fire Dynamics Simulator model
For validation purposes, sample input data for three validation cases can be found in the 'Input' folder. These inputs have been sourced from relevant literature, details of which are provided in the paper: Dubey, R. R., & Yaghoobian, N. (2024). A physics-based model of thermodynamically varying fuel moisture content for fire behavior prediction. Environmental Modelling & Software, 179, 106111.
