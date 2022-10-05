#### IOI Matlab/Python codes. These codes calculate the two-dataset and multi-dataset IOIs, the outlier indices for a given set of constraints and the Bayesian interpretation in the case when statistical errors also play a role. Inputs are summary statisitics in terms of the parameter means, uncertainties and correlation coefficients which can be obtained by, e.g., running GetDist to MCMC chains. A Jupyter notebook goes through all important steps.
 
 
#### For the jupyter notebook:

Steps:
1. Put all the .margestats and .corr files (obtained from GetDist) of the constraints of interest in one folder.
2. Run the jupyter notebook that calculate IOIs, the outlier indices and the Bayesian interpretation.


#### If the Bayesian interpretation is not needed (e.g., comparing different results between cosmological simulations), one can use **pyioi** or the Matlab script to calculate IOIs.

###### For **pyioi**:
python pyioi.py -r *Path/to/GetDist/outputs* -o *Output/directory* -p *Paramter1* *Parameter2* *...*

Steps: (In addition to the above 1)
e.g.: `python pyioi.py -r ./batch -o ./IOIouts/IOI.txt -p omegabh2 omegach2 theta logA ns`

###### For the Matlab Script:
A matlab script to calculate two-dataset and multi-dataset and multi-dataset IOIs, as well as all "outlier indices" for an arbitary number of contraints and in an arbitary parameter space. Outputs are stored in e.g., IOI_CMB.txt  

Steps: Steps: (In addition to the above 1)
2. Specify the directory that contains all the .margestats and .corr files.
3. Put the parameter names like below  
    e.g., H0 parameterization in LCDM model: Params = {'omegabh2','omegam','H0','sigma8','ns','tau'}   
    e.g., Theta parameterization in LCDM model: Params = {'omegabh2','omegach2','theta','logA','ns','tau'}   



#### 
#### For questions, contact weikanglin@sjtu.edu.cn
