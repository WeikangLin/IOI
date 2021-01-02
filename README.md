#### IOI Matlab/Python codes. These codes are to calculate the two-dataset and multi-dataset IOIs as well as the outlier indices for a given set of constraints. Those constraints are summary statisitics given by the parameter means, uncertainties and correlation coefficients which can be obtained by, e.g., running GetDist to MCMC chains.

#### For the Matlab Script:
A matlab script to calculate two-dataset and multi-dataset and multi-dataset IOIs, as well as all "outlier indices" for an arbitary number of contraints and in an arbitary parameter space. Outputs are stored in e.g., IOI_CMB.txt  

Steps:
1. Put all the .margestats and .corr files of the constraints of interest in one folder.
2. Specify the directory that contains all the .margestats and .corr files.
3. Put the parameter names like below  
    e.g., H0 parameterization in LCDM model: Params = {'omegabh2','omegam','H0','sigma8','ns','tau'}   
    e.g., Theta parameterization in LCDM model: Params = {'omegabh2','omegach2','theta','logA','ns','tau'}   
 
 
#### For the python script, usage of **pyioi**:

replace the 2nd and the 3rd steps above by:
python pyioi.py -r *Path/to/GetDist/outputs* -o *Output/directory* -p *Paramter1* *Parameter2* *...*

e.g.: `python pyioi.py -r ./batch -o ./IOIouts/IOI.txt -p omegabh2 omegach2 theta logA ns`



#### 
#### For questions, contact wlin23@ncsu.edu
