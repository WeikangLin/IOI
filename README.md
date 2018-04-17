# IOI_CosmoMC

A matlab script for CosmoMC (after running getdist) to calculate two-experiment IOIs for arbitary number of contraints and in an arbitary parameter space outputs in e.g., IOI_CMB.txt  

Steps:
1. Put all the .margestats and .corr files of the constraints of interest in one folder.
2. Specify the constraint files directory that contains all the .margestats and .corr files.
3. Put the parameter names like below  
   e.g., H0 parameterization in LCDM model: Params = {'omegabh2','omegam','H0','sigma8','ns','tau'};    
   e.g., Theta parameterization in LCDM model: Params = {'omegabh2','omegach2','theta','logA','ns','tau'};   
