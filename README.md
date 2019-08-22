#### IOI_CosmoMC

A matlab script for CosmoMC (after running getdist) to calculate two-dataset and multi-dataset IOIs, as well as all "outlier index" for an arbitary number of contraints and in an arbitary parameter space. Outputs are stored in e.g., IOI_CMB.txt  

Steps:
1. Put all the .margestats and .corr files of the constraints of interest in one folder.
2. Specify the directory that contains all the .margestats and .corr files.
3. Put the parameter names like below  
    e.g., H0 parameterization in LCDM model: Params = {'omegabh2','omegam','H0','sigma8','ns','tau'}   
    e.g., Theta parameterization in LCDM model: Params = {'omegabh2','omegach2','theta','logA','ns','tau'}   
    
Usage for `pyioi`:
python pyioi.py -r <Folder> -p <Paramter1> <Parameter2> ...

For questions, contact wlin23@ncsu.edu
