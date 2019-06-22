%%%%%% A script for CosmoMC (after running getdist)
%%%%%% Calculate two-experiment IOIs for arbitary number of contraints
%%%%%% and in an arbitary parameter space
%%%%%% outputs in e.g., IOI_CMB.txt
%%%%%% Steps:
%%%%%% 0. Put all the .margestats and .corr files of the constraints of
%%%%%%    interest in one folder.
%%%%%% 1. Specify the constraint files directory that contains all 
%%%%%%    the .margestats and .corr files.
%%%%%% 2. Put the parameter names below 
%%%%%%    e.g., H0 parameterization in LCDM model:
%%%%%%      Params = {'omegabh2','omegam','H0','sigma8','ns','tau'};
%%%%%%    e.g., Theta parameterization in LCDM model:
%%%%%%      Params = {'omegabh2','omegach2','theta','logA','ns','tau'}; 

clear;
Params = {'omegabh2','omegach2','theta','logA','ns'};
constraint_filedir = './batch/';

Outfile = 'IOIouts/IOI.txt';

margfiles = dir(fullfile(constraint_filedir, '*.margestats'));
Num_exp = length(margfiles);
if Num_exp<2
    ErrorMessage = sprintf(['Error: \n' ...
                'At least two constraints are required.\n']);
    disp(ErrorMessage);
    return        
end

ParamDim = length(Params);

exp_names = string(Num_exp);
for i=1:Num_exp
    exp_names(i) = erase(margfiles(i).name,'.margestats');
end

Message = sprintf([num2str(Num_exp) ' constraints, ' ...
    num2str(ParamDim) ' parameters.']);
disp(Message)



%%%%%% Find the common parameters
for i = 1:Num_exp
    fileID = fopen([constraint_filedir margfiles(i).name]);
    Marg_header = fgets(fileID);
    Marg_header = fgets(fileID);
    Marg_header = fgets(fileID);
    All_params = textscan(fileID,'%s %*[^\n]');
    All_params_without_star = All_params{1};
    for index_temp = 1:length(All_params_without_star)
         All_params_without_star(index_temp)= ... 
             erase(All_params_without_star(index_temp),"*");
    end
    if i==1
        Common_Params = All_params_without_star;
    end
    common_index = ismember(Common_Params, All_params_without_star);
    Common_Params = Common_Params(common_index);
end

Num_com_param = length(Common_Params);
txt = ['There are ' num2str(Num_com_param) ...
    ' common parameters (including derived),'...
    ' which are stored in variable Common_Params'];
disp(txt)



%%%%%% Extract mu and C from files
C = zeros(ParamDim,ParamDim,Num_exp);
mu = zeros(ParamDim,Num_exp);
sigma = zeros(ParamDim,Num_exp);
index = zeros(ParamDim,1);

delimiterIn = ' ';
headerlinesIn = 1;

for i = 1:Num_exp
    fileID = fopen([constraint_filedir margfiles(i).name]);
    Marg_header = fgets(fileID);
    Marg_header = fgets(fileID);
    Marg_header = fgets(fileID);
    meat = textscan(fileID,'%s %f %f %*[^\n]');
    for k = 1:ParamDim
        index(k) = 1;
        str=meat{1};
        NotFound = true;
        for ii=1:length(meat{1})
            if (strcmp(str{ii},Params{k}) == 1|strcmp(str{ii},[Params{k},'*']) == 1)
                index(k) = ii;
                NotFound = false;
            end
        end
        if NotFound == true
            ErrorMessage = sprintf(['Error: \n' ...
                Params{k} ' is not in experiment: ' exp_names{i}]);
            disp(ErrorMessage)
            return
        end
        mu(k,i) = meat{2}(index(k));    
        sigma(k,i) = meat{3}(index(k));   
    end
    fclose(fileID);
    
    Corr = importdata([constraint_filedir exp_names{i} '.corr']);
    offset = 0;
    for ii=1:length(Corr)
        if Corr(ii-offset,ii-offset)==0.0
           Corr(ii-offset,:) = [];
           Corr(:,ii-offset)  = [];
           offset = offset+1;
        end
    end   
    Corr_select = Corr([index],[index]);
    for n=1:ParamDim
        for m=1:ParamDim
            C(n,m,i) = sigma(n,i).*Corr_select(n,m)*sigma(m,i);
        end
    end
end


 
%%%%%% calculating all the two-experiment IOIs
IOI = zeros(Num_exp,Num_exp-1);
for i = 1: Num_exp
    for j = (i+1):Num_exp
        IOI(i,j) = 0.5*(mu(:,i)-mu(:,j))'*(C(:,:,i)+C(:,:,j))^-1*(mu(:,i)-mu(:,j));
        IOI(j,i) = IOI(i,j);
    end
end


%%%%%% Saving the results of the two-data IOIs in a matrix for
fileID = fopen(Outfile,'w');
fprintf(fileID,'%18s',' ');
fprintf(fileID,'%18s',exp_names{:});
for i=1:Num_exp
    fprintf(fileID,'\n');
    fprintf(fileID,'%18s',exp_names{i});
    fprintf(fileID,'%18.2f',round(IOI(i,:),2));
end
fprintf(fileID,'\n\n\n\n\n\n');



%%%%%% calculating the multi-experiment IOIs 

F_mul = zeros(ParamDim,ParamDim);
mu_mul_pre = zeros(ParamDim,1);
IOI_mul_pre = 0;
for i = 1: Num_exp
    F_mul = F_mul + C(:,:,i)^-1;
    mu_mul_pre = mu_mul_pre + C(:,:,i)^-1*mu(:,i);
    IOI_mul_pre = IOI_mul_pre + mu(:,i)'*C(:,:,i)^-1*mu(:,i);
end
C_mul = F_mul^-1;
mu_mul = C_mul*mu_mul_pre;
IOI_mul = (IOI_mul_pre - mu_mul'*F_mul*mu_mul)/Num_exp;

fprintf(fileID,'%25s', 'Multi-IOI');
fprintf(fileID,'%18.2f\n\n\n',round(IOI_mul,2));


%%%%%% calculate the multi-IOIs by successively taking out each expriment 
IOI_rmj = zeros(Num_exp,1);
Outlierness = zeros(Num_exp,1);
fprintf(fileID,'%25s%18s%18s\n', 'Removing','Remaining IOI','Outlierness');
if Num_exp > 2
    for j = 1: Num_exp
        F_rmj = F_mul - C(:,:,j)^-1;
        C_rmj = F_rmj^-1;
        mu_rmj = C_rmj*(mu_mul_pre - C(:,:,j)^-1*mu(:,j));
        IOI_rmj_pre = IOI_mul_pre - mu(:,j)'*C(:,:,j)^-1*mu(:,j);
        IOI_rmj(j) = (IOI_rmj_pre - mu_rmj'*F_rmj*mu_rmj)/(Num_exp-1); 
        Outlierness(j) = (Num_exp*IOI_mul-(Num_exp-1)*IOI_rmj(j))/2;
    end
    [IOI_rmj_sort, sort_index] = sort(IOI_rmj);
    for j = 1: Num_exp
        fprintf(fileID,'%25s', exp_names{sort_index(j)});
        fprintf(fileID,'%18.2f',round(IOI_rmj_sort(j),2));
        fprintf(fileID,'%18.2f\n',round(Outlierness(sort_index(j)),2));
    end
end


%%%%%% Finishing message:
FishingMessage = sprintf(['Finished: \n' ...
    'Two-experiment and multi-experiment IOIs have been saved in ' ...
    Outfile '.\n']);
disp(FishingMessage)
fclose(fileID);

