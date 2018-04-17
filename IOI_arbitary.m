%%%%%% two-experiment IOI for arbitary number of contraints
%%%%%% and in arbitary parameter space
%%%%%% outputs in e.g., IOI_CMB.txt
clear;

%%%%%% Put below all constraints to be compared:
exp_names = {'hiTT_lowTEB' ,'hiTE_lowTEB','hiEE_lowTEB'};

%%%%%% Put the parameter names below 
% e.g., H0 parameterization in LCDM model:
%Params = {'omegabh2','omegam','H0','sigma8','ns','tau'};
%%%%%% e.g., Theta parameterization in LCDM model:
Params = {'omegabh2','omegach2','theta','logA','ns','tau'};

ParamDim = length(Params);

Num_exp = length(exp_names);

delimiterIn = ' ';
headerlinesIn = 1;

C = zeros(ParamDim,ParamDim,Num_exp);
mu = zeros(ParamDim,Num_exp);
sigma = zeros(ParamDim,Num_exp);
index = zeros(ParamDim,1);

for i = 1:Num_exp
    Corr = importdata(['../batch/' exp_names{i} '.corr']);
    offset = 0;
    for ii=1:length(Corr)
        if Corr(ii-offset,1)==0.0
           Corr(ii-offset,:) = [];
           Corr(:,ii-offset)  = [];
           offset = offset+1;
        end
    end    
    for k = 1:ParamDim
        fileID = fopen(['../batch/' exp_names{i} '.margestats']);
        Marg_header = textscan(fileID,'%s',13);
        index(k) = 1;
        for j = 1:200
            meat = textscan(fileID,'%s %f %f%*[^\n]',1);
            str =meat{1};
            if (strcmp(str,Params{k}) == 1|strcmp(str,[Params{k},'*']) == 1)
                break
            end
            index(k) = index(k)+1;
        end
        mu(k,i) = meat{2};    
        sigma(k,i) = meat{3};  
        fclose(fileID);
    end
    Corr_select = Corr([index],[index]);
    for n=1:ParamDim
        for m=1:ParamDim
            C(n,m,i) = sigma(n,i).*Corr_select(n,m)*sigma(m,i);
        end
    end
 end


IOI = zeros(Num_exp,Num_exp-1);
for i = 1: Num_exp
    for j = (i+1):Num_exp
        IOI(i,j) = 0.5*(mu(:,i)-mu(:,j))'*(C(:,:,i)+C(:,:,j))^-1*(mu(:,i)-mu(:,j));
        IOI(j,i) = IOI(i,j);
    end
end

fileID = fopen('IOI_CMB.txt','w');
fprintf(fileID,'%18s',' ');
fprintf(fileID,'%18s',exp_names{:});
for i=1:Num_exp
    fprintf(fileID,'\n');
    fprintf(fileID,'%18s',exp_names{i});
    fprintf(fileID,'%18.2f',round(IOI(i,:),2));
end
fclose(fileID);