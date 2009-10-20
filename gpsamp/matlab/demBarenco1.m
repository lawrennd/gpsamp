%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load gene expression data for all replicals as well as the time points 
% (the expressions of each replica are currently assumed to have been obtained at the same
%  time points)
clear
load data/Barencodata.mat;
namefile = 'BarencoP53';
TimesG = times';
Genes(:,:,1) = y{1}'; 
Genes(:,:,2) = y{2}'; 
Genes(:,:,3) = y{3}';
GeneVars(:,:,1) = yvar{1}';
GeneVars(:,:,2) = yvar{2}';
GeneVars(:,:,3) = yvar{3}';
%%%%%%%%%%%%%%%%%%%%%% END LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% START Of USER SPECIFIED OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMBER OF TFS
NumOfTFs=1;
% Next we define the activation function that transforms the GP functions. 
% For multiple TFs see demDrosAntti1.m. The single TF case use 'exp' 
% plus 'michMenten'
% - INDIVIDUAL transform of the GP function
%   (this is implemetns the exp(f(t)))
modelOp.TFsingleAct = 'exp';% 'exp' or 'lin'
% - JOINT transform by passing through michMenten' assume that NumOfTFs=1)
modelOp.TFjointAct = 'michMenten';% 'lin' or 'michMenten' or 'sigmoid'
modelOp.TFjointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized
% If you have chosen modelOp.TFjointAct = 'michMenten', then you have 
% have to specify if the TFs acts as REPRESSOR or ACTIVATOR
% (you have the freedom to make differecne choice for each genes)
NumOfGenes = size(Genes,1);
if strcmp(modelOp.TFjointAct,'michMenten') == 1
    % set the TF to activate all genes
    modelOp.Net_X = ones(NumOfGenes,1); % 1 for activation, -1 for repression 
    modelOp.Net_Learn = 'no'; 
end

% maximum value for the gene delays
modelOp.tau_max = 0;

% CONSTRAINTS of the initial value of the TF at time t=0. 
% if 0, then the TF has concentration 0 at time t=0
% if 1, then the TF is free to take any value at time t=0
Constraints.Ft0 = 1;  % it must have NumOfTFs elements
% CONSTRAINTS of the initial conditions of each differential equation 
% for each gene  
% if 0, then the ODE initial cond for the jth gene is 0 at time t=0
% if 1, then the ODE initial cond is free to take any value at time t=0
Constraints.InitialConds = zeros(1,NumOfGenes);
% CONSTRAINTS coming from side information about which TFs do not regulate 
% certain genes. This means that certain values in the interaction 
% matrix W are constrained to be zero. 
% For the single TF case with Michaelis-Menten non-linearity 
% there is single weight per gene which should NOT be constrained 
% to be zero 
Constraints.W = ones(NumOfGenes,1);

% Noise Variances if you have Gene Variances from PUMA package
Gvariances = 'yes';% 'yes' or 'no'. Do you know the gene variance values
                   %(not stds) for each time point and gene and replica?. 
                   % If yes, store these variances into a matrix GeneVars 
                   % that has the same size the variable Genes.   

% MCMC OPTIONS (you can change these options if you like)   
trainOps.StoreEvery = 200; % store samples after burn in  every StoreEvery iterations
trainOps.Burnin = 10000;  % burn in time
trainOps.T = 100000; % sampling time
trainOps.Store = 0;  % store the results regularly in a file 
if trainOps.Store == 1
    % file name where to store the resuls
    ok = date; 
    trainOps.ResFile = [storeFile TFActivStruct.TFsingleAct TFActivStruct.TFjointAct TFActivStruct.TFjointActBin ok];
end 
trainOps.disp = 0; % display training by giving plot during sampling and providing 
                   % information about the acceptance rate                   
%%%%%%%%%%%% END Of USER SPECIFIED OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%  From old code
%TFActivStruct.TFsingleAct = 'exp';% or 'exp'
%TFActivStruct.TFjointAct = 'michMentenAct'; %'sigmoid';% or 'michMentenAct' or 'michMentenRepres'
%TFActivStruct.TFjointActBin = 0; % binarize the output of the joint activation function 
%                                 % to 0-1 (where the threshold is the 0.5)
%                                 % if 0 the outputs are not binarilzed
%                                 % if 1 the outputs are binarized
                                 
                                                                   
% DISCRETIZE the GP function
% (around 15 times more the number the discrete time points
% we have gene expressions) 
discr=15;
if (discr*size(TimesG,2)) > 200 
    discr = floor(200/size(TimesG,2)); 
end 
step = ((max(TimesG) - min(TimesG))/(size(TimesG(:),1)-1))/discr;
if step > (0.5*min(TimesG(2:end)-TimesG(1:end-1))/2)
    step = 0.5*min(TimesG(2:end)-TimesG(1:end-1))/2;
end
TimesF =[]; TimesF(1) = TimesG(1);
for j=1:size(TimesG(:),1)-1
    TimesF = [TimesF, ((TimesG(j)+step):step:TimesG(j+1))];
    if TimesF(end) ~= TimesG(j+1)
        TimesF = [TimesF, TimesG(j+1)];
    end
end  
% - Now discretize in [-tau_max,0) (the "delay" part of the TF) 
modelOp.sizTime = size(TimesF,2);
if abs(modelOp.tau_max) > 0
DelayP = -step:-step:-modelOp.tau_max; 
modelOp.tau_max = DelayP(end); 
TimesF = [DelayP(end:-1:1) TimesF];    
end

% CREATE the model
if strcmp(Gvariances,'yes')
model = modelCreate(modelOp, Constraints, Genes, NumOfTFs, TimesG, TimesF, GeneVars);
else
model = modelCreate(modelOp, Constraints, Genes, NumOfTFs, TimesG, TimesF);
end

%% store ground truth if available
model.GroundTr = 'no';
%if strcmp(model.GroundTr,'yes')
%model.GtrW0 = W0;
%model.GtrW = W;
%model.GtrKinetics = Kinetics;
%model.GtrX = Net_X;
%model.FF = F;
%end

% options for the adaptive phase in MCMC (you could change this, although not recommended) 
AdaptOps.T = 200;          AdaptOps.Burnin = 100;
AdaptOps.StoreEvery = 10;  AdaptOps.disp = 1;
AdaptOps.InitialNOfControlPnts = 3; AdaptOps.IncreaseSize = 1;
% adapt the proposal distribution
[model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFAdapt(model, Genes, TimesG, TimesF, AdaptOps);
% sample  
[model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFSample(model, PropDist, Genes, TimesG, TimesF, trainOps);

% Plot the results
plotCreate(model, samples, Genes, TimesG, TimesF, namefile, 0, GeneVars);

