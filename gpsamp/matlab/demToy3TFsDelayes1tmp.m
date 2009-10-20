

dataName = 'toy3TfsDelays';
expNo = 1;

% load gene expression data for all replicas as well as the time points 
% (the expressions of each replica are currently assumed to have been obtained at the same
%  time points)
[Genes, TimesG, GeneVars] = loadDataset(dataName);

% Set up model
options = gpsampOptions('multipleTF'); % multiple TF or single TF 
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = 100;
options.scale2var1 = 1; % scale data to have variance 1
%options.tieParam = 'tied';



%%%%%%%%%%%% START Of USER SPECIFIED OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMBER OF TFS
NumOfTFs=3;
% Next we define the activation function that transforms the GP functions. 
% For multiple TFs the GP functions are linearly mixed and then are passed 
% through a sigmoid function
% for multiple TFs USE ONLY 'lin' and 'sigmoid'. The rest options work only 
% for the single TF case (e.g. 'exp' plus 'michMenten')
% - INDIVIDUAL transform of the GP function
modelOp.TFsingleAct = 'loglogOnePlusExp'; %'exp' or 'lin' or 'loglogOnePlusExp'
% - ADDITIONAL/JOINT transform by passing through a sigmoid type of function.
% - ('michMenten' assume that NumOfTFs=1)
modelOp.TFjointAct = 'sigmoid';% 'lin' or 'michMenten' or 'sigmoid'
modelOp.TFjointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized

% If you have chosen modelOp.TFjointAct = 'michMenten', then you have 
% have to specify if the TFs acts as repressor or activator 
%if strcmp(modelOp.TFjointAct,'michMenten') == 1
%    modelOp.Net_X = 1% 1 for activation, -1 for repression  
%% initial connectivity of the network expressed by a NumOfGenes x NumOfTFs
%% matrix. Each element of the matrix can take 
%% 3 values (-1:repression, 0:non-regulation, 1:activation)
%% (This feature is only used when modelOp.TFjointAct = 'michMenten'
%   initX = randn(NumOfGenes,NumOfTFs);
%   initX(initX<=0)=-1;
%   initX(initX>0)=1;
%   modelOp.Net_X = initX; % randn(NumOfGenes,NumOfTFs);
%% if 'yes' the connectivity network will be sampled 
%   modelOp.Net_Learn = 'no'; 
%end

% maximum value for the gene delays   
modelOp.tau_max = 3;

% CONSTRAINTS of the initial value of the TF at time t=0. 
% if 0, then the TF has concentration 0 at time t=0
% if 1, then the TF is free to take any value at time t=0
Constraints.Ft0 = [0 0 0];  % it must have NumOfTFs elements
% CONSTRAINTS of the initial conditions of each differential equation 
% for each gene  
% if 0, then the ODE initial cond for the jth gene is 0 at time t=0
% if 1, then the ODE initial cond is free to take any value at time t=0
Constraints.InitialConds = ones(1,NumOfGenes);
% CONSTRAINTS coming from side information about which TFs do not regulate 
% certain genes. This means that certain values in the interaction 
% matrix W are constrained to be zero. 
% if w(i,j)=1 (for i gene and j TF) then the interaction is 
% allowed and the weight w(i,j) is learned
% if w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
% constrained to be zero. 
Constraints.W = Net_X; %ones(NumOfGenes,NumOfTFs);
%Constraints.W(:,3) = ones(NumOfGenes,1);

% Nosie Variances if you have Gene Variances from PUMA package
Gvariances = 'yes';% 'yes' or 'no'. Do you know the gene variance values
                   %(not stds) for each time point and gene and replica?. 
                   % If yes, store these variances into a matrix GeneVars 
                   % that has the same size the variable Genes.   

% MCMC OPTIONS (you can change these options if you like)   
trainOps.StoreEvery = 100; % store samples after burn in  every StoreEvery iterations
trainOps.Burnin = 5000;  % burn in time
trainOps.T = 50000; % sampling time
trainOps.Store = 0;  % store the results regularly in a file 
if trainOps.Store == 1
    % file name where to store the resuls
    ok = date; 
    trainOps.ResFile = [storeFile TFActivStruct.TFsingleAct TFActivStruct.TFjointAct TFActivStruct.TFjointActBin ok];
end 
trainOps.disp = 0; % display training by giving plot during sampling and providing 
                   % information about the acceptance rate                   
%%%%%%%%%%%% END Of USER SPECIFIED OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DISCRETIZE the GP function
% (around 10 times more the number the discrete time points
% we have gene expressions)

% TimesF discretizes the [-tau_max, T] where T = max(TimesG)
% - First discretize in [0,T] (where we have observed genes) 
discr=10;
% the GPs functions is going to be descretize in  discr*(size(TimesG,2)-1)+1
% points. This number of points must be a odd number 
if (discr*(size(TimesG,2)-1))+1 > 200 
   discr = floor(200/(size(TimesG,2)-1)) + 1; 
end 
if mod(discr*(size(TimesG,2)-1)+1,2) == 0
   discr = discr-1;
end 
step = ((max(TimesG) - min(TimesG))/(size(TimesG(:),1)-1))/discr;
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
GroundTr = 'yes';
if strcmp(GroundTr,'yes')
model.GroundTruth.W0 = W0;
model.GroundTruth.W = W;
model.GroundTruth.kinetics = Kinetics;
model.GroundTruth.X = Net_X;
%model.GroundTruth.F = F(:,31:end);
model.GroundTruth.F = F;
model.GroundTruth.Taus = LikParams.Taus; % delays
model.GroundTruth.sigmas = LikParams.sigmas; 
end

if 0
model.Likelihood.W = W;
model.Likelihood.W0 = W0;
end
if 0
model.Likelihood.kinetics = Kinetics;
end


% options for the adaptive phase in MCMC (you could change this, although not recommended) 
AdaptOps.T = 200;          AdaptOps.Burnin = 100;
AdaptOps.StoreEvery = 10;  AdaptOps.disp = 1;
AdaptOps.InitialNOfControlPnts = 3; AdaptOps.IncreaseSize = 1;
% adapt the proposal distribution
[model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFAdapt(model, Genes, TimesG, TimesF, AdaptOps);
% sample  
[model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFSample(model, PropDist, Genes, TimesG, TimesF, trainOps);

% Plot the results (where time starts from 1);
plotCreate(model, samples, Genes, TimesG+1, TimesF+1, nameFile, 0, GeneVars)