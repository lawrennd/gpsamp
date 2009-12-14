function model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options)
%model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options)
%
% Description: Creates the structure variable that contains all data, the
%              parameters, the priors over those parameters, the form 
%              of the likelihood etc
%
% Inputs: 
%     -- Genes:  NumOfGenes x NumOfTimes x Replicas matrix that stores the 
%            gene expressions for all genes in all times and replicas
%     -- GenesVar:  NumOfGenes x NumOfTimes x Replicas matrix
%            that specifies the variances of the gene expressions 
%            in the Gaussian likelihood (a separate variance is given 
%            for each element in the Genes). This is useful when 
%            GeneVars have been estimated previously from a package like 
%            PUMA (http://www.bioinf.manchester.ac.uk/resources/puma/).  
%            When this input argument is not given, then GeneVars is 
%            considered a random variable that needs to be inferred by MCMC. 
%            In such case a single variance is estimated for each gene, 
%            i.e. GeneVars(i,:,:) = sigma2_i * ones(NumOfTimes,Replicas)
%            and we sample sigma2_i
%     -- GenesTF:  NumOfTFs x NumOfTimes x Replicas matrix that stores the 
%            gene expressions for the TF genes (if given, otherwize is empty)
%     -- GenesTFVar:  NumOfTFs x NumOfTimes x Replicas matrix
%            that specifies the variances of the TF gene expressions similarly 
%            to GenesVar (if given, otherwize is empty).
%     -- TimesG: The time points where gene expressions are evaluated 
%     -- TimesF: The times where the GP functions (TFs) are evaluated.
%            Note that size(TimesF)>>size(TimesG).
%     -- options: Determines the functional form of the activation 
%            of the TFs in the differential equation. It contains 3 fields   
%            * singleAct = 'lin' or 'exp' or 'logOnePlusExp'. This specifies the
%                  individual activation of the GP functions   
%            * jointAct = 'genHill' or 'sigmoid' or 'michMenten' 
%                          (see the paper for definition)
%            * jointActBin = 0 or 1. Binarize the outputs of the joint
%                  activation function. If 0 the outputs are not binarized.
%                  If 1 the outputs are binarized.
%                  !IMPORTANT NOTE! singleAct = 'exp' should be used only 
%                  with the jointAct='michMentenAct' or 'michMentenRepres'
%                  but not for TFjointAct = 'sigmoid' (SEE DEMOS)  
%            * Constraints: Determines several constraints of the parameter during
%                   sampling. It contains 3 fields
%               - Ft0: 1 x NumOfTFs  (number of TFs; see below) binary vector that 
%                  detetmines is the initial value of the TF is 0 or not: 
%                  if Ft0(j)=0, then the j TF has concentration 0 at time t=0
%                  if Ft0(j)=1, then the j TF is free to take any value
%               - InitialConds: 1 x NumOfGenes binary vector that decides if the
%                  initial conditions of each differential equation for each gene  
%                  are constrained to be zero at the initial time (say t=0)  
%                  if InitialConds(j)=0, then the ODE initial cond for the jth 
%                  gene is 0. Otheriwise if is free to take any value at time t=0
%                - W:  NumOfGenes x NumOfTFs binary matrix that determines constraints 
%                  coming from side information about which TFs do not regulate 
%                  certain genes. This means that certain values in the interaction 
%                  matrix W are constrained to be zero. if w(i,j)=1 (for i gene and j TF) 
%                  then the interaction is allowed and the weight w(i,j) is learned
%                  If w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
%                  constrained to be zero. 
% 
%                      
% Outputs: 
%     -- model: The model stucture and parameters. In consists of 3 fields 
%            * Likelihood: The likelihood structure (functional form of the 
%                  activation of the TFs; see options input argument) 
%                  and parameters, e.g. kinetics and observation noise 
%            * GP: The GP prior type and kernel parameters. Only the 
%                  Gaussian or 'rbf' kernel function with unit variance 
%                  and unknown lengthscale 'ell^2' is currently supported:
%
%                    cov = exp[-0.5 * (x_i - x_j)^2/ell^2]
%
%                  The lengthscale ell^2 is encoded as 
%                  logth = 0.5*log(ell^2) = log(ell) and a Gaussian prior 
%                  is placed on the logth variable   
%
%            * prior: the prior distribution over all parameters 
%            * constraints: Specifies the status of several variables
%                  that can be constrained ('free') or 'fixed' during 
%                  sampling 
%      

% create the Likelihood model structrure
% differential equation model
model.Likelihood.type = 'diffeq';
% number of genes
[NumOfGenes NumOfTimes NumOfReplicas] = size(Genes);

% number of genes 
model.Likelihood.numGenes = NumOfGenes;
% number of times the gene expression is evaluated
model.Likelihood.numTimes = NumOfTimes;
% number of replicas per gene 
model.Likelihood.numReplicas = NumOfReplicas;
% number of transcription factors
model.Likelihood.numTFs = options.numTFs;

model.Likelihood.Genes = Genes;
% gene variances
if ~isempty(GenesVar)
   model.Likelihood.sigmas = GenesVar;
   model.constraints.sigmas = 'fixed';  % 'free' or 'fixed'
else
   model.Likelihood.sigmas = 0.05*ones(NumOfGenes, NumOfTimes, NumOfReplicas);
   model.constraints.sigmas = 'free';  % 'free' or 'fixed'
end


if ~isempty(GenesTF)
   % 
   model.Likelihood.GenesTF = GenesTF;
   if ~isempty(GenesTFVar)
      model.Likelihood.sigmasTF = GenesTFVar; 
      model.constraints.sigmasTF = 'fixed';  % 'free' or 'fixed'
   else
      model.Likelihood.sigmasTF = 0.05*ones(options.numTFs, NumOfTimes, NumOfReplicas);
      model.constraints.sigmasTF = 'free';  % 'free' or 'fixed'
   end   
   %
end


model.Likelihood.TimesG = TimesG;
model.Likelihood.TimesF = TimesF;
model.Likelihood.step = TimesF(2)-TimesF(1);

model.Likelihood.kineticsReadme = '1st column: basals, 2st decays, 3st sensitivities, 4st initial conds';
% decay rates
model.Likelihood.kinetics(:,1) = rand(NumOfGenes,1) + 0.5; %B: basal rates 
% sensitivity parameters
model.Likelihood.kinetics(:,2) = rand(NumOfGenes,1) + 0.5; %D: decay rates
% basal rates
model.Likelihood.kinetics(:,3) = rand(NumOfGenes,1) + 0.5; %S: sensitivities 
% initial conditions
model.Likelihood.kinetics(:,4) = model.Likelihood.kinetics(:,1)./model.Likelihood.kinetics(:,2);  % A: initial conditions
% Delays 
model.Likelihood.Tausreadme = 'taus are the delays in the ODEs'; 
model.Likelihood.Taus = zeros(1,NumOfGenes); 
model.Likelihood.Tausindex = (size(TimesF,2) - options.sizTime + 1)*ones(1,NumOfGenes); 
model.Likelihood.sizTime = options.sizTime;
model.Likelihood.startTime = (size(TimesF,2) - options.sizTime + 1);
model.Likelihood.tauMax = options.tauMax;

uu = model.Likelihood.TimesF(model.Likelihood.startTime:end);
[commonSlots, comInds] = intersect(uu,model.Likelihood.TimesG);
model.Likelihood.comInds = comInds; 

uu = model.Likelihood.TimesF;
[commonSlots, comIndsTF] = intersect(uu,model.Likelihood.TimesG);
model.Likelihood.comIndsTF = comIndsTF;

%model.Likelihood.numKins_perGene = 4;

%  additional kinetic parameters for the TF Genes
if ~isempty(GenesTF)
   %
   model.Likelihood.kineticsTFReadme = '1st column: decays, 2st sensitivities';
   model.Likelihood.kineticsTF(:,1) = 0.5*ones(options.numTFs,1);%rand(options.numTFs,1) + 0.5; %d: decay rates
   model.Likelihood.kineticsTF(:,2) = ones(options.numTFs,1);%rand(options.numTFs,1) + 0.5; %s: sensitivities
   %
end

% Gene - TFs interaction weights 
% (relevant only for multiple transcription factors)
model.Likelihood.W = 0.05*randn(NumOfGenes,options.numTFs);
% gene-specific bias terms for the Gene - TFs interaction 
% functions
model.Likelihood.W0 = 0.05*randn(NumOfGenes,1);


% single activation
model.Likelihood.singleAct = options.singleAct;
% valid values of jointAct are : 'sigmoid', 'michMenten'
model.Likelihood.jointAct = options.jointAct;
% binarize or not the joint activation function 
model.Likelihood.jointActBin = options.jointActBin;


% *** This feature is only used when 'michMenten' joint activation is
%     considered. When  the 'sigmoid' joint activation is used 
%     this feature remains inactive        
% ***
% -1: repression, 0 no regulation, 1 for activation 
% (default values are to assume that all TFs are activators) 
model.Likelihood.Net_Learn = 'no';
if strcmp(model.Likelihood.jointAct,'michMenten');
model.Likelihood.Net_X = options.Net_X;  
%model.Likelihood.Net_Learn = options.Net_Learn;
end

% this may be used to stored precomputed TF profiles 
% (for saving computations purposes)
model.Likelihood.TF = [];

% create the GP prior model 
timescale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
lengthscale = (max(TimesG(:))-min(TimesG(:)))/10;%(size(TimesG,2)-1);
%lengthscale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
for j=1:options.numTFs
   model.GP{j}.type = {'rbf','white'};
   model.GP{j}.TF = j;
   model.GP{j}.logtheta = [log(lengthscale) 0 0.5*log(1e-06)];
end
X = TimesF(:);
[n D] = size(X);
model.GP{1}.X2 = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);


% define the prior of the TFs regulation types 
% (repression, non-regulation, activation)
                       % repression, non-regulation, activation
model.prior.Net.type = 'discrete';
model.prior.Net.contraint = 'probability';
model.prior.Net.priorSpace = 'lin'; % it means NoTransform;
model.prior.Net.prob = [0.025 0.95 0.025];
model.prior.Net.readme = 'The 3 prior probabilities correspond to: repression, non-regulation, activation';


% prior for the kinetics parameters
model.prior.kinetics.type = 'normal';
%model.prior.kinetics.readme = 'kinetics are: decays, sensitivities, basals, initial conds'
model.prior.kinetics.contraint = 'positive';
model.prior.kinetics.priorSpace = 'log';
model.prior.kinetics.a = -0.5; % mean 
model.prior.kinetics.b = 2; % variance

model.prior.delays.type = 'discrete'; 
a = 1; 
AllTaus = model.Likelihood.tauMax:model.Likelihood.step:-eps;
AllTaus = [AllTaus, 0]; 
model.prior.delays.prob = exp(a*AllTaus)/(sum(exp(a*AllTaus)));

model.prior.weights.type = 'normal'; %'normal' or Laplace
if strcmp(model.Likelihood.jointAct,'michMenten')
model.prior.weights.constraint = 'positive';
model.prior.weights.priorSpace = 'log';
model.Likelihood.W = rand(NumOfGenes,options.numTFs)+0.1;
model.Likelihood.W0 = rand(NumOfGenes,1)+0.1;
else
model.prior.weights.constraint = 'real';
model.prior.weights.priorSpace = 'lin'; % it means NoTransform;   
end
model.prior.weights.mu = 0;
model.prior.weights.sigma2 = 1.5;

% prior for the gene specifc noise variances (the inverce of them) 
model.prior.invsigma2.type = 'gamma';
model.prior.invsigma2.constraint = 'positive';
model.prior.invsigma2.priorSpace = 'lin'; % it means NoTransform;
model.prior.invsigma2.a = 0.1;
model.prior.invsigma2.b = 0.01;

% prior for the lengthscales
model.prior.GPkernel.lenghtScale.type = 'normal';
model.prior.GPkernel.lenghtScale.constraint = 'positive';
model.prior.GPkernel.lenghtScale.priorSpace = 'log';
ok = 1.5*(max(TimesG(:))-min(TimesG(:)))/8;%(size(TimesG,2)-1);
model.prior.GPkernel.lenghtScale.a = 2*log(ok); % mean 
model.prior.GPkernel.lenghtScale.b = 2;         % variance
%model.prior.GPkernel.lenghtScale.a = 1;
%model.prior.GPkernel.lenghtScale.b = 1/(ok^2);
 

% initial value of the TFs
% if 0, then the TF has concentration 0 at time t=0
% if 1, then the TF is free to take any value at time t=0
model.constraints.Ft0 = options.constraints.Ft0;
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.Ft0 = 'fixed')
model.constraints.Ft0_value = 0; 

% initial condition of the differential equation 
% if 0, then the ODE initial cond for the jth gene is 0 at time t=0
% if 1, then the ODE initial cond is free to take any value at time t=0
model.constraints.InitialConds = options.constraints.initialConds; 
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.InitialConds = 0)
model.constraints.InitialConds_value = 0; 

%
model.constraints.geneTFsensitivity = options.constraints.geneTFsensitivity; 
model.constraints.geneTFsensitivity_value = 1; 

% constraints on the interaction weigths between TF and genes 
model.constraints.W = options.constraints.X;

% constraint on TF functions across different replicas  
model.constraints.replicas = options.constraints.replicas;
