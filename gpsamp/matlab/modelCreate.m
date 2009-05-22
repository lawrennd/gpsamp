function model = modelCreate(modelOp, Constraints, Genes, NumOfTFs, TimesG, TimesF, GeneVars)
%model = modelCreate(modelOp, Constraints, Genes, NumOfTFs, TimesG, TimesF, GeneVars)
%
% Description: Creates the structure variable that contains all the
%              parameters, the priors over those parameters, the form 
%              of the likelihood etc
%
% Inputs: 
%     -- modelOp: Determines the functional form of the activation 
%            of the TFs in the differential equation. It contains 3 fields   
%            * TFsingleAct = 'exp' or 'linear'. This specifies the
%                  individual activation of the GP functions   
%            * TFjointAct = 'sigmoid' or 'michMentenAct' or 'michMentenRepres'
%                          (see the paper for definition)
%            * TFjointActBin = 0 or 1. Binarize the outputs of the joint
%                  activation function. If 0 the outputs are not binarized.
%                  If 1 the outputs are binarized.
%                  !IMPORTANT NOTE! TFsingleAct = 'exp' should be used only 
%                  with the TFjointAct='michMentenAct' or 'michMentenRepres'
%                  but not for TFjointAct = 'sigmoid' (SEE DEMOS)  
%
%     -- Constraints: Determines several constraints of the parameter during
%            sampling. It contains 3 fields
%            * Ft0: 1 x NumOfTFs  (number of TFs; see below) binary vector that 
%                  detetmines is the initial value of the TF is 0 or not: 
%                  if Ft0(j)=0, then the j TF has concentration 0 at time t=0
%                  if Ft0(j)=1, then the j TF is free to take any value
%            * InitialConds: 1 x NumOfGenes binary vector that decides if the
%                  initial conditions of each differential equation for each gene  
%                  are constrained to be zero at the initial time (say t=0)  
%                  if InitialConds(j)=0, then the ODE initial cond for the jth 
%                  gene is 0. Otheriwise if is free to take any value at time t=0
%            * W:  NumOfGenes x NumOfTFs binary matrix that determines constraints 
%                  coming from side information about which TFs do not regulate 
%                  certain genes. This means that certain values in the interaction 
%                  matrix W are constrained to be zero. if w(i,j)=1 (for i gene and j TF) 
%                  then the interaction is allowed and the weight w(i,j) is learned
%                  If w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
%                  constrained to be zero. 
%
%     -- Genes:  NumOfGenes x NumOfTimes x Replicas matrix that stores the 
%            gene expressions for all genes in all times and replicas
%     -- NumOfTFs: A positive integer ( >=1) that is the number of TFs used 
%            in the differential equation 
%     -- TimesG: The time points where gene expressions are evaluated 
%     -- TimesF: The times where the GP functions (TFs) are evaluated.
%            Note that size(TimesF)>>size(TimesG).
%     -- GeneVars(OPTIONAL):  NumOfGenes x NumOfTimes x Replicas matrix
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
%                      
% Outputs: 
%     -- model: The model stucture and parameters. In consists of 3 fields 
%            * Likelihood: The likelihood structure (functional form of the 
%                  activation of the TFs; see modelOp input argument) 
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
%            * Constraints: Specifies the status of several variables
%                  that can be constrained ('free') or 'fixed' during 
%                  sampling 
%      

% create the Likelihood model structrure
% differential equation model
model.Likelihood.type = 'diffeq';
% number of genes
[NumOfGenes NumOfTimes NumOfReplicas] = size(Genes);

% number of genes 
model.Likelihood.NumOfGenes = NumOfGenes;
% number of times the gene expression is evaluated
model.Likelihood.NumOfTimes = NumOfTimes;
% number of replicas per gene 
model.Likelihood.NumOfReplicas = NumOfReplicas;
% number of transcription factors
model.Likelihood.NumOfTFs = NumOfTFs;

% decay rates
model.Likelihood.D = rand(1,NumOfGenes) + 0.5;
% sensitivity parameters
model.Likelihood.S = rand(1,NumOfGenes) + 0.5;
% basal rates
model.Likelihood.B = rand(1,NumOfGenes) + 0.5;
% initial conditions
model.Likelihood.A = model.Likelihood.B./model.Likelihood.D;
model.Likelihood.NumOfKins_perGene = 4;

% Gene - TFs interaction weights 
% (relevant only for multiple transcription factors)
model.Likelihood.W = 0.05*randn(NumOfGenes,NumOfTFs);
% gene-specific bias terms for the Gene - TFs interaction 
% functions
model.Likelihood.W0 = 0.05*randn(NumOfGenes,1);


% valid values of TFsingleAct are : 'lin', 'exp'
model.Likelihood.TFsingleAct = modelOp.TFsingleAct;
% valid values of TFjointAct are : 'sigmoid', 'michMenten'
model.Likelihood.TFjointAct = modelOp.TFjointAct;
% binarize or not the joint activation function 
model.Likelihood.TFjointActBin = modelOp.TFjointActBin;


% *** This feature is only used when 'michMenten' joint activation is
%     considered. When  the 'sigmoid' joint activation is used 
%     this feature remains inactive        
% ***
% -1: repression, 0 no regulation, 1 for activation 
% (default values are to assume that all TFs are activators) 
model.Likelihood.Net_Learn = 'no';
if strcmp(model.Likelihood.TFjointAct,'michMenten');
model.Likelihood.Net_X = modelOp.Net_X;  
model.Likelihood.Net_Learn = modelOp.Net_Learn;
end


% gene variances
if nargin == 7
   model.Likelihood.sigmas = GeneVars;
   model.Constraints.sigma2 = 'fixed';  % 'free' or 'fixed'
else
   model.Likelihood.sigmas = 0.05*ones(NumOfGenes, NumOfTimes, NumOfReplicas);
   model.Constraints.sigma2 = 'free';  % 'free' or 'fixed'
end

% define the prior of the TFs regulation types 
% (repression, non-regulation, activation)
                       % repression, non-regulation, activation
model.prior.Net.type = '3-valued multinomial';
model.prior.Net.contraint = 'probability';
model.prior.Net.priorSpace = 'lin'; % it means NoTransform;
model.prior.Net.prob = [0.025 0.95 0.025];
model.prior.Net.readme = 'The 3 prior probabilities correspond to: repression, non-regulation, activation';


% prior for the kinetcis parameters
model.prior.kinetics.type = 'normal';
model.prior.kinetics.contraint = 'positive';
model.prior.kinetics.priorSpace = 'log';
model.prior.kinetics.a = 0; % mean 
model.prior.kinetics.b = 2; % variance

model.prior.weights.type = 'normal'; %'normal' or Laplace
if strcmp(model.Likelihood.TFjointAct,'michMenten')
model.prior.weights.constraint = 'positive';
model.prior.weights.priorSpace = 'log';
model.Likelihood.W = rand(NumOfGenes,NumOfTFs)+0.1;
model.Likelihood.W0 = rand(NumOfGenes,1)+0.1;
else
model.prior.weights.constraint = 'real';
model.prior.weights.priorSpace = 'lin'; % it means NoTransform;   
end
model.prior.weights.mu = 0;
model.prior.weights.sigma2 = 1.5;

% prior for the gene sepcfic nosie variances (the invernce of them) 
 model.prior.invsigma2.type = 'gamma';
model.prior.invsigma2.constraint = 'positive';
model.prior.invsigma2.priorSpace = 'lin'; % it means NoTransform;
model.prior.invsigma2.a = 1;
model.prior.invsigma2.b = 0.1;

% prior for the lengthscales
model.prior.GPkernel.lenghtScale.type = 'normal';
model.prior.GPkernel.lenghtScale.constraint = 'positive';
model.prior.GPkernel.lenghtScale.priorSpace = 'log';
ok = 1.5*(max(TimesG(:))-min(TimesG(:)))/8;%(size(TimesG,2)-1);
model.prior.GPkernel.lenghtScale.a = 2*log(ok); % mean 
model.prior.GPkernel.lenghtScale.b = 2;         % variance
%model.prior.GPkernel.lenghtScale.a = 1;
%model.prior.GPkernel.lenghtScale.b = 1/(ok^2);
 
% create the GP prior model 
model.GP.type = 'rbf';
timescale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
lengthscale = (max(TimesG(:))-min(TimesG(:)))/10;%(size(TimesG,2)-1);
%lengthscale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
model.GP.logtheta = repmat([log(lengthscale) 0],NumOfTFs,1);
X = TimesF(:);
[n D] = size(X);
model.GP.X2 = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);

% initial value of the TFs
% if 0, then the TF has concentration 0 at time t=0
% if 1, then the TF is free to take any value at time t=0
model.Constraints.Ft0 = Constraints.Ft0;
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.Ft0 = 'fixed')
model.Constraints.Ft0_value = 0; 

% initial condition of the differential equation 
% if 0, then the ODE initial cond for the jth gene is 0 at time t=0
% if 1, then the ODE initial cond is free to take any value at time t=0
model.Constraints.InitialConds = Constraints.InitialConds; 
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.InitialConds = 0)
model.Constraints.InitialConds_value = 0; 
% constraints on the interaction weigths between TF and genes 
model.Constraints.W = Constraints.W;


