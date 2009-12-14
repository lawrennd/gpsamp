function [model PropDist samples accRates] = gpmtfTestGenesAdapt2(model, TFs, simMat,  AdaptOps)
%[model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
%
% Description: Sample the parameters of the Bayesian differential equation 
%              model so that to tune/adapt the proposal distribution (number 
%              of control points, variances of Gaussian proposals etc). 
%                ** This function should always be called before 
%                   gpTFSample. When the gpTFSample will be called the 
%                   proposal distribution (previously adapted)
%                   is kept fixed. **
%
% Inputs: 
%    -- model: the structure that contains the likelihood and GP
%              parameters as well as the priors for all these
%              quantities; write 'help model' for full details. 
%    -- PropDist: a stucture that defines the functional form of the 
%                 proposal distribution. It contains the following fields
%              * qF: contains precomputed variables needed for the 
%                    definition of the poposal distribution (that uses 
%                    control variables) over the GP latent functions 
%              * qKinVars: NumOfGenes x 4 matrix that contains the
%                    variances of the Gaussian proposal distributions
%                    used to sample the logarithm of the kinetic parameters 
%                    (represented in the log space; see modelCreate)
%                    for each gene. These Gaussian proposals have diagonal 
%                    covariance matrices so as each row corresponds to each 
%                    gene, 1 column to the basal rates (B), 2 column to 
%                    the decay rates (D), 3 column to the sensitivities(S)
%                    and 4 column to the initial condition parameter (A). 
%              * qWeigVars: NumOfGenes x M matrix that contains the
%                    variances of the Gaussian proposal distribution
%                    of all the parameters that exist in the activation 
%                    function of the latent GP function (that represent the 
%                    log of the TFs). Each row of the matrix corresponds to 
%                    each gene. The number M of parameters depends on the 
%                    number of TFs used and the functional form of the 
%                    activation function. When the sigmoid activation 
%                    function is used then  M=NumOfTFs+1. When the 
%                    Michaelis-Menten that is valid only for a single 
%                    TF (NumOfTFs=1), then M = 1; the single parameter 
%                    corresponds to the gamma Michaelis-Menten constant
%              * qLengScVars: NumOfTFs x 1 vector that contains all the
%                    variances of the Gaussian proposal distributions used to 
%                    sample the logarithm of the lengthscales of the NumOfTFs
%                    different rbf GP priors (see modelCreate) 
%    -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the
%               gene expressions for all genes, all times and 
%               replicas
%    -- TimesG: The time points where gene expression are evaluated 
%    -- TimesF: The times where the GP function are evaluated
%               TimesF>>TimesG
%    -- trainOps: User defined options about the burn-in and sampling
%               iterations, visualization features and others 
%               **see demos**
%
% Outputs: 
%    -- model: see above or write 'help model' for full details. 
%              The outputed model is updated to contain the parameters 
%              values of the final MCMC iteration
%    -- PropDist: as above. the precomputations in the PropDist.qF field
%              can be different (compared to the ones given at the input) 
%              due the update of the kernel lengthscale that determine 
%              the proposal over the latent GP functions
%    -- samples: the structure that contains the samples. In contains
%              
%              
%
%    -- accRateF: acceptance rates during sampling time (after burn-in) 
%                 for the TFs
%    -- accRateKin: >>  >> for the kinetic parameters 
%    -- accRateW: >>  >> for the interaction weigths and any other 
%              parameters that exist in the activation function of 
%              the TF (e.g. the gamma in the Michaels-Menten activation)  
%    -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%              kernel 

BurnInIters = AdaptOps.Burnin; 
Iters = AdaptOps.T; 

%%% SIZES 
[NumOfGenes SizG NumReplicas] = size(model.Likelihood.Genes);
% number of genes 
NumOfGenes = model.Likelihood.numGenes;
% number of times the gene expression is evaluated
SizG = model.Likelihood.numTimes;
% number of replicas per gene 
NumReplicas = model.Likelihood.numReplicas;
% number of transcription factors
NumOfTFs = model.Likelihood.numTFs;
SizKin = size(model.Likelihood.kinetics,2);
SizF = size(model.Likelihood.TimesF,2);


% initialize the transcription factors 
if strcmp(model.constraints.replicas,'free')
    %
    % this variable will never used (its jsut to cheat the fucntino compute likelihoods
    % which will use the precomputed TFs from the Fs )
    F = zeros(NumOfTFs, SizF, NumReplicas);
    for j=1:NumOfTFs
    for r=1:NumReplicas
         gPerm = randperm(size(TFs,2));       
         ch = gPerm(1); 
         model.Likelihood.TF(:,:,r) = TFs{ch}(:,:,r);
         TFindex(j,r) = ch;
    end
    end
    %
else
    %
    F = zeros(NumOfTFs, SizF);
    for j=1:NumOfTFs
        gPerm = randperm(size(TFs,2));       
        ch = gPerm(1); 
        model.Likelihood.TF = TFs{ch};
         TFindex(j) = ch;
    end
    %
end

% fake variable
model.F = F;

% index of the TF fucntion form the training samples 
model.TFindex = TFindex; 

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.05*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.05*ones(NumOfGenes,NumOfTFs+1);

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
epsilon = 0.1;

cnt = 0;
%
% do the adaption 
while 1
%
%  
   [model PropDist samples accRates] = gpmtfTestGenesSample2(model, TFs, simMat, PropDist, AdaptOps);
 
   accRateKin = accRates.Kin;
   accRateW = accRates.W;
   accRateF = accRates.F;
   model.TFindex
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   if AdaptOps.disp == 1
   
   fprintf(1,'Acceptance Rates for GP functions\n');
   disp(accRateF);
       
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
 
   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Average likelihood value %15.8f',mean(samples.LogL));
   fprintf(1,'\n');
   end
   fprintf(1,'------------------------------- \n',cnt+1);
   
   
   
   % if you got a good acceptance rate, then stop
   if  (min(accRateKin(:))>15) & (min(accRateW(:))>15) 
        disp('END OF ADAPTION: acceptance rates OK');
        %pause
        break;
   end
    
   cnt = cnt + 1;
   % do not allow more than 50 iterations when you adapt the proposal distribution
   if cnt == 50
       warning('END OF ADAPTION: acceptance rates were not all OK');
       break;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) + epsilon*PropDist.kin(j,:);
         if PropDist.kin(j,1) > qKinAbove 
             PropDist.kin(j,:) = qKinAbove*ones(1,SizKin);
         end
      end
      if accRateKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) - epsilon*PropDist.kin(j,:);    
         if PropDist.kin(j,1) < qKinBelow 
             PropDist.kin(j,:) = qKinBelow*ones(1,SizKin);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
         
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateW(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) + epsilon*PropDist.W(j,:);
         if PropDist.W(j,1) > qWabove 
             PropDist.W(j,:) = qWabove*ones(1,NumOfTFs+1);
         end
      end
      if accRateW(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) - epsilon*PropDist.W(j,:);    
         if PropDist.W(j,1) < qWbelow 
             PropDist.W(j,:) = qWbelow*ones(1,NumOfTFs+1);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   
%
%
end


