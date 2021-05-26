function [model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
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
SizKinTF = 2;
SizF = size(model.Likelihood.TimesF,2);

% if inialize model parameters if required
if AdaptOps.initParams == 1
%
% initialize the GP functions
if strcmp(model.constraints.replicas,'free')
    %
    F = zeros(NumOfTFs, SizF, NumReplicas);
    %
else
    %
    F = zeros(NumOfTFs, SizF);
    %
end


% if the initial values of the TFs is kept fixed, then set the first 
% function value to be zero (-5 si GP function is the log of the TF)
FixedFirst = zeros(1,NumOfTFs);
for j=1:NumOfTFs
if model.constraints.Ft0(j)==0  
   FixedFirst(j) = 1;
   if model.constraints.Ft0_value == 0
       if strcmp(model.Likelihood.singleAct,'lin') & strcmp(model.Likelihood.jointAct,'sigmoid') == 0
           model.u = zeros(NumOfTFs,NumReplicas);  
       else
           if strcmp(model.Likelihood.singleAct,'truncPositiveLinear')       
              model.u = exp(-10)*ones(NumOfTFs,NumReplicas);  
           else
              model.u = -5*ones(NumOfTFs,NumReplicas);
           end
       end
   end        
end
end

model.auxLikVar = ones(size(F));
if  strcmp(model.Likelihood.singleAct,'truncPositiveLinear')   | strcmp(model.Likelihood.singleAct,'lin')
for j=1:NumOfTFs
for r=1:NumReplicas
  model.auxLikVar(j,:,r) = var(model.Likelihood.GenesTF(j,:,r))*ones(1,size(F,2));
end
end
end

% initialiaze the imaginary data
model.Z = F + sqrt(model.auxLikVar).*randn(size(F)); 

% if the initial condition of the differential equations is zero, then set accordingly the
% corresponding kinetic parameter (initial condition)
for j=1:NumOfGenes
if model.constraints.InitialConds(j) == 0
   %if model.Constraints.A_value == 0
       model.Likelihood.kinetics(:,4) = model.Likelihood.kinetics(:,1)./model.Likelihood.kinetics(:,2);
   %end        
end
end

TimesF = model.Likelihood.TimesF;
model.F = F; 
n = SizF;

for j=1:NumOfTFs
  
  model.K(:,:,j) = covfunCompute(model.GP{j}, TimesF(:));   
  
  % condition on the first point being zero if needed
  if model.constraints.Ft0(j)==0  
     model.Kfu(:,j) = covfunCompute(model.GP{j}, TimesF(:), TimesF(1));   
     model.Kuu(j) = covfunCompute(model.GP{j}, TimesF(1), TimesF(1)); 
     % update prior
     %[cmu, cSigma, KInvK] =  gaussianFastConditional(zeros(1,111), model.K(:,:,j), 2:n,1); 
     model.K(:,:,j) = model.K(:,:,j) - (model.Kfu(:,j)*model.Kfu(:,j)')/model.Kuu(j);
     model.mu(:,j) = (model.Kfu(:,j)*model.u(j,1))/model.Kuu(j);
  end
  
  % update model structure
  [L,er]=jitterChol(model.K(:,:,j));
  L = L';
  model.LogDetK(j) = 2*sum(log(diag(L)));
  model.invL(:,:,j) = L\eye(n);
  
  X = TimesF(:);
  model.GP{j}.X2 = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);       
  %L = chol(model.K + diag(model.auxLikVar));    
  %tmp = L\eye(n);
  %invK = tmp*tmp';
  %model.invKsigmaK = invK*model.K;
  for r=1:NumReplicas
    model.invKsigmaK{r}(:,:,j) = (model.K(:,:,j) + diag(model.auxLikVar(j,:,r)))\model.K(:,:,j); 
    model.auxPostL{r}(:,:,j) = jitterChol( model.K(:,:,j) - model.K(:,:,j)*model.invKsigmaK{r}(:,:,j) );
  end
end

end % if for the initialization of model parameters


% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.5*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.5*ones(NumOfGenes,NumOfTFs+1);
PropDist.Kern = 0.5*ones(NumOfTFs, 2); % 2 possible parameters: kenrle variacne and kernel lengthscale


onlyPumaVar = 1; 
if model.Likelihood.noiseModel.active(2) > 0  
   onlyPumaVar = 0; 
end

% if there are puma variances for the TF genes, then use only those
onlyPumaVarTF = 1;
if model.Likelihood.noiseModel.active(1) == 0 & model.Likelihood.noiseModel.active(2) > 0
     onlyPumaVarTF = 0;
end


if onlyPumaVar == 0 
   % white nosie variance per gene possibly added to the PUMA variances
   PropDist.noiseModel = 0.5*ones(1, NumOfGenes); 
end


% additional proposal distribution for the TF kinetic parameters
if isfield(model.Likelihood,'GenesTF')  
  PropDist.TFkin = 0.5*ones(NumOfTFs, size(model.Likelihood.kineticsTF,2));
  if onlyPumaVar == 0
     PropDist.noiseModelTF = 0.5*ones(1, NumOfTFs); 
  end
  %
end

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
qNoiseMbelow = 0.000001;  qNoiseMabove = 2;
qKernBelow = 0.0001*PropDist.Kern(1); 
qKernAbove = 2;
epsilon = 0.1;
cnt = 0;
opt = 0.25;

%
% do the adaption 
while 1
%
%  
   %tic;
   [model PropDist samples accRates] = gpmtfSampleImdata(model, PropDist, AdaptOps);
  
   model.Likelihood.kineticsTF
   samples.lengthScale(:,end)
   %samples.LogL(end)
   %samples.LogLTF(end)
   
   %toc;
   accRateF = accRates.F;
   accRateKin = accRates.Kin;
   accRateW = accRates.W; 
   accRateNoiseM = accRates.noiseM;
   accRateKern = accRates.Kern;
   %
   if isfield(model.Likelihood,'GenesTF')
       accRateTFKin = accRates.TFKin;
       if onlyPumaVarTF == 0
          accRateNoiseMTF = accRates.noiseMTF;
       end
   end
   
   if AdaptOps.disp == 1
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   fprintf(1,'Acceptance Rates for GP functions\n');
   for jj=1:NumOfTFs 
          fprintf(1,'TF function #%2d (rows: #%2d replicas) \n',jj,NumReplicas);       
          disp(accRateF(jj,:));
   end    
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
   
   if isfield(model.Likelihood,'GenesTF')
      fprintf(1,'Acceptance Rates for kinetic parameters (per TF-gene))\n');
      disp(accRateTFKin);
   end

   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Acceptance Rates for kernel hyperparameters (per GP function)\n');
   disp(accRateKern);
  
   if onlyPumaVar == 0
         fprintf(1,'Acceptance rates for the noise parameters in the likelihood\n');
         disp(accRateNoiseM);
   end
   if isfield(model.Likelihood,'GenesTF') & (onlyPumaVarTF == 0)
         fprintf(1,'Acceptance rates for the noise parameters in the TF-Genes likelihood\n');
         disp(accRateNoiseMTF);
   end   
   
   
   fprintf(1,'Average likelihood value %15.8f',mean(samples.LogL));
   if isfield(model.Likelihood,'GenesTF')
       fprintf(1,' TFGenes LogL %15.8f\n',mean(samples.LogLTF));
   else
       fprintf(1,'\n');
   end
   
   fprintf(1,'------------------------------- \n',cnt+1);
   end   
       
   if (min(accRateF(:))>20) & (min(accRateKin(:))>15) & (min(accRateW(:))>15) & (min(accRateKern(:))>15)
   %
        disp('END OF ADAPTION: acceptance rates OK');
        %pause
        break;
   end
   
    
   cnt = cnt + 1;
   % do not allow more than 100 iterations when you adapt the proposal distribution
   if cnt == AdaptOps.maxIters
       warning('END OF ADAPTION: acceptance rates were not all OK');
       break;
   end
   
   %%%%%%%%%%%% ADAPT AUXILIARY LIKELIHOOD VARIANCE %%%%%%%%%%%
   for j=1:NumOfTFs
      for r=1:NumReplicas    
      if (accRateF(j,r) > (100*(opt+0.05))) | (accRateF(j,r) < (100*(opt-0.05)))
      %    
          model.auxLikVar(j,:,r) = model.auxLikVar(j,:,r) + (epsilon*((accRateF(j,r)/100 - opt)/opt))*model.auxLikVar(j,:,r);
          
          %[min(model.auxLikVar(j,:,r)) max(model.auxLikVar(j,:,r))] 
          %L = chol(model.K + diag(model.auxLikVar));  
          %tmp = L\eye(n);
          %invK = tmp*tmp';
          %model.invKsigmaK = invK*model.K;
          model.invKsigmaK{r}(:,:,j) = (model.K(:,:,j) + diag(model.auxLikVar(j,:,r)))\model.K(:,:,j);
          %model.invKsigmaK = (model.K + diag(model.auxLikVar.*model.Psi) )\model.K;
          model.auxPostL{r}(:,:,j) = jitterChol( model.K(:,:,j) - model.K(:,:,j)*model.invKsigmaK{r}(:,:,j) );  
      %   
      end
      end
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
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT TF-GENES KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%) 
   if isfield(model.Likelihood,'GenesTF')
      for j=1:NumOfTFs
      if accRateTFKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.TFkin(j,:) = PropDist.TFkin(j,:) + epsilon*PropDist.TFkin(j,:);
         if PropDist.TFkin(j,1) > qKinAbove 
             PropDist.TFkin(j,:) = qKinAbove*ones(1,SizKinTF);
         end
      end
      if accRateTFKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.TFkin(j,:) = PropDist.TFkin(j,:) - epsilon*PropDist.TFkin(j,:);    
         if PropDist.TFkin(j,1) < qKinBelow 
             PropDist.TFkin(j,:) = qKinBelow*ones(1,SizKinTF);
         end
         %
      end
       %
      end 
      % 
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT TF-GENES KINETICS PROPOSAL %%%%%%%%%%%%%%%%%%
   
   
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
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KERNEL HYPRS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfTFs
      if accRateKern(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.Kern(j, :) = PropDist.Kern(j, :) + epsilon*PropDist.Kern(j, :);
         if PropDist.Kern(j, 1) > qKernAbove
             PropDist.Kern(j, :) = qKernAbove*ones(1, size(PropDist.Kern, 2));
         end
      end
      if accRateKern(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.Kern(j, :) = PropDist.Kern(j, :) - epsilon*PropDist.Kern(j, :);    
         if PropDist.Kern(j, 1) < qKernBelow 
             PropDist.Kern(j, :) = qKernBelow*ones(1, size(PropDist.Kern, 2));
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT KERNEL HYPRS PROPOSAL %%%%%%%%%%%%%%%%
   
   if onlyPumaVar == 0 
    %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT NOISE-MODEL PROPOSAL %%%%%%%%%%%%%%%%
      for j=1:NumOfGenes
         if accRateNoiseM(j) > 35
            % incease the covariance to reduce the acceptance rate
            PropDist.noiseModel(j) = PropDist.noiseModel(j) + epsilon*PropDist.noiseModel(j);
            if PropDist.noiseModel(j) > qNoiseMabove 
               PropDist.noiseModel(j) = qNoiseMabove;
            end
         end
         if accRateNoiseM(j) < 15
            % decrease the covariance to incease the acceptance rate
            PropDist.noiseModel(j) = PropDist.noiseModel(j) - epsilon*PropDist.noiseModel(j);    
            if PropDist.noiseModel(j) < qNoiseMbelow 
               PropDist.noiseModel(j) = qNoiseMbelow;
            end
          %
         end
       %
      end
      %
      if isfield(model.Likelihood,'GenesTF') & (onlyPumaVarTF == 0) 
      for j=1:NumOfTFs
         if accRateNoiseMTF(j) > 35
            % incease the covariance to reduce the acceptance rate
            PropDist.noiseModelTF(j) = PropDist.noiseModelTF(j) + epsilon*PropDist.noiseModelTF(j);
            if PropDist.noiseModelTF(j) > qNoiseMabove 
               PropDist.noiseModelTF(j) = qNoiseMabove;
            end
         end
         if accRateNoiseMTF(j) < 15
            % decrease the covariance to incease the acceptance rate
            PropDist.noiseModelTF(j) = PropDist.noiseModelTF(j) - epsilon*PropDist.noiseModelTF(j);    
            if PropDist.noiseModelTF(j) < qNoiseMbelow 
               PropDist.noiseModelTF(j) = qNoiseMbelow;
            end
          %
         end
       %
      end   
      end
      %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT NOISE-MODEL PROPOSAL %%%%%%%%%%%%%%%%
   
%
%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function  visualization(model, F, Fnew, Z, j)
%
%
  trform = 'lin';
  subplot(1,2,1);
  FF = feval(model.Likelihood.singleAct,F);
  FFnew = feval(model.Likelihood.singleAct,Fnew);

  plot(model.Likelihood.TimesF, feval(trform,FF),'g','lineWidth',4);
  hold on;
  if isfield(model,'groundtr') == 1
    GrFtmp = feval(model.Likelihood.singleAct,model.groundtr.F(j,:));
    plot(model.Likelihood.TimesF, feval(trform,GrFtmp),'k','lineWidth',4);
  end
  pause(0.3);   
  set(gca,'FontSize',16);
  plot(model.Likelihood.TimesF, feval(model.Likelihood.singleAct, Fnew), '--b', 'lineWidth', 4); 
  title(j)
  hold off;
  
  subplot(1,2,2);
  plot(model.Likelihood.TimesF, Z, '+k', 'lineWidth', 2);
  pause(0.5);
  hold off;
  
