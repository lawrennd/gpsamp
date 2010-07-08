function [model PropDist samples accRates] = gpmtfSample(model, PropDist, trainOps)
% Description: Draw a set of samples from the Bayesian differential
%              equation model
%
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- trainOps: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   parameters values of the final MCMC iteration
%                   parameters as well as the priors
%         -- PropDist: as above. PropDist can be updated (compared to the input one) 
%                     due to the update of the kernel parameters that
%                     influence the proposal 
%         -- samples: the structure that contrains the samples 
%         -- accRates: acceptance rates 
%

%ppi = 0.5;
BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF; 
SizF = size(TimesF,2);
SizKin = size(model.Likelihood.kinetics,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);


NumOfTFs = model.Likelihood.numTFs;
istart = ones(NumOfTFs,1);
for j=1:NumOfTFs
if model.constraints.Ft0(j)==0
   % do not sample the first control point so that the function 
   % will be fixed at the time t=0
   istart(j)=2;
end
end

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

if isfield(model.Likelihood,'GenesTF')
    SizTFKin = size(model.Likelihood.kineticsTF,2);
end


% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.W.constraint,'positive')
    posw = 1;
end

onlyPumaVar = 1; 
if model.Likelihood.noiseModel.active(2) > 0
   onlyPumaVar = 0;
end

% if there are puma variances for the TF genes, then use only those 
onlyPumaVarTF = 1;
if model.Likelihood.noiseModel.active(1) == 0 & model.Likelihood.noiseModel.active(2) > 0
     onlyPumaVarTF = 0;
end
   

LikParams = model.Likelihood;


% the latent function values are also special parameters that appear in both 
% the likelihood and the GP prior
F = model.F; 
%F = model.groundtr.F;

n = SizF;

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
if strcmp(model.constraints.replicas,'free')
   for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   oldLogLik(r,:) = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
  
   % additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      oldLogLikTF(r,:) = gpmtfLogLikelihoodGeneTF(model.Likelihood, F(:,:,r), r, 1:NumOfTFs);
   end
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   oldLogLik(1,:) = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   % compute fast the additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      oldLogLikTF(1,:) = gpmtfLogLikelihoodGeneTF(model.Likelihood, F, 1, 1:NumOfTFs);
   end
   %
   % the predicted genes are the same for the remaining coupled replicas
   for r=2:NumOfReplicas
      % compute fast the likelihood terms for the remaining replicas  
      oldLogLik(r,:) = remainRepsLikelihood(LikParams,  PredictedGenes, r, 1:NumOfGenes);
      % compute fast the additional likelihood when you have observations for the TF genes  
      if isfield(model.Likelihood,'GenesTF')
         oldLogLikTF(r,:) = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, 1:NumOfTFs);
      end
   end
   %              
   %
end
%
% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);
if isfield(model.Likelihood,'GenesTF')
  LikkinTF = feval(TrspaceKin, LikParams.kineticsTF); 
  oldLogPriorKinTF = feval(lnpriorKin, LikkinTF, model.prior.kinetics);
end      
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.W.type,'pdf'];
TrspaceW = model.prior.W.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.W);
% log prior of the lengthscale lengthscale 
lnpriorLengSc = ['ln',model.prior.lengthScale.type,'pdf'];
for j=1:NumOfTFs
oldLogPriorLengSc(j) = feval(lnpriorLengSc, model.GP{j}.lengthScale, model.prior.lengthScale);
end

%%% IF there are observations for the TF-genes then you sample also
%%% the variance of the rbf kernel 
if isfield(model.Likelihood,'GenesTF')
   lnpriorNoiseRBFsigma2f = ['ln',model.prior.sigma2f.type,'pdf'];
   TrspaceNoiseRBFsigma2f = model.prior.sigma2f.priorSpace; 
   for j=1:NumOfTFs
       temp = feval(TrspaceNoiseRBFsigma2f, model.GP{j}.sigma2f);
       oldLogPriorNoiseRBFsigma2f(j) = feval(lnpriorNoiseRBFsigma2f, temp, model.prior.sigma2f);
   end
end

if onlyPumaVar == 0 
  %
  % evaluation of the noise variance prior in the likelihood 
  if model.Likelihood.noiseModel.active(2) == 1
  % 
     % white GP noise model is used
     lnpriorNoiseWHITE = ['ln',model.prior.sigma2.type,'pdf'];
     TrspaceNoiseWHITE = model.prior.sigma2.priorSpace; 
     temp = feval(TrspaceNoiseWHITE, LikParams.noiseModel.sigma2);
     oldLogPriorNoiseWHITE = feval(lnpriorNoiseWHITE, temp, model.prior.sigma2);
     
     if isfield(model.Likelihood,'GenesTF')
        lnpriorNoiseWHITETF = ['ln',model.prior.sigma2.type,'pdf'];
        TrspaceNoiseWHITETF = model.prior.sigma2.priorSpace; 
        temp = feval(TrspaceNoiseWHITETF, LikParams.noiseModel.sigma2_TF);        
        oldLogPriorNoiseWHITETF = feval(lnpriorNoiseWHITE, temp, model.prior.sigma2); 
     end
  %
  end
  %
end


cnt = 0;
acceptF = zeros(NumOfTFs,NumOfReplicas);
acceptKin = zeros(1,NumOfGenes);
acceptTFKin = zeros(1,NumOfTFs);
acceptW = zeros(1,NumOfGenes); 
acceptLengSc = zeros(1, NumOfTFs);  
accRateNoiseM = zeros(1, NumOfGenes); 
accRateNoiseMTF = zeros(1, NumOfTFs); 

n = size(model.F,2);

%
for it = 1:(BurnInIters + Iters) 
    % 
    
    % *
    % SAMPLE THE TFs USING CONTROL VARIABLES
    %F = model.groundtr.F;
    if 1 
    for j=1:NumOfTFs   
    for r=1:NumOfReplicas
         %
         %
         Z(j,:,r) = F(j,:,r) + sqrt(model.auxLikVar(j,:,r)).*randn(1,n);
          
         % STEP 2: M-H step
         if model.constraints.Ft0(j)==0  
           cmu = model.mu(:,j)' + (Z(j,:,r) - model.mu(:,j)')*model.invKsigmaK{r}(:,:,j);
         else
           cmu = Z(j,:,r) *model.invKsigmaK{r}(:,:,j); 
         end
         Fnew = gaussianFastSample(1, cmu, model.auxPostL{r}(:,:,j));
         
         FFnew = F(:,:,r);
         FFnew(j,:) = Fnew;
         
         if ~isfield(model.Likelihood,'GenesTF')
             % perform an evaluation of the likelihood p(Genes | F) 
             newLogLik = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);
       
             % Metropolis-Hastings to accept-reject the proposal
             [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(r,:),2), 0, 0);
         else
             % perform an evaluation of the likelihood p(Genes | F) 
             newLogLik = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);        
             newLogLikTF = gpmtfLogLikelihoodGeneTF(LikParams, FFnew, r, j); 
           
             % Metropolis-Hastings to accept-reject the proposal
             newL = sum(newLogLik(:)) + newLogLikTF;
             oldL = sum(oldLogLik(r,:),2) + oldLogLikTF(r,j);
             [accept, uprob] = metropolisHastings(newL, oldL, 0, 0);
         end
         
         % visualization 
         if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0)
         %    
             visualization(model, F(j,:,r), Fnew, Z(j,:,r), j);
         %    
         end
         %
         if (it > BurnInIters)
             acceptF(j,r) = acceptF(j,r) + accept; 
         end
         % update protein F
         if accept == 1
            F(j,:,r) = Fnew;
            oldLogLik(r,:) = newLogLik;
            if isfield(model.Likelihood,'GenesTF')      
                 oldLogLikTF(r,j) = newLogLikTF;
            end
            %   
         end
    end % num Replicas loop
    end % num TFs loop
    end % zero-one if
    % END SAMPLE THE TFs USING IMAGINARY DATA
    % *
    
    
    % *
    % SAMPLE KINETIC PARAMETERS FOR THE TF TRANSLATION MODEL
    if 1
    if isfield(model.Likelihood,'GenesTF')
    for j=1:NumOfTFs
        TFKineticsNew = randn(1,SizTFKin).*sqrt(PropDist.TFkin(j,:)) + log(LikParams.kineticsTF(j,:));
        TFKineticsNew(TFKineticsNew<-10) =-10; 
        TFKineticsNew(TFKineticsNew>10) = 10;
        TFKineticsNew = exp(TFKineticsNew); 
        %
        if model.constraints.geneTFsensitivity(j) == 0
        TFKineticsNew(2) = 1; 
        end
        
        LikParams1 = LikParams;
        LikParams1.kineticsTF(j,:)=TFKineticsNew; 
        newLogLik = []; 
        if strcmp(model.constraints.replicas,'free')   
           for r=1:NumOfReplicas
               % perform an evaluation of the likelihood p(GENES | TFs) 
               newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
               % 
           end
        else
           %
               % perform an evaluation of the likelihood p(GENES | TFs)
               % only for the first replica
               newLogLik(1,:) = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
               % computed faster the remaining likelihood terms  
               for r=2:NumOfReplicas
                   newLogLik(r,:) = remainRepsLikelihood(LikParams1, predgen, r, 1:NumOfGenes);
               end 
           %
        end
        %        
        Likkin = feval(TrspaceKin, TFKineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics);
        
        % Metropolis-Hastings to accept-reject the proposal
        newP = sum(newLogLik(:)) + sum(LogPriorKinNew(:));
        oldP = sum(oldLogLik(:)) + sum(oldLogPriorKinTF(j,:),2);
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        %[accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.kineticsTF(j,:) = TFKineticsNew;
           oldLogLik = newLogLik; 
           oldLogPriorKinTF(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptTFKin(j) = acceptTFKin(j) + accept;
        end
        %
    end
    end
    end
    % *
    % END SAMPLE KINETIC PARAMETERS FOR THE TF TRANSLATION MODEL
    
    % *
    % SAMPLE KINETIC PARAMETERS FOR EACH GENE ODE
    if 1 
    for j=1:NumOfGenes
        KineticsNew = randn(1,SizKin).*sqrt(PropDist.kin(j,:)) + log(LikParams.kinetics(j,:));
        KineticsNew(KineticsNew<-10) =-10; 
        KineticsNew(KineticsNew>10) = 10;
        KineticsNew = exp(KineticsNew); 
        %
        % set the initial condition to be zero
        % if you know that it should be zero
        if model.constraints.InitialConds(j) == 0
        KineticsNew(4) = KineticsNew(1)/KineticsNew(2); 
        end
    
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           %
           for r=1:NumOfReplicas
              % call the function only with j gene expressions  
              newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
            % call the function only for the first replica
            [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
            % computed faster the remaining likelihood terms  
            for r=2:NumOfReplicas
                newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
            end 
           %
        end
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorKin(j,:),2);
        newP = sum(newLogLik(:))+ sum(LogPriorKinNew(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        %[accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.kinetics(j,:) = KineticsNew;
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorKin(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptKin(j) = acceptKin(j) + accept;
        end
        %
    end
    %
    end % if 0
    % *
    % END SAMPLE KINETIC PARAMETERS FOR EACH GENE ODE
   
    % *
    % SAMPLE INTERACTION WEIGHTS FOR EACH GENE
    if 1  
    for j=1:NumOfGenes
        % 
        %  
        if posw == 1
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + log([LikParams.W(j,:), LikParams.W0(j)]);     
            Wnew = exp(Wnew);
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.constraints.W(j,:);
       
        LikParams1 = LikParams;
        LikParams1.W(j,:) = Wnew(1:NumOfTFs);
        LikParams1.W0(j)=Wnew(end);
      
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           for r=1:NumOfReplicas
               % call the function only with j gene expressions  
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
           % call the function only for the first replica
           [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
           end 
        end
        
        LikW = feval(TrspaceW, Wnew+eps);
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.W);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorW(j,:),2);
        newP = sum(newLogLik(:)) + sum(LogPriorWnew(:)); 
        %
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorW(j,:) = LogPriorWnew;
        end
        %
        if (it > BurnInIters) 
           acceptW(j) = acceptW(j) + accept;
        end
        %
    end
    end % if 0
    % *
    % END SAMPLE INTERACTION WEIGHTS FOR EACH GENE
    
    % * 
    % SAMPLE THE NOISE MODEL IN THE LIKELIHOOD of THE GENES 
    if onlyPumaVar == 0 
       %  only white noise (no PUMA variances)
       if model.Likelihood.noiseModel.active(1) == 0 &  model.Likelihood.noiseModel.active(2) == 1 
           %
           sumSquerrors1 = oldLogLik + 0.5*SizG*log(2*pi*repmat(LikParams.noiseModel.sigma2, NumOfReplicas, 1));
           sumSquerrors1 = -2*repmat(LikParams.noiseModel.sigma2, NumOfReplicas, 1).*sumSquerrors1;
           sumSquerrors = sum(sumSquerrors1,1);
                
           anew = model.prior.sigma2.a + 0.5*NumOfReplicas*SizG;
           bnew = model.prior.sigma2.b + 0.5*sumSquerrors;
           if exist('xgamrnd')
           newinvSigma2 = xgamrnd(anew, 1./bnew);
           else
           newinvSigma2 = gamrnd(anew, 1./bnew);
           end
           newSigma2 = 1./newinvSigma2;
         
           LikParams.noiseModel.sigma2 = newSigma2;
           
           okk = repmat(LikParams.noiseModel.sigma2, NumOfReplicas, 1); 
       
           oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
           % always accepted
           if (it > BurnInIters) 
                accRateNoiseM = accRateNoiseM + ones(1, NumOfGenes);
           end
           %
       elseif model.Likelihood.noiseModel.active(1) == 1 & model.Likelihood.noiseModel.active(2) == 1 
           %
           % there exist only white noise variance (per gene) that is added to the PUMA variance
           %
           % sample new variances from truncated Gaussians
           for j =1:NumOfGenes
               %
               newLogProp = 0;  % forward Hastings Q(sigma2_t+1 | sigma2_t)
               oldLogProp = 0;  % backward Hastings Q(sigma2_t| sigma2_t+1) 
               % sample from a truncated Gaussian proposal distribution
               while 1
                   newSigma2 = randn.*sqrt(PropDist.noiseModel(j)) + LikParams.noiseModel.sigma2(j);  
                   if newSigma2 > 0
                       break; 
                   end
               end
               %
                   
               LikParams1 = LikParams; 
               LikParams1.noiseModel.sigma2(j) = newSigma2;
 
               % evaluate the new log likelihood 
               newLogLik = zeros(1,NumOfReplicas);
               if strcmp(model.constraints.replicas,'free')
                   for r=1:NumOfReplicas
                       % call the function only with j gene expressions  
                       newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
                   end
               else
                   % call the function only for the first replica
                   [newLogLik(1), predgen]  = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
                   % compute faster the remaining likelihood terms  
                   for r=2:NumOfReplicas
                       newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
                   end 
               end
              
               
               % evaluation of the prior to the new white variance
               LikSigma2 = feval(TrspaceNoiseWHITE, newSigma2);
               LogPriorNewWHITE = feval(lnpriorNoiseWHITE, LikSigma2, model.prior.sigma2);
       
               % Metropolis-Hastings to accept-reject the proposal
               oldP = sum(oldLogLik(:,j),1) + oldLogPriorNoiseWHITE(j);
               newP = sum(newLogLik(:)) + LogPriorNewWHITE; 
         
               % compute the proposal distribution forward and backwards terms
               trprior.sigma2 = PropDist.noiseModel(j); 
               trprior.mu = LikParams.noiseModel.sigma2(j); 
        
               newLogProp = sum(lntruncNormalpdf(newSigma2, trprior)); 
            
               trprior.mu = newSigma2; 
               oldLogProp  = sum(lntruncNormalpdf(LikParams.noiseModel.sigma2(j), trprior)); 
               
               [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
               %
               if accept == 1
                  LikParams.noiseModel.sigma2(j) = newSigma2;
                  oldLogLik(:,j) = newLogLik(:); 
                  oldLogPriorNoiseWHITE(j) = LogPriorNewWHITE;
               end
               %
               if (it > BurnInIters) 
                  accRateNoiseM(j) = accRateNoiseM(j) + accept;
               end
               %
           end % for loop 
           %
       else
           %
           % do nothing (only PUMA variances are used whihc are fixed)
           %
       end % if-then-else statement 
       % 
       %
    end % if statement
    % END SAMPLE THE NOISE MODEL IN THE LIKELIHOOD of THE GENES 
    % *
   
    % * 
    % SAMPLE THE NOISE MODEL IN THE LIKELIHOOD OF THE !!TF!! GENES 
    if onlyPumaVarTF == 0 
       %  only white noise (no PUMA variances)
       if model.Likelihood.noiseModel.active(1) == 0 &  model.Likelihood.noiseModel.active(2) == 1 
           %
           sumSquerrors1 = oldLogLikTF + 0.5*SizG*log(2*pi*repmat(LikParams.noiseModel.sigma2_TF, NumOfReplicas, 1));
           sumSquerrors1 = -2*repmat(LikParams.noiseModel.sigma2_TF, NumOfReplicas,1).*sumSquerrors1; 
           sumSquerrors = sum(sumSquerrors1,1);
                
           anew = model.prior.sigma2.a + 0.5*NumOfReplicas*SizG;
           bnew = model.prior.sigma2.b + 0.5*sumSquerrors;
           if exist('xgamrnd')
           newinvSigma2 = xgamrnd(anew, 1./bnew);
           else
           newinvSigma2 = gamrnd(anew, 1./bnew);
           end
           newSigma2 = 1./newinvSigma2;
           
           LikParams.noiseModel.sigma2_TF = newSigma2;
           
           okk = repmat(LikParams.noiseModel.sigma2_TF, NumOfReplicas, 1); 
       
           oldLogLikTF = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
           % always accepted
           if (it > BurnInIters) 
                 accRateNoiseMTF = accRateNoiseMTF + ones(1, NumOfTFs);
           end
           %
       elseif model.Likelihood.noiseModel.active(1) == 1 & model.Likelihood.noiseModel.active(2) == 1 
           %
           % there exist only white noise variance (per gene) that is added to the PUMA variance
           %
           % sample new variances from truncated Gaussians
           for j =1:NumOfTFs
               %
               newLogProp = 0;  % forward Hastings Q(sigma2_t+1 | sigma2_t)
               oldLogProp = 0;  % backward Hastings Q(sigma2_t| sigma2_t+1) 
               % sample from a truncated Gaussian proposal distribution
               while 1
                   newSigma2 = randn.*sqrt(PropDist.noiseModelTF(j)) + LikParams.noiseModel.sigma2_TF(j);  
                   if newSigma2 > 0
                       break; 
                   end
               end 
               %
               
               LikParams1 = LikParams; 
               LikParams1.noiseModel.sigma2_TF(j) = newSigma2;
 
               % evaluate the new log likelihood 
               newLogLikTF = zeros(1,NumOfReplicas);
               if strcmp(model.constraints.replicas,'free')
                   for r=1:NumOfReplicas
                       % call the function only with j gene expressions  
                       newLogLikTF(r) = gpmtfLogLikelihoodGeneTF(LikParams1, F(:,:,r), r, j);
                   end
               else
                   % call the function only for the first replica
                   [newLogLikTF(1), predgen]  = gpmtfLogLikelihoodGeneTF(LikParams1, F, 1, j);
                   % compute faster the remaining likelihood terms  
                   for r=2:NumOfReplicas
                       newLogLikTF(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
                   end 
               end
            
               % evaluation of the prior to the new white variance
               LikSigma2 = feval(TrspaceNoiseWHITETF, newSigma2);
               LogPriorNewWHITETF = feval(lnpriorNoiseWHITETF, LikSigma2, model.prior.sigma2);
       
               % Metropolis-Hastings to accept-reject the proposal
               oldP = sum(oldLogLikTF(:,j),1) + oldLogPriorNoiseWHITETF(j);
               newP = sum(newLogLikTF(:)) + LogPriorNewWHITETF; 
         
               % compute the proposal distribution forward and backwards terms
               trprior.sigma2 = PropDist.noiseModelTF(j); 
               trprior.mu = LikParams.noiseModel.sigma2_TF(j); 
        
               newLogProp = sum(lntruncNormalpdf(newSigma2, trprior)); 
            
               trprior.mu = newSigma2; 
               oldLogProp  = sum(lntruncNormalpdf(LikParams.noiseModel.sigma2_TF(j), trprior)); 
               
               [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
               %
               if accept == 1
                  LikParams.noiseModel.sigma2_TF(j) = newSigma2;
                  oldLogLikTF(:,j) = newLogLikTF(:); 
                  oldLogPriorNoiseWHITETF(j) = LogPriorNewWHITETF;
               end
               %
               if (it > BurnInIters) 
                  accRateNoiseMTF(j) = accRateNoiseMTF(j) + accept;
               end
               %
           end % for loop 
           %
       else
           %
           % do nothing (only PUMA variances are used whihc are fixed)
           %
       end % if-then-else statement 
       % 
       %
    end % if statement
    % END SAMPLE THE NOISE MODEL IN THE LIKELIHOOD OF THE !!TF!! GENES 
    % *
    
    
    % * 
    % SAMPLE THE DELAYS IN THE GENE ODES 
    if 1 
    if model.Likelihood.tauMax < 0
    for j=1:NumOfGenes
        LikParams1 = LikParams;
        if LikParams1.Tausindex(j) == 1
            % propose to decrease tau
            LikParams1.Tausindex(j) = 2;
            LikParams1.Taus(j) = LikParams1.Taus(j) - LikParams1.step;
            logbias = log(0.5);
            %
        elseif LikParams1.Tausindex(j) == LikParams1.startTime
            % propose to increase tau 
            LikParams1.Tausindex(j) = LikParams1.startTime-1;
            LikParams1.Taus(j) = LikParams1.Taus(j) + LikParams1.step;
            logbias = -log(0.5);
            %
        else
            %
            % propose to decrease or increase with probability 0.5
            ud = round(rand); 
            ud(ud==0)=-1;
            logbias = 0;
            LikParams1.Tausindex(j) = LikParams1.Tausindex(j)-ud;
            LikParams1.Taus(j) = LikParams1.Taus(j) + ud*LikParams1.step;
        end
        
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           for r=1:NumOfReplicas
               % call the function only with j gene expressions  
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
           % call the function only for the first replica
           [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
           end  
        end
        
        oldP = log(model.prior.delays.prob(LikParams.Tausindex(j))) + sum(oldLogLik(:,j),1) - logbias;
        newP = log(model.prior.delays.prob(LikParams1.Tausindex(j))) + sum(newLogLik(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.Tausindex(j) = LikParams1.Tausindex(j);
           LikParams.Taus(j) = LikParams1.Taus(j);
           oldLogLik(:,j) = newLogLik(:);
        end
        %
        %
    end
    end
    end
    % END SAMPLE THE DELAYS IN THE GENE ODES 
    % *
    
    % * 
    % SAMPLE THE LENGTHSCALES OF THE GP FUNCTIONS 
    if 1 
    for j=1:NumOfTFs
        %
        % to samples the lengthscale you need to evaluate 
        newLogProp = 0;  % forward Hastings Q(kern_t+1 | kern_t)
        oldLogProp = 0;  % backward Hastings Q(kern_t| kenr_t+1) 
        % sample from a truncated Gaussian proposal distribution
        if isfield(model.Likelihood,'GenesTF')
            while 1
              newKern = randn(1, 2).*sqrt(PropDist.Kern(j, 1:2)) + [model.GP{j}.sigma2f, model.GP{j}.lengthScale];  
              if min(newKern) > 0
              %    
                 if strcmp(model.prior.lengthScale.type, 'uniform')
                    if newKern(2) >= model.prior.lengthScale.constraint(1)  & newKern(2) <= model.prior.lengthScale.constraint(2) 
                       break;
                    end
                 else
                    break;  
                 end
              %   
              end
            end 
        else    
            while 1
            %    
              tmpK = randn.*sqrt(PropDist.Kern(j, 2)) + model.GP{j}.lengthScale;  
              if min(tmpK) > 0
              %    
                 if strcmp(model.prior.lengthScale.type, 'uniform')
                    if tmpK >= model.prior.lengthScale.constraint(1)  & tmpK <= model.prior.lengthScale.constraint(2) 
                       break;
                    end
                 else
                    break;  
                 end
              %   
              end
            % 
            end 
            newKern = [model.GP{j}.sigma2f, tmpK];
        end        

        newK = newKern(1)*exp(-(0.5/newKern(2))*model.GP{j}.X2) + model.GP{j}.sigma2*eye(size(model.GP{j}.X2,1)); 
        
        
        % condition on the first point being zero if needed
        if model.constraints.Ft0(j)==0  
           newKfu = newKern(1)*exp(-(0.5/newKern(2))*(model.Likelihood.TimesF(:) - model.Likelihood.TimesF(1)).^2); 
           newKuu = newKern(1) + model.GP{j}.sigma2;
     
           % update prior
           newK = newK - (newKfu*newKfu')/newKuu;
           newmu = (newKfu*model.u(j,1))/newKuu;
        end
          
        % compute the Cholesky decomposition of the new K
        [newL, er]=jitterChol(newK);
        newL = newL';
        % evaluate the new log GP prior value 
        invnewL = newL\eye(SizF);
        newLogDetK = 2*sum(log(diag(newL)));
        
        newlogGP = - 0.5*NumOfReplicas*newLogDetK;
        oldlogGP = - 0.5*NumOfReplicas*model.LogDetK(j);
        if model.constraints.Ft0(j)==0  
           for r=1:NumOfReplicas     
              temp = invnewL*(F(j,:,r)' - newmu); 
              newlogGP = newlogGP - 0.5*temp'*temp;
              temp = model.invL(:,:,j)*(F(j,:,r)' - model.mu(:,j)); 
              oldlogGP = oldlogGP - 0.5*temp'*temp;
           end
        else
           for r=1:NumOfReplicas     
              temp = invnewL*(F(j,:,r)'); 
              newlogGP = newlogGP - 0.5*temp'*temp;
              temp = model.invL(:,:,j)*(F(j,:,r)'); 
              oldlogGP = oldlogGP - 0.5*temp'*temp;
           end 
        end
        
         
        LogPriorNewRBFsigma2fnew = feval(lnpriorNoiseRBFsigma2f, newKern(1), model.prior.sigma2f);
        LogPriorLengScnew = feval(lnpriorLengSc, newKern(2), model.prior.lengthScale);
        
        % Metropolis-Hastings to accept-reject the proposal
        oldlogGP = oldlogGP + oldLogPriorNoiseRBFsigma2f(j) + oldLogPriorLengSc(j);
        newlogGP = newlogGP + LogPriorNewRBFsigma2fnew + LogPriorLengScnew; 
       
        % compute the proposal distribution forward and backwards terms
        trprior.sigma2 = PropDist.Kern(j, 1:2); 
        trprior.mu = [model.GP{j}.sigma2f, model.GP{j}.lengthScale]; 
        
        newLogProp = sum(lntruncNormalpdf(newKern, trprior)); 
            
        trprior.mu = newKern; 
        oldLogProp  = sum(lntruncNormalpdf([model.GP{j}.sigma2f, model.GP{j}.lengthScale], trprior)); 
                                               
        [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, newLogProp, oldLogProp);
        
        %%%%%%%%%%%%%%%%  start accept/update proposal for the lengthscale %%%%%%%%%%%% 
        if accept == 1
           model.GP{j}.sigma2f = newKern(1); 
           model.GP{j}.lengthScale = newKern(2); 
           oldLogPriorNoiseRBFsigma2f(j) = LogPriorNewRBFsigma2fnew;
           oldLogPriorLengSc(j) = LogPriorLengScnew;
           model.K(:,:,j) = newK;
           
           if model.constraints.Ft0(j)==0  
              model.Kfu(:,j) = newKfu;
              model.Kuu(j) = newKuu;
              model.mu(:,j) = newmu;
           end
           
           model.LogDetK(j) = newLogDetK;
           model.invL(:,:,j) = invnewL;
           %L = chol(model.K + diag(model.auxLikVar));    
           %tmp = L\eye(n);
           %invK = tmp*tmp';
           %model.invKsigmaK = invK*model.K;
           for r=1:NumOfReplicas
              model.invKsigmaK{r}(:,:,j) = (model.K(:,:,j) + diag(model.auxLikVar(j,:,r)))\model.K(:,:,j); 
              model.auxPostL{r}(:,:,j) = jitterChol( model.K(:,:,j) - model.K(:,:,j)*model.invKsigmaK{r}(:,:,j) );
           end
        end
        % 
        if (it > BurnInIters) 
           acceptLengSc(j) = acceptLengSc(j) + accept;
        end
        %%%%%%%%%%%%%%%%%%%%%%% end accept/update proposal for the lengthscale %%%%%%%%%%%%%%%%
        %
    end
    end 
    % END SAMPLE THE LENGTHSCALES OF THE GP FUNCTIONS
    % *
    
    %
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.F{cnt} = F;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.W(:,:,cnt) = LikParams.W;
        samples.W0(:,cnt) = LikParams.W0;
        %
        if model.Likelihood.tauMax < 0
           samples.Taus(:,cnt) = LikParams.Taus(:);
           samples.Tausindex(:,cnt) = LikParams.Tausindex(:);
        end
        %
        if isfield(model.Likelihood,'GenesTF')
            samples.kineticsTF(:,:,cnt) = LikParams.kineticsTF;
            samples.LogLTF(cnt) = sum(oldLogLikTF(:));
            %samples.predGenesTF{cnt} = PredictedGenesTF;
            if onlyPumaVar == 0
                samples.sigma2_TF(:,cnt) = LikParams.noiseModel.sigma2_TF';
            end
        end
        if onlyPumaVar == 0
           samples.sigma2(:,cnt) = LikParams.noiseModel.sigma2';
        end
        %if netLearn == 1
        %    samples.NetX(:,:,cnt) = LikParams.Net_X;
        %end
        for jin=1:NumOfTFs
            samples.lengthScale(jin,cnt) = model.GP{jin}.lengthScale;
            samples.sigma2f(jin,cnt) = model.GP{jin}.sigma2f;
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        %
    end
    %
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.kinetics = LikParams.kinetics;
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;
model.Likelihood.Tausindex = LikParams.Tausindex;
model.Likelihood.Taus = LikParams.Taus;
if isfield(model.Likelihood,'GenesTF')
    model.Likelihood.kineticsTF = LikParams.kineticsTF;
end

if onlyPumaVar == 0
    model.Likelihood.noiseModel.sigma2 = LikParams.noiseModel.sigma2;
    if isfield(model.Likelihood,'GenesTF')
        model.Likelihood.noiseModel.sigma2_TF = LikParams.noiseModel.sigma2_TF;
    end
end

model.F = F;
accRates.F = (acceptF/Iters)*100; 
accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
accRates.Kern = (acceptLengSc/Iters)*100;
if isfield(model.Likelihood,'GenesTF')
   accRates.TFKin = (acceptTFKin/Iters)*100;
end

accRates.noiseM = (accRateNoiseM/Iters)*100;
if isfield(model.Likelihood,'GenesTF') & (onlyPumaVarTF == 0) 
    accRates.noiseMTF = (accRateNoiseMTF/Iters)*100;
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
  
