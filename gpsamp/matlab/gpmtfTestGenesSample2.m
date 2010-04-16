function [model PropDist samples accRates] = gpmtfTestGenesSample2(model, TFs, PropDist, trainOps)
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


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF; 
SizF = size(TimesF,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);
NumOfTFs = model.Likelihood.numTFs;

num_stored = floor(Iters/StoreEvery);
samples.TFindex = zeros(1, num_stored);
samples.kinetics = zeros(4, num_stored);
samples.W = zeros(NumOfTFs, num_stored);
samples.W0 = zeros(1, num_stored);
samples.LogL = zeros(1, num_stored);

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

onlyPumaVar = 1; 
if sum(model.Likelihood.noiseModel.active(2:3)) > 0
   onlyPumaVar = 0;
end


% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.W.constraint,'positive')
    posw = 1;
end


% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

n = SizF;
F = model.F;
TFindex = model.TFindex;


% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
if strcmp(model.constraints.replicas,'free')
   for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   oldLogLik(r,:) = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   oldLogLik(1,:) = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   % compute fast the additional likelihood when you have observations for the TF genes  
   %
   % the predicted genes are the same for the remaining coupled replicas
   for r=2:NumOfReplicas
      % compute fast the likelihood terms for the remaining replicas  
      oldLogLik(r,:) = remainRepsLikelihood(LikParams,  PredictedGenes, r, 1:NumOfGenes);
      % compute fast the additional likelihood when you have observations
      % for the TF genes  
   end
   %              
   %
end

% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);

% evaluation of the prior for the interaction bias
lnpriorW0 = ['ln',model.prior.W0.type,'pdf'];
TrspaceW0 = model.prior.W0.priorSpace; 
LikW0 = feval(TrspaceW0, LikParams.W0);
oldLogPriorW0 = feval(lnpriorW0, LikW0, model.prior.W0);

% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.W.type,'pdf'];
TrspaceW = model.prior.W.priorSpace; 
LikW = feval(TrspaceW, LikParams.W);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.W);

oldLogPriorNoiseWHITE = [];
if onlyPumaVar == 0 
    %
    % evaluation of the noise variance in the likelihood 
    if model.Likelihood.noiseModel.active(2) == 1
       % 
       % white GP noise model is used
       lnpriorNoiseWHITE = ['ln',model.prior.sigma2.type,'pdf'];
       TrspaceNoiseWHITE = model.prior.sigma2.priorSpace; 
       temp = feval(TrspaceNoiseWHITE, LikParams.noiseModel.sigma2);
       oldLogPriorNoiseWHITE = feval(lnpriorNoiseWHITE, temp, model.prior.sigma2);
       %
    else
       % rbf GP noise model is used
       lnpriorNoiseRBFsigma2f = ['ln',model.prior.sigma2f.type,'pdf'];
       TrspaceNoiseRBFsigma2f = model.prior.sigma2f.priorSpace; 
       temp = feval(TrspaceNoiseRBFsigma2f, LikParams.noiseModel.sigma2f);
       oldLogPriorNoiseRBFsigma2f = feval(lnpriorNoiseRBFsigma2f, temp, model.prior.sigma2f);
       
       lnpriorNoiseRBFlength = ['ln',model.prior.lengthScale.type,'pdf'];
       TrspaceNoiseRBFlength = model.prior.lengthScale.priorSpace; 
       temp = feval(TrspaceNoiseRBFlength, LikParams.noiseModel.lengthScale);
       oldLogPriorNoiseRBFlength = feval(lnpriorNoiseRBFlength, temp, model.prior.lengthScale);
       %
    end
    %
end

cnt = 0;

acceptF = 0; 
acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1, NumOfGenes);
if onlyPumaVar == 0 
   accRateNoiseM = zeros(1, NumOfGenes); 
end

numSamples = size(TFs,2);
%
for it = 1:(BurnInIters + Iters) 
    %

    % *
    % SAMPLE THE TFs FROM THE EMPIRICAL POSTERIOR DISTRIBUTION
    
    % choose a training sample 
    ch = round(rand*numSamples) + 1; 
    ch(ch>numSamples) = numSamples;
    
    LikParams1 = LikParams;
    % store the TF in the LikeParams to save computations 
    LikParams1.TF = TFs{ch};
    
    
    newLogProp = 0;  % forward Hastings Q(s_t+1 | s_t)
    oldLogProp = 0;  % backward Hastings Q(s_t| s_t+1) 
    
    newLogLik = zeros(NumOfReplicas, NumOfGenes);
    predgen = zeros(NumOfReplicas, model.Likelihood.sizTime);   
    if strcmp(model.constraints.replicas,'free') 
       % 
       for r=1:NumOfReplicas
       %   
           % perform an evaluation of the likelihood p(Genes | F) 
           newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
       %    
       end % num Replicas loop
       %
       
       % Metropolis-Hastings to accept-reject the proposal
       [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), newLogProp, oldLogProp);
           
       if (it > BurnInIters)
           acceptF = acceptF + accept; 
       end
    
       % update protein F
       if accept == 1
           %LikParams.TF(j,:,r) = TFs{ch}(j,:,r);
           LikParams.TF = TFs{ch};
           %TFindex(1:NumOfTFs,1:NumOfReplicas) = ch;
           TFindex = ch;
           oldLogLik = newLogLik;
        %   
       end
    %
    else
         % Replicas are coupled
         for j=1:NumOfTFs
             
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           gPerm = randperm(size(TFs,2));
           ch = gPerm(1); 
           
           newLogProp = 0;
           oldLogProp = 0; 
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           LikParams1.TF(j,:) = TFs{ch}(j,:);
             
           newLogLik = zeros(NumOfReplicas,NumOfGenes);
           % perform an evaluation of the likelihood p(Genes | F)      
           newLogLik(1,:) = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
              
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r,:) = remainRepsLikelihood(LikParams1,  predgen, r, 1:NumOfGenes);
           end                 
           % Metropolis-Hastings to accept-reject the proposal
           [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), newLogProp, oldLogProp);
          
           %
           if (it > BurnInIters)
               acceptF(j) = acceptF(j) + accept; 
           end
    
           % update protein F
           if accept == 1
               LikParams.TF(j,:) = TFs{ch}(j,:);
               TFindex(j) = ch;  
               oldLogLik = newLogLik;
            %   
           end
         end % num TFs loop
        %
    end % if end
    % END SAMPLE THE TFs FROM THE EMPIRICAL POSTERIOR DISTRIBUTION
    % *
    
    % *
    % SAMPLE KINETIC PARAMETERS 
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
        
        %%%%
        %KineticsNew = [eps 1 eps eps];
        %%%% 
       
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
    % END SAMPLE KINETIC PARAMETERS 
    % *
    
    % *
    % SAMPLE THE INTERACTION WEIGHTS
    for j=1:NumOfGenes
        % 
        %   
        newLogProp = 0;  % forward Hastings Q(w_t+1 | w_t)
        oldLogProp = 0;  % backward Hastings Q(w_t| w_t+1)  
        if posw == 1
            %
            % sample from a truncated Gaussian using rejection
            % sampling 
            while 1
                trW = randn(1,NumOfTFs).*sqrt(PropDist.W(j,1:NumOfTFs)) + LikParams.W(j,:);  
                if min(trW) > 0
                    break; 
                end
            end 
            Wnew(1:NumOfTFs) = trW;
            % sample also the bias which is allowed to be negative
            Wnew0 = randn.*sqrt(PropDist.W(j,NumOfTFs+1)) +  LikParams.W0(j);
            Wnew = [trW, Wnew0]; 
            
            trprior.sigma2 = PropDist.W(j,1:NumOfTFs); 
            trprior.mu = LikParams.W(j,:); 
        
            newLogProp = sum(lntruncNormalpdf(trW, trprior)); 
            
            trprior.mu = trW; 
            oldLogProp  = sum(lntruncNormalpdf(LikParams.W(j,:), trprior)); 
            
            %
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.constraints.W(j,:);
        
        if model.constraints.W0(j) == 0
           Wnew(NumOfTFs + 1) = 0;
        end
       
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
        
        
        % evaluation of the prior for the interaction bias 
        LikW0 = feval(TrspaceW0, Wnew(end));
        LogPriorWnew0 = feval(lnpriorW0, LikW0, model.prior.W0);
        % >>>  interaction weights
        LikW = feval(TrspaceW, Wnew(1:NumOfTFs));
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.W);
        
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorW(j,:),2) + oldLogPriorW0(j);
        newP = sum(newLogLik(:)) + sum(LogPriorWnew(:)) + LogPriorWnew0; 
         
        [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorW(j,:) = LogPriorWnew;
           oldLogPriorW0(j) = LogPriorWnew0;
        end
        %
        if (it > BurnInIters) 
           acceptW(j) = acceptW(j) + accept;
        end
        %
    end
    % END SAMPLE THE INTERACTION WEIGHTS
    % *
   
    % * 
    % SAMPLE THE NOISE MODEL IN THE LIKELIHOOD
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
                   newSigma2 = randn.*sqrt(PropDist.noiseModel(j,1)) + LikParams.noiseModel.sigma2(j);  
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
                   newLogLik(1) = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
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
         
               % compute the proposal distribution forwrd and backwards terms
               trprior.sigma2 = PropDist.noiseModel(j,1); 
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
           %
       elseif model.Likelihood.noiseModel.active(3) == 1 
           %
           % there exist only rbf GP noise  (per gene) that is added to the PUMA variance  
           
           % sample new  variance and lengthscale of the rbf kernel from truncated Gaussians
           for j =1:NumOfGenes
               %
               newLogProp = 0;  % forward Hastings Q(kern_t+1 | kern_t)
               oldLogProp = 0;  % backward Hastings Q(kern_t| kenr_t+1) 
               % sample from a truncated Gaussian proposal distribution
               while 1
                   newKern = randn(1,2).*sqrt(PropDist.noiseModel(j,2:3)) + [LikParams.noiseModel.sigma2f(j), LikParams.noiseModel.lengthScale(j)];  
                   if min(newKern) > 0
                       break; 
                   end
               end 
               
               LikParams1 = LikParams;
               LikParams1.noiseModel.sigma2f(j) = newKern(1); 
               LikParams1.noiseModel.lengthScale(j) = newKern(2); 
               
                        
               % Update the full covariance matrix for the new sample  
               if LikParams.noiseModel.active(1) == 1
               % add PUMA variances firstly  
                     for r=1:NumOfReplicas 
                         %
                         LikParams1.noiseModel.totalSigma{r}(:,:,j) = newKern(1)*exp(-(0.5/(newKern(2)^2)).*LikParams1.noiseModel.X2) ...
                                                                    + diag(LikParams1.noiseModel.pumaSigma2(j,:,r));
                         L = jitterChol( LikParams1.noiseModel.totalSigma{r}(:,:,j)  )'; % lower triangular 
                         LikParams1.noiseModel.InvL{r}(:,:,j) = L\eye(size(L,1)); 
                         LikParams1.noiseModel.LogDetSigma(j,r) = 2*sum(log(diag(L)));
                     end
                  %
               else % all replicas have the same covariance matrix
                  % 
                     LikParams1.noiseModel.totalSigma(:,:,j) = newKern(1)*exp(-(0.5/(newKern(2)^2)).*LikParams1.noiseModel.X2);
                     L =jitterChol( LikParams1.noiseModel.totalSigma(:,:,j)  )'; % lower triangular 
                     LikParams1.noiseModel.InvL(:,:,j) = L\eye(size(L,1)); 
                     LikParams1.noiseModel.LogDetSigma(j) = 2*sum(log(diag(L)));
                  %
               end
               
               % evaluate the new log likelihood 
               newLogLik = zeros(1,NumOfReplicas);
               if strcmp(model.constraints.replicas,'free')
                   for r=1:NumOfReplicas
                       % call the function only with j gene expressions  
                       newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
                   end
               else
                   % call the function only for the first replica
                   newLogLik(1) = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
                   % compute faster the remaining likelihood terms  
                   for r=2:NumOfReplicas
                       newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
                   end 
               end
              
               % evaluation of the prior for the new kernel variance
               LikSigma2f = feval(TrspaceNoiseRBFsigma2f, newKern(1));
               LogPriorNewRBFsigma2f = feval(lnpriorNoiseRBFsigma2f, LikSigma2f, model.prior.sigma2f);
               
               % evaluation of the prior for the new kernel variance
               LikLength = feval(TrspaceNoiseRBFlength, newKern(2));
               LogPriorNewRBFlength = feval(lnpriorNoiseRBFlength, LikLength, model.prior.lengthScale);
          
               % Metropolis-Hastings to accept-reject the proposal
               oldP = sum(oldLogLik(:,j),1) + oldLogPriorNoiseRBFsigma2f(j) + oldLogPriorNoiseRBFlength(j);
               newP = sum(newLogLik(:)) + LogPriorNewRBFsigma2f + LogPriorNewRBFlength;
               
               % compute the proposal distribution forward and backwards terms
               trprior.sigma2 = PropDist.noiseModel(j,2:3); 
               trprior.mu = [LikParams.noiseModel.sigma2f(j), LikParams.noiseModel.lengthScale(j)]; 
        
               newLogProp = sum(lntruncNormalpdf(newKern, trprior)); 
            
               trprior.mu = newKern; 
               oldLogProp  = sum(lntruncNormalpdf([LikParams.noiseModel.sigma2f(j), ...
                                                   LikParams.noiseModel.lengthScale(j)],  trprior)); 
               
               [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
               %
               if accept == 1
                  LikParams.noiseModel.sigma2f(j) = newKern(1);
                  LikParams.noiseModel.lengthScale(j) = newKern(2);
                  oldLogLik(:,j) = newLogLik(:); 
                  oldLogPriorNoiseRBFsigma2f(j) = LogPriorNewRBFsigma2f;
                  oldLogPriorNoiseRBFlength(j) = LogPriorNewRBFlength;
                  % Update kernel matrix and Cholesky decomposition
                  
                  % add PUMA variances firstly  
                  if LikParams.noiseModel.active(1) == 1
                  % 
                     for r=1:NumOfReplicas 
                         %
                         LikParams.noiseModel.totalSigma{r}(:,:,j) = LikParams1.noiseModel.totalSigma{r}(:,:,j); 
                         LikParams.noiseModel.InvL{r}(:,:,j) = LikParams1.noiseModel.InvL{r}(:,:,j); 
                         LikParams.noiseModel.LogDetSigma(j,r) = LikParams1.noiseModel.LogDetSigma(j,r);
                         %     
                     end
                  %
                  else
                  %               
                     LikParams.noiseModel.totalSigma(:,:,j) = LikParams1.noiseModel.totalSigma(:,:,j);
                     LikParams.noiseModel.InvL(:,:,j) = LikParams1.noiseModel.InvL(:,:,j);
                     LikParams.noiseModel.LogDetSigma(j) = LikParams1.noiseModel.LogDetSigma(j);
                     %
                  end
                  %
               end
               %
               if (it > BurnInIters) 
                  accRateNoiseM(j) = accRateNoiseM(j) + accept;
               end
               %
               %
           end % for loop
       else
           %
           % do nothing (only PUMA variances are used whihc are fixed)
           %
       end % if-then-else statement 
       % 
       %
    end % if statement
    % END SAMPLE THE NOISE MODEL IN THE LIKELIHOOD
    % *
    
    % *
    % KEEP SAMPLES AFTER BURN IN
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.TFindex(cnt) = TFindex;
        %samples.predGenes(:,:,cnt) = PredictedGenes;
        samples.kinetics(:,cnt) = LikParams.kinetics;
        samples.W(:,cnt) = LikParams.W;
        samples.W0(:,cnt) = LikParams.W0;
        if onlyPumaVar == 0
           if model.Likelihood.noiseModel.active(2) == 1 
              samples.sigma2(cnt) = LikParams.noiseModel.sigma2;
           end
           if model.Likelihood.noiseModel.active(3) == 1 
              samples.sigma2f(cnt) = LikParams.noiseModel.sigma2f;
              samples.lengthScale(cnt) = LikParams.noiseModel.lengthScale;
           end
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        samples.LogPrior(cnt) = sum(oldLogPriorW(:)) + sum(oldLogPriorW0(:)) + sum(oldLogPriorKin(:)); 
        if ~isempty(oldLogPriorNoiseWHITE)
           samples.LogPrior(cnt) = samples.LogPrior(cnt) + sum(oldLogPriorNoiseWHITE(:));
        end
        %
        %
    end
    % END KEEP SAMPLES AFTER BURN IN
    % *
    
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.kinetics = LikParams.kinetics;
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;
if onlyPumaVar == 0
   if model.Likelihood.noiseModel.active(2) == 1 
      model.Likelihood.noiseModel.sigma2 = LikParams.noiseModel.sigma2;
   end
   if model.Likelihood.noiseModel.active(3) == 1 
      model.Likelihood.noiseModel.sigma2f = LikParams.noiseModel.sigma2f;
      model.Likelihood.noiseModel.lengthScale = LikParams.noiseModel.lengthScale;
      model.Likelihood.noiseModel.totalSigma = LikParams.noiseModel.totalSigma;
      model.Likelihood.noiseModel.InvL = LikParams.noiseModel.InvL; 
      model.Likelihood.noiseModel.LogDetSigma = LikParams.noiseModel.LogDetSigma;
   end
end

model.TFindex = TFindex; 

accRates.F = (acceptF/Iters)*100;
accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
if onlyPumaVar == 0 & ~(model.Likelihood.noiseModel.active(1) == 0 &  model.Likelihood.noiseModel.active(2) == 1) 
    accRates.noiseM = (accRateNoiseM/Iters)*100;
else
    accRates.noiseM = 25;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikval = remainRepsLikelihood(LikParams,  PredictedGenes, r, Gindex)
%
    Gns = LikParams.Genes(:,:,r);
    loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,r)),2) ...
                - 0.5*sum(((Gns(Gindex,:) - PredictedGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,r),2);
    loglikval = loglikval';
% 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikvalTF = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, TFindex)
%    
  
  GnsTF = LikParams.GenesTF(:,:,r);
  loglikvalTF = - 0.5*sum(log(2*pi*LikParams.sigmasTF(TFindex,:,r)),2) ...
                       - 0.5*sum(((GnsTF(TFindex,:) - PredictedGenesTF(:,LikParams.comIndsTF)).^2)./LikParams.sigmasTF(TFindex,:,r),2);
  loglikvalTF = loglikvalTF';
%  
