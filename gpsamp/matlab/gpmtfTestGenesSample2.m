function [model PropDist samples accRates] = gpmtfTestGenesSample2(model, TFs, simMat, PropDist, trainOps)
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
samples.predGenes = zeros(NumOfReplicas, SizF, num_stored);
samples.kinetics = zeros(4, num_stored);
samples.Weights = zeros(NumOfTFs, num_stored);
samples.Weights0 = zeros(1, num_stored);
samples.LogL = zeros(1, num_stored);

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

% check if the observation noise is known/fixed
fixsigma2 = 0;
if strcmp(model.constraints.sigmas,'fixed')
    fixsigma2 = 1;
end

% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
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
   [oldLogLik(r,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
   PredictedGenes(r,:) = predgen;
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   [oldLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   PredictedGenes = predgen;
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
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);

% evaluation of the prior for the interaction bias
lnpriorW0 = ['ln',model.prior.weight0.type,'pdf'];
TrspaceW0 = model.prior.weight0.priorSpace; 
LikW0 = feval(TrspaceW0, LikParams.W0);
oldLogPriorW0 = feval(lnpriorW0, LikW0, model.prior.weight0);

% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, LikParams.W);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights);

cnt = 0;

%if strcmp(model.constraints.replicas,'free')
%   acceptF = zeros(NumOfTFs,NumOfReplicas);
%else
%   acceptF = zeros(NumOfTFs,1);   
%end
acceptF = 0; 

acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1,NumOfGenes);
numSamples = size(TFs,2);
%
for it = 1:(BurnInIters + Iters) 
    %
    % choose one sample for the TFS from the training set  
    
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
       for r=1:NumOfReplicas
       %  
       
       %for j=1:NumOfTFs
       %     
           
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           %ch = round(rand*numSamples) + 1; 
           %ch(ch>numSamples) = numSamples;
           
           %newLogProp = 0;  % forward Hastings Q(s_t+1 | s_t)
           %oldLogProp = 0;  % backward Hastings Q(s_t| s_t+1) 
           
           %% draw from the geometric distribution
           %kk = geornd(model.geo) + 1;
           %% propose the kkth nearest neighbor of the next TF
           %ch = simMat{j,r}(TFindex(j,r),kk);
           %
           % 
           %% forward Hastings Q(s_t+1 | s_t)
           %newLogProp = log(geopdf(kk-1,model.geo));
           %
           %% backward Hastings Q(s_t| s_t+1) 
           %bkk = find(simMat{j,r}(ch,:)==TFindex(j,r));
           %oldLogProp = log(geopdf(bkk-1,model.geo)); 
          
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           %LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           %LikParams1.TF(j,:,r) = TFs{ch}(j,:,r);
           
           % perform an evaluation of the likelihood p(Genes | F) 
           [newLogLik(r,:) predgen(r,:)] = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
       %end % num TFs loop
           
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
           PredictedGenes = predgen;
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
           [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
              
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
               PredictedGenes = predgen;
               oldLogLik = newLogLik;
            %   
           end
         end % num TFs loop
        %
    end % if end
    
    
    %
    % sample new kinetic parameters for each gene separately 
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
   
    % sample the interaction weights 
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
                if min(trW) >= 0
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
        LikW0 = feval(TrspaceW0, LikParams.W0);
        LogPriorWnew0 = feval(lnpriorW0, Wnew(end), model.prior.weight0);
        % >>>  interaction weights
        LikW = feval(TrspaceW, Wnew(1:NumOfTFs));
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights);
        
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
    
    % sample the noise of the likelihood when is free parameter.
    % The posterior is gamma so this step involves exact simualtion from
    % the gamma distribution
    if fixsigma2 == 0
       sumSquerrors1 = oldLogLik + 0.5*SizG*log(2*pi*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1));
       sumSquerrors1 = -2*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1).*sumSquerrors1;
       sumSquerrors = sum(sumSquerrors1,1);
                
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       LikParams.sigmas = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       % 
       okk = repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1);
       oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end
    
    
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.TFindex(cnt) = TFindex;
        %samples.predGenes(:,:,cnt) = PredictedGenes;
        samples.kinetics(:,cnt) = LikParams.kinetics;
        samples.Weights(:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
        if fixsigma2 == 0
           samples.sigmas(cnt) = LikParams.sigmas(1,1,1);
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

model.TFindex = TFindex; 

accRates.F = (acceptF/Iters)*100;
accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikval = remainRepsLikelihood(LikParams,  PredictedGenes, r, Gindex)
%
    Gns = LikParams.Genes(:,:,r);
    loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,r)),2)....
                - 0.5*sum(((Gns(Gindex,:) - PredictedGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,r),2);
    loglikval = loglikval';
% 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikvalTF = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, TFindex)
%    
  
  GnsTF = LikParams.GenesTF(:,:,r);
  loglikvalTF = - 0.5*sum(log(2*pi*LikParams.sigmasTF(TFindex,:,r)),2)....
                       - 0.5*sum(((GnsTF(TFindex,:) - PredictedGenesTF(:,LikParams.comIndsTF)).^2)./LikParams.sigmasTF(TFindex,:,r),2);
  loglikvalTF = loglikvalTF';
%  