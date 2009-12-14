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

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
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
   PredictedGenes(:,:,r) = predgen;
  
   % additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      [oldLogLikTF(r,:) predgen] = gpmtfLogLikelihoodGeneTF(model.Likelihood, F(:,:,r), r, 1:NumOfTFs);
      PredictedGenesTF(:,:,r) = predgen;
   end
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   [oldLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   PredictedGenes = predgen;
   % compute fast the additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      [oldLogLikTF(1,:) predgen] = gpmtfLogLikelihoodGeneTF(model.Likelihood, F, 1, 1:NumOfTFs);
      PredictedGenesTF = predgen;
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


% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]+eps);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights.mu, model.prior.weights.sigma2);

cnt = 0;

if strcmp(model.constraints.replicas,'free')
   acceptF = zeros(NumOfTFs,NumOfReplicas);
else
   acceptF = zeros(NumOfTFs,1);   
end

acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1,NumOfGenes);
%
for it = 1:(BurnInIters + Iters) 
    %
    % choose one sample for the TFS from the training set  
    for j=1:NumOfTFs
    if strcmp(model.constraints.replicas,'free') 
    %    
       for r=1:NumOfReplicas
       %     
           
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           gPerm = randperm(size(TFs,2));
           ch = gPerm(1); 
           
           newLogProp = 0;  % forward Hastings Q(s_t+1 | s_t)
           oldLogProp = 0;  % backward Hastings Q(s_t| s_t+1) 
           
           % draw from the geometric distribution
           kk = geornd(model.geo) + 1;
           % propose the kkth nearest neighbor of the next TF
           ch = simMat{j,r}(TFindex(j,r),kk);
           
            
           % forward Hastings Q(s_t+1 | s_t)
           newLogProp = log(geopdf(kk-1,model.geo));
           
           % backward Hastings Q(s_t| s_t+1) 
           bkk = find(simMat{j,r}(ch,:)==TFindex(j,r));
           oldLogProp = log(geopdf(bkk-1,model.geo)); 
           
           %[kk bkk]
           %[newLogProp oldLogProp]
           
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           LikParams1.TF(j,:,r) = TFs{ch}(j,:,r);
           
           if ~isfield(model.Likelihood,'GenesTF')
               % perform an evaluation of the likelihood p(Genes | F) 
               [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
       
               % Metropolis-Hastings to accept-reject the proposal
               [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(r,:),2), newLogProp, oldLogProp);
           else
               % perform an evaluation of the likelihood p(Genes | F) 
               [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);        
               [newLogLikTF predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams1, F(:,:,r), r, j); 
           
               % Metropolis-Hastings to accept-reject the proposal
               newL = sum(newLogLik(:)) + newLogLikTF;
               oldL = sum(oldLogLik(r,:),2) + oldLogLikTF(r,j);
               [accept, uprob] = metropolisHastings(newL, oldL, newLogProp, oldLogProp);
           end       
           %
           %[j r]
           %[sum(newLogLik(:))  sum(oldLogLik(r,:),2)]
           %pause
           
           if (it > BurnInIters)
               acceptF(j,r) = acceptF(j,r) + accept; 
           end
    
           % update protein F
           if accept == 1
               LikParams.TF(j,:,r) = TFs{ch}(j,:,r);
               TFindex(j,r) = ch;  
               PredictedGenes(:,:,r) = predgen;
               oldLogLik(r,:) = newLogLik;
               if isfield(model.Likelihood,'GenesTF')      
                  oldLogLikTF(r,j) = newLogLikTF;
                  PredictedGenesTF(j,:,r) = predgenTF;
               end
            %   
           end
       end % num Replicas loop
       %
    else
         % Replicas are coupled
           
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           gPerm = randperm(size(TFs,2));
           ch = gPerm(1); 
           
           newLogProp = 0;
           oldLogProp = 0; 
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           LikParams.TF(j,:) = TFs{ch}(j,:);
             
           newLogLik = zeros(NumOfReplicas,NumOfGenes);
           if ~isfield(model.Likelihood,'GenesTF')
               % perform an evaluation of the likelihood p(Genes | F)      
               [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
              
               % computed faster the remaining likelihood terms  
               for r=2:NumOfReplicas
                    newLogLik(r,:) = remainRepsLikelihood(LikParams1,  predgen, r, 1:NumOfGenes);
               end                 
               % Metropolis-Hastings to accept-reject the proposal
               [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), newLogProp, oldLogProp);
           else
               % perform an evaluation of the likelihood p(Genes | F)      
               [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);   
               [newLogLikTF(1) predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams, F, 1, j); 
           
               % computed faster the remaining likelihood terms  
               for r=2:NumOfReplicas
                   newLogLik(r,:) = remainRepsLikelihood(LikParams1,  predgen, r, 1:NumOfGenes);
                   newLogLikTF(r) = remainRepsLikelihoodTF(LikParams1, predgenTF, r, j);
               end    
             
               % Metropolis-Hastings to accept-reject the proposal
               newL = sum(newLogLik(:)) + sum(newLogLikTF(:));
               oldL = sum(oldLogLik(:)) + sum(oldLogLikTF(:,j),1);
               [accept, uprob] = metropolisHastings(newL, oldL, newLogProp, oldLogProp);
               %
           end       
       
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
         
               if isfield(model.Likelihood,'GenesTF')      
                  oldLogLikTF(:,j) = newLogLikTF(:);
                  PredictedGenesTF(j,:) = predgenTF;
               end
            %   
           end
        end % num Replicas loop
        %
    end % num TFs loop
    
    
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
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
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
        if posw == 1
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + log([LikParams.W(j,:), LikParams.W0(j)]+eps);     
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
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights.mu, model.prior.weights.sigma2);
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
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.TFindex = TFindex; 
        samples.predGenes{cnt} = PredictedGenes;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.Weights(:,:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
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