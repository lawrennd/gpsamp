function [model PropDist samples accRates] = gpmtfTestGenesSample(model, PropDist, trainOps)
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

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
for r=1:NumOfReplicas
  %
  % evaluate the likelihood 
  [oldLogLik(r,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
  PredictedGenes(:,:,r) = predgen;  
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
acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1,NumOfGenes); 
%
for it = 1:(BurnInIters + Iters) 
    %
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
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
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
      
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
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

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;

