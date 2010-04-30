function [model PropDist samples accRates] = gpmtfOnlyDecayModel(model, mcmcoptions)
%  ODE model with only the decay term and the initial condition without having
%  the TFs 
%  Full Bayesian inference is applied 

% adaption phase
[model PropDist samples accRates] = gpmtfOnlyDecayModelAdapt(model, mcmcoptions.adapt); 
% training/sampling phase
[model PropDist samples accRates] = gpmtfOnlyDecayModelSample(model, PropDist, mcmcoptions.train);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model PropDist samples accRates] = gpmtfOnlyDecayModelAdapt(model, AdaptOps)
%
%

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
SizKin = size(model.Likelihood.kinetics,2);

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.05*ones(NumOfGenes,SizKin);

onlyPumaVar = 1; 
if sum(model.Likelihood.noiseModel.active(2:3)) > 0  
   onlyPumaVar = 0;
end

if onlyPumaVar == 0  
   PropDist.noiseModel = 0.02*ones(NumOfGenes,1); 
end

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qNoiseMbelow = 0.000001;   qNoiseMabove = 2;
epsilon = 0.1;

cnt = 0;
%
% do the adaption 
minAccR = 18; 
maxAccR = 38; 

nextbreak = 0; 
while 1
%
%  
   [model PropDist samples accRates] = gpmtfOnlyDecayModelSample(model,  PropDist, AdaptOps);
 
   accRateKin = accRates.Kin;
   accRateNoiseM = accRates.noiseM;
   
   if AdaptOps.disp == 1
      % 
      fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
       
      fprintf(1,'Acceptance rates for kinetic parameters (per gene))\n');
      disp(accRateKin);
      
      if onlyPumaVar == 0
         fprintf(1,'Acceptance rates for the noise parameters in the likelihood\n');
         disp(accRateNoiseM);
      end
      fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
      fprintf(1,'------------------------------- \n',cnt+1);
      %
   end
   
   
   % if you got a good acceptance rate, then stop
   if (min(accRateKin(:))>minAccR)  & (min(accRateNoiseM(:))>minAccR)  & (max(accRateKin(:))<maxAccR)  & (max(accRateNoiseM(:))<maxAccR) 
      if nextbreak == 1
          disp('END OF ADAPTION: acceptance rates OK');
          break;
      else
          nextbreak = 1;
      end
   end
    
   cnt = cnt + 1;
   % do not allow more than 150 iterations when you adapt the proposal distribution
   if cnt == 150   
       disp('END OF ADAPTION: acceptance rates OK');
       break;
       %  
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
      if accRateKin(j) < 20
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
           
   if onlyPumaVar == 0 
      %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT NOISE-MODEL PROPOSAL %%%%%%%%%%%%%%%%
      % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
      for j=1:NumOfGenes
         if accRateNoiseM(j) > 35
            % incease the covariance to reduce the acceptance rate
            PropDist.noiseModel(j,:) = PropDist.noiseModel(j,:) + epsilon*PropDist.noiseModel(j,:);
            if PropDist.noiseModel(j,1) > qNoiseMabove 
               PropDist.noiseModel(j,:) = qNoiseMabove*ones(1, size(PropDist.noiseModel(j,:), 2));
            end
         end
         if accRateNoiseM(j) < 20
            % decrease the covariance to incease the acceptance rate
            PropDist.noiseModel(j,:) = PropDist.noiseModel(j,:) - epsilon*PropDist.noiseModel(j,:);    
            if PropDist.noiseModel(j,1) < qNoiseMbelow 
               PropDist.noiseModel(j,:) = qNoiseMbelow*ones(1, size(PropDist.noiseModel(j,:), 2));
            end
          %
         end
       %
      end
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT NOISE-MODEL PROPOSAL %%%%%%%%%%%%%%%%
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model PropDist samples accRates] = gpmtfOnlyDecayModelSample(model, PropDist, trainOps)
%
%

BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
[NumOfGenes SizG NumOfReplicas] = size(Genes);

num_stored = floor(Iters/StoreEvery);
samples.kinetics = zeros(3, num_stored);
samples.LogL = zeros(1, num_stored);

%
onlyPumaVar = 1; 
if sum(model.Likelihood.noiseModel.active(2:3)) > 0
   onlyPumaVar = 0;
end

% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
for r=1:NumOfReplicas
%
   for j=1:NumOfGenes
     %  
     B = LikParams.kinetics(j,1);
     D = LikParams.kinetics(j,2);
     A = LikParams.kinetics(j,3);
     PredGenes(j,:) = B/D  + (A - B/D)*exp(-TimesG*D);
     %
   end 
   
   % evaluate the likelihood  
   sigmas = zeros(size(NumOfGenes,1), LikParams.numTimes);
   if LikParams.noiseModel.active(1) == 1
       sigmas = LikParams.noiseModel.pumaSigma2(:, : ,r);
   end
   if LikParams.noiseModel.active(2) == 1
       sigmas = sigmas + repmat(LikParams.noiseModel.sigma2(:)', 1, LikParams.numTimes ); 
   end

   if isfield(LikParams, 'crValMask')    
       oldLogLik(r,:) = - 0.5*sum(log(2*pi*sigmas(:, LikParams.crValMask)),2)....
                    - 0.5*sum(((LikParams.Genes(:, LikParams.crValMask,r)...
                    - PredGenes(:,LikParams.crValMask)).^2)./sigmas(:,LikParams.crValMask),2);
   else 
       oldLogLik(r,:) = - 0.5*sum(log(2*pi*sigmas),2)....
                    - 0.5*sum(((LikParams.Genes(:,:,r) - PredGenes).^2)./sigmas,2);
   end
%     
end

%
% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);
%

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
    end    
    %
%
end

cnt = 0;
acceptKin = zeros(1,NumOfGenes);
if onlyPumaVar == 0 
   accRateNoiseM = zeros(1, NumOfGenes); 
end

%
for it = 1:(BurnInIters + Iters) 
    %
    
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
       
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = zeros(1,NumOfReplicas);
        %
        
        B = LikParams1.kinetics(j,1);
        D = LikParams1.kinetics(j,2);
        A = LikParams1.kinetics(j,3);
        PredGenes = B/D  + (A - B/D)*exp(-TimesG*D); 
        % evaluate the new log likelihood 
        newLogLik = zeros(1,NumOfReplicas);
        for r=1:NumOfReplicas
            % 
            % evaluate the likelihood  
            sigmas = zeros(size(NumOfGenes,1), LikParams.numTimes);
            if LikParams.noiseModel.active(1) == 1
               sigmas = LikParams.noiseModel.pumaSigma2(j, : , r);
            end
            if LikParams.noiseModel.active(2) == 1
               sigmas = sigmas + repmat(LikParams.noiseModel.sigma2(j), 1, LikParams.numTimes ); 
            end
            
            if isfield(LikParams, 'crValMask')    
               newLogLik(r) = - 0.5*sum(log(2*pi*sigmas(LikParams.crValMask)),2)....
                        - 0.5*sum(((LikParams.Genes(j, LikParams.crValMask, r)...
                        - PredGenes(LikParams.crValMask)).^2)./sigmas(LikParams.crValMask),2);
            else 
               newLogLik(r) = - 0.5*sum(log(2*pi*sigmas),2)...
                         - 0.5*sum(((LikParams.Genes(j,:,r) - PredGenes).^2)./sigmas,2);
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
	   if exist('xgamrnd'),
	     newinvSigma2 = xgamrnd(anew,1./bnew);
	   else
	     newinvSigma2 = gamrnd(anew,1./bnew);
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
 
               B = LikParams1.kinetics(j,1);
               D = LikParams1.kinetics(j,2);
               A = LikParams1.kinetics(j,3);
               PredGenes = B/D  + (A - B/D)*exp(-TimesG*D); 
               % evaluate the new log likelihood 
               newLogLik = zeros(1,NumOfReplicas);
               for r=1:NumOfReplicas
               % 
                   % evaluate the likelihood  
                   sigmas = zeros(size(NumOfGenes,1), LikParams.numTimes);
                   if LikParams.noiseModel.active(1) == 1
                      sigmas = LikParams1.noiseModel.pumaSigma2(j, : , r);
                   end
                   if LikParams.noiseModel.active(2) == 1
                      sigmas = sigmas + repmat(newSigma2, 1, LikParams1.numTimes ); 
                   end
                   
                   if isfield(LikParams, 'crValMask')    
                      newLogLik(r) = - 0.5*sum(log(2*pi*sigmas(LikParams.crValMask)),2)....
                            - 0.5*sum(((LikParams.Genes(j, LikParams.crValMask, r)...
                            - PredGenes(LikParams.crValMask)).^2)./sigmas(LikParams.crValMask),2);
                   else
                      newLogLik(r) = - 0.5*sum(log(2*pi*sigmas),2)...
                            - 0.5*sum(((LikParams.Genes(j,:,r) - PredGenes).^2)./sigmas,2);
                   end  
               %              
               end
            
               % evaluation of the prior to the new white variance
               LikSigma2 = feval(TrspaceNoiseWHITE, newSigma2);
               LogPriorNewWHITE = feval(lnpriorNoiseWHITE, LikSigma2, model.prior.sigma2);
               % Metropolis-Hastings to accept-reject the proposal
               oldP = sum(oldLogLik(:,j),1) + oldLogPriorNoiseWHITE(j);
               newP = sum(newLogLik(:)) + LogPriorNewWHITE; 
         
               % compute the proposal distribution forward and backwards terms
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
        samples.kinetics(:,cnt) = LikParams.kinetics(:);
        if onlyPumaVar == 0
           if model.Likelihood.noiseModel.active(2) == 1 
              samples.sigma2(cnt) = LikParams.noiseModel.sigma2;
           end
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        samples.LogPrior(cnt) = sum(oldLogPriorKin(:));
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
if onlyPumaVar == 0
   if model.Likelihood.noiseModel.active(2) == 1 
      model.Likelihood.noiseModel.sigma2 = LikParams.noiseModel.sigma2;
   end
end

accRates.Kin = (acceptKin/Iters)*100;
if onlyPumaVar == 0 & ~(model.Likelihood.noiseModel.active(1) == 0 &  model.Likelihood.noiseModel.active(2) == 1) 
    accRates.noiseM = (accRateNoiseM/Iters)*100;
else
    accRates.noiseM = 25;
end

