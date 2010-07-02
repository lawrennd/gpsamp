function [model PropDist samples accRates] = gpmtfBaselineMultTFModelTest(model, mcmcoptions)
%
%

AdaptOps = mcmcoptions.adapt;

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
%F = zeros(NumOfTFs, SizF, NumReplicas);
%% fake variable
%model.F = F;

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.05*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.05*ones(NumOfGenes,NumOfTFs+1);

onlyPumaVar = 1; 
if sum(model.Likelihood.noiseModel.active(2:3)) > 0  
   onlyPumaVar = 0;
end

if onlyPumaVar == 0 
   % they are (at most) three parameters in the likelihood noise model  
   % 1) white nosie variance, rbf variance, rbf lengthscale 
   PropDist.noiseModel = 0.02*ones(NumOfGenes,3); 
end

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
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
   [model PropDist samples accRates] = gpmtfBaselineMultTFModelSampleTest(model, PropDist, AdaptOps);
 
   %samples
   %Likelihood =  model.Likelihood; 
   %for r=1:NumReplicas
   %   Likelihood.kinetics = [3.3382e-04 3.2266e-04 0.7137 1.2491e-04]; % samples.kinetics;
   %   Likelihood.W = 44.2867; %samples.W; 
   %   Likelihood.W0 = 229.7094; % samples.W0;
   %   Likelihood.noiseModel.sigma2 = 1.9227; %samples.sigma2; 
   %   oldLogLik(r,:) = gpmtfLogLikelihoodGene(Likelihood, model.F(:,:,r), r, 1:NumOfGenes);   
   %end
   %
   %sum(sum(oldLogLik))
   %pause
   
   
   accRateKin = accRates.Kin;
   accRateW = accRates.W;
   accRateNoiseM = accRates.noiseM;
   
   if AdaptOps.disp == 1
      % 
      fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
       
      fprintf(1,'Acceptance rates for kinetic parameters (per gene))\n');
      disp(accRateKin);
 
      fprintf(1,'Acceptance rates for interaction weights (per gene)\n');
      disp(accRateW);
  
      if onlyPumaVar == 0
         fprintf(1,'Acceptance rates for the noise parameters in the likelihood\n');
         disp(accRateNoiseM);
      end
      fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
      fprintf(1,'------------------------------- \n',cnt+1);
      %
   end
   
   
   % if you got a good acceptance rate, then stop
   if (min(accRateKin(:))>minAccR) & (min(accRateW(:))>minAccR) & (min(accRateNoiseM(:))>minAccR)  & (max(accRateKin(:))<maxAccR) & (max(accRateW(:))<maxAccR)  & (max(accRateNoiseM(:))<maxAccR) 
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
      if accRateW(j) < 20
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
%
end

% run the sampler
[model PropDist samples accRates] = gpmtfBaselineMultTFModelSampleTest(model,  PropDist, mcmcoptions.train); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model PropDist samples accRates] = gpmtfBaselineMultTFModelSampleTest(model, PropDist, trainOps)
%
%
%
%
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


% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   oldLogLik(r,:) = gpmtfLogLikelihoodGene(model.Likelihood, model.F(:,:,r), r, 1:NumOfGenes);
   %     
end


cnt = 0;
acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1, NumOfGenes);
if onlyPumaVar == 0 
   accRateNoiseM = zeros(1, NumOfGenes); 
end

maxLogLik = -Inf; 

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
        for r=1:NumOfReplicas
            % call the function only with j gene expressions  
            newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
        end
 
        % Metropolis-Hastings to accept-reject the proposal
        oldP = model.Likelihood.invT*sum(oldLogLik(:,j),1);
        newP = model.Likelihood.invT*sum(newLogLik(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.kinetics(j,:) = KineticsNew;
           oldLogLik(:,j) = newLogLik(:);
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
        for r=1:NumOfReplicas
            % call the function only with j gene expressions  
            newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
        end
        
        % Metropolis-Hastings to accept-reject the proposal
        oldP = model.Likelihood.invT*sum(oldLogLik(:,j),1);
        newP = model.Likelihood.invT*sum(newLogLik(:)); 
         
        [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:);
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
               for r=1:NumOfReplicas
                   % call the function only with j gene expressions  
                   newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
               end
               
               % Metropolis-Hastings to accept-reject the proposal
               oldP = model.Likelihood.invT*sum(oldLogLik(:,j),1);
               newP = model.Likelihood.invT*sum(newLogLik(:)); 
         
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
               end
               %
               if (it > BurnInIters) 
                  accRateNoiseM(j) = accRateNoiseM(j) + accept;
               end
               %
           end % for loop 
       %
    end % if statement
    % END SAMPLE THE NOISE MODEL IN THE LIKELIHOOD
    % *
    
    % *
    % KEEP SAMPLES AFTER BURN IN
    if (it > BurnInIters)
        %
        cnt = cnt + 1;
        % store the parameters with the maximum likelihood value 
        if sum(oldLogLik(:)) > maxLogLik  
            samples.kinetics(:,:,1) = LikParams.kinetics;
            samples.W(:,:,1) = LikParams.W;
            samples.W0(:,1) = LikParams.W0;
            if onlyPumaVar == 0
              if model.Likelihood.noiseModel.active(2) == 1 
                samples.sigma2 = LikParams.noiseModel.sigma2;
              end
            end
            %if isfield(model.Likelihood,'GenesTF')
            %   samples.kineticsTF(:,:,1) = LikParams.kineticsTF;
            %end
            if onlyPumaVar == 0
              samples.sigma2(:,1) = LikParams.noiseModel.sigma2';
            end
            maxLogLik = sum(oldLogLik(:));  
            samples.LogL = maxLogLik;
        end
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
end

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
if onlyPumaVar == 0 & ~(model.Likelihood.noiseModel.active(1) == 0 &  model.Likelihood.noiseModel.active(2) == 1) 
    accRates.noiseM = (accRateNoiseM/Iters)*100;
else
    accRates.noiseM = 25;
end



