function [model PropDist samples accRates] = gpmtfBaselineMultTFModel(model, mcmcoptions)
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
SizKinTF = 2;


% if the initial condition of the differential equations is zero, then set accordingly the
% corresponding kinetic parameter (initial condition)
for j=1:NumOfGenes
if model.constraints.InitialConds(j) == 0
     model.Likelihood.kinetics(:,4) = model.Likelihood.kinetics(:,1)./model.Likelihood.kinetics(:,2);      
end
end

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.5*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.5*ones(NumOfGenes,NumOfTFs+1);


onlyPumaVar = 1; 
if model.Likelihood.noiseModel.active(2) > 0  
   onlyPumaVar = 0;
end

if onlyPumaVar == 0 
   % white nosie variance per gene possibly added to the PUMA variances
   PropDist.noiseModel = 0.5*ones(1, NumOfGenes); 
end


% additional proposal distribution for the TF kinetic parameters
if isfield(model.Likelihood,'GenesTF')    
  PropDist.TFkin = 0.5*ones(NumOfTFs,2);
end

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
qNoiseMbelow = 0.000001;  qNoiseMabove = 2;
epsilon = 0.1;
cnt = 0;
opt = 0.25;

% create the piece-wise linear function that passes through the Tf mRNA 
% values 
F = zeros( size(model.Likelihood.GenesTF, 1), size(model.Likelihood.TimesF,2), NumReplicas);
for R=1:NumReplicas
   for j=2:size(model.Likelihood.TimesG, 2)
      b = model.Likelihood.TimesG(j); 
      a = model.Likelihood.TimesG(j-1);
      x = a:model.Likelihood.step:b;
      F(:,model.Likelihood.comInds(j-1):model.Likelihood.comInds(j), R) = repmat( model.Likelihood.GenesTF(:,j,R), 1, size(x,2)).*repmat( (x-a)/(b-a), size(F,1), 1)  + repmat( model.Likelihood.GenesTF(:,j-1,R), 1, size(x,2)).*repmat( (b-x)/(b-a), size(F,1), 1); 
   end
end
model.F = F;

%for R=1:NumReplicas
%for j=1:size(F,1)
%  plot( model.Likelihood.TimesF, F(j,:,R),'b');
%  hold on;
%  plot( model.Likelihood.TimesG, model.Likelihood.GenesTF(j,:,R),'r+');
%  [j R]
%  pause;
%  hold off;
%end
%end


%
% do the adaption 
while 1
%
%  
   
   [model PropDist samples accRates] = gpmtfBaselineMultTFModelSample(model, PropDist, AdaptOps);
   
   %samples
   %Likelihood =  model.Likelihood; 
   %for r=1:NumReplicas
   %   Likelihood.kinetics = samples.kinetics;
   %   Likelihood.W = samples.W; 
   %   Likelihood.W0 = samples.W0;
   %   Likelihood.kineticsTF = samples.kineticsTF;
   %   oldLogLik(r,:) = gpmtfLogLikelihoodGene(Likelihood, model.F(:,:,r), r, 1:NumOfGenes);   
   %end
   
   accRateKin = accRates.Kin;
   accRateW = accRates.W; 
   accRateNoiseM = accRates.noiseM;
   %
   if isfield(model.Likelihood,'GenesTF')
       accRateTFKin = accRates.TFKin;
   end
   
   if AdaptOps.disp == 1
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
   
   if isfield(model.Likelihood,'GenesTF')
      fprintf(1,'Acceptance Rates for kinetic parameters (per TF-gene))\n');
      disp(accRateTFKin);
   end

   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
  
   if onlyPumaVar == 0
         fprintf(1,'Acceptance rates for the noise parameters in the likelihood\n');
         disp(accRateNoiseM);
   end
   fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
   fprintf(1,'------------------------------- \n',cnt+1);
   end   
       
   
   if  (min(accRateKin(:))>15) & (min(accRateW(:))>15) 
        disp('END OF ADAPTION: acceptance rates OK');
        break;
   end
   cnt = cnt + 1;
   % do not allow more than 100 iterations when you adapt the proposal distribution
   if cnt == 100
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
   end
  %
  %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT NOISE-MODEL PROPOSAL %%%%%%%%%%%%%%%%
%
end

% run the sampler
[model PropDist samples accRates] = gpmtfBaselineMultTFModelSample(model, PropDist, mcmcoptions.train); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model PropDist samples accRates] = gpmtfBaselineMultTFModelSample(model, PropDist, trainOps)
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
SizKin = size(model.Likelihood.kinetics,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);
NumOfTFs = model.Likelihood.numTFs;

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


LikParams = model.Likelihood;
% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
if strcmp(model.constraints.replicas,'free')
   for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   %oldLogLik(r,:) = gpmtfBaselineMultTFModelLogLikelihood(model.Likelihood, model.F(:,:,r), r, 1:NumOfGenes);
   % this assumes that the single activation function is linear
   oldLogLik(r,:) = gpmtfLogLikelihoodGene(model.Likelihood, model.F(:,:,r), r, 1:NumOfGenes);
   %     
   end
%   
end
%

cnt = 0;
maxLogLik = -Inf;
acceptKin = zeros(1,NumOfGenes);
acceptTFKin = zeros(1,NumOfTFs);
acceptW = zeros(1,NumOfGenes); 
accRateNoiseM = zeros(1, NumOfGenes); 
%
for it = 1:(BurnInIters + Iters) 
    % 
    % *
    % SAMPLE KINETIC PARAMETERS FOR THE TF TRANSLATION MODEL
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
               %newLogLik(r,:) = gpmtfBaselineMultTFModelLogLikelihood(LikParams1, model.F(:,:,r), r, 1:NumOfGenes);
               % this assumes that the single activation function is linear
               newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, 1:NumOfGenes);
               % 
           end
        end
        %        
        
        % Metropolis-Hastings to accept-reject the proposal
        newP = sum(newLogLik(:));
        oldP = sum(oldLogLik(:));
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.kineticsTF(j,:) = TFKineticsNew;
           oldLogLik = newLogLik;
        end
        %
        if (it > BurnInIters) 
           acceptTFKin(j) = acceptTFKin(j) + accept;
        end
        %
    end
    end
    % *
    % END SAMPLE KINETIC PARAMETERS FOR THE TF TRANSLATION MODEL
    
    % *
    % SAMPLE KINETIC PARAMETERS FOR EACH GENE ODE 
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
              %newLogLik(r) = gpmtfBaselineMultTFModelLogLikelihood(LikParams1, model.F(:,:,r), r, j);
              % this assumes that the single activation function is linear
              newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
           end
        end
        %        
        
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1);
        newP = sum(newLogLik(:));
        
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        %[accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
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
    % *
    % END SAMPLE KINETIC PARAMETERS FOR EACH GENE ODE
   
    % *
    % SAMPLE INTERACTION WEIGHTS FOR EACH GENE
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
               %newLogLik(r) = gpmtfBaselineMultTFModelLogLikelihood(LikParams1, model.F(:,:,r), r, j);
               % this assumes that the single activation function is linear
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
           end
        end
        
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1);
        newP = sum(newLogLik(:)); 
        %
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
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
    % *
    % END SAMPLE INTERACTION WEIGHTS FOR EACH GENE
    
    % * 
    % SAMPLE THE NOISE MODEL IN THE LIKELIHOOD of THE GENES 
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
                       %newLogLik(r) = gpmtfBaselineMultTFModelLogLikelihood(LikParams1, model.F(:,:,r), r, j);
                       % this assumes that the single activation function is linear
                       newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, model.F(:,:,r), r, j);
                   end
               end
       
               % Metropolis-Hastings to accept-reject the proposal
               oldP = sum(oldLogLik(:,j),1);
               newP = sum(newLogLik(:)); 
         
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
               end
               %
               if (it > BurnInIters) 
                  accRateNoiseM(j) = accRateNoiseM(j) + accept;
               end
               %
           end % for loop
       %
    end % if statement
    % END SAMPLE THE NOISE MODEL IN THE LIKELIHOOD of THE GENES 
    % *
 
    %
    % keep samples after burn in
    if (it > BurnInIters)  
        %
        cnt = cnt + 1;
        % store the parameters with the maximum likelihood value 
        if sum(oldLogLik(:)) > maxLogLik  
            samples.kinetics(:,:,1) = LikParams.kinetics;
            samples.W(:,:,1) = LikParams.W;
            samples.W0(:,1) = LikParams.W0;
            if isfield(model.Likelihood,'GenesTF')
               samples.kineticsTF(:,:,1) = LikParams.kineticsTF;
            end
            if onlyPumaVar == 0
              samples.sigma2(:,1) = LikParams.noiseModel.sigma2';
            end
            maxLogLik = sum(oldLogLik(:));  
            samples.LogL = maxLogLik;
        end
        %% store all log likelihoods jsut for reference 
        %samples.LogL(cnt) = sum(oldLogLik(:));
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
if isfield(model.Likelihood,'GenesTF')
    model.Likelihood.kineticsTF = LikParams.kineticsTF;
end

if onlyPumaVar == 0
    model.Likelihood.noiseModel.sigma2 = LikParams.noiseModel.sigma2;
end
accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
if isfield(model.Likelihood,'GenesTF')
   accRates.TFKin = (acceptTFKin/Iters)*100;
end
accRates.noiseM = (accRateNoiseM/Iters)*100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loglikval, PredGenes, PredTFs] = gpmtfBaselineMultTFModelLogLikelihood(LikParams, F, R, Gindex)
%
%

Ngenes = size(Gindex,2);
% Translational ODE model for the TFs  
uu = LikParams.TimesF;
Delta = uu(2)-uu(1); 
PredTFs = zeros(size(F));
for j=1:LikParams.numTFs
    %
    D = LikParams.kineticsTF(j,1);
    S = LikParams.kineticsTF(j,2);    
    % Trapezoid rule of numerical integration 
    expD = exp(uu*D);
    tmp = F(j,:).*expD; 
    expD = exp(-uu*D);
    PredTFs(j,2:end) = (0.5*Delta)*cumsum( ( tmp(1:end-1) + tmp(2:end)) );
    PredTFs(j,:) = S*(expD.*PredTFs(j,:));
    %
end

fx = jointactFunc(LikParams, PredTFs, Gindex);
uu = LikParams.TimesF(LikParams.startTime:end);
Delta = uu(2)-uu(1); 
PredGenes = zeros(Ngenes,size(uu,2));
for m=1:Ngenes
    j = Gindex(m);
    %
    B = LikParams.kinetics(j,1);
    D = LikParams.kinetics(j,2);
    S = LikParams.kinetics(j,3);   
    A = LikParams.kinetics(j,4);
  
    % Trapezoid rule of numerical integration
    %IntVals = exp(D*uu).*fx(m,:);
    ffx = exp(D*uu).*fx(m,LikParams.Tausindex(j):LikParams.Tausindex(j)+LikParams.sizTime-1);
    IntVals = zeros(size(ffx));
    IntVals(2:end) = .5 * Delta*cumsum(ffx(1:end-1) + ffx(2:end));
    %IntVals = Delta*cumtrapz(IntVals);
    %IntVals = IntVals(comInds);
    
    % Simpson rule of integration 
    %ffx = exp(D*uu).*fx(m,:); 
    %IntVals = ffx;
    %IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    %IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    %IntVals = cumsum(IntVals);
    %IntVals = (Delta/3)*IntVals;
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
    expD = exp(-uu*D);
    PredGenes(m,:) = B/D  + (A - B/D)*expD + S*(expD.*IntVals);
   %
end    

if LikParams.noiseModel.active(1) == 1
    sigmas = LikParams.noiseModel.pumaSigma2(Gindex, : , R);
else
    sigmas = zeros(size(Gindex(:),1), LikParams.numTimes);
end
if LikParams.noiseModel.active(2) == 1
    sigmas = sigmas + repmat(LikParams.noiseModel.sigma2(Gindex)', 1, LikParams.numTimes ); 
end
if isfield(LikParams, 'crValMask')    
   loglikval = - 0.5*sum(log(2*pi*sigmas(:, LikParams.crValMask)),2)....
               - 0.5*sum(((LikParams.Genes(Gindex, LikParams.crValMask,R)...
               - PredGenes(:,LikParams.comInds(LikParams.crValMask))).^2)./sigmas(:,LikParams.crValMask),2);
else 
     loglikval = - 0.5*sum(log(2*pi*sigmas),2) ...
               - 0.5*sum(((LikParams.Genes(Gindex,:,R) - PredGenes(:,LikParams.comInds)).^2)./sigmas,2);
end            
loglikval = loglikval(:)';   

