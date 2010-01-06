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

% check if the observation noise is known/fixed
fixsigma2 = 0;
if strcmp(model.constraints.sigmas,'fixed')
    fixsigma2 = 1;
end

% check if the observation noise is known/fixed for the TF genes
if isfield(model.Likelihood,'GenesTF')
fixsigma2TF = 0;  
SizTFKin = size(model.Likelihood.kineticsTF,2);
if strcmp(model.constraints.sigmasTF,'fixed')
    fixsigma2TF = 1;
end
end


% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
    posw = 1;
end

% Construct the Comp binary matrix when you learn the 
% structure with Gibbs sampling --> not included in this version
netLearn = 0; 

% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

% the latent function values are also special parameters that appear in both 
% the likelihood and the GP prior
F = model.F; 
%F = model.groundtr.F;

% store the control variables
Fu = model.Fu; % function values  
Xu = model.Xu; % locations  
M = model.M;
n = SizF;

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
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
%save ok PredictedGenesTF PredictedGenes;
%oldLogLikTF
%sum(oldLogLik(:))
%model.Likelihood.kineticsTF 
%model.groundtr.kineticsTF
%model.Likelihood.kinetics 
%model.groundtr.kinetics
%[model.Likelihood.Taus; model.groundtr.Taus]
%pause
% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);
if isfield(model.Likelihood,'GenesTF')
  LikkinTF = feval(TrspaceKin, LikParams.kineticsTF+eps); 
  oldLogPriorKinTF = feval(lnpriorKin, LikkinTF, model.prior.kinetics);
end      
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights);
% log prior of the lengthscale lengthscale 
lnpriorLengSc = ['ln',model.prior.GPkernel.lenghtScale.type,'pdf'];
%oldLogPriorLengSc = feval(lnpriorLengSc, exp(2*model.GP.logtheta(:,1)), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
for j=1:NumOfTFs
oldLogPriorLengSc(j) = feval(lnpriorLengSc, 2*model.GP{j}.logtheta(1), model.prior.GPkernel.lenghtScale);
end

cnt = 0;
for j=1:NumOfTFs
    %
    if strcmp(model.constraints.replicas,'free')
      acceptF{j} = zeros(M(j),NumOfReplicas);
    else
      acceptF{j} = zeros(M(j),1);   
    end
    %
end
acceptKin = zeros(1,NumOfGenes);
acceptTFKin = zeros(1,NumOfTFs);
acceptW = zeros(1,NumOfGenes); 
acceptLengSc = zeros(1,NumOfTFs); 
%
for it = 1:(BurnInIters + Iters) 
    %
    %
    %F = model.groundtr.F;
    % Sample the TFs -----------------------------------------
    if 1 
    for j=1:NumOfTFs
    if strcmp(model.constraints.replicas,'free')   
         for r=1:NumOfReplicas
         %
         %
         Fold = F(j,:,r);
         Fuold = Fu{j}(:,r);
         % iterate between control points one-at-a-time 
         for i=istart(j):M(j)
         %
             % sample new control point i  
             %if abs(ppi) >= 1
                Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],r);    
             %else
             %   % use the underlaxed schme to sample the control point
             %   mu = PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],r);
             %   Fui = mu  + ppi*(Fu{j}(i,r) -  mu)  + sqrt(1-ppi^2)*(randn.*sqrt(PropDist.qF{j}.ku(i)));
             %end    
             
             Funew = Fu{j}(:,r)';
             Funew(i) = Fui;
    
             % sample the remaining points 
             cmu = PropDist.qF{j}.cmuMinus + Funew*PropDist.qF{j}.KInvKu;
             Fnew = gaussianFastSample(1, cmu, PropDist.qF{j}.L);
   
             FFnew = F(:,:,r);
             FFnew(j,:) = Fnew;
   
             if ~isfield(model.Likelihood,'GenesTF')
             % perform an evaluation of the likelihood p(Genes | F) 
             [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);
       
             % Metropolis-Hastings to accept-reject the proposal
             [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(r,:),2), 0, 0);
             else
             % perform an evaluation of the likelihood p(Genes | F) 
             [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);        
             [newLogLikTF predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams, FFnew, r, j); 
           
             % Metropolis-Hastings to accept-reject the proposal
             newL = sum(newLogLik(:)) + newLogLikTF;
             oldL = sum(oldLogLik(r,:),2) + oldLogLikTF(r,j);
             [accept, uprob] = metropolisHastings(newL, oldL, 0, 0);
             end
             %
             if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
             % 
                 visualize(model, F, Fu, FFnew, Funew, i, j, r);
             %  
             end
       
             %
             if (it > BurnInIters)
                 acceptF{j}(i,r) = acceptF{j}(i,r) + accept; 
             end
    
             % update protein F
            if accept == 1
                 F(j,:,r) = Fnew;
                 Fu{j}(:,r) = Funew';
                 PredictedGenes(:,:,r) = predgen;
                 oldLogLik(r,:) = newLogLik;
         
                 if isfield(model.Likelihood,'GenesTF')      
                    oldLogLikTF(r,j) = newLogLikTF;
                    PredictedGenesTF(j,:,r) = predgenTF;
                 end
            %   
            end
         end % num of control points loop
         %
         end % num Replicas loop
    else
         % Replicas are coupled
         % ---------------------------------------------
         % iterate between control points one-at-a-time 
         Fold = F(j,:);
         Fuold = Fu{j}(:,1);
         for i=istart(j):M(j)
         %
             % sample new control point i  
             %if abs(ppi) >= 1
                Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],1);        
             %else
             %   % use the underlaxed schme to sample the control point
             %   mu = PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],1);
             %   Fui = mu  + ppi*(Fu{j}(i,1) -  mu)  + sqrt(1-ppi^2)*(randn.*sqrt(PropDist.qF{j}.ku(i)));
             %end    
       
             Funew = Fu{j}(:,1)';
             Funew(i) = Fui;
             
             % sample the remaining points 
             cmu = PropDist.qF{j}.cmuMinus + Funew*PropDist.qF{j}.KInvKu;
             Fnew = gaussianFastSample(1, cmu, PropDist.qF{j}.L);
             FFnew = F;
             FFnew(j,:) = Fnew;          
             
             newLogLik = zeros(NumOfReplicas,NumOfGenes);
             if ~isfield(model.Likelihood,'GenesTF')
                 % perform an evaluation of the likelihood p(Genes | F)      
                 [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, 1, 1:NumOfGenes);
              
                 % computed faster the remaining likelihood terms  
                 for r=2:NumOfReplicas
                      newLogLik(r,:) = remainRepsLikelihood(LikParams,  predgen, r, 1:NumOfGenes);
                 end                 
                 % Metropolis-Hastings to accept-reject the proposal
                 [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), 0, 0);
             else
                 % perform an evaluation of the likelihood p(Genes | F)      
                 [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, 1, 1:NumOfGenes);   
                 [newLogLikTF(1) predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams, FFnew, 1, j); 
           
                 % computed faster the remaining likelihood terms  
                 for r=2:NumOfReplicas
                     newLogLik(r,:) = remainRepsLikelihood(LikParams,  predgen, r, 1:NumOfGenes);
                     newLogLikTF(r) = remainRepsLikelihoodTF(LikParams, predgenTF, r, j);
                 end    
             
                 % Metropolis-Hastings to accept-reject the proposal
                 newL = sum(newLogLik(:)) + sum(newLogLikTF(:));
                 oldL = sum(oldLogLik(:)) + sum(oldLogLikTF(:,j),1);
                 [accept, uprob] = metropolisHastings(newL, oldL, 0, 0);
                 %
             end
             
             %
             if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
             % 
                 visualize(model, F, Fu, FFnew, Funew, i, j, 1);
             %  
             end
       
             %
             if (it > BurnInIters)
                 acceptF{j}(i) = acceptF{j}(i) + accept; 
             end
    
             % update protein F
            if accept == 1
                 F(j,:) = Fnew;
                 Fu{j} = Funew';
                 PredictedGenes = predgen;
                 oldLogLik = newLogLik;
         
                 if isfield(model.Likelihood,'GenesTF')      
                    oldLogLikTF(:,j) = newLogLikTF(:);
                    PredictedGenesTF(j,:) = predgenTF;
                 end
            %   
            end
         end % num of control points loop
         %
        end % num Replicas loop
        %
    end % num TFs loop
    end % zero-one if
   
    
    if 1
    % sample new kinetic parameters for each TF-gene  
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
               [newLogLik(1,:), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
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
    
    
    
    if 1
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
    
    
    
    if 1 
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
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights);
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
       %for r=1:NumOfReplicas
       %oldLogLik(r,:) = genesLogLikelihood(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %end  
       okk = repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1);
       oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end
    
    % sample the noise of the TF-Genes likelihood when is free parameter.
    if isfield(model.Likelihood,'GenesTF')
    if fixsigma2TF == 0
       sumSquerrors1 = oldLogLikTF + 0.5*SizG*log(2*pi*repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1));
       sumSquerrors1 = -2*repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1).*sumSquerrors1;
       sumSquerrors = sum(sumSquerrors1,1);
                
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       LikParams.sigmasTF = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       %for r=1:NumOfReplicas
       %oldLogLik(r,:) = genesLogLikelihood(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %end  
       okk = repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1);
       oldLogLikTF = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end
    end

    if 1 
    if model.Likelihood.tauMax < 0
    % sample the delay parameters in the ODEs 
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
   
   
    %
    %
    if 1
    % sample the lenghtscale of the TFs 
    for j=1:NumOfTFs
        %
        % to samples the lengthscale you need to evaluate 
        % the GP pior compute the new and the current/old GP prior 
        newlogEll = randn.*sqrt(PropDist.LengSc(j)) + 2*model.GP{j}.logtheta(1);
        %newlogEll = model.GP{j}.logtheta(1);
        newEll2 = exp(newlogEll);  
        newK = exp(-(0.5/newEll2)*model.GP{j}.X2) + ...
               exp(2*model.GP{j}.logtheta(end))*eye(size(model.GP{j}.X2,1)); 
        % compute the Cholesky decomposition of the new K
        [newL,er]=jitterChol(newK);
        newL = newL';
        % evaluate the new log GP prior value 
        invnewL = newL\eye(SizF+M(j));
        newLogDetK = 2*sum(log(diag(newL)));
        if strcmp(model.constraints.replicas,'free')
            %
            newlogGP = - 0.5*NumOfReplicas*newLogDetK;
            oldlogGP = - 0.5*NumOfReplicas*PropDist.qF{j}.LogDetK;
            for r=1:NumOfReplicas     
                temp = invnewL*([F(j,:,r), Fu{j}(:,r)']'); 
                newlogGP = newlogGP - 0.5*temp'*temp;
                temp = PropDist.qF{j}.invL*([F(j,:,r), Fu{j}(:,r)']'); 
                oldlogGP = oldlogGP - 0.5*temp'*temp;
            end
            %
        else
            newlogGP = - 0.5*newLogDetK;
            oldlogGP = - 0.5*PropDist.qF{j}.LogDetK;     
            temp = invnewL*([F(j,:), Fu{j}(:,1)']'); 
            newlogGP = newlogGP - 0.5*temp'*temp;
            temp = PropDist.qF{j}.invL*([F(j,:), Fu{j}(:,1)']'); 
            oldlogGP = oldlogGP - 0.5*temp'*temp;
        end
            
        LogPriorLengScnew = feval(lnpriorLengSc, newlogEll, model.prior.GPkernel.lenghtScale);
        % Metropolis-Hastings to accept-reject the proposal
        oldlogGP = oldlogGP + oldLogPriorLengSc(j);
        newlogGP = newlogGP + LogPriorLengScnew; 
        
        
        %YnewlogEll = 2*model.GP{j}.logtheta(1);
        %%newlogEll = model.GP{j}.logtheta(1);
        %YnewEll2 = exp(YnewlogEll);  
        %YnewK = kernCompute(model.GP{j}, [TimesF(:); Xu(:)]);
        %%ok     = exp(-(0.5/YnewEll2)*model.GP{1}.X2) + ...
        %%        exp(2*model.GP{j}.logtheta(end))*eye(size(model.GP{1}.X2,1));
            
           
        %% compute the Cholesky decomposition of the new K
        %[YnewL,er]=jitterChol(YnewK);
        %YnewL = YnewL';
        %% evaluate the new log GP prior value 
        %YinvnewL = YnewL\eye(SizF+M);
        %YnewLogDetK = 2*sum(log(diag(YnewL)));
        %YnewlogGP = - 0.5*NumOfReplicas*YnewLogDetK;
        %for r=1:NumOfReplicas     
        %   temp = YinvnewL*([F(j,:,r), Fu(j,:,r)]'); 
        %   YnewlogGP = YnewlogGP - 0.5*temp'*temp;
        %end
        %
        %YLogPriorLengScnew = feval(lnpriorLengSc, YnewlogEll, model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b); 
        %YnewlogGP = YnewlogGP + oldLogPriorLengSc(j);
        
        %if mod(it,100) == 0
        %[oldlogGP  newlogGP]
        %[oldLogPriorLengSc(j) LogPriorLengScnew]
        %[exp(2*model.GP{j}.logtheta(1)) newEll2]
        %end
        %
        
        [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, 0, 0);
        %%%%%%%%%%%%%%%%  start accept/update proposal for the lengthscale %%%%%%%%%%%% 
        if accept == 1
           U = n+1:n+M(j);
           model.GP{j}.logtheta(1) = newlogEll/2;
           oldLogPriorLengSc(j) = LogPriorLengScnew;
           %newEll2
           %pause
           PropDist.qF{j}.K = newK;
           PropDist.qF{j}.invL = invnewL; 
           PropDist.qF{j}.LogDetK = newLogDetK;
           
           % change the proposal distribution for the TF
           % compute the conditional GP prior given the control variables
           [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', newK, 1:n, U);
           [L,er]=jitterChol(cSigma);
           if er>0, L = real(sqrtm(cSigma)); end
           PropDist.qF{j}.cmuMinus = cmuMinus; 
           PropDist.qF{j}.cSigma = cSigma;
           PropDist.qF{j}.KInvKu = KInvKu;
           PropDist.qF{j}.L = L;
           for i=1:M(j)
           %  
              G = [1:i-1, i+1:M(j)];  
              [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
           %
           end
           PropDist.qF{j}.alpha = alpha;
           PropDist.qF{j}.ku = ku;
           PropDist.qF{j}.KInvK = KInvK;
           clear alpha  ku  KInvK;
        end
        % 
        if (it > BurnInIters) 
           acceptLengSc(j) = acceptLengSc(j) + accept;
        end
        %%%%%%%%%%%%%%%%%%%%%%% end accept/update proposal for the lengthscale %%%%%%%%%%%%%%%%
        %
    end
    end
    
    %
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.F{cnt} = F;
        samples.Fu{cnt} = Fu;
        %samples.predGenes{cnt} = PredictedGenes;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.Weights(:,:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
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
            if fixsigma2TF == 0
                samples.sigmasTF(:,:,cnt) = LikParams.sigmasTF(:,:,1);
            end
        end
        if fixsigma2 == 0
           samples.sigmas(:,:,cnt) = LikParams.sigmas(:,:,1);
        end
        %if netLearn == 1
        %    samples.NetX(:,:,cnt) = LikParams.Net_X;
        %end
        for jin=1:NumOfTFs
            samples.logthetas(jin,cnt) = model.GP{jin}.logtheta(1);
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        %save(trainOps.ResFile,'samples','model');
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
    model.Likelihood.sigmasTF = LikParams.sigmasTF;
end

model.Likelihood.sigmas = LikParams.sigmas;
if netLearn == 1
model.Likelihood.Net_x = LikParams.Net_X;
end
model.F = F;
model.Fu = Fu;
%
for j=1:NumOfTFs
    accRates.F{j} = (acceptF{j}/Iters)*100; 
    if istart == 2
       if strcmp(model.constraints.replicas,'free')  
          accRates.F{j}(1,:) = 100*ones(1,NumOfReplicas);
       else
          accRates.F{j}(1) = 100; 
       end
    end
end
%

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
accRates.LengSc = (acceptLengSc/Iters)*100;
if isfield(model.Likelihood,'GenesTF')
   accRates.TFKin = (acceptTFKin/Iters)*100;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function visualize(model,F,Fu,FFnew,Funew,i,j,r)
%
trform = 'lin'; % model.Likelihood.singleAct; % 'exp' or 'linear'
if strcmp(model.constraints.replicas,'free')
subplot(1,model.Likelihood.numReplicas,r);    
end
Futmp = feval(model.Likelihood.singleAct,Fu{j}(:,r)');
Ftmp = feval(model.Likelihood.singleAct,F(j,:,r));
Funewtmp = feval(model.Likelihood.singleAct,Funew);
FFnewtmp = feval(model.Likelihood.singleAct,FFnew);
%
plot(model.Xu{j}, feval(trform,Futmp),'or','MarkerSize', 14,'lineWidth', 3);
hold on;
plot(model.Likelihood.TimesF, feval(trform,Ftmp),'g','lineWidth',4);
if isfield(model,'groundtr') == 1
    GrFtmp = feval(model.Likelihood.singleAct,model.groundtr.F(j,:));
    plot(model.Likelihood.TimesF, feval(trform,GrFtmp),'k','lineWidth',4);
end
title(j);
pause(0.3);
plot(model.Xu{j}(i), feval(trform,Futmp(i)),'oy','MarkerSize', 14,'lineWidth', 3);
plot(model.Xu{j}(i), feval(trform,Funewtmp(i)), 'md','MarkerSize', 14, 'lineWidth',3);
plot(model.Likelihood.TimesF, feval(trform,FFnewtmp(j,:)), '--b', 'lineWidth', 4);
pause(0.5);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%% end of visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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