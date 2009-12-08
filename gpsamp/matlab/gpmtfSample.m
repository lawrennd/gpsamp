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
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
if isfield(model.Likelihood,'GenesTF')
  LikkinTF = feval(TrspaceKin, LikParams.kineticsTF+eps); 
  oldLogPriorKinTF = feval(lnpriorKin, LikkinTF, model.prior.kinetics.a, model.prior.kinetics.b);
end      
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]+eps);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights.mu, model.prior.weights.sigma2);
% log prior of the lengthscale lengthscale 
lnpriorLengSc = ['ln',model.prior.GPkernel.lenghtScale.type,'pdf'];
%oldLogPriorLengSc = feval(lnpriorLengSc, exp(2*model.GP.logtheta(:,1)), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
for j=1:NumOfTFs
oldLogPriorLengSc(j) = feval(lnpriorLengSc, 2*model.GP{j}.logtheta(1), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
end

cnt = 0;
for j=1:NumOfTFs
    acceptF{j} = zeros(M(j),NumOfReplicas);
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
    % sample first all the TFs
    if 1 
    for r=1:NumOfReplicas
    for j=1:NumOfTFs
       %
       %
       % iterate between control points one-at-a-time 
       Fold = F(j,:,r);
       Fuold = Fu{j}(:,r);
       for i=istart(j):M(j)
       %
       % sample new control point i    
       Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],r);    
       
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
    
       %%%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
       %  
       %
         trform = 'lin'; model.Likelihood.singleAct; % 'exp' or 'linear'
         subplot(1,NumOfReplicas,r);
         Futmp = feval(model.Likelihood.singleAct,Fu{j}(:,r)'); 
         Ftmp = feval(model.Likelihood.singleAct,F(j,:,r));
         Funewtmp = feval(model.Likelihood.singleAct,Funew); 
         FFnewtmp = feval(model.Likelihood.singleAct,FFnew);
         %
         plot(Xu{j}, feval(trform,Futmp),'or','MarkerSize', 14,'lineWidth', 3);
         hold on;
         plot(TimesF, feval(trform,Ftmp),'g','lineWidth',4);
         if isfield(model,'groundtr') == 1
         GrFtmp = feval(model.Likelihood.singleAct,model.groundtr.F(j,:));   
         plot(TimesF, feval(trform,GrFtmp),'k','lineWidth',4);
         end
         title(j);
         pause(0.3);
         plot(Xu{j}(i), feval(trform,Futmp(i)),'oy','MarkerSize', 14,'lineWidth', 3);
         plot(Xu{j}(i), feval(trform,Funewtmp(i)), 'md','MarkerSize', 14, 'lineWidth',3);
         plot(TimesF, feval(trform,FFnewtmp(j,:)), '--b', 'lineWidth', 4); 
         pause(0.5);
         hold off;
       %  
       %  
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%% end of visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end % num TFs loop
    end % num Replicas loop
    end % zero-one if
   
    
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
        
        %-- Barenco data -- %
        %if j == 4
        %   KineticsNew([2 3]) = [0.8 1];
        %end
        %KineticsNew = kinetics(j,:);
        %-- Barenco data -- %
  
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          %newLogLik(r) = genesLogLikelihood(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
          newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
        end
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
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
        for r=1:NumOfReplicas
          % perform an evaluation of the likelihood p(GENES | TFs) 
          newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
          % 
        end
        %        
        Likkin = feval(TrspaceKin, TFKineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
        
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
          %newLogLik(r) = genesLogLikelihood(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
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
    end % if 0
    
    %
    % learn binary network connections (this not working for the current version) 
    if netLearn == 1
        for j=1:NumOfGenes
           % Gibbs sampling 
           LikParams1 = LikParams;
           newLogLik = zeros(sComb,NumOfReplicas);
           for co=1:sComb
               LikParams1.Net_X(j,:) = Comb(co,:);
               for r=1:NumOfReplicas
                 newLogLik(co,r) = genesLogLikelihood(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
               end
           end
           % sample one combination
           SumnewLogLik = sum(newLogLik,2); 
           SumnewLogLik = SumnewLogLik - max(SumnewLogLik(:));
           prob = exp(SumnewLogLik);
           prob = prob/sum(prob(:));
           co = sample_discrete(prob, 1, 1);
           LikParams.Net_X(j,:) = Comb(co,:);
           oldLogLik(:,j) = newLogLik(co,:)'; 
        end
    end
    %  
    
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
        
        newLogLik = [];
        for r=1:NumOfReplicas
        % call the function only with j gene expressions  
        %newLogLik(r) = genesLogLikelihood(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
        newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
        end

        %[model.prior.delays.prob(LikParams.Tausindex(j)) model.prior.delays.prob(LikParams1.Tausindex(j))]
        %[LikParams.Tausindex(j) LikParams1.Tausindex(j)]
        
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
    
    %model.GroundTruth.Taus
    %LikParams.Taus
    %pause
   
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
        newlogGP = - 0.5*NumOfReplicas*newLogDetK;
        oldlogGP = - 0.5*NumOfReplicas*PropDist.qF{j}.LogDetK;
        for r=1:NumOfReplicas     
           temp = invnewL*([F(j,:,r), Fu{j}(:,r)']'); 
           newlogGP = newlogGP - 0.5*temp'*temp;
           temp = PropDist.qF{j}.invL*([F(j,:,r), Fu{j}(:,r)']'); 
           oldlogGP = oldlogGP - 0.5*temp'*temp;
        end
        LogPriorLengScnew = feval(lnpriorLengSc, newlogEll, model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
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
        samples.predGenes{cnt} = PredictedGenes;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.Weights(:,:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
        samples.Taus(:,cnt) = LikParams.Taus(:);
        samples.Tausindex(:,cnt) = LikParams.Tausindex(:);
        if isfield(model.Likelihood,'GenesTF')
            samples.kineticsTF(:,:,cnt) = LikParams.kineticsTF;
            samples.LogLTF(cnt) = sum(oldLogLikTF(:));
            samples.predGenesTF{cnt} = PredictedGenesTF;
            if fixsigma2TF == 0
                samples.sigmasTF(:,:,cnt) = LikParams.sigmasTF(:,:,1);
            end
        end
        if fixsigma2 == 0
           samples.sigmas(:,:,cnt) = LikParams.sigmas(:,:,1);
        end
        if netLearn == 1
            samples.NetX(:,:,cnt) = LikParams.Net_X;
        end
        for jin=1:NumOfTFs
            samples.logthetas(jin,:,cnt) = model.GP{jin}.logtheta;
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
    accRates.F{j}(1,:) = 100*ones(1,NumOfReplicas);
    end
end
%

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
accRates.LengSc = (acceptLengSc/Iters)*100;
if isfield(model.Likelihood,'GenesTF')
   accRates.TFKin = (acceptTFKin/Iters)*100;
end

