function [model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFSample(model, PropDist, Genes, TimesG, TimesF, trainOps)
% Description: Draw a set of samples from the Bayesian differential
%              equation model
%
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the gene expressions for 
%                    all genes for all times and replicas
%         -- TimesG: the time points where gene expression are evaluated 
%         -- TimesF: the times where the GP fucntion are evaluated
%                    TimesF>>TimesG
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
%         -- accRateF: acceptance rates (after burn-in) for the TFs
%         -- accRateKin: >>  >> for the kinetic parameters 
%         -- accRateW: >>  >> for the interaction weigths or any parameters 
%                      that is used in the activation fucntion of the TF (e.g.  
%                       the gamma in the Michaels-Menten activation)  
%         -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%                           kernel 
%


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;
SizF = size(TimesF,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);


NumOfTFs = model.Likelihood.NumOfTFs;
istart = ones(NumOfTFs,1);
for j=1:NumOfTFs
if model.Constraints.Ft0(j)==0
   % do not sample the first control point so that the function 
   % will be fixed at the time t=0
   istart(j)=2;
end
end

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.Constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

% check if the observation noise is known/fixed
fixsigma2 = 0;
if strcmp(model.Constraints.sigma2,'fixed')
    fixsigma2 = 1;
end

% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
    posw = 1;
end

netLearn = 0; 
% Construct the Comp binary matrix when you learn the 
% structure with Gibbs sampling --> not included in this version


% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams.sigmas = model.Likelihood.sigmas;
LikParams.kinetics(:,1) = model.Likelihood.B'; % basal rates
LikParams.kinetics(:,2) = model.Likelihood.D'; % decay rates
LikParams.kinetics(:,3) = model.Likelihood.S'; % sensitivities 
LikParams.kinetics(:,4) = model.Likelihood.A'; % initila conditions
SizKin = size(LikParams.kinetics,2);
LikParams.NumOfGenes = NumOfGenes;
%% Gene - TFs interaction weights 
LikParams.W = model.Likelihood.W.*model.Constraints.W;
%% gene bias term in the regulatory part
LikParams.W0 = model.Likelihood.W0;
LikParams.TFjointAct = model.Likelihood.TFjointAct;
LikParams.TFsingleAct = model.Likelihood.TFsingleAct;
LikParams.TFjointActBin = model.Likelihood.TFjointActBin;
if strcmp(model.Likelihood.TFjointAct,'michMenten')
   LikParams.Net_X = model.Likelihood.Net_X;
end
% the latent function values are also special parameters that appear in both 
% the likelihood and the GP prior
F = model.F; 

% store the control variables
Fu = model.Fu; % function values  
Xu = model.Xu; % locations  
M = size(Fu,2);
n = SizF;
U = n+1:n+M; 

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
for r=1:NumOfReplicas
  [oldLogLik(r,:) predgen] = logLTFdiffEquation(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
  PredictedGenes(:,:,r) = predgen;
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
% log prior of the lengthscale lengthscale 
lnpriorLengSc = ['ln',model.prior.GPkernel.lenghtScale.type,'pdf'];
%oldLogPriorLengSc = feval(lnpriorLengSc, exp(2*model.GP.logtheta(:,1)), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
oldLogPriorLengSc = feval(lnpriorLengSc, 2*model.GP.logtheta(:,1), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);

SizKin = size(LikParams.kinetics,2);
cnt = 0;
acceptF = zeros(NumOfTFs,M,NumOfReplicas);
acceptKin = zeros(1,NumOfGenes); 
acceptW = zeros(1,NumOfGenes); 
acceptLengSc = zeros(1,NumOfTFs); 
%
for it = 1:(BurnInIters + Iters) 
    %
    %
    % sample first all the TFs
    for r=1:NumOfReplicas
    for j=1:NumOfTFs
       %
       %
       % iterate between control points one-at-a-time 
       Fold = F(j,:,r);
       Fuold = Fu(j,:,r);
       for i=istart(j):M
       %
       % sample new control point i       
       Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu(j,[1:i-1, i+1:end],r)';    
   
       Funew = Fu(j,:,r);
       Funew(i) = Fui;
    
       % sample the remaining points 
       cmu = PropDist.qF{j}.cmuMinus + Funew*PropDist.qF{j}.KInvKu;
       Fnew = gaussianFastSample(1, cmu, PropDist.qF{j}.L);
   
       FFnew = F(:,:,r);
       FFnew(j,:) = Fnew;
    
       % perform an evaluation of the likelihood p(Genes | F) 
       [newLogLik predgen] = logLTFdiffEquation(LikParams, FFnew, r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       % Metropolis-Hastings to accept-reject the proposal
       [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(r,:),2), 0, 0);
    
       %%%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
       %         
         trform = 'exp';%model.Likelihood.TFsingleAct; % 'exp' or 'linear'
         subplot(1,NumOfReplicas,r);
         plot(Xu, feval(trform,Fu(j,:,r)),'or','MarkerSize', 14,'lineWidth', 3);
         hold on;
         plot(TimesF, feval(trform,F(j,:,r)),'g','lineWidth',4);
         if strcmp(model.GroundTr,'yes')==1
         plot(TimesF, feval(trform,model.FF(j,:)),'k','lineWidth',4);
         end
         title(j);
         pause(0.3);
         plot(Xu(i), feval(trform,Fu(j,i,r)),'oy','MarkerSize', 14,'lineWidth', 3); 
         plot(Xu(i), feval(trform,Funew(i)), 'md','MarkerSize', 14, 'lineWidth',3);
         plot(TimesF, feval(trform,FFnew(j,:)), '--b', 'lineWidth', 4); 
         pause(0.5);
         hold off;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%% end of visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       if (it > BurnInIters)
         acceptF(j,i,r) = acceptF(j,i,r) + accept; 
       end
    
       % update protein F
       if accept == 1
         F(j,:,r) = Fnew;
         Fu(j,:,r) = Funew;
         PredictedGenes(:,:,r) = predgen;
         oldLogLik(r,:)  = newLogLik;
       end
       end % num of control points loop
       %
    end % num Genes loop
    end % num Replicas loop
    
   
    % sample new kinetic parameters for each gene separately 
    for j=1:NumOfGenes
        KineticsNew = randn(1,SizKin).*sqrt(PropDist.qKinVars(j,:)) + log(LikParams.kinetics(j,:));
        KineticsNew(KineticsNew<-10) =-10; 
        KineticsNew(KineticsNew>10) = 10;
        KineticsNew = exp(KineticsNew); 
        %
        % set the initial condition to be zero
        % if you know that it should be zero
        if model.Constraints.InitialConds(j) == 0
        KineticsNew(4) = KineticsNew(1)/KineticsNew(2); 
        end
        % this makes the initial condition to be zero
        %KineticsNew(5) = 1;%kinetics(j,5); % dot not sample the gamma
        
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
          %newLogLik(r) = logLTFdiffEquation(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
          newLogLik(r) = logLTFdiffEquation(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
        end
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorKin(j,:),2);
        newP = sum(newLogLik(:)) + sum(LogPriorKinNew(:)); 
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
    
    % sample the interaction weights 
    for j=1:NumOfGenes
        % 
        %  
        if posw == 1
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.qWeigVars(j,:)) + log([LikParams.W(j,:), LikParams.W0(j)]+eps);     
            Wnew = exp(Wnew);
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.qWeigVars(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.Constraints.W(j,:);
       
        LikParams1 = LikParams;
        LikParams1.W(j,:) = Wnew(1:NumOfTFs);
        LikParams1.W0(j)=Wnew(end);
      
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          %newLogLik(r) = logLTFdiffEquation(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
          newLogLik(r) = logLTFdiffEquation(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
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
    % learn binary network connections (this not working for the current version) 
    if netLearn == 1
        for j=1:NumOfGenes
           % Gibbs sampling 
           LikParams1 = LikParams;
           newLogLik = zeros(sComb,NumOfReplicas);
           for co=1:sComb
               LikParams1.Net_X(j,:) = Comb(co,:);
               for r=1:NumOfReplicas
                 newLogLik(co,r) = logLTFdiffEquation(LikParams1, F(:,:,r), r, j, Genes(:,:,r), TimesG, TimesF);
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
                
       %for r=1:NumOfReplicas
       %[oldLogLik1(r,:) predgen] = logLTFdiffEquation(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %PredictedGenes(:,:,r) = predgen;
       %end
       %uu = TimesF; 
       %[commonSlots, comInds] = intersect(uu,TimesG);
       %sumSquerrors
       %sum(sum((Genes -  PredictedGenes(:,comInds)).^2))
       %pause
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       LikParams.sigmas = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       %for r=1:NumOfReplicas
       %oldLogLik(r,:) = logLTFdiffEquation(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %end  
       okk = repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1);
       oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end

    %
    %
    % sample the lenghtscale of the TFs 
    for j=1:NumOfTFs
        %
        % to samples the lengthscale you need to evaluate 
        % the GP pior compute the new and the current/old GP prior 
        newlogEll = randn.*sqrt(PropDist.qLengScVars(j)) + 2*model.GP.logtheta(j,1);
        %newlogEll = model.GP.logtheta(j,1);
        newEll2 = exp(newlogEll);  
        newK = exp(-(0.5/newEll2)*model.GP.X2);  
        % compute the Cholesky decomposition of the new K
        [newL,er]=jitterChol(newK);
        newL = newL';
        % evaluate the new log GP prior value 
        invnewL = newL\eye(SizF+M);
        newLogDetK = 2*sum(log(diag(newL)));
        newlogGP = - 0.5*NumOfReplicas*newLogDetK;
        oldlogGP = - 0.5*NumOfReplicas*PropDist.qF{j}.LogDetK;
        for r=1:NumOfReplicas     
           temp = invnewL*([F(j,:,r), Fu(j,:,r)]'); 
           newlogGP = newlogGP - 0.5*temp'*temp;
           temp = PropDist.qF{j}.invL*([F(j,:,r), Fu(j,:,r)]'); 
           oldlogGP = oldlogGP - 0.5*temp'*temp;
        end
        LogPriorLengScnew = feval(lnpriorLengSc, newlogEll, model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
        % Metropolis-Hastings to accept-reject the proposal
        oldlogGP = oldlogGP + oldLogPriorLengSc(j);
        newlogGP = newlogGP + LogPriorLengScnew; 
        %
        
        [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, 0, 0);
        %%%%%%%%%%%%%%%%  start accept/update proposal for the lengthscale %%%%%%%%%%%% 
        if accept == 1
           model.GP.logtheta(j,1) = newlogEll/2;
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
           for i=1:M
           %  
              G = [1:i-1, i+1:M];  
              [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
           %
           end
           PropDist.qF{j}.alpha = alpha;
           PropDist.qF{j}.ku = ku;
           PropDist.qF{j}.KInvK = KInvK;   
        end
        % 
        if (it > BurnInIters) 
           acceptLengSc(j) = acceptLengSc(j) + accept;
        end
        %%%%%%%%%%%%%%%%%%%%%%% end accept/update proposal for the lengthscale %%%%%%%%%%%%%%%%
        %
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
        if fixsigma2 == 0
        samples.sigma2(:,:,cnt) = LikParams.sigmas(:,:,1);
        end
        if netLearn == 1
            samples.NetX(:,:,cnt) = LikParams.Net_X;
        end
        samples.logthetas(:,:,cnt) = model.GP.logtheta;
        samples.LogL(cnt) = sum(oldLogLik(:));
        %save(trainOps.ResFile,'samples','model');
        %
    end
    %
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.B = LikParams.kinetics(:,1)'; % basal rates
model.Likelihood.D = LikParams.kinetics(:,2)'; % decay rates
model.Likelihood.S = LikParams.kinetics(:,3)'; % sensitivities 
model.Likelihood.A = LikParams.kinetics(:,4)'; % initial conditions
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;
if netLearn == 1
model.Likelihood.Net_x = LikParams.Net_X;
end
model.Likelihood.sigmas = LikParams.sigmas;
model.F = F;
model.Fu = Fu;
accRateF = (acceptF/Iters)*100; 
if istart == 2
accRateF(:,1,:) = 100*ones(size(accRateF(:,1,:)));
end
accRateKin = (acceptKin/Iters)*100;
accRateW = (acceptW/Iters)*100;
accRateLengSc = (acceptLengSc/Iters)*100;


