function [model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFAdapt(model, Genes, TimesG, TimesF, AdaptOps)
%[model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFAdapt(model, Genes, TimesG, TimesF, AdaptOps)
%
% Description: Sample the parameters of the Bayesian differential equation 
%              model so that to tune/adapt the proposal distribution (number 
%              of control points, variances of Gaussian proposals etc). 
%                ** This function should always be called before 
%                   gpTFSample. When the gpTFSample will be called the 
%                   proposal distribution (previously adapted)
%                   is kept fixed. **
%
% Inputs: 
%    -- model: the structure that contains the likelihood and GP
%              parameters as well as the priors for all these
%              quantities; write 'help model' for full details. 
%    -- PropDist: a stucture that defines the functional form of the 
%                 proposal distribution. It contains the following fields
%              * qF: contains precomputed variables needed for the 
%                    definition of the poposal distribution (that uses 
%                    control variables) over the GP latent functions 
%              * qKinVars: NumOfGenes x 4 matrix that contains the
%                    variances of the Gaussian proposal distributions
%                    used to sample the logarithm of the kinetic parameters 
%                    (represented in the log space; see modelCreate)
%                    for each gene. These Gaussian proposals have diagonal 
%                    covariance matrices so as each row corresponds to each 
%                    gene, 1 column to the basal rates (B), 2 column to 
%                    the decay rates (D), 3 column to the sensitivities(S)
%                    and 4 column to the initial condition parameter (A). 
%              * qWeigVars: NumOfGenes x M matrix that contains the
%                    variances of the Gaussian proposal distribution
%                    of all the parameters that exist in the activation 
%                    function of the latent GP function (that represent the 
%                    log of the TFs). Each row of the matrix corresponds to 
%                    each gene. The number M of parameters depends on the 
%                    number of TFs used and the functional form of the 
%                    activation function. When the sigmoid activation 
%                    function is used then  M=NumOfTFs+1. When the 
%                    Michaelis-Menten that is valid only for a single 
%                    TF (NumOfTFs=1), then M = 1; the single parameter 
%                    corresponds to the gamma Michaelis-Menten constant
%              * qLengScVars: NumOfTFs x 1 vector that contains all the
%                    variances of the Gaussian proposal distributions used to 
%                    sample the logarithm of the lengthscales of the NumOfTFs
%                    different rbf GP priors (see modelCreate) 
%    -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the
%               gene expressions for all genes, all times and 
%               replicas
%    -- TimesG: The time points where gene expression are evaluated 
%    -- TimesF: The times where the GP function are evaluated
%               TimesF>>TimesG
%    -- trainOps: User defined options about the burn-in and sampling
%               iterations, visualization features and others 
%               **see demos**
%
% Outputs: 
%    -- model: see above or write 'help model' for full details. 
%              The outputed model is updated to contain the parameters 
%              values of the final MCMC iteration
%    -- PropDist: as above. the precomputations in the PropDist.qF field
%              can be different (compared to the ones given at the input) 
%              due the update of the kernel lengthscale that determine 
%              the proposal over the latent GP functions
%    -- samples: the structure that contains the samples. In contains
%              
%              
%
%    -- accRateF: acceptance rates during sampling time (after burn-in) 
%                 for the TFs
%    -- accRateKin: >>  >> for the kinetic parameters 
%    -- accRateW: >>  >> for the interaction weigths and any other 
%              parameters that exist in the activation function of 
%              the TF (e.g. the gamma in the Michaels-Menten activation)  
%    -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%              kernel 

BurnInIters = AdaptOps.Burnin; 
Iters = AdaptOps.T; 

%%% SIZES stuff
[NumOfGenes SizG NumReplicas] = size(Genes);
% number of genes 
NumOfGenes = model.Likelihood.NumOfGenes;
% number of times the gene expression is evaluated
SizG = model.Likelihood.NumOfTimes;
% number of replicas per gene 
NumReplicas = model.Likelihood.NumOfReplicas;
% number of transcription factors
NumOfTFs = model.Likelihood.NumOfTFs;
SizKin = model.Likelihood.NumOfKins_perGene;
SizF = size(TimesF,2);

%initial number of control variables 
M = AdaptOps.InitialNOfControlPnts;

% initialize the transcription factors
F = zeros(NumOfTFs, SizF, NumReplicas);
% initialize the control variables 
Fu = zeros(NumOfTFs, M, NumReplicas);
for r=1:NumReplicas
ook = mean(Genes(:,:,r),1);       
F(:,:,r) = mean(ook(:))*ones(NumOfTFs,SizF);
Fu(:,:,r) = repmat(mean(Genes(:,1:M,r)),[NumOfTFs 1]); 
end



% if the initial values of the TFs is kept fixed, then set the first 
% control point to be zero (-5 even GP function is the log of the TF)
for j=1:NumOfTFs
if model.Constraints.Ft0(j)==0  
   if model.Constraints.Ft0_value == 0
       if strcmp(model.Likelihood.TFsingleAct,'lin') & strcmp(model.Likelihood.TFjointAct,'sigmoid') == 0
           Fu(:,1,:) = zeros(NumOfTFs,NumReplicas);
       elseif strcmp(model.Likelihood.TFsingleAct,'exp')
           Fu(:,1,:) = -5*ones(NumOfTFs,NumReplicas);     
       end
   end        
end
end

%% if the initial values of the TFs is kept fixed,  then set the first 
%% control point to be equal to this value  
%if strcmp(model.InitConds_Constraints.Ft0,'fixed')==1  % 'free' or 'fixed' 
%   if model.InitConds_Constraints.Ft0_value == 0
%       Fu(:,1,:) = -5*ones(NumOfTFs,NumReplicas);       
%   elseif model.InitConds_Constraints.Ft0_value > 0
%       Fu(:,1,:) = log(model.InitConds_Constraints.Ft0_value)*ones(NumOfTFs,NumReplicas);
%   else
%       disp('Negative values for the TFs are not allowed');
%   end        
%end
%Fu(:,1,:) = model.FF(:,1);

% if the initial condition of the differential equations is zero, then set accordingly the 
% corresponding kinetic parameter (initial condition)
for j=1:NumOfGenes
if model.Constraints.InitialConds(j) == 0
   %if model.Constraints.A_value == 0
       model.Likelihood.A = model.Likelihood.B./model.Likelihood.D;
   %end        
end
end

% Initial input locations of the control points
% (placed in a regular grid)
step = (max(TimesF) - min(TimesF))/(M-1);
Xu = TimesF(1):step:TimesF(end);
Xu = Xu(1:M);

% Initial proposal distribution  for the GP latent function
n = SizF;
U = n+1:n+M; 
for j=1:NumOfTFs 
   PropDist.qF{j}.m = zeros(n+M,1);
   PropDist.qF{j}.K = covfuncCompute(model.GP.logtheta(j,:), [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]); 
   %PropDist.qF{1}.K = PropDist.qF{1}.K + 0.1*eye(size(PropDist.qF{1}.K)); 
   L=jitterChol(PropDist.qF{j}.K)';
   PropDist.qF{j}.invL = L\eye(n+M); 
   PropDist.qF{j}.LogDetK = 2*sum(log(diag(L)));
      
   % compute the conditional GP prior given the control variables
   [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, 1:n, U);
   [L,er]=jitterChol(cSigma);
   if er>0, L = real(sqrtm(cSigma)); end
   for r=1:NumReplicas
   cmu = cmuMinus + Fu(j,:,r)*KInvKu;
   F(j,:,r) = gaussianFastSample(1, cmu, L);
   end
   PropDist.qF{j}.cmuMinus = cmuMinus;
   PropDist.qF{j}.cSigma = cSigma;
   PropDist.qF{j}.KInvKu = KInvKu;
   PropDist.qF{j}.L = L;
   % compute all the conditional variances for the control Points
   for i=1:M
   % 
       G = [1:i-1, i+1:M];  
       [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
   %
   end 
   %
   PropDist.qF{j}.alpha = alpha;
   PropDist.qF{j}.ku = ku;
   PropDist.qF{j}.KInvK = KInvK;
end % end NumOfTFs loop


% precompution
X = [TimesF(:); Xu(:)];
model.GP.X2 = -2*X*X' + repmat(sum(X.*X,2)',n+M,1) + repmat(sum(X.*X,2),1,n+M);

% Initial proposal distribution for the kinetic parameters,
% interection weights and the lengthscale of the GPs 
PropDist.qKinVars = ones(NumOfGenes,SizKin);
PropDist.qWeigVars = 0.5*ones(NumOfGenes,NumOfTFs+1);
PropDist.qLengScVars = 0.1*(1/model.prior.GPkernel.lenghtScale.b)*ones(1,NumOfTFs);
% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.0001; qKinAbove = 2;
qWbelow = 0.0001;   qWabove = 2;
qLengScBelow = 0.001*PropDist.qLengScVars(1); 
qLengScAbove = 2*PropDist.qLengScVars(1);
epsilon = 0.1;

model.F = F;
model.Fu = Fu;
cnt = 0;
%
%
% do the adaption 
while 1
%
%  
   model.Xu = Xu;
   [model PropDist samples accRateF accRateKin accRateW accRateLengSc] = gpTFSample(model, PropDist, Genes, TimesG, TimesF, AdaptOps);
   
   if AdaptOps.disp == 1
   fprintf(1,'------ ADAPTION STEP #%2d, Number of Control Points %2d ------ \n',cnt+1,M);    
   fprintf(1,'Acceptance Rates for GP functions\n');
   for rr=1:NumReplicas
   fprintf(1,' * Replica #%2d (rows diffrent TFs, columns control points) \n',rr);       
   disp(accRateF(:,:,rr));
   end    
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Acceptance Rates for lengthscales (per GP function)\n');
   disp(accRateLengSc);
   fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
   end
   fprintf(1,'------------------------------- \n',cnt+1);
   
   % if you got a good acceptance rate, then stop
   if ((min(accRateF(:)) > ((0.2/M)*100)) & (min(accRateKin(:))>15) & (min(accRateW(:))>15) & (min(accRateLengSc(:))>15))
        disp('END OF ADAPTION: acceptance rates OK');
        %pause
        break;
   end
   
   
   cnt = cnt + 1;
   % do not allow more than 50 iterations when you adapt the proposal distribution
   if cnt == 50
       warning('END OF ADAPTION: acceptance rates were not all OK');
       break;
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START ADAPTING CONTROL VARIABLES %%%%%%%%%%%%%%%%
   % adapt the proposal distribution over control variables
   % by adding one control variable  
   if (min(accRateF(:)) < ((0.2/M)*100)) & (mod(cnt,2)==0) 
   M = M+1;
   step = (max(TimesF) - min(TimesF))/(M-1);
   Xu = TimesF(1):step:TimesF(end);
   Xu = Xu(1:M);
   % initialize the control variables given the current F
   model.Fu = zeros(NumOfTFs,M,NumReplicas);
   U = n+1:n+M;
   for j=1:NumOfTFs
      PropDist.qF{j}.m = zeros(n+M,1);
      PropDist.qF{j}.K = covfuncCompute(model.GP.logtheta(j,:), [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]);     
      L=jitterChol(PropDist.qF{j}.K)';
      PropDist.qF{j}.invL = L\eye(n+M); 
      PropDist.qF{j}.LogDetK = 2*sum(log(diag(L)));
      
      [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, U, 1:n);
      [L,er]=jitterChol(cSigma);
      if er>0, L = real(sqrtm(cSigma)); end
      for r=1:NumReplicas
      cmu = cmuMinus + model.F(j,:,r)*KInvKu;
      model.Fu(j,:,r) = gaussianFastSample(1, cmu, L);
      end
      
      % compute the conditional GP prior given the control variables
      [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, 1:n, U);
      [L,er]=jitterChol(cSigma);
      if er>0, L = real(sqrtm(cSigma)); end
      PropDist.qF{j}.cmuMinus = cmuMinus; 
      PropDist.qF{j}.cSigma = cSigma;
      PropDist.qF{j}.KInvKu = KInvKu;
      PropDist.qF{j}.L = L;
      clear alpha ku KInvK;
      for i=1:M
      %  
         G = [1:i-1, i+1:M];  
         [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
      %
      end
      PropDist.qF{j}.alpha = alpha;
      PropDist.qF{j}.ku = ku;
      PropDist.qF{j}.KInvK = KInvK; 
      X = [TimesF(:); Xu(:)];
      model.GP.X2 = -2*X*X' + repmat(sum(X.*X,2)',n+M,1) + repmat(sum(X.*X,2),1,n+M);
   end
   %pause
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT CONTROL VARIABLES %%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.qKinVars(j,:) = PropDist.qKinVars(j,:) + epsilon*PropDist.qKinVars(j,:);
         if PropDist.qKinVars(j,1) > qKinAbove 
             PropDist.qKinVars(j,:) = qKinAbove*ones(1,SizKin);
         end
      end
      if accRateKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.qKinVars(j,:) = PropDist.qKinVars(j,:) - epsilon*PropDist.qKinVars(j,:);    
         if PropDist.qKinVars(j,1) < qKinBelow 
             PropDist.qKinVars(j,:) = qKinBelow*ones(1,SizKin);
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
         PropDist.qWeigVars(j,:) = PropDist.qWeigVars(j,:) + epsilon*PropDist.qWeigVars(j,:);
         if PropDist.qWeigVars(j,1) > qWabove 
             PropDist.qWeigVars(j,:) = qWabove*ones(1,NumOfTFs+1);
         end
      end
      if accRateW(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.qWeigVars(j,:) = PropDist.qWeigVars(j,:) - epsilon*PropDist.qWeigVars(j,:);    
         if PropDist.qWeigVars(j,1) < qWbelow 
             PropDist.qWeigVars(j,:) = qWbelow*ones(1,NumOfTFs+1);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT LENGTHSCALES PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfTFs
      if accRateLengSc(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.qLengScVars(j) = PropDist.qLengScVars(j) + epsilon*PropDist.qLengScVars(j);
         if PropDist.qLengScVars(j) > qLengScAbove
             PropDist.qLengScVars(j) = qLengScAbove;
         end
      end
      if accRateLengSc(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.qLengScVars(j) = PropDist.qLengScVars(j) - epsilon*PropDist.qLengScVars(j);    
         if PropDist.qLengScVars(j) < qLengScBelow 
             PropDist.qLengScVars(j) = qLengScBelow;
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT LENGTHSCALES PROPOSAL %%%%%%%%%%%%%%%%
%
%
end

