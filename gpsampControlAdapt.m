function [model PropDist samples accRates] = gpsampControlAdapt(model, ops)
% adaptive MCMc to specify the number and location of control points 
% as well as proposal distribution for kernel hyperparameters 
%
%

BurnInIters = ops.Burnin; 
Iters = ops.T; 
% initial number of control variables 
M = ops.initialNumContrPnts;
if M < 3, 
    M = 3; 
end
Init = ops.initialNumContrPnts;
IncM = ops.incrNumContrBy;

[n D] = size(model.X);
perm = randperm(n);
Xu =[];
X = model.X;


% initial input location of the control variables
if strcmp(model.constraints.kernHyper, 'fixed') 
 %
 % keep adding control points till the trace falls under the 10% 
 model.K = kernCompute(model.GP, X); 
   while 1 
      if isempty(Xu) 
         Xu = X(perm(1:M),:);
      else
         Xu = [Xu; X(perm(M+1:M+IncM),:)];
      end
      M = size(Xu,1);
      [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta);
      Xu = reshape(Xu,M,D);
      if f(end) < (0.1*sum(diag(model.K)))
          break;
      end
   end
   %     
else
% 
  Xu = X(perm(1:M),:);
  [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta); 
  Xu = reshape(Xu,M,D);
%
end

M = size(Xu,1); 
Fu = zeros(1, M);

% proposal distribution for F
U = n+1:n+M;  
PropDist.qF.m = zeros(n+M,1);
PropDist.qF.K = kernCompute(model.GP, [X; Xu]); 
L=jitterChol(PropDist.qF.K)';
PropDist.qF.invL = L\eye(n+M); 
PropDist.qF.LogDetK = 2*sum(log(diag(L)));      
% compute the conditional GP prior given the control variables
[cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, 1:n, U);
[L,er]=jitterChol(cSigma);
if er>0, L = real(sqrtm(cSigma)); end
cmu = cmuMinus + Fu*KInvKu;
F = gaussianFastSample(1, cmu, L);
model.F = F;
PropDist.qF.cmuMinus = cmuMinus;
PropDist.qF.cSigma = cSigma;
PropDist.qF.KInvKu = KInvKu;
PropDist.qF.L = L;
% compute all the conditional variances for the control Points
for i=1:M
% 
   G = [1:i-1, i+1:M];
   [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
%
end 
%
PropDist.qF.alpha = alpha;
PropDist.qF.ku = ku;
PropDist.qF.KInvK = KInvK;

% proposal for the kernel hyperparameters
PropDist.kern = 0.2*(1/model.prior.kernParams.b)*ones(1,model.GP.nParams);
if model.Likelihood.nParams > 0
  PropDist.lik = 0.2*(1/model.prior.likParams.b)*ones(1,model.Likelihood.nParams);
end

nextbreak = 0;
epsilon = 0.1;
cnt = 0;
model.F = F;
model.Fu = Fu;
model.Xu = Xu;

% do the adaption  
while 1
    %
    model.Xu = Xu;
    [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, ops);
    
    accRateF = accRates.F;
    accRateKern = accRates.kern;
    accRateLik = accRates.lik;
 
    fprintf(1,'------ ADAPTION STEP #%2d, Number of Control Points %2d ------ \n',cnt+1,M); 
    if ops.disp == 1
       fprintf(1,'Acceptance Rates for GP function\n');
       fprintf(1,' (Rates per control point) \n');       
       disp(accRateF);    
       if ~strcmp(model.constraints.kernHyper, 'fixed')
           fprintf(1,'Acceptance Rates for kernel hyperparameters\n');
           disp(accRateKern);
       end
       if ~strcmp(model.constraints.likHyper, 'fixed') & (model.Likelihood.nParams > 0)
           fprintf(1,'Acceptance Rates for likelihood hyperparameters\n');
           disp(accRateLik);
       end
       fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
    end

    % if you got a descent acceptance rate, then stop
    if ((min(accRateF(:)) > ((0.2/M)*100)) &  (accRateKern>15) & (accRateLik>15) )
        if nextbreak == 1
           disp('END OF ADAPTION: acceptance rates OK');
           break; 
        else
           nextbreak = 1;
        end
    end
    
    cnt = cnt + 1;
    % do not allow more than 80 iterations when you adapt the proposal distribution
    if cnt == 80
        warning('END OF ADAPTION: acceptance rates were not all OK');
        break;
    end

    %  increase the number of control points (add controls in every second iteration)
    if (mod(cnt,2)==0) & (min(accRateF(:)) < ((0.2/M)*100)) 
       %
       %
       Xu = [Xu; X(perm(M+1:M+IncM),:)]; 
       M = M + IncM;
       
       % optimize over control input location  by minimizing the trace term
       [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta);
       Xu = reshape(Xu,M,D);
       if D == 1
           Xu = sort(Xu);
       end
       model.Xu = Xu;
       
       % initialize the new control variables given the current F
       model.Fu = zeros(1,M);
       U = n+1:n+M;
       % update proposal distribution 
       PropDist.qF.m = zeros(n+M,1);
       PropDist.qF.K = kernCompute(model.GP, [X; Xu]);     
       L=jitterChol(PropDist.qF.K)';
       PropDist.qF.invL = L\eye(n+M); 
       PropDist.qF.LogDetK = 2*sum(log(diag(L)));
      
       [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, U, 1:n);
       [L,er]=jitterChol(cSigma);
       if er>0, L = real(sqrtm(cSigma)); end
       cmu = cmuMinus + model.F*KInvKu;
       model.Fu = gaussianFastSample(1, cmu, L);
 
       % compute the conditional GP prior given the control variables
       [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, 1:n, U);
       [L,er]=jitterChol(cSigma);
       if er>0, L = real(sqrtm(cSigma)); end
       PropDist.qF.cmuMinus = cmuMinus; 
       PropDist.qF.cSigma = cSigma;
       PropDist.qF.KInvKu = KInvKu;
       PropDist.qF.L = L;
       clear alpha ku KInvK;
       for i=1:M
       %  
         G = [1:i-1, i+1:M];  
         [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
       %
       end
       PropDist.qF.alpha = alpha;
       PropDist.qF.ku = ku;
       PropDist.qF.KInvK = KInvK; 
    end
   
    
    % adapt likelihood hyperparameters (everything apart from F)
    if ~strcmp(model.constraints.likHyper, 'fixed') & (model.Likelihood.nParams > 0)
       if accRateLik > 35
        % incease the covariance to reduce the acceptance rate
        PropDist.lik = PropDist.lik + epsilon*PropDist.lik;
       end
       if accRateLik < 15
        % decrease the covariance to incease the acceptance rate
        PropDist.lik = PropDist.lik - epsilon*PropDist.lik;    
        %
       end
    end
    
    
    % adapt kernel hyperparameters proposal 
    if ~strcmp(model.constraints.kernHyper, 'fixed') 
       if accRateKern > 35
        % incease the covariance to reduce the acceptance rate
        PropDist.kern = PropDist.kern + epsilon*PropDist.kern;
       end
       if accRateKern < 15
        % decrease the covariance to incease the acceptance rate
        PropDist.kern = PropDist.kern - epsilon*PropDist.kern;    
        %
       end
    end
    %
    %
end
