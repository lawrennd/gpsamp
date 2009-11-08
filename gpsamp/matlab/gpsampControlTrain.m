function [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, trainOps)
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
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F = zeros(num_stored, n);
samples.Fu = zeros(num_stored, n);
samples.LogL = zeros(1, num_stored);
X = model.X;
Y = model.y;
F = model.F; 
Fu = model.Fu; % function values  
Xu = model.Xu; % control locations  
M = size(Fu,2);
n = size(X,1);
U = n+1:n+M; 

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:)); 
%if noise == 1
%vars = exp(model.Likelihood.logtheta);
%oldLogLik1 = -0.5*n*log(2*pi*vars) - (0.5/vars)*sum((Y-F(:)).^2);
%elseif noise == 2
%   [pp oldLogLik] = cumGauss(Y,F(:));
%   oldLogLik = sum(oldLogLik(:));
%end


% evaluation of the log prior for the kernel hyperparameters
if strcmp(model.constraints.kernHyper, 'free') 
   lnpriorK = ['ln', model.prior.kernParams.type,'pdf'];
   oldLogPriorK = feval(lnpriorK, model.GP.logtheta, model.prior.kernParams.a, model.prior.kernParams.b);
end

% evaluation of the log prior for the likelihood hyperparameters
if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
   lnpriorLik = ['ln', model.prior.likParams.type,'pdf'];
   oldLogPriorLik = feval(lnpriorLik, model.Likelihood.logtheta, model.prior.likParams.a, model.prior.likParams.b);       
end     


cnt = 0;
acceptF = zeros(1,M);
acceptK = 0;
acceptL = 0;

for it = 1:(BurnInIters + Iters) 
    %
     % sample the control points one-at-a-time 
    Fold = F;
    Fuold =Fu;
    for i=1:M
    %
       % sample the i control point given the rest (Gibbs-like)       
       Fui = randn.*sqrt(PropDist.qF.ku(i)) + PropDist.qF.KInvK(i,:)*Fu([1:i-1, i+1:end])';    
    
       Funew = Fu;
       Funew(i) = Fui;
    
       % resample the whole function
       cmu = PropDist.qF.cmuMinus + Funew*PropDist.qF.KInvKu;
       Fnew = gaussianFastSample(1, cmu, PropDist.qF.L);
    
       % perform an evaluation of the likelihood p(Y | F) 
       newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
       newLogLik = sum(newLogLik(:)); 
       %
       
       % Metropolis-Hastings to accept-reject the proposal
       newP = newLogLik;
       oldP = oldLogLik;
       [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
    
       % visualization
       if trainOps.disp & (mod(it,100) == 0) 
       %[newLogLik oldLogLik]
          if (D==1) 
             trform = 'lin';
             if strcmp(model.Likelihood.type,'Poisson') 
             trform = 'exp';
             end
             plot(X, Y, '+k', 'lineWidth', 2);
             hold on;
             if strcmp(model.Likelihood.type,'Sigmoid') | strcmp(model.Likelihood.type,'Probit')
                 plot(X, zeros(size(X,1)), 'k:')
             end
             plot(X, feval(trform, F), 'g', 'lineWidth', 4);
             %plot(X, F,'or','MarkerSize', 14,'lineWidth', 3);
             pause(0.2);
             plot(Xu, feval(trform, Fu),'or','MarkerSize', 14,'lineWidth', 3);
             plot(Xu(i), feval(trform, Fu(i)),'oy','MarkerSize', 14,'lineWidth', 3);
             %legend(hh,'Current state','Control points');
             set(gca,'FontSize',16);
             plot(Xu(i), feval(trform, Funew(i)), 'md','MarkerSize', 14, 'lineWidth',3);
             plot(X, feval(trform, Fnew), '--b', 'lineWidth', 4); 
             pause(0.5)
             %plot(Xu, Funew, 'md','MarkerSize', 14, 'lineWidth', 3);
             %legend(hh,'Current state','Control points','Proposed state');
             set(gca,'FontSize',16);
             clear hh;
             hold off;
          end
       %
       end
    
       if (it > BurnInIters) 
          acceptF(i) = acceptF(i) + accept;
       end
    
       if accept == 1
          F = Fnew;
          Fu = Funew;
          oldLogLik = newLogLik;
       end
       % 
    end % end control point loop 
    
    
    % sample gp likelihood parameters 
    if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
       % 
       newlogLik = randn.*sqrt(PropDist.lik) + model.Likelihood.logtheta;
       
       Lik1 = model.Likelihood;
       Lik1.logtheta = newlogLik;
       % perform an evaluation of the likelihood p(Y | F) 
       newLogLik = loglikHandle(Lik1, Y, F(:));
       newLogLik = sum(newLogLik(:)); 
       
       newLogPriorLik = feval(lnpriorLik, newlogLik, model.prior.likParams.a, model.prior.likParams.b);
       
       % Metropolis-Hastings to accept-reject the proposal
       oldlogP = oldLogLik + sum(oldLogPriorLik(:));
       newlogP = newLogLik + sum(newLogPriorLik(:)); 
       %
       [accept, uprob] = metropolisHastings(newlogP, oldlogP, 0, 0);
       if accept == 1
        %
          model.Likelihood.logtheta = newlogLik;
          oldLogPriorLik = newLogPriorLik;
          oldLogLik = newLogLik;
        %
       end
        %
       if (it > BurnInIters) 
           acceptL = acceptL + accept;
       end
       %
    end
    
    
    % sample the hyperparameters 
    if strcmp(model.constraints.kernHyper, 'free')
       % 
       newlogK = randn.*sqrt(PropDist.kern) + model.GP.logtheta;
       GPtmp = model.GP; 
       GPtmp.logtheta = newlogK;
       newK = kernCompute(GPtmp, [X; Xu]);    
       [newL,er]=jitterChol(newK);
       newL = newL';
       % evaluate the new log GP prior value 
       invnewL = newL\eye(n+M);
       newLogDetK = 2*sum(log(diag(newL)));
      
       newlogGP = - 0.5*newLogDetK;
       oldlogGP = - 0.5*PropDist.qF.LogDetK;
       temp = invnewL*[F, Fu]'; 
       newlogGP = newlogGP - 0.5*temp'*temp;
       temp = PropDist.qF.invL*[F, Fu]'; 
       oldlogGP = oldlogGP - 0.5*temp'*temp;
       newLogPriorK = feval(lnpriorK, newlogK, model.prior.kernParams.a, model.prior.kernParams.b);
       
       % Metropolis-Hastings to accept-reject the proposal
       oldlogGP = oldlogGP + sum(oldLogPriorK(:));
       newlogGP = newlogGP + sum(newLogPriorK(:)); 
       %
       [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, 0, 0);
       if accept == 1
        %
          model.GP.logtheta = newlogK;
          oldLogPriorK = newLogPriorK;
          % update proposal for F
          PropDist.qF.K = newK;
          PropDist.qF.invL = invnewL; 
          PropDist.qF.LogDetK = newLogDetK;
          [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', newK, 1:n, U);
          [L,er]=jitterChol(cSigma);
          if er>0, L = real(sqrtm(cSigma)); end
          PropDist.qF.cmuMinus = cmuMinus; 
          PropDist.qF.cSigma = cSigma;
          PropDist.qF.KInvKu = KInvKu;
          PropDist.qF.L = L;
          for i=1:M
          %  
             G = [1:i-1, i+1:M];  
             [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
          %
          end
          PropDist.qF.alpha = alpha;
          PropDist.qF.ku = ku;
          PropDist.qF.KInvK = KInvK;  
          %
        end
        % 
        %[accept acceptH oldlogGP]
        
        if (it > BurnInIters) 
           acceptK = acceptK + accept;
        end
        %
    end
  
    % keep samples after burn in 
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
    %
        cnt = cnt + 1;
        samples.F(cnt,:) = F;
        samples.Fu(cnt,:) = F;
        if strcmp(model.constraints.kernHyper, 'free')
           samples.kernLogtheta(cnt,:) = model.GP.logtheta;    
        end 
        if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
           samples.likLogtheta(cnt,:) = model.Likelihood.logtheta;  
        end
        samples.LogL(cnt) = oldLogLik;
    %
    end
    %        
end

%
model.F = F;
model.Fu = Fu;
accRates.F = (acceptF/Iters)*100;
if strcmp(model.constraints.kernHyper, 'free')
   accRates.kern = (acceptK/Iters)*100;
else
   accRates.kern = 100;
end

if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
   accRates.lik = (acceptL/Iters)*100;
else
   accRates.lik = 100;
end
