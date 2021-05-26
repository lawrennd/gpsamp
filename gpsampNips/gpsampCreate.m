function model = gpsampCreate(y, X, options)
%
%

model.type = 'gpmodel';
model.Likelihood.type = options.Likelihood;

[numData D] = size(X);

model.y = y; 
model.X = X;
model.numData = numData;
model.D = D;

switch model.Likelihood.type
    case 'Gaussian' % standard regression
         %
         model.Likelihood.nParams = 1; % parameters (excluding gp function F)
         model.Likelihood.logtheta = log(0.05); % log(sigma2) 
         %
    case 'Probit'  % binary classifcation
         %
         model.Likelihood.nParams = 0;
         %
    case 'Sigmoid' % binary classification
         %
         model.Likelihood.nParams = 0;
         %
    case 'Poisson' % for counts data      
         %
         model.Likelihood.nParams = 0;  
         %
    case 'ODE'
         % %%
end     

model.constraints.kernHyper = options.constraints.kernHyper;
model.constraints.likHyper = options.constraints.likHyper;

model.GP.type = 'rbf';
% kernel hyperparameters 
model.GP.logtheta = [2*log((max(X) - min(X))*0.2) 0];
model.GP.nParams = D+1;

% prior over the likelihood parameters  
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.likParams.type = 'normal';
model.prior.likParams.constraint = 'positive';
model.prior.likParams.priorSpace = 'log';
model.prior.likParams.a = 0; % mean 
model.prior.likParams.b = 2; % variance
 

% prior over GP kernel hyperparameters 
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.kernParams.type = 'normal';
model.prior.kernParams.constraint = 'positive';
model.prior.kernParams.priorSpace = 'log';
model.prior.kernParams.a = 0; % mean 
model.prior.kernParams.b = 2; % variance
 
% GP latent function values needed to define the likelihood
model.F = zeros(1,model.numData);
