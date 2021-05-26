function options = gpsampOptions(gplik) 
%
%

options.kern = 'rbf';
options.constraints.kernHyper = 'free'; % 'free' or 'fixed'

switch gplik 
    case 'regression'
        % Gaussian likelihood 
        options.Likelihood = 'Gaussian';
    case 'classification' % binary classification 
        % Probit likelihood
        options.Likelihood = 'Sigmoid';
    case 'poissonRegr' % binary classification 
        % Probit likelihood
        options.Likelihood = 'Poisson';
    case 'ODE' 
        % not included 
end

options.constraints.likHyper = 'free'; % 'free' or 'fixed'

        
