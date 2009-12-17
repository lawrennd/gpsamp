function logpdf = lnnormalpdf(x, prior)
% natural logarithm of the normal density
%
%

logpdf = -0.5*log(2*pi*prior.sigma2) - (0.5/prior.sigma2).*((x-prior.mu).^2);