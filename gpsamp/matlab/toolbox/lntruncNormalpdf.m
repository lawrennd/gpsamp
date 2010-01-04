function logpdf = lntruncNormalpdf(x, prior)
% natural logarithm of half truncated  normal density (x > mu)
%
%

sigma  = sqrt(prior.sigma2);

z = (x-prior.mu)./sigma;

logpdf = -0.5*log(2*pi*prior.sigma2) - 0.5.*(z.^2); 

den = 1 - normcdf(0, prior.mu, sigma); 

logpdf = logpdf - log(den); 