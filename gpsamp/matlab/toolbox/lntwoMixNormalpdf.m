function logpdf = lntwoMixNormalpdf(x, prior)
% natural logarithm of the normal density
%
%

%logpdf = -0.5*log(2*pi*prior.sigma2) - (0.5/prior.sigma2).*((x-prior.mu).^2);
%ok = exp(logpdf);
g1 = (1/sqrt(2*pi*prior.sigma2))*exp(- (0.5/prior.sigma2).*((x-prior.mu).^2));
g2 = (1/sqrt(2*pi*prior.spikeSigma2))*exp(- (0.5/prior.spikeSigma2).*((x-prior.spikeMu).^2));

%[ok; g1; g2]
%pause

logpdf = log((1-prior.pis).*g1 + prior.pis.*g2); 