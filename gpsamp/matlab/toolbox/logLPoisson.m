function out = logLPoisson(lik, Y, F)
%function out = logLPoisson(lik, Y,F)
%
%Description: Poisson GP log-likelihood useful for counts data 
%

out = Y(:).*F(:) - exp(F(:)) - gammaln(Y(:)+1); 