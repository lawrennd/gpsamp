function out = logLPoisson(lik, Y, F)
%function out = logLPoisson(lik, Y,F)
%
%Description: Poisson GP log-likelihood  
%

out = Y.*F - exp(F) - gammaln(Y+1); 