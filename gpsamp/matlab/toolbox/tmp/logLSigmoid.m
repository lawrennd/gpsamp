function out = logLSigmoid(lik, Y, F)
%
%

YF= Y.*F; 
out = 1./(1 + exp(-YF));