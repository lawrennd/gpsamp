function logpdf = lnspikeTruncNormalpdf(x, prior)
%  A mixture model with half truncated Gaussian (x>=mu) and a spike Gaussian 
%

% !!!! This will only work if prior.mu = prior.spikeMu = 0 !!!!


if (min(x) >= prior.mu) &  (min(x) >= prior.spikeMu) 
   g1 = (1/sqrt(2*pi*prior.sigma2))*exp(- (0.5/prior.sigma2).*((x-prior.mu).^2));
   g1 = 2*g1;

   g2 = (1/sqrt(2*pi*prior.spikeSigma2))*exp(- (0.5/prior.spikeSigma2).*((x-prior.spikeMu).^2));

   g2 = 2*g2;
   logpdf = log((1-prior.pis).*g1 + prior.pis.*g2);
else
    disp('x cannot be smaller than zero')
end