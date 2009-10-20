function logpdf = lnnormalpdf(x,mu,sigma2)
% natural logarithm of the normal density
%
%

logpdf = -0.5*log(2*pi*sigma2) - (0.5/sigma2).*((x-mu).^2);