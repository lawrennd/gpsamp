function logpdf = lninvGammapdf(x,prior)
%
%

a = prior.a;
b = prior.b;

logpdf = (-a-1)*log(x) - b./x + a*log(b) - gammaln(a);