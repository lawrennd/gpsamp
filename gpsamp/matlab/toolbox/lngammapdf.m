function logpdf = lngammapdf(x,prior)
%
%

%z = x ./ b;
%u = (a - 1) .* log(z) - z - gammaln(a);
%y = exp(u) ./ b;

a = prior.a;
b = prior.b;

logpdf = (a-1)*log(x) - b*x + a*log(b) - gammaln(a);