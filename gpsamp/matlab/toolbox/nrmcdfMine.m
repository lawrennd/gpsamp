function p = nrmcdfMine(x,mu,sigma,pcov,alpha)
%
%


z = (x-mu) ./ sigma;

p = 0.5 * erfc(-z ./ sqrt(2));
%p = 0.5 * erf(-z ./ sqrt(2));