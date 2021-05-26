function Knm = covfuncCompute(GPprior, X, Xu)
%function Knm = covfuncCompute(logtheta, X, Xu)
%
%Description:  It computes the covariance function between 
%              two set of inputs points: X and Xu.  
%

%
%Supported covariance functions:  RBF and ARD kernel.   

jitter = exp(2*GPprior.logtheta(end));

[n D] = size(X);
logtheta = GPprior.logtheta(:);
sigmaf = exp(logtheta(D+1));
X = X ./ repmat(exp(logtheta(1:D))',n,1);

if nargin == 3
   [m,D] = size(Xu);   
   Xu = Xu ./ repmat(exp(logtheta(1:D))',m,1);
   %
   Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
   Knm = sigmaf*exp(-0.5*Knm');
else
   Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   Knm = sigmaf*exp(-0.5*Knm') + jitter*eye(n); 
end
