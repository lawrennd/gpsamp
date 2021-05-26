function Knm = covfuncCompute(GPprior, X, Xu)
%function Knm = covfuncCompute(logtheta, X, Xu)
%
%Description:  It computes the covariance function between 
%              two set of inputs points: X and Xu.  
%

%
%Supported covariance functions:  RBF and ARD kernel with added white noise   

jitter = GPprior.sigma2;

[n D] = size(X);
leng = sqrt(GPprior.lengthScale);
sigma2f = GPprior.sigma2f;
X = X ./ repmat(leng(:)',n,1);

if nargin == 3
   [m,D] = size(Xu);   
   Xu = Xu ./ repmat(leng(:)', m, 1);
   %
   Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
   Knm = sigma2f*exp(-0.5*Knm');
else
   Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   Knm = sigma2f*exp(-0.5*Knm') + jitter*eye(n); 
end
