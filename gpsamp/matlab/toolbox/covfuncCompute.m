function Knm = covfuncCompute(logtheta, X, Xu)
%function Knm = covfuncCompute(logtheta, X, Xu)
%
%Description:  It computes the covariance function between 
%              two set of inputs points: X and Xu.  
%
%Supported covariance functions:  ARD kernel,   

[n D] = size(X);
[m,D] = size(Xu);

%logtheta = covfunParams.logtheta;  

sigmaf = exp(2*logtheta(D+1));
X = X ./ repmat(exp(logtheta(1:D))',n,1);
Xu = Xu ./ repmat(exp(logtheta(1:D))',m,1);
%
Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
Knm = sigmaf*exp(-0.5*Knm');


