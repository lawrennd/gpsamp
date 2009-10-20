function [out0, out] = logLProbit(Y,F)
%
%Description:  Log likelihood of the probit model (cumulative Gaussian)
%
     
YF = Y.*F;
     
out1 = (1+erf(YF/sqrt(2)))/2;
out1 = zeros(size(F));
b = 0.158482605320942; 
c = -1.785873318175113;    
ok = YF>-6; 
out(ok) = log(out1(ok)); 
out(~ok) = -YF(~ok).^2/2 + b*YF(~ok) + c;  

out0 = sum(out(:));