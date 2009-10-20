function PredTFs = gpmtfComputeTFODE(LikParams, F, TFindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- TFindex: Indices of the TF-genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- PredTFs: Predicted TF functions at time in TimesF 
%
% Notes: Numerical integration was done using the Trapezoid rule


Ntfs = size(TFindex,2);

% apply single activation for the GP functions  
fx = singleactFunc(LikParams.singleAct,F);

%
uu = LikParams.TimesF;
Delta = uu(2)-uu(1); 
    
for m=1:Ntfs
    j = TFindex(m);
    % 
    D = LikParams.kineticsTF(j,1);
    S = LikParams.kineticsTF(j,2);    
    
    % Trapezoid rule of numerical integration
    IntVals = exp(D*uu).*fx(m,:);
    IntVals = Delta*cumtrapz(IntVals);
    
    % Simpson rule of integration 
    %ffx = exp(D*uu).*fx(m,:); 
    %IntVals = ffx;
    %IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    %IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    %IntVals = cumsum(IntVals);
    %IntVals = (Delta/3)*IntVals;
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
    
    expD = exp(-uu*D);
    PredTFs(m,:) = S*(expD.*IntVals);
    %
end
 