function PredTFs = gpmtfComputeTF(LikParams, F, TFindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- TFindex: Indices of the TF-genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- PredTFs: Predicted TF functions at time in TimesF 
%



if ~isfield(LikParams,'GenesTF')
    %
    % apply single activation for the GP functions  
    PredTFs = singleactFunc(LikParams.singleAct,F);
    %
else
    %
    % Translational ODE model for the TFs  
    PredTFs = gpmtfComputeTFODE(LikParams, F, TFindex);
   %
end