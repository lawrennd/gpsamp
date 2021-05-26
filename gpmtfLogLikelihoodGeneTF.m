function [loglikvalTF, PredGenesTF] = gpmtfLogLikelihoodGeneTF(LikParams, F, R, TFindex)
% function [loglikval, PredictedGenes] = logLTFdiffEquation(LikParams, F, R, Gindex, Genes, TimesG, TimesF)
%
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- Gindex: Indicets the gens for which you compute the log
%           likelihood      
%        -- Genes: Expression of all gene for the replica R
%        -- TimesG: Times you have gene expression measurements 
%        -- TimesF: Times you you have discretize the TFs 
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredictedGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 
%


GenesTF = LikParams.GenesTF(:,:,R);
PredGenesTF = singleactFunc(LikParams.singleAct, F(TFindex,:));


if LikParams.noiseModel.active(1) == 1
    sigmas = LikParams.noiseModel.pumaSigma2_TF(TFindex, : , R);
else
    sigmas = zeros(size(TFindex(:),1), LikParams.numTimes);
end
    
%
if LikParams.noiseModel.active(2) == 1
    sigmas = sigmas + repmat(LikParams.noiseModel.sigma2_TF(TFindex)', 1, LikParams.numTimes ); 
end

loglikvalTF = - 0.5*sum(log(2*pi*sigmas),2) ...
               - 0.5*sum(((GenesTF(TFindex,:) - PredGenesTF(:,LikParams.comIndsTF)).^2)./sigmas,2);

loglikvalTF = loglikvalTF';
%sum(abs(loglikvalTF - loglikvalTF1))
%pause
