function [loglikvalTF, PredGenesTF, totalmse] = gpmtfLogLikelihoodGeneTF(LikParams, F, R, TFindex)
% function [loglikval, PredictedGenes, totalmse] = logLTFdiffEquation(LikParams, F, R, Gindex, Genes, TimesG, TimesF)
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

uu = LikParams.TimesF;
[commonSlots, comInds] = intersect(uu,LikParams.TimesG);

for m=1:size(TFindex,2)
    %
    j = TFindex(m);
    loglikvalTF(m) = - 0.5*sum(log(2*pi*LikParams.sigmasTF(j,:,R)),2)....
                   - 0.5*sum(((GenesTF(j,:) - PredGenesTF(m,comInds)).^2)./LikParams.sigmasTF(j,:,R),2);
               
    totalmse(m) = sum(sum((GenesTF(j,:) -  PredGenesTF(m,comInds)).^2)); 
    
    %[GenesTF(j,:);PredGenesTF(m,comInds)]
    %plot(LikParams.TimesG, LikParams.GenesTF(j,:),'r');
    %hold on; 
    %plot(LikParams.TimesF, PredGenesTF(m,:));
    %disp('in TF')
    %pause
    %hold off
    %
end
