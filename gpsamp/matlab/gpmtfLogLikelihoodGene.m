function [loglikval, PredGenes, totalmse] = gpmtfLogLikelihoodGene(LikParams, F, R, Gindex)
%function [loglikval, PredGenes, totalmse] = gpmtfLogLikelihoodgene(LikParams, F, R, Gindex)
%
% Description: Computes the log likelihood of the linear differential
%              equation model with multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R (from whihc TFs are constructed) 
%        -- R: the replica for which you compute the log likelihood  
%        -- Gindex: Indices of the genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredictedGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 


Genes = LikParams.Genes(:,:,R);
PredGenes = gpmtfComputeGeneODE(LikParams, F, Gindex);

uu = LikParams.TimesF(LikParams.startTime:end);
[commonSlots, comInds] = intersect(uu,LikParams.TimesG);

%
for m=1:size(Gindex,2)
    %
    j = Gindex(m);
    loglikval(m) = - 0.5*sum(log(2*pi*LikParams.sigmas(j,:,R)),2)....
                   - 0.5*sum(((Genes(j,:) - PredGenes(m,comInds)).^2)./LikParams.sigmas(j,:,R),2);
               
    totalmse(m) = sum(sum((Genes(j,:) -  PredGenes(m,comInds)).^2)); 
    
    %[Genes(j,:);PredGenes(m,comInds)]
    %plot(LikParams.TimesG, LikParams.Genes(j,:),'r');
    %hold on; 
    %plot(LikParams.TimesF(31:end), PredGenes(m,:));
    %pause
    %hold off
    %
end
%totalmse
%pause