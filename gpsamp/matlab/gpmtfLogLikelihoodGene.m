function [loglikval, PredGenes] = gpmtfLogLikelihoodGene(LikParams, F, R, Gindex)
%function [loglikval, PredGenes] = gpmtfLogLikelihoodgene(LikParams, F, R, Gindex)
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


PredGenes = gpmtfComputeGeneODE(LikParams, F, R, Gindex);

%
%for m=1:size(Gindex,2)
%    %
%    j = Gindex(m);
%    loglikval(m) = - 0.5*sum(log(2*pi*LikParams.sigmas(j,:,R)),2)....
%                   - 0.5*sum(((Genes(j,:) - PredGenes(m,LikParams.comInds)).^2)./LikParams.sigmas(j,:,R),2);
%               
%    %totalmse(m) = sum(sum((Genes(j,:) -  PredGenes(m,LikParams.comInds)).^2)); 
%    
%    %[Genes(j,:);PredGenes(m,comInds)]
%    %plot(LikParams.TimesG, LikParams.Genes(j,:),'r');
%    %hold on; 
%    %plot(LikParams.TimesF(31:end), PredGenes(m,:));
%    %pause
%    %hold off
%    %
%end


%loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,R)),2)....
%                   - 0.5*sum(((Genes(Gindex,:) - PredGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,R),2);

% REMINDER about what the active variable mean
% active(1) -> puma variances are used
% active(2) -> white noise variances are used (possibly added to the puma variances if also active(1)=1)
% active(3) -> rbf GP is used (possibly added to the puma variances if also active(1)=1)

% if the rbf GP noise model is used then you have to deal with full covariances
if LikParams.noiseModel.active(3) == 1
   % 
   %
   for m=1:size(Gindex(:), 1)
   %
   %
      j = Gindex(m); 
      
      diff = (LikParams.Genes(j,:,R) - PredGenes(m, LikParams.comInds))';
      % note that if they are puma variances, then these have been added to the covariances elsewhere
      % and the whole covariacne matrix is already inverted
      if LikParams.noiseModel.active(1) == 0
         c = LikParams.noiseModel.InvL(:,:,j)*diff;
         logDet = LikParams.noiseModel.LogDetSigma(j);
      else
         c = LikParams.noiseModel.InvL{R}(:,:,j)*diff;  
         logDet = LikParams.noiseModel.LogDetSigma(j,R); 
      end
      
      loglikval(m) = - (0.5*LikParams.numTimes)*log(2*pi) - 0.5*logDet - 0.5*(c'*c);
      %
   end
   %
else
   %
   %
   if LikParams.noiseModel.active(1) == 1
       sigmas = LikParams.noiseModel.pumaSigma2(Gindex, : , R);
   else
       sigmas = zeros(size(Gindex(:),1), LikParams.numTimes);
   end
    
   %
   if LikParams.noiseModel.active(2) == 1
       sigmas = sigmas + repmat(LikParams.noiseModel.sigma2(Gindex)', 1, LikParams.numTimes ); 
   end
    
   loglikval = - 0.5*sum(log(2*pi*sigmas),2)....
               - 0.5*sum(((LikParams.Genes(Gindex,:,R) - PredGenes(:,LikParams.comInds)).^2)./sigmas,2);
             
   %
   %
end

loglikval = loglikval(:)';                 