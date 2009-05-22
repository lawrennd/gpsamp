function [loglikval, PredictedGenes, totalmse] = logLTFdiffEquation(LikParams, F, R, Gindex, Genes, TimesG, TimesF)
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
% Notes: Numerical integration was done using the composite Simpson rule 


Ngenes = size(Gindex,2);

ff = F; 
uu = TimesF;
Delta = uu(2)-uu(1); 
[commonSlots, comInds] = intersect(uu,TimesG);

% compute the joint activation function of the TFs 
% i.e. g(f_1(u),...,f_M(u);w_j) 
fx = TFactivFun(LikParams,ff,Gindex);

% genes are independent given the TF proteins 
for m=1:Ngenes
    j = Gindex(m);
    %
    B = LikParams.kinetics(j,1);
    D = LikParams.kinetics(j,2);
    S = LikParams.kinetics(j,3);   
    A = LikParams.kinetics(j,4);
        
    % Trapezoid rule of numerical integration
    %IntVals = exp(D*uu).*fx;
    %IntVals = Delta*cumtrapz(IntVals);
    %IntVals = IntVals(comInds);
    %IntVals1 = IntVals;
    %Ngenes
    
    % Simpson rule of integration 
    ffx = exp(D*uu).*fx(m,:); 
    IntVals = ffx;
    IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    IntVals = cumsum(IntVals);
    
    %IntVals = (Delta/3)*IntVals(comInds);
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(comInds(1:end-1));
    IntVals = (Delta/3)*IntVals;
    IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
   
    %ok = exp(-TimesG*D); 
    ok = exp(-TimesF*D);
    
    PredictedGenes(m,:) = B/D  + (A - B/D)*ok + S*(ok.*IntVals);

    %loglikval(j) = - 0.5*sum(log(2*pi*vars(j,:)),2)....
    %               - 0.5*sum(((Genes(j,:) -  PredictedGenes(j,:)).^2)./vars(j,:),2);
    loglikval(m) = - 0.5*sum(log(2*pi*LikParams.sigmas(j,:,R)),2)....
                   - 0.5*sum(((Genes(j,:) -  PredictedGenes(m,comInds)).^2)./LikParams.sigmas(j,:,R),2);
               
    totalmse(m) = sum(sum((Genes(j,:) -  PredictedGenes(m,comInds)).^2));            
    %           
    %
end

