function PredGenes = gpmtfComputeGeneODE(LikParams, F, R, Gindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- Gindex: Indicets the gens for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 
%
% Notes: Numerical integration was done using the Trapezoid rule 


Ngenes = size(Gindex,2);

% This useful to speed up computation when TF are fixed (therefore precomputed once)
if ~isempty(LikParams.TF) 
   F = LikParams.TF(:,:,R);
else
   F = gpmtfComputeTF(LikParams, F, 1:LikParams.numTFs);
end
%

% compute the joint activation function of the TFs 
% i.e. g(f_1(u-tau_j),...,f_M(u-tau_j);w_j) 
fx = jointactFunc(LikParams,F,Gindex);
%fx = jointactFunc2(LikParams,F,Gindex);

uu = LikParams.TimesF(LikParams.startTime:end);
Delta = uu(2)-uu(1); 

PredGenes = zeros(Ngenes,size(uu,2));

for m=1:Ngenes
    j = Gindex(m);
    %
    B = LikParams.kinetics(j,1);
    D = LikParams.kinetics(j,2);
    S = LikParams.kinetics(j,3);   
    A = LikParams.kinetics(j,4);
  
    % Trapezoid rule of numerical integration
    %IntVals = exp(D*uu).*fx(m,:);
    ffx = exp(D*uu).*fx(m,LikParams.Tausindex(j):LikParams.Tausindex(j)+LikParams.sizTime-1);
    IntVals = zeros(size(ffx));
    IntVals(2:end) = .5 * Delta*cumsum(ffx(1:end-1) + ffx(2:end));
    %IntVals = Delta*cumtrapz(IntVals);
    %IntVals = IntVals(comInds);
    
    % Simpson rule of integration 
    %ffx = exp(D*uu).*fx(m,:); 
    %IntVals = ffx;
    %IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    %IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    %IntVals = cumsum(IntVals);
    %IntVals = (Delta/3)*IntVals;
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
    
    expD = exp(-uu*D);
    PredGenes(m,:) = B/D  + (A - B/D)*expD + S*(expD.*IntVals);
   %
end    
