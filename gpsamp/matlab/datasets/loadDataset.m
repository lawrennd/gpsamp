function [Genes, TimesG, GeneVars, TFGenes, TFGeneVars, GroundTruth] = loadDataset(dname)
%
% Load the gene expression dataset 

GroundTruth = [];
switch dname 
    case 'toy3TfsDelays' 
        %
        load ToyData3TFDelays04_Sep_09;
        Genes = Y;
        GeneVars = LikParams.sigmas;
        TFGenes = [];
        TFGeneVars = [];
        GroundTruth.F = F;
        GroundTruth.kinetics = Kinetics;
        GroundTruth.W = W;
        GroundTruth.W0 = W0;
        GroundTruth.X = Net_X;
        GroundTruth.Taus = LikParams.Taus; 
        GroundTruth.sigmas = LikParams.sigmas;  
    case 'toy3TFDelaysTFGenes' 
        %
        load ToyData3TFDelays18_Oct_09; 
        Genes = LikParams.Genes;
        TimesG = LikParams.TimesG;
        GeneVars = LikParams.sigmas;
        TFGenes = LikParams.GenesTF;
        TFGeneVars = LikParams.sigmasTF; 
        
        % GP fucntion
        GroundTruth.F = LikParams.F;
        % TF function 
        GroundTruth.TF = LikParams.TF;
        GroundTruth.kinetics = LikParams.kinetics;
        GroundTruth.kineticsTF = LikParams.kineticsTF;
        GroundTruth.W = LikParams.W;
        GroundTruth.W0 = LikParams.W0;
        GroundTruth.X = LikParams.Net_X;
        GroundTruth.Taus = LikParams.Taus; 
        GroundTruth.Tausindex = LikParams.Tausindex;
        GroundTruth.sigmas = LikParams.sigmas;
        GroundTruth.sigmasTF = LikParams.sigmasTF; 
        %
    case 'p53Barenco5Genes'
        %
        load Barencodata;
        TimesG = [0 2 4 6 8 10 12];
        Genes(:,:,1) = y{1}'; 
        Genes(:,:,2) = y{2}'; 
        Genes(:,:,3) = y{3}';
        GeneVars(:,:,1) = yvar{1}';
        GeneVars(:,:,2) = yvar{2}';
        GeneVars(:,:,3) = yvar{3}';
        TFGenes = [];
        TFGeneVars = [];
        %
    case 'ecoliRogersetal14Genes'
        %
        load ecoliNormalisedData; 
        TimesG = times';
        Genes(:,:,1) = y{1}'; 
        TimesG = times'; 
        TFGenes = [];
        TFGeneVars = []; 
        GroundTruth =[];
        %
    case 'p53Barenco50Genes'
        %
        load dataBarenco50Genes;
        %Genes = Genes; 
        %TimesG = TimesG; 
        GeneVars = GenesVars;
        TFGenes = [];
        TFGeneVars = []; 
        GroundTruth =[];
    case 'drosAntti'
        %
        load dros_multitf_example_data;
        Genes = sampledata.mean;
        Genes = reshape(Genes,8,12,3);
        GeneVars = sampledata.var;
        GeneVars = reshape(GeneVars,8,12,3);
        TimesG = 0:(NumOfTimes-1);
        TFGenes = [];
        TFGeneVars = []; 
        GroundTruth =[];
        %
end

[NumOfGenes NumOfTimes NumOfReplicas] = size(Genes);
if (size(TimesG,2) ~= NumOfTimes) 
    disp('Error: the number of time points must be consistent with the Gene expressions.');
end


