function options = gpmtfOptions(Genes,numTFs) 
%
%
%


NumOfGenes = size(Genes,1);

if numTFs == 1
        %
        % NUMBER of TFs
        options.numTFs = 1;
        %
        % INDIVIDUAL transform of the GP function
        % --- either  'logOnePlusExp' or 'loglogOnePlusExp' or 'lin' or 'exp'
        % --- 'exp' or 'logOnePlusExp' must combine with TFjointAct= 'lin' or 'michMenten'
        % --- 'loglogOnePlusExp' or 'lin' with TFjointAct='sigmoid'
        options.singleAct = 'logOnePlusExp';
        %
        % - JOINT transform by passing through a sigmoid type of function.
        options.jointAct = 'michMenten';
        
        options.Net_X = ones(NumOfGenes,1);
        %        
        % Make the output BINARY at the end or not 
        options.jointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized
        %                         
        % maximum positve value for the gene DELAYS
        % -- 0 means that no delays are allowed 
        options.tauMax = 0;
        %
        % CONSTRAINTS of the initial value of the TF at time t=0. 
        % if 0, then the TF has concentration 0 at time t=0
        % if 1, then the TF is free to take any value at time t=0
        options.constraints.Ft0 = 0; 
        %
        % CONSTRAINTS of the initial conditions of each differential equation 
        % for each gene  
        % if 0, then the ODE initial cond for the jth gene is 0 at time t=0
        % if 1, then the ODE initial cond is free to take any value at time t=0
        options.constraints.initialConds = ones(1,NumOfGenes); 
        options.constraints.geneTFsensitivity = 0;
        %
        %
        options.constraints.X = ones(NumOfGenes,1);
        options.constraints.replicas = 'free'; % 'free' or 'coupled'
        %
else
        %
        % NUMBER of TFs
        options.numTFs = numTFs;
        %
        % Define the activation function that transforms the GP functions. 
        % For multiple TFs the GP functions are linearly mixed and then are passed 
        % through a sigmoid function
        % INDIVIDUAL transform of the GP function
        % --- either  'loglogOnePlusExp' or 'lin' or 'exp'
        % --- 'exp' must combined with TFjointAct= 'lin'
        % --- the rest with TFjointAct='sigmoid'
        options.singleAct = 'logOnePlusExp';
        %
        % - JOINT transform by passing through a sigmoid type of function.
        options.jointAct = 'genHill';
        
        % Make the output BINARY at the end or not 
        options.jointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized
        %                         
        % maximum positve value for the gene DELAYS
        % -- 0 means that no delays are allowed 
        options.tauMax = 0;
        
        % use or not the two mixture (one spike and one borad) for the interaction weights 
        options.spikePriorW = 'no'; 
        %
        % CONSTRAINTS of the initial value of the TF at time t=0. 
        % if 0, then the TF has concentration 0 at time t=0
        % if 1, then the TF is free to take any value at time t=0
        options.constraints.Ft0 = zeros(1,options.numTFs); 
        %
        % CONSTRAINTS of the initial conditions of each differential equation 
        % for each gene  
        % if 0, then the ODE initial cond for the jth gene is 0 at time t=0
        % if 1, then the ODE initial cond is free to take any value at time t=0
        options.constraints.initialConds = ones(1,NumOfGenes); 
        %
        options.constraints.geneTFsensitivity = zeros(1,numTFs);
        %
        % CONSTRAINTS coming from side information about which TFs do not regulate 
        % certain genes. This means that certain values in the interaction 
        % matrix W are constrained to be zero. 
        % if w(i,j)=1 (for i gene and j TF) then the interaction is 
        % allowed and the weight w(i,j) is learned
        % if w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
        % constrained to be zero. 
        options.constraints.X = ones(NumOfGenes,numTFs);
        options.constraints.replicas = 'free'; % 'free' or 'coupled'
        %
end
