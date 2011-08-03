

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

load /local/mtitsias/michalis/mlprojects/gpsamp/matlab/datasets/drosophila_data;
load /local/mtitsias/michalis/mlprojects/gpsamp/matlab/datasets/testset;
genesAndChip = importdata('/local/mtitsias/michalis/mlprojects/gpsamp/matlab/datasets/eileen_nature_training_set.txt'); 

if ~isfield(drosexp, 'fbgns'),
  drosexp.fbgns = drosexp.genes;
end

remainder = 1;
modulus = 1;
normST = 0;
positivConstr = 1; 


comb = [0 0 0 0 0;
            1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;
	        1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	        0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	        0 0 1 1 0; 0 0 1 0 1;
	        0 0 0 1 1;
	        1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1;
	        1 0 1 1 0; 1 0 1 0 1;
	        1 0 0 1 1;
	        0 1 1 1 0; 0 1 1 0 1;
	        0 1 0 1 1;
	        0 0 1 1 1;
	        1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;
	        1 1 1 1 1];

testindices = remainder:modulus:length(testset.indices);
indices = testset.indices(testindices);
numGenes = length(indices);
mygenes = drosexp.genes(indices);
Genes = drosexp.fitmean(indices, :);
GenesVar = drosexp.fitvar(indices, :);
fbgns = genesAndChip.textdata(2:end,1);
 
% Normalization
for j=1:size(Genes,1)
    G = Genes(j,:);
    G = G(:);
    sc = mean(G.^2, 1);  
    Genes(j,:) = Genes(j,:)/sqrt(sc); 
    GenesVar(j,:) = GenesVar(j,:)/sc;
end    
Genes = reshape(Genes, numGenes, 12, 3);
GenesVar = reshape(GenesVar,numGenes, 12, 3);
       
% collect the genes for the TFs
if isstruct(drosTF.fbgns),
  fbgnsTF = struct2cell(drosTF.fbgns);
else
  fbgnsTF = drosTF.fbgns; 
end
GenesTF = [];
GenesTFVar = [];
for i=1:size(fbgnsTF(:),1) 
   %
   I = find(strcmp(fbgnsTF(i), drosexp.fbgns));
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);
   % select one probe for the geneTF 
   [val j] = max(prScore);
   GenesTF = [GenesTF; drosexp.fitmean(I(j), :)];
   GenesTFVar = [GenesTFVar; drosexp.fitvar(I(j), :)];
   %
end
numTFs = 5;
for j=1:size(GenesTF,1)
   G = GenesTF(j,:);
   G = G(:);
   sc = mean(G.^2, 1);
   G =GenesTF(j,:)/sqrt(sc); 
   GenesTF(j,:) = G;
   GenesTFVar(j,:) = GenesTFVar(j,:)/sc;
end
GenesTF = reshape(GenesTF, numTFs, 12, 3);
GenesTFVar = reshape(GenesTFVar,numTFs, 12, 3);
  

TimesG = 0:11;
numGenes = size(Genes, 1);

% number of fold in cross-validation
K = 12; 

predErrors1 = zeros(numGenes, size(comb,1)); 
predErrors2 = zeros(numGenes, size(comb,1)); 

vv = zeros(numGenes,1); 

% process each target gene separately
for j=1:numGenes
    TestGene = squeeze(Genes(j,:,:))';
    TestGeneVar = squeeze(GenesVar(j,:,:))';
    
  
    %if mod(j,200) == 0
    %end
    
    vv(j) = var(TestGene(:));
    
    
    while 1
        perm1 = randperm(K); 
        perm2 = randperm(K); 
        perm3 = randperm(K);
        ind12 = find((perm1 - perm2)==0); 
        ind13 = find((perm1 - perm3)==0); 
        ind23 = find((perm2 - perm3)==0); 
        if (isempty(ind12) == 1) &  (isempty(ind13) == 1)  &  (isempty(ind23) == 1) 
            break;
        end
    end
   
    
    % cross validate
    for i=1:K
    %
        TrY = [TestGene(1,perm1([1:(i-1),(i+1):end])), ...
                  TestGene(2,perm2([1:(i-1),(i+1):end])), TestGene(3,perm3([1:(i-1),(i+1):end]))]';
        TsY = [TestGene(1,perm1(i)), TestGene(2,perm2(i)), TestGene(3,perm3(i))]'; 
        
        TrYVar = [TestGeneVar(1,perm1([1:(i-1),(i+1):end])), ...
                      TestGeneVar(2,perm2([1:(i-1),(i+1):end])), TestGeneVar(3,perm3([1:(i-1),(i+1):end]))]';
        TsYVar = [TestGeneVar(1,perm1(i)), TestGeneVar(2,perm2(i)), TestGeneVar(3,perm3(i))]'; 

        TrX = [GenesTF(:,perm1([1:(i-1),(i+1):end]),1), ...
                  GenesTF(:,perm2([1:(i-1),(i+1):end]),2), GenesTF(:,perm3([1:(i-1),(i+1):end]),3)]';     
        TsX = [GenesTF(:,perm1(i),1), GenesTF(:,perm2(i),2), GenesTF(:,perm3(i),3)]';
        

        % Train all possible 2^5 models 
        for c=1:size(comb,1)
          
            % model selection mask
            mask = repmat(comb(c,:), size(TrX,1),1);
            maskTs = repmat(comb(c,:), size(TsX,1),1);
            XX= [ones(11*3,1), TrX.*mask];   
            jit = 1e-7;
            if positivConstr == 1
                H = XX'*XX;
                H = H +  jit*eye(size(H,1));
                f = -XX'*TrY; 
                A = -eye(size(f,1));
                cc = zeros(size(f,1),1);
                % Do not constrain the bias parameter
                cc(1) = 1000;          
                
                b1 = quadprog(H,f,A,cc);
                predErrors1(j,c) = predErrors1(j,c) + sum( (TsY - [ones(3,1),(maskTs.*TsX)]*b1).^2); 
                predErrors2(j,c) = predErrors2(j,c) + sum( (TsY - [ones(3,1),(maskTs.*TsX)]*b1).^2)/var(TsY(:));     
            else
                b1 = ridge(TrY,TrX.*mask,0,0); 
                b2 = ( XX'*(diag(1./TrYVar)*TrY));
                L = chol(XX'*diag(1./TrYVar)*XX +  jit*eye(size(b2,1)));
                b2 = L\(L'\b2);
                % predict the held-out data
                predErrors1(j,c) = predErrors1(j,c) + sum( (TsY - [ones(3,1),(maskTs.*TsX)]*b1).^2);
                predErrors2(j,c) = predErrors2(j,c) + sum( (TsY - [ones(3,1),(maskTs.*TsX)]*b2).^2);
            end
            
        end
    end
end

%save  results/res_linearPositive.mat   predErrors1; 


% standardized mean-squared error 
%N = 3*12; 
%predErrors1norm = predErrors1./(repmat(vv, 1, 32)*N); 
%predErrors2norm = predErrors2./(repmat(vv, 1, 32)*N); 

