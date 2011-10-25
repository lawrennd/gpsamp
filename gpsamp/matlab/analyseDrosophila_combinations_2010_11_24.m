
% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

%tfnames = {'tin', 'bin', 'twi', 'bap', 'Mef2'};
tfnames = {'TIN', 'BIN', 'TWI', 'BAP', 'MEF2'};

% USER-specified: Sizes of the ranking sets 
T = [20, 50, 100, 150, 200];

% USER-specified:: Directory whwre you store figures
ddir = 'figures/';
printPlot = 1; % 0 means not printing

only32 = 1;

fontSize = 7;

% models (as stored  in the results variables; see below) 
% correspidng to 5 TFs being active/inactive 
combConstr = [0 0 0 0 0;
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

combUnconstr = [0 0 0 0 0;
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

baselinecomb = [0 0 0 0 0;
		1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0 ; 0 0 0 1 0; 0 0 0 0 1;
		1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
		0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
		0 0 1 1 0; 0 0 1 0 1;
		0 0 0 1 1];

    
if flag == 1     
    % load the results
    results_b = sortResults(load('results/multitf8c_2010-12-14_summary.mat'));
    %results_b = sortResults(load('results/multitf8b_2010-12-06_summary.mat'));

    % load the 15 baseline models (zerhoth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];
else 
    % load the unconstrained summaries 
    results_b = sortResults(load('results/multitf7a_2010-07-13_summary.mat'));
    
    % !! For Baseline uses the constrained models (if we train the uncosntrained ones, modify the code below) !!
    % load the 15 baseline models (zerhoth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];
    
    % the code below refers to the combConstr variable, so re-define it
    combConstr = combUnconstr; 
end


% You need to exclude the 92 training genes 
genesAndChip = importdata('datasets/eileen_nature_training_set.txt'); 
Trfbgns = genesAndChip.textdata(2:end,1);
mask = ones(size(results_b.genes,1),1); 
for i=1:size(results_b.genes,1)
    for j=1:size(Trfbgns,1)
         if strcmp(results_b.genes(i), Trfbgns(j)) 
             mask(i) = 0; 
         end
    end
end
results_b.marlls = results_b.marlls(mask==1,:);
results_b.genes = results_b.genes(mask==1); 
baseline_a.marlls = baseline_a.marlls(mask==1,:);
baseline_a.genes = baseline_a.genes(mask==1);

numTFs = size(combConstr,2);
M = drosMakeValidationMatrix(chipdistances, results_b.genes, 2000);
% mask out the nans
mask = ones(size(results_b.genes,1),1); 
for i=1:size(results_b.genes,1)
    if isnan(sum( M(i,:) ))
       mask(i) = 0;
    end
end
%regression_a.marlls = regression_a.marlls(mask==1,:);
%regression_a.genes = regression_a.genes(mask==1); 
results_b.marlls = results_b.marlls(mask==1,:);
results_b.genes = results_b.genes(mask==1); 
baseline_a.marlls = baseline_a.marlls(mask==1,:);
baseline_a.genes = baseline_a.genes(mask==1);

% load inferelator results
load inferelator_models.mat; 

singles = singles(mask==1,:);
doubles = doubles(mask==1,:);
M = M(mask==1,:);
numGenes = size(M,1); 



% expand doubles
dd = zeros(numGenes,10);
% twi & Mef2 
dd(:,9) = doubles(:,1);
% bin & bap 
dd(:,6) = doubles(:,2);
% tin & bap 
dd(:,3) = doubles(:,3);
% bap & Mef2
dd(:,10) = doubles(:,4);
% tin & Mef2
dd(:,4) = doubles(:,5);
doubles = dd;
% compute also Antti's doubles  
pairs =  [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2  5; 3 4; 3 5; 4 5];
for i=1:numGenes
   for k=1:10
       doubles(i,k) = max(abs(dd(i,k)), min( abs(singles(i,pairs(k,1))),  abs(singles(i,pairs(k,2))) ) );
   end
end


% indices of all 32 models 
ind32 = 1:size(results_b.marlls, 2); 
% indices that correspond to 16 models (at most 2 TFs) 
ind16 = find(sum(combConstr,2)<=2);
% index of the zeroth model 
ind0 = find(sum(combConstr,2)==0);
% computation of pair-probabilities using set of 32 hypotheses (models)
J_Pair32 = {};
% computation of pair-probabilities using set of 16 hypotheses (models)
J_Pair16 = {};
% computation of pair-probabilities using set of 4 hypotheses (models)
J_Pair4 = {};
% computation of pair-probabilities using by combining heurisrtically 
% the single TF models
J_PairFromSingleTF = {};
cnt = 0;
%
for k=1:numTFs
  for g=(k+1):numTFs
  %   
  cnt = cnt + 1;
  
  % find all the indices inside comb where the TFs "k" and "g" are active 
  indPair = find(combConstr(:,k)==1 & combConstr(:,g)==1);
   
  % posterior probability of the TF-pair under the 32 hypotheses case
  [foo, J_Pair32{cnt}] = sort(logsumexp(results_b.marlls(:,  indPair), 2) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind32, indPair)), 2), 'descend');
    
  % find the *single* index inside comb (restricted to 16 hypotheses) where 
  % the TFs "k" and "g" are active 
  indPairSingle = find(combConstr(:,k)==1 & combConstr(:,g)==1 & sum(combConstr,2)==2 );
 
  % posterior probability of the TF-pair under the 16 hypotheses case
  [foo, J_Pair16{cnt}] = sort(results_b.marlls(:,  indPairSingle) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind16, indPairSingle)), 2), 'descend');
  
  % find the indices inside comb restricted to 4 models 
  %(0 0; k 0; 0 g; k g )
  ind4 = find(  sum(combConstr(:, setdiff(1:numTFs, [k g]) ), 2)==0  );
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  [foo, J_Pair4{cnt}] = sort(results_b.marlls(:, indPairSingle) - ...
  			             logsumexp(results_b.marlls(:, setdiff(ind4, indPairSingle)), 2), 'descend');
          
  % posterior probability of the TF-pair by combining *Heuristically*
  % the single TFs models 
  indk = find(combConstr(:,k)==1  &  sum(combConstr,2)==1);
  indg = find(combConstr(:,g)==1  &  sum(combConstr,2)==1);
  
  [foo, J_PairFromSingleTF{cnt}] = sort(results_b.marlls(:, indk) + results_b.marlls(:, indg)  - ...
  			             results_b.marlls(:, ind0) - logsumexp(results_b.marlls(:, [ind0  indk indg]),2), 'descend');        

  % baseline maximum likelihood model
  indPSinglebase = find(baselinecomb(:,k)==1 & baselinecomb(:,g)==1 & sum(baselinecomb,2)==2 ); 
  [foo, J_Pairbase{cnt}] = sort(baseline_a.marlls(:, indPSinglebase) - baseline_a.marlls(:, 1),'descend');  
  
  % Inferelator 
  [foo, J_PairInfer{cnt}] = sort(dd(:,cnt),'descend');  
  
  %             
  end
end


% computation of three-link probabilities using set of 32 hypotheses (models)
J_Triple32 = {};
% computation of three link probabilities using set of 8 hypotheses (models)
J_Triple8 = {};
J_TripleFromSingleTF = {};
J_Triplebase = {};
cnt = 0;
%
for k=1:numTFs
  for g=(k+1):numTFs
  for f=(g+1):numTFs    
  %   
  cnt = cnt + 1;
  
  % find all the indices inside comb where the TFs "k", "g" and "f" are active 
  indTriple = find(combConstr(:,k)==1 & combConstr(:,g)==1 & combConstr(:,f)==1);
   
  % posterior probability of the TF-pair under the 32 hypotheses case
  [foo, J_Triple32{cnt}] = sort(logsumexp(results_b.marlls(:,  indTriple), 2) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind32, indTriple)), 2), 'descend');
      
  % find the indices inside comb restricted to 8 models 
  %(combination of k,g,f )
  ind8 = find(  sum(combConstr(:, setdiff(1:numTFs, [k g f]) ), 2)==0  );
  
  % find the *single* index inside comb (restricted to 8 hypotheses) where 
  % the TFs "k", "g" and "f"  are active 
  indTripleSingle = find(combConstr(:,k)==1 & combConstr(:,g)==1 & combConstr(:,f)==1 & sum(combConstr,2)==3 );
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  [foo, J_Triple8{cnt}] = sort(results_b.marlls(:, indTripleSingle) - ...
  			             logsumexp(results_b.marlls(:, setdiff(ind8, indTripleSingle)), 2), 'descend');
        
  % posterior probability of the TF-pair by combining *Heuristically*
  % the single TFs models 
  indk = find(combConstr(:,k)==1  &  sum(combConstr,2)==1);
  indg = find(combConstr(:,g)==1  &  sum(combConstr,2)==1);
  indf = find(combConstr(:,g)==1  &  sum(combConstr,2)==1);
  tkg = results_b.marlls(:, indk) + results_b.marlls(:, indg); 
  tkf = results_b.marlls(:, indk) + results_b.marlls(:, indf); 
  tgf = results_b.marlls(:, indg) + results_b.marlls(:, indf); 
  tk0 = results_b.marlls(:, indk) + results_b.marlls(:, ind0); 
  tg0 = results_b.marlls(:, indg) + results_b.marlls(:, ind0); 
  tf0 = results_b.marlls(:, indf) + results_b.marlls(:, ind0); 
  t00 = results_b.marlls(:, ind0) + results_b.marlls(:, ind0); 
  [foo, J_TripleFromSingleTF{cnt}] = sort(results_b.marlls(:, indk) + results_b.marlls(:, indg) + results_b.marlls(:, indf) - ...
  			             results_b.marlls(:, ind0) - ...
                         logsumexp([tkg tkf tgf tk0 tg0 tf0 t00],2), 'descend');        
                         
  % baseline maximum likelihood model
  %indPSinglebase = find(baselinecomb(:,k)==1 & baselinecomb(:,g)==1 & baselinecomb(:,f)==1 & sum(baselinecomb,2)==3 ); 
  %[foo, J_Triplebase{cnt}] = sort(baseline_a.marlls(:, indPSinglebase) - baseline_a.marlls(:, 1),'descend');  
  %             
  end
  end
end



J_indiv32 = {};
J_indiv16 = {};
J_indiv2 = {};
J_indbase = {};
J_indInfer = {};
comb16 = combConstr(ind16,:);
for k=1:numTFs,
%
  indSingle = find(combConstr(:,k)==1);
  % single TF probability computed from 32 models 
  [foo, J_indiv32{k}] = sort(logsumexp(results_b.marlls(:,indSingle),2) - ...
			            logsumexp(results_b.marlls(:,setdiff(ind32, indSingle)),2),'descend');
    
  % single TF probability computed from 16 models                
  [foo, J_indiv16{k}] = sort(logsumexp(results_b.marlls(:,comb16(:, k)==1),2) - ...
			            logsumexp(results_b.marlls(:,comb16(:, k)==0),2),'descend');
      
  % find the *single* index inside comb where only "k" is active 
  indIndiv = find(combConstr(:,k)==1 & sum(combConstr,2)==1 );     
  
  % single TF probability computed from 2 models             
  [foo, J_indiv2{k}] = sort(results_b.marlls(:, indIndiv) - results_b.marlls(:, ind0),'descend');
  
  % baseline maximum likelihood evaluation         
  [foo, J_indbase{k}] = sort(baseline_a.marlls(:, k+1) - baseline_a.marlls(:, 1), 'descend');
  
   % posterior of the TF-link being active using inferelator   
  [foo, J_indInfer{k}] = sort(abs(singles(:,k)), 'descend');   
%  
end


% 1 PLOT ---------------------------------------------------------------
% K TOP LISTS OF POSITIVE PREDICTION OF SINGLE  LINKS
% START  ---------------------------------------------------------------
h1 = figure; 
set(gca, 'FontSize', fontSize);
for k=1:numTFs,
  figure(h1);
  subplot(2, 3, k);
  set(gca, 'FontSize', fontSize);
  if only32 == 0
    drosPlotAccuracyBars({J_indiv32{k}, J_indiv16{k}, J_indiv2{k}, J_indbase{k}}, M(:, k), T);
  else
    drosPlotAccuracyBars({J_indiv32{k},  J_indiv2{k}, J_indbase{k}, J_indInfer{k}}, M(:, k), T);  
  end   
  title(sprintf('%s', tfnames{k}));
end
subplot(2, 3, 6);
set(gca, 'FontSize', fontSize);
if only32 == 0
bar(rand(4));
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-32', 'Posterior-16',  'Posterior-2', 'ML-Baseline', 'Inferelator', 'Random');
else
bar(rand(4));
hold on
plot([0 1], [0 1], 'k--')
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-32', 'Posterior-2', 'ML-Baseline', 'Inferelator', 'Random');  
end
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 18 9])
% 1 PLOT ---------------------------------------------------------------
% K TOP LISTS OF POSITIVE PREDICTION OF SINGLE  LINKS
% END    ---------------------------------------------------------------


% 2 PLOT ---------------------------------------------------------------
% K TOP LISTS OF POSITIVE PREDICTION OF PAIR-TF-LINKS
% START  ---------------------------------------------------------------
h2 = figure; 
set(gca, 'FontSize', fontSize);
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
  cnt = cnt + 1;    
  figure(h2);
  subplot(4, 3, cnt);
  set(gca, 'FontSize', fontSize);
  %drosPlotAccuracyBars({J_joint32{cnt}, J_joint16{cnt}, J_joint4{cnt},  J_joint2{cnt}, J_jointbase{cnt}}, prod(M(:, [k g]), 2), T);
  if only32 == 0
  drosPlotAccuracyBars({J_Pair32{cnt}, J_Pair16{cnt}, J_Pair4{cnt},  J_PairFromSingleTF{cnt}, J_Pairbase{cnt}}, prod(M(:, [k g]), 2), T);
  else
  drosPlotAccuracyBars({J_Pair32{cnt}, J_Pair4{cnt},  J_PairFromSingleTF{cnt}, J_Pairbase{cnt}, J_PairInfer{cnt}}, prod(M(:, [k g]), 2), T);   
  end
  title(sprintf('%s & %s', tfnames{k}, tfnames{g}));
  end
end
figure(h2);
subplot(4, 3, 12);
set(gca, 'FontSize', fontSize);
if only32 == 0
bar(rand(5));
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-32', 'Posterior-16', 'Posterior-4', 'Posterior from single-TF models', 'ML-Baseline', 'Inferelator', 'Random');
else
bar(rand(5));
hold on
plot([0 1], [0 1], 'k--')
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-32', 'Posterior-4', 'Posterior from single-TF models', 'ML-Baseline', 'Inferelator', 'Random');   
end
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 18 12])
% 2 PLOT ---------------------------------------------------------------
% K TOP LISTS OF POSITIVE PREDICTION OF PAIR-TF-LINKS
% END    ---------------------------------------------------------------


% % 3 PLOT ---------------------------------------------------------------
% % K TOP LISTS OF POSITIVE PREDICTION OF TRIPLE-TF-LINKS
% % START  ---------------------------------------------------------------
% h3 = figure; 
% set(gca, 'FontSize', fontSize);
% cnt = 0;
% for k=1:numTFs
%   for g=(k+1):numTFs
%   for f=(g+1):numTFs  
%   cnt = cnt + 1;    
%   figure(h3);
%   subplot(3, 5, cnt);
%   set(gca, 'FontSize', fontSize);
%   drosPlotAccuracyBars({J_Triple32{cnt}, J_Triple8{cnt}, J_TripleFromSingleTF{cnt}}, prod(M(:, [k g f]), 2), T);
%   title(sprintf('%s + %s + %s', tfnames{k}, tfnames{g}, tfnames{f}));
%   end
%   end
% end
% figure(h3);
% subplot(3, 5, cnt+4);
% set(gca, 'FontSize', fontSize);
% bar(rand(3));
% axis([-10 -9 -10 -9]);
% axis off;
% legend('Posterior-32', 'Posterior-8', 'Posterior from single-TF models');
% set(gcf, 'PaperUnits', 'centimeters')
% set(gcf, 'PaperPosition', [0, 0, 15, 15])
% % 3 PLOT ---------------------------------------------------------------
% % K TOP LISTS OF POSITIVE PREDICTION OF TRIPLE-TF-LINKS
% % START  ---------------------------------------------------------------


% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end
dd = date;
if printPlot 
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'PairOfTFs_only32_', num2str(only32) property '.eps']);
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'SingleTFs_only32_', num2str(only32) property '.eps']); 
   %print(h3, '-depsc2', [ddir 'drosophilaBars_' 'TripleTFs_only32_', num2str(only32) property '.eps']); 
end

