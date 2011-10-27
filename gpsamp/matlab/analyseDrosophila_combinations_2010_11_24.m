
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


analyseDrosophila_constants;
if (~exist('posteriors')),
  [posteriors, linkMargPosteriors, linkPairPosteriors, priors, M, post_genes] = analyseDrosophila_loadResults;
else
  fprintf('Re-using old posteriors, clear ''posteriors'' to force reload/recompute\n');
end



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
  
  % posterior probability of the TF-pair under the 32 hypotheses case
  [foo, J_Pair32{cnt}] = sort(linkPairPosteriors{1}(:, cnt), 'descend');
    
  % posterior probability of the TF-pair under the 16 hypotheses case
  [foo, J_Pair16{cnt}] = sort(linkPairPosteriors{3}(:, cnt), 'descend');
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  [foo, J_Pair4{cnt}] = sort(linkPairPosteriors{5}(:, cnt), 'descend');
          
  % posterior probability of the TF-pair by combining *Heuristically*
  % the single TFs models 
  [foo, J_PairFromSingleTF{cnt}] = sort(linkPairPosteriors{9}(:, cnt), 'descend');

  % baseline maximum likelihood model
  [foo, J_Pairbase{cnt}] = sort(linkPairPosteriors{7}(:, cnt),'descend');  
  
  % Inferelator 
  [foo, J_PairInfer{cnt}] = sort(linkPairPosteriors{8}(:, cnt),'descend');  
  
  %             
  end
end



J_indiv32 = {};
J_indiv16 = {};
J_indiv2 = {};
J_indbase = {};
J_indInfer = {};
for k=1:numTFs,
%
  % single TF probability computed from 32 models 
  [foo, J_indiv32{k}] = sort(linkMargPosteriors{1}(:, k),'descend');
    
  % single TF probability computed from 16 models                
  [foo, J_indiv16{k}] = sort(linkMargPosteriors{3}(:, k),'descend');
      
  % single TF probability computed from 2 models             
  [foo, J_indiv2{k}] = sort(linkMargPosteriors{5}(:, k),'descend');
  
  % baseline maximum likelihood evaluation         
  [foo, J_indbase{k}] = sort(linkMargPosteriors{7}(:, k), 'descend');
  
   % posterior of the TF-link being active using inferelator   
  [foo, J_indInfer{k}] = sort(linkMargPosteriors{8}(:, k), 'descend');
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

