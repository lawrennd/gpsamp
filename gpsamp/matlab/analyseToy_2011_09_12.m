%results{1} = sortResults(load('results/multitfToyOneCond6a_2010-07-09_summary.mat'));
results{1} = sortResults(load('results/multitfToyOneCond7a_2011-06-27_summary.mat'));
%results{2} = sortResults(load('results/multitfToyTwoConds6a_2010-07-09_summary.mat'));
results{2} = sortResults(load('results/multitfToyOtherCond7a_2011-09-12_summary.mat'));
results{3} = sortResults(load('results/multitfToyTwoConds7a_2011-06-27_summary.mat'));

load datasets/toy4TFs28_June_10.mat

tcomb = [0 0 0; 1 0 0; 0 1 0; 0 0 1;
	 1 1 0; 1 0 1; 0 1 1;
	 1 1 1];
gstyles = {'m--', 'g--'};
styles = {'r', 'b', 'g', 'r-.', 'b-.'};
styles2 = {'ro-.', 'bx-', 'gv-'};
LinWidth = 3; 
FONTSIZE=12;
% USER-specified:: Directory where you store figures
ddir = 'figures/';
printPlot = 1; % 0 means not printing
plotAll = 0;

NUMCONDS = 3;

%TFnames = {'first', 'second', 'third'};
TFnames = {'ANT', 'BEE', 'CAR'};
T = [20, 50, 100, 150, 200 500];
numTFs = size(tcomb,2);
% indices of all 32 models 
ind8 = 1:size(results{1}.marlls, 2); 
% index of the zeroth model 
ind0 = find(sum(tcomb,2)==0);


% 1 PLOT ---------------------------------------------------------------
% ROC CURVES FOR GLOBAL PREDICTION OF SINGLE  LINKS
% START  ---------------------------------------------------------------
auc = {};
M = Net(results{1}.genes, 1:3);
linkprobs8 = {};
linkprobs2 = {};
I8 = {};
I2 = {};
val8 = {};
val2 = {};
ind0 = 1;
h1 = figure;
for k=1:NUMCONDS, % two conditions 
  for l=1:3,
    % posterior link probabilites using all possible 8 models  
    linkprobs8{k}(:, l) = logsumexp(results{k}.marlls(:, tcomb(:, l)==1), 2) - ...
	logsumexp(results{k}.marlls(:, tcomb(:, l)==0), 2);
    % posterior link probabilites using 2 models  
    % find the *single* index inside comb where only "k" is active 
    indIndiv = find(tcomb(:,l)==1 & sum(tcomb,2)==1 );
    linkprobs2{k}(:, l) = results{k}.marlls(:, indIndiv) -  results{k}.marlls(:, ind0);    
  end
  [foo, I8{k}] = sort(linkprobs8{k}(:), 'descend');
  %[foo, I2{k}] = sort(linkprobs2{k}(:), 'descend');
  val8{k} = M(I8{k});
  auc8{1}(k) = drosPlotROC(val8{k}', sum(val8{k}), length(val8{k}), styles{k}, 'LineWidth',LinWidth);
  hold on
  %val2{k} = M(I2{k});
  %auc2{1}(k) = drosPlotROC(val2{k}', sum(val2{k}), length(val2{k}), styles{k+2}, 'LineWidth',LinWidth);  
  set(gca, 'FontSize', FONTSIZE);
end
plot([0,1], [0,1], 'k:', 'LineWidth',LinWidth);
set(gca, 'FontSize', FONTSIZE);
legend(sprintf('Cond. 1\nAUC=%.2f', auc8{1}(1)),...
       sprintf('Cond. 2\nAUC=%.2f', auc8{1}(2)),...
       sprintf('Two conds\nAUC=%.2f', auc8{1}(3)), 'Location', 'SouthEast');
set(gca, 'FontSize', FONTSIZE);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 10 7]);
title('Any');
hold off;
% 1 PLOT ---------------------------------------------------------------
% ROC CURVES FOR GLOBAL PREDICTION OF SINGLE  LINKS
% END    ---------------------------------------------------------------
%


% 2 PLOT ---------------------------------------------------------------
% SEPARATE ROC CURVES FOR PREDICTION OF EACH LINK
% START  ---------------------------------------------------------------
sepAuc = {};
sepVal8 = {};
sepVal2 = {};
for l=1:3,
  hh(l) = figure;  
  for k=1:NUMCONDS, % two conditions
     [foo, II8] = sort(linkprobs8{k}(:,l), 'descend');
     [foo, II2] = sort(linkprobs2{k}(:,l), 'descend');
     sepVal8 = M(II8,l);
     %sepVal2 = M(II2,l);
     sepAuc8{1}(k) = drosPlotROC(sepVal8', sum(sepVal8), length(sepVal8), styles{k}, 'LineWidth',LinWidth);
     hold on;
     %sepAuc2{1}(k) = drosPlotROC(sepVal2', sum(sepVal2), length(sepVal2), styles{k+2}, 'LineWidth',LinWidth);
     set(gca, 'FontSize', FONTSIZE);
  end
  plot([0,1], [0,1], 'k:', 'LineWidth',LinWidth);
  set(gca, 'FontSize', FONTSIZE);
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperSize', [20 20])
  set(gcf, 'PaperPosition', [0 0 10 7])
  legend(sprintf('Cond. 1\nAUC=%.2f', sepAuc8{1}(1)),...
         sprintf('Cond. 2\nAUC=%.2f', sepAuc8{1}(2)),...
         sprintf('Two conds \nAUC=%.2f', sepAuc8{1}(3)),'Location', 'SouthEast');
  %tt = ['TF' num2str(l)]; 
  title(sprintf('%s', TFnames{l}));
  %title(tt);
  hold off;
end
% 2 PLOT ---------------------------------------------------------------
% SEPARATE ROC CURVES FOR PREDICTION OF EACH LINKS
% END    ---------------------------------------------------------------
%


% 2b PLOT ---------------------------------------------------------------
% SEPARATE ROC CURVES FOR PREDICTION OF PAIR LINKS
% START  ---------------------------------------------------------------
sepPairAuc8 = {};
sepPairVal8 = {};
sepPairAucSingle = {};
sepPairValSingle = {};
cnt = 0;
for k=1:numTFs
for g=(k+1):numTFs
    cnt = cnt + 1;
    hPair(cnt) = figure;  
    for cond=1:NUMCONDS, % two conditions 
      % posterior pair-link probabilites using all possible 8 models 
      % find all the indices inside comb where the TFs "k" and "g" are active 
      indPair = find(tcomb(:,k)==1 & tcomb(:,g)==1);
   
      % posterior probability of the TF-pair under the 8 hypotheses case
      pairLinkprobs8{cond}(:,cnt) = logsumexp(results{cond}.marlls(:,  indPair), 2) - ...
			              logsumexp(results{cond}.marlls(:, setdiff(ind8, indPair)), 2);
          
      [foo, I8] = sort(pairLinkprobs8{cond}(:,cnt), 'descend');  
      sepPairVal8 = prod(M(I8,[k,g]),2);
      sepPairAuc8{1}(cond) = drosPlotROC(sepPairVal8', sum(sepPairVal8), length(sepPairVal8), styles{cond}, 'LineWidth',LinWidth);
      hold on;
      set(gca, 'FontSize', FONTSIZE);
      
      if plotAll == 1
      % posterior probability of the TF-pair by combining *Heuristically*
      % the single TFs models 
      indk = find(tcomb(:,k)==1  &  sum(tcomb,2)==1);
      indg = find(tcomb(:,g)==1  &  sum(tcomb,2)==1);
  
      pairLinkprobsFromSingleTF8{cond}(:,cnt) = results{cond}.marlls(:, indk) + results{cond}.marlls(:, indg)  - ...
  			             results{cond}.marlls(:, ind0) - logsumexp(results{cond}.marlls(:, [ind0  indk indg]),2);   
 
      [foo, I8] = sort(pairLinkprobsFromSingleTF8{cond}(:,cnt), 'descend');  
      sepPairValSingle = prod(M(I8,[k,g]),2);
      sepPairAucSingle{1}(cond) = drosPlotROC(sepPairValSingle', sum(sepPairValSingle), length(sepPairValSingle), styles{cond+2}, 'LineWidth',LinWidth);          
      end
      
    end
    if plotAll == 1
    sepPairAucSingle{1}
    end
    plot([0,1], [0,1], 'k:', 'LineWidth',LinWidth);
    set(gca, 'FontSize', FONTSIZE);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [20 20])
    set(gcf, 'PaperPosition', [0 0 10 7])
    
    if plotAll == 1 
      legend(sprintf('One cond\nAUC=%.2f', sepPairAuc8{1}(1)),...
         sprintf('Two conds \nAUC=%.2f', sepPairAuc8{1}(2)),...
         sprintf('One cond,\nsingle-TF models\nAUC=%.2f', sepPairAucSingle{1}(1)),...
         sprintf('Two conds,\nsingle-TF models\nAUC=%.2f', sepPairAucSingle{1}(2)), 'Location', 'SouthEast');
    else
      legend(sprintf('Cond. 1\nAUC=%.2f', sepPairAuc8{1}(1)),...
         sprintf('Cond. 2\nAUC=%.2f', sepPairAuc8{1}(2)),...
         sprintf('Two conds \nAUC=%.2f', sepPairAuc8{1}(3)),'Location', 'SouthEast');
    end
    title(sprintf('%s & %s', TFnames{k}, TFnames{g}));
    hold off;
end
end
% 2b PLOT ---------------------------------------------------------------
% SEPARATE ROC CURVES FOR PREDICTION OF PAIR LINKS
% END    ---------------------------------------------------------------
%


% 2c PLOT ---------------------------------------------------------------
% GLOBAL ROC CURVE FOR PREDICTION OF PAIR LINKS
% START  ---------------------------------------------------------------
sepPairAuc = {};
sepPairVal8 = [];
sepPairAucSingle = {};
sepPairValSingle = [];
cnt = 0;  
hPairGlob = figure;  
hold on;
for cond=1:NUMCONDS, % two conditions   
    [foo, I8] = sort(pairLinkprobs8{cond}(:), 'descend');  
    for j=1:size(I8,1)
        gene = mod(I8(j), 1000);
        gene(gene==0)=100;          
        TFpair = floor(I8(j)/1000) + 1;
        %TFpair(TFpair==4)=3;
        %[gene TFpair I8(j)]
        %pause
        if TFpair == 1 
           sepPairVal8(j) = prod(M(gene, [1, 2]),2); 
        elseif TFpair == 2
           sepPairVal8(j) = prod(M(gene, [1, 3]),2); 
        elseif TFpair == 3
           sepPairVal8(j) = prod(M(gene, [2, 3]),2); 
        end
    end
    sepPairAuc8{1}(cond) = drosPlotROC(sepPairVal8, sum(sepPairVal8), length(sepPairVal8), styles{cond}, 'LineWidth',LinWidth);
    hold on;
    set(gca, 'FontSize', FONTSIZE);                 
    
    if plotAll == 1
    [foo, I8] = sort(pairLinkprobsFromSingleTF8{cond}(:), 'descend');  
    for j=1:size(I8,1)
        gene = mod(I8(j), 1000);
        gene(gene==0)=100;          
        TFpair = floor(I8(j)/1000) + 1;
        %TFpair(TFpair==4)=3;
        %[gene TFpair I8(j)]
        %pause
        if TFpair == 1 
           sepPairValSingle(j) = prod(M(gene, [1, 2]),2); 
        elseif TFpair == 2
           sepPairValSingle(j) = prod(M(gene, [1, 3]),2); 
        elseif TFpair == 3
           sepPairValSingle(j) = prod(M(gene, [2, 3]),2); 
        end
    end
    sepPairAucSingle{1}(cond) = drosPlotROC(sepPairValSingle, sum(sepPairValSingle), length(sepPairValSingle), styles{cond+2}, 'LineWidth',LinWidth);    
    end
end
if plotAll == 1
    sepPairAucSingle{1}
end
plot([0,1], [0,1], 'k:', 'LineWidth',LinWidth);
set(gca, 'FontSize', FONTSIZE);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 10 7]);
if plotAll == 1
legend(sprintf('One cond\nAUC=%.2f', sepPairAuc8{1}(1)),...
       sprintf('Two conds \nAUC=%.2f', sepPairAuc8{1}(2)),...
       sprintf('One cond, \nsingle-TF models\nAUC=%.2f', sepPairAucSingle{1}(1)),...
       sprintf('Two conds, \nsingle-TF models\nAUC=%.2f', sepPairAucSingle{1}(2)), 'Location', 'SouthEast');   
else
legend(sprintf('Cond. 1\nAUC=%.2f', sepPairAuc8{1}(1)),...
       sprintf('Cond. 2\nAUC=%.2f', sepPairAuc8{1}(2)),...
       sprintf('Two conds \nAUC=%.2f', sepPairAuc8{1}(3)),'Location', 'SouthEast');
end
box on;
title('Any pair');
hold off;
% 2c PLOT ---------------------------------------------------------------
% GLOBAL ROC CURVE FOR PREDICTION OF PAIR LINKS
% END    ---------------------------------------------------------------
%




% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING OF TOP K LISTS OF GENES
% START  ---------------------------------------------------------------
modelprobs = {};
modelinds = {};
TG = 50:50:1000;
h2 = figure;
for k=1:NUMCONDS,
  %subplot(1, 2, k);
  [foo, modelinds{k}] = max(results{k}.marlls, [], 2);
  J = sub2ind(size(results{k}.marlls), 1:length(results{k}.genes), modelinds{k}');
  modelprobs{k} = results{k}.marlls(J') - ...
	logsumexp(results{k}.marlls, 2);

  val{k} = all(M == tcomb(modelinds{k}, :), 2);
  [foo, I2{k}] = sort(modelprobs{k}, 'descend');
  mean(val{k}(I2{k}));
  averates = cumsum(val{k}(I2{k})) ./ (1:1000)';
  plot(TG, 100*averates(TG), styles2{k}, 'LineWidth',LinWidth);
  set(gca, 'FontSize', FONTSIZE);
  hold on
end
plot([50,1000], (100/8)*[1,1], 'k:', 'LineWidth',LinWidth)
xlabel('# of most confident genes');
ylabel('Accuracy of MAP predictions (%)');
set(gca, 'FontSize', FONTSIZE);
legend('Cond. 1', ...
       'Cond. 2', ...
       'Two conds')
       %'Random', 'Location', 'EastOutside')
set(gca, 'FontSize', FONTSIZE);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 10 7])
axis([50 1000 0 100]);
hold off;
% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING OF TOP K LISTS OF GENES
% END    ---------------------------------------------------------------
%


% computation of pair-probabilities using set of 8 hypotheses (models)
J_Pair8 = {};
% computation of pair-probabilities using set of 4 hypotheses (models)
J_Pair4 = {};
% computation of pair-probabilities using by *combining heurisrtically* 
% the single TF models
J_PairFromSingleTF = {};
cnt = 0;
for cond=1:NUMCONDS
for k=1:numTFs
  for g=(k+1):numTFs
  %   
  cnt = cnt + 1;
  % find all the indices inside comb where the TFs "k" and "g" are active 
  indPair = find(tcomb(:,k)==1 & tcomb(:,g)==1);
   
  % posterior probability of the TF-pair under the 8 hypotheses case
  [foo, J_Pair8{cnt}] = sort(logsumexp(results{cond}.marlls(:,  indPair), 2) - ...
			              logsumexp(results{cond}.marlls(:, setdiff(ind8, indPair)), 2), 'descend');
    
  % find the *single* index inside comb where 
  % the TFs "k" and "g" are active 
  indPairSingle = find(tcomb(:,k)==1 & tcomb(:,g)==1 & sum(tcomb,2)==2 );
 
  % find the indices inside comb restricted to 4 models 
  %(0 0; k 0; 0 g; k g )
  ind4 = find(  sum(tcomb(:, setdiff(1:numTFs, [k g]) ), 2)==0  );
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  [foo, J_Pair4{cnt}] = sort(results{cond}.marlls(:, indPairSingle) - ...
  			             logsumexp(results{cond}.marlls(:, setdiff(ind4, indPairSingle)), 2), 'descend');
          
  % posterior probability of the TF-pair by combining *Heuristically*
  % the single TFs models 
  indk = find(tcomb(:,k)==1  &  sum(tcomb,2)==1);
  indg = find(tcomb(:,g)==1  &  sum(tcomb,2)==1);
  
  [foo, J_PairFromSingleTF{cnt}] = sort(results{cond}.marlls(:, indk) + results{cond}.marlls(:, indg)  - ...
  			             results{cond}.marlls(:, ind0) - logsumexp(results{cond}.marlls(:, [ind0  indk indg]),2), 'descend');        

  %% baseline maximum likelihood model
  %indPSinglebase = find(baselinecomb(:,k)==1 & baselinecomb(:,g)==1 & sum(baselinecomb,2)==2 ); 
  %[foo, J_Pairbase{cnt}] = sort(baseline_a.marlls(:, indPSinglebase) - baseline_a.marlls(:, 1),'descend');  
  %             
  end
end
end

% computation of single-link probabilities 
J_indiv8 = {};
J_indiv2 = {};
%J_indbase = {};
cnt = 0; 
for cond=1:NUMCONDS
for k=1:numTFs,
%
  cnt = cnt + 1;
  indSingle = find(tcomb(:,k)==1);
  % single TF probability computed from 8 models 
  [foo, J_indiv8{cnt}] = sort(logsumexp(results{cond}.marlls(:,indSingle),2) - ...
			            logsumexp(results{cond}.marlls(:,setdiff(ind8, indSingle)),2),'descend');
         
  % find the *single* index inside comb where only "k" is active 
  indIndiv = find(tcomb(:,k)==1 & sum(tcomb,2)==1 );     
  
  % single TF probability computed from 2 models             
  [foo, J_indiv2{cnt}] = sort(results{cond}.marlls(:, indIndiv) - results{cond}.marlls(:, ind0),'descend');
  
  %% baseline maximum likelihood evaluation         
  %[foo, J_indbase{k}] = sort(baseline_a.marlls(:, k+1) - baseline_a.marlls(:, 1), 'descend');
%  
end
end


% 4 PLOT ---------------------------------------------------------------
% RANKING OF TOP K LISTS OF GENES FOR SINGLE TF LINKS
% START    -------------------------------------------------------------
h3 = figure; 
cnt = 0;
for cond=1:NUMCONDS 
for k=1:numTFs,
  cnt = cnt + 1;  
  figure(h3);
  subplot(3, 4, cnt);
  %drosPlotAccuracyBars({J_indiv8{k},  J_indiv2{k}, J_indbase{k}}, M(:, k), T);
  drosPlotAccuracyBars({J_indiv8{cnt},  J_indiv2{cnt}}, M(:, k), T);
  if cond == 1
  title(sprintf('%s: %d cond', TFnames{k}, cond));
  else
  title(sprintf('%s: %d conds', TFnames{k}, cond));   
  end  
end
end
figure(h3);
subplot(3, 4, cnt+2);
bar(rand(2));
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-8', 'Posterior-2');
%legend('Posterior-32', 'Posterior-16',  'Posterior-2', 'Baseline');
% 4 PLOT ---------------------------------------------------------------
% RANKING OF TOP K LISTS OF GENES FOR SINGLE TF LINKS
% END    ---------------------------------------------------------------


% 5 PLOT ---------------------------------------------------------------
% RANKING OF TOP K LISTS OF GENES FOR  LINK-PAIRS
% START    -------------------------------------------------------------
h4 = figure; 
cnt = 0;
for cond=1:NUMCONDS
for k=1:numTFs 
for g=(k+1):numTFs
  cnt = cnt + 1;    
  figure(h4);
  subplot(3, 4, cnt);
  %drosPlotAccuracyBars({J_joint32{cnt}, J_joint16{cnt}, J_joint4{cnt},  J_joint2{cnt}, J_jointbase{cnt}}, prod(M(:, [k g]), 2), T);
  drosPlotAccuracyBars({J_Pair8{cnt}, J_Pair4{cnt},  J_PairFromSingleTF{cnt}}, prod(M(:, [k g]), 2), T);
  if cond == 1
  title(sprintf('%s & %s: %d cond', TFnames{k}, TFnames{g}, cond));
  else
  title(sprintf('%s & %s: %d conds', TFnames{k}, TFnames{g}, cond));   
  end   
  end
end
end
figure(h4);
subplot(3, 4, cnt+2);
bar(rand(3));
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-8', 'Posterior-4',  'Posterior from single-TF models');
% 5 PLOT ---------------------------------------------------------------
% RANKING OF TOP K LISTS OF GENES FOR  LINK-PAIRS
% END    ---------------------------------------------------------------


if printPlot
   % plot the ground truth TF
   %for j=1:3
   %   figure; 
   %   hold on;  
   %   plot(TimesF, exp(TFs(j,:,1)), 'b'); 
   %   plot(TimesF, exp(TFs(j,:,2)), 'r'); 
   %   plot(TimesF, exp(TFs(4,:,1)), 'k-.'); 
   %end
   print(hPairGlob, '-depsc2', [ddir 'toy_pairGlobal_roc' '.eps']);
   print(h1, '-depsc2', [ddir 'toy_link_roc' '.eps']);
   print(h2, '-depsc2', [ddir 'toy_map_accuracy' '.eps']);
   % plot the seprate roc for each TF 
   for l=1:3
       print(hh(l), '-depsc2', [ddir 'toy_link_roc_TF' num2str(l) '.eps']);
       print(hPair(l) , '-depsc2', [ddir 'toy_pair_roc_TF' num2str(l) '.eps']);
   end
   print(h3, '-depsc2', [ddir 'toy_singleLinks_bars' '.eps']);
   print(h4, '-depsc2', [ddir 'toy_pairsOfTFs_bars' '.eps']);

   % print also plots from the training phase
   %load datasets/toy4TFs28_June_10.mat;
   %load demtoy_dataTwoConds109-Jul-2010.mat 
   %gpmtfPlot(model, samples, 'toyTwoConds', {'TF1', 'TF2', 'TF3'}, 1, 'figures/')
   % 
   %load datasets/toy4TFs28_June_10.mat;
   %load demtoy_dataOneCond108-Jul-2010.mat
   %gpmtfPlot(model, samples, 'toyOneCond', {'TF1', 'TF2', 'TF3'}, 1, 'figures/')  
end