% Demos used to generate the results in the paper

% Toy data data  
demToyTrain_3TFs30GenesOneCondImData1;    % Training phase:  running with one condition 
demToyTest_3TFsOneCondAllModels1;               % Test phase:  running with one condition 
demToyTrain_3TFs30GenesTwoCondsImData1;  % Training phase:  running with two conditions 
demToyTest_3TFsTwoCondsAllModels1;            % Test phase:  running with two conditions 

% Drosophila data
demDrosophilaTrain_5TFs92GenesImData1;   % runs with 92 genes  
demDrosophilaTrain_5TFs92GenesImData6;   % filters out 25 genes and re-runs using 25 genes 
demDrosophilaTest_5TFsReallyAllModels4;     % Test/sceeening 
% Drosophila data: baseline model  
demDrosophilaTest_5TFsAllModelsBaseLine1;  % Test/screening 


% Drosophila data: Regression with positivity constraints and 
% cross validation accross all 32 models 
demDrosophilaRegression;