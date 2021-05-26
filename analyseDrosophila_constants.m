% USER-specified: Sizes of the ranking sets for the first plot
T1 = [200, 400, 800, 1600, 3200, 6003];
% Sizes of the ranking sets for the second plot
T2 = [50, 100, 200, 400, 800, 1600, 3200];

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

% indices of all 32 models 
ind32 = 1:size(combConstr, 1);
% indices that correspond to 16 models (at most 2 TFs) 
ind16 = find(sum(combConstr,2)<=2); 
% index of the zeroth model 
ind0 = find(sum(combConstr,2)==0);
comb16 = combConstr(ind16,:);

numTFs = size(combConstr,2);

pairs = zeros(numTFs*(numTFs-1)/2, 2);
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
    %   
    cnt = cnt + 1;
    pairs(cnt,:) = [k g];
  end
end
