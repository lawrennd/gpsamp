
load dataBarencoOption_0_Genes_50.mat


Genes = zeros(50,7,3);
GeneVars = zeros(50,7,3);
for r=1:3 
for j=1:50
   Genes(j,:,r) = data{r}.Ytrain{j};
   GenesVars(j,:,r) = data{r}.Yvar{j};
end
end
TimesG = data{1}.Xtrain';

save dataBarenco50Genes data Genes GenesVars TimesG; 
