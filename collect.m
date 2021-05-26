 

 cnt = 0;
 g=1:12;
 name = 'multitf5a_2010-04-22_m12_r';
 for j=1:size(g,2),
     name1 = [name num2str(g(j)) '.mat'];
     load(['../honkela/scratch/multitf5_results/' name1]);
     for i=1:size(testGene,1)
        if isfield(testGene{i,4},'kinetics')
           cnt = cnt + 1;
           Kin(cnt,:) = median(testGene{i,4}.kinetics,2)';
        end
     end
end
  

