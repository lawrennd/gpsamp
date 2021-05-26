function writePosteriors(genes, post, comb)
%

fid = fopen('posteriors.txt','wt'); 
fprintf(fid, 'geneID\t');
for j=1:size(comb,1)
    for i=1:size(comb,2);
        fprintf(fid,'%d',comb(j,i));
    end
    if j < size(comb,1)
      fprintf(fid,'\t'); 
    else
      fprintf(fid,'\n');
    end
end

for j=1:size(post,1)
    fprintf(fid,'%s\t',genes{j}); 
    for i=1:size(post,2)
        fprintf(fid,'%f\t',post(j,i));
    end
    fprintf(fid,'\n');
end
