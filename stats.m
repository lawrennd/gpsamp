


for n=1:size(testGene,2)
    %
    for j=1:5
        ok = testGene{n}.Weights(j,:);
        ok(ok>=-0.0316 & ok<=0.0316) = 0; 
        ok(ok~=0)=1;
        prob(n,j) = sum(ok)/3000;
    end
    %
end