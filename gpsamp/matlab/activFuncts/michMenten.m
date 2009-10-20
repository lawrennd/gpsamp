function fx = michMenten(ff, W, NetX)
% function fx = michMenten(x, W, netX)
% 
% Description: Impements the Michaelis-Menten multiple TF activation 
%              The single TF Michaelis-Menten equations are obtained as special cases
%
% Inputs: 
%    -- ff:This is NumoFTs x NumOfTimes size matrix that contains  
%           that stores TFs values.
%    -- W: This a NumOfGenes x NumOfTFs matrix that contains 
%           connectivity weights that take NON-NEGATIVE values. 
%           Each element W(i,j) stores the interaction weight 
%           between the i TF and the j gene. 
%    -- netX: This a NumOfGenes x NumOfTFs matrix that describes the 
%           network connectivity. Each element netX(i,j) can take three
%           values, (-1, 0, 1). When netX(j,i) = -1, then the i TF acts 
%           as repressor for the j gene. When netX(j,i) = 0, the i TF does not 
%           interact with the j gene, and when netX(j,i) = 1, the i TF 
%           activates the j gene. 
%
% Output: fx: This a  NumOfGenes x NumOfTimes matrix that computes the 
%          mulitple Tf Michaelis-Menten function:
%            
%             fx(j,t) =  (\sum_{i in R} w_ji+\sum_{i in A} w_ji f_i(t))/...
%                         (1 + \sum_{i in R} w_ji f(t)  + \sum_{i in A} w_ji f_i(t))
%
%          where 'A' is the set of TF activator for the j gene and
%          'R' is the set of repressors for the j gene. 
%                        
% Notes: When we use only one TF, then we obtain the M-M single TF equation. 
%        In that case the netX is NumOfGenes x 1 vector with -1s and 1s.
%        E.g. for one gene j netX(j) = -1, the equation becomes the repressor
%        M-M given by
%             fx(j,t) =  w_j /(1 + w_j f(t))
%        where you can extract the M-M constant form gamma_j =1/w_j. 


SizF = size(ff,2);

% find all repressors 
R = (netX==-1);
% find all activators 
A = (netX==1);
xp = (W.*A)*ff;
xpR = W.*R;

% combine the TFs
fx = (xp + repmat(sum(xpR,2),[1 SizF]))./(1 + xp + xpR*ff);


