function [options, TimesF] = gpmtfDiscretize(TimesG, options)
% DISCRETIZE the GP function
% (around discr times more the number the discrete time points
% we have gene expressions)


% TimesF discretizes the [-tau_max, T] where T = max(TimesG)
% - First discretize in [0,T] (where we have observed genes) 

% the GPs functions is going to be descretize in  discr*(size(TimesG,2)-1)+1
% points. This number of points must be a odd number 
discr=10;

if (discr*(size(TimesG,2)-1))+1 > 200 
   discr = floor(200/(size(TimesG,2)-1)) + 1; 
end 
if mod(discr*(size(TimesG,2)-1)+1,2) == 0
   discr = discr-1;
end 
step = ((max(TimesG) - min(TimesG))/(size(TimesG(:),1)-1))/discr;
TimesF =[]; TimesF(1) = TimesG(1);
for j=1:size(TimesG(:),1)-1
   TimesF = [TimesF, ((TimesG(j)+step):step:TimesG(j+1))];
   if TimesF(end) ~= TimesG(j+1)
      TimesF = [TimesF, TimesG(j+1)];
   end
end

% - Now discretize in [-tau_max,0) (the "delay" part of the TF) 
options.sizTime = size(TimesF,2);
if options.tauMax > 0
DelayP = -step:-step:-options.tauMax; 
options.tauMax = DelayP(end); 
TimesF = [DelayP(end:-1:1) TimesF];    
end
