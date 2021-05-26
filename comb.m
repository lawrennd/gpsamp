%if size(comb,1) == 3
  %
  % if comb(c,1) == 1
  %   %
  %   TFset = 3; 
  %   if comb(c,2) == 1
  %      TFset = [3 5];
  %   end
  %   %
  % else
  %   TFset = 5; 
  % end
  %
  %elseif size(comb,1) == 7
  % %
  % if comb(c,1) == 1
  %    TFset = 1; 
  %   if comb(c,2) == 1
  %      TFset = [1 3];
  %      if comb(c,3) == 1
  %           TFset = [1 3 5];
  %      end
  %   else
  %      if comb(c,3) == 1
  %           TFset = [1 5];
  %      end
  %   end
  %   elseif comb(c,2) == 1
  %        TFset = 3;
  %        if comb(c,3) == 1
  %            TFset = [3 5];
  %        end
  %   else   
  %      TFset = 5; 
  % end
  % %
  %end