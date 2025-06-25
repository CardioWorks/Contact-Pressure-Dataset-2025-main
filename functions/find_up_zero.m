%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        find_up_zero
% Designed by Yongbo
% Last update date:2017-07-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_point] = find_up_zero(num_point,num_start,num_end,data)
 for j = num_start : num_end
        if data(j) < 0
            j;
            break;
        end
    end
%find the up slope zero, first find the negative value as start
    
    for i = j : num_end 
        if data(i)*data(i+1) < 0
           num_point=i;
           break;
        end    
    end
    num_point;
end