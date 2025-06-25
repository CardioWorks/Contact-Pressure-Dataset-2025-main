%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        find_down_zero
% Designed by Yongbo
% Last update date:2017-07-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_point] = find_down_zero(num_point,num_start,num_end,data)
    for j = num_start : num_end
        if data(j) > 0
            break;
        end
    end
%find the down slope zero, first find the positive value as start
    
    for i = j : num_end 
        if data(i)*data(i+1) < 0
           num_point=i;
           break;
        end    
    end
end