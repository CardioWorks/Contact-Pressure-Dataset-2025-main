%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        find_max
% Designed by Yongbo
% Last update date:2017-07-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_point]=find_max(num_point,num_start,num_end,data)
    temp = data(num_start);   
    for i = num_start:num_end
       if temp <= data(i)
          temp = data(i);
          num_point = i;
       end
    end
end