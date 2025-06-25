%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        find_peaks
% Designed by Yongbo
% Last update date:2017-07-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxv,maxl] = find_peaks(data,ratio,distance)    

[num_max]=find_max(0,1,floor(length(data)/2),data);
threshold_height = ratio*data(num_max);

[maxv,maxl] = findpeaks(data,'MinPeakHeight',threshold_height,'MinPeakDistance',distance);
end