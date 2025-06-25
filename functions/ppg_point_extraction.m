function [error_code,Loc_PPG,Loc_VPG,Loc_APG] = ppg_point_extraction(error_code,ppg,vpg,apg,third_block)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        ppg_point_extraction.m
% Identify and extract the feature points from PPG, VPG and APG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input:
%error_code----Extracted PPG feature points or not.(0:YES, nonZero:NO)
%sample_time----The sample time of PPG signal
%ppg----PPG waveform of two adjacent heartbeat cycles
%vpg----VPG waveform of two adjacent heartbeat cycles
%apg----APG waveform of two adjacent heartbeat cycles

%% Output:
%error_code----Updated error code for calculation.(0:YES, 1:NO)
%PPG_Loc----The location of PPG feature points
%VPG_Loc----The location of VPG feature points
%APG_Loc----The location of APG feature points

%% Notes: ppg features points location
% They are respectively as below:
%【O】【S】【N】【D】【O_next】in PPG
%【w】【y】【z】【w_next】 in VPG
%【a】【b】【c】【d】【e】【b2】 in APG

%%%%%%
max_num = length(ppg)-1;

%% Error Report Mechanism
% Give each point a initial value(nonZero)
%and the value can be used to check the pass_or_fail status of PPG feature points extraction
num_O = 1;
num_S = 2;
num_N = 3;
num_D = 4;
num_O_next = 5;

num_w = 11;
num_y = 12;
num_z = 13;
num_w_next = 14;

num_a = 21;
num_b = 22;
num_c = 23;
num_d = 24;
num_e = 25;
num_a_next = 26;
num_b2 = 27;

% PPG_Loc = zeros(5);
% VPG_Loc = zeros(4);
% APG_Loc = zeros(6);

%% STEP1----[a: error code = 21]Find 'a' location from APG: num_a
if error_code == 0
    [a_peak,a_peak_loc] = find_peaks(apg, 0.6);       %%%%0.6
    mean_a_peak = mean(a_peak(2:length(a_peak)));
    
    if (a_peak(1)- mean_a_peak)/mean_a_peak > 0.4
        num_a = a_peak_loc(2);
    else
        num_a = find_max(num_a,101,fix(max_num/2),apg);
    end
    
    error_code = ppg_error_report(num_a);
end

% [a_peak,a_peak_loc] = find_peaks(apg,0.5);
% mean_a_peak = mean(a_peak(2:length(a_peak)));
    n = length(a_peak_loc);
%     if l == 7
%         Loc_PPG = struct('Loc_30mmHg',zeros(5,n-3),'Loc_40mmHg',zeros(5,n-3),'Loc_50mmHg',zeros(5,n-3),...
%         'Loc_60mmHg',zeros(5,n-3),'Loc_70mmHg',zeros(5,n-3),'Loc_80mmHg',zeros(5,n-3));
%         Loc_VPG = struct('Loc_30mmHg',zeros(4,n-3),'Loc_40mmHg',zeros(4,n-3),'Loc_50mmHg',zeros(4,n-3),...
%         'Loc_60mmHg',zeros(4,n-3),'Loc_70mmHg',zeros(4,n-3),'Loc_80mmHg',zeros(4,n-3));
%         Loc_APG = struct('Loc_30mmHg',zeros(6,n-3),'Loc_40mmHg',zeros(6,n-3),'Loc_50mmHg',zeros(6,n-3),...
%         'Loc_60mmHg',zeros(6,n-3),'Loc_70mmHg',zeros(6,n-3),'Loc_80mmHg',zeros(6,n-3));
%     end


for s = 2:1:n-2                                                       % a 点迭代次数，25s内循环
    num_a = a_peak_loc(1,s);
    
    
    %% STEP2----[O: error code = 1]Find 'O' location from VPG(zero-crossing point：right to left)：num_O
    if error_code == 0
        for i = num_a : -1 : 2
            if vpg(i)*vpg(i-1) < 0
                num_O = i;
                break;
            end
        end
        
        error_code = ppg_error_report(num_O);
    end
    
    %% STEP3----[w: error code = 11]Find 'w' location from APG(zero-crossing point):num_w
    if error_code == 0
        num_w = find_down_zero(num_w,num_a,max_num,apg);
        error_code = ppg_error_report(num_w);
    end
    
    %% STEP4----[S: error code = 2]Find 'S' location from VPG(zero-crossing point):num_S
    if error_code == 0
        num_S = find_down_zero(num_S,num_w,max_num,vpg);
        error_code = ppg_error_report(num_S);
    end
    
    % %% STEP5----[y: error code = 12]Find 'y' location from VPG(local minimum value):num_y
    % if error_code == 0
    %     num_y = find_min(num_y,num_S,num_S+100,vpg);
    %     error_code = ppg_error_report(num_y);
    % end
    
    %% STEP6----[b: error code = 22]Find 'b' location from third derivate(zero-crossing point):num_b
    if error_code == 0
        num_b = find_up_zero(num_b,num_w,max_num,third_block);
        error_code = ppg_error_report(num_b);
    end
    
    %% STEP7----[c: error code = 23]Find 'c' and 'd' location from third derivate(zero-crossing point):num_c and num_d
    if error_code == 0
        num_c = find_down_zero(num_c,num_b,max_num,third_block);
        num_d = find_up_zero(num_d,num_c,max_num,third_block);
        
        %adjust one or two peak of third_block in b--c interval
        [pks,locs] = findpeaks(third_block(num_b:num_c));
        
        locs = locs + num_b;
        if length(locs) ~= 1
            num_c = find_min(num_c,locs(1),locs(2),third_block);
            num_d = fix((num_c+locs(2))/2);
        end
        
        error_code = ppg_error_report(num_c);
    end
    
    %% STEP5----[y: error code = 12]Find 'y' location from VPG(local minimum value):num_y
    if error_code == 0
        num_y = find_min(num_y,num_b,num_c,vpg);
        error_code = ppg_error_report(num_y);
    end
    
    
    %% STEP8----[e: error code = 25]Find 'e' location from third derivate(zero-crossing point): num_e
    if error_code == 0
        num_e = find_down_zero(num_e,num_d,max_num,third_block);
        error_code = ppg_error_report(num_e);
    end
    
    %% STEP9----[z: error code = 13]Find 'z' location from APG(zero-crossing point):num_z
    % if error_code == 0
    %     num_z = find_down_zero(num_z,num_e,max_num,apg);
    %     error_code = ppg_error_report(num_z);
    % end
    % if error_code == 0
    %     num_z = find_max(num_z,num_e,num_O_next,vpg);
    %     error_code = ppg_error_report(num_z);
    % end
    
    %% STEP10----[N: error code = 3]Find 'N' location from VPG(zero-crossing point):num_N
    % if error_code == 0
    %     if vpg(num_z) > 0    %above 0
    %         num_N = find_up_zero(num_N,num_e,num_z,vpg);
    %     else
    %         num_N =num_e + fix((num_z - num_e)/2);
    %     end
    %
    %     error_code = ppg_error_report(num_N);
    % end
    %
    %     num_D = find_down_zero(num_D,num_b,max_num,third_block);
    %     num_N = find_up_zero(num_N,num_D,max_num,third_block);
    %
    %     %adjust one or two peak of third_block in b--c interval
    %     [pks,locs] = findpeaks(third_block(num_b:num_D));
    %
    %     locs = locs + num_b;
    %     if length(locs) ~= 1
    %         num_D = find_min(num_D,locs(1),locs(2),third_block);
    %         num_N = fix((num_D+locs(2))/2);
    %     end
    
    
    %% STEP11----[D: error code = 4]Find 'D' location from VPG(zero-crossing point):num_D
    % if error_code == 0
    %     if vpg(num_z) > 0    %above 0
    %         num_D = find_down_zero(num_D,num_z,max_num,vpg);
    %     else
    %         num_D =num_z + fix((num_z - num_e)/2);
    %     end
    %
    %     error_code = ppg_error_report(num_D);
    % end
    
    %% STEP12----[a_next and w_next]Find 'a_next' and 'w_next' :num_a_next, num_w_next
    if error_code == 0
        for k =1:1:length(a_peak)
            if a_peak_loc(k) == num_a && length(a_peak) > 1 && k ~= length(a_peak)
                num_a_next = a_peak_loc(k+1);
                break;
            end
        end
        
        error_code = ppg_error_report(num_a_next);
    end
    
    if error_code == 0
        num_w_next = find_down_zero(num_w_next,num_a_next,max_num,apg);
        error_code = ppg_error_report(num_w_next);
    end
    
    %% STEP13----[O_next: error code = 5]Find 'O_next' location from VPG(zero-crossing point):num_O_next
    % %search 'O_next' from right to left
    if error_code == 0
        for i = num_a_next : -1 : 2
            if vpg(i)*vpg(i-1) < 0
                num_O_next = i;
                break;
            end
        end
        
        error_code = ppg_error_report(num_O_next);
    end
    
    if error_code == 0
        num_z = find_max(num_z,num_e+20,num_O_next-30,vpg);
        error_code = ppg_error_report(num_z);
    end
    
    if error_code == 0
        if vpg(num_z) > 0    %above 0
            num_N = find_up_zero(num_N,num_e,num_z,vpg);
        else
            num_N =num_e + fix((num_z - num_e)/2);
        end
        if num_N < fix(num_e+num_z)/2
            num_N = num_e + fix((num_z - num_e)/2);
        end
        error_code = ppg_error_report(num_N);
    end
    
    if error_code == 0
        if vpg(num_z) > 0    %above 0
            num_D = find_down_zero(num_D,num_z,max_num,vpg);
        else
            num_D =num_z + fix((num_z - num_e)/2);
        end
        if num_D > num_O_next
            num_D = num_z + fix((num_O_next - num_z)/2);
        end
        error_code = ppg_error_report(num_D);
    end
    
    
    
    %% STEP14----[b2: error code = 19]Find 'b2' from PPG:num_b2
    % if error_code == 0
    %     num_b2 = num_d;
    %     for j= num_S : num_N
    %         if abs(ppg(j) - ppg(num_b)) < 0.02
    %             num_b2=j;
    %             break;
    %         end
    %         num_b2 = num_d;
    %     end
    %
    %     error_code = ppg_error_report(num_b2);
    % end
    if error_code == 0
        num_b2 = find_min(num_b2,num_O_next,num_O_next+200,apg);
        error_code = ppg_error_report(num_b2);
    end
    
    
    %% PPG Feature Point Location
    if error_code == 0
%         if l == 7
%             Loc_PPG.Loc_30mmHg(1,s-1) = num_O;         Loc_PPG.Loc_30mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_30mmHg(3,s-1) = num_N;         Loc_PPG.Loc_30mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_30mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_30mmHg(1,s-1) = num_w;         Loc_VPG.Loc_30mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_30mmHg(3,s-1) = num_z;         Loc_VPG.Loc_30mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_30mmHg(1,s-1) = num_a;         Loc_APG.Loc_30mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_30mmHg(3,s-1) = num_c;         Loc_APG.Loc_30mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_30mmHg(5,s-1) = num_e;         Loc_APG.Loc_30mmHg(6,s-1) = num_b2;
            
%         elseif l == 9
%             Loc_PPG.Loc_40mmHg(1,s-1) = num_O;         Loc_PPG.Loc_40mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_40mmHg(3,s-1) = num_N;         Loc_PPG.Loc_40mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_40mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_40mmHg(1,s-1) = num_w;         Loc_VPG.Loc_40mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_40mmHg(3,s-1) = num_z;         Loc_VPG.Loc_40mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_40mmHg(1,s-1) = num_a;         Loc_APG.Loc_40mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_40mmHg(3,s-1) = num_c;         Loc_APG.Loc_40mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_40mmHg(5,s-1) = num_e;         Loc_APG.Loc_40mmHg(6,s-1) = num_b2;
            
%         elseif l == 11
%             Loc_PPG.Loc_50mmHg(1,s-1) = num_O;         Loc_PPG.Loc_50mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_50mmHg(3,s-1) = num_N;         Loc_PPG.Loc_50mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_50mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_50mmHg(1,s-1) = num_w;         Loc_VPG.Loc_50mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_50mmHg(3,s-1) = num_z;         Loc_VPG.Loc_50mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_50mmHg(1,s-1) = num_a;         Loc_APG.Loc_50mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_50mmHg(3,s-1) = num_c;         Loc_APG.Loc_50mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_50mmHg(5,s-1) = num_e;         Loc_APG.Loc_50mmHg(6,s-1) = num_b2;
            
%         elseif l == 13
%             Loc_PPG.Loc_60mmHg(1,s-1) = num_O;         Loc_PPG.Loc_60mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_60mmHg(3,s-1) = num_N;         Loc_PPG.Loc_60mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_60mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_60mmHg(1,s-1) = num_w;         Loc_VPG.Loc_60mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_60mmHg(3,s-1) = num_z;         Loc_VPG.Loc_60mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_60mmHg(1,s-1) = num_a;         Loc_APG.Loc_60mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_60mmHg(3,s-1) = num_c;         Loc_APG.Loc_60mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_60mmHg(5,s-1) = num_e;         Loc_APG.Loc_60mmHg(6,s-1) = num_b2;
%             
%         elseif l == 15
%             Loc_PPG.Loc_70mmHg(1,s-1) = num_O;         Loc_PPG.Loc_70mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_70mmHg(3,s-1) = num_N;         Loc_PPG.Loc_70mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_70mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_70mmHg(1,s-1) = num_w;         Loc_VPG.Loc_70mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_70mmHg(3,s-1) = num_z;         Loc_VPG.Loc_70mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_70mmHg(1,s-1) = num_a;         Loc_APG.Loc_70mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_70mmHg(3,s-1) = num_c;         Loc_APG.Loc_70mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_70mmHg(5,s-1) = num_e;         Loc_APG.Loc_70mmHg(6,s-1) = num_b2;
            
%         else
%             Loc_PPG.Loc_80mmHg(1,s-1) = num_O;         Loc_PPG.Loc_80mmHg(2,s-1) = num_S;
%             Loc_PPG.Loc_80mmHg(3,s-1) = num_N;         Loc_PPG.Loc_80mmHg(4,s-1) = num_D;
%             Loc_PPG.Loc_80mmHg(5,s-1) = num_O_next;
%             
%             Loc_VPG.Loc_80mmHg(1,s-1) = num_w;         Loc_VPG.Loc_80mmHg(2,s-1) = num_y;
%             Loc_VPG.Loc_80mmHg(3,s-1) = num_z;         Loc_VPG.Loc_80mmHg(4,s-1) = num_w_next;
%             
%             Loc_APG.Loc_80mmHg(1,s-1) = num_a;         Loc_APG.Loc_80mmHg(2,s-1) = num_b;
%             Loc_APG.Loc_80mmHg(3,s-1) = num_c;         Loc_APG.Loc_80mmHg(4,s-1) = num_d;
%             Loc_APG.Loc_80mmHg(5,s-1) = num_e;         Loc_APG.Loc_80mmHg(6,s-1) = num_b2;
%         end

            Loc_PPG(1,s-1) = num_O;         Loc_PPG(2,s-1) = num_S;
            Loc_PPG(3,s-1) = num_N;         Loc_PPG(4,s-1) = num_D;
            Loc_PPG(5,s-1) = num_O_next;
            
            Loc_VPG(1,s-1) = num_w;         Loc_VPG(2,s-1) = num_y;
            Loc_VPG(3,s-1) = num_z;         Loc_VPG(4,s-1) = num_w_next;
            
            Loc_APG(1,s-1) = num_a;         Loc_APG(2,s-1) = num_b;
            Loc_APG(3,s-1) = num_c;         Loc_APG(4,s-1) = num_d;
            Loc_APG(5,s-1) = num_e;         Loc_APG(6,s-1) = num_b2;

    end
end
    

end

function flag_search = ppg_error_report(initial_value)
if initial_value > 0 && initial_value < 30
    flag_search = initial_value;
else
    flag_search = 0;
end
end

function [num_point] = find_up_zero(num_point,num_start,num_end,data)
for j = num_start : num_end
    if data(j) < 0
        j;
        break;
    end
end
%find the up slope zero, find the negative value as start firstly.

for i = j : num_end
    if data(i)*data(i+1) < 0
        num_point=i;
        break;
    end
end
num_point;
end

function [num_point] = find_down_zero(num_point,num_start,num_end,data)
for j = num_start : num_end
    if data(j) > 0
        break;
    end
end
%find the down slope zero, first find the positive value as start firstly.

for i = j : num_end
    if data(i)*data(i+1) < 0
        num_point=i;
        break;
    end
end
end

function [num_point]=find_max(num_point,num_start,num_end,data)
temp = data(num_start);
for i = num_start:num_end
    if temp <= data(i)
        temp = data(i);
        num_point = i;
    end
end
end

function [num_point] = find_min(num_point,num_start,num_end,data)
temp = data(num_start);
for i = num_start:num_end
    if temp >= data(i)
        temp = data(i);
        num_point = i;
    end
end
end

function [maxv,maxl] = find_peaks(data,ratio)
[num_max]=find_max(0,1,floor(length(data)/2),data);
threshold_height = ratio*data(num_max)*0.8;

[maxv,maxl] = findpeaks(data,'MinPeakHeight',threshold_height,'MinPeakDistance',180);
end