%% %%%%%%%%%%%%%%------Perfusion index calculation for PPG data----%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%-----Authors Ziyi Wang and Mohamed Elgendi----%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%-------Last updated: 24 June, 2025------%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%------STEP ONE---Calculate Amplitude & Average----%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd));
cd '\Contact Pressure-PPG (CP-PPG) Dataset'
cp_vars = {'cp_30mmHg', 'cp_40mmHg', 'cp_50mmHg', 'cp_60mmHg', 'cp_70mmHg', 'cp_80mmHg'};
for i = 1:1:142
    %%%%%%%%%%%-------Read data-------%%%%%%%%%%%%%
    disp(i)
    filename = [num2str(i), '.mat'];
    PPG_data = load(filename);
    for j = 1:1:6
        varname = cp_vars{j};
        pressure = varname(4:end);
        
        if  isfield(PPG_data, varname) && ~isempty(PPG_data.(varname))
            %%%%%%%%%%%-------Preprocessing-------%%%%%%%%%%%%%
            PPG_raw = PPG_data.(varname)(:, 2);   %% Inrared light wavelength
            PPG_fill = filloutliers(PPG_raw, "previous", "grubbs");
            PPG_fill = PPG_fill(~isnan(PPG_fill));
            filter_type = 1;     filter_order = 6;     f_sample = 400;     fL = 0.1;     fH = 10;
            [PPG_filt,~] = filter_home_bandpass(PPG_fill, filter_type, filter_order, f_sample, fL, fH);
            PPG_filt = flipud(PPG_filt);
            PPG_nor = mapminmax( (PPG_filt)', 0, 1);
            
            VPG = diff( PPG_nor );   APG = diff(VPG);   PPG_third = diff(APG);
            loc_varname = ['Loc_PPG_', pressure];
            
            try
                [~, loc_result, ~, ~] = ppg_point_extraction(0, PPG_nor, VPG, APG, PPG_third);
            catch ME
                
                fprintf('Error processing ID%d, %s: %s\n', i, varname, ME.message);
                continue; 
            end
            
            assignin('base', loc_varname, loc_result);
            amplitude_var = ['Amplitude_', pressure];
            average_var = ['Average_', pressure];
            
            current_loc = eval(loc_varname);
            beats = size(current_loc, 2);
            for beat = 1:1:beats
                %%%%%%%%%%%-------Calculate amplitude-------%%%%%%%%%%%%%
                eval([amplitude_var, '(i, beat) = PPG_nor(1, current_loc(2,beat)) - PPG_nor(1, current_loc(1,beat));']);
                
                %%%%%%%%%%%%-------Calculate average-------%%%%%%%%%%%%%%
                PPG_flip = flipud(PPG_fill);
                PPG_nor_b = mapminmax(PPG_flip', 0, 1);
                eval([average_var, '(i, beat) = mean(PPG_nor_b(1, current_loc(1,beat):current_loc(5,beat)));']);
            end
        end
    end
    
end

%% %%%%%%%%%%%%----------STEP TWO---Calculate PI----------%%%%%%%%%%%%%%%%
cp_vars = {'30mmHg', '40mmHg', '50mmHg', '60mmHg', '70mmHg', '80mmHg'};
for p = 1:1:6
    pressure = cp_vars{p};
    
    amplitude_var = ['Amplitude_', num2str(pressure)];
    average_var = ['Average_', num2str(pressure)];
    perfusion_var = ['Perfusion_', num2str(pressure)];
    trimmed_var = ['Perfusion_', num2str(pressure), '_trimmed'];
    
    for i = 1:142 
        row_value = eval([amplitude_var, '(i, :)']);
        non_zero_idx = find(row_value ~= 0);
        valid_amps = row_value(non_zero_idx);
        valid_avgs = eval([average_var, '(i, non_zero_idx)']);
        %%%%%%%%%%%%---------Calculate PI--------%%%%%%%%%%%%%%
        valid_perf = valid_amps ./ valid_avgs;
        eval([perfusion_var, '(i, 1:length(valid_perf)) = valid_perf;']);
        %%%%%%%%%%%%-------Trimmed average-------%%%%%%%%%%%%%%
        current_data = eval([perfusion_var, '(i, :)']);
        current_data = current_data(find(current_data ~= 0));
        data_ascend = sort(current_data);
        b = length(current_data);
        a = floor(0.15*b);
        data_ascend(1:a) = 0;
        data_ascend(b-a+1:end) = 0;
        c = sum(data_ascend)/(b-2*a);
        eval([trimmed_var, '(i) = c;']);
    end
end

%% %%%%%%-----STEP THREE---Calculate PI after trimmed average-----%%%%%%%%%
Perfusion_30mmHg_trimmed(2, :) = Perfusion_40mmHg_trimmed;
Perfusion_30mmHg_trimmed(3, :) = Perfusion_50mmHg_trimmed;
Perfusion_30mmHg_trimmed(4, :) = Perfusion_60mmHg_trimmed;
Perfusion_30mmHg_trimmed(5, :) = Perfusion_70mmHg_trimmed;
Perfusion_30mmHg_trimmed(6, :) = Perfusion_80mmHg_trimmed;

for i = 1:1:6
        DATA = Perfusion_30mmHg_trimmed(i,:);
        finiteIndices = isfinite(DATA);
        DATA = DATA(finiteIndices);
        DATA_ascend = sort(DATA);
        b = length(DATA);
        a = floor(0.15*b);
        DATA_ascend(1:a) = 0;
        DATA_ascend(b-a+1:end) = 0;
        c = sum(DATA_ascend)/(b-2*a);
        Perfusion_trimmed(i,1) = c;
end

%% %%%%%%--------STEP FOUR---Calculate b/a feature--------%%%%%%%%%
%%%%%%%%%%%-------Amplitude of num_b divide num_a-------%%%%%%%%%%%%%
b_a_feature = APG( 1, Loc_APG(2,beat) ) / APG( 1, Loc_APG(1,beat) );     