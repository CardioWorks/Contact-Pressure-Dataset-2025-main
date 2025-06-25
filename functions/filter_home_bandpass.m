%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        filter_home_bandpass.m
% Function: some FIR and IIR bandpass filter
% Last update date:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Notes of Parameter Settting
% filter_type------------filter type
% order------------------filter order(level-Wavelet or point--Fir)
% raw_data---------------raw signal
% sample_freq------------sample frequency
% fL---------------------StopFrequency1
% fH---------------------StopFrequency2
% For more information, please search and read the matlab help file.

function [filter_data,filter_SOS] = filter_home_bandpass(raw_data,filter_type,order,sample_freq,fL,fH)
Fn = sample_freq/2;

switch filter_type
    case 1
        [A,B,C,D] = butter(order,[fL fH]/Fn);
        [filter_SOS,g] = ss2sos(A,B,C,D);
        filter_data = filtfilt(filter_SOS,g,raw_data);
    case 2
        [A,B,C,D] = cheby1(order,0.1,[fL fH]/Fn);
        [filter_SOS,g] = ss2sos(A,B,C,D);
        filter_data = filtfilt(filter_SOS,g,raw_data);
    case 3
        [A,B,C,D] = cheby2(order,20,[fL fH]/Fn);
        [filter_SOS,g] = ss2sos(A,B,C,D);
        filter_data = filtfilt(filter_SOS,g,raw_data);
    case 4
        [A,B,C,D] = ellip(order,0.1,30,[fL fH]/Fn);
        [filter_SOS,g] = ss2sos(A,B,C,D);
        filter_data = filtfilt(filter_SOS,g,raw_data);
    case 5
        d = fir1(order,[fL fH]/Fn,'bandpass');
        filter_data = filtfilt(d,1,raw_data);
        filter_SOS = d;
    case 6
        d = designfilt('bandpassfir','FilterOrder',order,'StopbandFrequency1',fL-0.2,'PassbandFrequency1',fL,...
            'PassbandFrequency2',fH,'StopbandFrequency2',fH+2,'DesignMethod','ls','SampleRate',sample_freq);
        filter_data = filtfilt(d,raw_data);
        filter_SOS = d;
    case 7
        filter_data = smooth(raw_data,order);
        filter_SOS = 0;  %No meaning
    case 8
        filter_data = medfilt1(raw_data,order);
        filter_SOS = 0;  %No meaning
    case 9
        filter_data= wden(raw_data,'modwtsqtwolog','s','mln',order,'db2'); %Wavelet level : order
        filter_SOS = 0;  %No meaning
end

end