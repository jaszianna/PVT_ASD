clear all; clc;
% addpath('d:\Program files\eeglab14_1_2b\');
% eeglab
% MAKE SURE THAT DATA SEGMENTS ARE ORGANIZED IN ROWS INSTEAD OF COLUMNS (because of hanning windowing)!!!
% maindir = 'e:\work\KOKI_mouse_sleep\';
psd_dir = 'h:\work\KOKI_mouse_sleep\ÁBRA_20221111\Inhibited SwiChR\PRE\'; %[maindir];
fPath_out = 'h:\work\KOKI_mouse_sleep\ÁBRA_20221111\Inhibited SwiChR\PRE\';
fs = 2000; % sampling frequency
% fs = 3000; % sampling frequency
fft_length = 4; % fft data length in seconds
zero_padding = 0; % zero padding length in seconds
overlap = 0.5; % 50% overlap
fqs=[2 48]; % relevan frequencies for plotting and linear fitting
freq_p = [0:0.25:125]; % if you want to use periodogram function!
fPath = strcat(psd_dir,'\');
List = dir(fullfile(fPath,'*.mat*')); % list all mat files in the PSD_DIR

for loop = 1:length(List)
    psd_all = [];
    EEG =load([fPath,List(loop).name]);
    fns = fieldnames(EEG);
    num_segments = numel(fieldnames(EEG));
    c = 1; % loop counting (for the ovetlapping sub-segments)
    
    for segment = 1:num_segments % PSD is calculated for each segment in the EEG structure corresponding to the same recording
        fft_data = EEG.(fns{segment});
        m = 1; % start datapoint
        
        while m  <= (length(fft_data)-(fft_length*fs)) % segment length minus fft window

            data_fft = fft_data(m:m+fft_length*fs-1);
            padding = zeros(1,zero_padding*fs); % initialize zero padding
            hann_win = hann(length(data_fft(1,:)))'; % initialize hanning window
            data = [(data_fft.*hann_win),padding]; 

            N = length(data);
            xdf = fft(data);
            xdft = xdf(1:N/2+1); %deleting the mirrored second half of the FFT
            psdx = (1/(fs*N)) * abs(xdft).^2; % compute the square of the absolute value
            psdx(2:end-1) = 2*psdx(2:end-1); % add the second half
            freq = 0:fs/length(data):fs/2;
        
%             psdx_pg = periodogram(data,rectwin(length(data)),freq_p,fs); % calculate power by periodogram
            
            psd_all(c,:) = psdx;
%             psd_pg(c,:) = psdx_pg; % periodogram
            
            clear psdx;
            c = c+1;
            %             end
            
            m = m+(fft_length*fs*overlap);
        end
        
    end
    
    
a =find(freq==fqs(1)); b = find(freq==fqs(2));
% c =find(freq_p==fqs(1)); d = find(freq_p==fqs(2)); % for periodogram
% function
rf = freq(a:b);
% rf_p = freq_p(c:d); % for periodogram
psd_data = squeeze(mean(psd_all(:,:),1));
% psd_pg_data = squeeze(mean(psd_pg(:,:),1)); % for periodogram
figure(loop)
plot(log(rf),log(psd_data(a:b)),'b');
xticks([log(4) log(6) log(8) log(10) log(12)  log(14) log(16) log(20) log(30) log(40) log(50)]);
xticklabels({'4','6','8','10','12','14','16','20','30','40','50'});
xlabel('Frequency (Hz)'), ylabel('power(log) (\muV^2)');
title('NREM power spectum');
% hold on
% plot(log(rf_p),log(psd_pg_data(c:d)),'r'); % for periodogram
% hold off

% figure(2)
% plot((rf),log(psd_data(a:b)),'b');
saveas(gcf, [fPath_out, 'psd_',List(loop).name(1:end-4),'.tif']);

% save([fPath,'psd.dat'],'psd_data','-ascii');
dlmwrite([fPath_out,cell2mat(fns),'psd_NREM.dat'], psd_data', 'delimiter','\t');
clear psd_all; 
clear psd_data;
close all;
end

clear all;

% 
