%%This code transform the output from the SpykingCircus Clustering Program to Matlab variables

% read "spikes" matlab variable
phyfolder = uigetdir;
folder = uigetdir;
sampRate = 30000;
%amplitudes = double(readNPY([phyfolder filesep 'amplitudes.npy']));
spiketimes = double(readNPY([phyfolder filesep 'spike_times.npy']));
phycluster = double(readNPY([phyfolder filesep 'spike_clusters.npy']));

clusterfile = [phyfolder filesep 'cluster_info.tsv'];
tsvclusters = readtable(clusterfile,'FileType','text');
goodclusters = tsvclusters.cluster_id(cellfun(@(x) contains(x,'good'),tsvclusters.group));
goodch = tsvclusters.ch(cellfun(@(x) contains(x,'good'),tsvclusters.group));
%goodclusters = tsvclusters.cluster_id(ismember(tsvclusters.ch,[28 29 30 31])); % ide azok a csatornak kerulnek, amik a tetrodehoz tartoznak

spikestemp = cell(length(goodclusters),1);

for c = 1:length(goodclusters)
    spikestemp{c,1} = spiketimes(phycluster == goodclusters(c));  
end

%spikes = cell2table(spikestemp,'RowNames',arrayfun(@(x) ['clu_' num2str(x)], goodclusters,'UniformOutput',false),'VariableNames',{'spikes'});

wfhalftime = 0.001; % 1 ms
wfhalftime = round(wfhalftime*sampRate);    
waveform_ALL= NaN(wfhalftime*2+3,4*length(goodclusters));



% channelpool = unique(goodch + 1); % all channels that contain good unit
% for i = 1:length(channelpool)
%     if mod(channelpool(i),4)==0 
%         
%        channelpool = [channelpool;channelpool(i)-1;channelpool(i)-2;channelpool(i)-3];
%        
%     elseif mod(goodchannel,4)==3
%         
%         channelpool = [channelpool;channelpool(i)-1;channelpool(i)-2;channelpool(i)+1];
%         
%     elseif mod(goodchannel,4)==2
%         
%         channelpool = [channelpool;channelpool(i)-1;channelpool(i)+2;channelpool(i)+1];
%         
%     elseif mod(goodchannel,4)==1
%         
%         channelpool = [channelpool;channelpool(i)+1;channelpool(i)+2;channelpool(i)+3];
%         
%     end
% end 
% channelpool  = unique(channelpool);
% 
% eegfilename0 ='110_CH1.continuous';
% init_eeg = load_open_ephys_data([folder filesep eegfilename0]);
% eegmatrix = NaN(length(channelpool),length(init_eeg)/2); % filelength: eeg channel data points
% 
% for k = 1:length(channelpool)
%     eeglong = load_open_ephys_data([folder filesep '110_CH' num2str(channelpool(k))  '.continuous']);
%     eegmatrix(k,:) = eeglong(1:length(eegmatrix));
% end
% 
% eegmatrix = transpose(eegmatrix);
% eegmatrix = [reshape(channelpool,1,[]) ; eegmatrix];

num_of_col_ALL = 1;

for u = 1:size(spikestemp,1)
    waveform = NaN(wfhalftime*2+2,4);
    actual = spikestemp{u,1};
    goodchannel = goodch(u) + 1; % for .continuous files
    
    if mod(goodchannel,4)==0 
        
       startchannel = goodchannel-3;
       
    elseif mod(goodchannel,4)==3
        
        startchannel = goodchannel-2;
        
    elseif mod(goodchannel,4)==2
        
        startchannel = goodchannel-1;
        
    elseif mod(goodchannel,4)==1
        
        startchannel = goodchannel;
        
    end
    
    
    num_of_col = 1;
    for channel = startchannel:1:startchannel+3
%   eeg channel:
    eegfilename = ['110_CH' num2str(channel) '.continuous'];
    eeg = load_open_ephys_data([folder filesep eegfilename]);
%     eeg = eegmatrix(2:end,(eegmatrix(1,:) == channel)); 
    wfchannel = filtfilt(fir1(2048, 250/(sampRate/2), 'high'),1,double(eeg)); % high-pass filter >250Hz
    
    % spike triggers:
    wfnum = round(length(actual)/10); % modify if needed, now: 1/10 of all spike num
    randomindex = randi(length(actual),wfnum,1);
    trigger = actual(randomindex);
    
    % window definition:
    trigger(trigger < (wfhalftime+1)) = [];
    trigger(trigger > length(eeg)-wfhalftime) = [];
    
    wfmatrix = NaN(2*wfhalftime+1,length(trigger));
    for w = 1:length(trigger)
        wfmatrix(:,w) = wfchannel((trigger(w)-wfhalftime):(trigger(w)+wfhalftime));
    end
    waveform(2:end,num_of_col) = mean(wfmatrix,2);
    num_of_col= num_of_col+1;
    end
    waveform = [startchannel:startchannel+3 ; waveform];
    waveform_ALL(:,num_of_col_ALL:num_of_col_ALL+3) = waveform;
    num_of_col_ALL = num_of_col_ALL+4;
end