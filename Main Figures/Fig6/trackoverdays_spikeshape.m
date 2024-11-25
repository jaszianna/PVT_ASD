%%This is a modified code of Balázs Hangya, we computed similarity scores with this code

function trackoverdays_spikeshape(cellids)
%TRACKOVERDAYS   Determine if the same neuron was recorded in different sessions.
%   TRACKOVERDAYS(CELLIDS) estimates the probability that the two clusters,
%   recorded on different days but on the same tetrode, with 150 microns 
%   based on track reconstruction, belong to the same neuron. The method is
%   based on Fraser and Schwartz (2012). Three similarity scores are used:
%   maximum waveform crosscorrelation (Fisher-transformed), autocorrelogram
%   correlation (Fisher-transformed) and absolute log difference of firing
%   rates. The scores are tested against a bootstrap distribution of pairs
%   of cells recorded in different animals. p=0.05 is used for the first
%   two and p=0.1 is used for the third scores (pairs have to pass all
%   three to be considered the same neuron). The firing rate criterion is
%   more laxed than the other two, because firing rates can change between
%   and even within sessions. Thelo crosscorrelation criterion was omitted
%   from the Fraser paper, since neighbors can change when the electrodes
%   are moved.
%
%   Reference:
%   Fraser GW, Schwartz AB (2012) Recording from the same neurons 
%   chronically in motor cortex, J Neurophys, 107: 1970 –1978.

%   Balazs Hangya, 4-Feb-2021
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu 

% Initialize bootstrap null distribution
NumCells = numel(cellids);  % number of neurons
prms = nchoosek(1:NumCells,2);  % all possible pairs
rprms = prms(randperm(size(prms,1)),:);  % randomize order

mx = 100;   % maximum of pairs considered, to save CPU time

c1 = cellids(rprms(1:mx,1));  % candidates for the first neurons of pairs
c2 = cellids(rprms(1:mx,2));  % candidates for the second neurons of pairs
r1 = getvalue('ratnum',c1);  % animal ID of the second neurons;runanalysis(@cellid2vals)
r2 = getvalue('ratnum',c2);  % animal ID of the second neurons
diffmouseinx = find(r1-r2~=0);  % when the pair is NOT from the same mouse

bno = 50;  % bootstrap sample size
PairsOfCells = [c1(diffmouseinx(1:bno))'...
    c2(diffmouseinx(1:bno))'];   % take a bootstrap sample; pairs are NOT from the same mouse

% Bootstrap null distribution
[FWR, FAR] = deal(nan(1,bno));
%[FWR, FAR, FFR] = deal(nan(1,bno));
for iP = 1:bno
    disp(iP/bno)
    cellid1 = PairsOfCells{iP,1};  % first cellid in the pair
    cellid2 = PairsOfCells{iP,2};  % second cellid in the pair
    [FWR(iP)] = simscores(cellid1,cellid2);   % Fisher-transformed waveform and ACG correlation
    %[FWR(iP), FAR(iP), FFR(iP)] = simscores(cellid1,cellid2);   % Fisher-transformed waveform and ACG correlation
    close all
end
% keyboard
% load('simscore_distributions3.mat')

% Critical values
%cFWR = 0.20;
cFWR = prctile(FWR,5);   % critical value corresponding to the upper 0.05
%cFAR = prctile(FAR,90);
%cFFR = prctile(FFR,10);

% Potential duplicates
[Duplicates, Tested_Pairs] = deal({});
Tested_Pairs_Scores = [];
for iC = 1:NumCells
    disp(iC)
    cellid = cellids{iC};
    % [animalID, sessionID, Tetrode, Unit] = cellid2tags(cellid);
    RatID = getvalue('ratnum',cellid);  % animal ID
    datenum = getvalue('ses',cellid);  % recording date
    tetrode = getvalue('tetrode',cellid);  % recording tetrode
    %DV = getvalue('DVpos',cellid);  % recording depth
    alldatenum = getvalue('ses',cellids);  % recording date for all cells
    allRatIDs = getvalue('ratnum',cellids);  % animal ID for all cells
    alltetrode = getvalue('tetrode',cellids);  % recording tetrode for all cells
    %allDV = getvalue('DVpos',cellids);
    samemouseinx = allRatIDs == RatID;  % indices of cells from the same mouse
    diffsessioninx = alldatenum > datenum;  % indices of cells from later sessions
    sametetrodeinx = tetrode == alltetrode;  % indices of cells recorded on the same tetrode
    potential_duplicates = cellids(samemouseinx&sametetrodeinx&diffsessioninx);
    %potential_duplicates = cellids(samemouseinx&sametetrodeinx&...
    %    diffsessioninx&allDV-DV>=0&allDV-DV<=150);  % cell in the same mouse, on the same tetrode but different session,
    % within 150 um but not above the neuron
    
    % Test the similarity scores of potential duplicates
    NumDup = length(potential_duplicates);  % number of potential duplicates
    for iD = 1:NumDup   % loop through potential duplicates
        cellidd = potential_duplicates{iD};
        [dFWR] = simscores(cellid,cellidd);   % Fisher-transformed waveform and ACG correlation
        %[dFWR, dFAR, dFFR] = simscores(cellid,cellidd);   % Fisher-transformed waveform and ACG correlation
        close all
        Tested_Pairs = [Tested_Pairs; {cellid cellidd}]; %#ok<AGROW>  % store similarity score values
        Tested_Pairs_Scores = [Tested_Pairs_Scores; dFWR]; %#ok<AGROW>  % store similarity score values
        %Tested_Pairs_Scores = [Tested_Pairs_Scores; dFWR, dFAR, dFFR]; %#ok<AGROW>  % store similarity score values
        %if dFWR > cFWR && dFAR > cFAR && dFFR < cFFR  % if all similarity scores are above critical value (below for FR similarity)
        if dFWR < cFWR % if all similarity scores are above critical value (below for FR similarity)
        Duplicates = [Duplicates; {cellid cellidd}];   %#ok<AGROW> % add to list of 'duplicates'
        
        end
        save('tested_pairs','Tested_Pairs','Tested_Pairs_Scores','Duplicates')
        
    end
end

keyboard

% delinx = zeros(1,size(Duplicates,1));
% for k = 1:size(Duplicates,1)
%     datenum1 = getvalue('DateNum',Duplicates(k,1));  % recording date for the first cell
%     datenum2 = getvalue('DateNum',Duplicates(k,2));  % recording date for the second cell
%     if datenum2 < datenum1  % pairs may be counted twice
%         delinx(k) = 1;
%     end
% end
% Duplicates(logical(delinx),:) = [];   % remove duplicated duplicates

% delinx = zeros(1,size(Duplicates,1));
% for k = 1:size(Duplicates,1)
%     fr1 = getvalue('baseline_FR',Duplicates(k,1));  % FR for the first cell
%     fr2 = getvalue('baseline_FR',Duplicates(k,2));  % FR for the second cell
%     dFFR = abs(log(fr1)-log(fr2));  % log difference (assumed log-normal FR distribution)
%     if dFFR > cFFR  % drop pairs with FR difference above 10% bootstrap critical value
%         delinx(k) = 1;
%     end
% end
% Duplicates(logical(delinx),:) = [];   % remove cellx with large FR difference

%Cells2Exclude = unique(Duplicates(:,2));

% -------------------------------------------------------------------------
function [FWR] = simscores(cellid1,cellid2)
% [FWR, FAR, FFR] = simscores(cellid1,cellid2)
% Waveform correlation


path1 = cellid2fnames(cellid1,'tfile') ;
path2 = cellid2fnames(cellid2,'tfile') ;


tetrode1 = getvalue('tetrode',cellid1);  % recording tetrode
tetrode2 = getvalue('tetrode',cellid2);  % recording tetrode
unit1 = getvalue('unit',cellid1);  % recording tetrode
unit2 = getvalue('unit',cellid2);  % recording tetrode
Waveform1 = load([ path1(1:end-2),'_waveforms']);
Waveform2 = load([ path2(1:end-2),'_waveforms']);

wave1 = Waveform1.(['TT',int2str(tetrode1),'_',int2str(unit1)]);
wave2 = Waveform2.(['TT',int2str(tetrode2),'_',int2str(unit2)]);
twave1 = wave1;  % transpose
twave2 = wave2;
catwave1 = twave1(:)';  % concatenate channels
catwave2 = twave2(:)';
ncatwave1 = catwave1 ./ max(abs(catwave1));   % normalize to a maximum of 1 (remove the effect of proportional amplitude changes)
ncatwave2 = catwave2 ./ max(abs(catwave2));
FWR = sum((ncatwave1-ncatwave2).^2);
%[xr, lags] = xcorr(ncatwave1,ncatwave2);  % cross-correlogram of waveforms
%hv = floor(size(wave1,1)/2);  % half of the waveform 'length' (about half ms)
%WR = max(xr(lags>=-hv&lags<=hv));  % waveform correlation: maximal cross-correlation within half-wave times from 0 lag
%WR = WR / length(xr);  % normalize between -1 and 1
%FWR = atanh(WR);  % Fisher's Z-transform (approximately normal distribution)

% Autocorrelogram
% set(0,'DefaultFigureVisible','off');
% 
%  ac1 = acg(cellid1,0.1,'dt',0.005,'minspikeno',0,'maxspikeno',100000);  % autocorrelogram: 100 ms window, 5 ms resolution
%  [ac2, aclags] = acg(cellid2,0.1,'dt',0.005,'minspikeno',0,'maxspikeno',100000);
%  acc = corrcoef(ac1(aclags>0),ac2(aclags>0));  % Pearson's correlation
%  AR = acc(1,2);  % ACG correlation
%  FAR = atanh(AR);  % Fisher's Z-transform (approximately normal distribution)
%   set(0,'DefaultFigureVisible','on');
% % 
% % Baseline firing rate
% fr1 = getvalue('baseline_FR',cellid1);  % baseline firing rate
% fr2 = getvalue('baseline_FR',cellid2);
% FFR = abs(log(fr1)-log(fr2));  % log difference (assumed log-normal FR distribution)

% function waveform = extractWaveforms(cellid)
% 
% folder = cellid2fnames(cellid,'Session'); % modify
% sampRate = 30000;
% 
% % read "spikes" matlab variable
% phyfolder = [folder, filesep, 'Continuous_Datasc.GUI'];
% %amplitudes = double(readNPY([phyfolder filesep 'amplitudes.npy']));
% spiketimes = double(readNPY([phyfolder filesep 'spike_times.npy']));
% phycluster = double(readNPY([phyfolder filesep 'spike_clusters.npy']));
% 
% clusterfile = [phyfolder filesep 'cluster_info.tsv'];
% tsvclusters = readtable(clusterfile,'FileType','text');
% goodclusters = tsvclusters.cluster_id(cellfun(@(x) contains(x,'good'),tsvclusters.group));
% %goodclusters = tsvclusters.cluster_id(ismember(tsvclusters.ch,[28 29 30 31])); % ide azok a csatornak kerulnek, amik a tetrodehoz tartoznak
% 
% spikestemp = cell(length(goodclusters),1);
% 
% for c = 1:length(goodclusters)
%     spikestemp{c,1} = spiketimes(phycluster == goodclusters(c));  
% end
% 
% spikes = cell2table(spikestemp,'RowNames',arrayfun(@(x) ['clu_' num2str(x)], goodclusters,'UniformOutput',false),'VariableNames',{'spikes'});
% 
% wfhalftime = 0.001; % 1 ms
% wfhalftime = round(wfhalftime*sampRate);    
% waveform = NaN(wfhalftime*2+1,4);
% for u = 1:size(spikes,1)
%     actual = spikestemp{u,1};
%     goodchannel = tsvclusters.ch(goodclusters(u),1) + 1; % for .continuous files
%     startchannel = goodchannel - mod(goodchannel,4) + 1;
%     
%     for channel = startchannel:1:startchannel+3
% %   eeg channel:
%     eegfilename = ['110_CH' num2str(channel) '.continuous'];
%     eeg = load_open_ephys_data([folder filesep eegfilename]);
%     wfchannel = filtfilt(fir1(2048, 250/(sampRate/2), 'high'),1,double(eeg)); % high-pass filter >250Hz
%     
%     % spike triggers:
%     wfnum = round(length(actual)/10); % modify if needed, now: 1/10 of all spike num
%     randomindex = randi(length(actual),wfnum,1);
%     trigger = actual(randomindex);
%     
%     % window definition:
%     trigger(trigger < (wfhalftime+1)) = [];
%     trigger(trigger > length(eeg)-wfhalftime) = [];
%     
%     wfmatrix = NaN(2*wfhalftime+1,length(trigger));
%     for w = 1:length(trigger)
%         wfmatrix(:,w) = wfchannel((trigger(w)-wfhalftime):(trigger(w)+wfhalftime));
%     end
%     waveform(:,channel) = mean(wfmatrix,2);
%     end
% end

