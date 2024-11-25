%%This code was created to compute firing rates, burst/high frequency activity, burst (high frequency activity) ratio and autocorrellograms in every states (wake, nest, sleep)

function FiringRate_HighFrequencyActivity()
%inputs:
path='X:\_Projects_2\Midline_anxiety\DATA\_Tetrode\ANALÍZIS\Cellbase_data\A113\220325w\';
resultspath='X:\_Projects_2\Midline_anxiety\DATA\_Tetrode\ANALÍZIS\Cellbase_data';
session='SD_A113_220325';
 
awake_off1=3600; %in seconds
nest_on1=3600; %in seconds
nest_off1=3800; %in seconds
sleep_on1=3800; %in seconds
sleep_off1=4100; % in seconds
awake_off2=1125; %in seconds
nest_on2=1125; %in seconds
nest_off2=1779;%in seconds
sleep_on2=1779; %in seconds
sleep_off2=2369; % in seconds


%------------------------------------------
acg_wind=400; %in ms
burst_maxlag=10; % in ms
baseline_window=[-300,-200]; %in ms
acg_res=5; % in ms

fileList = dir([path,session,'\','TT*.mat']);
[~, timestamps, ~] = load_open_ephys_data([path,session,'\','110_CH12.continuous']);
starttime=timestamps(1);
[~, stims, ~] = load_open_ephys_data([path,session,'\','all_channels.events']);
stims=stims(1:2:end)-starttime;
borders=find(diff(stims)>1);
borders=[0;borders;length(stims)];


for i=1:length(fileList)
    %for i=1:1
    
    TS=load([path,session,'\',fileList(i).name]);
    TS=TS.TS;
    TS=(TS-TS(1))/10000;
%     for q=1:length(borders)-1
%         TS(TS>(stims(borders(q)+1)) & TS<(stims(borders(q+1))))=[];
%     end
    endtext=strfind(fileList(i).name,'_');
    Tetrode=str2double(fileList(i).name(3:endtext-1));
    if Tetrode<9
        awake_off=awake_off1;
        nest_on=nest_on1;
        nest_off=nest_off1;
        sleep_on=sleep_on1;
        sleep_off=sleep_off1;
    elseif Tetrode>8
        awake_off=awake_off2;
        nest_on=nest_on2;
        nest_off=nest_off2;
        sleep_on=sleep_on2;
        sleep_off=sleep_off2;
    end
    awake_filt=TS<nest_on;
    if ~isnan(nest_on)
        nesting_filt=TS>nest_on & TS<nest_off;
    end
    sleep_filt=TS>sleep_on & TS<sleep_off;
    firing_rate_awake=sum(awake_filt)/awake_off; %Hz
    if ~isnan(nest_on)
        firing_rate_nesting=sum(nesting_filt)/(nest_off-nest_on); %Hz
    else
        firing_rate_nesting=0;
    end
    firing_rate_sleep=sum(sleep_filt)/(sleep_off-sleep_on); %Hz
    figure
    bar([firing_rate_awake,firing_rate_nesting,firing_rate_sleep])
    title([session,' ',fileList(i).name(1:end-4),' firing rate changes'])
    ylabel('Firing rate (Hz)')
    set(gca,'XTickLabel',str2mat('Awake','Nesting','Sleep'))
    savefig([path,session,'\',session,'_',fileList(i).name(1:end-4),'_FR.fig'])
    saveas(gcf,[path,session,'\',session,'_',fileList(i).name(1:end-4),'_FR.png'])
    
    
    G=figure
    ax1=subplot(1,3,1)
    x=diff(TS(awake_filt))*1000;
    [~,edges] = histcounts(log10(x));
    histogram(x,10.^edges)
    set(gca, 'xscale','log')
    xlim([10^-1,10^4])
    ylabel('Number of spikes')
    xlabel('ISI (ms)')
    title('awake')
    if ~isnan(nest_on)
        ax2=subplot(1,3,2)
        x=diff(TS(nesting_filt))*1000;
        [~,edges] = histcounts(log10(x));
        histogram(x,10.^edges)
        set(gca,'xscale','log')
        ylabel('Number of spikes')
        xlabel('ISI (ms)')
        title('nesting')
    end
    ax3=subplot(1,3,3)
    x=diff(TS(sleep_filt))*1000;
    [~,edges] = histcounts(log10(x));
    histogram(x,10.^edges)
    set(gca,'xscale','log')
    ylabel('Number of spikes')
    xlabel('ISI (ms)')
    title('sleep')
    linkaxes([ax1 ax2 ax3],'x')
    %suptitle([session,' ',fileList(i).name(1:end-4),' ISIs'])
    maximize_figure(G)
    savefig([path,session,'\',session,'_',fileList(i).name(1:end-4),'_ISI.fig'])
    saveas(gcf,[path,session,'\',session,'_',fileList(i).name(1:end-4),'_ISI.png'])
    
    G=figure
    subplot(1,3,1)
    [~,~,Burst_index_A]=acorr(TS(awake_filt),TS(awake_filt),acg_wind,acg_res,burst_maxlag,baseline_window);
    str=sprintf('Awake, Burst-index=%0.3f',Burst_index_A);
    title(str);
    xlabel('lag (ms)')
    subplot(1,3,2)
    [~,~,Burst_index_N]=acorr(TS(nesting_filt),TS(nesting_filt),acg_wind,acg_res,burst_maxlag,baseline_window);
    str=sprintf('Nesting, Burst-index=%0.3f',Burst_index_N);
    title(str);
    xlabel('lag (ms)')
    subplot(1,3,3)
    [~,~,Burst_index_S]=acorr(TS(sleep_filt),TS(sleep_filt),acg_wind,acg_res,burst_maxlag,baseline_window);
    str=sprintf('Sleep, Burst-index=%0.3f',Burst_index_S);
    title(str);
    xlabel('lag (ms)')
    %suptitle([session,' ',fileList(i).name(1:end-4),' ACGs'])
    maximize_figure(G);
    savefig([path,session,'\',session,'_',fileList(i).name(1:end-4),'_ACG.fig'])
    saveas(gcf,[path,session,'\',session,'_',fileList(i).name(1:end-4),'_ACG.png'])
    
   
    [awake_burst_num,awake_burst_ratio] = burst_detection(TS(awake_filt),burst_maxlag,awake_off);
    [nesting_burst_num,nesting_burst_ratio] = burst_detection(TS(nesting_filt),burst_maxlag,nest_off-nest_on);
    [sleep_burst_num,sleep_burst_ratio] = burst_detection(TS(sleep_filt),burst_maxlag,sleep_off-sleep_on);
   
    
    R_isexist = dir([resultspath,'\','Results_A113_2.mat']);
    
    if isempty(R_isexist)
        Results(1).ID=[session,'_',fileList(i).name];
    else
        Results=load([resultspath,'\','Results_A113_2.mat']);
        Results=Results.Results;
    end
    R_index=find(strcmp({Results.ID}, [session,'_',fileList(i).name])==1);
    if isempty(R_index)
        R_index=length(Results)+1;
        Results(R_index).ID=[session,'_',fileList(i).name];
    end
    Results(R_index).firing_rate_awake=firing_rate_awake;
    Results(R_index).firing_rate_nesting=firing_rate_nesting;
    Results(R_index).firing_rate_sleep=firing_rate_sleep;
    Results(R_index).burst_index_awake=Burst_index_A;
    Results(R_index).burst_index_nesting=Burst_index_N;
    Results(R_index).burst_index_sleep=Burst_index_S;
    Results(R_index).awake_burst_num=awake_burst_num;
    Results(R_index).nesting_burst_num=nesting_burst_num;
    Results(R_index).sleep_burst_num=sleep_burst_num;
    Results(R_index).awake_burst_ratio=awake_burst_ratio;
    Results(R_index).nesting_burst_ratio=nesting_burst_ratio;
    Results(R_index).sleep_burst_ratio=sleep_burst_ratio;
    
    
    save([resultspath,'\','Results_A113_2.mat'],'Results');
    
end








function [ccr, lags,Burst_index] = acorr(ncc1,ncc2,wn,res,burst_maxlag,baseline_window)

% Calculate spike times in milliseconds
sr = 1000;
nc1 = ncc1 * sr;
mn1 = nc1(1);  % only relative spike times count; avoid out of memory
nc1 = nc1 - mn1;
nc1(nc1<0.5*res) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
wn2 = wn / 1000;    % window size in seconds

% Calculate spike times in milliseconds
nc2 = ncc2 * sr;
mn2 = nc2(1);  % only relative spike times count; avoid out of memory
nc2 = nc2 - mn2;
nc2(nc2<0.5*res) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39

% Auto-correlogram
zunit1 = zeros(1,round(nc1(end)/res)+5);
zunit1(round(nc1/res)) = 1;
zunit2 = zeros(1,round(nc2(end)/res)+5);
zunit2(round(nc2/res)) = 1;
[ccr, lags] = xcorr(zunit1,zunit2,wn2*sr/res);     % 1->2; window: -wn ms - wn ms
ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
lags(length(lags)/2+0.5) = [];
lags = lags * res;   % in ms

% Plot
%H1 = figure;
bar(lags,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])
hold on
burst_inx=length(ccr)/2-fix(burst_maxlag/res):length(ccr)/2;
int_burst=sum(ccr(burst_inx))/length(burst_inx);
bar(lags(burst_inx),ccr(burst_inx),'FaceColor','red')
[~,b_wind(1)]=min(abs(baseline_window(1)-lags));
[~,b_wind(2)]=min(abs(baseline_window(2)-lags));
int_baseline=sum(ccr((b_wind(1)):(b_wind(2))))/(b_wind(2)-b_wind(1)+1);
bar(lags(b_wind(1):b_wind(2)),ccr(b_wind(1):b_wind(2)),'FaceColor','blue')
Burst_index=int_burst/int_baseline;

function [burstNumNorm,Burst_ratio]=burst_detection(unitAct,maxBurstL,time_norm)
%

ISI = diff(unitAct);
burstInd = ISI<maxBurstL/1000;
burstWindow = diff([0;burstInd]);
burstNumNorm = sum(burstWindow==1) / time_norm;

    Bursts = [];
    w1=find(burstWindow==1);
    w2=find(burstWindow==-1);
    for b= 1:length(w1)-1
        Bursts = [Bursts, unitAct(w1(b):w2(b))'];
    end
    Burst_ratio = length(Bursts) / length(unitAct);
    