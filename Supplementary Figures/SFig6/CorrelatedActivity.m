%%This code is a modified version of the code of Balázs Hangya's lab
%%This code was created to compute correlated activity of the given unit

function CorrelatedActivity()
%inputs:
path='X:\_Projects_2\Midline_anxiety\DATA\_Tetrode\ANALÍZIS\Cellbase_data\A090\210303w\';
resultspath='X:\_Projects_2\Midline_anxiety\DATA\_Tetrode\ANALÍZIS\Cellbase_data';
session='SD_A090_A091_210303';
 
awake_off1=2642; %in seconds
nest_on1=2642; %in seconds
nest_off1=3155; %in seconds
sleep_on1=3155; %in seconds
sleep_off1=3755; % in seconds
awake_off2=126; %in seconds
nest_on2=126; %in seconds
nest_off2=366;%in seconds
sleep_on2=366; %in seconds
sleep_off2=966; % in seconds


%------------------------------------------
ccg_wind=200; %in ms
ccg_res=1; % in ms
coincidence_bins = 5; % 
shuff_min = 500;     % in ms
shuff_max = 5000; %in ms
cutoff = 0; % 4 Hz
num_shuff = 5000;
corr_trsh =30; %in ms (treshold for correlation peak width


fileList = dir([path,session,'\','TT*.mat']);
[~, timestamps, ~] = load_open_ephys_data([path,session,'\','110_CH19.continuous']);
starttime=timestamps(1);
[~, stims, ~] = load_open_ephys_data([path,session,'\','all_channels.events']);
stims=stims(1:2:end)-starttime;
borders=find(diff(stims)>1);
borders=[0;borders;length(stims)];

mkdir([path,session,'\ccgs'])

jk=1;
for i=1:length(fileList)
%for i=1:1
   
    TS=load([path,session,'\',fileList(i).name]);
    TS=TS.TS;
    TS=(TS)/10000-starttime;
    for q=1:length(borders)-1
        TS(TS>(stims(borders(q)+1)) & TS<(stims(borders(q+1))))=[];
    end
    endtext=strfind(fileList(i).name,'_');
    dottext=strfind(fileList(i).name,'.');
    Tetrode=str2double(fileList(i).name(3:endtext-1));
    cellnum=str2double(fileList(i).name(endtext+1:dottext-1));
    listofpairs=[];
    if Tetrode<9
    awake_off=awake_off1;
    nest_on=nest_on1;
    nest_off=nest_off1;
    sleep_on=sleep_on1;
    sleep_off=sleep_off1;
    for j=1:length(fileList)
        endtext_j=strfind(fileList(j).name,'_');
        Tetrode_j=str2double(fileList(j).name(3:endtext_j-1));
        dottext_j=strfind(fileList(j).name,'.');
        cellnum_j=str2double(fileList(j).name(endtext_j+1:dottext_j-1));
        
        if ((Tetrode_j > Tetrode) | ((Tetrode_j == Tetrode) & (cellnum_j > cellnum))) & Tetrode_j<9
        listofpairs=[listofpairs, j];
        end
    end
    elseif Tetrode>8
    awake_off=awake_off2;
    nest_on=nest_on2;
    nest_off=nest_off2;
    sleep_on=sleep_on2;
    sleep_off=sleep_off2;
    for j=1:length(fileList)
        endtext_j=strfind(fileList(j).name,'_');
        Tetrode_j=str2double(fileList(j).name(3:endtext_j-1));
        dottext_j=strfind(fileList(j).name,'.');
        cellnum_j=str2double(fileList(j).name(endtext_j+1:dottext_j-1));
        
        if ((Tetrode_j > Tetrode) | ((Tetrode_j == Tetrode) & (cellnum_j > cellnum))) & Tetrode_j>8
        listofpairs=[listofpairs, j];
        end
    end
    end
    awake_filt=TS<nest_on;
    if ~isnan(nest_on)
    nesting_filt=TS>nest_on & TS<nest_off;
    end
    sleep_filt=TS>sleep_on & TS<sleep_off;
   
    for k=1:length(listofpairs)
        
        TS_k=load([path,session,'\',fileList(listofpairs(k)).name]);
        TS_k=TS_k.TS;
        TS_k=(TS_k)/10000-starttime;
        for q=1:length(borders)-1
            TS_k(TS_k>(stims(borders(q)+1)) & TS_k<(stims(borders(q+1))))=[];
        end
        
        G=figure
        subplot(1,3,1)
        
        [iscorr_awake,ccr_awake,~]=ccorr(TS(TS<nest_on),TS_k(TS_k<nest_on),ccg_wind,ccg_res,shuff_min,shuff_max,cutoff,num_shuff,corr_trsh);
        title('awake')
        xlabel('lag (ms)')
        subplot(1,3,2)
        [iscorr_nest,ccr_nest,~]=ccorr(TS(TS>nest_on & TS<nest_off),TS_k(TS_k>nest_on & TS_k<nest_off),ccg_wind,ccg_res,shuff_min,shuff_max,cutoff,num_shuff,corr_trsh);
        title('nesting')
        xlabel('lag (ms)')
        subplot(1,3,3)
        [iscorr_sleep,ccr_sleep,~]=ccorr(TS(TS>sleep_on & TS<sleep_off),TS_k(TS_k>sleep_on & TS_k<sleep_off),ccg_wind,ccg_res,shuff_min,shuff_max,cutoff,num_shuff,corr_trsh);
        title('sleep')
        xlabel('lag (ms)')
        
        
        %suptitle([session,' ',fileList(i).name(1:end-4) '&' ,fileList(listofpairs(k)).name(1:end-4),' CCGs'])
        maximize_figure(G);
        savefig([path,session,'\ccgs\',session,'_',fileList(i).name(1:end-4),'&',fileList(listofpairs(k)).name(1:end-4),'_CCG.fig'])
        saveas(gcf,[path,session,'\ccgs\',session,'_',fileList(i).name(1:end-4),'&',fileList(listofpairs(k)).name(1:end-4),'_CCG.png'])
        saveas(gcf,[path,session,'\ccgs\',session,'_',fileList(i).name(1:end-4),'&',fileList(listofpairs(k)).name(1:end-4),'_CCG.svg'])
        
        CCG_results(jk).cell1=fileList(i).name(1:end-4);
        CCG_results(jk).cell2=fileList(listofpairs(k)).name(1:end-4);
        CCG_results(jk).iscorr_awake=iscorr_awake;
        CCG_results(jk).iscorr_nest=iscorr_nest;
        CCG_results(jk).iscorr_sleep=iscorr_sleep;
        CCG_results(jk).coincidence_awake=sum(ccr_awake(ccg_wind+1-coincidence_bins:ccg_wind+1+coincidence_bins))/awake_off;
        CCG_results(jk).coincidence_nest=sum(ccr_nest(ccg_wind+1-coincidence_bins:ccg_wind+1+coincidence_bins))/(nest_off-nest_on);
        CCG_results(jk).coincidence_sleep=sum(ccr_sleep(ccg_wind+1-coincidence_bins:ccg_wind+1+coincidence_bins))/(sleep_off-sleep_on);
        CCG_results(jk).ccr_awake=ccr_awake;
        CCG_results(jk).ccr_nest=ccr_nest;
        CCG_results(jk).ccr_sleep=ccr_sleep;
        
        jk=jk+1;
        
    
    
    end
end
save([path,session,'\ccgs\CCG_results'],'CCG_results');
close all       
    

function [iscorr,ccr, lags] = ccorr(ncc1,ncc2,wn,res,shuff_min,shuff_max,cutoff,num_shuff,corr_trsh) 


% Calculate spike times in milliseconds
sr = 1000;
nc1 = ncc1 * sr;
mn1 = nc1(1);  % only relative spike times count; avoid out of memory
nc1 = nc1 - mn1;
nc1(nc1<0.5*res) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
wn2 = wn / 1000;    % window size in seconds

% Calculate spike times in milliseconds
nc2 = ncc2 * sr;
nc2 = nc2 - mn1;
nc2(nc2<0.5*res) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39

% Auto-correlogram
czunit1 = zeros(1,round(nc1(end)/res)+5);
zunit1(round(nc1/res)) = 1;
zunit2 = zeros(1,round(nc2(end)/res)+5);
zunit2(round(nc2/res)) = 1;
% [ccr, lags] = xcorr(zunit1,zunit2,wn2*sr/res);     % 1->2; window: -wn ms - wn ms
% %ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
% %lags(length(lags)/2+0.5) = [];
% lags = lags * res;   % in ms
% 
% % Plot
% %H1 = figure;
% bar(lags,ccr,'FaceColor','black')
% set(gca,'XLim',[-wn wn])



shuff_min = shuff_min * sr / 1000;
shuff_max = shuff_max * sr / 1000;

shf = round(rand(1,shuff_max)*num_shuff+shuff_min);
[ccr,lags,rccg] = xcorr_wrand_filter(sr,cutoff,zunit1,zunit2,wn2*sr/res,shf);     % 1->2; window: -wn ms - wn ms
% ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
% if isequal(ncc1,ncc2)
%     ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
%     rccg(:,end) = [];
% end

% Plot

time = linspace(-wn,wn,length(ccr));
bar(time,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])


%plot(time,ccr/sum(ccr)*length(ccr));

% Plot confidence interval
hold on
[lc nc] = size(rccg);
ptc = ceil(lc*0.0005);
upr = zeros(1,nc);
lwr = zeros(1,nc);
for k = 1:nc
    sts = sort(rccg(:,k),'ascend');
    upr(k) = sts(end-ptc);
    lwr(k) = sts(ptc);
end
plot(time,upr,'Color',[0.7 0.7 0.7])
plot(time,lwr,'Color',[0.7 0.7 0.7])

j=0;
k=0;
for i=1:length(ccr)
    if ccr(i)>upr(i)
        j=j+1;
    else
        j=0;    
    end
    
    if j>k
        k=j;
    end
end

if k*res>corr_trsh
    iscorr=1;
else
    iscorr=0;
end

