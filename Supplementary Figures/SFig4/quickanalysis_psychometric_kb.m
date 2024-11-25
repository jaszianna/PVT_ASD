%%This is a modified code of Balázs Hangya's lab to generate Supplementary Figure 4

function quickanalysis_psychometric_kb(animalNO,sessionID,sessionspec,protocoltag)
%based on: quickanalysis_pavlovian2_p(animalNO,sessionID,sessionspec,protocoltag)
%QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC) performs the analysis
%   for a session specified by the first two input arguments. SESSIONSPEC
%   should be a 1x3 logical array indicating the need for behavior,
%   recording and stimulation analysis.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC,PROTOCOLTAG) accepts a
%   PROTOCOLTAG argument to allow calls to trial event conversion programs
%   for different versions of the behavioral protocol.
%
%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.

% Input argument check
error(nargchk(0,4,nargin))

% Behavior, recording or both
if nargin < 4
    protocoltag = '';
end
if nargin < 3
    isbeh = 1;
    isrec = 1;
    isstim = 0;
else
    isbeh = sessionspec(1);
    isrec = sessionspec(2);
    isstim = sessionspec(3);
end

% Animal, session
if nargin < 2
    sessionID = 'error';
end
if nargin < 1
    animalID2 = 'error';
    animalID = 'error';
else
    %animalID2 = [num2str(animalNO)]; %ez nem kell
    animalID = [num2str(animalNO)];
end

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
%convert_events(fullpth);
% Stop if error
dbstop if error

% % Directories
% global DATAPATH
% resdir = [DATAPATH 'C:\Users\balazsfi.diana\Documents\MATLAB\Cellbase\CRT5\170412a'];
% if ~isdir(resdir)
%     mkdir(resdir)
% end

global uniqueSound
uniqueSound=zeros(1,5);

%% Create trial events structure
if isbeh
    if isempty(protocoltag)
        TE = solo2trialevents_psychometric_kb(([fullpth animalID sessionID '.mat']),1); %1: I want to save; based on solo2trialevents_auditory_cuedoutcome
    %TE = solo2trialevents4_auditory_gonogo([fullpth 'data_@auditory_gonogo_balazs_' animalID2 '_' sessionID '.mat']);
%     else
%          TE = solo2trialevents_3csrtt_extra_bd(([fullpth animalID sessionID '.mat']),1);
% %         evalstr = ['TE = solo2trialevents4_pavlovian2_Sanchari_' protocoltag '([fullpth ''data_@auditory_pavlovian2_' protocoltag '_Sanchari_'' animalID2 ''_'' sessionID ''.mat'']);'];
% %         eval(evalstr)
%     end
%     if isrec
        MakeTrialEvents2_psychometric_kb(fullpth)  % synchronize
    end
end

%% Update CellBase
if isrec
    
    addnewcells('dir',[animalID filesep sessionID]);
    %addnewcells('dir',[animalID filesep sessionID ]);
    %addnewcells('dir','CRT5\170405a\');
    cellids = findcell('mouse',animalID,'session',sessionID);
    disp(cellids)
    
end

%% Align spikes to defined behaviors
%if isbeh && isrec

% Prealign spikes for trial events
if isempty(protocoltag)
    problem_behav_cellid = [];
    for iC = 1:length(cellids),
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_psychometric_kb,'filetype','event','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_behav_cellid = [problem_behav_cellid cellid];
        end
    end
else
    problem_behav_cellid = [];
    for iC = 1:length(cellids),
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_psychometric_kb,'filetype','event','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_behav_cellid = [problem_behav_cellid cellid];
        end
    end
end


% %Is predictive?
% for k = 1:length(cellids)
%     H = figure;
%         viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent', 'StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
%     viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#TrialType','window',[-5 5])
%     maximize_figure(H)
% 
%     cellidt = cellids{k};
%     cellidt(cellidt=='.') = '_';
%     fnm = [resdir cellidt '_IPD.jpg'];   % save
%     saveas(H,fnm)
%     close(H)
%      end

mkdir(fullpth,'picturespertrials')
mkdir(fullpth,'stim')
mkdir(fullpth,'stimb')

%% Figures per trials
% 
 if isstim
   %lightpsth(fullpth)
 end
 
  if isbeh && isrec

      


for k = 1:length(cellids)
    G = figure;
    pause(0.01)
    viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Type','window',[-3 3])
    maximize_figure(G)
    
    cellidt = cellids{k};
    cellidt(cellidt=='.') = '_';
    fnm = [fullpth 'picturespertrials' '\' cellidt '_DF.jpg'];   % save
    saveas(G,fnm)
    close(G)
end

for k = 1:length(cellids)
    G = figure;
    pause(0.01)
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Type3','window',[-3 3])
    maximize_figure(G)
    
    cellidt = cellids{k};
    cellidt(cellidt=='.') = '_';
    fnm = [fullpth 'picturespertrials' '\' cellidt '_SO.jpg'];   % save
    saveas(G,fnm)
   
    close(G)
end


end
%% Light effects
if isrec && isstim
    
    % Create stimulus events
    %MakeStimEvents2_psychometric_kb_2(fullpth,7,'BurstStartNttl',4)
    MakeStimEvents2_psychometric_kb_2(fullpth,7,'BurstStartNttl',4)
%     Prealign spikes to stimulus events
    problem_stim_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0)
            %prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stimb','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_stim_cellid = [problem_stim_cellid cellid];
        end
    end
    
    % View light-triggered raster and PSTH
    TrigEvent = 'PulseOn';
    SEvent = 'BurstOff';
    win = [-4 4]; %%gatlos
    %win = [-0.2 0.5]; %%serkentos
     %parts = 'all'; %% ez ha nem kell kiszedni (serkentos)
    parts = '#StimSerieNum2'; %% ez szedi ki a pulzusokat a stimulusok kozotti elso nagy szunet elott
    %parts = '#StimSerieNum';
    dt = 0.01; %%gatlos
    %dt = 0.001; %% serkentes
    sigma = 0.02; %%gatlos
    %sigma = 0.002; %%serkentes
    PSTHstd = 'on';
    ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    for iCell = 1:length(cellids)
        for stimtype={'stim'}
        cellid = cellids(iCell);
        H = figure;
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
            'FigureNum',H,'eventtype',stimtype,'window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
            'EventMarkerWidth',0,'PlotZeroLine','off')
        maximize_figure(H)
        
        cellidt = cellid{1};
        cellidt(cellidt=='.') = '_';
        fnm = [fullpth stimtype{1} '\' cellidt '_optotagging.jpg'];   % save
        saveas(H,fnm)
        print(H,[fullpth stimtype{1} '\' cellidt '_optotagging'],'-dsvg');
        close(H)
        
         win = [-1,3];
         parts = '#StimSerieNum2';
        try
         [psth spsth spsth_se tags spt] = ultimate_psth(cellid{1},'stim','PulseOn',win,'dt',0.01,'sigma',0.01,'parts',parts);
         spt = spt{1};
         for i = 1:win(2)-win(1)
            S(i)= sum(sum(spt(:,(i-1)*(1/dt)+1:(i)*(1/dt)+1)))    
         end
        H = figure;
        fnm = [fullpth stimtype{1} '\' cellidt '_optotaggingbar.jpg'];   % save
        
        bar(S)
        saveas(H,fnm)
        close(H)
        end
        end
    end
end

% Cluster quality
global CheckCluster_DrawAllWaves
if isrec
    CheckCluster_DrawAllWaves = true;
    BatchSessionClust(fullpth)
    CheckCluster_DrawAllWaves = false;
end