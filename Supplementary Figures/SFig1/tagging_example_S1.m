%%This code generates the representative example of Supplementary Figure 1

cellid = 'A099_210602w_4.2';
% View light-triggered raster and PSTH
TrigEvent = 'PulseOn';
SEvent = 'BurstOff';
win = [-1 15]; %%gatlos
%win = [-0.03 0.03]; %%serkentos
parts = 'all'; %% ez ha nem kell kiszedni (serkentos)
%parts = '#StimSerieNum2'; %% ez szedi ki a pulzusokat a stimulusok kozotti elso nagy szunet elott
%parts = '#StimSerieNum';
dt = 0.01; %%gatlos
%dt = 0.001; %% serkentes
%sigma = 0.02; %%gatlos
sigma = 0.001; %%serkentes
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
stimtype = 'stim';
H = figure;
viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',H,'eventtype',stimtype,'window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','off')
maximize_figure(H)