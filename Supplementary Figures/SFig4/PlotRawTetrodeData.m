%%This code plot the raw LFP data of the given tetrode and marks with some colour the spikes belong to the given neurons

figure
numberofcells = 3;
[data, timestamps, ~]=load_open_ephys_data(['D:\Data\Tetrode\A070\2021-01-06_15-40-36\110_CH',mat2str(1),'.continuous']);       
    for k = 1:numberofcells
    load(['D:\Data\Tetrode\A070\2021-01-06_15-40-36\Continuous_Data\Continuous_Datasc.GUI\TT1_',mat2str(k),'.mat']);      
    APs{k} = nan(length(data),1);
    TS2{k} = round((TS-timestamps(1)*10000)*3);
    sr = 30000;
    end
    
    
for j = 1:4
    [data, timestamps, ~]=load_open_ephys_data(['D:\Data\Tetrode\A070\2021-01-06_15-40-36\Continuous_Data\110_CH',mat2str(0+j),'.continuous']);   
    [b,a] = butter(3,150/(sr/2),'high');
    data2 = filter(b,a,data);
    for k=1:numberofcells
    for i = 1:length(TS2{k})
        APs{k}(TS2{k}(i)-6:TS2{k}(i)+6) = data2(TS2{k}(i)-6:TS2{k}(i)+6);
    end
    end
    g(j)=subplot(4,1,j)
    plot(data2)
    hold on
    for k=1:numberofcells
        plot(APs{k})
    end
end
linkaxes(g);