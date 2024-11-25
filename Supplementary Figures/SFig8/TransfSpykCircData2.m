%%This code transform the output from the SpykingCircus Clustering Program to Matlab variables

amplitudes = double(readNPY([phyfolder filesep 'amplitudes.npy']));

spiketimes = double(readNPY([phyfolder filesep 'spike_times.npy']));

phycluster = double(readNPY([phyfolder filesep 'spike_clusters.npy']));

spiketempl = double(readNPY([phyfolder filesep 'spike_templates.npy']));

templates = double(readNPY([phyfolder filesep 'templates.npy']));

 

clusterfile = [phyfolder filesep 'cluster_info.tsv'];

tsvclusters = readtable(clusterfile,'FileType','text');

goodclusters = tsvclusters.cluster_id(cellfun(@(x) contains(x,'good'),tsvclusters.group));

goodclusters_phyChannel = tsvclusters.ch(cellfun(@(x) contains(x,'good'),tsvclusters.group))+1;

goodclusters_depth = tsvclusters.depth(cellfun(@(x) contains(x,'good'),tsvclusters.group));

goodclusters_amp = tsvclusters.amp(cellfun(@(x) contains(x,'good'),tsvclusters.group));

 

phyclust = cell(length(filenames),1);

phytimes = cell(length(filenames),1);

phytemplates = cell(length(filenames),1);

for k = 1:length(filenames)

    phyclust{k,1} = phycluster(spiketimes>=filestart(k) & spiketimes<=fileend(k));

    phytimes{k,1} = spiketimes(spiketimes>=filestart(k) & spiketimes<=fileend(k));

    phytemplates{k,1} = spiketempl(spiketimes>=filestart(k) & spiketimes<=fileend(k));

    if k>1

        phytimes{k,1} = phytimes{k,1}-fileend(k-1);

    end

end

 

spikestemp = cell(length(filenames),length(goodclusters));

wavetemp = cell(length(filenames),length(goodclusters));

for f = 1:length(filenames)

    for c = 1:length(goodclusters)

        spikestemp{f,c} = phytimes{f,1}(phyclust{f,1} == goodclusters(c));

        template = mode(phytemplates{f,1}(phyclust{f,1} == goodclusters(c)));

        if ~isnan(template)

            wavetemp{f,c} = templates(template+1,:,goodclusters_phyChannel(c));

        else

            wavetemp{f,c} = NaN(1,61);

        end

    end

end

spikes = cell2table(spikestemp,'VariableNames',arrayfun(@(x) ['clu_' num2str(x)], goodclusters,'UniformOutput',false),'RowNames',filenames);

waveforms = cell2table(wavetemp,'VariableNames',arrayfun(@(x) ['clu_' num2str(x)], goodclusters,'UniformOutput',false),'RowNames',filenames);