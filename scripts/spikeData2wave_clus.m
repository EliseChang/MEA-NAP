% Directories
homeDir = ("D:\MATLAB\MEA-NAP");
cd(homeDir)
spikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData09Nov2022\1_SpikeDetection\1A_SpikeDetectedData';
addpath(spikeDir)

outputDir = 'E:\spikeSorting\inputFiles'; % for wave_clus: 'E:\spikeSorting\inputFiles'
addpath(outputDir)

% Get file names
metadataSpreadsheet = 'MEC.xlsx'; % file name
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
addpath(outputDir)
xlSheet = 'WT_KO_14-35';
xlRange = 'A49:A55';
[~,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1); % name of sample

sr = 25000; % sampling rate
channelsN = 60;
methods = {'bior1p5','bior1p3','db2','thr5'};

% wave_clus
for n = 1:length(samples)

    disp(samples{n})

    mkdir(fullfile(outputDir,samples{n})); % folder to store each channel file

    spikeWaveforms = load(fullfile(spikeDir,strcat(samples{n},'_spikes')),"-mat","spikeWaveforms").spikeWaveforms;
    spikeTimes = load(fullfile(spikeDir,strcat(samples{n},'_spikes')),"-mat","spikeTimes").spikeTimes; % in s
    mergedSpikeTimes = cell(1, channelsN);
    for ch = 1:channelsN

        % Collate spike waveforms, merge spike times and convert times from s to ms
        waveforms = spikeWaveforms{1, ch};
        times = spikeTimes{1, ch};
        spikes = zeros(1, 51);
        index = 0;
        for m = 1:numel(methods)
            method = methods{m};
            spikes = union(spikes, waveforms.(method), "rows", "stable");
            index = union(index, times.(method), "stable");
        end
        spikes(1,:) = [];
        index(1) = [];
        index = index*1000;
        clear waveforms times
        
        % Check spikes and index are same length
        assert(size(spikes,1) == length(index), "Number of spike waveforms and times are not equal.")
        
        % for wave_clus
        chFileName = strcat(samples{n},"_channel_",num2str(ch),".mat");
        save(fullfile(outputDir,samples{n},chFileName),"spikes","index","sr")
        

%         % for Python to convert to Neo
%         mergedSpikeTimes{1, ch} = index;
%         clear spikes index
    end
    
%     allSpikeTimes = cell2mat(mergedSpikeTimes);
%     totalSpikesN = numel(allSpikeTimes);
%     spikeData = zeros(spikesN, 2);
%     lastIndex = 0;
%     for ch = 1:channelsN
%         chSpikesN = numel(mergedSpikeTimes{1, ch});
%         spikeData(lastIndex:lastIndex + chSpikesN, 2) = ch;
%     end
%     spikeData(:, 1) = allSpikeTimes;
%     writematrix(spikeData, fullfile(outputDir, strcat(samples{n}, '.csv')))
    
end