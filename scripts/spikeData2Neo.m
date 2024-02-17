% Directories
homeDir = ("D:\MATLAB\MEA-NAP");
cd(homeDir)
spikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData09Nov2022\1_SpikeDetection\1A_SpikeDetectedData';
addpath(spikeDir)

outputDir = 'C:\Users\elise\Python\ReservoirComputing\data\spike_times';
addpath(outputDir)

% Get file names
metadataSpreadsheet = 'MEC.xlsx'; % file name
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
addpath(outputDir)
xlSheet = 'Sheet1';
xlRange = 'A2:A129';
[~,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1); % name of sample

channelsN = 60;
methods = {'bior1p5','bior1p3','db2','thr5'};

% wave_clus
for n = 1:length(samples)

    disp(samples{n})
    spikeTimes = load(fullfile(spikeDir,strcat(samples{n},'_spikes')),"-mat","spikeTimes").spikeTimes; % in s
    mergedSpikeTimes = cell(channelsN, 1);
    for ch = 1:channelsN

        % Collate spike waveforms, merge spike times and convert times from s to ms
        times = spikeTimes{1, ch};
        mergedTimes = [];
        for m = 1:numel(methods)
            method = methods{m};
            methodTimes = round(times.(method) * 1000); % convert to ms and round to nearest ms -- spikes detected by
                                                        % different methods within 1 ms of each other will be considered
                                                        % identical
            mergedTimes = union(mergedTimes, methodTimes, "stable");
        end
        mergedSpikeTimes{ch, 1} = mergedTimes;
        clear times mergedTimes
    end
    
    allSpikeTimes = cell2mat(mergedSpikeTimes);
    totalSpikesN = numel(allSpikeTimes);
    spikeData = zeros(totalSpikesN + 1, 2);
    lastIndex = 0;
    for ch = 1:channelsN
        chSpikesN = numel(mergedSpikeTimes{ch, 1});
        % fill concatenated vector of spike times with channel/unit index
        spikeData(lastIndex + 1:lastIndex + 1 + chSpikesN, 2) = ch - 1; % for Python 0-indexing
        lastIndex = lastIndex + chSpikesN;
    end
    spikeData(end, :) = []; % because of indexing when filling, extra row added -- delete
    spikeData(:, 1) = allSpikeTimes;
    writematrix(spikeData, fullfile(outputDir, strcat(samples{n}, '.csv')))
    
end