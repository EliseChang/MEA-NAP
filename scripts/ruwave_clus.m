% Directories
homeDir = ("D:\MATLAB\MEA-NAP");
cd(homeDir)
fileDir = 'E:\spikeSorting\inputFiles';
addpath(genpath(fileDir))
addpath(genpath('\Functions\wave_clus'))
csvOutputDir = 'C:\Users\elise\Python\ReservoirComputing\data\connectivity';
matOutputDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData05Jan2023';

% Get all file names
metadataSpreadsheet = 'MEC.xlsx'; % file name
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
xlSheet = 'WT_KO_14-35';
xlRange = 'A2:C55';
[~,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1); % name of sample

unitsN = zeros(1, length(samples));

% Batch run channels
for n = 1:length(samples)
    sampleFolder = samples{n};
    cd(fullfile(fileDir, sampleFolder))

    % Make .txt file containing .mat files for each channel
    inputFiles = cell(60, 1);
    outputFiles = cell(60, 1);
    for c = 1:60
        inputFiles{c} = strcat(sampleFolder, '_channel_', num2str(c), '.mat');
        outputFiles{c} = strcat('times_', sampleFolder, '_channel_', num2str(c), '.mat');
    end
    writecell(inputFiles, 'files.txt')
    
    % Run spike sorting
    Do_clustering('files.txt', 'parallel', true, 'make_plots', false, 'save_spikes', false)
    cd(homeDir)

    % Get spike times
    currClus = 1;
    spikeTimes = cell(1, currClus);
    firingRates = [];
    channelClusterIds = [];
    for c = 1:60
        cluster_class = load(outputFiles{c}, 'cluster_class').cluster_class;
        clusClass = cluster_class(:,1);
        clusTimes = cluster_class(:,2);
        clusN = length(unique(clusClass));
        if any(clusClass == 0)
            clusIDs = [1:clusN] - 1;
        else
            clusIDs = 1:clusN;
        end
        
        for clus = clusIDs
    
            timesStruct = struct;
            times = sort(clusTimes(find(clusClass == clus)) / 1000); % convert from ms to s
            rate = length(times) / 600;
            timesStruct.mergedSorted = times;
            spikeTimes{1, currClus} = timesStruct;
            firingRates = [firingRates, rate];
            currClus = currClus + 1;
            
        end
        channelClusterIds = [channelClusterIds, repelem(c, clusN)];
    end

    nUnits = currClus - 1;    
    
    % Get adjacency matrix
    [adjM, ~] = adjM_thr_parallel(spikeTimes, 'mergedSorted', 10, 0.05, 25e3,...
                600, 200);
    density = density_und(adjM);
    nDisconnectedNodes = sum(sum(adjM, 1) == 0);
    if nDisconnectedNodes
        idxDiconnectedNodes = find(sum(adjM, 1) == 0);
    else
        idxDiconnectedNodes = NaN;
    end

for n = 1:length(samples)

    sampleFolder = samples{n};

    load(fullfile('D:\MATLAB\MEA-NAP\outputs\OutputData05Jan2023', sampleFolder));
    
    % Calculate firing rates
    firingRates = zeros(1, nUnits);
    for u = 1:nUnits
        firingRates(u) = length(spikeTimes{u}.mergedSorted) / 600;    
    end
    
    % Add inhibitory weights
    eiRatio = 1;
    inhibStrength = 2;
    adjMInhib = adjM;
    inhibUnits = 1:floor(nUnits / eiRatio); % Create 1 inhibitory unit for every 4 excitatory units
    inhibUnitsIdx = inhibUnits + nUnits; % index of inhibitory unit in the new matrix
    inhibPoolIdx = repelem(inhibUnits, eiRatio); % Assign excitatory units to inhibitory pools
    inhibPoolIdx(end+1:nUnits) = inhibUnits(end); % Assign any leftover units to the final pool
    % Get mean weight and activity for each unit and pool
    avgUnitWeights = mean(adjM, 1); % Will set E -> I weights as average for each unit
    unitActivities = avgUnitWeights .* firingRates; % Weight activity according to average weight for each unit
    poolActivities = accumarray(inhibPoolIdx', unitActivities'); % Summed activity in each pool
    maxPoolActivity = max(poolActivities) * max(adjM, [], 'all'); % Ensures max I -> E weight = max E -> E weight * I strength
    
    for i = inhibUnits(end)
        unitIdx = inhibUnitsIdx(i);
        pool = find(inhibPoolIdx == i);
        poolSize = numel(pool);
        adjMInhib(pool, unitIdx) = avgUnitWeights(pool);
        adjMInhib(unitIdx, pool) = - poolActivities(i) / maxPoolActivity * inhibStrength;
        adjMInhib(unitIdx, inhibUnitsIdx) = poolActivities(i) / maxPoolActivity;
        adjMInhib(unitIdx, unitIdx) = 0;
    end

    writematrix(adjMInhib, fullfile(csvOutputDir, strcat(sampleFolder, '_spike_sorted.csv')))
%     writematrix(channelClusterIds, fullfile(csvOutputDir, strcat(sampleFolder, '_spike_sorted_unit_ids.csv')))
%     save(fullfile(matOutputDir, sampleFolder), "spikeTimes", "channelClusterIds", "adjM", "adjMInhib", "nUnits", "density")

end
