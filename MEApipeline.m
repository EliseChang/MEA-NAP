% Process data from MEA recordings of 2D and 3D neuronal cultures
% Created: RC Feord, May 2021
% Authors: T Sit, RC Feord, AWE Dunn, J Chabros and other members of the Synaptic and Network Development (SAND) Group
%% USER INPUT REQUIRED FOR THIS SECTION
% In this section all modifiable parameters of the analysis are defined.
% No subsequent section requires user input.
% Please refer to the documentation for guidance on parameter choice here:
% https://analysis-pipeline.readthedocs.io/en/latest/pipeline-steps.html#pipeline-settings
% clear all 
close all
restoredefaultpath
% Directories
HomeDir = '/Users/elise/MEA-NAP'; % Where the MEA-NAP (MEA Network Analysis Pipeline) code is located
Params.outputDataFolder = fullfile(HomeDir, 'outputs');   % Where to save the output data, leave as '' if same as HomeDir 
rawData = fullfile(HomeDir, 'mat', 'MeCP2DIV240FebMar2024Orgs');  % path to raw data .mat files
% Added following global parameters to process stimulation data
Params.stim = true; Params.rawDataFieldName = 'postSALPA0_500ms'; % if true, spike detection will require file containing 'stimDat' struct with data matrix stored in 'fieldName', otherwise 'dat' used
if ~Params.stim
    Params.rawDataFieldName = '';
end
Params.priorAnalysisPath = '';  % path to prev analysis, leave as [''] if no prior anlaysis
spikeDetectedData = ''; % fullfile(Params.priorAnalysisPath, '1_SpikeDetection','1A_SpikeDetectedData'); % path to spike-detected data, leave as '' if no previously detected spike data

% Input and output filetype
spreadsheet_file_type = 'excel'; % 'csv' or 'excel'
spreadsheet_filename = '/Users/elise/MEA-RC/metadata/MeCP2DIV240FebMar2024Orgs.xlsx'; 
sheet = 'MEA-RC'; % 'BaselineMEA-NAP'  specify excel sheet
xlRange = 'A2:D11'; % specify range on the sheet (e.g., 'A2:C7' would analyse the first 6 files)
csvRange = [2, Inf]; % read the data in the range [StartRow EndRow], e.g. [2 Inf] means start reading data from row 2
Params.output_spreadsheet_file_type = 'csv';  % .xlsx or .csv

% Analysis step settings
Params.priorAnalysisDate = ''; % prior analysis date in format given in output data folder e.g., '27Sep2021'
Params.priorAnalysis = 0; % use previously analysed data? 1 = yes, 0 = no
Params.startAnalysisStep = 1; % if Params.priorAnalysis=0, default is to start with spike detection
Params.optionalStepsToRun = {''}; % include 'generateCSV' to generate csv for rawData folder
                                  % include 'Stats' to look at feature
                                  % correlation and classification across groups
                                  % include 'combineDIVplots' to combine plots across DIVs

% Spike detection settings
detectSpikes = 1; % run spike detection? % 1 = yes, 0 = no
Params.runSpikeCheckOnPrevSpikeData = 0; % whether to run spike detection check without spike detection 
Params.fs = 25000; % Sampling frequency you selected when acquiring data, e.g., MCS: 25000, Axion: 12500;
Params.dSampF = 25000; % down sampling factor for spike detection check, e.g., can set the same as the sampling frequency
Params.potentialDifferenceUnit = 'uV';  % Unit for voltage signal, e.g., MCS: uV, Axion: V 
Params.channelLayout = 'MCS60';  % 'MCS60' (for MEA2100), 'Axion64' (for 6-well plates), or 'MCS60old' (for MEA6100)
Params.thresholds = {'4', '5'}; % standard deviation multiplier threshold(s), eg. {'2.5', '3.5', '4.5'}
Params.wnameList = {'bior1.5', 'bior1.3', 'db2'}; % wavelet methods to use e.g., {'bior1.5', 'bior1.3', 'mea'}; 
Params.costList = -0.12;
Params.SpikesMethod = 'mergedAll';  % wavelet methods, e.g., 'bior1p5', or 'mergedAll', or 'mergedWavelet'

% Grounding
Params.grdElecsIncluded = 1; % grounded electrodes are specified for each recording in 4th column of the spreadsheet as comma+space (, )-separated lists
Params.electrodesToGroundPerRecordingUseName = 1;  % 1 = electrode IDs given, 0 = electrodes indices given

% Functional connectivity inference settings
Params.FuncConLagval = [10, 25, 50]; % set the different lag values (in ms), default to [10, 15, 25]
Params.TruncRec = 0; % truncate recording? 1 = yes, 0 = no
Params.TruncLength = 120; % length of truncated recordings (in seconds)
Params.adjMtype = 'weighted'; % 'weighted' or 'binary'

% Connectivity matrix thresholding settings
Params.ProbThreshRepNum = 200; % probabilistic thresholding number of repeats (recommend at least 180) 
Params.ProbThreshTail = 0.05; % probabilistic thresholding percentile threshold 
Params.ProbThreshPlotChecks = 1; % randomly sample recordings to plot probabilistic thresholding check, 1 = yes, 0 = no
Params.ProbThreshPlotChecksN = 5; % number of random checks to plot

% Node cartography settings 
Params.autoSetCartographyBoudariesPerLag = 1;  % whether to fit separate boundaries per lag value
Params.cartographyLagVal = [10, 25, 50]; % lag value (ms) to use to calculate PC-Z distribution (only applies if Params.autoSetCartographyBoudariesPerLag = 0)
Params.autoSetCartographyBoundaries = 1;  % whether to automatically determine bounds for hubs or use custom ones

% Statistics and machine learning settings 
Params.classificationTarget = 'AgeDiv';  % which property of the recordings to classify 
Params.classification_models = {'linearSVM', 'kNN', 'fforwardNN', 'decisionTree', 'LDA'};
Params.regression_models = {'svmRegressor', 'regressionTree', 'ridgeRegression', 'fforwardNN'};

% Plot settings
Params.figExt = {'.png'};  % supported options are '.fig', '.png', and '.svg'
Params.fullSVG = 1;  % whether to insist svg even with plots with large number of elements
Params.showOneFig = 1;  % otherwise, 0 = pipeline shows plots as it runs, 1: supress plots

%% GUI / Tutorial mode settings 
Params.guiMode = 1;   % GUI mode? 1 = on, 0 = off
if Params.guiMode == 1
    runPipelineApp
    spikeDetectedData = Params.spikeDetectedData;
end 

%% Paths 
% add all relevant folders to path
cd(HomeDir)
addpath(genpath('Functions'))
addpath('Images')


%% END OF USER REQUIRED INPUT SECTION
% The rest of the MEApipeline.m runs automatically. Do not change after this line
% unless you are an expert user.
% Define output folder names
formatOut = 'ddmmmyyyy'; 
Params.Date = datestr(now,formatOut); 
clear formatOut

AdvancedSettings

if Params.runSpikeCheckOnPrevSpikeData
    fprintf(['You specified to run spike detection check on previously extracted spikes, \n', ... 
            'so I will skip over the spike detection step \n'])
    detectSpikes = 0;
end 

Params.detectSpikes = detectSpikes;  % As a record of option selection

% Allow starting from a subset of steps 
if length(Params.startAnalysisStep) > 1
    Params.startAnalysisSubStep = Params.startAnalysisStep(2);
    Params.startAnalysisStep = str2num(Params.startAnalysisStep(1));
else
    if isstr(Params.startAnalysisStep)
        Params.startAnalysisStep = str2num(Params.startAnalysisStep);
    end
    Params.startAnalysisSubStep = 'ALL';
end 

%% Optional step : generate csv 
if any(strcmp(Params.optionalStepsToRun,'generateCSV')) 
    fprintf('Generating CSV with given rawData folder \n')
    mat_file_list = dir(fullfile(rawData, '*mat'));
    name_list = {mat_file_list.name}';
    name_without_ext = {};
    div = {};
    for filenum = 1:length(name_list)
        name_without_ext{filenum} = name_list{filenum}(1:end-4);
        div{filenum} = name_list{filenum}((end-5):end-4);
    end 
    name = name_without_ext'; 
    div = div';
    name_table = table([name, div]);
    writetable(name_table, spreadsheet_filename)
end 

%% setup - additional setup
setUpSpreadSheet  % import metadata from spreadsheet
[~,Params.GrpNm] = findgroups(ExpGrp);
[~,Params.DivNm] = findgroups(ExpDIV);

% create output data folder if doesn't exist
CreateOutputFolders(HomeDir, Params.outputDataFolder, Params.Date, Params.GrpNm, Params)

% Set up one figure handle to save all the figures
oneFigureHandle = NaN;
oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);
    
% plot electrode layout 
plotElectrodeLayout(Params.outputDataFolder, Params, oneFigureHandle)

% export parameters to csv file
outputDataWDatePath = fullfile(Params.outputDataFolder, strcat('OutputData',Params.Date));
ParamsTableSavePath = fullfile(outputDataWDatePath, strcat('Parameters_',Params.Date,'.csv'));
writetable(struct2table(Params,'AsArray',true), ParamsTableSavePath)

% save metadata
metaDataSaveFolder = fullfile(outputDataWDatePath, 'ExperimentMatFiles');
for ExN = 1:length(ExpName)
    Info.FN = ExpName(ExN);
    Info.DIV = num2cell(ExpDIV(ExN));
    Info.Grp = ExpGrp(ExN);
    InfoSavePath = fullfile(metaDataSaveFolder, strcat(char(Info.FN),'_',Params.Date,'.mat'));
    save(InfoSavePath,'Info')
end

% create a random sample for checking the probabilistic thresholding
if Params.ProbThreshPlotChecks == 1
    Params.randRepCheckExN = randi([1 length(ExpName)],1,Params.ProbThreshPlotChecksN);
    Params.randRepCheckLag = Params.FuncConLagval(randi([1 length(Params.FuncConLagval)],1,Params.ProbThreshPlotChecksN));
    Params.randRepCheckP = [Params.randRepCheckExN;Params.randRepCheckLag];
end

%% Step 1 - spike detection

if ((Params.priorAnalysis == 0) || (Params.runSpikeCheckOnPrevSpikeData)) && (Params.startAnalysisStep == 1) 

    if (detectSpikes == 1) || (Params.runSpikeCheckOnPrevSpikeData)
        if iscell(rawData)
            for pathIdx = 1:length(rawData)
                addpath(rawData{pathIdx});
            end 
        else 
            addpath(rawData)
        end 
    else
        addpath(spikeDetectedData)
    end
    
    % If voltage data is stored in a structure, save the results in a folder named
    % after the field used in this structure
    if ~isempty(Params.rawDataFieldName)
        savePath = fullfile(Params.outputDataFolder, ...
                        strcat('OutputData', Params.Date), ...
                        '1_SpikeDetection', '1A_SpikeDetectedData', Params.rawDataFieldName);
        if ~isfolder(savePath)
            mkdir(savePath)
        end
    else
        savePath = fullfile(Params.outputDataFolder, ...
                        strcat('OutputData', Params.Date), ...
                        '1_SpikeDetection', '1A_SpikeDetectedData');
    end
    
    % Run spike detection
    if detectSpikes == 1
        batchDetectSpikes(rawData, savePath, option, ExpName, Params);
    end
    
    % Specify where ExperimentMatFiles are stored
    experimentMatFileFolder = fullfile(Params.outputDataFolder, ...
           strcat('OutputData',Params.Date), 'ExperimentMatFiles');

    % Plot spike detection results 
    for  ExN = 1:length(ExpName)

        if Params.runSpikeCheckOnPrevSpikeData
            spikeDetectedDataOutputFolder = spikeDetectedData;
        else
            spikeDetectedDataOutputFolder = savePath;
        end 

        spikeFilePath = fullfile(spikeDetectedDataOutputFolder, strcat(char(ExpName(ExN)),'_spikes.mat'));
        load(spikeFilePath,'spikeTimes','spikeDetectionResult','spikeWaveforms')

        experimentMatFilePath = fullfile(experimentMatFileFolder, ...
            strcat(char(ExpName(ExN)),'_',Params.Date,'.mat'));
        load(experimentMatFilePath,'Info')
    
        if ~isempty(Params.rawDataFieldName)
            spikeDetectionCheckGrpFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), '1_SpikeDetection', '1B_SpikeDetectionChecks', char(Info.Grp),...
                Params.rawDataFieldName);
        else
            spikeDetectionCheckGrpFolder = fullfile(Params.outputDataFolder, ...
            strcat('OutputData',Params.Date), '1_SpikeDetection', '1B_SpikeDetectionChecks', char(Info.Grp));
        end

        FN = char(Info.FN);
        spikeDetectionCheckFNFolder = fullfile(spikeDetectionCheckGrpFolder, FN);

        if ~isfolder(spikeDetectionCheckFNFolder)
            mkdir(spikeDetectionCheckFNFolder)
        end 

        plotSpikeDetectionChecks(spikeTimes, spikeDetectionResult, spikeWaveforms, Info, Params, spikeDetectionCheckFNFolder, oneFigureHandle)

        % Check whether there are no spikes at all in the recording 
        checkIfAnySpikes(spikeTimes, ExpName{ExN});

        clear spikeTimes spikeDetectionResult channels spikeWaveforms Info
    end

end

%% Step 2 - neuronal activity

if Params.priorAnalysis==0 || Params.priorAnalysis==1 && Params.startAnalysisStep<3
    fprintf('Running step 2 of MEA-NAP: neuronal activity \n')
    experimentMatFolderPath = fullfile(Params.outputDataFolder, ...
        strcat('OutputData',Params.Date), 'ExperimentMatFiles');
    priorAnalysisExperimentMatFolderPath = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');

    for ExN = 1:length(ExpName)
        experimentMatFname = strcat(char(ExpName(ExN)),'_',Params.Date,'.mat');
        experimentMatFpath = fullfile(experimentMatFolderPath, experimentMatFname);
        load(experimentMatFpath, 'Info')

        % extract spike matrix, spikes times and associated info
        disp(char(Info.FN))
        if Params.priorAnalysis==1 && Params.startAnalysisStep==2
            spikeDetectedDataFolder = spikeDetectedData;
        else
            if detectSpikes == 1
                if ~isempty(Params.rawDataFieldName)
                    spikeDetectedDataFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData', Params.Date), '1_SpikeDetection', ...
                    '1A_SpikeDetectedData',Params.rawDataFieldName);
                else
                    spikeDetectedDataFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData', Params.Date), '1_SpikeDetection', ...
                    '1A_SpikeDetectedData');
                end
            else
                spikeDetectedDataFolder = spikeDetectedData;
            end
        end

        channelLayout =  Params.channelLayout; %PerRecording{ExN};

        [spikeMatrix,spikeTimes,~,Info] = formatSpikeTimes(...
        char(Info.FN), Params, Info, spikeDetectedDataFolder, channelLayout);
    
        % Load LFP and save ground electrodes
        load(fullfile(spikeDetectedDataFolder, strcat(char(Info.FN), '_spikes.mat')), 'LFP', 'spikeDetectionResult') 
        Info.grdElecs = spikeDetectionResult.params.groundElecs;

        % initial run through to establish max values for scaling
        spikeFreqMax(ExN) = prctile((downSampleSum(full(spikeMatrix), Info.duration_s*10)),95,'all');
        infoFnFilePath = fullfile(experimentMatFolderPath, ...
                          strcat(char(Info.FN),'_',Params.Date,'.mat'));
        save(infoFnFilePath, 'Info', 'spikeTimes', 'spikeMatrix', 'LFP')
        clear spikeTimes spikeMatrix
    end

    % extract and plot neuronal activity

    % Set up one figure handle to save all the figures
    oneFigureHandle = NaN;
    oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);

    disp('Electrophysiological properties')

    spikeFreqMax = max(spikeFreqMax);

    for  ExN = 1:length(ExpName)

        disp(ExpName{ExN})
        
        experimentMatFname = strcat(char(ExpName(ExN)),'_',Params.Date,'.mat'); 
        experimentMatFpath = fullfile(experimentMatFolderPath, experimentMatFname);
        load(experimentMatFpath,'spikeTimes','spikeMatrix','Info','LFP');

        % get firing rates and burst characterisation
        Ephys = firingRatesBursts(spikeMatrix,Params,Info);

        idvNeuronalAnalysisGrpFolder = fullfile(Params.outputDataFolder, ...
            strcat('OutputData',Params.Date), '2_NeuronalActivity', ...
            '2A_IndividualNeuronalAnalysis', char(Info.Grp));

        if ~isfolder(idvNeuronalAnalysisGrpFolder)
            mkdir(idvNeuronalAnalysisGrpFolder)
        end 

        idvNeuronalAnalysisFNFolder = fullfile(idvNeuronalAnalysisGrpFolder, char(Info.FN));
        if ~isfolder(idvNeuronalAnalysisFNFolder)
            mkdir(idvNeuronalAnalysisFNFolder)
        end

        % generate and save raster plot
        % rasterPlot(char(Info.FN),spikeMatrix,Params,spikeFreqMax, idvNeuronalAnalysisFNFolder, oneFigureHandle)
        rasterPlotLFP(char(Info.FN),spikeMatrix,LFP,Params,idvNeuronalAnalysisFNFolder, oneFigureHandle)
        % electrode heat maps
        coords = Params.coords; %{ExN}
        % Electrode firing rate distribution plots
        electrodeHeatMaps(Info, Ephys, spikeFreqMax, Params, coords, idvNeuronalAnalysisFNFolder, oneFigureHandle)
        firingRateElectrodeDistribution(char(Info.FN), Ephys, Params, ... 
            Info, idvNeuronalAnalysisFNFolder, 'histogram', oneFigureHandle)
        firingRateElectrodeDistribution(char(Info.FN), Ephys, Params, ... 
            Info, idvNeuronalAnalysisFNFolder, 'halfViolin', oneFigureHandle)

        infoFnFilePath = fullfile(experimentMatFolderPath, ...
                          strcat(char(Info.FN),'_',Params.Date,'.mat'));
        save(infoFnFilePath,'Info','Params','spikeTimes','Ephys','-v7.3')

        clear spikeTimes spikeMatrix LFP

    end

    % create combined plots across groups/ages -- note: plots FR bounds
    PlotEphysStats(ExpName,Params,HomeDir, oneFigureHandle)
    saveEphysStats(ExpName, Params)
    cd(HomeDir)

end


%% Step 3 - functional connectivity, generate adjacency matrices

if Params.priorAnalysis==0 || Params.priorAnalysis==1 && Params.startAnalysisStep<4

    disp('generating adjacency matrices')
    
    % Set up one figure handle to save all the figures
    oneFigureHandle = NaN;
    oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);

    for  ExN = 1:length(ExpName)

        if Params.priorAnalysis==1 && Params.startAnalysisStep==3
            priorAnalysisExpMatFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
            spikeDataFname = strcat(char(ExpName(ExN)),'_',Params.priorAnalysisDate,'.mat');
            spikeDataFpath = fullfile(priorAnalysisExpMatFolder, spikeDataFname);
            load(spikeDataFpath, 'spikeTimes', 'Ephys', 'Info')
        else
            ExpMatFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), 'ExperimentMatFiles');
            spikeDataFname = strcat(char(ExpName(ExN)),'_',Params.Date,'.mat');
            spikeDataFpath = fullfile(ExpMatFolder, spikeDataFname);
            load(spikeDataFpath, 'Info', 'spikeTimes', 'Ephys') % Params
        end

        disp(char(Info.FN))

        adjMs = generateAdjMs(spikeTimes, ExN, Params, Info, oneFigureHandle);

        ExpMatFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), 'ExperimentMatFiles');
        infoFnFname = strcat(char(Info.FN),'_',Params.Date,'.mat');
        infoFnFilePath = fullfile(ExpMatFolder, infoFnFname);
        save(infoFnFilePath, 'Info', 'spikeTimes', 'Ephys', 'adjMs') % 'Params',
    end

end

%% Step 4 - network activity

if Params.priorAnalysis==0 || Params.priorAnalysis==1 && Params.startAnalysisStep<=4
    
    % Set up one figure handle to save all the figures
    oneFigureHandle = NaN;
    oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);
    
    % Step 4 Analysis step
    if strcmp(Params.startAnalysisSubStep, 'ALL')
        for  ExN = 1:length(ExpName) 

            if Params.priorAnalysis==1 && Params.startAnalysisStep==4
                priorAnalysisExpMatFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
                spikeDataFname = strcat(char(ExpName(ExN)),'_',Params.priorAnalysisDate,'.mat');
                spikeDataFpath = fullfile(priorAnalysisExpMatFolder, spikeDataFname);
                load(spikeDataFpath, 'spikeTimes', 'Ephys','adjMs','Info')
            else
                ExpMatFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData',Params.Date), 'ExperimentMatFiles');
                spikeDataFname = strcat(char(ExpName(ExN)),'_',Params.Date,'.mat');
                spikeDataFpath = fullfile(ExpMatFolder, spikeDataFname);
                load(spikeDataFpath, 'Info', 'spikeTimes', 'Ephys','adjMs') %'Params', 
            end

            disp(char(Info.FN))

            idvNetworkAnalysisGrpFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), '4_NetworkActivity', ...
                '4A_IndividualNetworkAnalysis', char(Info.Grp));

            idvNetworkAnalysisFNFolder = fullfile(idvNetworkAnalysisGrpFolder, char(Info.FN));
            if ~isfolder(idvNetworkAnalysisFNFolder)
                mkdir(idvNetworkAnalysisFNFolder)
            end 

            if Params.priorAnalysis == 1
                if isempty(spikeDetectedData)
                    spikeDetectedDataFolder = fullfile(Params.outputDataFolder, ...
                        strcat('OutputData', Params.Date), '1_SpikeDetection', ...
                        '1A_SpikeDetectedData');
                else 
                    spikeDetectedDataFolder = spikeDetectedData;
                end 
            else
                spikeDetectedDataFolder = fullfile(Params.outputDataFolder, ...
                        strcat('OutputData', Params.Date), '1_SpikeDetection', ...
                        '1A_SpikeDetectedData');
            end 

            channelLayout = Params.channelLayout; % sPerRecording{ExN};
            [spikeMatrix, spikeTimes, Params, Info] = formatSpikeTimes(char(Info.FN), ...
                Params, Info, spikeDetectedDataFolder, channelLayout);

            Params.networkActivityFolder = idvNetworkAnalysisFNFolder;

            coords = Params.coords; %{ExN};
            channels = Params.channels; %{ExN};
            NetMet = ExtractNetMet(adjMs, spikeTimes, ...
                Params.FuncConLagval, Info,HomeDir,Params, spikeMatrix, coords, channels, oneFigureHandle);

            ExpMatFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData',Params.Date), 'ExperimentMatFiles');
            infoFnFname = strcat(char(Info.FN),'_',Params.Date,'.mat');
            infoFnFilePath = fullfile(ExpMatFolder, infoFnFname);

            save(infoFnFilePath, 'Info', 'Params', 'spikeTimes', 'Ephys', 'adjMs','NetMet', '-append')

            clear adjMs

        end

        % save and export network data to spreadsheet
        saveNetMet(ExpName, Params, HomeDir)
 
    
        % Set up one figure handle to save all the figures
        oneFigureHandle = NaN;
        oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);


    % Aggregate all files and run density analysis to determine boundaries
    % for node cartography
    if Params.autoSetCartographyBoundaries
        usePriorNetMet = 0;  % set to 0 by default
        if Params.priorAnalysis==1 && usePriorNetMet
            experimentMatFileFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
            % cd(fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles'));   
            fig_folder = fullfile(Params.priorAnalysisPath, ...
                '4_NetworkActivity', '4B_GroupComparisons', '7_DensityLandscape');
        else
            experimentMatFileFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData', Params.Date), 'ExperimentMatFiles');
            % cd(fullfile(strcat('OutputData', Params.Date), 'ExperimentMatFiles'));  
            fig_folder = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), ...
                '4_NetworkActivity', '4B_GroupComparisons', '7_DensityLandscape');
        end 

        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end 

        ExpList = dir(fullfile(experimentMatFileFolder, '*.mat'));
        add_fig_info = '';

        if Params.autoSetCartographyBoudariesPerLag
            for lag_val = Params.FuncConLagval
                [hubBoundaryWMdDeg, periPartCoef, proHubpartCoef, nonHubconnectorPartCoef, connectorHubPartCoef] = ...
                TrialLandscapeDensity(ExpList, fig_folder, add_fig_info, lag_val, oneFigureHandle);

                if isnan(hubBoundaryWMdDeg)
                    hubBoundaryWMdDeg = Params.hubBoundaryWMdDeg;
                    periPartCoef = Params.periPartCoef;
                    proHubpartCoef = Params.proHubpartCoef;
                    nonHubconnectorPartCoef = Params.nonHubconnectorPartCoef; 
                    connectorHubPartCoef = Params.connectorHubPartCoef;
                end 

                Params.(strcat('hubBoundaryWMdDeg', sprintf('_%.fmsLag', lag_val))) = hubBoundaryWMdDeg;
                Params.(strcat('periPartCoef', sprintf('_%.fmsLag', lag_val))) = periPartCoef;
                Params.(strcat('proHubpartCoef', sprintf('_%.fmsLag', lag_val))) = proHubpartCoef;
                Params.(strcat('nonHubconnectorPartCoef', sprintf('_%.fmsLag', lag_val))) = nonHubconnectorPartCoef;
                Params.(strcat('connectorHubPartCoef', sprintf('_%.fmsLag', lag_val))) = connectorHubPartCoef;
            end 

        else 
            lagValIdx = 1;
            lag_val = Params.FuncConLagval;
            [hubBoundaryWMdDeg, periPartCoef, proHubpartCoef, nonHubconnectorPartCoef, connectorHubPartCoef] = ...
                TrialLandscapeDensity(ExpList, fig_folder, add_fig_info, lag_val(lagValIdx), oneFigureHandle);
            Params.hubBoundaryWMdDeg = hubBoundaryWMdDeg;
            Params.periPartCoef = periPartCoef;
            Params.proHubpartCoef = proHubpartCoef;
            Params.nonHubconnectorPartCoef = nonHubconnectorPartCoef;
            Params.connectorHubPartCoef = connectorHubPartCoef;
        end 

        % save the newly set boundaries to the Params struct
        experimentMatFileFolderToSaveTo = fullfile(Params.outputDataFolder, ...
                strcat('OutputData', Params.Date), 'ExperimentMatFiles');
        for nFile = 1:length(ExpList)
            FN = ExpList(nFile).name;
            FNPath = fullfile(experimentMatFileFolderToSaveTo, FN);
            save(FNPath, 'Params', '-append')
        end 


    end 
        % Run through each file to do node cartography analysis 
        for ExN = 1:length(ExpName) 
            if Params.priorAnalysis==1 && Params.startAnalysisStep==4 && usePriorNetMet
                experimentMatFileFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
                experimentMatFilePath = fullfile(experimentMatFileFolder, strcat(char(ExpName(ExN)),'_',Params.priorAnalysisDate,'.mat'));
                expData = load(experimentMatFilePath, 'spikeTimes','Ephys','adjMs','Info', 'NetMet');
            else
                experimentMatFileFolder = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), 'ExperimentMatFiles');
                experimentMatFilePath = fullfile(experimentMatFileFolder, strcat(char(ExpName(ExN)),'_',Params.Date,'.mat'));
                expData = load(experimentMatFilePath,'Info','Params', 'spikeTimes','Ephys','adjMs', 'NetMet');
            end
            fileNameFolder = fullfile(Params.outputDataFolder, strcat('OutputData',Params.Date), ...
                                      '4_NetworkActivity', '4A_IndividualNetworkAnalysis', ...
                                      char(expData.Info.Grp), char(expData.Info.FN));
            originalCoords = Params.coords; %{ExN};
            originalChannels = Params.channels; %{ExN};
            Params.ExpNameGroupUseCoord = 1;
            NetMet = calNodeCartography(expData.adjMs, Params, expData.NetMet, expData.Info, originalCoords, originalChannels, ...
            HomeDir, fileNameFolder, oneFigureHandle);
            % save NetMet now that we have node cartography data as well
            experimentMatFileFolderToSaveTo = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), 'ExperimentMatFiles');
            experimentMatFilePathToSaveTo = fullfile(experimentMatFileFolderToSaveTo, strcat(char(expData.Info.FN),'_',Params.Date,'.mat'));

            spikeTimes = expData.spikeTimes;
            Ephys = expData.Ephys;
            adjMs = expData.adjMs;
            Info = expData.Info;
            save(experimentMatFilePathToSaveTo,'Info','Params','spikeTimes','Ephys','adjMs','NetMet')
        end
        
    end 
    
    
    %%% 4B: create combined plots %%%
    if strcmp(Params.startAnalysisSubStep, 'ALL') || strcmp(Params.startAnalysisSubStep, 'B')
        
        if strcmp(Params.startAnalysisSubStep, 'ALL')
            experimentMatFileFolder = fullfile(Params.outputDataFolder, ... 
                strcat('OutputData',Params.Date), 'ExperimentMatFiles');
        elseif strcmp(Params.startAnalysisSubStep, 'B')
            experimentMatFileFolder =  fullfile(Params.outputDataFolder, ... 
                strcat('OutputData', Params.priorAnalysisDate), 'ExperimentMatFiles');
        end
        % (everything except node cartography)
        PlotNetMet(ExpName, Params, experimentMatFileFolder, oneFigureHandle)

        % Plot node cartography metrics across all recordings     
        NetMetricsE = Params.networkLevelNetMetToPlot;
        NetMetricsC = Params.unitLevelNetMetToPlot;
        experimentMatFileFolderToSaveTo = fullfile(Params.outputDataFolder, ...
            strcat('OutputData', Params.Date), 'ExperimentMatFiles');
        combinedData = combineExpNetworkData(ExpName, Params, NetMetricsE, ...
            NetMetricsC, experimentMatFileFolder, experimentMatFileFolderToSaveTo);
        figFolder = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), ...
            '4_NetworkActivity', '4B_GroupComparisons', '6_NodeCartographyByLag');
        plotNetMetNodeCartography(combinedData, ExpName,Params, HomeDir, figFolder, oneFigureHandle)
    end 

    
    %%%%% 4A: individual plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4A Part 1: General network metrics
    % Make network plots with shared colorbar and edge weight widths etc.
    if strcmp(Params.startAnalysisSubStep, 'ALL') || strcmp(Params.startAnalysisSubStep, 'A')
        outputDataDateFolder = fullfile(Params.outputDataFolder, ...
            strcat('OutputData', Params.Date));
        
        % Search current analysis folder for csv, if none found, then use
        % the one from the previously analysed data
        spreadsheetFname = strcat('NetworkActivity_RecordingLevel', '.csv');
        if isfile(fullfile(outputDataDateFolder, spreadsheetFname))
            minMax = findMinMaxNetMetTable(outputDataDateFolder, Params);
        else
            minMax = findMinMaxNetMetTable(Params.priorAnalysisPath, Params);
        end
        
        minMax.EW = [0.1, 1];
        Params.metricsMinMax = minMax;
        Params.useMinMaxBoundsForPlots = 1;
        Params.sideBySideBoundPlots = 1;
        for ExN = 1:length(ExpName) 

            % Make figure handle per recording 
            % Set up one figure handle to save all the figures
            oneFigureHandle = NaN;
            oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);

            disp(ExpName(ExN))
            % load NetMet 
            if isfile(fullfile(outputDataDateFolder, spreadsheetFname))
                experimentMatFileFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData', Params.Date), 'ExperimentMatFiles');
                experimentMatFilePath = fullfile(experimentMatFileFolder, ...
                    strcat(char(ExpName(ExN)),'_',Params.Date,'.mat'));
            else
                experimentMatFileFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
                ExpFPathSearchName = dir(fullfile(experimentMatFileFolder, [char(ExpName(ExN)), '*.mat'])).name;
                experimentMatFilePath = fullfile(experimentMatFileFolder, ExpFPathSearchName);
            end

            expData = load(experimentMatFilePath);
            idvNetworkAnalysisGrpFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), '4_NetworkActivity', ...
                '4A_IndividualNetworkAnalysis', char(expData.Info.Grp));

            idvNetworkAnalysisFNFolder = fullfile(idvNetworkAnalysisGrpFolder, char(expData.Info.FN));
            if ~isfolder(idvNetworkAnalysisFNFolder)
                mkdir(idvNetworkAnalysisFNFolder)
            end 

            Params.networkActivityFolder = idvNetworkAnalysisFNFolder;

            originalCoords = Params.coords; %{ExN};
            originalChannels = Params.channels; %{ExN};
            PlotIndvNetMet(expData, Params, expData.Info, originalCoords, originalChannels,  oneFigureHandle)

            if Params.showOneFig
                clf(oneFigureHandle)
            else
                close all 
            end 

            % Testing: Always close all the figures
            close all
            if strcmp(Params.verboseLevel, 'High')
                getMemoryUsage
                numTotalGraphicObjects = length(findall(groot));
                fprintf(sprintf('Total number of graphic objects: %.f \n', numTotalGraphicObjects))
            end 

        end

        % 4A Part 2: NMF individual plot
        if Params.includeNMFcomponents
            % Plot NMF 
            if isfile(fullfile(outputDataDateFolder, spreadsheetFname))
                experimentMatFolder = fullfile(Params.outputDataFolder, ...
                    strcat('OutputData', Params.Date), 'ExperimentMatFiles');
            else
                experimentMatFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
            end
            
            plotSaveFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date), '4_NetworkActivity', ...
                '4A_IndividualNetworkAnalysis');
            plotNMF(experimentMatFolder, plotSaveFolder, Params)
        end 

        % 4A Part 3: Node cartography individual plot
        % Plot node cartography plots using either custom bounds or
        % automatically determined bounds

        % Group the ExpNames by their file identity, to anchor coordinates to
        % the last DIV
        for ExN = 1:length(ExpName)
            stringParts = strsplit(ExpName{ExN}, '_');
            stringPartWithDIV_index = find(contains(stringParts, 'DIV'));
            stringPartWithoutDIV_index = find(~contains(stringParts, 'DIV'));
            ExpNamesWithoutDIV{ExN} = strjoin(stringParts(stringPartWithoutDIV_index), '_');
            divPerExN(ExN) = str2num(stringParts{stringPartWithDIV_index}(4:end));  % first three characters are always 'DIV'
        end
        ExpNameGroup = findgroups(ExpNamesWithoutDIV);
        ExpNameGroupUseCoord = zeros(length(ExpNamesWithoutDIV), 1);
        for groupNum = unique(ExpNameGroup)
            subsetIndex = find(ExpNameGroup == groupNum);
            divMaxIndex = find(divPerExN(subsetIndex) == max(divPerExN(subsetIndex)));
            ExpNameGroupUseCoord(subsetIndex(divMaxIndex)) = 1; 
        end 

        % Loop through each group 
        for groupNum = unique(ExpNameGroup)
            subsetIndex = find(ExpNameGroup == groupNum);
            % Put the recording to be used as the anchor as the first recording
            % to analyse
            [~, sort_index] = sort(ExpNameGroupUseCoord(subsetIndex), 'descend');
            ExNorderToAnalyse = subsetIndex(sort_index);
            for ExN = ExNorderToAnalyse
                Params.ExpNameGroupUseCoord = ExpNameGroupUseCoord(ExN);
                if Params.priorAnalysis==1 && Params.startAnalysisStep==4 && usePriorNetMet
                    experimentMatFileFolder = fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles');
                    experimentMatFilePath = fullfile(experimentMatFileFolder, strcat(char(ExpName(ExN)),'_',Params.priorAnalysisDate,'.mat'));
                    expData = load(experimentMatFilePath, 'spikeTimes','Ephys','adjMs','Info', 'NetMet');
                else
                    experimentMatFileFolder = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), 'ExperimentMatFiles');
                    experimentMatFilePath = fullfile(experimentMatFileFolder, strcat(char(ExpName(ExN)),'_',Params.Date,'.mat'));
                    expData = load(experimentMatFilePath,'Info','Params', 'spikeTimes','Ephys','adjMs', 'NetMet');
                end

                disp(char(expData.Info.FN))

                fileNameFolder = fullfile(Params.outputDataFolder, strcat('OutputData',Params.Date), ...
                                      '4_NetworkActivity', '4A_IndividualNetworkAnalysis', ...
                                      char(expData.Info.Grp), char(expData.Info.FN));

                originalCoords = Params.coords; %{ExN};
                originalChannels = Params.channels; %{ExN};

                % apply the anchoring index used in previous recording
                if Params.ExpNameGroupUseCoord == 0
                    for lagval = Params.FuncConLagval
                         expData.NetMet.(sprintf('AnchoredReorderingIndex%.fmslag', lagval)) = tempNetMet.(sprintf('AnchoredReorderingIndex%.fmslag', lagval)); 
                    end
                end

                % note here the Params from expData is not used because there
                % is a variable I want to preserve between runs... may come
                % up with a better solution down the line...
                % TODO: remove the analysis step in the function, already
                % dealt with earlier
                NetMet = plotNodeCartography(expData.adjMs, Params, expData.NetMet, expData.Info, originalCoords, originalChannels, ...
                      HomeDir, fileNameFolder, oneFigureHandle);
                % plotNodeCartographyProportions(expData.NetMet, Params.FuncConLagval, char(expData.Info.FN), ...
                % Params, fileNameFolder, oneFigureHandle)
                % save NetMet now that we have node cartography data as well
                experimentMatFileFolderToSaveTo = fullfile(Params.outputDataFolder, strcat('OutputData', Params.Date), 'ExperimentMatFiles');
                experimentMatFilePathToSaveTo = fullfile(experimentMatFileFolderToSaveTo, strcat(char(expData.Info.FN),'_',Params.Date,'.mat'));

                spikeTimes = expData.spikeTimes;
                Ephys = expData.Ephys;
                adjMs = expData.adjMs;
                Info = expData.Info;
                save(experimentMatFilePathToSaveTo,'Info','Params','spikeTimes','Ephys','adjMs','NetMet')

                % save the current in use reordering index 
                for lagval = Params.FuncConLagval
                    tempNetMet.(sprintf('AnchoredReorderingIndex%.fmslag', lagval)) = NetMet.(sprintf('AnchoredReorderingIndex%.fmslag', lagval)); 
                end
            end
        end
    end 
    

end

%% Optional step: Run density landscape to determine the boundaries for the node cartography 
if any(strcmp(Params.optionalStepsToRun,'getDensityLandscape')) 
    cd(fullfile(Params.priorAnalysisPath, 'ExperimentMatFiles'));
    
    fig_folder = fullfile(Params.priorAnalysisPath, '4_NetworkActivity', ...
        '4B_GroupComparisons', '7_DensityLandscape');
    if ~isfolder(fig_folder)
        mkdir(fig_folder)
    end 
    
    % loop through multiple DIVs
    for DIV = [14, 17, 21, 24, 28]
        ExpList = dir(sprintf('*DIV%.f*.mat', DIV));
        add_fig_info = strcat('DIV', num2str(DIV));
        [hubBoundaryWMdDeg, periPartCoef, proHubpartCoef, nonHubconnectorPartCoef, connectorHubPartCoef] ...
            = TrialLandscapeDensity(ExpList, fig_folder, add_fig_info, Params.cartographyLagVal);
    end 
end 

%% Optional step: statistics and classification of genotype / ages
% Statistics and machine learning settings 
% Params.classificationTarget = 'AgeDiv';  % which property of the recordings to classify 
Params.classification_models = {'linearSVM', 'kNN', 'decisionTree', 'LDA'}; % 'fforwardNN'
Params.regression_models = {'svmRegressor', 'regressionTree', 'ridgeRegression', 'fforwardNN'};
Params.statsRandomSeed = 1;

if any(strcmp(Params.optionalStepsToRun,'Stats'))
    if Params.showOneFig
        if ~isfield(Params, 'oneFigure')
            Params.oneFigure = figure;
        end 
    end 
    
    if Params.priorAnalysis && Params.startAnalysisStep >= 5
        statsDataFolder = Params.priorAnalysisPath;
    else
        statsDataFolder = fullfile(Params.outputDataFolder, ...
                strcat('OutputData',Params.Date));
    end 
    
    nodeLevelFile = fullfile(statsDataFolder, 'NetworkActivity_NodeLevel.csv');
    nodeLevelData = readtable(nodeLevelFile);
    
    recordingLevelFile = fullfile(statsDataFolder, 'NetworkActivity_RecordingLevel.csv');
    recordingLevelData = readtable(recordingLevelFile);
    
    %{
    if sum(ismember(recordingLevelData.Properties.VariableNames, 'AgeDiv')) == 0 
        if sum(ismember(recordingLevelData.Properties.VariableNames, 'DIV')) == 1
            recordingLevelData.AgeDiv = recordingLevelData.DIV;
        end 
    end 

    if sum(ismember(recordingLevelData.Properties.VariableNames, 'eGrp')) == 0 
        if sum(ismember(recordingLevelData.Properties.VariableNames, 'Grp')) == 1
            recordingLevelData.eGrp = recordingLevelData.Grp;
        end 
    end 
    %}

    % Traditional Statistics
    statsTable = doStats(0, recordingLevelData, Params);
    statsTableSavePath = fullfile(statsDataFolder, 'stats.csv');
    writetable(statsTable, statsTableSavePath);
    oneFigureHandle = NaN;
    oneFigureHandle = checkOneFigureHandle(Params, oneFigureHandle);
    plotSaveFolder = fullfile(statsDataFolder, '5_Stats');
    if ~isfolder(plotSaveFolder)
        mkdir(plotSaveFolder)
    end 
    plotStats(statsTable, plotSaveFolder, Params, oneFigureHandle)
    
    % Classification and Regression 
    for lag_val = Params.FuncConLagval
        plotSaveFolder = fullfile(statsDataFolder, '5_Stats', sprintf('%.fmsLag', lag_val));
        if ~isfolder(plotSaveFolder)
            mkdir(plotSaveFolder)
        end 
        featureCorrelation(nodeLevelData, recordingLevelData, Params, lag_val, plotSaveFolder);
        doLDA(recordingLevelData, Params, lag_val);
        doClassification(recordingLevelData, Params, lag_val, plotSaveFolder, oneFigureHandle);
    end 


end 

%% Optional step : combine plots across DIVs
if any(strcmp(Params.optionalStepsToRun,'combineDIVplots'))
    if Params.priorAnalysis == 1
        featureFolder = fullfile(Params.priorAnalysisPath, '4_NetworkActivity', '4A_IndividualNetworkAnalysis');
    else
        featureFolder = fullfile(Params.outputDataFolder, ['OutputData' Params.Date], '4_NetworkActivity', '4A_IndividualNetworkAnalysis');
    end 
    featureFolderSearch = dir(featureFolder);
    dirFlags = [featureFolderSearch.isdir];
    folderNames = {featureFolderSearch.name};
    groupFolders = folderNames(dirFlags);
    groupFolders = groupFolders(~ismember(groupFolders, {'.', '..'}));
    combinedPlotFolder = fullfile(Params.outputDataFolder, ['OutputData' Params.Date], ...
        '4_NetworkActivity', '4B_GroupComparisons', '8_CombinedPlotsByDiv');
    if 1 - isfolder(combinedPlotFolder)
        mkdir(combinedPlotFolder)
    end 
    
    for grpNameIdx = 1:length(Params.GrpNm)
        combinedPlotGroupFolder = fullfile(combinedPlotFolder, Params.GrpNm{grpNameIdx});
        if 1 - isfolder(combinedPlotGroupFolder)
            mkdir(combinedPlotGroupFolder)
        end 
    end 
    
    Params.includeIdvScaledPlotsInCombinedPlots = 1;
    Params.plotNames = {'3_scaled_MEA_NetworkPlotNodedegreeBetweennesscentrality.png', ...
                        '4_scaled_MEA_NetworkPlotNodedegreeParticipationcoefficient.png', ...
                        '5_scaled_MEA_NetworkPlotNodestrengthLocalefficiency.png', ...
                        '7_scaled_MEA_NetworkPlotNodedegreeAveragecontrollability.png', ... 
                        '8_scaled_MEA_NetworkPlotNodedegreeModalcontrollability.png', ...
                        '2_scaled_MEA_NetworkPlot.png', ...
                        };


    for nGroupFolder = 1:length(groupFolders)

        % get the recording folders 
        groupFolder = fullfile(featureFolder, groupFolders{nGroupFolder});
        groupFolderSearch = dir(groupFolder);
        dirFlags = [groupFolderSearch.isdir];
        folderNames = {groupFolderSearch.name};
        recordingFolders = folderNames(dirFlags);
        recordingFolders = recordingFolders(~ismember(recordingFolders, {'.', '..'}));

        % get the recording name excluding DIV 
        numRecordings = length(recordingFolders);
        recordingNames = cell(numRecordings, 1);
        for recordingIdx = 1:numRecordings
            recordingNameParts = split(recordingFolders{recordingIdx}, '_');
            recordingNames(recordingIdx) = join(recordingNameParts(1:end-1), '_');
        end 

        uniqueRecordings = unique(recordingNames);

        for uniqueRecordingIdx = 1:length(uniqueRecordings)

            recordingName = uniqueRecordings{uniqueRecordingIdx};
            recordingDIVfoldersSearch = dir(fullfile(groupFolder, sprintf('%s*', recordingName)));
            dirFlags = [recordingDIVfoldersSearch.isdir];
            recordingDIVfoldersSearchNames = {recordingDIVfoldersSearch.name};
            recordingDIVfolders = recordingDIVfoldersSearchNames(dirFlags);
            recordingDIVfolders = recordingDIVfolders(~ismember(recordingDIVfolders, {'.', '..'}));

            recordingDIVfolderFullPath = fullfile(groupFolder, recordingDIVfolders{1});
            recordingDIVfoldersSearch = dir(recordingDIVfolderFullPath);
            dirFlags = [recordingDIVfoldersSearch.isdir];
            recordingDIVfoldersSearchNames = {recordingDIVfoldersSearch.name};
            lagFolders = recordingDIVfoldersSearchNames(dirFlags);
            lagFolders = lagFolders(~ismember(lagFolders, {'.', '..'}));

            for lagIdx = 1:length(lagFolders)
                for plotNameIdx = 1:length(Params.plotNames)
                    % make the list of plot paths to combine
                    plotName = Params.plotNames{plotNameIdx};
                    plotPathsToCombine = cell(length(recordingDIVfolders), 1);

                    for divIdx = 1:length(recordingDIVfolders)
                        plotPathsToCombine{divIdx} = fullfile(...
                        groupFolder, recordingDIVfolders{divIdx}, ...
                        lagFolders{lagIdx}, plotName);
                    end 

                    % save the plot in 4B
                    outputFolder = fullfile(combinedPlotFolder, ...
                        groupFolders{nGroupFolder}, recordingName, ...
                        lagFolders{lagIdx});
                    if 1 - isdir(outputFolder)
                        mkdir(outputFolder)
                    end 
                    outputFilePath = fullfile(outputFolder,  plotName(1:end-4)); 
                    %  outputFilePath = fullfile(recordingDIVfolderFullPath, ... 
                    %     sprintf('combined_%s_%s', lagFolders{lagIdx}, plotName(1:end-4)));

                    combinePlots(plotPathsToCombine, outputFilePath, Params)

                end 
            end 

        end

    end
end

%% Optional Step: compare pre-post TTX spike activity 
if any(strcmp(Params.optionalStepsToRun,'comparePrePostTTX')) 
    % see find_best_spike_result.m for explanation of the parameters
    Params.prePostTTX.max_tolerable_spikes_in_TTX_abs = 100; 
    Params.prePostTTX.max_tolerable_spikes_in_grounded_abs = 100;
    Params.prePostTTX.max_tolerable_spikes_in_TTX_per_s = 1; 
    Params.prePostTTX.max_tolerable_spikes_in_grounded_per_s = 1;
    Params.prePostTTX.start_time = 0;
    Params.prePostTTX.default_end_time = 600;  
    Params.prePostTTX.sampling_rate = 1;  
    Params.prePostTTX.threshold_ground_electrode_name = 15;
    Params.prePostTTX.default_grounded_electrode_name = 15;
    Params.prePostTTX.min_spike_per_electrode_to_be_active = 0.5;
    Params.prePostTTX.wavelet_to_search = {'mea', 'bior1p5'};
    Params.prePostTTX.use_TTX_to_tune_L_param = 0;
    Params.prePostTTX.spike_time_unit = 'frame'; 
    Params.prePostTTX.custom_backup_param_to_use = []; 
    Params.prePostTTX.regularisation_param = 10;
    
    
    % Get spike detection result folder
    spike_folder = strcat(HomeDir,'/OutputData',Params.Date,'/1_SpikeDetection/1A_SpikeDetectedData/');
    spike_folder(strfind(spike_folder,'\'))='/';
    
    pre_post_ttx_plot_folder = fullfile(HomeDir, 'OutputData', ... 
        Params.Date, '1_SpikeDetection', '1C_prePostTTXcomparison'); 
    
    find_best_spike_result(spike_folder, pre_post_ttx_plot_folder, Params)
end 



