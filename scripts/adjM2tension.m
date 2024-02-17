% Directories
srcDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData11Nov2022\ExperimentMatFiles'; addpath(srcDir)
srcDate = '11Nov2022';
outputDir = 'C:\Users\elise\Python\ReservoirComputing\data\connectivity'; addpath(outputDir)

% Get file names
metadataSpreadsheet = 'MEC.xlsx'; % file name
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
addpath(outputDir)
xlSheet = 'Sheet1';
xlRange = 'A2:A129';
[~,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1); % name of sample

lags = [10]; %  25, 50

for l = lags

    for n = 1:length(samples)
        load(fullfile(srcDir, strcat(samples{n}, "_", srcDate)), "-mat", "adjMs")
        adjM = adjMs.(strcat("adjM", num2str(l), "mslag"));
%         mat2np(adjM, fullfile(outputDir, strcat(ExpName{ExN}, "_", num2str(l), "mslag")),...
%             'float64')
        % TODO: add inhibitory nodes
        writematrix(adjM, fullfile(outputDir,strcat(samples{n}, '.csv')))
    end

end
