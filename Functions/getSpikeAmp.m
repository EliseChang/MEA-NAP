function [minAmp, maxAmp, meanAmp, minSlope, maxSlope] = getSpikeAmp(spikeWaveforms)

% Description:
%   Computes spike waveform properties for (useful for slope-based online
%   spike detection with MC_Rack)
% Inputs (required):
%   spikeWaveforms: spikes x frames matrix of voltage data for spikes detected on
%   a single electrode
%   
% Outputs:
    % minAmp: minimum peak-trough/trough-peak voltage amplitude across all spikes
    % maxAmp: maximum peak-trough/trough-peak voltage amplitude across all spikes
    % meanAmp: peak-trough voltage amplitude averaged across all spikes
    % minSlope: minimum peak-trough (downward deflection only) i.e. shallowest slope across all spikes
    % maxSlope: maximum peak-trough (downward deflection only) i.e. steepest slope across all spikes

% Note: commented lines from previous version which took the whole
% spikeWaveforms structure output from batch spike detection and output
% min. and max. values across all electrodes

% % Create arrays to store values for each channel
% chMinAmp = zeros(1,60);
% chMinSlope = zeros(1,60);
% chMaxSlope = zeros(1,60);
% for ch = 1:60
%     spikes = spikeWaveforms{1, ch}.(method);

[~, allSpikesMinIdx] = min(spikeWaveforms, [], 2); % find position of negative peaks
inclSpikes = spikeWaveforms(allSpikesMinIdx == 25, :); % the 'clean' spike waveforms have negative peak centred on the midpoint of the cutout window, 
                                                %  so just consider these for now
if ~isempty(inclSpikes)

    [spikeMinAmp,spikeMinIdx] = min(inclSpikes, [], 2);
    [spikeMaxAmp,~] = max(inclSpikes, [], 2); % note that positive peak is allowed before or after negative peak for calculating amplitude

    amps = spikeMaxAmp - spikeMinAmp;
    meanAmp = mean(amps);
    minAmp = min(amps(amps~=0));
    maxAmp = max(amps);
    % chMinAmp(1, ch) = minAmp
    % chMaxAmp(1, ch) = maxAmp
    % clear minAmp maxAmp
    
    preSpikes = inclSpikes(:,1:25); % find the max. value and index before the spike minimum
    [~,spikeMaxIdx] = max(preSpikes, [], 2); % note that only positive peak *before* negative peak is allowed for calculating slope
    slopeDur = (spikeMaxIdx - spikeMinIdx) * 40; % time between spike max. and min. in microseconds
    slopes = amps ./ slopeDur; % slope in microvolts / microseconds
    minSlope = min(slopes(slopes~=0));
    maxSlope = max(slopes);
    % chMinSlope(1, ch) = minSlope
    % chMaxSlope(1, ch) = maxSlope
    % clear minSlope maxSlope

else
    minAmp = NaN;
    maxAmp = NaN;
    meanAmp = NaN;
    minSlope = NaN;
    maxSlope = NaN;
end

clear allSpikesMinIdx inclSpikes spikeMinAmp spikeMinIdx preSpikes spikeMaxAmp spikeMaxIdx amps slopeDur slopes

% minAmp = min(chMinAmp(chMinAmp ~= 0));
% minSlope = min(chMinSlope(chMinSlope ~= 0));
% maxSlope = max(chMaxSlope(chMaxSlope ~= 0));
% clear chMinAmp chMinSlope chMaxSlope

end