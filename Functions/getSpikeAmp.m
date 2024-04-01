function [minAmp, maxAmp, meanAmp, minSlope, maxSlope] = getSpikeAmp(spikeWaveforms)

% TODO: add documentation
% function originally created to estimate spike waveform properties for
% slope-based online spike detection with MC_Rack

% % Create arrays to store values for each channel
% chMinAmp = zeros(1,60);
% chMinSlope = zeros(1,60);
% chMaxSlope = zeros(1,60);

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
    % chMinAmp(1, ch) = prctile(amps, 10); % set the minimum amplitude as the 10th percentile
    
    preSpikes = inclSpikes(:,1:25); % find the max. value and index before the spike minimum
    [~,spikeMaxIdx] = max(preSpikes, [], 2); % note that only positive peak *before* negative peak is allowed for calculating slope
    slopeDur = (spikeMaxIdx - spikeMinIdx) * 40; % time between spike max. and min. in microseconds
    slopes = amps ./ slopeDur; % slope in microvolts / microseconds
    minSlope = min(slopes(slopes~=0));
    maxSlope = max(slopes);
    % chMinSlope(1, ch) = prctile(slopes, 10);
    % chMaxSlope(1, ch) = prctile(slopes, 90);

else
    minAmp = NaN;
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