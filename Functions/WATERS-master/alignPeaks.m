function [newSpikeTimes, spikeWaveforms, spikeFreeTrace] = alignPeaks(spikeTimes, unit, fs, trace, win,...
    artifactFlg, varargin)

% Description:
%   Aligns spikes by negative peaks and removes artifacts by amplitude

% INPUT:
%   spikeTimes: vector containing spike times
%   unit: unit in which spike times are given ('ms' | 's' | 'frames')
%   fs: sampling frequency (Hz)
%   trace: [n x 1] filtered voltage trace
%   win: [scalar] window around the spike in [frames]; 
%        this is the width of the bin that is used to search for the peak 
%        c.f. with waveform_width, which is the half width of the aligned
%        spike (so the full width will be [-waveform_width peak
%        +waveform_wdith])
%   artifactFlg: [logical] flag for artifact removal; 1 to remove artifacts, 0 otherwise
%
% Optional arguments (only used in post-hoc artifact removal)
%   varargin{1} = minPeakThrMultiplier;
%   varargin{2} = maxPeakThrMultiplier;
%   varargin{3} = posPeakThrMultiplier;

% OUTPUT:
%   spikeTimes: [#spikes x 1] new spike times aligned to the negative amplitude peak
%   spikeWaveforms: [51 x #spikes] waveforms of the detected spikes
%   spikeFreeTrace: [#frames-#spikes*win x 1] filtered voltage trace with spike cut-outs removed

% Author:
%   Jeremy Chabros, University of Cambridge, 2020
%   email: jjc80@cam.ac.uk
%   github.com/jeremi-chabros

if strcmp(unit, 'ms')
    spikeFrames = round(spikeTimes*(fs/1000));
elseif strcmp(unit, 's')
    spikeFrames = round(spikeTimes*fs);
elseif strcmp(unit,'frames')
    spikeFrames = spikeTimes;
end

waveform_width = 25;

% Obtain thresholds for artifact removal
threshold = median(trace) - median(abs(trace))/0.6745;

% TODO: should be a user option to use multiplier or absolute threshold
% Comment out to use the multiplier
% if artifactFlg
%     minPeakThr = threshold * varargin{1};
%     maxPeakThr = -threshold * varargin{2};
%     posPeakThr = -threshold * varargin{3};
% end

% Uses absolute threshold in microvolts
if artifactFlg
    minPeakThr = varargin{1}; % e.g. -7 uV
    maxPeakThr = varargin{2}; % e.g. -100 uV
    posPeakThr = varargin{3}; % % e.g. 100 uV
end

sFr = zeros(length(spikeFrames),1);
spikeWaveforms = zeros(length(spikeFrames),waveform_width*2+1);
spikeFreeTrace = trace;
for i = 1:length(spikeFrames)
    
    if spikeFrames(i)+win < length(trace)-1 && spikeFrames(i)-win > 1
        
        % Look into a window around the spike
        bin = trace(spikeFrames(i)-win:spikeFrames(i)+win);
        spikeFreeTrace(spikeFrames(i)-win:spikeFrames(i)+win) = NaN;
        negativePeak = min(bin);
        positivePeak = max(bin);
        pos = find(bin == negativePeak);
        
        % Remove artifacts and assign new timestamps
        if artifactFlg
            if (negativePeak < minPeakThr) && (positivePeak < posPeakThr) && (negativePeak > maxPeakThr)
                newSpikeFrame = spikeFrames(i)+pos-win;
                if newSpikeFrame+waveform_width < length(trace) && newSpikeFrame-waveform_width > 1
                    waveform = trace(newSpikeFrame-waveform_width:newSpikeFrame+waveform_width);
                    sFr(i) = newSpikeFrame;
                    spikeWaveforms(i, :) = waveform;
                end
            end
        else
            newSpikeFrame = spikeFrames(i)+pos-win;
            if newSpikeFrame+waveform_width < length(trace) && newSpikeFrame-waveform_width > 1
                waveform = trace(newSpikeFrame-waveform_width:newSpikeFrame+waveform_width);
                sFr(i) = newSpikeFrame;
                spikeWaveforms(i, :) = waveform;
            end
        end
    end
end

% Pre-allocation & logical indexing made it a lot faster
% than using (end+1) indexing in the loop above
% Convert back to original units
if strcmp(unit, 'ms')
    newSpikeTimes = sFr/(fs/1000);
elseif strcmp(unit, 's')
    newSpikeTimes = sFr/fs;
elseif strcmp(unit, 'frames')
    newSpikeTimes = sFr;
end
spikeFreeTrace = rmmissing(spikeFreeTrace); % remove NaN values
end

