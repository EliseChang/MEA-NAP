function IFR = getInstFR(spikeTimes, bin)

% spikeTimes: cell array output from pipeline after merging
% bin: width of bin in seconds

% IFR: instaneous firing rate

channelsN = length(spikeTimes);

allSpikeTimes = [];
for ch = 1:channelsN
    allSpikeTimes = [allSpikeTimes, spikeTimes{1,ch}.merged];
end

edges = 0:1:600;
spikeCounts = histcounts(allSpikeTimes,edges) / channelsN;
plot(spikeCounts)

t_samp = 0.01; % Sample in bins of 0.01 s
[fr, tt] = hist(allSpikeTimes, [0:t_samp:600]);
[b, a] = butter(1, t_samp); % Prepare the filter
ffr = filtfilt(b, a, fr / t_samp); % Calculate smoothed FR.
figure
plot(tt, ffr);