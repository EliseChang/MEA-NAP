function spikeMatrix = spikeTimesToSpikeMatrix(spikeTimes, fs, duration)

end_time = duration;
bin_edges = 0:1/fs:end_time;

numUnits = length(spikeTimes);
num_bins = length(bin_edges) - 1;
spikeMatrix = zeros(num_bins, numUnits);

for unit = 1:numUnits
    channel_spike_times = spikeTimes{unit};
    spike_vector = histcounts(channel_spike_times, bin_edges);
    spikeMatrix(:, unit) = spike_vector;
end 


end 