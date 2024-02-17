function block = saveSpikeTimesNeo(spikeTimes, endT, outputDir, fileName)

channelsN = length(spikeTimes);
methods = {'bior1p5','bior1p3','db2','thr5'};

% Merge spikes
mergedSpikeTimes = cell(1,channelsN);
for ch = 1:channelsN
    chSpikes = [];
    for m = 1:numel(methods)
        method = methods{m};
        chSpikes = union(chSpikes, spikeTimes{1,ch}.(method), "stable");
        mergedSpikeTimes{1,ch} = chSpikes;
    end
end

% TODO: add single units after sorting

% Block
block = struct();
block.segments = cell(1);
block.name = '10-minute baseline recording';

% Segment
s = 1;
seg = struct();
seg.name = strcat('segment ',num2str(s));

% ChannelIndex
seg.channelindexes = { };
for ch = 1:channelsN
    chx = struct();
    chx.name = strcat('channel ', num2str(ch));
    chx.index = ch;
    
    % Unit
    chx.units = { };
    u = 1; % For un-spike sorted data
    unt = struct();
    unt.name = strcat('unit ', num2str(u));
    unt.spiketrains = { };
    
    % SpikeTrain
    spikes = mergedSpikeTimes{1, ch}; % For un-spike sorted data
    spktr = struct();
    spktr.times = spikes;
    spktr.times_units = 'sec'; % change back to ms after spike sorting
    spktr.t_start = 0;
    spktr.t_start_units = 'sec';
    spktr.t_stop = endT;
    spktr.t_stop_units = 'sec';

    % Add spiketrain to unit
    unt.spiketrains{1} = spktr;

    % Add unit to channelindex
    chx.units{u} = unt;

    % Add channelindex to segment
    seg.channelindexes{ch} = chx;

end
spktr = struct();
spktr.times = spikeTimes;
spktr.times_units = 'ms';
spktr.t_start = 0;
spktr.t_start_units = 'ms';
spktr.t_stop = endT;
spktr.t_stop_units = 'ms';

% Add segment to block
block.segments{1} = seg;

save(fullfile(outputDir,strcat(fileName,'.mat')), 'block', '-v7')

end