function t = plotSpikeWaveforms(channel, incSpikes, spikeTimes, spikeWaveforms, origTrace, fs, methods, figPos)

%     [~, unique_idx, ~] = mergeSpikes(spike_times{channel},'all'); yt
%     unique spikes for each method -- for now, just plot all spikes

figure('Position',figPos)
t = tiledlayout(2, ceil(length(methods)/2), 'TileSpacing','compact');
title(t, ['Electrode ', num2str(channel)])
for i = 1:length(methods)

    method = methods{i};
    mask = incSpikes{channel}.(method);
    allWaveforms = spikeWaveforms{channel}.(method);
    plotWaveforms = allWaveforms(mask,:); % plot only waveforms for spikes specified in incSpikes
    
    if ~isempty(plotWaveforms)
        plotOrigTrace = zeros(size(plotWaveforms,1),51); % check this??
        plotSpikeTimes = spikeTimes{channel}.(method);
        plotSpikeFrames = int64(plotSpikeTimes * fs);
        for f = 1:length(plotSpikeFrames)
            centreFrame = plotSpikeFrames(f);
            if centreFrame < 25
                plotOrigTrace(f,:) = origTrace(1:51, channel);
            elseif centreFrame > length(origTrace)-25
                plotOrigTrace(f,:) = origTrace(end-50:end, channel);
            else
                plotOrigTrace(f,:) = origTrace(centreFrame-25:centreFrame+25, channel);
            end
        end
    end
%     if ~strcmp(method, 'all')
%         spk_method = find(unique_idx == i);
%         spk_waves_method = spikeWaveforms{channel}.(method);
%     end
%         % convert spike waveform to the appropriate units as well 
%         if isstring(Params.potentialDifferenceUnit)
%             if strcmp(Params.potentialDifferenceUnit, 'V')
%                 spk_waves_method = spk_waves_method .* 10^6;
%             elseif strcmp(Params.potentialDifferenceUnit, 'mV')
%                 spk_waves_method = spk_waves_method .* 10^3;
%             elseif strcmp(Params.potentialDifferenceUnit, 'uV')
%                 spk_waves_method = spk_waves_method;
%             end
%         else 
%             % convert to V by provided multiplication factor, then convert to uV
%             % for plotting
%             spk_waves_method = spk_waves_method .* Params.potentialDifferenceUnit .* 10^6;
%         end 

%     if size(spk_waves_method,2) > 1000
%         spk_waves_method = spk_waves_method(:, round(linspace(1,length(spk_waves_method),1000)));
%     end

    nexttile
    if exist('plotOrigTrace','var')
        plot(plotOrigTrace', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    end
    hold on
    plot(plotWaveforms', 'linewidth', 1.0, 'color', [0 0.4470 0.7410])
    plot(mean(plotWaveforms), 'linewidth', 1.5, 'color', [0 0 0])
    title(method)
    box off;
    axis tight
    pbaspect([1,2,1]);
%         ylim([-6*std(trace) 5*std(trace)]);
    ylabel('Voltage [\muV]')
    set(gca, 'xcolor', 'none');
    aesthetics

    clear method mask allWaveforms plotWaveforms plotOrigTrace

end
