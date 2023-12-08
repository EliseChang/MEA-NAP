function t = plotTrialHeatmap(figPos,stimProt,psthData,psthBin,fs)

figure('Position', figPos)
t = tiledlayout(1, 2, 'TileSpacing','Compact');
maskA = logical(stimProt);
maskB = logical(~stimProt);

% Pattern A trials
trialsA = psthData(maskA,:);
[~,trialOrder] = sort(sum(trialsA,2),'descend');
trialsSorted = trialsA(trialOrder,:);
h(1) = nexttile(t);
imagesc(trialsSorted)
xticklabels(xticks*psthBin/fs*1e3)
aesthetics
clear trialsSorted trialOrder

% Pattern B trials
trialsB = psthData(maskB,:);
[~,trialOrder] = sort(sum(trialsB,2),'descend');
trialsSorted = trialsB(trialOrder,:);
h(2) = nexttile(t);
imagesc(trialsSorted)
xticklabels(xticks*psthBin/fs*1e3)
aesthetics
clear trialsSorted trialOrder

linkaxes([h(1) h(2)])
xlabel(t,"Time post-stimulus (ms)")
ylabel(t, "Trials")
cmin = min([trialsA,trialsB],[],'all');
cmax = max([trialsA,trialsB],[],'all');
set(h, 'Colormap', parula, 'CLim', [cmin cmax])
cb = colorbar();
cb.TickDirection = 'out';
cb.Box = 'off';
cb.Layout.Tile = 'east';
ylabel(cb, 'Spike count per 5 ms')

end