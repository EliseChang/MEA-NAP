function electrodeHeatMaps(Info, Ephys, spikeFreqMax, Params, coords, figFolder, oneFigureHandle)
% Plots the firing rate of each node / electrode with a circle representing 
% the spatial location of the electrode / node, and the color representing 
% the firing rate (spikes/s)
% 
% Parameters
% -----------
% FN : char
%     name of the recording file (including extention)
%     eg. 'HP_tc165_DIV28_24Feb2022.mat'
% spikeMatrix : T x N matrix
% channels : 1 x N vector 
%     the integer ID of each channel 
% spikeFreqMax : float
%     the maximum mean firing rate to include in the heatmap that is scaled 
%     the same way across the entire batch of recordings, to allow for visual
%     comparison between recordings
% Params : struct
%    the following fields are used 
%    coords : N x 2 matrix 
%         the first column is the x-coordinate of each node 
%         the second column is the y-coordinate of each node
%    fs : int 
%         the sampling frequency of the recording
%    figMat : bool 
%         whether to save figure in .fig format
%    figPng : bool 
%         whether to save figure in .png format
%    figEps : bool 
%         whether to save figure in .eps format
% 
% Returns 
% -------
% None 

%% Get file and ephys variables

FN = char(Info.FN);
channels = Info.channels;
allFR = Ephys.FR;
groundElecs = Info.grdElecs;
inactiveElecs = setdiff(Ephys.inactiveElecs, groundElecs);

% Plot settings for inactive and ground electrodes
groundElecColor = 'black';
inactiveElecColor = 'red';
%% plot
p = [50 100 1150 570];

if Params.showOneFig
    if isgraphics(oneFigureHandle)
        set(oneFigureHandle, 'OuterPosition', p);
    else 
        oneFigureHandle = figure;
        set(oneFigureHandle, 'OuterPosition', p);
    end 
else 
    F1 = figure;
    F1.OuterPosition = p;
end 

tiledlayout(1,2)
aesthetics; axis off; hold on

%% coordinates

% Perform transpose if not column vector 
if size(channels, 1) == 1
    channels = channels'; 
end 

% xc = Params.coords(:,1);
% yc = Params.coords(:,2);

xc = coords(:, 1);
yc = coords(:, 2); 

%% plot electrodes
% TODO: I think a lot of rectangle coloring can be simplified
mycolours = colormap;
numCbarTicks = 5;


%% Left electrode plot (scaled to individual recording)
nexttile
uniqueXc = sort(unique(xc));
nodeScaleF = 2/3; 

% makeMinSpikeCountZero = 1;
% 
% if makeMinSpikeCountZero == 1
    minSpikeCountToPlot = 0;
% else 
%     minSpikeCountToPlot = min(allFR);
% end 

for i = 1:length(allFR)

    pos = [xc(i)-(0.5*nodeScaleF) yc(i)-(0.5*nodeScaleF) nodeScaleF nodeScaleF];
            
            if ismember(i, groundElecs) % (allFR(i) - minSpikeCountToPlot) / (prctile(allFR,95,'all') - minSpikeCountToPlot) <= 0
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
                    'w','EdgeColor',groundElecColor)
%                 rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
%                     mycolours(ceil(length(mycolours)*((allFR(i)- minSpikeCountToPlot)/(prctile(allFR,99,'all')-minSpikeCountToPlot))+0.00001),1:3),'EdgeColor','w','LineWidth',0.1)
            elseif ismember(i, inactiveElecs) % (allFR(i) - minSpikeCountToPlot) / (prctile(allFR,95,'all') - minSpikeCountToPlot) <= 0
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
                    'w','EdgeColor',inactiveElecColor)
            else
                try
                    colorToUse = mycolours(ceil(length(mycolours) * ((allFR(i) - minSpikeCountToPlot)/(prctile(allFR,99,'all')-minSpikeCountToPlot))),1:3);
                    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',colorToUse,'EdgeColor','w','LineWidth',0.1)
                catch
                    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',mycolours(length(mycolours),1:3),'EdgeColor','w','LineWidth',0.1)
                end
            end
        
    if Params.includeChannelNumberInPlots 
        text(pos(1) + 0.5 * nodeScaleF, pos(2) + 0.5 * nodeScaleF, ...
            sprintf('%.f', channels(i)), 'HorizontalAlignment','center')
    end 
end
ylim([min(yc) - 1, max(yc) + 1])
xlim([min(xc) - 1, max(xc) + 1])
axis off

cb = colorbar;
cb.Box = 'off';
cb.Ticks = linspace(0, 1, numCbarTicks);

tickLabels = cell(numCbarTicks, 1);
for nTick = 1:numCbarTicks
    if nTick == 1
        tickLabels{nTick} = num2str(minSpikeCountToPlot);
    else 
        tickLabels{nTick} = num2str(round((nTick-1) / numCbarTicks * prctile(allFR,99,'all'),2));
    end 
end 

cb.TickLabels = tickLabels;


cb.TickDirection = 'out';
cb.Label.String = 'mean firing rate (Hz)';
title({strcat(regexprep(FN,'_','','emptymatch'),' Electrode heatmap scaled to recording'),' '});

%% Right electrode plot (scaled to all recordings)

nexttile
for i = 1:length(allFR)
    pos = [xc(i)-(0.5*nodeScaleF) yc(i)-(0.5*nodeScaleF) nodeScaleF nodeScaleF];
    if ismember(i, groundElecs) % (allFR(i) - minSpikeCountToPlot) / (prctile(allFR,95,'all') - minSpikeCountToPlot) <= 0
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
            'w','EdgeColor',groundElecColor);
%                 rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
%                     mycolours(ceil(length(mycolours)*((allFR(i)- minSpikeCountToPlot)/(prctile(allFR,99,'all')-minSpikeCountToPlot))+0.00001),1:3),'EdgeColor','w','LineWidth',0.1)
    elseif ismember(i, inactiveElecs) % (allFR(i) - minSpikeCountToPlot) / (prctile(allFR,95,'all') - minSpikeCountToPlot) <= 0
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
            'w','EdgeColor',inactiveElecColor);
    else
        
        try
            colorToUse = mycolours(ceil(length(mycolours) * ((allFR(i) - minSpikeCountToPlot)/(prctile(allFR,99,'all')-minSpikeCountToPlot))),1:3);
            rectangle('Position',pos,'Curvature',[1 1],'FaceColor',colorToUse,'EdgeColor','w','LineWidth',0.1)
        catch
            rectangle('Position',pos,'Curvature',[1 1],'FaceColor',mycolours(length(mycolours),1:3),'EdgeColor','w','LineWidth',0.1)
        end
    end
end
ylim([min(yc)-1 max(yc)+1])
xlim([min(xc)-1 max(xc)+1])
axis off

cb = colorbar;
cb.Box = 'off';
cb.Ticks = linspace(0, 1, numCbarTicks);

tickLabels = cell(numCbarTicks, 1);
for nTick = 1:numCbarTicks
    if nTick == 1
        tickLabels{nTick} = num2str(minSpikeCountToPlot);
    else 
        tickLabels{nTick} = num2str(round((nTick-1) / numCbarTicks * spikeFreqMax, 2));
    end 
end 

cb.TickLabels = tickLabels;
cb.TickDirection = 'out';
cb.Label.String = 'mean firing rate (Hz)';
title({strcat(regexprep(FN,'_','','emptymatch'),' Electrode heatmap scaled to entire data batch'),' '});

% dummy plot for legend
hold on
plot(NaN, NaN, 'Color',groundElecColor)
if ~isempty(inactiveElecs)
    plot(NaN, NaN, 'Color',inactiveElecColor)
    legend({'Ground', 'Inactive'})
else
    legend({'Ground'})
end
legend('Location','best', 'Orientation','horizontal')
legend("boxoff")
hold off

% save figure
figName = 'Heatmap';
figPath = fullfile(figFolder, figName);

if Params.showOneFig
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle);
else
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG)
end 

if Params.showOneFig
    clf(oneFigureHandle, 'reset')
else 
    close(F1);
end 

end