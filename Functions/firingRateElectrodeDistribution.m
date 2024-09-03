function firingRateElectrodeDistribution(File, Ephys, Params, Info, figFolder, style, oneFigureHandle)
%
% Plots the firing rate distribution across electrodes
% Parameters
% ----------
% File : 
% Ephys : 
% Params : 
% Info : 
% 
% Returns
% -------
% None

%% Create a half violin plot of the firing rate for individual electrodes
FR = Ephys.FR;
FR(Info.grdElecs) = [];
max_ephys_fr = max(FR);
max_ephys_fr = max([FR, 0.1]);  % ensures a minimum of 0.1

p = [50 50 500 600];
set(0, 'DefaultFigurePosition', p)

if Params.showOneFig
    if isgraphics(oneFigureHandle)
        set(oneFigureHandle, 'OuterPosition', p);
    else 
        oneFigureHandle = figure;
        set(oneFigureHandle, 'OuterPosition', p);
    end 
else 
    F1 = figure;
end 

if strcmp(style, 'halfViolin')

    HalfViolinPlot(Ephys.FR,1,[0.5 0.5 0.5], Params.kdeHeight, Params.kdeWidthForOnePoint)
    xlim([0.5 1.5])
    xticks([])
    xlabel(strcat('age',num2str(cell2mat(Info.DIV))))
    aesthetics
    ylim([0 max_ephys_fr+max_ephys_fr*0.15])
    ylabel('mean firing rate per electrode (Hz)')
    figName = '1_FiringRateByElectrodeHalfViolin';

elseif strcmp(style, 'histogram')
    
    h = histogram(Ephys.FR, 'Normalization','percentage');
    xticks(h.BinEdges)
    xlabel('mean firing rate per electrode (Hz)')
    ylabel('percentage')
    aesthetics
    figName = '1_FiringRateByElectrodeHistogram';

end
title({strcat(regexprep(File,'_','','emptymatch')),' '});
ax = gca;
ax.TitleFontSizeMultiplier = 0.7;

%% save the figure
figPath = fullfile(figFolder, figName);

saveas(oneFigureHandle, figPath)
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