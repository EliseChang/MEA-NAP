function plotNetMetNodeCartography(combinedData, ExpName, Params,HomeDir, figFolder)
% Plots network node cartography of individual recordings
% Parameters
% ----------
% combinedData : struct
%     structure with fields of the different groups (eg. genotypes:
%     combinedData.WT, combinedData.KO, ...)
% ExpName : str
% Params : structure
% HomeDir : str 
% 
% Returns
% -------
% 
%

%% groups and DIV

Grps = Params.GrpNm;
AgeDiv = Params.DivNm;

if strcmp(char(Grps{1}),'HET') && strcmp(char(Grps{2}),'KO') && strcmp(char(Grps{3}),'WT')
   clear Grps
   Grps{1} = 'WT'; Grps{2} = 'HET'; Grps{3} = 'KO';
end

%% Plotting 

c1 = [0.8 0.902 0.310]; % light green
c2 = [0.580 0.706 0.278]; % medium green
c3 = [0.369 0.435 0.122]; % dark green
c4 = [0.2 0.729 0.949]; % light blue
c5 = [0.078 0.424 0.835]; % medium blue
c6 = [0.016 0.235 0.498]; % dark blue

eMet = {'NCpn1','NCpn2','NCpn3','NCpn4','NCpn5','NCpn6'}; 

p = [100 100 1200 800]; % this can be ammended accordingly
set(0, 'DefaultFigurePosition', p)
if isfield(Params, 'oneFigure')
    set(Params.oneFigure, 'Position', p);
end 

for l = 1:length(Params.FuncConLagval)
    if ~isfield(Params, 'oneFigure')
        F1 = figure;
    end 
    xt = 1:length(AgeDiv);
    for g = 1:length(Grps)
        h(g) = subplot(length(Grps),1,g);
        eGrp = cell2mat(Grps(g));
        for n = 1:length(eMet)
            eMeti = char(eMet(n));
            for d = 1:length(AgeDiv)
                eDiv = strcat('TP',num2str(d));
                %VNe = strcat(eGrp,'.',eDiv,'.',eMeti);
                %eval(['DatTemp = ' VNe ';']);
                DatTemp = combinedData.(eGrp).(eDiv).(eMeti); 
                % Tim: temp fix to make zero vector of DIV is empty 
                if isempty(DatTemp)
                    DatTemp = zeros(1, length(Params.FuncConLagval));
                end 
                meanTP(d)= nanmean(DatTemp(:,l));
                stdTP(d)= nanstd(DatTemp(:,l));
                xtlabtext{d} = num2str(AgeDiv(d));
            end
            eval(['c = c' num2str(n) ';']);
            UpperStd = meanTP+stdTP; % upper std line
            LowerStd = meanTP-stdTP; % lower std line
            Xf =[xt,fliplr(xt)]; % create continuous x value array for plotting
            Yf =[UpperStd,fliplr(LowerStd)]; % create y values for out and then back
            h1 = fill(Xf,Yf,c,'edgecolor','none');
            % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.
            set(h1,'facealpha',0.3)
            hold on
            eval(['y' num2str(n) '= plot(xt,meanTP,''Color'',c,''LineWidth'',3);']);   % TODO : fix this to not use eval
            xticks(xt)
            xticklabels(xtlabtext)
            xlabel('Age')
            ylabel('node cartography')
            clear DatTemp meanTP meanSTD ValStd UpperStd LowerStd
        end
        linkaxes(h,'xy')
        h(1).XLim = [min(xt)-0.5 max(xt)+0.5];
        title(eGrp)
        aesthetics
        set(gca,'TickDir','out');
        legend([y1 y2, y3, y4, y5, y6],'proportion peripheral nodes','proportion non-hub connectors', ... 
            'proportion non-hub kinless nodes','proportion provincial hubs','proportion connector hubs','proportion kinless hubs','Location','eastoutside')
        legend Box off
    end

    % Export figure
    figName = strcat(['NodeCartography', num2str(Params.FuncConLagval(l)), 'mslag']);
    figPath = fullfile(figFolder, figName);
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG)

    % Close figure or clear the one shared figures
    if ~isfield(Params, 'oneFigure')
        close(gcf)
    else
        set(0, 'CurrentFigure', Params.oneFigure);
        clf reset
    end 
end




end 