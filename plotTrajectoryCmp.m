function [ f ] = plotTrajectoryCmp( varargin )
%PLOTRUN Plot two trajectories together
%   Expected arguments: trajectory1 trajectory2 label1 label2...

    f = figure;

    runsToCmp = nargin / 2;

    numberOfPlots = 1 + length(varargin{1,1}.compoundNames);

    xsize = 3;
    ysize = ceil(numberOfPlots/3.0);

    plotSymbols = {'o', '-', 'x', '-.'};

    providedLegend = {varargin{1,(runsToCmp + 1):nargin}};
    myLegend = {};

    subplot(ysize, xsize, 1)

    % Biomass plots

    for j = 1:runsToCmp
        for k = 1:size(varargin{1,j}.biomass, 2)
            plot(varargin{1,j}.time, varargin{1,j}.biomass(:, k) - varargin{1,j}.biomass(1, k), plotSymbols{j})
            myLegend{length(myLegend) + 1} = strcat(providedLegend{j}, {' '}, varargin{1,j}.modelNames{k});
            hold on
        end
    end

    xlabel('Time (hours)')
    ylabel('Concentration offset (gDW/L)')
    title('Biomass')
    legend([myLegend{:}], 'Interpreter', 'none')
    hold off

    for i = 1:(numberOfPlots - 1)
        subplot(ysize, xsize, i + 1)
        for j = 1:runsToCmp
            plot(varargin{1,j}.time, varargin{1,j}.compounds(:, i), plotSymbols{j})
            hold on
        end
        xlabel('Time (hours)')
        ylabel('Concentration (mmol/L)')
        title(varargin{1,1}.compoundNames(i), 'Interpreter', 'none')
        hold off
    end

end

