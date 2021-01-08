function [ f ] = plotExchangeFluxes(trajectory)
%PLOTRUN Plot the simulated trajectory
%   Detailed explanation goes here
% 
    f = figure;

    % prep time coordinate: fluxes are valid between iteration time points

    for i = 1:(length(trajectory.time) - 1)
        time(i) = (trajectory.time(i) + trajectory.time(i+1)) / 2.0;
    end

    xsize = 2;
    ysize = length(trajectory.FBA);

    panelCounter = 1;
    for i = 1:ysize
        subplot(ysize, xsize, panelCounter)
        panelCounter = panelCounter + 1;

        % plot coupled reactions
        ids = trajectory.FBA(i).coupledReactions.ReacID;

        counter = 1; finalIds = -1; myLegend = {};
        for j = 1:length(ids)
            if sum(abs(trajectory.FBA(i).fluxes(:,ids(j)))) ~= 0.0
                finalIds(counter) = ids(j);
                myLegend(counter) = trajectory.FBA(i).coupledReactions{int2str(finalIds(counter)), 'ReacName'};
                myLegend(counter) = cellstr(strjoin([myLegend(counter) 'Sns' int2str(trajectory.FBA(i).coupledReactions{int2str(finalIds(counter)), 'SecretionSense'})]));
                counter = counter + 1;
            end
        end
        
        if finalIds(1) ~= -1
            plot(time, trajectory.FBA(i).fluxes(:,finalIds))

            xlabel('Time (hours)')
            ylabel('Flux (mmol/gDW/h)')
            title({char(trajectory.modelNames(i)), 'Coupled (non-zero) fluxes'}, 'Interpreter', 'none')
            legend(myLegend, 'Interpreter', 'none'); % no super/subscripting
        end
        subplot(ysize, xsize, panelCounter)
        panelCounter = panelCounter + 1;

        % plot exchange reactions which are not coupled
        ids = setdiff(trajectory.FBA(i).exchangeReactions.ReacID, trajectory.FBA(i).coupledReactions.ReacID);

        counter = 1; finalIds = -1; myLegend = {};
        for j = 1:length(ids)
            if sum(abs(trajectory.FBA(i).fluxes(:,ids(j)))) ~= 0.0
                finalIds(counter) = ids(j);
                myLegend(counter) = trajectory.FBA(i).exchangeReactions{int2str(finalIds(counter)), 'ReacName'};
                myLegend(counter) = cellstr(strjoin([myLegend(counter) 'Sns' int2str(trajectory.FBA(i).exchangeReactions{int2str(finalIds(counter)), 'SecretionSense'})]));
                counter = counter + 1;
            end
        end
        
        if finalIds(1) ~= -1
            plot(time, trajectory.FBA(i).fluxes(:,finalIds))
            xlabel('Time (hours)')
            ylabel('Flux (mmol/gDW/h)')
            title({char(trajectory.modelNames(i)), 'Uncoupled (non-zero) exchange fluxes'}, 'Interpreter', 'none')
            legend(myLegend, 'Interpreter', 'none'); % no super/subscripting
        end
    end
end

