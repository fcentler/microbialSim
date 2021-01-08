function [f] = plotTopBiomassCompounds(trajectory, maxNumberOfBiomass, maxNumberOfCompounds)
%PLOTTOPBIOMASSCOMPOUNDS Summary of this function goes here
%   Detailed explanation goes here
f = figure

subplot(1,2,1)
%xlim([0 1])
numberOfData = min(maxNumberOfBiomass, size(trajectory.biomass, 2));
colors = colorcube(numberOfData);
[dummy, idx] = sort(trajectory.biomass(end,:), 'descend');
for i=1:numberOfData
            hold all
            h(i)=plot(trajectory.time,trajectory.biomass(:,idx(i)),'Color',colors(i,:));
end
l = legend(trajectory.modelNames(idx(1:numberOfData)), 'Interpreter', 'none', 'Location', 'southoutside');
title('Biomass')
xlabel('Time (hours)')
ylabel('Biomass concentration (gDW/L)')

subplot(1,2,2)
%xlim([0 1])
numberOfData = min(maxNumberOfCompounds, size(trajectory.compounds, 2));
colors = colorcube(numberOfData);
[dummy, idx] = sort(trajectory.compounds(end,:), 'descend');
for i=1:numberOfData
            hold all
            h(i)=plot(trajectory.time,trajectory.compounds(:,idx(i)),'Color',colors(i,:));
end
l = legend(trajectory.compoundNames(idx(1:numberOfData)), 'Interpreter', 'none', 'Location', 'southoutside');
title('Compounds')
xlabel('Time (hours)')
ylabel('Concentration (mmol/L)')

end

