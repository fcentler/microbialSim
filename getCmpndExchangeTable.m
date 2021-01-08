function [exchangeTable, reportTime] = getCmpndExchangeTable(trajectory, time)
%GETCMPEXCHANGEMATRIX Summary of this function goes here
%   Detailed explanation goes here

%% get index to desired report time
if trajectory.time(end-1) < time %"-1" as for the last index, there is no corresponding flux!
    exchangeTable = 0;
    reportTime = -1;
    return
end
i = 1;
while trajectory.time(i) < time
    i = i + 1;
end
if i > 1 && (time - trajectory.time(i-1)) < (trajectory.time(i) - time)
    i = i - 1;
end
index = i;
reportTime = trajectory.time(index);

% truncate too long model names (Matlab restriction)
for i = 1:length(trajectory.modelNames)
	if length(trajectory.modelNames{i}) > 63
        trajectory.modelNames{i} = trajectory.modelNames{i}(1:63);
        trajectory.modelNames{i}(63) = 'x';
    end
end

%% set up exchange matrix, entries have unit mM/h/L, negative values indicate uptake
exchangeMatrix = zeros(length(trajectory.compoundNames), length(trajectory.modelNames));
for i = 1:length(trajectory.modelNames)
    exchangeMatrix(:, i) = mapExchangeToReactorCompounds(getFluxes(trajectory.FBA(i).fluxes(index, :), trajectory.FBA(i).coupledReactions.ReacID) * trajectory.biomass(index, i) .* transpose(trajectory.FBA(i).coupledReactions.SecretionSense), trajectory.FBA(i).reactorCompoundIDs, length(trajectory.compoundNames));
end

%% prepare exchange table
exchangeTable = array2table(exchangeMatrix, 'VariableNames', trajectory.modelNames, 'RowNames', trajectory.compoundNames);

end