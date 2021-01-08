function [filteredExchangeTable] = filterCmpndExchangeTable(exchangeTable, mode, threshold)
%FILTERCMPNDEXCHANGEMATRIX Filter table of metabolite exchange according to
%specific criteria
%   mode == 0: remove compounds with zero entries only
%   mode == 1: true exchange
%   mode == 2: only uptake/competition
%   mode == 3: only production
% threshold indicates a value, below which fluxes are neglected

rowsToKeep = [];

tmp = table2array(exchangeTable);
tmp(abs(tmp) < threshold) = 0;
tmp2 = array2table(tmp, 'VariableNames', exchangeTable.Properties.VariableNames, 'RowNames', exchangeTable.Properties.RowNames);
exchangeTable = tmp2;

%% iterate across table rows (== compounds)
for i = 1:size(exchangeTable, 1)
    keepIt = false;
    switch mode
        case 0
            if sum(abs(exchangeTable{i,:})) ~= 0.0
              keepIt = true;  
            end
        case 1
            if min(exchangeTable{i,:}) < 0 && max(exchangeTable{i,:}) > 0
                keepIt = true;
            end
        case 2
            if min(exchangeTable{i,:}) < 0 && max(exchangeTable{i,:}) <= 0
                keepIt = true;
            end
        case 3
            if min(exchangeTable{i,:}) >= 0 && max(exchangeTable{i,:}) > 0
                keepIt = true;
            end
        otherwise
            error('Unknown mode.');
    end
    if keepIt
        rowsToKeep(end+1) = i;
    end
end


%% set up table
filteredExchangeTable = exchangeTable(rowsToKeep, :);

%% don't include species with zero values only

for i = size(filteredExchangeTable, 2):-1:1
    if sum(abs(filteredExchangeTable{:,i})) == 0
        filteredExchangeTable(:,i) = [];
    end
end

%writetable(filteredExchangeTable, 'filename.txt', 'WriteRowNames', true, 'Delimiter', '\t')
