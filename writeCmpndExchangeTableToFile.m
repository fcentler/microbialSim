function [exchangeMatrix] = writeCmpndExchangeTableToFile(exchangeTable, filename)
%WRITECMPNDEXCHANGETABLETOFILE Summary of this function goes here
%   Detailed explanation goes here

dimension = sum(size(exchangeTable));

exchangeMatrix = zeros(dimension, dimension);

for i = 1:size(exchangeTable, 1)
    for j = 1:size(exchangeTable, 2)
        if exchangeTable{i, j} < 0
            exchangeMatrix(i, j) = -exchangeTable{i,j};
        else
            exchangeMatrix(j + size(exchangeTable, 1), i + size(exchangeTable, 2)) = exchangeTable{i,j};
        end
    end
end

stringsToReplace = ['[', ']'];
replacementStrings = ['_', '_'];
cmpNames = exchangeTable.Properties.RowNames;
spcsNames = exchangeTable.Properties.VariableNames;
% leading numbers not allowed as variable name in table
cmpNames = regexprep(cmpNames, '^([0-9])', 'C_$1');
for i = 1:length(stringsToReplace)
    cmpNames = strrep(cmpNames, stringsToReplace(i), replacementStrings(i));
    spcsNames = strrep(spcsNames, stringsToReplace(i), replacementStrings(i));
end
spcsNames = transpose(spcsNames);
newRowNames = cat(1, cmpNames, spcsNames);
newColNames = cat(1, spcsNames, cmpNames);

%% positive integer values required for circos plots
%exchangeMatrix = exchangeMatrix / min(exchangeMatrix(exchangeMatrix > 0));

fullExchangeTable = array2table(exchangeMatrix, 'VariableNames', cellstr(newColNames), 'RowNames', cellstr(newRowNames));

writetable(fullExchangeTable, strcat(filename, '.txt'), 'WriteRowNames', true, 'Delimiter', '\t')

end

