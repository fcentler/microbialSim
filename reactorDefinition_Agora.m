function [ reactor ] = reactorDefinition_Agora( numberOfSpecies, externalCompounds, dietFile)
%reactorDefinition Summary of this function goes here
%   Detailed explanation goes here

% define the reactor
    reactor.volume = 1.0;	% liter

    reactor.flowRate = 0.0;	% liter/hour

    numberOfCompounds = size(externalCompounds, 1);

    reactor.compounds = externalCompounds';

    reactor.compoundsInit = 0.01 * ones(1, numberOfCompounds);

    reactor.biomassInit = 0.1 * ones(1, numberOfSpecies);

    reactor.compoundsInflow = 0.01 * ones(1, numberOfCompounds);

    reactor.biomassInflow = 0.0 * ones(1, numberOfSpecies);

    if ~isempty(dietFile)
        dietFile = string(dietFile);
        disp(strcat("Loading diet file '", dietFile, "'"))
        dietFileData = readtable(dietFile, "FileType","text",'Delimiter', '\t');
        % remove "EX_" in names
        dietFileData.Reaction = erase(dietFileData.Reaction, 'EX_');
        % reset compoundsInflow (maybe optional?)
        %reactor.compoundsInflow = 0.0 * ones(1, numberOfCompounds);
        for i = 1:size(dietFileData, 1)
            cmpIndex = find(strcmp(externalCompounds, dietFileData{i,"Reaction"}));
            if isempty(cmpIndex)
                disp(strcat("Diet compound '", dietFileData{i,"Reaction"}, "' not utilized by microbiome."))
            else % convert 'mmol/human/day' to inflow concentration
                reactor.compoundsInflow(cmpIndex) = dietFileData{i,"FluxValue"} / 24.0 / reactor.flowRate;
            end
        end
        %reset compoundsInit (maybe optional?)
        reactor.compoundsInit = reactor.compoundsInflow;
    end

end

