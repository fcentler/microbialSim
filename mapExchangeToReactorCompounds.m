function [ dcompound ] = mapExchangeToReactorCompounds( dcompoundRAW, reactorMapping, numOfReactorCompounds )
%MAPEXCHANGETOREACTORMETABOLITES Summary of this function goes here
%   Detailed explanation goes here
%   dcompound(i,:) = mapExchangeToReactorCompounds(dcompoundRAWE, metabolicModel(i).reactorCompoundIDs);

    dcompound = zeros(1,numOfReactorCompounds);

    for i = 1:length(reactorMapping)
        dcompound(reactorMapping(i)) = dcompoundRAW(i);
    end

end

