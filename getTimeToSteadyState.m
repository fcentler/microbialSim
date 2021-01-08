function [ deltaTime ] = getTimeToSteadyState(metabolicModel, compoundID, metabConc, deltaCompound, deltaFlow, slope, biomass, reactor, solverPars)
%GETTIMETOSTEADYSTATE Summary of this function goes here
%   Detailed explanation goes here

% get production rate & record consuming species

consumers = [];
productionRate = 0;
consumptionRate = 0;

for i = 1:length(metabolicModel)
    if deltaCompound(i) > 0
        productionRate = productionRate + deltaCompound(i);
    else if deltaCompound(i) < 0
            consumers(end + 1) = i;
            consumptionRate = consumptionRate + deltaCompound(i);
        end
    end
end

% for each consuming species, get index of coupled reaction

reacIndexStr = cell(1,1);

for i = 1:length(consumers)
    tmp = find(~(metabolicModel(consumers(i)).reactorCompoundIDs - compoundID)); % find index at which position value compoundID is located (finding the index of the coupled reaction which refers to the reactor compound)
    reacIndexStr{i} = int2str(metabolicModel(consumers(i)).coupledReactions{tmp, 'ReacID'});
end

% get steady state concentration

if productionRate + max(deltaFlow, 0) == 0 % only consuming processes present
    targetConcentration = 0;
else % consumption & producing processes present
    topBound = metabConc;
    lowBound = 0;
    consumptionRate = - productionRate - 1000; % just to ensure to enter the while loop
    oldConsumptionRate = Inf;
    while abs(productionRate + consumptionRate) > solverPars.SteadyStateAccuracy
        targetConcentration = (topBound + lowBound) / 2;
        consumptionRate = 0;
        for i = 1:length(consumers)
            consRate = getMaxUptakeFlux(metabolicModel(consumers(i)), reacIndexStr(i), targetConcentration);
            consRate = consRate * biomass(consumers(i));
            consumptionRate = consumptionRate + consRate;
        end
        consumptionRate = consumptionRate + considerFlow(reactor.compoundsInflow(compoundID), targetConcentration, reactor);
        if productionRate + consumptionRate < 0
            topBound = targetConcentration;
        else
            lowBound = targetConcentration;
        end
        if oldConsumptionRate == consumptionRate % rate does not depend on concentration, e.g. if Km = 0!
            targetConcentration = 0;
            break
        else
            oldConsumptionRate = consumptionRate;
        end
    end
end
% compute deltaTime at which steady state concentration is reached
deltaTime = (targetConcentration - metabConc) / slope;
end
