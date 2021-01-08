function [ dydt ] = getDyWrapper(y, metabolicModel, reactor, solverPars, workerModels)
%GETDYWRAPPER This is a wrapper function for getDy which is used when the
% Matlab ODE solver are used
%   Detailed explanation goes here

global ODEresult

numberOfSpecies = length(metabolicModel);

biomass = y(1:numberOfSpecies);
compounds = y((numberOfSpecies+1):end);

if solverPars.parallel == 1
    for i = 1:length(ODEresult.metabolicModel)
        myODEresult_metabolicModel_lastFlux{i} = ODEresult.metabolicModel(i).lastFlux;
    end
    spmd
        for i = 1:length(workerModels)
            workerModels(i).lastFlux = myODEresult_metabolicModel_lastFlux{labindex+(i-1)*numlabs};
        end
    end
else
    for i = 1:numberOfSpecies
        metabolicModel(i).lastFlux = ODEresult.metabolicModel(i).lastFlux;
    end
end

if solverPars.parallel == 1
    %ticBytes(gcp);
    spmd
        for i = 1:length(workerModels)
            [mydbio(i), mydcompoundRAW{i}, myflux{i}, mysuccess(i)] = getDy(workerModels(i), transpose(compounds), biomass(labindex + numlabs * (i-1)), solverPars);
        end
    end
    %tocBytes(gcp)
    for i = 1:length(mydbio)
        dbio(i:length(mydbio):length(metabolicModel)) = mydbio{i};
        dcompoundRAW(i:length(mydbio):length(metabolicModel)) = mydcompoundRAW{i};
        flux(i:length(mydbio):length(metabolicModel)) = myflux{i};
        success(i:length(mydbio):length(metabolicModel)) = mysuccess{i};
    end
else
    for i = 1:numberOfSpecies
        [ dbio(i), dcompoundRAW{i}, flux{i}, success(i) ] = getDy(metabolicModel(i), transpose(compounds), biomass(i), solverPars);
    end
end

for i = 1:numberOfSpecies
    dcompound(i,:) = mapExchangeToReactorCompounds(dcompoundRAW{i}, metabolicModel(i).reactorCompoundIDs, length(reactor.compounds));
    % consider flow (biomasses)
    dbio(i) = dbio(i) + considerFlow(reactor.biomassInflow(i), biomass(i), reactor);
end

ODEresult.ODElastFlux = flux;

% consider flow (compounds)
if numberOfSpecies > 1
    myDcompound = sum(dcompound);
else
    myDcompound = dcompound;
end

dcompound = myDcompound + considerFlow(reactor.compoundsInflow, compounds', reactor);

dydt = transpose([dbio dcompound]);

end
