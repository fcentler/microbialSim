function [success] = writeStateToFile(recordFilename, n, result, metabolicModel, stronglyConsumedFlag, isLastRecord, solverPars, workerModels)
%WRITESTATETOFILE Summary of this function goes here
%   Detailed explanation goes here
    recordResult.time = result.time(n);
    recordResult.compounds = result.compounds(n,:);
    recordResult.biomass = result.biomass(n,:);
    recordResult.stronglyConsumedFlag = stronglyConsumedFlag;
    if solverPars.parallel == 1
        for i = 1:numel(workerModels)
            tmpWorker = workerModels{i};
            for j = 1:numel(tmpWorker)
                metabolicModel(i+numel(workerModels)*(j-1)).lastFlux = tmpWorker(j).lastFlux;
                metabolicModel(i+numel(workerModels)*(j-1)).previousFlux = tmpWorker(j).previousFlux;
            end
        end
    end
    
    for m = 1:length(metabolicModel)
        if exist('result.FBA')
            recordResult.FBA(m).fluxes = result.FBA(m).fluxes(n,:);
        end
        if isLastRecord
            recordResult.metabolicModel(m).lastFlux = metabolicModel(m).lastFlux;
        else
            recordResult.metabolicModel(m).lastFlux = metabolicModel(m).previousFlux;
        end
    end
    save(recordFilename, 'recordResult');
    clear recordResult
    success = 1;
end
