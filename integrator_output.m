function [ status ] = integrator_output(t, y, flag, metabolicModel, solverPars)
%INTEGRATOR_OUTPUT Collecting results when using Matlab ODE solver
%   Detailed explanation goes here

    global ODEresult

    switch flag
        case 'init'
            t = t(1);     % first time this function is called, t is [0 tmax](?)
            disp('Starting.');
            status = 0;
            return
        case 'done'
            t = solverPars.tend; % last time this function is called, t is empty
            disp('Done.');
            status = 0;
            return
    end

    for i = 1:length(metabolicModel)
        numberOfReversibleReactions = length(metabolicModel(i).reversibleReacIDs);
        switch metabolicModel(i).FBAsolver
            case 1
                numberOfReactions(i) = length(metabolicModel(i).CNAmodel.reacID) - numberOfReversibleReactions;
            case 2
                numberOfReactions(i) = length(metabolicModel(i).COBRAmodel.rxns) - numberOfReversibleReactions;
        end
    end
    flux = ODEresult.ODElastFlux;
    
    n = ODEresult.counter;
    ODEresult.counter = ODEresult.counter + 1;

    for i = 1:length(metabolicModel)
        ODEresult.metabolicModel(i).lastFlux = flux{i};
        if isfield(metabolicModel(i), 'reversibleReacIDs')
            ODEresult.FBA(i).fluxes(n, :) = flux{i}(1:numberOfReactions(i));
            % add reverse reactions
            for j = 1:length(metabolicModel(i).reversibleReacIDs)
                ODEresult.FBA(i).fluxes(n, metabolicModel(i).reversibleReacIDs(j)) = ODEresult.FBA(i).fluxes(n, metabolicModel(i).reversibleReacIDs(j)) + flux{i}(numberOfReactions(i)+j);
            end
        else
            ODEresult.FBA(i).fluxes(n, :) = flux{i};
        end
        ODEresult.mu(n, i) = ODEresult.FBA(i).fluxes(n, metabolicModel(i).biomassReac);
    end

    fprintf('\nTime %f (%2.1f%% done)\n\n', t(1), t(1) / solverPars.tend * 100);
    
    for k = 2:size(t, 2)
        n = ODEresult.counter;
        ODEresult.counter = ODEresult.counter + 1;

        % to have values at each time step ...
        for i = 1:length(metabolicModel)
            ODEresult.FBA(i).fluxes(n, :) = ODEresult.FBA(i).fluxes(n-1, :);
            ODEresult.mu(n, i) = ODEresult.mu(n-1, i);
        end

        % ... but more correct would be this, but requiring the filtering of
        % negative values during plotting.
%         for i = 1:length(metabolicModel)
%             ODEresult.FBA(i).fluxes(n, :) = -ones(1, size(ODEresult.FBA(i).fluxes, 2));
%             ODEresult.mu(n, i) = -1;
%         end
        
        
        fprintf('\n(Time %f (%2.1f%% done))\n\n', t(k), t(k) / solverPars.tend * 100);
    end
    
    status = 0;
end

