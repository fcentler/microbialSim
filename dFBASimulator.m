function [ result ] = dFBASimulator( metabolicModel, reactor, solverPars, workerModels)
%DFBASIMULATOR Running a dynamic FBA simulation with µbialSim
%   Running the dFBA Simulation. For solverPars.solverType == 0, the step-wise iteration is
%   used, solverPars.solverType == 1 uses the Matlab ODE solver.

% PREPARE DATA

global ODEresult

initFile = solverPars.readInitialStateFrom;

if isempty(initFile)
    startTime = 0;
    stronglyConsumedFlag=ones(length(reactor.compounds),1)*0; % this flags compounds that in one time step were more consumed than produced leading to time step size reduction to avoid negative concentrations
else
    myInit = load(initFile);
    if solverPars.readInitialStateResetTime
        startTime = 0;
    else
        startTime = myInit.recordResult.time;
    end
    if solverPars.readInitialStateSelection > 0
        stronglyConsumedFlag = myInit.recordResult.stronglyConsumedFlag;
        if solverPars.parallel == 1
            for i = 1:length(myInit.recordResult.metabolicModel)
                myInit_recordResult_metabolicModel_lastFlux{i} = myInit.recordResult.metabolicModel(i).lastFlux;
            end
            spmd
                for i = 1:length(workerModels)
                    workerModels(i).lastFlux = myInit_recordResult_metabolicModel_lastFlux{labindex+(i-1)*numlabs};
                end
            end
        else
            for m = 1:length(metabolicModel)
                metabolicModel(m).lastFlux = myInit.recordResult.metabolicModel(m).lastFlux;
            end
        end
    else
        stronglyConsumedFlag=ones(length(reactor.compounds),1)*0;
    end
end

timeSteps = floor((solverPars.tend - startTime) / solverPars.timeStepSize);
numberOfCompounds = length(reactor.compounds);

if length(reactor.biomassInit) ~= length(reactor.biomassInflow) || length(metabolicModel) ~= length(reactor.biomassInit)
    error('Number of biomasses, inflow, and FBA models mismatches!');
end

% iterate over all FBA models

for i = 1:length(metabolicModel)
    numberOfReversibleReactions = length(metabolicModel(i).reversibleReacIDs);
    switch metabolicModel(i).FBAsolver
        case 1
            numberOfReactions(i) = length(metabolicModel(i).CNAmodel.reacID) - numberOfReversibleReactions;
        case 2
            numberOfReactions(i) = length(metabolicModel(i).COBRAmodel.rxns) - numberOfReversibleReactions;
            % Flux reporting will be done based on number of reactions in
            % original model.
    end
    result.FBA(i).fluxes = zeros(timeSteps, numberOfReactions(i));
    if isempty(initFile) || solverPars.readInitialStateSelection < 1
        metabolicModel(i).lastFlux=[];
    end
end

if solverPars.parallel == true
    spmd
        for i = 1:length(workerModels)
            workerModels(i).lastFlux=[];
        end
    end
end

if isempty(initFile)
    result.time(1) = 0;
    result.compounds(1,:) = reactor.compoundsInit;
    result.biomass(1,:) = reactor.biomassInit;
else
   if solverPars.readInitialStateResetTime
       result.time(1) = 0;
   else 
       result.time(1) = myInit.recordResult.time;
   end
   if solverPars.readInitialStateSelection ~= 1
       if solverPars.readInitialStateSelection == 2
           result.compounds(1,:) = myInit.recordResult.compounds;
       else % a change of community composition might lead to a different number / order of chemical compounds, careful mapping required
           result.compounds(1,:) = reactor.compoundsInit;
	   for i = 1:length(myInit.recordResult.compounds)
	       found = 0;
	       for j = 1:length(reactor.compounds)
		   if strcmp(myInit.recordResult.compoundNames{i}, reactor.compounds(j))
		       result.compounds(1,j) = myInit.recordResult.compounds(i);
                       found = 1;
		       break
		   end
	       end
	       if found ~= 1
		       error("Compound name from init file not found in current reactor!")
	       end
	   end

       end
   else
       result.compounds(1,:) = reactor.compoundsInit;
   end
   if solverPars.readInitialStateSelection > 0
       result.biomass(1,:) = myInit.recordResult.biomass;
   else
       result.biomass(1,:) = reactor.biomassInit;
   end
end

if (solverPars.recordLimitingFluxes == 1)
    result.limitingFluxes = cell2table(cell(0,9), 'VariableNames', {'ModelID', 'ModelName', 'Time', 'FluxID', 'FluxName', 'FluxValue', 'IsFixed', 'Compound', 'ExchangeCoupledFlux'});
end    

% DO SIMULATION

if solverPars.solverType == 1
    ODEresult = result;
    ODEresult.counter = 1;
    for i = 1:length(metabolicModel)
        ODEresult.metabolicModel(i).lastFlux=[];
    end

    options = odeset('OutputFcn', @(t,y,flag) integrator_output(t,y,flag,metabolicModel,solverPars), 'RelTol', solverPars.relTol, 'AbsTol', solverPars.absTol);
    %ODE Simulation
    if solverPars.nonNegative == 1
        options = odeset(options, 'NonNegative', 1:1:(length(metabolicModel)+numberOfCompounds));
    end
    switch solverPars.solver
        case 1
            [fba.t,fba.y] = ode45(@(t,y) getDyWrapper(y,metabolicModel,reactor, solverPars, workerModels), [0:solverPars.timeStepSize:solverPars.tend], [reactor.biomassInit reactor.compoundsInit], options);
        case 2
            [fba.t,fba.y] = ode15s(@(t,y) getDyWrapper(y,metabolicModel,reactor, solverPars, workerModels), [0:solverPars.timeStepSize:solverPars.tend], [reactor.biomassInit reactor.compoundsInit], options);
        otherwise
            error('Unknown solver.')
    end
    result = ODEresult;
    result.biomass = fba.y(:,1:length(metabolicModel));
    result.compounds = fba.y(:,(length(metabolicModel)+1:end));
    result.time(2:(timeSteps+1)) = linspace(solverPars.timeStepSize, solverPars.tend, timeSteps)';
else
    n = 1;

    myTimeStepSize = solverPars.timeStepSize;

    solveFullFBA = 1;

    while result.time(n) < solverPars.tend
        % get FBA predictions for all models
        result.compounds(n+1, :) = result.compounds(n, :);
        result.biomass(n+1, :) = result.biomass(n, :);
        if solveFullFBA == 1
            clearvars dbio dcompound flux success 
            if solverPars.parallel == 1
                %ticBytes(gcp);
                myResult.compounds = result.compounds(n, :);
                myResult.biomass = result.biomass(n, :);
                spmd
                    for i = 1:length(workerModels)
                        [mydbio(i), mydcompoundRAW{i}, myflux{i}, mysuccess(i), mylimitedFluxes{i}] = getDy(workerModels(i), myResult.compounds, myResult.biomass(labindex + numlabs * (i-1)), solverPars);
                    end
                end
                %tocBytes(gcp)
                for i = 1:length(mydbio)
                    dbio(i:length(mydbio):length(metabolicModel)) = mydbio{i};
                    dcompoundRAW(i:length(mydbio):length(metabolicModel)) = mydcompoundRAW{i};
                    flux(i:length(mydbio):length(metabolicModel)) = myflux{i};
                    success(i:length(mydbio):length(metabolicModel)) = mysuccess{i};
                    limitedFluxes(i:length(mydbio):length(metabolicModel)) = mylimitedFluxes{i};
                end
            else
                for i = 1:length(metabolicModel)
                    [dbio(i), dcompoundRAW{i}, flux{i}, success(i), limitedFluxes{i}] = getDy(metabolicModel(i), result.compounds(n, :), result.biomass(n, i), solverPars);
                end
            end
            
            for i = 1:length(metabolicModel)
                dcompound(i,:) = mapExchangeToReactorCompounds(dcompoundRAW{i}, metabolicModel(i).reactorCompoundIDs, length(reactor.compounds));
            end

            if solverPars.parallel == 1
                spmd
                    for i = 1:length(workerModels)
                        workerModels(i).previousFlux = workerModels(i).lastFlux;
                        workerModels(i).lastFlux = flux{labindex+(i-1)*numlabs};
                    end
                end
            else
                for i = 1:length(metabolicModel)
                    metabolicModel(i).previousFlux = metabolicModel(i).lastFlux;
                    metabolicModel(i).lastFlux = flux{i};
                end
            end
            if length(metabolicModel) > 1
                myDcompound = sum(dcompound);
            else
                myDcompound = dcompound;
            end

            dbioTotal = dbio + considerFlow(reactor.biomassInflow, result.biomass(n, :), reactor);
            dcompoundFlow = considerFlow(reactor.compoundsInflow, result.compounds(n, :), reactor);
            dcompoundTotal = myDcompound + dcompoundFlow;

            % record limiting fluxes
            if(solverPars.recordLimitingFluxes == 1)
                for i = 1:length(metabolicModel)
                    for j = 1:size(limitedFluxes{i}, 1)
                    %'ModelID', 'ModelName', 'Time', 'FluxID', 'FluxName',
                    %'FluxValue', 'IsFixed', 'Compound', 'Exchange/Coupled Flux'
                        newEntry = {};
                        newEntry{1,1} = i;
                        newEntry{1,2} = metabolicModel(i).modelName;
                        newEntry{1,3} = result.time(n);
                        newEntry{1,4} = limitedFluxes{i}(j,1);
                        switch metabolicModel(i).FBAsolver
                        case 1
                            newEntry{1,5} = metabolicModel(i).CNAmodel.reacID(limitedFluxes{i}(j, 1));
                        case 2
                            newEntry{1,5} = metabolicModel(i).COBRAmodel.rxnNames(limitedFluxes{i}(j, 1));
                        end
                        newEntry{1,6} = limitedFluxes{i}(j,2);
                        newEntry{1,7} = limitedFluxes{i}(j,3);
						if ismember(newEntry{1,4}, metabolicModel(i).exchangeReactions.ReacID)
						    newEntry{1,9} = 1;
						else
						    newEntry{1,9} = 0;
						end
                        myIndex = find(newEntry{1,4} == metabolicModel(i).coupledReactions.ReacID);
                        if length(myIndex) == 1 
                            newEntry{1,8} = reactor.compounds(metabolicModel(i).reactorCompoundIDs(myIndex));
                            newEntry{1,9} = 2;
						else
                            newEntry{1,8} = '-';
                        end
						
                        result.limitingFluxes = [result.limitingFluxes;newEntry];
                    end
                end
            end
        end
        
        %Update reactor with changes due to biotic activity and flow
        result.biomass(n+1, :) = result.biomass(n+1, :) + dbioTotal * myTimeStepSize;
        result.compounds(n+1, :) = result.compounds(n+1, :) + dcompoundTotal * myTimeStepSize;

        % check if time step size is ok or needs to be reduced due to negative
        % concentrations, or as exchanged, strongly consumed compounds
        % deviate too much

        if solveFullFBA == 1 && solverPars.augmentedEuler == 1
            % do negative biomass concentrations occur?
            for i = 1:length(metabolicModel)
                if result.biomass(n+1, i) <= 0 && result.biomass(n, i) ~= 0
                    deltaTimeToZero = (0 - result.biomass(n, i)) / ((result.biomass(n+1, i) - result.biomass(n, i))/solverPars.timeStepSize);
                    if abs(result.biomass(n+1, i)) < solverPars.myBioAccuracy
                        requiredTimeStepSize = deltaTimeToZero;
                    else
                        requiredTimeStepSize = deltaTimeToZero / solverPars.biomassReductionFactor;
                    end
                    if requiredTimeStepSize < myTimeStepSize
                        myTimeStepSize = requiredTimeStepSize;
                    end
                end
            end
            
            % do negative compound concentrations occur?
            for i = 1:numberOfCompounds
                if result.compounds(n+1, i) <= 0 && result.compounds(n, i) ~= 0
                    if result.compounds(n+1, i) < 0 && result.compounds(n, i) <= 0
                        error('In TimeStepAdjustment: previous concenctration was zero or negative! How was uptake possible? Something is wrong. Aborting')
                    end
                    deltaConcentration = result.compounds(n+1, i) - result.compounds(n, i);
                    if result.compounds(n+1, i) < 0 && deltaConcentration == 0
                        error('In TimeStepAdjustment: deltaConcentration is zero. I am confused. Aborting')
                    end
                    requiredTimeStepSize = getTimeToSteadyState(metabolicModel, i, result.compounds(n, i), dcompound(:,i), dcompoundFlow(i), deltaConcentration/solverPars.timeStepSize, result.biomass(n,:), reactor, solverPars);
                    if requiredTimeStepSize < myTimeStepSize
                        myTimeStepSize = requiredTimeStepSize;
                    end
                    stronglyConsumedFlag(i) = 1;
                else
                % avoid fluctuating concentrations for compounds produced and
                % strongly consumed
                    if stronglyConsumedFlag(i) == 1 && result.compounds(n, i) ~= 0
                       if abs(result.compounds(n+1,i) - result.compounds(n, i)) > result.compounds(n,i)/100 * solverPars.maxDeviation;
                           requiredTimeStepSize = result.compounds(n, i) / 100.0 * solverPars.maxDeviation / abs(result.compounds(n+1,i)-result.compounds(n,i)) * solverPars.timeStepSize;
                            if requiredTimeStepSize < myTimeStepSize
                                myTimeStepSize = requiredTimeStepSize;
                            end
% Uncommenting the next two lines speeds up computation, but leads to non-smooth compound trajectories. Runtime and accuracy will be in the middle between
%     setting solverPars.maxDeviation to "5.0" (or other finite values) (slow computation, high accurracy) and "inf" (fast computation, low accurracy).    
%                       else
%                           stronglyConsumedFlag(i) = 0;
                       end
                    end
                end
                if stronglyConsumedFlag(i) == 1
                    if min(dcompound(:,i)) >= 0 || (dcompoundFlow(i) <= 0 && max(dcompound(:,i)) <= 0)
                        stronglyConsumedFlag(i) = 0;
                    end
                end
            end
            
            if myTimeStepSize ~= solverPars.timeStepSize
		if solverPars.logLevel > 1
                    fprintf('Adjusting time step size')
		end
                solveFullFBA = 0;
                continue
                % recalculate results
            end
        end
        
	if solverPars.augmentedEuler == 1
		% for small biomass and compound concentrations: set to zero
		for i = 1:numberOfCompounds
		    if abs(result.compounds(n+1, i)) < solverPars.myAccuracy
			result.compounds(n+1, i) = 0.0;
		    end
		end
		for i = 1:length(metabolicModel)
		    if abs(result.biomass(n+1, i)) < solverPars.myBioAccuracy
			result.biomass(n+1, i) = 0.0;
		    end
		end
	else
		% setting negative concentrations to zero
		for i = 1:numberOfCompounds
			result.compounds(n+1, i) = max(result.compounds(n+1, i), 0);
		end
		for i = 1:length(metabolicModel)
			result.biomass(n+1, i) = max(result.biomass(n+1, i), 0);
		end
		
	end	
        
        % STORE FLUX DATA
        for i = 1:length(metabolicModel)
            if isfield(metabolicModel(i), 'reversibleReacIDs')
                result.FBA(i).fluxes(n, :) = flux{i}(1:numberOfReactions(i));
                % add reverse reactions
                for j = 1:length(metabolicModel(i).reversibleReacIDs)
                    result.FBA(i).fluxes(n, metabolicModel(i).reversibleReacIDs(j)) = result.FBA(i).fluxes(n, metabolicModel(i).reversibleReacIDs(j)) + flux{i}(numberOfReactions(i)+j);
                end
            else
                result.FBA(i).fluxes(n, :) = flux{i};
            end
            %store current specific growth rate µ
            result.mu(n, i) = result.FBA(i).fluxes(n, metabolicModel(i).biomassReac);
        end

        result.time(n+1) = result.time(n) + myTimeStepSize;
        if solverPars.recording == 1
            recordFilename = strcat(solverPars.trajectoryFile, '_', int2str(n),'.mat');
            writeStateToFile(recordFilename, n, result, metabolicModel, stronglyConsumedFlag, 0, solverPars, workerModels);
        end
        n = n + 1;
        solveFullFBA = 1;
        myTimeStepSize = solverPars.timeStepSize;
	switch solverPars.logLevel
	    case 1
                if (floor((result.time(n)-result.time(1)) / (solverPars.tend-result.time(1)) * 100) - floor((result.time(n-1)-result.time(1)) / (solverPars.tend-result.time(1)) * 100)) >= 1
                    fprintf('Time %f (%2.1f%% done)\n', result.time(n), (result.time(n)-result.time(1)) / (solverPars.tend-result.time(1)) * 100);
		        end
	    case 2
                fprintf('Time %f (%2.1f%% done)\n', result.time(n), (result.time(n)-result.time(1)) / (solverPars.tend-result.time(1)) * 100);
        otherwise
            error("unkown logLeve, aborting");
	end
    end
    if solverPars.recording == 1 % record final step, there are no fluxes then!
        recordFilename = strcat(solverPars.trajectoryFile, '_', int2str(n),'.mat');
        writeStateToFile(recordFilename, n, result, metabolicModel, stronglyConsumedFlag, 1, solverPars, workerModels);
    end

    % record restart file
    result.compoundNames = reactor.compounds;
    result.modelNames = {};
    for i = 1:length(metabolicModel)
        if ischar(metabolicModel(i).modelName)
            metabolicModel(i).modelName = cellstr(metabolicModel(i).modelName);
        end
        result.modelNames{i} = metabolicModel(i).modelName;
        % add reaction names (only tested for COBRAToolbox); expecting all
        % reversible reactions at end of reaction list!
        if (solverPars.FBAsolver == 2)
            result.FBA(i).fluxes = array2table(result.FBA(i).fluxes, "VariableNames", metabolicModel(i).COBRAmodel.rxns(1:(length(metabolicModel(i).COBRAmodel.rxns)-length(metabolicModel(i).reversibleReacIDs))));
        end
    end
    result.modelNames = [result.modelNames{:}];
    recordFilename = strcat(solverPars.trajectoryFile, '_restartInit.mat');
    writeStateToFile(recordFilename, n, result, metabolicModel, stronglyConsumedFlag, 1, solverPars, workerModels);
end

result.time = transpose(result.time);

for i = 1:length(metabolicModel)
    result.FBA(i).biomassReac = metabolicModel(i).biomassReac;
    result.FBA(i).ngamReac = metabolicModel(i).ngamReac;
    result.FBA(i).exchangeReactions = metabolicModel(i).exchangeReactions;
    result.FBA(i).coupledReactions = metabolicModel(i).coupledReactions;
    result.FBA(i).reactorCompoundIDs = metabolicModel(i).reactorCompoundIDs;
    result.FBA(i).modelName = metabolicModel(i).modelName;
end

% do mass balance for all exchange fluxes
if solverPars.doMassBalance == 1
    if solverPars.parallel == true
        modelBiomass = result.biomass;
        modelTime = result.time;
%        ticBytes(gcp);
        modelFluxes = Composite();
        modelBiomass = Composite();
        for i = 1:numel(modelFluxes)
            modelFluxes{i} = cat(2, result.FBA(i:numel(modelFluxes):length(metabolicModel)).fluxes);
            modelBiomass{i} = result.biomass(:,i:numel(modelFluxes):length(metabolicModel));
        end
        spmd
            fluxStartIndex = 1;
            for i = 1:length(workerModels)
                fluxEndIndex = fluxStartIndex - 1 + length(workerModels(i).lastFlux) - length(workerModels(i).reversibleReacIDs);
                for j = 1:height(workerModels(i).exchangeReactions)
                    myMassBalance(i, j) = 0;
                    for k = 1:size(modelFluxes, 1)
                        myMassBalance(i, j) = myMassBalance(i, j) + modelFluxes(k, fluxStartIndex - 1 + workerModels(i).exchangeReactions{j, 'ReacID'}) * modelBiomass(k, i) * (modelTime(k+1)-modelTime(k));                
                    end
                end
                abs_sum(i) = sum(sum(abs(modelFluxes(:,fluxStartIndex:fluxEndIndex))));
                dev_sum(i) = 0;
                for j = 2:size(modelFluxes, 1)
                    dev_sum(i) = dev_sum(i) + sum(abs(modelFluxes(j,fluxStartIndex:fluxEndIndex) - modelFluxes(j-1,fluxStartIndex:fluxEndIndex)));
                end
                fluxStartIndex = fluxEndIndex + 1;
            end
        end
 %       tocBytes(gcp);
        for i = 1:length(myMassBalance)
            tmp = abs_sum{i};
            tmp2 = dev_sum{i};
            k = 1;
            for j = i:length(myMassBalance):length(metabolicModel)
                result.FBA(j).flux_sum = tmp(k);
                result.FBA(j).flux_dev_sum = tmp2(k);
                k = k + 1;
            end
            myMB = myMassBalance{i};
            myIdx = 1;
            for j = i:length(myMassBalance):length(metabolicModel)
                result.FBA(j).exchangeReactions = addvars(result.FBA(j).exchangeReactions, transpose(myMB(myIdx, 1:height(result.FBA(j).exchangeReactions))), 'NewVariableNames', 'MassBalance');
                myIdx = myIdx + 1;
            end
        end
    else
        for i = 1:length(metabolicModel)    % for all models
            result.FBA(i).exchangeReactions.MassBalance = zeros(height(result.FBA(i).exchangeReactions), 1); % init with zeros
            for j = 1:height(result.FBA(i).exchangeReactions)   % for all exchange reactions
                for k = 1:size(result.FBA(i).fluxes, 1)  % for all time steps / recorded rates
                    result.FBA(i).exchangeReactions{j, 'MassBalance'} = result.FBA(i).exchangeReactions{j, 'MassBalance'} + result.FBA(i).fluxes(k, result.FBA(i).exchangeReactions{j, 'ReacID'}) * result.biomass(k, i) * (result.time(k+1)-result.time(k));
                end
            end

            % compute absolute fluxes & flux difference
            result.FBA(i).flux_sum = sum(sum(abs(result.FBA(i).fluxes)));
            result.FBA(i).flux_dev_sum = 0;
            for j = 2:size(result.FBA(i).fluxes,1)
                result.FBA(i).flux_dev_sum = result.FBA(i).flux_dev_sum + sum(abs(result.FBA(i).fluxes(j,:) - result.FBA(i).fluxes(j-1,:)));
            end
        end
    end
    for i = 1:length(metabolicModel)    % for all models
        % copy data to coupledReactions table
        result.FBA(i).coupledReactions.MassBalance = zeros(height(result.FBA(i).coupledReactions), 1); % init with zeros
        for j = 1:height(result.FBA(i).coupledReactions)
            result.FBA(i).coupledReactions{j, 'MassBalance'} = result.FBA(i).exchangeReactions{int2str(result.FBA(i).coupledReactions{j, 'ReacID'}), 'MassBalance'};
        end

        % Add final flux info to reaction tables
        result.FBA(i).exchangeReactions.finalFlux = getFluxes(result.FBA(i).fluxes(size(result.FBA(i).fluxes, 1),:), result.FBA(i).exchangeReactions.ReacID)';
        result.FBA(i).coupledReactions.finalFlux = getFluxes(result.FBA(i).fluxes(size(result.FBA(i).fluxes, 1),:), result.FBA(i).coupledReactions.ReacID)';

        fluxLimits = getLimitedFluxes(result.FBA(i).fluxes(size(result.FBA(i).fluxes, 1),:), metabolicModel(i), solverPars, result.compounds(size(result.FBA(i).fluxes, 1),:));

        result.FBA(i).exchangeReactions.finalFluxLimit = fluxLimits(result.FBA(i).exchangeReactions.ReacID)';
        result.FBA(i).coupledReactions.finalFluxLimit = fluxLimits(result.FBA(i).coupledReactions.ReacID)';
    end
end
result.solverPars = solverPars;
if solverPars.solverType == 1
    result.compoundNames = reactor.compounds;
    result.modelNames = {};
    for i = 1:length(metabolicModel)
        if ischar(metabolicModel(i).modelName)
            metabolicModel(i).modelName = cellstr(metabolicModel(i).modelName);
        end
        result.modelNames{i} = metabolicModel(i).modelName;
        % add reaction names (only tested for COBRAToolbox); expecting all
        % reversible reactions at end of reaction list!
        if (solverPars.FBAsolver == 2)
            result.FBA(i).fluxes = array2table(result.FBA(i).fluxes, "VariableNames", metabolicModel(i).COBRAmodel.rxns(1:(length(metabolicModel(i).COBRAmodel.rxns)-length(metabolicModel(i).reversibleReacIDs))));
        end
    end
    result.modelNames = [result.modelNames{:}];
end

% compute final derivatives
for i = 1:length(metabolicModel)
    result.finalGradientX(i) = (result.biomass(end, i) - result.biomass(end-1, i)) / (result.time(end) - result.time(end-1)); % gDW/l/h
end
for i = 1:numberOfCompounds
    result.finalGradientC(i) = (result.compounds(end, i) - result.compounds(end-1, i)) / (result.time(end) - result.time(end-1)); % mM/h
end

if solverPars.parallel == 1
    p = gcp('nocreate');
    result.numberOfWorkers = p.NumWorkers;
end
end

