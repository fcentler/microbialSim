function [ dbio, dcompound, flux, success, resultLimitingFluxes ] = getDy(metabolicModel, compounds, biomass, solverPars)
%GETDY Compute the current rate of change for biomass and reactor compounds
%based on the current reactor state
%   Detailed info

    % Check if species still alive
    if biomass == 0
	if solverPars.logLevel > 1
            fprintf('Species died out.');
	end
        switch metabolicModel.FBAsolver
            case 1
                flux = zeros(length(metabolicModel.CNAmodel.objFunc), 1);
            case 2
                flux = zeros(length(metabolicModel.COBRAmodel.c), 1);
        end
        success = 1;
        dbio = 0.0;
        dcompound = zeros(1, length(metabolicModel.coupledReactions.ReacID));
        resultLimitingFluxes = '';
        return
    end

    % SET UPTAKE LIMITS
    for i = 1:length(metabolicModel.coupledReactions.ReacID)
        limit = getMaxUptakeFlux(metabolicModel, i, compounds(metabolicModel.reactorCompoundIDs(i)));
        if metabolicModel.coupledReactions{i, 'SecretionSense'} > 0
            switch metabolicModel.FBAsolver
                case 1
                    metabolicModel.CNAmodel.reacMin(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
                case 2
                    metabolicModel.COBRAmodel.lb(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
            end
        else
            switch metabolicModel.FBAsolver
                case 1
                    metabolicModel.CNAmodel.reacMax(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
                case 2
                    metabolicModel.COBRAmodel.ub(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
            end
        end
    end
  
    % FBA
    switch metabolicModel.FBAsolver
        case 1
            [flux, success, status] = CNAoptimizeFlux(metabolicModel.CNAmodel, metabolicModel.CNAconstraints);
            if success
                flux = validateFlux(flux, metabolicModel.CNAmodel.reacMin, metabolicModel.CNAmodel.reacMax);
            end
        case 2
            solution = optimizeCbModel(metabolicModel.COBRAmodel, 'max', 0, 1);
            
            flux = solution.x;
            success = solution.stat == 1;
            if success
                flux = validateFlux(flux, metabolicModel.COBRAmodel.lb, metabolicModel.COBRAmodel.ub);
            end
        otherwise
            error('In getDy: unknown FBA solver')
    end
    
    % not a simulation, but just called by estimateVmax.m
    if biomass == -1
        dbio = 0;
        dcompound = 0;
        return
    end
    
    if success == false || flux(metabolicModel.biomassReac) < solverPars.minimalGrowth
	if solverPars.logLevel > 1
            fprintf('No growth conditions.');
	end
        switch metabolicModel.FBAsolver
            case 1
                flux = zeros(length(metabolicModel.CNAmodel.objFunc), 1);
            case 2
                flux = zeros(length(metabolicModel.COBRAmodel.c), 1);
        end
    else
        
% perform parsimonious FBA, fix growth rate and minimize internal fluxes,
% however not allowing reactions to switch direction
        originalSum=(sum(abs(flux)));
        if solverPars.dopFBA
            % record current constraints etc.
            switch metabolicModel.FBAsolver
                case 1
                    currentObj = metabolicModel.CNAmodel.objFunc;
                    currentConstraints = metabolicModel.CNAconstraints;
                    currentLB = metabolicModel.CNAmodel.reacMin;
                    currentUB = metabolicModel.CNAmodel.reacMax;
                    numberOfReactions = length(metabolicModel.CNAmodel.objFunc);
                case 2
                    currentObj = metabolicModel.COBRAmodel.c;
                    currentLB = metabolicModel.COBRAmodel.lb;
                    currentUB = metabolicModel.COBRAmodel.ub;
                    numberOfReactions = length(metabolicModel.COBRAmodel.c);
            end
            % prep objective function    
            for i = 1:numberOfReactions
                if (~ismember(i, metabolicModel.exchangeReactions.ReacID))
                    if currentLB(i) * currentUB(i) < 0
                        error('pFBA only works if (internal) reactions are all irreversible.')
                    end
                    if max(currentLB(i), currentUB(i)) > 0 && currentLB(i) ~= currentUB(i)
                        switch metabolicModel.FBAsolver
                            case 1
                                metabolicModel.CNAmodel.objFunc(i) = 1;
                            case 2
                                metabolicModel.COBRAmodel.c(i) = 1;
                        end
                    else
                        if min(currentLB(i), currentUB(i)) < 0 && currentLB(i) ~= currentUB(i)
                            switch metabolicModel.FBAsolver
                                case 1
                                    metabolicModel.CNAmodel.objFunc(i) = -1;
                                case 2
                                    metabolicModel.COBRAmodel.c(i) = -1;
                            end
                        else % flux restricted to zero or fixed value
                            switch metabolicModel.FBAsolver % ignore reactions which carry no flux at optimal growth
                                case 1
                                    metabolicModel.CNAconstraints(i) = flux(i);
                                    metabolicModel.CNAmodel.objFunc(i) = 0;
                                case 2
                                    metabolicModel.COBRAmodel.lb(i) = flux(i);
                                    metabolicModel.COBRAmodel.ub(i) = flux(i);
                                    metabolicModel.COBRAmodel.c(i) = 0; 
                            end
                        end
                    end
                end
            end
            
            % fix biomass reaction, do simulation, and reset objective
            % function
            for k = 1:numberOfReactions
                if currentObj(k) ~= 0
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAconstraints(k) = flux(k);
                            metabolicModel.CNAmodel.objFunc(k) = 0;
                        case 2
                            metabolicModel.COBRAmodel.lb(k) = flux(k);
                            metabolicModel.COBRAmodel.ub(k) = flux(k);
                            metabolicModel.COBRAmodel.c(k) = 0;
                    end
                end
            end
            switch metabolicModel.FBAsolver
                case 1  
                        [flux, success, status] = CNAoptimizeFlux(metabolicModel.CNAmodel, metabolicModel.CNAconstraints);
                        if success == false
                            flux = zeros(length(metabolicModel.CNAmodel.objFunc), 1);
			    if solverPars.logLevel > 1
                                fprintf('Whoa, solution in minimization did not work out. Setting to no growth conditions.')
			    end
                        end
                        metabolicModel.CNAmodel.objFunc = currentObj;
                        metabolicModel.CNAconstraints = currentConstraints;
                case 2
                        solution = optimizeCbModel(metabolicModel.COBRAmodel, 'min', 0, 1);
                        if (solution.stat ~= 1)
                            flux = zeros(length(metabolicModel.COBRAmodel.c), 1);
			    if solverPars.logLevel > 1
                                fprintf('Whoa, solution in minimization did not work out. Setting to no growth conditions.')
			    end
                        else
                            flux = solution.x;
                        end
                        metabolicModel.COBRAmodel.c = currentObj;
                        metabolicModel.COBRAmodel.lb = currentLB;
                        metabolicModel.COBRAmodel.ub = currentUB;
            end
            my_diff = originalSum - sum(abs(flux));
            if my_diff > 0
		if solverPars.logLevel > 1
                    fprintf('pFBA improvement (%f)\n', my_diff)
		end
            else
		if solverPars.logLevel > 1
                    fprintf('pFBA failure (%f)\n', my_diff)
		end
            end
            flux = validateFlux(flux, currentLB, currentUB);
        end
        
        % Minimize deviation to previous flux
        
        if (solverPars.doMin2PrevFlux == 1 && sum(abs(metabolicModel.lastFlux)) ~= 0) % previous flux available, try to minimize deviation
            originalDev = sum(abs(flux-metabolicModel.lastFlux));
            % retain current state
            switch metabolicModel.FBAsolver
                case 1
                    currentM = metabolicModel.CNAmodel.stoichMat;
                    currentUB = metabolicModel.CNAmodel.reacMax;
                    currentLB = metabolicModel.CNAmodel.reacMin;
                    currentObj = metabolicModel.CNAmodel.objFunc;
                    currentConstraints = metabolicModel.CNAconstraints;
                    currentNumr = metabolicModel.CNAmodel.numr;
                    % CellNetAnalyzer does not allow to define right hand
                    % side; hence we need to do this on the left side
                    numberOfSpecies = length(metabolicModel.CNAmodel.specID);
                    metabolicModel.CNAmodel.reacMax = [metabolicModel.CNAmodel.reacMax; 1000*ones(numberOfSpecies,1)];
                    metabolicModel.CNAmodel.reacMin = [metabolicModel.CNAmodel.reacMin; -1000*ones(numberOfSpecies,1)];
                    numberOfReactions = length(metabolicModel.CNAmodel.objFunc);
                    metabolicModel.CNAmodel.objFunc = [metabolicModel.CNAmodel.objFunc; 0*ones(numberOfSpecies,1)];
                    bconstraints = -metabolicModel.CNAmodel.stoichMat * metabolicModel.lastFlux;
                    metabolicModel.CNAmodel.stoichMat = [metabolicModel.CNAmodel.stoichMat -eye(numberOfSpecies)];
                    metabolicModel.CNAconstraints = [metabolicModel.CNAconstraints; bconstraints];
                    metabolicModel.CNAmodel.numr = metabolicModel.CNAmodel.numr + numberOfSpecies;
                case 2
                    currentB = metabolicModel.COBRAmodel.b;
                    currentUB = metabolicModel.COBRAmodel.ub;
                    currentLB = metabolicModel.COBRAmodel.lb;
                    currentObj = metabolicModel.COBRAmodel.c;
                    metabolicModel.COBRAmodel.b = -metabolicModel.COBRAmodel.S * metabolicModel.lastFlux;
                    numberOfReactions = length(metabolicModel.COBRAmodel.c);
            end
            % define upper and lower limits
            for i = 1:numberOfReactions
                if (~ismember(i, metabolicModel.exchangeReactions.ReacID))
                    fluxDiff = flux(i) - metabolicModel.lastFlux(i);
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAmodel.reacMin(i) = flux(i) - 2 * abs(fluxDiff);
                            metabolicModel.CNAmodel.reacMax(i) = flux(i) + 2 * abs(fluxDiff);
                        case 2
                            metabolicModel.COBRAmodel.lb(i) = flux(i) - 2 * abs(fluxDiff);
                            metabolicModel.COBRAmodel.ub(i) = flux(i) + 2 * abs(fluxDiff);
                    end
                    % maintain sign of delta x during optimization!
                    if flux(i) > metabolicModel.lastFlux(i)
                        switch metabolicModel.FBAsolver
                            case 1
                                metabolicModel.CNAmodel.reacMin(i) = metabolicModel.lastFlux(i);
                                metabolicModel.CNAmodel.objFunc(i) = 1;
                            case 2
                                metabolicModel.COBRAmodel.lb(i) = metabolicModel.lastFlux(i);
                                metabolicModel.COBRAmodel.c(i) = 1;
                        end
                    else
                        if flux(i) < metabolicModel.lastFlux(i)
                            switch metabolicModel.FBAsolver
                                case 1
                                    metabolicModel.CNAmodel.reacMax(i) = metabolicModel.lastFlux(i);
                                    metabolicModel.CNAmodel.objFunc(i) = -1;
                                 case 2
                                    metabolicModel.COBRAmodel.ub(i) = metabolicModel.lastFlux(i);
                                    metabolicModel.COBRAmodel.c(i) = -1;
                            end
                        else
                            % flux was identical, hard to know in which
                            % direction we should allow some slack ... just
                            % ignore this reaction.
                            switch metabolicModel.FBAsolver
                                case 1
                                    if isnan(metabolicModel.CNAconstraints(i))
                                        %metabolicModel.CNAconstraints(i) = flux(i); % fixing does not work unfortunately
                                        metabolicModel.CNAmodel.reacMin(i) = flux(i) - solverPars.fluxTolerance;
                                        metabolicModel.CNAmodel.reacMax(i) = flux(i) + solverPars.fluxTolerance;
                                        metabolicModel.CNAmodel.objFunc(i) = 0;
                                    end
                                case 2
                                    metabolicModel.COBRAmodel.lb(i) = flux(i) - solverPars.fluxTolerance;
                                    metabolicModel.COBRAmodel.ub(i) = flux(i) + solverPars.fluxTolerance;
                                    metabolicModel.COBRAmodel.c(i) = 0;
                             end
                         end
                    end
                    reacIsFixed = 0;
                    switch metabolicModel.FBAsolver
                        case 1
                            reacIsFixed = ~isnan(metabolicModel.CNAconstraints(i));
                        case 2
                            reacIsFixed = currentLB(i) == currentUB(i);
                    end
                    if reacIsFixed
                        switch metabolicModel.FBAsolver
                            case 1
                                metabolicModel.CNAmodel.objFunc(i) = 0;
                                metabolicModel.CNAconstraints(i) = 0;
                                if metabolicModel.CNAmodel.reacMin(i) == metabolicModel.CNAmodel.reacMax(i)
                                    metabolicModel.CNAmodel.reacMin(i) = -1;
                                    metabolicModel.CNAmodel.reacMax(i) = +1;
                                end
                            case 2
                                metabolicModel.COBRAmodel.lb(i) = 0;
                                metabolicModel.COBRAmodel.ub(i) = 0;
                                metabolicModel.COBRAmodel.c(i) = 0;
                        end
                    else
                        switch metabolicModel.FBAsolver
                            case 1
                                myLB = metabolicModel.CNAmodel.reacMin(i);
                                myUB = metabolicModel.CNAmodel.reacMax(i);
                            case 2
                                myLB = metabolicModel.COBRAmodel.lb(i);
                                myUB = metabolicModel.COBRAmodel.ub(i);
                        end
                        if currentLB(i) > myLB
                            switch metabolicModel.FBAsolver
                                case 1
                                    metabolicModel.CNAmodel.reacMin(i) = currentLB(i);
                                case 2
                                    metabolicModel.COBRAmodel.lb(i) = currentLB(i);
                            end
                        end
                        if currentUB(i) < myUB
                            switch metabolicModel.FBAsolver
                                case 1
                                    metabolicModel.CNAmodel.reacMax(i) = currentUB(i);
                                case 2
                                    metabolicModel.COBRAmodel.ub(i) = currentUB(i);
                            end
                        end
                        % convert to diff
                        switch metabolicModel.FBAsolver
                            case 1
                                metabolicModel.CNAmodel.reacMin(i) = metabolicModel.CNAmodel.reacMin(i) - metabolicModel.lastFlux(i);
                                metabolicModel.CNAmodel.reacMax(i) = metabolicModel.CNAmodel.reacMax(i) - metabolicModel.lastFlux(i);
                            case 2
                                metabolicModel.COBRAmodel.lb(i) = metabolicModel.COBRAmodel.lb(i) - metabolicModel.lastFlux(i);
                                metabolicModel.COBRAmodel.ub(i) = metabolicModel.COBRAmodel.ub(i) - metabolicModel.lastFlux(i);
                        end
                    end
                else % fix exchange reactions
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAconstraints(i) = flux(i) - metabolicModel.lastFlux(i);
                            metabolicModel.CNAmodel.objFunc(i) = 0;
                        case 2
                            metabolicModel.COBRAmodel.lb(i) = flux(i) - metabolicModel.lastFlux(i);
                            metabolicModel.COBRAmodel.ub(i) = flux(i) - metabolicModel.lastFlux(i);
                            metabolicModel.COBRAmodel.c(i) = 0; 
                    end
                end
            end
            for k = 1:numberOfReactions
                if currentObj(k) ~= 0
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAconstraints(k) = flux(k) - metabolicModel.lastFlux(k);
                            metabolicModel.CNAmodel.objFunc(k) = 0;
                        case 2
                            metabolicModel.COBRAmodel.lb(k) = flux(k) - metabolicModel.lastFlux(k);
                            metabolicModel.COBRAmodel.ub(k) = flux(k) - metabolicModel.lastFlux(k);
                            metabolicModel.COBRAmodel.c(k) = 0;
                    end
                end
            end
            
            switch metabolicModel.FBAsolver
                case 1
                    [myflux, mysuccess, mystatus] = CNAoptimizeFlux(metabolicModel.CNAmodel, metabolicModel.CNAconstraints);
                case 2
                    solution = optimizeCbModel(metabolicModel.COBRAmodel, 'min', 0, 1);
                    mysuccess = solution.stat == 1;
            end
                        
            succeeding = 1;
            if (~mysuccess)
                if metabolicModel.FBAsolver == 2
		    if solverPars.logLevel > 1
                        fprintf('Minimization to previous flux did not work out. Relaxing restrictions ...')
		    end
                end
                
                tryValue = solverPars.minRelaxValue;
                switch metabolicModel.FBAsolver
                    case 1
                        myUB = metabolicModel.CNAmodel.reacMax;
                        myLB = metabolicModel.CNAmodel.reacMin;
                    case 2
                        myUB = metabolicModel.COBRAmodel.ub;
                        myLB = metabolicModel.COBRAmodel.lb;
                end
                while (~mysuccess && tryValue < solverPars.maxRelaxValue && metabolicModel.FBAsolver == 2)
		    if solverPars.logLevel > 1	
                        fprintf('Trying %d\n', tryValue)
		    end
                    for i = 1:numberOfReactions
                        if myLB(i) ~= myUB(i)
                            if myLB(i) == 0
                                switch metabolicModel.FBAsolver
                                    case 1
                                        metabolicModel.CNAmodel.reacMin(i) = -tryValue;
                                    case 2
                                        metabolicModel.COBRAmodel.lb(i) = -tryValue;
                                end
                            else
                                if myUB(i) == 0
                                    switch metabolicModel.FBAsolver
                                        case 1
                                            metabolicModel.CNAmodel.reacMax(i) = tryValue;
                                       case 2
                                            metabolicModel.COBRAmodel.ub(i) = tryValue;
                                    end
                                end
                            end
                        end
                    end
                    switch metabolicModel.FBAsolver
                        case 1
                            [myflux, mysuccess, mystatus] = CNAoptimizeFlux(metabolicModel.CNAmodel, metabolicModel.CNAconstraints);
                        case 2
                            solution = optimizeCbModel(metabolicModel.COBRAmodel, 'min', 0, 1);
                            mysuccess = solution.stat == 1;
                    end
                    tryValue = tryValue * 10;
                end
                if ~mysuccess
		    if solverPars.logLevel > 1
                        fprintf('Failing, using original flux.')
		    end
                    succeeding = 0;
                else
		    if solverPars.logLevel > 1
                        fprintf('Success')
		    end
                end
            end

            if succeeding == 1
                switch metabolicModel.FBAsolver
                    case 1
                        myflux = validateFlux(myflux, metabolicModel.CNAmodel.reacMin, metabolicModel.CNAmodel.reacMax);
                        flux = myflux(1:length(metabolicModel.lastFlux)) + metabolicModel.lastFlux;
                    case 2
                        solution.x = validateFlux(solution.x, metabolicModel.COBRAmodel.lb, metabolicModel.COBRAmodel.ub);
                        flux = solution.x + metabolicModel.lastFlux;
                end
                
                my_dev = originalDev - sum(abs(flux-metabolicModel.lastFlux));
                if my_dev > 0
		    if solverPars.logLevel > 1	
                        fprintf('devPrevFlux improvement (%f)\n', my_dev)
		    end
                else
		    if solverPars.logLevel > 1	
                        fprintf('devPrevFlux failure (%f)\n', my_dev)
		    end
                end
            end
            % restore original model
            switch metabolicModel.FBAsolver
                case 1
                    metabolicModel.CNAmodel.stoichMat = currentM;
                    metabolicModel.CNAmodel.reacMax = currentUB;
                    metabolicModel.CNAmodel.reacMin = currentLB;
                    metabolicModel.CNAmodel.objFunc = currentObj;
                    metabolicModel.CNAconstraints = currentConstraints;
                    metabolicModel.CNAmodel.numr = currentNumr;
                case 2
                    metabolicModel.COBRAmodel.b = currentB;
                    metabolicModel.COBRAmodel.c = currentObj;
                    metabolicModel.COBRAmodel.lb = currentLB;
                    metabolicModel.COBRAmodel.ub = currentUB;
            end
        end
    end
    
    dbio = flux(metabolicModel.biomassReac) * biomass; % g/l/hour
    
    for i = 1:length(metabolicModel.coupledReactions.ReacID)
        dcompound(i) = flux(metabolicModel.coupledReactions.ReacID(i)) * metabolicModel.coupledReactions.SecretionSense(i) * biomass;
%         if dcompound(i) > 0 % don't feed back production to reactor to test what happens without compound exchange!
%             dcompound(i) = 0;
%         end
    end
    
% recording limiting fluxes
if solverPars.recordLimitingFluxes == 1
    limitingFluxes = getLimitedFluxes(flux, metabolicModel, solverPars)';
    fixedFluxes = find(limitingFluxes == 100);
    limitedFluxes = find(abs(limitingFluxes) == 1);
    resultLimitingFluxes = [fixedFluxes; limitedFluxes];
    resultLimitingFluxes(:, 2) = [flux(fixedFluxes); flux(limitedFluxes)];
    resultLimitingFluxes(:, 3) = [ones(length(fixedFluxes), 1); zeros(length(limitedFluxes), 1)];
else
    resultLimitingFluxes = '';
end
end

