function [vmax, flux, exchangeFluxes] = estimateVmax(metabolicModel, reactor, solverPars, coupledReactionID, targetMu)
%ESTIMATEVMAX Given a growth rate, this function computes which uptake flux
%of the growth limiting substrate gives rise to this growth rate, assuming
%that all other uptake fluxes are set to reasonable / problem-specific
%values. If 
%   Detailed explanation goes here

orig.vmax = metabolicModel.coupledReactions.vmax;
orig.ks = metabolicModel.coupledReactions.ks;

numberOfReactions = length(metabolicModel.coupledReactions.vmax);

%% report on growth rate under conditions specified in model
metabolicModel.coupledReactions.ks = zeros(numberOfReactions, 1);
[ dbio, dcompound, flux, success ] = getDy(metabolicModel, ones(1,length(reactor.compounds)), -1, solverPars);
if success == true
    fprintf("Growth rate under original model uptake conditions = %d 1/h.\n", flux(metabolicModel.biomassReac))
else
    fprintf("Original conditions do not allow for growth.\n");
end

%% check if selected compound is required for growth; if not: adjust all
% flux boundaries

growthLimiting = true;

metabolicModel.coupledReactions.vmax(coupledReactionID) = 0;
[ dbio, dcompound, flux, success ] = getDy(metabolicModel, ones(1,length(reactor.compounds)), -1, solverPars);
if success == true
    growthLimiting = false;
    fprintf("Selected compound not required for growth, adjusting all flux limits together.\n");
    fprintf("Check if uncoupled (and coupled) exchange reactions still contain\ngrowth-relevant fluxes (such as c-sources).\n")
else
    fprintf("Selected compound is required for growth.\n");
end

if growthLimiting == true
    stepSize = 10;
    subStepSize = 1;
else
    stepSize = 0.1;
    subStepSize = 0.01;
    switch metabolicModel.FBAsolver
        case 1
            orig.lb = metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID);
            orig.ub = metabolicModel.CNAmodel.reacMax(metabolicModel.exchangeReactions.ReacID);
        case 2
            orig.lb = metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID);
            orig.ub = metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID);
        otherwise
            error('In estimateVmax: unknown FBA solver')
    end
end


%% increase vmax until growth rate is larger then target growth rate
currentMu = 0;
if growthLimiting == true
    vmax = 0;
else
    vmax = stepSize/1000.0 - stepSize; % factor applied to all uptake fluxes
end
while currentMu < targetMu && vmax < 1000.0
    vmax = vmax + stepSize;
    if growthLimiting == true
        metabolicModel.coupledReactions.vmax(coupledReactionID) = vmax;
    else
        for i = 1:length(metabolicModel.exchangeReactions.ReacID)
            if metabolicModel.exchangeReactions.SecretionSense(i) == 1
                switch metabolicModel.FBAsolver
                    case 1
                        metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * vmax;
                    case 2
                        metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * vmax;
                    otherwise
                        error('In estimateVmax: unknown FBA solver')
                end
            else
                if metabolicModel.exchangeReactions.SecretionSense(i) == -1
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAmodel.reacMax(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * vmax;
                        case 2
                            metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * vmax;
                        otherwise
                            error('In estimateVmax: unknown FBA solver')
                    end
                end
            end
        end
        metabolicModel.coupledReactions.vmax = orig.vmax * vmax;
    end
    [ dbio, dcompound, flux, success ] = getDy(metabolicModel, ones(1,length(reactor.compounds)), -1, solverPars);
    if success == true
        currentMu = flux(metabolicModel.biomassReac);
    else
        currentMu = 0;
    end
end
if currentMu < targetMu
    fprintf('Target µ (%d) not reachable, just reached %d 1/h\n', targetMu, currentMu);
    vmax = -1;
    return
end

%% establish lower bound with growth
lowBound = vmax - stepSize - subStepSize;
solvable = false;
while solvable == false
    lowBound = lowBound + subStepSize;
    if growthLimiting == true
        metabolicModel.coupledReactions.vmax(coupledReactionID) = lowBound;
    else
        for i = 1:length(metabolicModel.exchangeReactions.ReacID)
            if metabolicModel.exchangeReactions.SecretionSense(i) == 1
                switch metabolicModel.FBAsolver
                    case 1
                        metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * lowBound;
                    case 2
                        metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * lowBound;
                    otherwise
                        error('In estimateVmax: unknown FBA solver')
                end
            else
                if metabolicModel.exchangeReactions.SecretionSense(i) == -1
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAmodel.reacMax(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * lowBound;
                        case 2
                            metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * lowBound;
                        otherwise
                            error('In estimateVmax: unknown FBA solver')
                    end
                end
            end
        end
        metabolicModel.coupledReactions.vmax = orig.vmax * lowBound;
    end
    [ dbio, dcompound, flux, success ] = getDy(metabolicModel, ones(1,length(reactor.compounds)), -1, solverPars);
    if success == true
        solvable = true;
    end
end

if flux(metabolicModel.biomassReac) > targetMu
    fprintf('Failed to find proper lower bound. Chose smaller subStepSize\n');
    vmax = -1;
    return
end

%% now approach target mu
highBound = vmax;
while abs(currentMu - targetMu) > solverPars.minimalGrowth
    vmax = (lowBound + highBound) / 2.0;
    if growthLimiting == true
        metabolicModel.coupledReactions.vmax(coupledReactionID) = vmax;
    else
        for i = 1:length(metabolicModel.exchangeReactions.ReacID)
            if metabolicModel.exchangeReactions.SecretionSense(i) == 1
                switch metabolicModel.FBAsolver
                    case 1
                        metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * vmax;
                    case 2
                        metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID(i)) = orig.lb(i) * vmax;
                    otherwise
                        error('In estimateVmax: unknown FBA solver')
                end
            else
                if metabolicModel.exchangeReactions.SecretionSense(i) == -1
                    switch metabolicModel.FBAsolver
                        case 1
                            metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * vmax;
                        case 2
                            metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID(i)) = orig.ub(i) * vmax;
                        otherwise
                            error('In estimateVmax: unknown FBA solver')
                    end
                end
            end
        end
        metabolicModel.coupledReactions.vmax = orig.vmax * vmax;
    end
    [ dbio, dcompound, flux, success ] = getDy(metabolicModel, ones(1,length(reactor.compounds)), -1, solverPars);
    if success == true
        currentMu = flux(metabolicModel.biomassReac);
        if currentMu < targetMu
            lowBound = vmax;
        else
            highBound = vmax;
        end
    else
        fprintf('Fail in approaching target µ!\n');
        vmax = -1;
        return
    end
end
[exchangeFluxes, limitedFluxesIDs] = getExchangeFluxes(metabolicModel, flux, solverPars, ones(1,length(reactor.compounds)));

if growthLimiting == true
    fprintf('For getting µ = %d 1/h, choose vmax = %d mmol/gDW/h \nfor the %i th coupled reaction.\n', targetMu, vmax, coupledReactionID);
else
    fprintf('Applying a factor of %d to all uptake limits results in µ = %d 1/h.\n', vmax, targetMu);
    fprintf('vmax values for coupled reactions:\n');
    for i = 1:numberOfReactions
        fprintf('%i: %d\n', i, metabolicModel.coupledReactions.vmax(i));
    end
end

fprintf('\nFluxes at limit:\n');
for i = 1:length(limitedFluxesIDs)
    switch metabolicModel.FBAsolver
        case 1
            fprintf("ID %i (%s): %d ", limitedFluxesIDs(i), metabolicModel.CNAmodel.reacID(limitedFluxesIDs(i),:), flux(limitedFluxesIDs(i)));
        case 2
            fprintf("ID %i (%s): %d ", limitedFluxesIDs(i), metabolicModel.COBRAmodel.rxnNames{limitedFluxesIDs(i)}, flux(limitedFluxesIDs(i)));
        otherwise
            error('In estimateVmax: unknown FBA solver')
    end
    if ismember(limitedFluxesIDs(i), metabolicModel.exchangeReactions.ReacID)
        fprintf("(Exchange reaction, ");
        if ismember(limitedFluxesIDs(i), metabolicModel.coupledReactions.ReacID)
            fprintf("coupled)");
        else
            fprintf("uncoupled)");
        end
    end
    fprintf('\n');
end
fprintf('\n');

%likely not really necessary; just for peace of mind
metabolicModel.coupledReactions.vmax = orig.vmax;
metabolicModel.coupledReactions.ks = orig.ks;
if growthLimiting == false
    switch metabolicModel.FBAsolver
        case 1
            metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID) = orig.lb;
            metabolicModel.CNAmodel.reacMax(metabolicModel.exchangeReactions.ReacID) = orig.ub;
        case 2
            metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID) = orig.lb;
            metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID) = orig.ub;
        otherwise
            error('In estimateVmax: unknown FBA solver')
    end
end
end

