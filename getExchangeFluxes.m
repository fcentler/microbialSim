function [exchangeFluxTable, limitedFluxesIDs] = getExchangeFluxes(metabolicModel, flux, solverPars, varargin)
%UNTITLED Report on the exchange fluxes given a flux distribution, mark
%coupled reactions, mark reactions at limit
%   Detailed explanation goes here

if nargin > 4
    error('Too many arguments in getExchangeFluxes()!');
end

exchangeFluxTable = metabolicModel.exchangeReactions;

% mapping back reversible reactions
switch metabolicModel.FBAsolver
    case 1
        numberOfReactions = length(metabolicModel.CNAmodel.objFunc);
    case 2
        numberOfReactions = length(metabolicModel.COBRAmodel.c);
    otherwise
        error('In getExchangeFluxes: unknown FBA solver')
end
if isfield(metabolicModel, 'reversibleReacIDs')
    numberOfReactions = numberOfReactions - length(metabolicModel.reversibleReacIDs);
    for i = 1:length(metabolicModel.reversibleReacIDs)
        flux(metabolicModel.reversibleReacIDs(i)) = flux(metabolicModel.reversibleReacIDs(i)) + flux(numberOfReactions + i);
    end
end
fluxes = flux(1:numberOfReactions);

FluxValues = getFluxes(fluxes, metabolicModel.exchangeReactions.ReacID);

if nargin == 4
    compounds = varargin{1};
    [limitedFluxes, evalLB, evalUB] = getLimitedFluxes(fluxes, metabolicModel, solverPars, compounds);
else
    [limitedFluxes, evalLB, evalUB] = getLimitedFluxes(fluxes, metabolicModel, solverPars);
end
limitedFluxesIDs = find(limitedFluxes);

ProductionConsumption = FluxValues .* metabolicModel.exchangeReactions.SecretionSense;
AtFluxLimit = limitedFluxes(metabolicModel.exchangeReactions.ReacID)';
coupledReac = zeros(1, length(fluxes))';
coupledReac(metabolicModel.coupledReactions.ReacID) = 1;
IsCoupled = coupledReac(metabolicModel.exchangeReactions.ReacID);

exchangeFluxTable = addvars(exchangeFluxTable, FluxValues);
exchangeFluxTable = addvars(exchangeFluxTable, ProductionConsumption);
exchangeFluxTable = addvars(exchangeFluxTable, AtFluxLimit);
exchangeFluxTable = addvars(exchangeFluxTable, IsCoupled);

exchangeFluxTable = addvars(exchangeFluxTable, evalLB(metabolicModel.exchangeReactions.ReacID), 'NewVariableNames', 'evalLB');
exchangeFluxTable = addvars(exchangeFluxTable, evalUB(metabolicModel.exchangeReactions.ReacID), 'NewVariableNames', 'evalUB');

switch metabolicModel.FBAsolver
    case 1
        lb = metabolicModel.CNAmodel.reacMin(metabolicModel.exchangeReactions.ReacID);
        ub = metabolicModel.CNAmodel.reacMax(metabolicModel.exchangeReactions.ReacID);
    case 2
        lb = metabolicModel.COBRAmodel.lb(metabolicModel.exchangeReactions.ReacID);
        ub = metabolicModel.COBRAmodel.ub(metabolicModel.exchangeReactions.ReacID);
    otherwise
        error('In getExchangeFluxes: unknown FBA solver')
end

exchangeFluxTable = addvars(exchangeFluxTable, lb);
exchangeFluxTable = addvars(exchangeFluxTable, ub);

end

