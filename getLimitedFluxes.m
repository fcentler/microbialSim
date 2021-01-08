function [ result, lb, ub ] = getLimitedFluxes( flux, metabolicModel, solverPars, varargin )
%GETLIMITEDFLUXES Identify reactions, for which the current flux is identical
%to a set flux boundary
%   Fluxes at lower limit are indicated by -1, at upper limit by +1, fixed
%   reactions are marked as 100.

if nargin > 4
    error('Too many arguments in getLimitedFluxes()!');
end

switch metabolicModel.FBAsolver
    case 1
        lb = metabolicModel.CNAmodel.reacMin;
        ub = metabolicModel.CNAmodel.reacMax;
    case 2
        lb = metabolicModel.COBRAmodel.lb;
        ub = metabolicModel.COBRAmodel.ub;
    otherwise
         error('getLimitedFluxes: Unkown model type, aborting')
end

% map back to reversible reactions

numberOfReactions = length(lb);

if isfield(metabolicModel, 'reversibleReacIDs')
    numberOfReactions = numberOfReactions - length(metabolicModel.reversibleReacIDs);
    for i = 1:length(metabolicModel.reversibleReacIDs)
        lb(metabolicModel.reversibleReacIDs(i)) = lb(numberOfReactions + i);
    end
end

% get current uptake limits
if nargin == 4
    compounds = varargin{1};
    for i = 1:length(metabolicModel.coupledReactions.ReacID)
        limit = getMaxUptakeFlux(metabolicModel, i, compounds(metabolicModel.reactorCompoundIDs(i)));
        if metabolicModel.coupledReactions{i, 'SecretionSense'} > 0
            lb(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
        else
            ub(metabolicModel.coupledReactions{i, 'ReacID'}) = limit;
        end
    end
end 
% for fixed CNA reactions

if metabolicModel.FBAsolver == 1
    for i = 1:numberOfReactions
        if ~isnan(metabolicModel.CNAconstraints(i))
            lb = metabolicModel.CNAconstraints(i);
            ub = metablicModel.CNAconstraints(i);
        end
    end
end

% find non-zero reactions at the upper or lower limit
result = [];
for i = 1:numberOfReactions
    if abs(flux(i)) > solverPars.fluxTolerance % only check non-zero fluxes
        if abs(flux(i) - lb(i)) < solverPars.fluxTolerance
            result(i) = -1;
        else
            if abs(flux(i) - ub(i)) < solverPars.fluxTolerance
                result(i) = 1;
            else
                result(i) = 0;
            end
        end
        if result(i) == -1 && lb(i) == ub(i)
            result(i) = 100;    % indicating a fixed reaction
        end
    else
        result(i) = 0;
    end
end
end

