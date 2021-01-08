function [ model ] = convertToOnlyIrrevReac( model )
%CONVERTTOONLYIRREVREAC Modifying a model such that it only contains
%irreversible reactions
%   For reversible reactions, these reactions are converted to only
%   proceeding in the forward direction. The backward reaction is appended
%   to the list of reactions. Exchange reactions are not modified.

% Collecting the ids of originally reversible reactions
model.reversibleReacIDs = [];

switch model.FBAsolver
    case 1
        numberOfReactions = length(model.CNAmodel.objFunc);
    case 2
        numberOfReactions = length(model.COBRAmodel.c);
    otherwise
        error('In convertToOnlyIrrecReac: unknown FBA solver')
end

for i = 1:numberOfReactions
    % don't modify exchange reactions
    if (~ismember(i, model.exchangeReactions.ReacID))
    % check whether reversible
    switch model.FBAsolver
        case 1
            isReversible = sign(model.CNAmodel.reacMin(i)) * sign(model.CNAmodel.reacMax(i)) < 0;
        case 2
            isReversible = sign(model.COBRAmodel.lb(i)) * sign(model.COBRAmodel.ub(i)) < 0;
        otherwise
            error('In convertToOnlyIrrecReac: unknown FBA solver')
    end
    if isReversible
        model.reversibleReacIDs(length(model.reversibleReacIDs) + 1) = i;
        % copy reaction in stoichiometric matrix and bounds and
        % optimization coefficient
        switch model.FBAsolver
            case 1
                model.CNAmodel.stoichMat(:,size(model.CNAmodel.stoichMat, 2) + 1) = model.CNAmodel.stoichMat(:,i);
                model.CNAmodel.reacMin(end + 1) = model.CNAmodel.reacMin(i);
                model.CNAmodel.reacMax(end + 1) = model.CNAmodel.reacMax(i);
                model.CNAmodel.objFunc(end + 1) = model.CNAmodel.objFunc(i);
                model.CNAconstraints(end + 1) = model.CNAconstraints(i);
                model.CNAmodel.numr = model.CNAmodel.numr + 1;
                % fix directions & optimization coefficient
                model.CNAmodel.reacMin(i) = 0;
                model.CNAmodel.reacMax(length(model.CNAmodel.reacMin)) = 0;
                if (model.CNAmodel.objFunc(length(model.CNAmodel.objFunc)) ~=0)
                    error('Case not considered yet!')
                end
                if ~isnan(model.CNAconstraints(i))
                    if model.CNAconstraints(i) >= 0
                        model.CNAconstraints(length(model.CNAconstraints)) = 0;
                    else
                        model.CNAconstraints(i) = 0;
                    end
                end
                newName = strcat(model.CNAmodel.reacID(i,:),'-reverse');
                newName = strcat(newName, repmat('_',1,length(model.CNAmodel.reacID(1,:))-length(newName)));
                model.CNAmodel.reacID(end + 1,:) = newName;
            case 2
                model.COBRAmodel.S(:,size(model.COBRAmodel.S, 2) + 1) = model.COBRAmodel.S(:,i);
                model.COBRAmodel.lb(end + 1) = model.COBRAmodel.lb(i);
                model.COBRAmodel.ub(end + 1) = model.COBRAmodel.ub(i);
                model.COBRAmodel.c(end + 1) = model.COBRAmodel.c(i);
                % fix directions & optimization coefficient
                model.COBRAmodel.lb(i) = 0;
                model.COBRAmodel.ub(length(model.COBRAmodel.lb)) = 0;
                if (model.COBRAmodel.c(length(model.COBRAmodel.c)) ~=0)
                    error('Case not considered yet!')
                end
                model.COBRAmodel.rxns(end + 1) = strcat(model.COBRAmodel.rxns(i),'-reverse');
                model.COBRAmodel.rxnNames(end + 1) = strcat(model.COBRAmodel.rxnNames(i),'-reverse');

            otherwise
                error('In convertToOnlyIrrecReac: unknown FBA solver')
        end
    end
    end
end
end

