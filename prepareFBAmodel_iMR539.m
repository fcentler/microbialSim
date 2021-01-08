function [ model ] = prepareFBAmodel_iSfu648Model(filename, FBAsolver)
%prepareFBAModel
%   Load FBA model and, if necessary, adapt flux boundaries. Functions for
%   doing this vary between CellNetAnalyzer and COBRA Toolbox, hence the
%   switch over FBAsolver. Additionally, reaction indices for the NGAM and
%   biomass reaction must be provided, and exchange reactions must be
%   defined that link compounds in the reactor to cellular uptake and
%   secretion. It must be indicated what a positive/negative flux means:
%   uptake or secretion. Values provided for NGAM & uptake kinetics are
%   only placeholders here and will be overwritten by function 
%   parametrize*Model().

ngamReac=366; % stoichiometry of phosphor is unconventional
switch FBAsolver
    case 1
    	error("not implemented");
    case 2
        COBRAmodel = readCbModel(filename);
    otherwise
        error('FBA code wrong');
end

% biomass reaction
biomassReac = 475;

% exchange reactions
IDs = [biomassReac 476:515 531 540:541 544 547 598 614 616 640:642 652:654 674:675];
switch FBAsolver
    case 1
        ReactionIDName = CNAmodel.reacID(IDs, :);
    case 2
        %ReactionIDName = COBRAmodel.rxns(IDs, :);
        ReactionIDName = COBRAmodel.rxnNames(IDs, :);
end

%prohibit acetate uptake
COBRAmodel.lb(531) = 0;

% reaction sense: 1 means: positive flux is metabolite excretion, -1 means: negative flux is excretion
Sense = ones(length(IDs), 1);

exchangeReactions = table(IDs', cellstr(ReactionIDName), Sense, 'VariableNames',{'ReacID' 'ReacName' 'SecretionSense'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% SUBSTRATES / PRODUCTS to couple to medium
IDs = [498 515 514]; % 
% DUMMY VALUES, UPDATED in PARAMETRIZATION FUNCTION!
vmax = [10 10 0];
ks = [0.1 0.1 0];
inhibitedCoupledReactions = [0 0 0]; % never allow uptake through these reactions
switch FBAsolver
    case 1
        ReacNames = CNAmodel.reacID(IDs, :);
    case 2
        %ReacNames = COBRAmodel.rxns(IDs, :);
        ReacNames = COBRAmodel.rxnNames(IDs, :);
end
coupledReactions = table(IDs', cellstr(ReacNames), exchangeReactions.SecretionSense(strtrim(cellstr(num2str(IDs'))')), vmax', ks', inhibitedCoupledReactions' ,'VariableNames',{'ReacID' 'ReacName' 'SecretionSense' 'vmax' 'ks' 'inhibited'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% couple model fluxes to reactor compounds
reactorCompoundIDs = [1 2 3]; % Position: as in coupled reactions, number: refers to position in reactor compounds

% prepare return values
switch FBAsolver
    case 1
        model.CNAmodel = CNAmodel;
        model.CNAconstraints = CNAconstraints;
    case 2
        model.COBRAmodel = COBRAmodel;
end
model.biomassReac = biomassReac;
model.ngamReac = ngamReac;
model.exchangeReactions = exchangeReactions;
model.coupledReactions = coupledReactions;
model.reactorCompoundIDs = reactorCompoundIDs;
model.FBAsolver = FBAsolver;
model.modelName = 'iMR539';
