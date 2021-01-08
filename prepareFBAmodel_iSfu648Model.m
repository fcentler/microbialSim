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

ngamReac=1;
switch FBAsolver
    case 1
CNAmodel = CNAloadNetwork({filename, 1}, 1, 1);

%non growth-dependend ATP demand
CNAconstraints = nan(length(CNAmodel.reacID), 1);
%will be changed in parametrizeGonnerman function
CNAconstraints(ngamReac) = 2.0;

%l-cys not required for growth
CNAmodel.reacMin(221) = 0;
%gly-ald not required for growth
CNAmodel.reacMin(231) = 0;

    case 2
        COBRAmodel = readCbModel(filename);
        % Apply settings according to Reed Excel-File, S2, for Co-culture
        % Growth of Fumaroxidans
        
        % cytFDH + cytH2ase = SDH implemented via additional
        % pseudo-metabolite
        
        for i = [29,28,24,25,812]
            COBRAmodel.lb(i) = 0;
            COBRAmodel.ub(i) = 0;
        end
        
        %Forward only
        for i = {'FRD', 'SDH', 'R00199', 'R00229', 'R00253', 'R00315', 'R00509', 'R00512', 'R00586', 'R00669', 'R00920', 'R01058', 'R01197', 'R01353', 'R01359', 'R01512', 'R01777', 'R03005', 'R03316', 'R03596', 'R08060'}
            COBRAmodel = changeRxnBounds(COBRAmodel,i, 0, 'l');
        end
        
        %Backward only
        for i = {'R00206', 'R00214', 'R00351', 'R01724', 'R02569', 'R03313'}
            COBRAmodel = changeRxnBounds(COBRAmodel,i, 0, 'u');
        end
        
        
    otherwise
        error('FBA code wrong');
end


% biomass reaction
biomassReac = 2;

% exchange reactions
IDs = [7:21 biomassReac];
switch FBAsolver
    case 1
        ReactionIDName = CNAmodel.reacID(IDs, :);
    case 2
        ReactionIDName = COBRAmodel.rxns(IDs, :);
end

% reaction sense: 1 means: positive flux is metabolite excretion, -1 means: negative flux is excretion
Sense = ones(length(IDs), 1);

exchangeReactions = table(IDs', cellstr(ReactionIDName), Sense, 'VariableNames',{'ReacID' 'ReacName' 'SecretionSense'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% SUBSTRATES / PRODUCTS to couple to medium
IDs = [17 11 16 12 9 19 13]; % prop, ac, fum, succ, CO2, H2, for, these all must be included in the exchange reactions defined above!
% DUMMY VALUES, UPDATED in PARAMETRIZATION FUNCTION!
vmax = [100  0  100 0 0 0 0];
ks = [0 0 0 0 0 0 0];
inhibitedCoupledReactions = [0 0 0 0 0 0 0]; % never allow uptake through these reactions
switch FBAsolver
    case 1
        ReacNames = CNAmodel.reacID(IDs, :);
    case 2
        ReacNames = COBRAmodel.rxns(IDs, :);
end
coupledReactions = table(IDs', cellstr(ReacNames), exchangeReactions.SecretionSense(strtrim(cellstr(num2str(IDs'))')), vmax', ks', inhibitedCoupledReactions' ,'VariableNames',{'ReacID' 'ReacName' 'SecretionSense' 'vmax' 'ks' 'inhibited'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% couple model fluxes to reactor compounds
reactorCompoundIDs = [1 2 3 4 5 6 7]; % Position: as in coupled reactions, number: refers to position in reactor compounds

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
model.modelName = 'iSfu648';
