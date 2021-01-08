function [ model ] = prepareFBAmodel_iMhu428Model(filename, FBAsolver)
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

switch FBAsolver
    case 1
        error('Model not configured to work with CNA!');
        CNAmodel = CNAloadNetwork({filename, 1}, 1, 1);
        % non growth-dependend ATP demand (value will be overwritten in
        % parametrize*Model().
        CNAconstraints = nan(length(CNAmodel.reacID), 1);
        CNAconstraints(ngamReac) = 2.0; % dummy value
        % l-cys not required for growth
        CNAmodel.reacMin(221) = 0;
        % gly-ald not required for growth
        CNAmodel.reacMin(231) = 0;
    case 2
        COBRAmodel = readCbModel(filename);
        % allow cell to export all metabolites coupled by exchange
        % reactions
        for i = 234:274
            COBRAmodel.ub(i) = 1000;
        end
         
        % Apply settings according to SI S2 (Excel-File), for co-culture
        % growth of iMhu428 (Hamilton et al., 2015)
        
        % allow for h2s, h2, formate, and phosphate input
        COBRAmodel = changeRxnBounds(COBRAmodel, 'H2ST', 1000, 'u');
        COBRAmodel = changeRxnBounds(COBRAmodel, 'H2TD', -1000, 'l');
        COBRAmodel = changeRxnBounds(COBRAmodel, 'FORT', -1000, 'l');
        COBRAmodel = changeRxnBounds(COBRAmodel, 'PIT', 1000, 'u');
        
        % restrict these reactions to forward only
        for i = {'5HBCR', 'ALDD20X', 'GK1', 'NAKT_1', 'PGM', 'URIDK2R', 'ADCPS1', 'ASNS1', 'GLNS', 'NCCT', 'PIABC', 'XPPT', 'ADK1', 'ATPS1', 'GLUDC', 'NIT_N1P4', 'PPK2', 'ZNABC2', 'ADK3', 'BTNABC', 'HKT', 'NNAT', 'PPKR', 'ADK4', 'CD2ABC1', 'MAN1PT', 'NTD11', 'PROABC', 'ADSS', 'CF3SA', 'MHPGLUT', 'NTD4', 'SULABC', 'AIRC2', 'DCTPD', 'NAABC', 'NTD9', 'TMDK1'}
            COBRAmodel = changeRxnBounds(COBRAmodel,i, 0, 'l');
            COBRAmodel = changeRxnBounds(COBRAmodel,i, 1000, 'u');
        end
        
        % restrict these reactions to backward only
        for i = {'MDH', 'MDHY', 'PPDK'}
            COBRAmodel = changeRxnBounds(COBRAmodel,i, 0, 'u');
            COBRAmodel = changeRxnBounds(COBRAmodel,i, -1000, 'l');
        end
        
        % allow uptake for these compounds
        for i = {'EX_AC_E', 'EX_CO2_E', 'EX_NH4_E', 'EX_H2S_E', 'EX_H2_E', 'EX_FOR_E','EX_NI2_E', 'EX_PI_E', 'EX_H2O_E', 'EX_H_E', 'EX_COBALT2_E'}%, 'EX_HCO3_E'}
            COBRAmodel = changeRxnBounds(COBRAmodel,i, -1000, 'l');
        end
        % fix direction of FMFTSPFT
            COBRAmodel = changeRxnBounds(COBRAmodel, 'FMFTSPFT', 0, 'l');
            COBRAmodel = changeRxnBounds(COBRAmodel, 'FMFTSPFT', 1000, 'u');
    otherwise
        error('FBA code wrong');
end

% index of NGAM reaction
ngamReac=108; 
% biomass reaction
biomassReac = 549;
% exchange reactions; name all, not only those which will be coupled to
% reactor compounds
IDs = [234:274 486 571 biomassReac];
switch FBAsolver
    case 1
        ReactionIDName = CNAmodel.reacID(IDs, :);
    case 2
        ReactionIDName = COBRAmodel.rxns(IDs, :);
end

% Define the reaction sense. A 1 indicates that a positive flux of the
% reaction indicates metabolite excretion, a -1, that a negative flux
% indicates excretion.

Sense = ones(length(IDs), 1);
Sense(42) = -1; % nitrogen diffusion, reaction 486
Sense(43) = -1; % inorganig phosphate diffusion, reaction 571

exchangeReactions = table(IDs', cellstr(ReactionIDName), Sense, 'VariableNames',{'ReacID' 'ReacName' 'SecretionSense'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% Define reactions that will be coupled to reactor compounds (substrate
% and/or products). All these reactions must have been indicated as
% exchange reactions above.
IDs = [237 245 255 249 244]; % acetate, co2, h2, Formate, methane
% Define uptake kinetics. These dummy values will be overwritten in
% function parametrize*model().
vmax = [100 100 100 100 0];
ks = [0 0 0 0 0];
% Never allow uptake through these reactions (set value to 1).
inhibitedCoupledReactions = [0 0 0 0 0]; 

switch FBAsolver
    case 1
        ReacNames = CNAmodel.reacID(IDs, :);
    case 2
        ReacNames = COBRAmodel.rxns(IDs, :);
end

coupledReactions = table(IDs', cellstr(ReacNames), exchangeReactions.SecretionSense(strtrim(cellstr(num2str(IDs'))')), vmax', ks', inhibitedCoupledReactions' ,'VariableNames',{'ReacID' 'ReacName' 'SecretionSense' 'vmax' 'ks' 'inhibited'}, 'RowNames', strtrim(cellstr(num2str(IDs'))'));

% Couple coupled reactions to reactor compounds. The position in the array
% is the same as in "IDs" of the coupled reactions, the number indicates
% the position in the order of the reactor compounds. For example, a 5 in
% the second position indicates that the second coupled reaction describes
% uptake and/or secretion of the compound at the 5th position in the
% reactor compound array.
reactorCompoundIDs = [2 5 6 7 8];

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
model.modelName = 'iMhu428';
