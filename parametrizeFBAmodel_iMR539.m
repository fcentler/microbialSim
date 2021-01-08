function [ model ] = parametrizeFBAmodel_iSfu648Model(model)
%PARAMETRIZE Setting the parameters (NGAM and uptake kinetics)for the model
%   Detailed explanation goes here

%% NGAM

NGAMValue = 5.1176; % mmolATP/gDW/h

switch model.FBAsolver
    case 1
        model.CNAconstraints(model.ngamReac) = NGAMValue;
    case 2
        model.COBRAmodel.lb(model.ngamReac) = NGAMValue;
        model.COBRAmodel.ub(model.ngamReac) = NGAMValue;
    otherwise
       display('parametrizeModel(): unkown type selected')
end

%% UPTAKE KINETICS
%CO2, H2, Acetate, Methane
model.coupledReactions.vmax = [1000 189.3 0]';
model.coupledReactions.ks = [0 4.375e-4 0]';

