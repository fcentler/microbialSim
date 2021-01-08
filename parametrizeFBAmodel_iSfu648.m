function [ model ] = parametrizeFBAmodel_iSfu648Model(model)
%PARAMETRIZE Setting the parameters (NGAM and uptake kinetics)for the model
%   Detailed explanation goes here

%% NGAM

NGAMValue = 3.36/24.0; % mmolATP/gDW/h

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
%prop, ac, fum, CO2, H2, for
model.coupledReactions.vmax = [1.1738 0 0 0 0 0 0]';
model.coupledReactions.ks = [2.7 0 0 0 0 0 0]';

