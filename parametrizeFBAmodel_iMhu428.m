function [ model ] = parametrizeFBAmodel_iMhu428(model)
%PARAMETRIZE Setting the parameters (NGAM and uptake kinetics)for the model
%   Detailed explanation goes here

%% NGAM

NGAMValue = 0.6/24.0; % mmolATP/gDW/h

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
% acetate, co2, h2, Formate, methane
model.coupledReactions.vmax = [0 100.0 27.6 0 0]';
model.coupledReactions.ks = [0 0 0.006 0 0]';

