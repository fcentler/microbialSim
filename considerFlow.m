function [ dflow ] = considerFlow(inflowConc, current, reactor)
%CONSIDERFLOW Compute the change of reactor concentrations due to
%in-/outflow
%   The current flow rate is computed based on set flow rate and reactor
%   volume, and the rate of change of reactor concentrations due to flow
%   conditions is then computed and returned. The function is called for
%   biomass concentrations and compound concentrations separately.

        flowrate = reactor.flowRate / reactor.volume;
        dflow = flowrate * (inflowConc - current);

end

