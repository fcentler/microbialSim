function [ maxFlux ] = getMaxUptakeFlux(metabolicModel, i, metabConcentration)
%GETMAXUPTAKEFLUX Compute current uptake limits based on Monod kinetics.
%   Detailed explanation goes here

        if metabolicModel.coupledReactions{i, 'inhibited'} == 1
            maxFlux = 0.0;
        else
            maxFlux = metabolicModel.coupledReactions{i, 'vmax'};
            maxFlux = maxFlux * -1 * metabolicModel.coupledReactions{i, 'SecretionSense'};
            maxFlux = maxFlux * metabConcentration;
            if maxFlux ~= 0   % avoid division by zero if conc is zero and ks is zero
                maxFlux = maxFlux / (metabolicModel.coupledReactions{i, 'ks'} + metabConcentration);
            end
        end

end

