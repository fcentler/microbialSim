function [ flux ] = validateFlux( flux_in, currentLB, currentUB )
%VALIDATEFLUX Summary of this function goes here
%   Detailed explanation goes here
    flux = flux_in;
    % sanity check, sometimes fluxes are slightly off
    for k = 1:length(flux)
        if flux(k) < currentLB(k)
    %                         disp(k);
    %                         disp(flux(k));
            flux(k) = currentLB(k);
    %                         disp('set to');
    %                         disp(flux(k));
        else
            if flux(k) > currentUB(k)
    %                             disp(k);
    %                             disp(flux(k));
                flux(k) = currentUB(k);
    %                             disp('set to');
    %                             disp(flux(k));
            end                        
        end
    end

end

