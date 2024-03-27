function [ f ] = plotTrajectoryDynamics( trajectory, reactor )
%PLOTRUN Plot Dynamics of a simulated trajectory, should work with full and simplified trajectory files (*_dynamics)
%   Expected arguments: filename or trajectory + reactor

if nargin == 1
    name = trajectory;
    load(trajectory);
    h=figure('visible', 'off');
else
    h=figure;
end

subplot(2,2,1)
plot(trajectory.time, trajectory.biomass)
legend(trajectory.modelNames, 'Interpreter', 'none')
title('Biomass')
xlabel('Time (hours)')
ylabel('Concentration (gDW/L)')
subplot(2,2,2)
plot(trajectory.time, trajectory.compounds)
legend(trajectory.compoundNames, 'Interpreter', 'none')
title('Compounds')
xlabel('Time (hours)')
ylabel('Concentration (mmol/L)')
subplot(2,2,3)
plot(trajectory.time(2:end), trajectory.mu)
%legend(trajectory.modelNames, 'Interpreter', 'none')
title('Specific growth rate µ')
xlabel('Time (hours)')
ylabel('µ (1/h)')
if reactor.flowRate ~= 0
    hold on
    yline(reactor.flowRate, '--')
    hold off
end

if nargin == 1
    saveas(h, strcat(name, '.png'))
end

end

