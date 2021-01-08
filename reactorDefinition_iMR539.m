function [ reactor ] = reactorDefinition_iMR539( input_args )
%reactorDefinition Summary of this function goes here
%   Detailed explanation goes here

% define the reactor
reactor.volume = 1.0;	% liter
reactor.flowRate = 0*0.1;	% liter/hour
reactor.compounds = ...
{
'CO2', ...
'H2', ...
'Methane' ...
};

reactor.compoundsInit = ... %mM
[
1.0 ...
1e-2 ...
0.0 ...
];

reactor.biomassInit = ... % gDW / l
[
1e-4 ...
];


reactor.compoundsInflow = ... %mM
[
0.0 ...
0.0 ...
0.0 ...
];

reactor.biomassInflow = ... %g / l
[
0.0 ...
];

end

