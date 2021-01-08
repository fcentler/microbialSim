function [ trajectory ] = microbialSimMain(scenarioID)
%microbialSimMain v1.1.1
%   This is the main routine of the dFBA simulator �bialSim. Different
%   simulation scenarios can be set up by configuring "scenarioID"s. Out of
%   the box, �bialSim can simulate three batch cultures:
%   scenarioID == 1: a Methanococcus maripaludis single culture (Richards et al., 2016)
%   scenarioID == 2: a syntrophic methanogenesis biculture (Hamilton et al., 2015)
%   scenarioID == 3: 8-species human gut consortium (SIHUMIx) with models taken from the
%   collection of 773 human gut microbiome species (Magnusdottir et al., 2017)

%% PARAMETER
solverPars.FBAsolver = 2; % Select FBA solver, set to 1 for CellNetAnalyzer
%   and to 2 for COBRA Toolbox
solverPars.tend = 1; % Simulation end time (in hours)
solverPars.timeStepSize=0.002; % Default simulation step size (in hours)
solverPars.saveLoadedModelsToFile = 0; % Set to 1 to save loaded SBML
%   models in Matlab format for speeding up subsequent simulation runs.
%   Models are stored in "loadedModels<date><time>.mat".
solverPars.readInitialStateFrom = ''; % Provide trajectory file name to
%   continue a previous simulation run ("simulatedTrajectory_*.mat_restartInit.mat").
solverPars.parallel = 0; % Set to 1 to compute FBA models in parallel using
%   Matlab's spmd environment
solverPars.maxWorkers = 12; % maximum number of workers to recruit for
%   parallel pool
solverPars.dopFBA = 1; % Minimize sum of total flux in each FBA computation
%   (pFBA). Set to 0 to deactivate.
solverPars.doMin2PrevFlux = 1; % Minimize deviation to flux distribution of
%   previous time step for each model. Set to 0 to deactivate.
solverPars.fluxTolerance = 1e-6; % During minimization of deviation to
%   previous flux, allow for this deviation for unaltered fluxes. Also
%   used as the tolerance when reporting limiting fluxes.
solverPars.minimalGrowth = 1e-6; % If predicted growth rate is below this
%   threshold, assume that no growth is possible
solverPars.solverType = 0; % Set to 0 for augmented forward Euler method,
%   1 for direct approach (using Matlab's ODE solver)
solverPars.doMassBalance = 1; % Set to 0 to avoid calculation of mass balances
%   for all exchange fluxes at the end of the simulation (can take a long
%   time)

% Only used if augmented forward Euler method is used:
solverPars.myAccuracy = 1e-15; % Compound concentrations below this value
%   are evaluated to be zero.
solverPars.myBioAccuracy = 1e-15; % Biomass concentrations below this 
%   vaue are evaluated to be zero.
solverPars.SteadyStateAccuracy = 1e-5; % Allow for this deviation between 
%   production and consumption rate when computing steady state during
%   time step size reduction
solverPars.maxDeviation = 5.0; % For compounds of high demand, reduce
%   time step size so that compound change is not bigger than this value
%   (specified in %). This helps to avoid oscillatory behavior. To disable
%   this feature and speed up simulation, set to "inf"
solverPars.biomassReductionFactor = 2.0; % If negative biomass concentrations
%   occur, reduce time step size such that at new time step, biomass
%   becomes old biomass divided by this factor.
solverPars.recordLimitingFluxes = 1; % record during the simulation which
%   fluxes are at their boundary
solverPars.recording = 0; % Set to 1 to save reactor state at each iteration
%   to separate files.

% Only used if COBRA is used:
solverPars.minRelaxValue = 1e-11; % If minimization of deviation to
%   previous flux fails, relax restrictions, starting with this value
%   (allow small zero-boundary violations by this amount) ...
solverPars.maxRelaxValue = 1; % ... and stop at this value.

% Only used for direct approach:
solverPars.solver = 2; % Set to 1 for ode45 solver, to 2 for ode15s solver.
solverPars.nonNegative = 0; % Set to 1 to activate NonNegative option for
%   solver.
solverPars.relTol = 1e-9; % Set relative tolerance
solverPars.absTol = 1e-9; % Set absolute tolerance

% Filenames
fileSuffix = strrep(char(datetime('now')),':', '_');
modelFile = strcat('loadedModels_', fileSuffix, '.mat');
solverPars.trajectoryFile = strcat('simulatedTrajectory_', fileSuffix, '.mat');

for name = {modelFile, solverPars.trajectoryFile}
    if exist(char(name), 'file') == 2
        error(strcat('Filename ', char(name), ' already exists. Aborting.'))
    end
end
%% FBA TOOL BOX SETUP

% initialize FBA solver
switch solverPars.FBAsolver
    case 1
        returnValue = initCNA;
    case 2
    if 2 ~= exist('optimizeCbModel')
        currDir = pwd;
        % ADD PATH TO YOUR COBRA INSTALLATION HERE (REPLACE OR APPEND AT END OF LIST)
        placesToLookForCobraToolbox = {
            'W:/cobratoolbox', ...
            'C:/Users/centlerf/cobratoolbox', ...
            'Y:/Home/centlerf/Projects/MetabolicModeling/FBATools/cobratoolbox', ...
            '/data/cobratoolbox' ...
            '../cobratoolbox' ...
            '/Users/flori/cobratoolbox'
            };

        for i = 1:length(placesToLookForCobraToolbox)
            if 7 == exist(char(placesToLookForCobraToolbox(i)), 'dir') % directory exists
                break
            end
            if i == length(placesToLookForCobraToolbox)
                error('COBRA Toolbox not found. Aborting.')
            end
        end
        
        cd(char(placesToLookForCobraToolbox(i)))
        initCobraToolbox(false); % no updating
        cd(currDir);
    else
        changeCobraSolver('glpk', 'LP');
    end

    otherwise
        error('FBA solver type unkown!')
end

%% SCENARIO SPECIFIC SETUP
clearvars models
workerModels = 0;

switch scenarioID
    case 1 % M. maripaluids monoculture 
        solverPars.tend = 1.0;
        solverPars.timeStepSize = 0.002;
        solverPars.solverType = 0; % 0 dFBA or 1 for ODE
        
        model = prepareFBAmodel_iMR539('models/MODELRICHARDS 2016.xml', solverPars.FBAsolver);
        model = convertToOnlyIrrevReac(model);
        model = parametrizeFBAmodel_iMR539(model);
        models(1) = model;

        reactor = reactorDefinition_iMR539;
        
        %Examplary usage of function estimateVmax:
        %[myVmax, myFlux, myExchangeFluxes] = estimateVmax(models(1), reactor, solverPars, 2, 0.0875);
        %trajectory = myExchangeFluxes;
        %return
        
    case 2 % syntrophic propionate degrading binary consortium
        solverPars.tend = 1.0;
        solverPars.timeStepSize = 0.002;
        solverPars.solverType = 0; % 0 dFBA or 1 ODE
        
        iSfu648Model = prepareFBAmodel_iSfu648Model('models/S5_Dataset_SDH_Coupling.xml', solverPars.FBAsolver);
        iSfu648Model = convertToOnlyIrrevReac(iSfu648Model);
        iSfu648Model = parametrizeFBAmodel_iSfu648(iSfu648Model);
        iMhu428Model = prepareFBAmodel_iMhu428Model('models/S3_Dataset.xml', solverPars.FBAsolver);
        iMhu428Model = convertToOnlyIrrevReac(iMhu428Model);
        iMhu428Model = parametrizeFBAmodel_iMhu428(iMhu428Model);
        models(1) = iSfu648Model;
        models(2) = iMhu428Model;
        
        reactor = reactorDefinition_synProp;
    case 3 % 773 gut microbiome
        solverPars.tend = 1.0;
        solverPars.timeStepSize = 0.002;

        % consider all species
        %speciesToConsider = 1:773;

        % consider only selected species
        %speciesToConsider = [1, 2, 3];
        speciesToConsider = [46, 125, 164, 178, 245, 270, 362, 482];
        
        [models, externalCompounds] = prepareFBAmodel_773Agora(speciesToConsider, 'models/AGORA-1.01-Western-Diet/', solverPars.FBAsolver);
        solverPars.saveLoadedModelsToFile = 1; % look for "loadedModels_*.mat" in current directory, can be used in subsequent runs for speed-up (see next line)
        %load('./loadedModels_26-Mar-2019 11_34_00.mat');

        solverPars.readInitialStateFrom = '';

        solverPars.parallel = 0;
        solverPars.recording = 0;
        
        % set KS to nonzero
        ksvalue = 0.01;
        for i = 1:length(models)
            for j = 1:length(models(i).coupledReactions.ks)
                models(i).coupledReactions.ks(j) = ksvalue;
            end
        end

        reactor = reactorDefinition_773Agora(length(models), externalCompounds);

    otherwise
        error('Unknown model scenario. Aborting')
end

if solverPars.parallel == 1
    myPool = parpool(min(length(models), solverPars.maxWorkers));

    global CBT_LP_SOLVER
    solver = CBT_LP_SOLVER;
    parfor i = 1:length(models)
        changeCobraSolver(solver, 'LP', 1, -1);
    end
    
    % each worker gets his share of FBA models
    % (memory requirements could be halfed if models were directly created
    % within workers)
    workerModels = Composite();
    for i = 1:numel(workerModels)
        workerModels{i} = models(i:numel(workerModels):length(models));
    end
end

if solverPars.saveLoadedModelsToFile == 1
    switch scenarioID
        case 3
            save(modelFile, 'models', 'externalCompounds', '-v7.3'); % for > 2GB
        otherwise
            save(modelFile, 'models', '-v7.3'); % for > 2GB
    end
end

tic
if solverPars.parallel == 1
    trajectory = dFBASimulator(models, reactor, solverPars.solverType, solverPars, workerModels);
else
    trajectory = dFBASimulator(models, reactor, solverPars.solverType, solverPars, 0);
end
elapsedSeconds = toc;
elapsedHours = elapsedSeconds / 60 / 60;

save(solverPars.trajectoryFile, 'trajectory', 'reactor', 'elapsedHours', '-v7.3');

%% plotting
% Quick-n-dirty-Plot
figure
subplot(2,1,1)
plot(trajectory.time, trajectory.biomass)
legend(trajectory.modelNames, 'Interpreter', 'none')
title('Biomass')
xlabel('Time (hours)')
ylabel('Concentration (gDW/L)')
subplot(2,1,2)
plot(trajectory.time, trajectory.compounds)
legend(trajectory.compoundNames, 'Interpreter', 'none')
title('Compounds')
xlabel('Time (hours)')
ylabel('Concentration (mmol/L)')

plotTrajectoryCmp(trajectory, 'Simulation')
plotExchangeFluxes(trajectory)

if solverPars.parallel == 1
   %delete(gcp('nocreate'))
   delete(myPool);
end

end
