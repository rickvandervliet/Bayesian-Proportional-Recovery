load('Data/CombinedTransform.mat')

id = repmat(1:size(FMt,2),8,1);
selectNonNan = ~isnan(FMt);
FMt = FMt(selectNonNan);
tt = tt(selectNonNan);
id = id(selectNonNan);
NSubjects = numel(unique(id));
NMeasurements = numel(FMt);
NChains = 1;
doparallel = 0;
NBurnin  = 10^4; 
NSamples  = 10^4;
gpalpha = 1.6;

%% Overfitted model
numGroups = 10;

dataStruct = struct('y', FMt, ...
    'NSubjects', NSubjects,...
    'NMeasurements', NMeasurements,...
    'NGroups', numGroups,...
    't', tt./7,...
    'gpalpha', gpalpha*ones(numGroups,1),...
    'id',id);

init0(1:NChains) = struct('alphaL', cell(1, 1),'rL', cell(1, 1),'gp', cell(1, 1),...
    'tauShift', cell(1, 1),'yp', cell(1, 1),'g', cell(1, 1),'alphap', cell(1, 1),'alpham', cell(1, 1));
for i=1:NChains
    init0(i).alphaL = 0.5*ones(NSubjects,1);
    init0(i).g = round(numGroups*rand(NSubjects,1)+0.5);
    init0(i).rL = zeros(numGroups,1);
    init0(i).tauShift = 4*ones(numGroups,1);
    init0(i).alpham = -2*ones(numGroups,1);
    init0(i).alphap = 0.1*ones(numGroups,1);
    init0(i).gp = (1./numGroups)*ones(numGroups,1); 
    init0(i).yp = 10^-2;  
end;  

[samples, stats, structArray] = matjags( ...
dataStruct, ...                     % Observed data   
fullfile(pwd, 'proportionalRecovery.txt'), ...    % File that contains model definition
init0, ...                          % Initial values for latent variables
'doparallel' , doparallel, ...      % Parallelization flag
'nchains', NChains,...              % Number of MCMC chains
'nburnin', NBurnin,...              % Number of burnin steps
'nsamples', NSamples, ...           % Number of samples to extract
'thin', 1, ...                      % Thinning parameter
'dic', 0, ...                       % Do the DIC?
'monitorparams', {'yp','g','gpalpha','gp',...
'tau', 'alpha','r','alpham','alphap'}, ...   % List of latent variables to monitor
'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
'cleanup' , 1);                    % clean up of temporary files?
save(['Output/NumberOfGroupsIs' int2str(numGroups) '.mat'], 'samples', 'stats', 'structArray')
processProportionalRecoveryModel(samples)
clear samples stats structArray init0 dataStruct S