function output = clusterBugs(set,fmLim,mod)
    FMt = set.FMtsub;
    tt = set.ttsub;
    id = set.idsub;
    
    NSubjects = id(end);
    NMeasurements = numel(FMt);
    NChains = 5;
    NBurnin  = 10^3;
    NSamples  = 10^3;
    doparallel = 1;
    
    dataStruct = struct('y', FMt, ...
        'NSubjects', NSubjects,...
        'NMeasurements', NMeasurements,...
        't', tt,...
        'id',id,...
        'alpham', squeeze(mod.alpham),...
        'alphap', squeeze(mod.alphap),...
        'r', squeeze(mod.r),...
        'yp', mod.yp,...
        'gp', mod.gp.*normcdf(fmLim,mod.alpham,1./sqrt(mod.alphap)),...
        'tau',squeeze(mod.tau));
    
    init0(1:NChains) = struct('alphaL', cell(1, 1), 'g', cell(1, 1));
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(NSubjects,1);
        init0(i).g = round(numel(mod.clust)*rand(NSubjects,1)+0.5);
    end;  

    [samples, ~, ~] = matjags( ...
    dataStruct, ...                     % Observed data   
    fullfile(pwd, 'proportionalRecoveryCluster.txt'), ...    % File that contains model definition
    init0, ...                          % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', NChains,...              % Number of MCMC chains
    'nburnin', NBurnin,...              % Number of burnin steps
    'nsamples', NSamples, ...           % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 0, ...                       % Do the DIC?
    'monitorparams', {'g'}, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1);                    % clean up of temporary files?
    g = reshape(samples.g,[],NSubjects);
    output(1,:) = mode(g);
    output(2,:) = mode(changem(g,mod.clust,1:numel(mod.clust)));
end