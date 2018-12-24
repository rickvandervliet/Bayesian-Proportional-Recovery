function output = predictionBugs(FMt,tt,id,ttfinal,ttinitial,FMfinal,FMinitial,mod)

    NSubjects = id(end);
    NMeasurements = numel(FMt);
    NChains = 1;
    NBurnin  = 10^3;
    NSamples  = 10^4;
    doparallel = 0;
    
    dataStruct = struct('y', FMt, ...
        'NSubjects', NSubjects,...
        'NMeasurements', NMeasurements,...
        't', tt./7,...
        'id',id,...
        'alpham', squeeze(mod.alpham),...
        'alphap', squeeze(mod.alphap),...
        'r', squeeze(mod.r),...
        'yp', mod.yp,...
        'gp', squeeze(mod.gp),...
        'tau',squeeze(mod.tau),...
        'ttfinal',ttfinal./7,...
        'ttinitial',ttinitial./7);
    
    init0(1:NChains) = struct('alphaL', cell(1, 1), 'g', cell(1, 1));
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(NSubjects,1);
        init0(i).g = round(numel(mod.clust)*rand(NSubjects,1)+0.5);
    end;  

    [samples, stats, ~] = matjags( ...
    dataStruct, ...                     % Observed data   
    fullfile(pwd, 'proportionalRecoveryPrediction.txt'), ...    % File that contains model definition
    init0, ...                          % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', NChains,...              % Number of MCMC chains
    'nburnin', NBurnin,...              % Number of burnin steps
    'nsamples', NSamples, ...           % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 0, ...                       % Do the DIC?
    'monitorparams', {'yfinalpred','ydiffpred','g'}, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1);                    % clean up of temporary files?
    
    output(1,:) = stats.mean.yfinalpred';
    output(2,:) = abs(FMfinal-stats.mean.yfinalpred');
    output(3,:) = stats.mean.ydiffpred';
    output(4,:) = abs((FMfinal-FMinitial)-stats.mean.ydiffpred');
    output(5,:) = mode(reshape(samples.g,[],NSubjects));
    output(6,:) = mode(changem(squeeze(samples.g),mod.clust,1:5));
end