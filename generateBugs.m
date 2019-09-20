function outputSet = generateBugs(inputSet,trEff,fmLim,NSamples,NMult,mod,simNoise,tauFixed)

    NSubjects = inputSet.id(end);
    NMeasurements = numel(inputSet.FMt);
    NChains = 1;
    NBurnin  = 100;
    doparallel = 0;
    gp = mod.gp.*normcdf(fmLim,mod.alpham,1./sqrt(mod.alphap));
    gp = gp./(repmat(sum(gp,1),5,1));
    
    dataStruct = struct('y', inputSet.FMt, ...
        'trGr', inputSet.trGr,...
        'trEff', trEff,...
        'fmLim', fmLim,...
        'NSubjects', NSubjects,...
        'NMeasurements', NMeasurements,...
        't', inputSet.tt,...
        'id',inputSet.id,...
        'tStart',inputSet.tStart,...
        'alpham', mod.alpham,...
        'alphap', mod.alphap,...
        'r', mod.r,...
        'yp', mod.yp,...
        'gp', gp,...
        'tau', mod.tau,...
        'tauFixed', tauFixed,...
        'NPostSamp',size(mod.r,2),...
        'simMeasNoise',simNoise==1|simNoise==3,...
        'simTreatNoise',simNoise==2|simNoise==3);
    
    init0(1:NChains) = struct('alphaL', cell(1, 1), 'g', cell(1, 1));
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(NSubjects,1);
        init0(i).g = round(size(mod.r,1)*rand(NSubjects,1)+0.5);
    end;  
    
    if tauFixed > 0
        [samples, ~, ~] = matjags( ...
        dataStruct, ...                     % Observed data   
        fullfile(pwd, 'proportionalRecoveryGenerationTauFixed.txt'), ...    % File that contains model definition
        init0, ...                          % Initial values for latent variables
        'doparallel' , doparallel, ...      % Parallelization flag
        'nchains', NChains,...              % Number of MCMC chains
        'nburnin', NBurnin,...              % Number of burnin steps
        'nsamples', NSamples*NMult, ...           % Number of samples to extract
        'thin', 1, ...                      % Thinning parameter
        'dic', 0, ...                       % Do the DIC?
        'monitorparams', {'y','g','tNoise','tStartNoise'}, ...   % List of latent variables to monitor
        'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
        'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
        'cleanup' , 1);      
    else    
        [samples, ~, ~] = matjags( ...
        dataStruct, ...                     % Observed data   
        fullfile(pwd, 'proportionalRecoveryGeneration.txt'), ...    % File that contains model definition
        init0, ...                          % Initial values for latent variables
        'doparallel' , doparallel, ...      % Parallelization flag
        'nchains', NChains,...              % Number of MCMC chains
        'nburnin', NBurnin,...              % Number of burnin steps
        'nsamples', NSamples*NMult, ...           % Number of samples to extract
        'thin', 1, ...                      % Thinning parameter
        'dic', 0, ...                       % Do the DIC?
        'monitorparams', {'y','g','tNoise','tStartNoise'}, ...   % List of latent variables to monitor
        'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
        'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
        'cleanup' , 1);                    % clean up of temporary files?
    end
    outputSet.FMt = reshape(permute(samples.y,[3,2,1]),[],NSamples);
    outputSet.FMt(outputSet.FMt<0) = 0;
    outputSet.FMt(outputSet.FMt>66) = 66;
    outputSet.tt = reshape(permute(samples.tNoise,[3,2,1]),[],NSamples);
    outputSet.ttOrg = repmat(inputSet.tt,size(outputSet.tt,1)/numel(inputSet.tt),1);
    outputSet.tStart = reshape(permute(samples.tStartNoise,[3,2,1]),[],NSamples);
    outputSet.g = reshape(permute(samples.g,[3,2,1]),[],NSamples);
    outputSet.id = reshape(repmat((1:size(outputSet.g,1)),size(outputSet.tt,1)./size(outputSet.g,1),1),[],1);
    outputSet.trGr = reshape(repmat((0:1)',size(outputSet.g,1)/2,1),[],1);
    outputSet.c = changem(outputSet.g,mod.clust,1:numel(mod.clust));
end