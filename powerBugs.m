function [trEff,plEff,p,gagr,cagr,ys] = powerBugs(set,mod,analysis)
    NSubjects = set.idsub(end);
    NMeasurements = size(set.FMtsub,1);
    NChains = 4;
    NBurnin  = 10^4;
    NSamples  = 2.5*10^3;
    doparallel = 1;
    c = changem(set.gsub,mod.clust,1:5);
    
    dataStruct = struct('y', set.FMtsub, ...
        'gorg', set.gsub,...
        'corg', c,...
        'trGr', set.trGrsub,...
        'tStart',set.tStartsub,...
        'NSubjects', NSubjects,...
        'NMeasurements', NMeasurements,...
        't', set.ttsub,...
        'id',set.idsub,...
        'alpham', mod.alpham(:,1),...
        'alphap', mod.alphap(:,1),...
        'r', mod.r(:,1),...
        'gp', mod.gp(:,1),...
        'tau', mod.tau(:,1),...
        'clust',mod.clust,...
        'NClusters',numel(unique(mod.clust)));
    init0(1:NChains) = struct('alphaL', cell(1, 1), 'g', cell(1, 1));
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(NSubjects,1);
        init0(i).g = round(size(mod.r,1)*rand(NSubjects,1)+0.5);
    end;  
    
    if analysis==1 || analysis ==3
        %% Single treatment effect
        [~, stats, ~] = matjags( ...
        dataStruct, ...                     % Observed data   
        fullfile(pwd, 'proportionalRecoveryPower.txt'), ...    % File that contains model definition
        init0, ...                          % Initial values for latent variables
        'doparallel' , doparallel, ...      % Parallelization flag
        'nchains', NChains,...              % Number of MCMC chains
        'nburnin', NBurnin,...              % Number of burnin steps
        'nsamples', NSamples, ...           % Number of samples to extract
        'thin', 1, ...                      % Thinning parameter
        'dic', 0, ...                       % Do the DIC?
        'monitorparams', {'trEff','plEff','gagr','cagr','p','ys'}, ...   % List of latent variables to monitor
        'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
        'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
        'cleanup' , 1);                    % clean up of temporary files?

        p(1) = stats.mean.p>0.975;
        trEff(1) = stats.mean.trEff;
        plEff(1) = stats.mean.plEff;
        gagr(1) = stats.mean.gagr;
        cagr(1) = stats.mean.cagr;
        ys(1) = stats.mean.ys;
    else
        p(1) = NaN;
        trEff(1) = NaN;
        plEff(1) = NaN;
        gagr(1) = NaN;
        cagr(1) = NaN;
        ys(1) = NaN;
    end;
    
    if analysis==2 || analysis==3
        %% Multiple treatment effect
        [~, stats, ~] = matjags( ...
        dataStruct, ...                     % Observed data   
        fullfile(pwd, 'proportionalRecoveryPowerClusters.txt'), ...    % File that contains model definition
        init0, ...                          % Initial values for latent variables
        'doparallel' , doparallel, ...      % Parallelization flag
        'nchains', NChains,...              % Number of MCMC chains
        'nburnin', NBurnin,...              % Number of burnin steps
        'nsamples', NSamples, ...           % Number of samples to extract
        'thin', 1, ...                      % Thinning parameter
        'dic', 0, ...                       % Do the DIC?
        'monitorparams', {'trEff','plEff','ys','gagr','cagr','p'}, ...   % List of latent variables to monitor
        'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
        'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
        'cleanup' , 1);                    % clean up of temporary files?
        p(2:4) = stats.mean.p>0.975;
        trEff(2:4) = stats.mean.trEff;
        plEff(2) = stats.mean.plEff;
        gagr(2) = stats.mean.gagr;
        cagr(2) = stats.mean.cagr;
        ys(2) = stats.mean.ys;
    end;
end