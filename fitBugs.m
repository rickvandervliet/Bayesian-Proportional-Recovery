function params = fitBugs(set,NPostSamples,knnmdl,mod)
    FMt = set.FMt;
    tt = set.tt;
    id = set.id;
    NIHSS = set.NIHSS;
    bamford = set.bamford;
    sa = set.sa;
    fe = set.fe;
    
    NSubjects = id(end);
    NMeasurements = numel(FMt);
    NChains = 1;
    NBurnin  = 10^3;
    NSamples  = 10^3;
    doparallel = 0;
    numGroups = 5;
    
    dataStruct = struct('y', FMt, ...
        'NSubjects', NSubjects,...
        'NMeasurements', NMeasurements,...
        'NGroups', numGroups,...
        't', tt,...
        'gpalpha', 1.6*ones(numGroups,1),...
        'id',id);
    
    init0(1:NChains) = struct('alphaL', cell(1, 1),'rL', cell(1, 1),'gp', cell(1, 1),...
        'tauShift', cell(1, 1),'yp', cell(1, 1),'g', cell(1, 1),'alpham', cell(1, 1),'alphap', cell(1, 1));
    
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(NSubjects,1);
        init0(i).g = round(numGroups*rand(NSubjects,1)+0.5);
        
        init0(i).rL = log(mod.rETI(:,1)./(1-mod.rETI(:,1)));
        init0(i).tauShift = mod.tauETI(:,1)-1/7;
        init0(i).yp = mod.yp(1);
        init0(i).alpham = mod.alphamETI(:,1);  
        init0(i).alphap = mod.alphapETI(:,1);
        init0(i).gp = mod.gpETI(:,1);  
    end;    

    [samples, stats, ~] = matjags( ...
    dataStruct, ...                     % Observed data   
    fullfile(pwd, 'proportionalRecovery.txt'), ...    % File that contains model definition
    init0, ...                          % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', NChains,...              % Number of MCMC chains
    'nburnin', NBurnin,...              % Number of burnin steps
    'nsamples', NSamples, ...           % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 0, ...                       % Do the DIC?
    'monitorparams', {'yp','g','gp',...
    'tau','r','alpham','alphap'}, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1);                    % clean up of temporary files?
    
    randInd = randperm(NSamples,NPostSamples);
    params = struct('yp', cell(1, NPostSamples+1),'gp', cell(1, NPostSamples+1),...
        'tau', cell(1, NPostSamples+1),'r', cell(1, NPostSamples+1),'alpham', cell(1, NPostSamples+1),...
        'alphap', cell(1, NPostSamples+1),'clust', cell(1, NPostSamples+1));
    
    params(1).yp = stats.mean.yp;
    params(1).gp = stats.mean.gp;
    params(1).tau = stats.mean.tau;
    params(1).r = stats.mean.r;
    params(1).alpham = stats.mean.alpham;
    params(1).alphap = stats.mean.alphap;
    params(1).clust = predict(knnmdl,[params(1).r', params(1).tau', params(1).alpham', params(1).alphap']);
    g = median(squeeze(samples.g(1,:,:)));
    for gg=1:numGroups
        NIHSSgamma = fitdist(NIHSS(g==gg)','gamma');
        params(1).NIHSSgamma.a(gg) = NIHSSgamma.a;
        params(1).NIHSSgamma.b(gg) = 1./NIHSSgamma.b;
        params(1).bamfordp(gg,:) = [mean(bamford(g==gg)==1),mean(bamford(g==gg)==2),mean(bamford(g==gg)==3)];
        params(1).fep(gg) = mean(fe(g==gg));
        params(1).sap(gg) = mean(sa(g==gg));
    end;  
        
    for ii=1:NPostSamples
        params(ii+1).yp = samples.yp(1,randInd(ii));
        params(ii+1).gp = samples.gp(1,randInd(ii),:);
        params(ii+1).tau = samples.tau(1,randInd(ii),:);
        params(ii+1).r = samples.r(1,randInd(ii),:);
        params(ii+1).alpham = samples.alpham(1,randInd(ii),:);
        params(ii+1).alphap = samples.alphap(1,randInd(ii),:);
        g = squeeze(samples.g(1,randInd(ii),:));
        for gg=1:numGroups
            NIHSSgamma = fitdist(NIHSS(g==gg)','gamma');
            params(ii+1).NIHSSgamma.a(gg) = NIHSSgamma.a;
            params(ii+1).NIHSSgamma.b(gg) = 1./NIHSSgamma.b;
            params(ii+1).bamfordp(gg,:) = [mean(bamford(g==gg)==1),mean(bamford(g==gg)==2),mean(bamford(g==gg)==3)];
            params(ii+1).fep(gg) = mean(fe(g==gg));
            params(ii+1).sap(gg) = mean(sa(g==gg));
        end;    
        params(ii+1).clust = params(1).clust;
    end;
    
end