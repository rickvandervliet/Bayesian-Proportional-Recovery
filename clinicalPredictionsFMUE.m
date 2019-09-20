function [ttpred,FMtpred,cmode,cp,FMtend] = clinicalPredictionsFMUE(tt,FMt)
    %% Fugl-Meyer predictions
    % Input: dataset.FMt - Fugl-Meyer scores
    % dataset.tt - Measurement times in weeks
    
    numOrgMeas = numel(tt);
    ttpred = (max(round(tt))+1):24;
    FMt = [reshape(FMt,1,[]), NaN(size(ttpred))];
    tt = [reshape(tt,1,[]), ttpred];
        
    NMeasurements = numel(FMt);
    NChains = 1;
    NBurnin  = 10^4;
    NSamples  = 5*10^4;
    doparallel = 0;
    
    alpham = [-3.1602   -2.0922   -2.7862   -1.3081    0.0201];
    alphap = [3.8376    0.2340    0.1394    0.1326    0.1774];
    tau = [5.3114   10.1417    9.8122    2.6722    1.1822];
    r = [0.0884    0.4591    0.8634    0.8883    0.9309];
    gp = [0.2674    0.1395    0.1118    0.1783    0.3031];
    yp = 0.0673;
    
    dataStruct = struct('y', FMt, ...
        'NMeasurements', NMeasurements,...
        'NGroups', numel(r),...
        't', tt,...
        'alpham', alpham,...
        'alphap', alphap,...
        'r', r,...
        'yp', yp,...
        'gp', gp,...
        'tau', tau);
    
    init0(1:NChains) = struct('alphaL', cell(1, 1), 'g', cell(1, 1));
    for i=1:NChains
        init0(i).alphaL = 0.5*ones(1);
        init0(i).g = round(numel(r)*rand(1)+0.5);
    end;  

    [samples, ~, ~] = matjags( ...
    dataStruct, ...                     % Observed data   
    fullfile(pwd, 'clinicalPredictionsFMUE.txt'), ...    % File that contains model definition
    init0, ...                          % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', NChains,...              % Number of MCMC chains
    'nburnin', NBurnin,...              % Number of burnin steps
    'nsamples', NSamples, ...           % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 0, ...                       % Do the DIC?
    'monitorparams', {'FMUE','g'}, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 0 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1);                    % clean up of temporary files?
    
    FMt = squeeze(samples.FMUE);
    FMtpred = determineETI(FMt(:,numOrgMeas+1:end),0.05);
    FMtend = FMt(:,end);
    c = changem(samples.g,[1 2 2 3 3],1:5);
    cmode = mode(changem(samples.g,[1 2 2 3 3],1:5));
    cp = 100*sum(c==cmode)./numel(c);
end