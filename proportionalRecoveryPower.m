load('Output/NumberOfGroupsIs10.mat')
load('Output/ResultsNumberOfGroupsIs10.mat')

minSize = 10;
stepSize = 10;
NSamples = 500;
NClassicalSamples = 10^4;
mult = 25;

maxSize = 140;
steps = minSize:stepSize:maxSize;
NSteps = numel(steps);

classicalMaxSize = 700;
classicalSteps = minSize:stepSize:classicalMaxSize;
NClassicalSteps = numel(classicalSteps);

trEff = 4.25;
fmLim = log((52.5/66)/(1-52.5/66));
fu(1).tt = [4 12 26];
fu(2).tt = 26;
cl(1).tt = 1;
cl(2).tt = 2;
cl(3).tt = [1 2];
powerBugsAnalysis = [3 0 0 1 1 0;1 1 1 1 1 1];
NPBA = size(powerBugsAnalysis,2);
t = [1 2 4 12 26]';
NFU = size(fu,2);
NCL = size(cl,2);

nonClustOutput(1).trEffEst = NaN(NSamples,NSteps,4,NFU,NPBA);
nonClustOutput(1).plEffEst = NaN(NSamples,NSteps,2,NFU,NPBA);
nonClustOutput(1).gagrEst = NaN(NSamples,NSteps,2,NFU,NPBA);
nonClustOutput(1).cagrEst = NaN(NSamples,NSteps,2,NFU,NPBA);
nonClustOutput(1).power = NaN(NSamples,NSteps,4,NFU,NPBA);
nonClustOutput(1).ys = NaN(NSamples,NSteps,2,NFU,NPBA);

clustOutput(1).trEffEst = NaN(NSamples,NSteps,3,4,NCL,NFU);
clustOutput(1).plEffEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).gAcc = NaN(NSamples,NSteps);
clustOutput(1).cAcc = NaN(NSamples,NSteps,NSteps);
clustOutput(1).gagrEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).cagrEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).missrate = NaN(NSamples,3,NCL);
clustOutput(1).ppv = NaN(NSamples,3,NCL);
clustOutput(1).power = NaN(NSamples,NSteps,3,4,NCL,NFU);

classicalOutput(1).power = NaN(NSamples,NClassicalSteps,4);

params.r = mod.rETI(:,1);
params.tau = mod.tauETI(:,1);
params.alpham = mod.alphamETI(:,1);
params.alphap = mod.alphapETI(:,1);
params.alphas = mod.alphasETI(:,1);
params.yp = mod.yp(1);
params.clust = mod.clust;
params.gp = mod.gpETI(:,1);

generateSet.id = [ones(numel(t),1);2*ones(numel(t),1)];
generateSet.tt = repmat(t,2,1);
generateSet.FMt = NaN(size(generateSet(1).tt));
generateSet.trGr = [0; 1];
generateSet.tStart = repmat(3,2,1);

tic
paramsAlt = params;
paramsAlt.yp = paramsAlt.yp*1./(1.25)^2;

dataSet(1) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,params,0,0);
dataSet(2) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,params,1,0);
dataSet(3) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,params,3,0);
dataSet(4) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,params,0,2);
dataSet(5) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,params,0,10);
dataSet(6) = generateBugs(generateSet,trEff,fmLim,NSamples,mult*maxSize/2,paramsAlt,0,0);
dataSet(7) = generateBugs(generateSet,trEff,fmLim,NClassicalSamples,classicalMaxSize/2,params,0,0);

fprintf('Data generation time=%3.1f\n',toc)

%% Classical analysis
for ii=1:NClassicalSamples
    for n=1:NClassicalSteps
        nn = classicalSteps(n);
        gr(1).subs = find(dataSet(7).trGr==0);
        gr(2).subs = find(dataSet(7).trGr==1);
        randSub = sort([gr(1).subs(randperm(numel(gr(1).subs),nn/2)); gr(2).subs(randperm(numel(gr(2).subs),nn/2))]);
        randInd = ismember(dataSet(7).id,randSub)&ismember(dataSet(7).tt(:,ii),[cl(1).tt fu(1).tt(end)]);
        FMtsub = reshape(dataSet(7).FMt(randInd,ii),2,[]);
        trGrsub = dataSet(7).trGr(randSub);
        deltaFM = diff(FMtsub);
        endpointFM = FMtsub(2,:);
        classicalOutput(1).power(ii,n,1) = ranksum(endpointFM(trGrsub==0),endpointFM(trGrsub==1))<0.05;
        classicalOutput(1).power(ii,n,2) = ranksum(deltaFM(trGrsub==0),deltaFM(trGrsub==1))<0.05;
        classicalOutput(1).power(ii,n,3) = ttest2(endpointFM(trGrsub==0),endpointFM(trGrsub==1));
        classicalOutput(1).power(ii,n,4) = ttest2(deltaFM(trGrsub==0),deltaFM(trGrsub==1));
        classicalOutput(1).std(ii,n,1) = std(endpointFM(trGrsub==0|trGrsub==1));
        classicalOutput(1).std(ii,n,2) = std(deltaFM(trGrsub==0|trGrsub==1));
    end
end

%% Longitudinal model analysis
for ii=1:NSamples
    display(['Sample=' int2str(ii)])
  
    %Do not cluster
    for n=1:NSteps
        nn = steps(n);

        for ss=1:NFU
            tic
            for zz=1:NPBA
                if powerBugsAnalysis(ss,zz)>0
                    gr(1).subs = find(dataSet(zz).trGr==0);
                    gr(2).subs = find(dataSet(zz).trGr==1);
                    randSub = sort([gr(1).subs(randperm(numel(gr(1).subs),nn/2)); gr(2).subs(randperm(numel(gr(2).subs),nn/2))]);
                    randInd = ismember(dataSet(zz).id,randSub)&ismember(dataSet(zz).ttOrg,[cl(1).tt fu(ss).tt]);
                    dataSet(zz).FMtsub = dataSet(zz).FMt(randInd,ii);
                    dataSet(zz).ttsub = dataSet(zz).tt(randInd,ii);
                    dataSet(zz).idsub = changem(dataSet(zz).id(randInd),1:nn,unique(randSub));
                    dataSet(zz).gsub = dataSet(zz).g(randSub);
                    dataSet(zz).csub = dataSet(zz).c(randSub);
                    dataSet(zz).trGrsub = dataSet(zz).trGr(randSub);
                    dataSet(zz).tStartsub = dataSet(zz).tStart(randSub,ii);

                    %Noise combinations
                    try
                        [nonClustOutput.trEffEst(ii,n,1:4,ss,zz),nonClustOutput.plEffEst(ii,n,1:2,ss,zz),nonClustOutput.power(ii,n,1:4,ss,zz),...
                        nonClustOutput.gagrEst(ii,n,1:2,ss,zz),nonClustOutput.cagrEst(ii,n,1:2,ss,zz),nonClustOutput.ys(ii,n,1:2,ss,zz)] = powerBugs(dataSet(zz),params,powerBugsAnalysis(ss,zz));
                    catch
                        disp('error')
                    end
                end
            end
            fprintf('No clustering Sample=%3.0f Subjects=%3.0f Follow-up=%1.0f Time=%3.1f\n',[ii nn ss toc])
        end
    end
    
    %Cluster
    for mm=1:NCL
        dataSet(1).FMtsub = dataSet(1).FMt(ismember(dataSet(1).ttOrg,cl(mm).tt),ii);
        dataSet(1).ttsub = dataSet(1).tt(ismember(dataSet(1).ttOrg,cl(mm).tt),ii);
        dataSet(1).idsub = dataSet(1).id(ismember(dataSet(1).ttOrg,cl(mm).tt));
        dataSet(1).gsubs = dataSet(1).g(unique(dataSet(1).idsub),ii);
        dataSet(1).csubs = dataSet(1).c(unique(dataSet(1).idsub),ii);
        
        clusterOutput = clusterBugs(dataSet(1),fmLim,params);
        clustOutput(1).gAcc(ii,mm) = sum(dataSet(1).gsubs==clusterOutput(1,:)')./numel(dataSet(1).gsubs);
        clustOutput(1).cAcc(ii,mm) = sum(dataSet(1).csubs==clusterOutput(2,:)')./numel(dataSet(1).csubs);
        for cc=1:3
            clustOutput(1).missrate(ii,cc,mm) = sum(clusterOutput(2,dataSet(1).csubs==cc)~=cc)./sum(dataSet(1).csubs==cc);
            clustOutput(1).ppv(ii,cc,mm) = sum(dataSet(1).csubs(clusterOutput(2,:)==cc)==cc)./sum(clusterOutput(2,:)==cc);
            gr(1).subs = find(dataSet(1).trGr==0 & clusterOutput(2,:)'==cc);
            gr(2).subs = find(dataSet(1).trGr==1 & clusterOutput(2,:)'==cc);
                            
            for n=1:NSteps
                nn = steps(n);
                if numel(gr(1).subs)>=nn/2&&numel(gr(2).subs)>=nn/2
                    for ss=1:NFU
                        tic
                        randSub = sort([gr(1).subs(randperm(numel(gr(1).subs),nn/2));gr(2).subs(randperm(numel(gr(2).subs),nn/2))]);
                        randInd = ismember(dataSet(1).id,randSub)&ismember(dataSet(1).ttOrg,[cl(mm).tt fu(ss).tt]);
                        dataSet(1).FMtsub = dataSet(1).FMt(randInd,ii);
                        dataSet(1).ttsub = dataSet(1).tt(randInd,ii);
                        dataSet(1).idsub = changem(dataSet(1).id(randInd),1:nn,unique(randSub));
                        dataSet(1).gsub = dataSet(1).g(randSub,ii);
                        dataSet(1).csub = dataSet(1).c(randSub,ii);
                        dataSet(1).trGrsub = dataSet(1).trGr(randSub);
                        dataSet(1).tStartsub = dataSet(1).tStart(randSub);
                        try
                            [clustOutput(1).trEffEst(ii,n,cc,1:4,mm,ss),clustOutput(1).plEffEst(ii,n,cc,1:2,mm,ss),clustOutput(1).power(ii,n,cc,1:4,mm,ss),...
                            clustOutput(1).gagrEst(ii,n,cc,1:2,mm,ss),clustOutput(1).cagrEst(ii,n,cc,1:2,mm,ss),clustOutput(1).ys(ii,n,cc,1:2,mm,ss)] = powerBugs(dataSet(1),params,2);
                        catch
                            disp('error')
                        end
                        fprintf('Clustering Sample=%3.0f Cluster=%3.0f Subjects=%3.0f Follow-up=%3.0f Time=%3.1f\n',[ii cc nn ss toc])
                    end
                end
            end
        end
    end
end
save('Output/proportionalRecoveryPowerResults.mat','clustOutput','nonClustOutput','classicalOutput')