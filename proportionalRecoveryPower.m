load('Output/NumberOfGroupsIs10.mat')
load('Output/ResultsNumberOfGroupsIs10.mat')
minSize = 400;
maxSize = 400;
stepSize = 10;
NSamples = 1;
mult = 1;
steps = minSize:stepSize:maxSize;
NSteps = numel(steps);

trEff = 0;
fmLim = log((66/66)/(1-66/66));
fu(1).tt = [3 5 8 12 26];
fu(2).tt = 26;
cl(1).tt = 1;
cl(2).tt = 3;
t = [0 3/7 1 2 3 8 12 26]';
NFU = size(fu,2);
NCL = size(cl,2);

nonClustOutput(1).trEffEst = NaN(NSamples,NSteps,4,NFU,3);
nonClustOutput(1).plEffEst = NaN(NSamples,NSteps,2,NFU,3);
nonClustOutput(1).gagrEst = NaN(NSamples,NSteps,2,NFU,3);
nonClustOutput(1).cagrEst = NaN(NSamples,NSteps,2,NFU,3);
nonClustOutput(1).power = NaN(NSamples,NSteps,4,NFU,3);
nonClustOutput(1).ys = NaN(NSamples,NSteps,2,NFU,3);

clustOutput(1).trEffEst = NaN(NSamples,NSteps,3,4,NCL,NFU);
clustOutput(1).plEffEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).gAcc = NaN(NSamples,NSteps);
clustOutput(1).cAcc = NaN(NSamples,NSteps,NSteps);
clustOutput(1).gagrEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).cagrEst = NaN(NSamples,NSteps,3,2,NCL,NFU);
clustOutput(1).missrate = NaN(NSamples,NSteps,3);
clustOutput(1).ppv = NaN(NSamples,NSteps,3);
clustOutput(1).power = NaN(NSamples,NSteps,3,4,NCL,NFU);

params.r = mod.rETI(:,1);
params.tau = mod.tauETI(:,1);
params.alpham = mod.alphamETI(:,1);
params.alphap = mod.alphapETI(:,1);
params.alphas = mod.alphasETI(:,1);
params.yp = mod.yp(1);
params.clust = mod.clust;
params.gp = mod.gpETI(:,1);

set(3).id = reshape(repmat(1:maxSize*mult,numel(t),1),[],1);
set(3).tt = reshape(repmat(t,1,maxSize*mult),[],1);
set(3).FMt = NaN(size(set.tt));
set(3).trGr = [zeros(maxSize*mult/2,1); ones(maxSize*mult/2,1)];
set(3).tStart = 3*ones(maxSize*mult,1);
set(3).FMtsub = set.FMt(randInd,ii);
set(3).ttsub = NaN;
set(3).idsub = NaN;
set(3).gsub = NaN;
set(3).trGrsub = NaN;
set(3).tStartsub = NaN;
                
numPostSamp = size(samples.r,2);
samp = randperm(numPostSamp,1000);
randomParams.r = squeeze(samples.r(1,samp,mod.grInc))';
randomParams.tau = squeeze(samples.tau(1,samp,mod.grInc))';
randomParams.alpham = squeeze(samples.alpham(1,samp,mod.grInc))';
randomParams.alphap = squeeze(samples.alphap(1,samp,mod.grInc))';
randomParams.yp = squeeze(samples.yp(samp));
randomParams.gp = squeeze(samples.gp(1,samp,mod.grInc))';
params.r=repmat(params.r,1,2)
params.tau=repmat(params.tau,1,2)
params.alpham=repmat(params.alpham,1,2)
params.alphap=repmat(params.alphap,1,2)
params.alphas=repmat(params.alphas,1,2)
params.gp=repmat(params.gp,1,2)
params.yp=repmat(params.yp,1,2)
tic
[set(1).FMt, set(1).tt, set(1).g] = generateBugs(set,trEff,fmLim,NSamples,params,0);
%[set(2).FMt, set(2).tt, set(2).g] = generateBugs(set,trEff,fmLim,NSamples,randomParams,1);
%[set(3).FMt, set(3).tt, set(3).g] = generateBugs(set,trEff,fmLim,NSamples,randomParams,3);
toc
set(1).c = changem(set(1).g,params.clust,1:numel(params.clust));

for ii=1:NSamples
    display(['Sample=' int2str(ii)])
  
    %Do not cluster
    for n=1:NSteps
        nn = steps(n);
        for ss=1:NFU
            tic
            for zz=1:3
                randSub = [randperm(maxSize*mult/2,nn/2), (maxSize*mult)/2+randperm(maxSize*mult/2,nn/2)];
                randInd = ismember(set.id,randSub)&ismember(set.tt,[cl(1).tt fu(ss).tt]);
                set(zz).FMtsub = set.FMt(randInd,ii);
                set(zz).ttsub = set.tt(randInd,ii);
                set(zz).idsub = changem(set.id(randInd),1:nn,unique(randSub));
                set(zz).gsub = set.g(randSub);
                set(zz).trGrsub = set.trGr(randSub);
                set(zz).tStartsub = set.tStart(randSub);

                %Noise combinations
                [nonClustOutput.trEffEst(ii,n,1:4,ss,zz),nonClustOutput.plEffEst(ii,n,1:2,ss,zz),nonClustOutput.power(ii,n,1:4,ss,zz),...
                nonClustOutput.gagrEst(ii,n,1:2,ss,zz),nonClustOutput.cagrEst(ii,n,1:2,ss,zz),nonClustOutput.ys(ii,n,1:2,ss,1)] = powerBugs(set(zz),params,3);
            end
            fprintf('No clustering Sample=%3.0f Subjects=%3.0f Follow-up=%1.0f Time=%3.1f\n',[ii nn ss toc])
        end;
    end;
    
    %Cluster
    for mm=1:NCL
        set(1).FMtsub = set(1).FMt(ismember(set(1).tt,cl(mm).tt),ii);
        set(1).ttsub = set(1).tt(ismember(set(1).tt,cl(mm).tt));
        set(1).idsub = set(1).id(ismember(set(1).tt,cl(mm).tt));
        set(1).gsub = set(1).g(set(1).idsub,ii);
        set(1).csub = set(1).c(set(1).idsub,ii);
        
        clusterOutput = clusterBugs(set(1),fmLim,params);
        clustOutput(1).gAcc(ii,mm) = sum(set(1).gsub==clusterOutput(1,:)')./numel(set(1).gsub);
        clustOutput(1).cAcc(ii,mm) = sum(set(1).csub==clusterOutput(2,:)')./numel(set(1).csub);
        for cc=1:3
            clustOutput(1).missrate(ii,cc,mm) = sum(clusterOutput(2,set(1).csub==cc)~=cc)./sum(set(1).csub==cc);
            clustOutput(1).ppv(ii,cc,mm) = sum(set(1).csub(clusterOutput(2,:)==cc)==cc)./sum(clusterOutput(2,:)==cc);
            plGr = find(clusterOutput(2,1:maxSize*mult/2)==mm);
            trGr = maxSize*mult/2+find(clusterOutput(2,maxSize*mult/2+1:end)==mm);
            
            for n=1:NSteps
                nn = steps(n);
                if numel(plGr)>=nn/2&&numel(trGr)>=nn/2
                    for ss=1:NFU
                        tic
                        randSub = [plGr(randperm(numel(plGr),nn/2)),trGr(randperm(numel(trGr),nn/2))];
                        randInd = ismember(set(1).id,randSub)&ismember(set(1).tt,[mm fu(ss).tt]);
                        set(1).FMtsub = set(1).FMt(randInd,ii);
                        set(1).ttsub = set(1).tt(randInd,ii);
                        set(1).idsub = changem(set(1).id(randInd),1:nn,unique(randSub));
                        set(1).gsub = set(1).g(randSub,ii);
                        set(1).trGrsub = set(1).trGr(randSub);
                        set(1).tStartsub = set(1).tStart(randSub);
                        [clustOutput(1).trEffEst(ii,n,cc,1:4,mm,ss),clustOutput(1).plEffEst(ii,n,cc,1:2,mm,ss),clustOutput(1).power(ii,n,cc,1:4,mm,ss),...
                        clustOutput(1).gagrEst(ii,n,cc,1:2,mm,ss),clustOutput(1).cagrEst(ii,n,cc,1:2,mm,ss),clustOutput(1).ys(ii,n,cc,1:2,mm,ss)] = powerBugs(set(1),params,2);
                        fprintf('No clustering Sample=%3.0f Cluster=%3.0f Subjects=%3.0f Follow-up=%3.0f Time=%3.1f\n',[ii cc nn ss toc])
                    end;
                end;
            end;
        end;
    end;
end;
save('Output/proportionalRecoveryPowerResults.mat','clustOutput','nonClustOutput')