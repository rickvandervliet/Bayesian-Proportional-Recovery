load('Data/CombinedTransform.mat')
load('Output/ResultsNumberOfGroupsIs10.mat')
load('Output/clusterEstimate.mat')

timePoints = 1:12;
NPostSamples = 10;
numFolds = size(dem(2).FMt,2);
numTime = numel(timePoints);
output = NaN(NPostSamples+1,size(dem(2).FMt,2),numTime,6);
 
id = repmat(1:size(dem(2).FMt,2),8,1);
randInd = randperm(size(dem(2).FMt,2));
FMts = dem(2).FMt(:,randInd);
ids = id(:,randInd);
tts = dem(2).tt(:,randInd)./7;
selectNonNan = ~isnan(FMts);

knnmdl = fitcknn([mod.rETI(:,1), mod.tauETI(:,1), mod.alphamETI(:,1), mod.alphapETI(:,1)],[1 2 2 3 3],'Standardize',1);
fits = struct('params', cell(1, numFolds));

for nn=73:numFolds
    fprintf('Num=%i \n',nn)
    set(1).subs = round((nn-1).*(size(FMts,2)./numFolds))+1: round(nn*size(FMts,2)./numFolds);
    set(2).subs = (1:size(FMts,2));
    set(2).subs(ismember(set(2).subs,set(1).subs)) = [];
    
    for k=1:2
        set(k).selectNonNan = selectNonNan(:,set(k).subs);
        set(k).FMt = FMts(:,set(k).subs);
        set(k).FMt = set(k).FMt(set(k).selectNonNan);
        set(k).tt = tts(:,set(k).subs);
        set(k).tt = set(k).tt(set(k).selectNonNan);
        set(k).id = id(:,set(k).subs);
        set(k).id = set(k).id(set(k).selectNonNan);
        set(k).id = changem(set(k).id,1:numel(unique(set(k).id)),unique(set(k).id));
        set(k).orgid = ids(1,set(k).subs);
        [~, set(k).uniqueidLast] = unique(set(k).id,'last');
        [~, set(k).uniqueidFirst] = unique(set(k).id,'first');
        set(k).FMfinal = set(k).FMt(set(k).uniqueidLast);
        set(k).ttfinal = set(k).tt(set(k).uniqueidLast);
        set(k).FMinitial = set(k).FMt(set(k).uniqueidFirst);
        set(k).ttinitial = set(k).tt(set(k).uniqueidFirst);
        set(k).NIHSS = dem(2).NIHSS(set(k).orgid)+1;
        set(k).bamford = dem(2).bamford(set(k).orgid);
        set(k).sa = dem(2).SA(set(k).orgid);
        set(k).fe = dem(2).FE(set(k).orgid);
    end;

    fits(nn).params = fitBugs(set(2),NPostSamples,knnmdl,mod);
    
    for t=1:numTime
        time = timePoints(t);
        selectEarly = set(1).tt<=time;
        set(1).FMtsub = set(1).FMt(selectEarly);
        set(1).ttsub = set(1).tt(selectEarly);
        set(1).idsub = set(1).id(selectEarly);
        set(1).ttinitialsub = set(1).ttinitial(unique(set(1).idsub));
        set(1).ttfinalsub = set(1).ttfinal(unique(set(1).idsub));
        set(1).FMinitialsub = set(1).FMinitial(unique(set(1).idsub));
        set(1).FMfinalsub = set(1).FMfinal(unique(set(1).idsub));
        set(1).NIHSSsub = set(1).NIHSS(unique(set(1).idsub));
        set(1).bamfordsub = set(1).bamford(unique(set(1).idsub));
        set(1).sasub = set(1).sa(unique(set(1).idsub));
        set(1).fesub = set(1).fe(unique(set(1).idsub));
        
        set(1).idsub = changem(set(1).idsub,1:numel(unique(set(1).idsub)),unique(set(1).idsub));
        set(1).orgidsub = set(1).orgid(unique(set(1).idsub));
        
        if ~isempty(set(1).FMtsub)
            for rr=1:NPostSamples+1
                tic
                outputTemp = predictionBugs(set(1),fits(nn).params(rr));
                toc
                output(rr,set(1).orgidsub,t,:) = permute(outputTemp,[4,2,3,1]);
            end;
        end;
    end;
end;

FMtnn = dem(2).FMt(~isnan(dem(2).FMt));
idnn = id(~isnan(dem(2).FMt));
[~, uniqueidLast] = unique(idnn,'last');
[~, uniqueidFirst] = unique(idnn,'first');
FMinitial = FMtnn(uniqueidFirst);
FMfinal = FMtnn(uniqueidLast);
FMdiff = FMfinal - FMinitial;

results.correlation = NaN(NPostSamples,numTime,2);
results.correlationETI = NaN(3,numTime,2);
results.missrate = NaN(NPostSamples,numTime,3);
results.missrateETI = NaN(3,numTime,3);
results.ppv = NaN(NPostSamples,numTime,3);
results.ppvETI = NaN(3,numTime,3);
results.accuracy = NaN(NPostSamples,numTime);
results.accuracyETI = NaN(3,numTime);

for tt=1:numTime
    for rr=1:NPostSamples
        outputNonNan = ~isnan(output(rr+1,:,tt,1));
        clustETINonNan = clusterEstimate(outputNonNan);
        results.correlation(rr,tt,1) = corr(FMfinal(outputNonNan),output(rr+1,outputNonNan,tt,1)');
        results.correlation(rr,tt,2) = corr(FMdiff(outputNonNan),output(rr+1,outputNonNan,tt,3)');
        
        results.accuracy(rr,tt) = sum(output(rr+1,outputNonNan,tt,6) == clustETINonNan)./numel(clustETINonNan);
        for cc=1:3
            results.missrate(rr,tt,cc) = sum(output(rr+1,clusterEstimate==cc&outputNonNan,tt,6)~=cc)./sum(clustETINonNan==cc);
            results.ppv(rr,tt,cc) = sum(clustETINonNan(output(rr+1,outputNonNan,tt,6)==cc)==cc)./sum(output(rr+1,outputNonNan,tt,6)==cc);
        end;
    end;
    results.correlationETI(:,tt,1) = determineETI(results.correlation(:,tt,1),0.05);
    results.correlationETI(:,tt,2) = determineETI(results.correlation(:,tt,2),0.05);
    results.accuracyETI(:,tt) = determineETI(results.accuracy(:,tt),0.05);
    for cc=1:3
        results.missrateETI(:,tt,cc) = determineETI(results.missrate(:,tt,cc),0.05);
        results.ppvETI(:,tt,cc) = determineETI(results.ppv(:,tt,cc),0.05);
    end;
end;

save('Output/CrossValidationNumberOfGroupsIs10Init.mat','output','timePoints','fits','results')