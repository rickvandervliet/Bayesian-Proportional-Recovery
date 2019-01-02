load('Data/CombinedTransform.mat')
load('Output/ResultsNumberOfGroupsIs5.mat')

timePoints = 1:12;
NPostSamples = 100;
numFolds = size(FMt,2);
numTime = numel(timePoints);
output = NaN(NPostSamples+1,size(FMt,2),numTime,6);

id = repmat(1:size(FMt,2),8,1);
randInd = randperm(size(FMt,2));
FMts = FMt(:,randInd);
ids = id(:,randInd);
tts = tt(:,randInd);
selectNonNan = ~isnan(FMts);

knnmdl = fitcknn([mod.rETI(:,1), mod.tauETI(:,1), mod.alphamETI(:,1), mod.alphapETI(:,1)],[1 2 2 3 3],'Standardize',1);
fits = struct('params', cell(1, numFolds));

for nn=1:numFolds
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
    end;

    fits(nn).params = fitBugs(set(2).FMt,set(2).tt,set(2).id,NPostSamples,knnmdl,mod);
    
    for t=1:numTime
        time = timePoints(t);
        selectEarly = set(1).tt<time.*7;
        FMtsub = set(1).FMt(selectEarly);
        ttsub = set(1).tt(selectEarly);
        idsub = set(1).id(selectEarly);
        ttinitialsub = set(1).ttinitial(unique(idsub));
        ttfinalsub = set(1).ttfinal(unique(idsub));
        FMinitialsub = set(1).FMinitial(unique(idsub));
        FMfinalsub = set(1).FMfinal(unique(idsub));
        idsub = changem(idsub,1:numel(unique(idsub)),unique(idsub));
        orgidsub = set(1).orgid(unique(idsub));
        
        for rr=1:NPostSamples+1
            if ~isempty(FMtsub)
                outputTemp = predictionBugs(FMtsub,ttsub,idsub,ttfinalsub,ttinitialsub,FMfinalsub,FMinitialsub,fits(nn).params(rr));
                output(rr,orgidsub,t,:) = permute(outputTemp,[4,2,3,1]);
            end;
        end;
    end;
end;

FMtnn = FMt(~isnan(FMt));
idnn = id(~isnan(FMt));
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
        outputNonNan = ~isnan(output(rr,:,tt,1));
        clustETINonNan = mod.clustETI(outputNonNan);
        results.correlation(rr,tt,1) = corr(FMfinal(outputNonNan),output(rr+1,outputNonNan,tt,1)');
        results.correlation(rr,tt,2) = corr(FMdiff(outputNonNan),output(rr+1,outputNonNan,tt,3)');
        
        results.accuracy(rr,tt) = sum(output(rr+1,outputNonNan,tt,6) == clustETINonNan)./numel(clustETINonNan);
        for cc=1:3
            results.missrate(rr,tt,cc) = sum(output(rr+1,(mod.clustETI==cc)&outputNonNan,tt,6)~=cc)./sum((mod.clustETI==cc)&outputNonNan);
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

save('Output/CrossValidationNumberOfGroupsIs5.mat','output','timePoints','fits','results')