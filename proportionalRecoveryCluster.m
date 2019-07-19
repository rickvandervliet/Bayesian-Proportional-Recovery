load('Data/CombinedTransform.mat')
load('Output/ResultsNumberOfGroupsIs10.mat')
numReps = 100;
overlap = ismember(dem(2).subs,dem(1).subs);
clusterEstimate(overlap) = mod.clustETI;
FMt = dem(2).FMt(:,~overlap);
tt = dem(2).tt(:,~overlap);

id = repmat(1:size(FMt,2),8,1);
selectNonNan = ~isnan(FMt);
FMt = FMt(selectNonNan);
tt = tt(selectNonNan);
id = id(selectNonNan);

load('Output/NumberOfGroupsIs10.mat')
randInd = randperm(5*10^4,numReps);

r = squeeze(samples.r(1,randInd,mod.grInc));
tau = squeeze(samples.tau(1,randInd,mod.grInc));
alpham = squeeze(samples.alpham(1,randInd,mod.grInc));
alphap = squeeze(samples.alphap(1,randInd,mod.grInc));
gp = squeeze(samples.gp(1,randInd,mod.grInc));
yp = squeeze(samples.yp(1,randInd));
params.clust = mod.clust;

clustersNonOverlap = NaN(numReps,sum(~overlap));
set.FMtsub = FMt;
set.ttsub = tt./7;
set.idsub = id;

for ii=1:numReps
    params.r = r(ii,:);
    params.tau = tau(ii,:);
    params.alpham = alpham(ii,:);
    params.alphap = alphap(ii,:);
    params.gp = gp(ii,:);
    params.yp = yp(ii);
    outputTemp = clusterBugs(set,params);
    clustersNonOverlap(ii,:) = outputTemp(2,:);
end;    
clusterEstimate(~overlap) = median(clustersNonOverlap);
save('Output/clusterEstimate.mat','clusterEstimate')