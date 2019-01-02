function processProportionalRecoveryModel(samples)

numGroups = size(samples.r,3);
yp = squeeze(samples.yp);
r = squeeze(samples.r);
tau = squeeze(samples.tau);
alpham = squeeze(samples.alpham);
alphap = squeeze(samples.alphap);
alphas = 1./sqrt(alphap);
alpha = squeeze(samples.alpha);
gp = squeeze(samples.gp);
g = squeeze(samples.g);
numSubs = size(squeeze(samples.alpha),2);

mod.rETI = NaN(numGroups,3);
mod.tauETI = NaN(numGroups,3);
mod.alphamETI = NaN(numGroups,3);
mod.alphapETI = NaN(numGroups,3);
mod.alphasETI = NaN(numGroups,3);
mod.gpETI = NaN(numGroups,3);
mod.alphaETI = NaN(numSubs,3);
mod.nonEmptygroups = nonEmptyClasses(g,0.05);
mod.gETI = mode(g);

for ii=1:numSubs
    mod.alphaETI(ii,:) = determineETI(alpha(:,ii),0.05);
end;

for ii=1:numGroups
    mod.rETI(ii,:) = determineETI(r(:,ii),0.05);
    mod.tauETI(ii,:) = determineETI(tau(:,ii),0.05);
    mod.alphamETI(ii,:) = determineETI(alpham(:,ii),0.05);
    mod.alphapETI(ii,:) = determineETI(alphap(:,ii),0.05);
    mod.alphasETI(ii,:) = determineETI(alphas(:,ii),0.05);
    mod.gpETI(ii,:) = determineETI(gp(:,ii),0.05);
end;

[~,sortIndices]=sort(mod.rETI(:,1));
mod.rETI = mod.rETI(sortIndices,:);
mod.tauETI = mod.tauETI(sortIndices,:);
mod.alphamETI = mod.alphamETI(sortIndices,:);
mod.alphapETI = mod.alphapETI(sortIndices,:);
mod.alphasETI = mod.alphasETI(sortIndices,:);
mod.gpETI = mod.gpETI(sortIndices,:);
mod.gETI = changem(mod.gETI,1:numGroups,sortIndices);
mod.ys  = determineETI(1./sqrt(yp),0.05);
mod.yp  = determineETI(yp,0.05);

if numGroups==5
    mod.clustETI = mode(changem(g,[1 2 2 3 3],sortIndices));
end;    

save(['Output/ResultsNumberOfGroupsIs' int2str(numGroups) '.mat'], 'mod')