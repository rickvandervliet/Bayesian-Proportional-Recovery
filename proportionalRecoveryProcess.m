function proportionalRecoveryProcess(samples,dem,cutoff)

orgNumGroups = size(samples.r,3);
gp = squeeze(samples.gp);
grInc = find(mean(gp)>=cutoff);
r = squeeze(samples.r(:,:,grInc));
[~,sortInd]=sort(mean(r));
grInc = grInc(sortInd);
r = r(:,sortInd);

gp = gp(:,grInc);
gp = gp./repmat(sum(gp,2),1,size(gp,2));
yp = squeeze(samples.yp);

tau = squeeze(samples.tau(:,:,grInc));
alpham = squeeze(samples.alpham(:,:,grInc));
alphap = squeeze(samples.alphap(:,:,grInc));
alphas = 1./sqrt(alphap);
alpha = squeeze(samples.alpha);
g = squeeze(samples.g);
g = changem(g,1:5,grInc);
numSubs = size(squeeze(samples.alpha),2);

mod.numGroups = numel(grInc);
mod.grInc = grInc;
mod.rETI = NaN(mod.numGroups,3);
mod.tauETI = NaN(mod.numGroups,3);
mod.alphamETI = NaN(mod.numGroups,3);
mod.alphapETI = NaN(mod.numGroups,3);
mod.alphasETI = NaN(mod.numGroups,3);
mod.gpETI = NaN(mod.numGroups,3);
mod.alphaETI = NaN(numSubs,3);
mod.gETI = mode(g);

for ii=1:numSubs
    mod.alphaETI(ii,:) = determineETI(alpha(:,ii),0.05);
end;

sub.age = NaN(mod.numGroups,3);
sub.rtpa = NaN(mod.numGroups,1);
sub.gender = NaN(mod.numGroups,1);
sub.dominant = NaN(mod.numGroups,1);
sub.preferred = NaN(mod.numGroups,1);
sub.sa = NaN(mod.numGroups,3);
sub.fe = NaN(mod.numGroups,3);
sub.sens = NaN(mod.numGroups,3);
sub.negl = NaN(mod.numGroups,3);
sub.nihss = NaN(mod.numGroups,3);
sub.bamford = NaN(mod.numGroups,3);
sub.mi = NaN(mod.numGroups,3);
sub.numSub = NaN(mod.numGroups,3);
age = repmat(dem.age,size(g,1),1);
gender = repmat(dem.gender,size(g,1),1);
dominant = repmat(dem.dominant,size(g,1),1);
preferred = repmat(dem.preferred,size(g,1),1);
sa = repmat(dem.SA,size(g,1),1);
fe = repmat(dem.FE,size(g,1),1);
sens = repmat(dem.SENS,size(g,1),1);
negl = repmat(dem.NEGL,size(g,1),1);
nihss = repmat(dem.NIHSS,size(g,1),1);
mi = repmat(dem.MI,size(g,1),1);
bamford = repmat(dem.bamford,size(g,1),1);

for ii=1:mod.numGroups
    mod.rETI(ii,:) = determineETI(r(:,ii),0.05);
    mod.tauETI(ii,:) = determineETI(tau(:,ii),0.05);
    mod.alphamETI(ii,:) = determineETI(alpham(:,ii),0.05);
    mod.alphapETI(ii,:) = determineETI(alphap(:,ii),0.05);
    mod.alphasETI(ii,:) = determineETI(alphas(:,ii),0.05);
    mod.gpETI(ii,:) = determineETI(gp(:,ii),0.05);
    sub.numSub(ii,:) = determineETI(sum(g==ii,2),0.05);
    sub.age(ii,:) = determineETI(age(g==ii),0.05);
    rtpa = repmat(dem.rtpa,size(g,1),1);
    rtpa(g~=ii) = NaN;
    
    sub.rtpa(ii,:) = sum(rtpa(g==ii)==1)./(numel(~isnan(rtpa(g==ii))));
    sub.gender(ii,:) = sum(gender(g==ii)==1)./(numel(~isnan(gender(g==ii))));
    sub.dominant(ii,:) = sum(dominant(g==ii)==1)./(numel(~isnan(dominant(g==ii))));
    sub.preferred(ii,:) = sum(preferred(g==ii)==1)./(numel(~isnan(preferred(g==ii))));
    sub.sa(ii,:) = sum(sa(g==ii)>0)./(numel(~isnan(sa(g==ii))));
    sub.fe(ii,:) = sum(fe(g==ii)>0)./(numel(~isnan(fe(g==ii))));
    sub.sens(ii,:) = [sum(sens(g==ii)==0),sum(sens(g==ii)==1),sum(sens(g==ii)==2)]./(numel(~isnan(sens(g==ii))));
    sub.negl(ii,:) = [sum(negl(g==ii)==0),sum(negl(g==ii)==1),sum(negl(g==ii)==2)]./(numel(~isnan(negl(g==ii))));
    sub.nihss(ii,:) = determineETI(nihss(g==ii),0.05);
    sub.mi(ii,:) = determineETI(mi(g==ii),0.05);
    sub.bamford(ii,:) = [sum(bamford(g==ii)==1),sum(bamford(g==ii)==2),sum(bamford(g==ii)==3)]./(numel(~isnan(bamford(g==ii))));
    
end;

mod.ys  = determineETI(1./sqrt(yp),0.05);
mod.yp  = determineETI(yp,0.05);

if mod.numGroups==5
    mod.clust = [1 2 2 3 3];
    mod.clustETI = mode(changem(g,mod.clust,1:5));
end;    

save(['Output/ResultsNumberOfGroupsIs' int2str(orgNumGroups) '.mat'], 'mod','sub')