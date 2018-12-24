function processProportionalRecoveryModel(samples,dem)

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

nsamples = size(g,1);
sub.age = NaN(numGroups,3);
ages = NaN(nsamples,1);
sub.rtpa = NaN(numGroups,3);
rtpas = NaN(nsamples,1);
sub.gender = NaN(numGroups,3);
genders = NaN(nsamples,1);
sub.dominant = NaN(numGroups,3);
dominants = NaN(nsamples,1);
sub.preferred = NaN(numGroups,3);
preferreds = NaN(nsamples,1);
sub.sa = NaN(numGroups,3);
sas = NaN(nsamples,1);
sub.fe = NaN(numGroups,3);
fes = NaN(nsamples,1);
sub.sens = NaN(numGroups,3);
senss = NaN(nsamples,1);
sub.negl = NaN(numGroups,3);
negls = NaN(nsamples,1);
sub.nihss = NaN(numGroups,3);
nihsss = NaN(nsamples,1);
sub.bamford = NaN(numGroups,3,3);
bamfords = NaN(nsamples,3);
sub.mi = NaN(numGroups,3);
mis = NaN(nsamples,1);

for ii=1:numGroups
    mod.rETI(ii,:) = determineETI(r(:,ii),0.05);
    mod.tauETI(ii,:) = determineETI(tau(:,ii),0.05);
    mod.alphamETI(ii,:) = determineETI(alpham(:,ii),0.05);
    mod.alphapETI(ii,:) = determineETI(alphap(:,ii),0.05);
    mod.alphasETI(ii,:) = determineETI(alphas(:,ii),0.05);
    mod.gpETI(ii,:) = determineETI(gp(:,ii),0.05);
    
    for s=1:nsamples
        ages(s) = nanmean(dem.age(g(s,:)==ii));
        rtpas(s) = nanmean(dem.rtpa(g(s,:)==ii));
        genders(s) = nanmean(dem.gender(g(s,:)==ii));
        dominants(s) = nanmean(dem.dominant(g(s,:)==ii));
        preferreds(s) = nanmean(dem.preferred(g(s,:)==ii));
        sas(s) = nanmean(dem.SA(g(s,:)==ii));
        fes(s) = nanmean(dem.FE(g(s,:)==ii));        
        senss(s) = nanmean(dem.SENS(g(s,:)==ii));
        negls(s) = nanmean(dem.NEGL(g(s,:)==ii));
        nihsss(s) = nanmean(dem.NIHSS(g(s,:)==ii));
        mis(s) = nanmean(dem.MI(g(s,:)==ii));
        for kk=1:3
            bamfords(s,kk) = nanmean(dem.bamford(g(s,:)==ii)==kk);
        end;
    end;
    sub.age(ii,:) = determineETI(ages,0.05);
    sub.rtpa(ii,:) = determineETI(rtpas,0.05);
    sub.gender(ii,:) = determineETI(genders,0.05);
    sub.dominant(ii,:) = determineETI(dominants,0.05);
    sub.preferred(ii,:) = determineETI(preferreds,0.05);
    sub.sa(ii,:) = determineETI(sas,0.05);
    sub.fe(ii,:) = determineETI(fes,0.05);
    sub.sens(ii,:) = determineETI(senss,0.05);
    sub.negl(ii,:) = determineETI(negls,0.05);
    sub.nihss(ii,:) = determineETI(nihsss,0.05);
    sub.mi(ii,:) = determineETI(mis,0.05);
    for kk=1:3
    	sub.bamford(ii,:,kk) = determineETI(bamfords(:,kk),0.05);
	end;
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

sub.age = sub.age(sortIndices,:);
sub.rtpa = sub.rtpa(sortIndices,:);
sub.gender = sub.gender(sortIndices,:);
sub.dominant = sub.dominant(sortIndices,:);
sub.preferred = sub.preferred(sortIndices,:);
sub.sa = sub.sa(sortIndices,:);
sub.fe = sub.fe(sortIndices,:);
sub.sens = sub.sens(sortIndices,:);
sub.negl = sub.negl(sortIndices,:);
sub.nihss = sub.nihss(sortIndices,:);
sub.bamford = sub.bamford(sortIndices,:,:);
sub.mi = sub.mi(sortIndices,:);

save(['Output/ResultsNumberOfGroupsIs' int2str(numGroups) '.mat'], 'mod','sub')