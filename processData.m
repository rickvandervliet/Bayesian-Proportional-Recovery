%% Combined dataset
combined = importdata('Data/Long format met NIHSS sens en neglect.xlsx');

% Time data
t = reshape(combined.data(:,16),13,[]);
t(t<0)=NaN;
FM = reshape(combined.data(:,13),13,[]);
FM(FM<0)=NaN;

MI = reshape(combined.data(:,10),13,[]);
MI(MI<0)=NaN;
[~,SAind]=max(~isnan(MI));
MI = MI(sub2ind(size(MI), SAind, 1:size(MI,2)));

SA = reshape(combined.data(:,9),13,[]);
SA(SA<0)=NaN;
[~,SAind]=max(~isnan(SA));
SA = SA(sub2ind(size(SA), SAind, 1:size(SA,2)));
SA(~isnan(SA)) = SA(~isnan(SA)) >= 9;

FE = reshape(combined.data(:,12),13,[]);
FE(FE<0)=NaN;
[~,FEind]=max(~isnan(FE));
FE = FE(sub2ind(size(FE), FEind, 1:size(FE,2)));
FE(~isnan(FE)) = FE(~isnan(FE)) >= 1;

SENS = reshape(combined.data(:,17),13,[]);
SENS(SENS<0)=NaN;
[~,SENSind]=max(~isnan(SENS));
SENS = SENS(sub2ind(size(SENS), SENSind, 1:size(SENS,2)));

NEGL = reshape(combined.data(:,18),13,[]);
NEGL(NEGL<0)=NaN;
[~,NEGLind]=max(~isnan(NEGL));
NEGL = NEGL(sub2ind(size(NEGL), NEGLind, 1:size(NEGL,2)));

NIHSS = reshape(combined.data(:,15),13,[]);
NIHSS(NIHSS<0)=NaN;
[~,NIHSSind]=max(~isnan(NIHSS));
NIHSS = NIHSS(sub2ind(size(NIHSS), NIHSSind, 1:size(NIHSS,2)));

% Predictors
age = reshape(combined.data(:,2),13,[]);
gender = reshape(combined.data(:,3),13,[]);
bamford = reshape(combined.data(:,4),13,[]);
rtpa = reshape(combined.data(:,5),13,[]);
affected = reshape(combined.data(:,6),13,[]);
preferred = reshape(combined.data(:,7),13,[]);
cirstot = reshape(combined.data(:,8),13,[]);

tt = NaN(8,479);
FMt = NaN(8,479);
deltaFM = NaN(1,479);
deltat = NaN(1,479);
numMeas = NaN(1,479);
initialt = NaN(1,479);

for ii=1:479
    ind = find(~isnan(FM(:,ii)) & ~isnan(t(:,ii)));
    FMt(1:numel(ind),ii) = FM(ind,ii);
    tt(1:numel(ind),ii) = t(ind,ii);
    
    numMeas(ii) =  numel(ind);
    
    if numel(ind)>1
        deltaFM(ii) = min(repmat(FM(ind(end),ii),size(FM(ind,ii)))-FM(ind,ii));
        deltat(ii) = (max(tt(:,ii))-min(tt(:,ii)));
    end
end
dem(1).subs = find(deltat>=12*7 & numMeas>=2);

gender = changem(gender,[1 0],[1 2]);
rtpa(rtpa>1) = NaN;
affected = changem(affected,[1 0],[1 2]);
preferred = changem(preferred,[1 0],[1 2]);
preferred(preferred==3) = NaN;
dominant = +(preferred==affected);
dominant(isnan(preferred) | isnan(affected)) = NaN;

dem(1).FMt = FMt(:,dem(1).subs);
dem(1).tt = tt(:,dem(1).subs);
dem(1).age = age(1,dem(1).subs);
dem(1).gender = gender(1,dem(1).subs);
dem(1).bamford = bamford(1,dem(1).subs);
dem(1).rtpa = rtpa(1,dem(1).subs);
dem(1).affected = affected(1,dem(1).subs);
dem(1).preferred = preferred(1,dem(1).subs);
dem(1).dominant = dominant(1,dem(1).subs);
dem(1).cirstot = cirstot(1,dem(1).subs);
dem(1).SA = SA(1,dem(1).subs);
dem(1).MI = MI(1,dem(1).subs);
dem(1).FE = FE(1,dem(1).subs);
dem(1).SENS = SENS(1,dem(1).subs);
dem(1).NEGL = NEGL(1,dem(1).subs);
dem(1).NIHSS = NIHSS(1,dem(1).subs);

save('Data/CombinedTransform.mat','dem')