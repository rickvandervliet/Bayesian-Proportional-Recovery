%% Coloring
c = linspecer(5);
colors=c;
darkcolors=brighten(c,-0.7);
lightcolors=brighten(c,0.25);

%% Figure 1. Mixture model results
figure(1)
set(gcf, 'color','white','Units','centimeters','Position', [26 2 23 13])

load('Data/CombinedTransform.mat')
load('Output/ResultsNumberOfGroupsIs5.mat')

t = [3/7 2 3 4 5 8 12 26];
groups = unique(mod.gETI);
numGroups = numel(groups);

x = -10:0.1:10;
xL = 66./(1+exp(-x));
FMi = NaN(numGroups,numel(x));
for nn=1:numGroups
    FMi(nn,:) = normpdf(x,mod.alphamETI(nn,1),mod.alphasETI(nn,1));
end

subplot(2,4,[1:2 5:6])
for nn=1:numGroups
    FMsub = FMt(:,mod.gETI==numGroups+1-nn);
    tsub = tt(:,mod.gETI==numGroups+1-nn);
    plot(tsub./7, FMsub, 'color', colors(numGroups+1-nn,:),'linewidth',0.1); hold on;
end;
for nn=1:numGroups
    alpha = 66./(1+exp(-mod.alphamETI(numGroups+1-nn,1)));
    FMmean = alpha + mod.rETI(numGroups+1-nn,1)*(66-alpha) * (1- exp(-(0:31)/mod.tauETI(numGroups+1-nn,1)));
    plot((0:31), FMmean,'linewidth', 3, 'color',darkcolors(numGroups+1-nn,:)); hold on;
    plot(31, FMmean(end),'marker', '>','color',darkcolors(numGroups+1-nn,:),'markerfacecolor',darkcolors(numGroups+1-nn,:));
    text(32,FMmean(end),sprintf('%0.2f',mod.rETI(numGroups+1-nn,1)),'fontweight','bold','fontsize',10,'backgroundcolor','white')
end;
for nn=1:numGroups
    alpha = 66./(1+exp(-mod.alphamETI(numGroups+1-nn,1)));
    plot(repmat(mod.tauETI(numGroups+1-nn,1),2,1), [-2 alpha+(66-alpha).*0.63.*mod.rETI(numGroups+1-nn)],'linewidth', 2, 'color',darkcolors(numGroups+1-nn,:));
    plot(mod.tauETI(numGroups+1-nn,1), alpha+(66-alpha).*0.63.*mod.rETI(numGroups+1-nn),'Marker', 'o', 'color',darkcolors(numGroups+1-nn,:), 'markerfacecolor',darkcolors(numGroups+1-nn,:));
    plot(mod.tauETI(numGroups+1-nn,1), -1,'Marker', 'v', 'color',darkcolors(numGroups+1-nn,:), 'markerfacecolor',darkcolors(numGroups+1-nn,:));
end;   
xlim([-1 35])
ylim([-2 66])
box off
xlabel('Time since stroke (w)','fontweight','bold','fontsize',10)
ylabel('FMA-UE','fontweight','bold','fontsize',10)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10)

subplot(2,4,3)
for nn=1:numGroups
    bar(nn,mod.gpETI(nn,1),'facecolor', colors(nn,:)); hold on;
end;
errorbar(1:numGroups,mod.gpETI(:,1),mod.gpETI(:,1)-mod.gpETI(:,2),mod.gpETI(:,3)-mod.gpETI(:,1),'linestyle','none','color','black')
xlabel('Groups (#)','fontweight','bold','fontsize',10)
ylabel('Probability','fontweight','bold','fontsize',10)
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',1:5)
title('p_k','fontweight','bold','fontsize',10)
xlim([0 6])

subplot(2,4,4)
for nn=1:numGroups
    bar(nn,mod.rETI(nn,1),'facecolor', colors(nn,:)); hold on;
end;   
errorbar(1:numGroups,mod.rETI(:,1),mod.rETI(:,1)-mod.rETI(:,2),mod.rETI(:,3)-mod.rETI(:,1),'linestyle','none','color','black')
xlabel('Groups (#)','fontweight','bold','fontsize',10)
ylabel('Recovery coefficient','fontweight','bold','fontsize',10)
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',1:5,'YTick',0:0.25:1)
title('r_k','fontweight','bold','fontsize',10)
xlim([0 6])

subplot(2,4,7)
for nn=1:numGroups
    bar(nn,mod.tauETI(nn,1),'facecolor', colors(nn,:)); hold on;
end;
errorbar(1:numGroups,mod.tauETI(:,1),mod.tauETI(:,1)-mod.tauETI(:,2),mod.tauETI(:,3)-mod.tauETI(:,1),'linestyle','none','color','black')
xlabel('Groups (#)','fontweight','bold','fontsize',10)
ylabel('Time since stroke (w)','fontweight','bold','fontsize',10)
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',1:5)
title('\tau_k','fontweight','bold','fontsize',10)
xlim([0 6])

subplot(2,4,8)
for nn=1:numGroups
    plot(xL,FMi(nn,:),'color',colors(nn,:),'linewidth',2); hold on;
end
xlim([0 66])
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',0:10:66,'YTick',0:0.1:0.6)
ylabel('Probability','fontweight','bold','fontsize',10)
title('Initial FMA-UE','fontweight','bold','fontsize',10)
xlabel('Initial FMA-UE','fontweight','bold','fontsize',10)
xlim([0 66])
ylim([0 0.65])

%% Figure 2
figure(2)
set(gcf, 'color','white','Units','centimeters','Position', [26 2 23 13])
load('Output/CrossValidationNumberOfGroupsIs5.mat')

sbp(1)=subplot(2,3,1);
crossValidationEndpoint = reshape(output(2:end,:,:,2),[],size(output,3));
boxplot(crossValidationEndpoint,1:12,'Symbol','','color','k')
ylim([0 35])
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12,'XTickLabel',2:2:12,'YTick',5:5:30)
set(findobj(sbp(1),'-regexp','Tag','\w*Whisker'),'LineStyle','-')
set(findobj(sbp(1),'-regexp','Tag','\w*Box'),'Color',colors(1,:))

sbp(2)=subplot(2,3,2);
crossValidationRecovery = reshape(output(2:end,:,:,4),[],size(output,3));
boxplot(crossValidationRecovery,'Symbol','','color','k')
ylim([0 35])
box off
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12,'XTickLabel',2:2:12,'YTick',5:5:30)
set(findobj(sbp(2),'-regexp','Tag','\w*Whisker'),'LineStyle','-')
set(findobj(sbp(2),'-regexp','Tag','\w*Box'),'Color',colors(2,:))

subplot(2,3,3)
errorbar(timePoints,results.correlationETI(1,:,1),results.correlationETI(1,:,1)-results.correlationETI(2,:,1),results.correlationETI(3,:,1)-results.correlationETI(1,:,1),'linewidth',1,'color',colors(1,:),'markerfacecolor','w','marker','o'); hold on;
errorbar(timePoints,results.correlationETI(1,:,2),results.correlationETI(1,:,2)-results.correlationETI(2,:,2),results.correlationETI(3,:,2)-results.correlationETI(1,:,2),'linewidth',1,'color',colors(2,:),'markerfacecolor','w','marker','v');
ylabel('Correlation','fontweight','bold','fontsize',10)
xlabel('Time since stroke (w)','fontweight','bold','fontsize',10)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12)
box off
xlim([0 12.5])
ylim([0.7 1.05])
legend('Endpoint','Recovery', 'Location', 'NorthWest')
legend boxoff

subplot(2,3,4)
errorbar(timePoints,results.accuracyETI(1,:),results.accuracyETI(1,:)-results.accuracyETI(2,:),results.accuracyETI(3,:)-results.accuracyETI(1,:),'linewidth',1,'color','k','markerfacecolor','w','marker','o'); hold on;
ylabel('Accuracy','fontweight','bold','fontsize',10)
xlabel('Time since stroke (w)','fontweight','bold','fontsize',10)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12)
box off
xlim([0 12.5])
ylim([0.75 0.95])
legend('Overall','Location', 'SouthEast');
legend boxoff

subplot(2,3,5)
errorbar(timePoints,results.ppvETI(1,:,1),results.ppvETI(1,:,1)-results.ppvETI(2,:,1),results.ppvETI(3,:,1)-results.ppvETI(1,:,1),'linewidth',1,'color',colors(3,:),'markerfacecolor','w','marker','o'); hold on;
errorbar(timePoints-0.1,results.ppvETI(1,:,2),results.ppvETI(1,:,2)-results.ppvETI(2,:,2),results.ppvETI(3,:,2)-results.ppvETI(1,:,2),'linewidth',1,'color',colors(4,:),'markerfacecolor','w','marker','v');
errorbar(timePoints+0.1,results.ppvETI(1,:,3),results.ppvETI(1,:,3)-results.ppvETI(2,:,3),results.ppvETI(3,:,3)-results.ppvETI(1,:,3),'linewidth',1,'color',colors(5,:),'markerfacecolor','w','marker','^');
ylabel('Positive Predictive Value','fontweight','bold','fontsize',10)
xlabel('Time since stroke (w)','fontweight','bold','fontsize',10)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12)
box off
xlim([0 12.5])
ylim([0.4 1.0])
legend('Poor', 'Moderate', 'Excellent','Location', 'SouthEast');
legend boxoff

subplot(2,3,6)
errorbar(timePoints,results.missrateETI(1,:,1),results.missrateETI(1,:,1)-results.missrateETI(2,:,1),results.missrateETI(3,:,1)-results.missrateETI(1,:,1),'linewidth',1,'color',colors(3,:),'markerfacecolor','w','marker','o'); hold on;
errorbar(timePoints-0.1,results.missrateETI(1,:,2),results.missrateETI(1,:,2)-results.missrateETI(2,:,2),results.missrateETI(3,:,2)-results.missrateETI(1,:,2),'linewidth',1,'color',colors(4,:),'markerfacecolor','w','marker','v');
errorbar(timePoints+0.1,results.missrateETI(1,:,3),results.missrateETI(1,:,3)-results.missrateETI(2,:,3),results.missrateETI(3,:,3)-results.missrateETI(1,:,3),'linewidth',1,'color',colors(5,:),'markerfacecolor','w','marker','^');
ylabel('Miss Rate','fontweight','bold','fontsize',10)
xlabel('Time since stroke (w)','fontweight','bold','fontsize',10)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10,'XTick',2:2:12)
box off
xlim([0 12.5])
ylim([0 0.8])

%% Tables
load('Data/CombinedTransform.mat')
load('Output/ResultsNumberOfGroupsIs5.mat')
clear table
table(1,1:5) = {'1','2','3','4','5'};
for ii=1:5
    table{2,ii} = sprintf('%0.2f [%0.2f %0.2f]',mod.gpETI(ii,1:3));
    table{3,ii} = sprintf('%0.2f [%0.2f %0.2f]',mod.rETI(ii,1:3));
    table{4,ii} = sprintf('%2.1f [%2.1f %2.1f]',mod.tauETI(ii,1:3));
    table{5,ii} = sprintf('%2.1f [%2.1f %2.1f]',mod.alphamETI(ii,1:3));
    table{6,ii} = sprintf('%2.1f [%2.1f %2.1f]',mod.alphasETI(ii,1:3));
end;
xlswrite('Paper/Table1.xlsx',table)
clear table

table(1,1:6) = {'Total','1','2','3','4','5'};
table{2,1} = sprintf('%2.1f',nanmean(dem.age));
table{3,1} = sprintf('%2.1f',100*nanmean(dem.gender));
table{4,1} = sprintf('%2.1f',100*nanmean(dem.preferred));
table{5,1} = sprintf('%2.1f',100*nanmean(dem.dominant));
table{6,1} = sprintf('%2.1f',100*nanmean(dem.bamford==1));
table{7,1} = sprintf('%2.1f',100*nanmean(dem.bamford==2));
table{8,1} = sprintf('%2.1f',100*nanmean(dem.bamford==3));
table{9,1} = sprintf('%2.1f',100*nanmean(dem.rtpa));
table{10,1} = sprintf('%2.1f',nanmean(dem.NIHSS));
table{11,1} = sprintf('%2.2f',nanmean(dem.SENS));
table{12,1} = sprintf('%2.2f',nanmean(dem.NEGL));
table{13,1} = sprintf('%2.2f',nanmean(dem.MI));
table{14,1} = sprintf('%2.2f',nanmean(dem.SA));
table{15,1} = sprintf('%2.2f',nanmean(dem.FE));

for ii=1:max(mod.gETI)
    table{2,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',sub.age(ii,:));
    table{3,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.gender(ii,:));
    table{4,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.preferred(ii,:));
    table{5,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.dominant(ii,:));
    table{6,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.bamford(ii,:,1));
    table{7,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.bamford(ii,:,2));
    table{8,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.bamford(ii,:,3));
    table{9,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',100*sub.rtpa(ii,:));
    table{10,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',sub.nihss(ii,:));
    table{11,ii+1} = sprintf('%2.2f [%2.2f %2.2f]',sub.sens(ii,:));
    table{12,ii+1} = sprintf('%2.2f [%2.2f %2.2f]',sub.negl(ii,:));
    table{13,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',sub.mi(ii,:));
    table{14,ii+1} = sprintf('%2.1f [%2.1f %2.1f]',sub.sa(ii,:));
    table{15,ii+1} = sprintf('%2.2f [%2.2f %2.2f]',sub.fe(ii,:));
end;
xlswrite('Paper/Table 2.xlsx',table)