function [results,figHandles] = externalValidation(id,tt,FMt,timePoints)
%externalValidation
%   id: patient id [integer]
%   tt: time in weeks 
%   FMt: Fugl-Meyer score
%   timePoints: time points in weeks for making predictions

    %Model parameters
    params.alpham = [-3.1602   -2.0922   -2.7862   -1.3081    0.0201];
    params.alphap = [3.8376    0.2340    0.1394    0.1326    0.1774];
    params.tau = [5.3114   10.1417    9.8122    2.6722    1.1822];
    params.r = [0.0884    0.4591    0.8634    0.8883    0.9309];
    params.gp = [0.2674    0.1395    0.1118    0.1783    0.3031];
    params.yp = 0.0673;
    params.clust = [1 2 2 3 3];
    
    %Parameter processing
    selectNonNan = ~isnan(FMt);
    id = id(selectNonNan);
    tt = tt(selectNonNan);
    FMt = FMt(selectNonNan);
    id = changem(id,1:numel(unique(id)),unique(id));
    [~, uniqueidLast] = unique(id,'last');
    [~, uniqueidFirst] = unique(id,'first');
    FMfinal = FMt(uniqueidLast);
    ttfinal = tt(uniqueidLast);
    FMinitial = FMt(uniqueidFirst);
    ttinitial = tt(uniqueidFirst);
    FMdiff = FMfinal - FMinitial;
    
    timePointsExt = [timePoints(:); max(tt)];
    numTime = numel(timePoints);
    output = NaN(numel(unique(id)),numTime,6);
    for t=1:(numTime+1)
        time = timePointsExt(t);
        select = tt<=time;
        dataset(1).FMtsub = FMt(select);
        dataset(1).ttsub = tt(select);
        dataset(1).idsub = id(select);
        dataset(1).ttinitialsub = ttinitial(unique(dataset(1).idsub));
        dataset(1).ttfinalsub = ttfinal(unique(dataset(1).idsub));
        dataset(1).FMinitialsub = FMinitial(unique(dataset(1).idsub));
        dataset(1).FMfinalsub = FMfinal(unique(dataset(1).idsub));
        dataset(1).orgidsub = unique(dataset(1).idsub);
        dataset(1).idsub = changem(dataset(1).idsub,1:numel(dataset(1).orgidsub),dataset(1).orgidsub);

        if ~isempty(dataset(1).FMtsub)
            outputTemp = predictionBugs(dataset(1),params);
            output(dataset(1).orgidsub,t,:) = permute(outputTemp,[2,3,1]);
        end
    end

    results.correlation = NaN(numTime,2);
    results.missrate = NaN(numTime,3);
    results.ppv = NaN(numTime,3);
    results.accuracy = NaN(numTime,1);

    for t=1:numTime
        outputNonNan = ~isnan(output(:,t,1));
        clustETINonNan = output(outputNonNan,3,6);
        results.correlation(t,1) = corr(FMfinal(outputNonNan),output(outputNonNan,t,1));
        results.correlation(t,2) = corr(FMdiff(outputNonNan),output(outputNonNan,t,3));

        results.accuracy(t) = sum(output(outputNonNan,t,6) == clustETINonNan)./numel(clustETINonNan);
        for cc=1:3
            results.missrate(t,cc) = sum(output(output(:,3,6)==cc&outputNonNan,t,6)~=cc)./sum(clustETINonNan==cc);
            results.ppv(t,cc) = sum(clustETINonNan(output(outputNonNan,t,6)==cc)==cc)./sum(output(outputNonNan,t,6)==cc);
        end;
    end;
    
    %% Plot figure
    c = linspecer(5);
    colors=c;
    figHandles.figure = figure(1);
    set(gcf, 'color','white','Units','centimeters','Position', [26    2   23   15])

    figHandles.subplot(1) = subplot(3,3,1);
    numPat = NaN(numel(timePoints),1);
    for ii=timePoints
        numPat(ii) = numel(unique(id(tt<=ii)));
    end;
    plot(timePoints,numPat,'linewidth',1,'color','k','markerfacecolor','w','marker','o')
    ylabel('Patients (#)','fontweight','bold','fontsize',10)
    box off
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    
    figHandles.subplot(2) = subplot(3,3,2);
    numMeas = NaN(numel(timePoints),3);
    for ii=timePoints
        numMeasTemp = sum(repmat(id(tt<=ii),1,max(id))==repmat(1:max(id),sum(tt<=ii),1));
        numMeasTemp(numMeasTemp==0) = NaN;
        numMeas(ii,:) = [nanmedian(numMeasTemp),prctile(numMeasTemp,5),prctile(numMeasTemp,95)];
    end;
    errorbar(timePoints,numMeas(:,1),numMeas(:,2)-numMeas(:,1),numMeas(:,1)-numMeas(:,3),'linewidth',1,'color','k','markerfacecolor','w','marker','o')
    ylabel('Measurements (#)','fontweight','bold','fontsize',10)
    box off
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    set(findobj(figHandles.subplot(2),'-regexp','Tag','\w*Whisker'),'LineStyle','-')
    set(findobj(figHandles.subplot(2),'-regexp','Tag','\w*Box'),'Color',colors(2,:))

    potRecs = NaN(numel(timePoints),max(id));
    figHandles.subplot(3) = subplot(3,3,3);
    for ii=timePoints
        for ss=1:max(id)
            ttsel = tt(id==ss);
            FMtsel = FMt(id==ss);
            potRecTemp = FMtsel(end)-FMtsel(find(ttsel<=ii,1,'last'));
            if ~isempty(potRecTemp)
                potRecs(ii,ss) = potRecTemp;
            end;
        end;
    end;
    boxplot(potRecs',timePoints,'Symbol','','color','k')
    box off
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    set(findobj(figHandles.subplot(3),'-regexp','Tag','\w*Whisker'),'LineStyle','-')

    figHandles.subplot(4) = subplot(3,3,4);
    boxplot(output(:,1:numTime,2),timePoints,'Symbol','','color','k')
    ylim([0 31])
    box off
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    set(findobj(figHandles.subplot(4),'-regexp','Tag','\w*Whisker'),'LineStyle','-')
    set(findobj(figHandles.subplot(4),'-regexp','Tag','\w*Box'),'Color',colors(1,:))

    figHandles.subplot(5) = subplot(3,3,5);
    boxplot(output(:,1:numTime,4),'Symbol','','color','k')
    ylim([0 31])
    box off
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    set(findobj(figHandles.subplot(5),'-regexp','Tag','\w*Whisker'),'LineStyle','-')
    set(findobj(figHandles.subplot(5),'-regexp','Tag','\w*Box'),'Color',colors(2,:))

    figHandles.subplot(6) = subplot(3,3,6);
    plot(timePoints,results.correlation(:,1),'linewidth',1,'color',colors(1,:),'markerfacecolor','w','marker','o'); hold on;
    plot(timePoints,results.correlation(:,2),'linewidth',1,'color',colors(2,:),'markerfacecolor','w','marker','v');
    ylabel('Correlation','fontweight','bold','fontsize',10)
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    box off

    figHandles.subplot(7) = subplot(3,3,7);
    plot(timePoints,results.accuracy,'linewidth',1,'color','k','markerfacecolor','w','marker','o'); hold on;
    ylabel('Accuracy','fontweight','bold','fontsize',10)
    xlabel('Time post stroke (w)','fontweight','bold','fontsize',10)
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    box off

    figHandles.subplot(8) = subplot(3,3,8);
    plot(timePoints,results.ppv(:,1),'linewidth',1,'color',colors(5,:),'markerfacecolor','w','marker','o'); hold on;
    plot(timePoints-0.1,results.ppv(:,2),'linewidth',1,'color',colors(4,:),'markerfacecolor','w','marker','v');
    plot(timePoints+0.1,results.ppv(:,3),'linewidth',1,'color',colors(3,:),'markerfacecolor','w','marker','^');
    ylabel('Positive Predictive Value','fontweight','bold','fontsize',10)
    xlabel('Time post stroke (w)','fontweight','bold','fontsize',10)
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    box off

    figHandles.subplot(9) = subplot(3,3,9);
    plot(timePoints,results.missrate(:,1),'linewidth',1,'color',colors(5,:),'markerfacecolor','w','marker','o'); hold on;
    plot(timePoints-0.1,results.missrate(:,2),'linewidth',1,'color',colors(4,:),'markerfacecolor','w','marker','v');
    plot(timePoints+0.1,results.missrate(:,3),'linewidth',1,'color',colors(3,:),'markerfacecolor','w','marker','^');
    ylabel('Miss Rate','fontweight','bold','fontsize',10)
    xlabel('Time post stroke (w)','fontweight','bold','fontsize',10)
    set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
    box off
end

