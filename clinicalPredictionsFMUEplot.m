function clinicalPredictionsFMUEplot(tt,FMt,ttpred,FMtpred,c,cp,FMtend)
    recoveryClusters = {'Poor','Moderate','Excellent'};
    colors = [0.6769,0.4447,0.7114;1.0000,0.5984,0.2000;0.4416,0.7490,0.4322;];
    
    subplot(1,2,1)
    plot(tt,FMt,'color','black','marker','o','linestyle','none','markerfacecolor','white'); hold on;
    shadedErrorBar(ttpred',FMtpred(:,1)',abs(FMtpred(:,[3 2])'-repmat(FMtpred(:,1)',2,1)),'lineProps',{'color',colors(c,:)}); 
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',10,'linewidth',2)
    box off
    xlabel('Time since stroke (weeks)','fontweight','bold','fontsize',10)
    ylabel('Fugl-Meyer upper extremity','fontweight','bold','fontsize',10)
    ylim([0 67])
    
    subplot(1,2,2)
    histogram(FMtend,'Normalization','probability','facecolor',colors(c,:),'edgecolor','none','numbins',50)
    xlabel('Fugl-Meyer upper extremity','fontweight','bold','fontsize',10)
    ylabel('Probability','fontweight','bold','fontsize',10)
    set(gca,'fontweight','bold','fontsize',10,'linewidth',2)
    box off
    xlim([-1 67])
    
	suptitle([recoveryClusters{c} ' recovery cluster (' int2str(cp) '%)'])
end

