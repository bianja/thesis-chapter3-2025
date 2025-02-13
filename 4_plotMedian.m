clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\3_Temperature_TN_correctX_median_BaselineCorrectionFit.mat')

% select TN type
for TN = 1 : 2
    if TN == 1
        Data = AllAni(find([AllAni.TNtype] == 1));
    elseif TN == 2
        Data = AllAni(find([AllAni.TNtype] == 2));
    else
        Data = AllAni;
    end

    % plot data (median)
    plotval = [1 3 9 45 76 98 160 202 328 408 622 664 832 1000];
    for k = 1 : 2
        if k == 1
            dir = 'fw';
        else
            dir = 'bw';
        end
        xval = Data(1).xvelo.w8;
        figure(2+TN)
        hold on
        for i = 1 : size(Data,2)
            maxval = max([Data(i).R01.yfreq.mean.translation.fw.w8' Data(i).R01.yfreq.mean.translation.bw.w8' Data(i).R02.yfreq.mean.translation.fw.w8' Data(i).R02.yfreq.mean.translation.bw.w8']);
            yvalCold(i,:) = Data(i).R01.yfreq.mean.translation.(char(dir)).w8/maxval;
            yvalWarm(i,:) = Data(i).R02.yfreq.mean.translation.(char(dir)).w8/maxval;
            subplot(2,3,i)
            plot(xval(plotval),yvalCold(i,plotval),'b-','LineWidth',k)
            hold on
            plot(xval(plotval),yvalWarm(i,plotval),'r-','LineWidth',k)
            plot([10^1 10^4],[0 0],'--','Color',[.5 .5 .5 .5])
            set(gca,'xscale','log','Box','on','XTick',[10 100 1000],'XTickLabel',{'10','100','1000'})
            axis square
            title(num2str(Data(i).Animal))
        end
        figure(0+TN)
        title(['median of TN',num2str(TN)])
        hold on       
        plot([10^1 10^4],[0 0],'--','Color',[.5 .5 .5 .5])
        plot_distribution_prctile(xval(plotval),yvalCold(:,plotval),'Color',[18/255 136/255 215/255],'LineWidth',1.5,'Alpha',0.15,'Prctile',[50])
        plot_distribution_prctile(xval(plotval),yvalWarm(:,plotval),'Color',[170/255 32/255 55/255],'LineWidth',1.5,'Alpha',0.15,'Prctile',[50])
         axis square
        set(gcf,'Position',[730 300 490 390])
        set(gca,'xscale','log','Box','on','XTick',[10 100 1000],'XTickLabel',{'10','100','1000'})
        xlim([10^1.32 10^3.45])
        xlim([xval(plotval(1)) xval(plotval(end))])
        ylim([-.6 1])
    end
end
figure(1)
print('Median_TC_TN1','-dsvg','-r300','-painters')
savefig('Median_TC_TN1.fig')
figure(2)
print('Median_TC_TN2','-dsvg','-r300','-painters')
savefig('Median_TC_TN2.fig')