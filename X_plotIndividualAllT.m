clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\1_Temperature_TN_correctX.mat')

% vectors to set first and last recording per animal and the position of
% the same temperature conditions
first = [2 4 1 3 1 1 1 1 1 1];
last = [3 9 6 6 2 2 6 6 6 9];
width = [2 4 8];
pw = [1 1.5 2 2.5 3 3.5 4 4.5];

%% baseline correction
for i = 1 : size(AllAni,2)
    for j = first(i) : last(i)
        for w = 1 : 3
            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.(['w',num2str(width(w))]) = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.(['w',num2str(width(w))]) - AllAni(i).(['R0',num2str(j)]).background.sum(1);%rawmean.translation.fw.(['w',num2str(width(w))]);  
            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.(['w',num2str(width(w))]) = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.(['w',num2str(width(w))]) - AllAni(i).(['R0',num2str(j)]).background.sum(1);%rawmean.translation.bw.(['w',num2str(width(w))]);  
        end
    end
end

%% plot individual units for all temperature conditions and repetitions
figure
for i = 1 : size(AllAni,2)-1
    subplot(3,3,i)
    hold on
    c = 1;
    for j = first(i) : last(i)
        plot(AllAni(i).xvelo.w8,AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(1:length(AllAni(i).xvelo.w8)),'k-','LineWidth',pw(c))
        plot(AllAni(i).xvelo.w8,AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(1:length(AllAni(i).xvelo.w8)),'k-','LineWidth',pw(c))
        c = c + 1;
    end
    box on
    set(gca,'XScale','log','XTick',[10 100 1000],'XTickLabel',{'10','100','1000'})
    axis square
    plot([10^0 10^3.5],[0 0],'--','Color',[.5 .5 .5 .5],'LineWidth',2)
    xlim([12 3000])
    title(['A',num2str(AllAni(i).Animal),', U',num2str(AllAni(i).UnitNbr),', TN',num2str(AllAni(i).TNtype)])
end
set(gcf,'Position',[100 100 850 890])

print('Individual_TC','-dsvg','-r300','-painters')
savefig('Individual_TC.fig')