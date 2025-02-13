function [] = plotFreqT3M(Data, type, dir, d, b)
% data is struct with all units to analyse(e.g. data.unit01; data.unit02; ..)

% calculate number of units to plot and set numbers of subplots
% d = 6;
units = size(struct2table(Data),1);
% if rem(units,d) == 0 % even 
%     b = units/d;
% else % uneven
%     b = (units+d-1)/d;
% end
% b = 6;
% figure
for k = 1 : units
    subplot(b,d,k)
    hold on
%     T1
    meanT1 = median([Data(k).R01.yfreq.mean.(char(type)).(char(dir)).w8'; Data(k).R02.yfreq.mean.(char(type)).(char(dir)).w8'; Data(k).R03.yfreq.mean.(char(type)).(char(dir)).w8']);
    meanT2 = median([Data(k).R04.yfreq.mean.(char(type)).(char(dir)).w8'; Data(k).R05.yfreq.mean.(char(type)).(char(dir)).w8'; Data(k).R06.yfreq.mean.(char(type)).(char(dir)).w8']);
    % meanT2 = median([Data(k).R07.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R08.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R09.yfreq.mean.(char(type)).(char(dir)).w8]);

    plot(Data(k).xvelo.w8, meanT1,'-', 'Color', [158/255 1/255 66/255], 'LineWidth', 2)
    plot(Data(k).xvelo.w8, meanT2,'-', 'Color', [69/255 117/255 180/255], 'LineWidth', 2)
    % plot(Data(k).xvelo.w8, meanT3,'-', 'Color', [255/255 156/255 0 0.5], 'LineWidth', 2)
    
    % bg = [Data(k).R01.background.sum(1) Data(k).R02.background.sum(1) Data(k).R03.background.sum(1) Data(k).R04.background.sum(1) Data(k).R05.background.sum(1) Data(k).R06.background.sum(1)];
    % bg1 = [Data(k).R01.background.sum(1) Data(k).R02.background.sum(1) Data(k).R03.background.sum(1)];
    % bg2 = [Data(k).R04.background.sum(1) Data(k).R05.background.sum(1) Data(k).R06.background.sum(1)];
    % bg2 = [Data(k).R07.background.sum(1) Data(k).R08.background.sum(1) Data(k).R09.background.sum(1)];
    set(gca,'XScale','log')
    lim = get(gca,'xlim');
    set(gca,'xtick',[10^0 10^1 10^2 10^3],'xticklabels',{'1','10','100','1000'},'Box','on')
    line([10^-2 10^4],[0 0],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5 0.5])
%     line([10^-2 10^3],[median(bg) median(bg)],'LineStyle','--','LineWidth',1,'Color',[0.5 0.5 0.5])  
    % line([10^-2 10^4],[median(bg1) median(bg1)],'LineStyle','--','LineWidth',1,'Color',[158/255 1/255 66/255 0.5])
    % line([10^-2 10^4],[median(bg2) median(bg2)],'LineStyle','--','LineWidth',1,'Color',[69/255 117/255 180/255 0.5])
    % line([10^-2 10^3],[median(bg3) median(bg3)],'LineStyle','--','LineWidth',1,'Color',[255/255 156/255 0 0.5])
    xlim([10^0.9 10^3.5])
end




