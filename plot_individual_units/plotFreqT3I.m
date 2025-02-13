function [] = plotFreqT3I(Data, type, dir, d, b)
% data is struct with all units to analyse(e.g. data.unit01; data.unit02; ..)

units = size(struct2table(Data),1);
for k = 1 : units
    subplot(b,d,k)
    hold on
    meanT1 = mean([Data(k).R01.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R02.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R03.yfreq.mean.(char(type)).(char(dir)).w8]);
    meanT2 = mean([Data(k).R04.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R05.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R06.yfreq.mean.(char(type)).(char(dir)).w8]);
    % meanT3 = median([Data(k).R07.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R08.yfreq.mean.(char(type)).(char(dir)).w8; Data(k).R09.yfreq.mean.(char(type)).(char(dir)).w8]);
    
    
    plot(Data(k).xvelo.w8, Data(k).R01.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)
    plot(Data(k).xvelo.w8, Data(k).R02.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)
    plot(Data(k).xvelo.w8, Data(k).R03.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)

    plot(Data(k).xvelo.w8, Data(k).R04.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [69/255 117/255 180/255 0.5], 'LineWidth', 1.5)
    plot(Data(k).xvelo.w8, Data(k).R05.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [69/255 117/255 180/255 0.5], 'LineWidth', 1.5)
    plot(Data(k).xvelo.w8, Data(k).R06.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [69/255 117/255 180/255 0.5], 'LineWidth', 1.5)

    % plot(Data(k).xvelo.w8, Data(k).R07.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [255/255 156/255 0 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R08.yfreq.mean.(char(type)).(char(dir)).w8,'--', 'Color', [255/255 156/255 0 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R09.yfreq.mean.(char(type)).(char(dir)).w8,':', 'Color', [255/255 156/255 0 0.5], 'LineWidth', 1.5)
    % 
    % bg1 = [Data(k).R01.background.sum(1) Data(k).R02.background.sum(1) Data(k).R03.background.sum(1)];
    % bg2 = [Data(k).R04.background.sum(1) Data(k).R05.background.sum(1) Data(k).R06.background.sum(1)];
    % bg3 = [Data(k).R07.background.sum(1) Data(k).R08.background.sum(1) Data(k).R09.background.sum(1)];


    % plot(Data(k).xvelo.w8, Data(k).R01.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R02.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [69/255 117/255 180/255 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R03.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)

    % plot(Data(k).xvelo.w8, Data(k).R04.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [154/255 1/255 172/255 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R05.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [103/255 174/255 223/255 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R06.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [69/255 117/255 180/255 0.5], 'LineWidth', 1.5)
    % plot(Data(k).xvelo.w8, Data(k).R07.yfreq.mean.(char(type)).(char(dir)).w8,'-', 'Color', [158/255 1/255 66/255 0.5], 'LineWidth', 1.5)

    % bg1 = [Data(k).R01.background.sum(1) Data(k).R01.background.sum(1) Data(k).R01.background.sum(1)];
    % bg2 = [Data(k).R02.background.sum(1) Data(k).R02.background.sum(1) Data(k).R02.background.sum(1)];
    % bg3 = [Data(k).R03.background.sum(1) Data(k).R03.background.sum(1) Data(k).R03.background.sum(1)]; 
    % bg4 = [Data(k).R04.background.sum(1) Data(k).R04.background.sum(1) Data(k).R04.background.sum(1)]; 
    % bg5 = [Data(k).R05.background.sum(1) Data(k).R05.background.sum(1) Data(k).R05.background.sum(1)]; 
    % bg6 = [Data(k).R06.background.sum(1) Data(k).R06.background.sum(1) Data(k).R06.background.sum(1)];    
    
    
    set(gca,'XScale','log')
    lim = get(gca,'xlim');
    set(gca,'xlim',[10^1 10^3.7],'xtick',[10^-1 10^0 10^1 10^2],'xticklabels',{'0.1','1','10','100'},'Box','on')
    % line([10^-2 10^3],[median(bg1) median(bg1)],'LineStyle','--','LineWidth',1,'Color',[158/255 1/255 66/255 0.5])
    % line([10^-2 10^3],[median(bg2) median(bg2)],'LineStyle','--','LineWidth',1,'Color',[69/255 117/255 180/255 0.5])
    % line([10^-2 10^3],[median(bg3) median(bg3)],'LineStyle','--','LineWidth',1,'Color',[158/255 1/255 66/255 0.5])
    % line([10^-2 10^3],[median(bg4) median(bg4)],'LineStyle','--','LineWidth',1,'Color',[154/255 1/255 172/255 0.5])
    % line([10^-2 10^3],[median(bg5) median(bg5)],'LineStyle','--','LineWidth',1,'Color',[103/255 174/255 223/255 0.5])
    % line([10^-2 10^3],[median(bg6) median(bg6)],'LineStyle','--','LineWidth',1,'Color',[69/255 117/255 180/255 0.5])
end




