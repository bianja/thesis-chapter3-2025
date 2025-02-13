
close all; clear all;
load('Temperature_TN.mat')
dir = 'fw';
Data = AllAni;
AllAni = Data(8:9);
% x = [0.5 1 5 10 20 40 60 80 100 120];
% 
% x2 = [5.625 11.25 56.25 112.5 225 450 675 900 1125 1350]; x2 = x2 * 0.0889;
% x4 = [8.4375 16.875 84.375 168.75 337.5 675 1012.5 1350 1687.5 2025]; x4 = x4 * 0.0444;
% x8 = [14.0625 28.125 140.625 281.25 562.5 1125 1687.5 2250 2812.5 3375]; x8 = x8 * 0.0222;
% v2 = [5.625 11.25 56.25 112.5 225 450 675 900 1125 1350];
% v4 = [8.4375 16.875 84.375 168.75 337.5 675 1012.5 1350 1687.5 2025];
% v8 = [14.0625 28.125 140.625 281.25 562.5 1125 1687.5 2250 2812.5 3375];

%% TF - plot frequency - this is oanly for w8
% AllAni(1).xvelo.w8 = AllAni(1).xvelo.w8(1:7);
% AllAni(2).xvelo.w8 = AllAni(2).xvelo.w8(1:7);
close all
figure % tr, bw, freq
plotFreqT3I(AllAni,'translation','bw',2,1);
plotFreqT3M(AllAni,'translation','bw',2,1);
subplot(1,2,1)
ylim([-100 150])
subplot(1,2,2)
ylim([-20 60])
set(gcf,'Position',[711 336 560 213])
print('TN_163_U1U5_btf','-dsvg','-r300','-painters')
savefig('TN_163_U1U5_btf.fig')

figure % tr, fw, freq
plotFreqT3I(AllAni,'translation','fw',2,1);
plotFreqT3M(AllAni,'translation','fw',2,1);
subplot(1,2,1)
ylim([-100 150])
subplot(1,2,2)
ylim([-20 60])
set(gcf,'Position',[711 336 560 213])
print('TN_163_U1U5_ftb','-dsvg','-r300','-painters')
savefig('TN_163_U1U5_ftb.fig')

%% substract bg
% % AllAni = Data;
% for i = 1 : size(AllAni,2)
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.bw.w2 AllAni(i).R02.yfreq.mean.translation.bw.w2 AllAni(i).R03.yfreq.mean.translation.bw.w2]);
%     AllAni(i).R01.yfreq.mean.translation.bw.w2 = AllAni(i).R01.yfreq.mean.translation.bw.w2/maxval;
%     AllAni(i).R02.yfreq.mean.translation.bw.w2 = AllAni(i).R02.yfreq.mean.translation.bw.w2/maxval;
%     AllAni(i).R03.yfreq.mean.translation.bw.w2 = AllAni(i).R03.yfreq.mean.translation.bw.w2/maxval;
% 
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.bw.w4 AllAni(i).R02.yfreq.mean.translation.bw.w4 AllAni(i).R03.yfreq.mean.translation.bw.w4]);
%     AllAni(i).R01.yfreq.mean.translation.bw.w4 = AllAni(i).R01.yfreq.mean.translation.bw.w4/maxval;
%     AllAni(i).R02.yfreq.mean.translation.bw.w4 = AllAni(i).R02.yfreq.mean.translation.bw.w4/maxval;
%     AllAni(i).R03.yfreq.mean.translation.bw.w4 = AllAni(i).R03.yfreq.mean.translation.bw.w4/maxval;
% 
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.bw.w8 AllAni(i).R02.yfreq.mean.translation.bw.w8 AllAni(i).R03.yfreq.mean.translation.bw.w8]);
%     AllAni(i).R01.yfreq.mean.translation.bw.w8 = AllAni(i).R01.yfreq.mean.translation.bw.w8/maxval;
%     AllAni(i).R02.yfreq.mean.translation.bw.w8 = AllAni(i).R02.yfreq.mean.translation.bw.w8/maxval;
%     AllAni(i).R03.yfreq.mean.translation.bw.w8 = AllAni(i).R03.yfreq.mean.translation.bw.w8/maxval;
% 
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.fw.w2 AllAni(i).R02.yfreq.mean.translation.fw.w2 AllAni(i).R03.yfreq.mean.translation.fw.w2]);
%     AllAni(i).R01.yfreq.mean.translation.fw.w2 = AllAni(i).R01.yfreq.mean.translation.fw.w2/maxval;
%     AllAni(i).R02.yfreq.mean.translation.fw.w2 = AllAni(i).R02.yfreq.mean.translation.fw.w2/maxval;
%     AllAni(i).R03.yfreq.mean.translation.fw.w2 = AllAni(i).R03.yfreq.mean.translation.fw.w2/maxval;
% 
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.fw.w4 AllAni(i).R02.yfreq.mean.translation.fw.w4 AllAni(i).R03.yfreq.mean.translation.fw.w4]);
%     AllAni(i).R01.yfreq.mean.translation.fw.w4 = AllAni(i).R01.yfreq.mean.translation.fw.w4/maxval;
%     AllAni(i).R02.yfreq.mean.translation.fw.w4 = AllAni(i).R02.yfreq.mean.translation.fw.w4/maxval;
%     AllAni(i).R03.yfreq.mean.translation.fw.w4 = AllAni(i).R03.yfreq.mean.translation.fw.w4/maxval;
% 
%     maxval = max([AllAni(i).R01.yfreq.mean.translation.fw.w8 AllAni(i).R02.yfreq.mean.translation.fw.w8 AllAni(i).R03.yfreq.mean.translation.fw.w8]);
%     AllAni(i).R01.yfreq.mean.translation.fw.w8 = AllAni(i).R01.yfreq.mean.translation.fw.w8/maxval;
%     AllAni(i).R02.yfreq.mean.translation.fw.w8 = AllAni(i).R02.yfreq.mean.translation.fw.w8/maxval;
%     AllAni(i).R03.yfreq.mean.translation.fw.w8 = AllAni(i).R03.yfreq.mean.translation.fw.w8/maxval;
% end

%% mean of all ani
for i = 1 : size(AllAni, 2)
    maxval = max([AllAni(i).R01.yfreq.mean.translation.(char(dir)).w8 AllAni(i).R02.yfreq.mean.translation.(char(dir)).w8 AllAni(i).R03.yfreq.mean.translation.(char(dir)).w8]);
    T1(i,:) = AllAni(i).R01.yfreq.mean.translation.(char(dir)).w8/maxval;
    T2(i,:) = AllAni(i).R02.yfreq.mean.translation.(char(dir)).w8/maxval;
    T3(i,:) = AllAni(i).R03.yfreq.mean.translation.(char(dir)).w8/maxval;
end


x = [.5 1 5 10 20 40 60 80 100 120];
figure
for i = 1 : size(T1,1)
    
    subplot(1,3,1)
    plot(x,T1(i,:),'-','Color',[.7 .7 .7],'LineWidth',1.5)
    set(gca,'XScale','log','Box','on')
    xlabel('temporal frequency (Hz)')
    ylabel('norm spike rate')
    hold on
    plot(x,median(T1),'-','Color',[.3 .3 .3],'LineWidth',2.5)
    
    subplot(1,3,2)
    plot(x,T2(i,:),'-','Color',[.7 .7 .7],'LineWidth',1.5)
    set(gca,'XScale','log','Box','on')
    xlabel('temporal frequency (Hz)')
    ylabel('norm spike rate')
    hold on
    plot(x,median(T2),'-','Color',[.3 .3 .3],'LineWidth',2.5)
    
    subplot(1,3,3)
    plot(x,T3(i,:),'-','Color',[.7 .7 .7],'LineWidth',1.5)
    set(gca,'XScale','log','Box','on')
    xlabel('temporal frequency (Hz)')
    ylabel('norm spike rate')
    hold on
    plot(x,median(T3),'-','Color',[.3 .3 .3],'LineWidth',2.5)
end
addpath shaded_plots
figure
plot_distribution_prctile(x,T1,'Color',[5/255 48/255 97/255],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
hold on
plot_distribution_prctile(x,T2,'Color',[158/255 1/255 66/255],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
plot_distribution_prctile(x,T3,'Color',[69/255 117/255 180/255],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
set(gca,'XScale','log','Box','on')
xlabel('temporal frequency (Hz)')
ylabel('norm spike rate')

% print('median_firstData_ftb','-depsc','-r300','-tiff','-painters')
% savefig('median_firstData_ftb.fig')

%% gain
% ani = [1 2 3 4 5 9 10 11];
for i = 1 : size(AllAni, 2)
    % find max
    maxval = [max(AllAni(i).R01.yfreq.mean.translation.(char(dir)).w8) max(AllAni(i).R02.yfreq.mean.translation.(char(dir)).w8) max(AllAni(i).R03.yfreq.mean.translation.(char(dir)).w8)];
    T1(i,:) = AllAni(i).R01.yfreq.mean.translation.(char(dir)).w8/max(maxval);
    T2(i,:) = AllAni(i).R02.yfreq.mean.translation.(char(dir)).w8/max(maxval);
    T3(i,:) = AllAni(i).R03.yfreq.mean.translation.(char(dir)).w8/max(maxval);
    % for gain calculation
    temp = max(T1(i,:));
    g1(i) = mean(temp);
    temp = max(T2(i,:));
    g2(i) = mean(temp);
    temp = max(T3(i,:));
    g3(i) = mean(temp);
end

% gain
clear temp gainT1 gainT2 gainC table6
for i = 1 : size(AllAni, 2)
    gainT1(i) = g2(i)/g1(i);
    gainT2(i) = g2(i)/g3(i);
    gainC(i) = g3(i)/g1(i);
end

table6 = [gainT1' gainT2' gainC'];
Boxplot_B(table6,3,15,[.7 .3 .3; .5 .3 .3; .3 .3 .3],{'T01','T02','T01'},{'1','2','1'})
set(gcf,'Position',[1082, 375, 251, 325])
ylabel('gain')


for i = 1 : size(AllAni, 2)
    S1(i,:) = (T2(i,:))./(T1(i,:));
    S2(i,:) = (T2(i,:))./(T3(i,:));
    S3(i,:) = (T3(i,:))./(T1(i,:));
    S1s = std(T2)./std(T1);
    S2s = std(T2)./std(T3);
    S3s = std(T3)./std(T1);
    for j = 1 : 10
        if S1(i,j) == Inf
            S1(i,j) = NaN;
        end
        if S2(i,j) == Inf
            S2(i,j) = NaN;
        end
        if S3(i,j) == Inf
            S3(i,j) = NaN;
        end
    end
end

figure
hold on
plot_distribution_prctile(x,S1,'Color',[.7 .3 .3],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
plot_distribution_prctile(x,S2,'Color',[.5 .3 .3],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
plot_distribution_prctile(x,S3,'Color',[.3 .3 .3],'LineWidth',2.5,'Alpha',0.15,'Prctile',[25])
set(gca,'XScale','log','Box','on')
xlabel('temporal frequency (Hz)')
ylabel('gain')




