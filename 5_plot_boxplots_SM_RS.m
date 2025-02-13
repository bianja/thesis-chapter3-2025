clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\3_Temperature_TN_correctX_median_BaselineCorrectionFit.mat')
plotval = [1 3 9 45 76 98 160 202 328 408 622 664 832 1000];
xval = AllAni(1).xvelo.w8;
x = xval(plotval);

SMC = [140 561 1686 449 499 449 899 449 899 NaN];
SMW = [140 561 1686 449 899 449 899 449 899 NaN];

for i = 1 : size(AllAni,2)
    if AllAni(i).TNtype == 1 % TN1, exc in btf
        ExcC(i,:) = AllAni(i).R01.yfreq.mean.translation.bw.w8(plotval);
        ExcW(i,:) = AllAni(i).R02.yfreq.mean.translation.bw.w8(plotval);
        ind(i) = 1;
    else % TN2, exc in ftb
        ExcC(i,:) = AllAni(i).R01.yfreq.mean.translation.fw.w8(plotval);
        ExcW(i,:) = AllAni(i).R02.yfreq.mean.translation.fw.w8(plotval);
        ind(i) = 2;
    end
end

for i = 1 : size(AllAni,2)
    AllAni(i).SMC = x(ExcC(i,:) == max(ExcC(i,:)));
    AllAni(i).SMW = x(ExcW(i,:) == max(ExcW(i,:)));
    AllAni(i).RSC = max(ExcC(i,:));
    AllAni(i).RSW = max(ExcW(i,:));
end

% Ani 161, correct for warm condition
AllAni(7).SMW = 899;
AllAni(7).RSW = AllAni(7).R02.yfreq.mean.translation.fw.w8(328);

%% boxplots
for i = 1 : size(AllAni,2)
    if AllAni(i).TNtype ~= 0
        SM(i,1:2) = [AllAni(i).SMC AllAni(i).SMW];
        RS(i,1:2) = [AllAni(i).RSC AllAni(i).RSW];
    end
end

% Sensitivity maximum
Boxplot_B(SM,2,15,[18/255 136/255 215/255; 170/255 32/255 55/255],{'RT','EV'},[1 2])
x = ones(size(SM,1),1).*(1+(rand(size(SM,1),1)*2-1)/10)+0;
scatter(x,SM(:,1),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(SM,1),1).*(1+(rand(size(SM,1),1)*2-1)/10)+1;
scatter(x,SM(:,2),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
set(gcf, 'Position',[875 255 375 370])
ylim([0 2000])
xlim([0 3])
xlabel('Head Temperature condition')
ylabel('sensitivity maximum (° s^-^1)')
print('SM_bp','-dsvg','-r300','-painters')
savefig('SM_bp.fig')

figure
hold on
for i = 1 : 9
    if AllAni(i).TNtype == 1
        plot(SM(i,:),'-c')
    else
        plot(SM(i,:),'-m')
    end
end
plot(mean(SM),'-k','LineWidth',2)
set(gcf, 'Position',[875 255 375 370])
xlim([0 3])
ylim([0 2000])
box on
set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
xlabel('Head Temperature condition')
ylabel('sensitivity maximum (° s^-^1)')
print('SM_lp','-dsvg','-r300','-painters')
savefig('SM_lp.fig')

% Response strength
Boxplot_B(RS,2,15,[18/255 136/255 215/255; 170/255 32/255 55/255],{'RT','EV'},[1 2])
x = ones(size(RS,1),1).*(1+(rand(size(RS,1),1)*2-1)/10)+0;
scatter(x,RS(:,1),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(RS,1),1).*(1+(rand(size(RS,1),1)*2-1)/10)+1;
scatter(x,RS(:,2),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
set(gcf, 'Position',[500 255 375 370])
ylim([0 350])
xlim([0 3])
xlabel('Head Temperature condition')
ylabel('response strength (spikes s^-^1)')
print('RS_bp','-dsvg','-r300','-painters')
savefig('RS_bp.fig')

figure
hold on
for i = 1 : 9
    if AllAni(i).TNtype == 1
        plot(RS(i,:),'-c')
    else
        plot(RS(i,:),'-m')
    end
end
plot(mean(RS),'-k','LineWidth',2)
set(gcf, 'Position',[500 255 375 370])
xlim([0 3])
ylim([0 350])
box on
set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
xlabel('Head Temperature condition')
ylabel('response strength (spikes s^-^1)')
print('RS_lp','-dsvg','-r300','-painters')
savefig('RS_lp.fig')

%% baseline
for i = 1 : size(AllAni,2)
    if AllAni(i).TNtype ~= 0
        bC(i) = AllAni(i).R01.background.sum(1);
        bW(i) = AllAni(i).R02.background.sum(1);
    end
end

figure
hold on
for i = 1 : 9
    if AllAni(i).TNtype == 1
        plot([bC(i) bW(i)],'-c')
    else
        plot([bC(i) bW(i)],'-m')
    end
end
plot([mean(bC) mean(bW)],'-k','LineWidth',2)
set(gcf, 'Position',[500 255 375 370])
xlim([0 3])
ylim([0 120])
box on
set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
xlabel('Head Temperature condition')
ylabel('baseline activity (spikes s^-^1)')
print('Baseline_lp','-dsvg','-r300','-painters')
savefig('Baseline_lp.fig')

% Q10 calculation
for i = 1 : size(SM,1)
    AllAni(i).Q10SM = (SM(i,2)/SM(i,1))^(10/(AllAni(1).R02.Temp-AllAni(1).R01.Temp));
    AllAni(i).Q10RS = (RS(i,2)/RS(i,1))^(10/(AllAni(1).R02.Temp-AllAni(1).R01.Temp));
    AllAni(i).Q10B = (bW(i)/bC(i))^(10/(AllAni(1).R02.Temp-AllAni(1).R01.Temp));
end

% Q10 plots
for i = 1 : size(AllAni,2)
    if AllAni(i).TNtype ~= 0
        Q10(i,1:3) = [AllAni(i).Q10SM AllAni(i).Q10RS AllAni(i).Q10B];
    end
end

Boxplot_B(Q10(:,1:2),2,15,[120/255 120/255 120/255; 120/255 120/255 120/255],{'SM','RS'},[1 2])
temp = find([AllAni.TNtype] == 1);
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+0;
scatter(x,Q10(temp,1),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+1;
scatter(x,Q10(temp,2),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
temp = find([AllAni.TNtype] == 2);
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+0;
scatter(x,Q10(temp,1),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+1;
scatter(x,Q10(temp,2),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
plot([0 5],[1 1],'--','Color',[.5 .5 .5 .5],'LineWidth',1.5)
set(gcf, 'Position',[500 255 375 370])
ylim([0 3])
xlim([0 4])
ylabel('Q10')
print('SM_Q10_bp','-dsvg','-r300','-painters')
savefig('SM_Q10_bp.fig')

Boxplot_B([Q10(:,3) Q10(:,3)],2,15,[120/255 120/255 120/255; 120/255 120/255 120/255],{'B','B'},[1 2])
temp = find([AllAni.TNtype] == 1);
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+0;
scatter(x,Q10(temp,3),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
temp = find([AllAni.TNtype] == 2);
x = ones(size(temp,2),1).*(1+(rand(size(temp,2),1)*2-1)/10)+0;
scatter(x,Q10(temp,3),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
plot([0 5],[1 1],'--','Color',[.5 .5 .5 .5],'LineWidth',1.5)
set(gcf, 'Position',[500 255 375 370])
ylim([0 25])
xlim([0 4])
ylabel('Q10')
print('Baseline_Q10_bp','-dsvg','-r300','-painters')
savefig('Baseline_Q10_bp.fig')


