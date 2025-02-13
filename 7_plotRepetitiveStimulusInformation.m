clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\4_Temperature_TN_correctX_median_BaselineCorrectionFit_RepetitiveStim.mat')
Data = AllAni;

% TN = find([Data(2:end-1).TNtype] == 1);
% AllAni = Data(TN);

%% baseline
Data = AllAni;
fv = 4;
a = [1 1 1 1 1 1 1];
b = [1 1 1 1 1 1 1];
% b = [6 1 1 6 2 2 2];
a = [2 1 1 2 2 2 2];
b = [2 1 1 2 2 2 2];

% a = a(TN);
% b = b(TN);

close all
for k = fv : size(Data,2)-1
    tempT = []; tempS = []; tempRSbtf = []; tempRSbtfT = []; tempRSftb = []; tempRSftbT = [];
    for i = a(k-fv+1) : b(k-fv+1) 
        tempT = [tempT Data(k).(char(['C0',num2str(i)])).bg.T(1:end-2)];
        tempS = [tempS Data(k).(char(['C0',num2str(i)])).bg.spikes(1:end-2)];
        tempRSbtf = [tempRSbtf Data(k).(char(['C0',num2str(i)])).btf.spikes];
        tempRSbtfT = [tempRSbtfT Data(k).(char(['C0',num2str(i)])).btf.T];
        tempRSftb = [tempRSftb Data(k).(char(['C0',num2str(i)])).ftb.spikes];
        tempRSftbT = [tempRSftbT Data(k).(char(['C0',num2str(i)])).ftb.T];
    end
    
    clear ttempS ttempT
    c = 1;
    for i = 1 : 2 : length(tempT)
        ttempT(c) = (tempT(i)+tempT(i+1))/2;
        ttempS(c) = (tempS(i)+tempS(i+1))/2;
        c = c + 1;
    end
    clear tempT tempS
    tempT = ttempT;
    tempS = ttempS;

    table2(k,1) = mean(tempS(find(tempT == min(tempT)))); % min bg
    table2(k,2) = mean(tempS(find(tempT == max(tempT)))); % max bg

    figure(80)
    subplot(3,3,k-fv+1)
    hold on
    axis square
    box on
    plot(tempT,tempS,'.k')
    [P,S] = polyfit(tempT,tempS,1);
    rsq(k) = S.rsquared;
    yfit = P(1)*tempT+P(2);
    slope(k) = P(1);
    ycut(k) = P(2);
    plot(tempT,yfit,'r-.','LineWidth',1,'Color','k')
    xlim([22 34])

    % Steigungsgerade für ftb und btf stimuli
    [Pftb,Sftb] = polyfit(tempRSftbT,tempRSftb,1);
    B26ftb(k) = Pftb(1)*26+Pftb(2);
    B28ftb(k) = Pftb(1)*28+Pftb(2);
    B30ftb(k) = Pftb(1)*30+Pftb(2);
    B32ftb(k) = Pftb(1)*32+Pftb(2);
    B34ftb(k) = Pftb(1)*34+Pftb(2);
    [Pbtf,Sbtf] = polyfit(tempRSbtfT,tempRSbtf,1);
    B26btf(k) = Pbtf(1)*26+Pbtf(2);
    B28btf(k) = Pbtf(1)*28+Pbtf(2);
    B30btf(k) = Pbtf(1)*30+Pbtf(2);
    B32btf(k) = Pbtf(1)*32+Pbtf(2);
    B34btf(k) = Pbtf(1)*34+Pbtf(2);
    slopeftb(k) = Pftb(1);
    slopebtf(k) = Pbtf(1);
    rsqftb(k) = Sftb.rsquared;
    rsqbtf(k) = Sbtf.rsquared;
    
    % plot die Steigungsgeraden für ftb und btf
    yfitftb = Pftb(1)*tempRSftbT+Pftb(2);
    yfitbtf = Pbtf(1)*tempRSbtfT+Pbtf(2);
    if Data(k).TNtype == 1
        plot(tempRSbtfT,tempRSbtf,'.m')
        plot(tempRSbtfT,yfitbtf,'r-.','LineWidth',1,'Color','m')
    elseif Data(k).TNtype == 2
        plot(tempRSftbT,tempRSftb,'.c')
        plot(tempRSftbT,yfitftb,'r-.','LineWidth',1,'Color','c')
    end
    if k == 23
        legend('baseline','','ftb','','btf')
    end

    % figure(85)
    % hold on
    % axis square
    % box on
    % plot(tempT,yfit,'-','LineWidth',1.5,'Color',[.3 .3 .3])
    tempc = corrcoef(tempT,tempS);
    corres(k) = tempc(1,2);
    corrp(k) = tempc(1,1);
    B26(k) = P(1)*26+P(2);
    B28(k) = P(1)*28+P(2);
    B30(k) = P(1)*30+P(2);
    B32(k) = P(1)*32+P(2);
    B34(k) = P(1)*34+P(2);
    % ylim([0 300])
    % xlim([22 34])
end

% repetitive
for k = fv : size(Data,2)-1
    tempTf = []; tempSf = []; tempTb = []; tempSb = [];
    for i = a(k-fv+1) : b(k-fv+1) % Data(k).first : Data(k).last-1
        tempTf = [tempTf Data(k).(char(['C0',num2str(i)])).ftb.T];
        tempSf = [tempSf Data(k).(char(['C0',num2str(i)])).ftb.spikes];
        tempTb = [tempTb Data(k).(char(['C0',num2str(i)])).btf.T];
        tempSb = [tempSb Data(k).(char(['C0',num2str(i)])).btf.spikes];
    end
    P = polyfit(tempTf,tempSf,1);
    yfitf = P(1)*tempTf+P(2);
    P = polyfit(tempTb,tempSb,1);
    yfitb = P(1)*tempTb+P(2);

    table1(k,1) = mean(tempSf(find(tempTf == min(tempTf)))); % minftb
    table1(k,2) = mean(tempSf(find(tempTf == max(tempTf)))); % maxftb
    table1(k,3) = mean(tempSb(find(tempTb == min(tempTb)))); % minbtf
    table1(k,4) = mean(tempSb(find(tempTb == max(tempTb)))); % maxbtf
end
table1(find(table1 == 0)) = NaN;
figure(80)
set(gcf,'Position',[10 70 900 900])
print('TN_scatter_slope_baseline_individual','-dsvg','-r300','-painters')
savefig('TN_scatter_slope_baseline_individual.fig')
% figure(85)
% set(gcf,'Position',[1320 250 320 330])
% print([char(savename),'_scatter_slope_baseline'],'-dsvg','-r300','-painters')
% savefig([char(savename),'_scatter_slope_baseline.fig'])

% für Steigung m ist 1 der höchstmögliche Wert (unter der Annahme, dass wir
% von 24°C auf 34°C eine maximale Steigung von 0 auf 1 haben). In dem Fall
% ist 0.1 der höchstmögliche Wert, um das zu vereinfachen rechne ich alle
% Werte mal 10, sodass der mögliche Wertebereich von 0 bis 1 ist
slope(find(slope == 0)) = NaN;
slopeftb(find(slopeftb == 0)) = NaN;
slopebtf(find(slopebtf == 0)) = NaN;
corres(find(corres == 0)) = NaN;
rsq(find(rsq == 0)) = NaN;
rsqftb(find(rsqftb == 0)) = NaN;
rsqbtf(find(rsqbtf == 0)) = NaN;

% combine TN types
slopem(1:9) = NaN;
slopem(find([AllAni.TNtype] == 1)) = slopebtf([AllAni.TNtype] == 1);
slopem(find([AllAni.TNtype] == 2)) = slopeftb([AllAni.TNtype] == 2);
rsqm(1:9) = NaN;
rsqm(find([AllAni.TNtype] == 1)) = rsqbtf([AllAni.TNtype] == 1);
rsqm(find([AllAni.TNtype] == 2)) = rsqftb([AllAni.TNtype] == 2);

Boxplot_B([slope' slopem'],2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'baseline','stimulus'},[1 2])
ppos = find([AllAni.TNtype] == 1);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25);
scatter(x,slope(ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+1;
scatter(x,slopebtf(ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
ppos = find([AllAni.TNtype] == 2);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25);
scatter(x,slope(ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+1;
scatter(x,slopeftb(ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
xlim([.5 3.5])
ylim([-2 12])
ylabel('slope')
% set(gca,'XTick',[1 2],'XTickLabel',{'baseline','stimulus'})
set(gcf,'Position',[1020 250 320 330])
print('TN_bp_slope_baseline','-dsvg','-r300','-painters')
savefig('TN_bp_slope_baseline.fig')

Boxplot_B([rsq' rsqm'],2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'baseline','stimlus'},[1 2])
ppos = find([AllAni.TNtype] == 1);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+0;
scatter(x,rsq(ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+1;
scatter(x,rsqbtf(ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
ppos = find([AllAni.TNtype] == 2);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+0;
scatter(x,rsq(ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/25)+1;
scatter(x,rsqftb(ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
xlim([.5 3.5])
ylim([0 1])
ylabel('R^2')
% set(gca,'XTick',[1 2],'XTickLabel',{'baseline','stimulus'})
set(gcf,'Position',[1020 250 320 330])
print('TN_bp_rsq_baseline','-dsvg','-r300','-painters')
savefig('TN_bp_rsq_baseline.fig')

% LP für slope und r^2
figure
hold on
ppos = find([AllAni.TNtype] == 1);
plot([1 2],[rsq(ppos); rsqbtf(ppos)],'LineWidth',1.2,'Color','c')
ppos = find([AllAni.TNtype] == 2);
plot([1 2],[rsq(ppos); rsqftb(ppos)],'LineWidth',1.2,'Color','m')
plot(1:2,nanmedian([rsq' rsqm']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('R^2')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([0 1])
box on
print('TN_lp_rsq_baseline','-dsvg','-r300','-painters')
savefig('TN_lp_rsq_baseline.fig')
figure
hold on
ppos = find([AllAni.TNtype] == 1);
plot([1 2],[slope(ppos); slopebtf(ppos)],'LineWidth',1.2,'Color','c')
ppos = find([AllAni.TNtype] == 2);
plot([1 2],[slope(ppos); slopeftb(ppos)],'LineWidth',1.2,'Color','m')
plot(1:2,nanmedian([slope' slopem']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('m')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([-2 12])
box on
print('TN_lp_slope_baseline','-dsvg','-r300','-painters')
savefig('TN_lp_slope_baseline.fig')


%% Q10 für baseline
Q10.B = (B34./B26).^(10/(34-26));
Q10.Bftb = (B34ftb./B26ftb).^(10/(34-26));
Q10.Bbtf = (B34btf./B26btf).^(10/(34-26));

Q10.m(1:9) = NaN;
Q10.m(find([AllAni.TNtype] == 1)) = Q10.Bbtf([AllAni.TNtype] == 1);
Q10.m(find([AllAni.TNtype] == 2)) = Q10.Bftb([AllAni.TNtype] == 2);

Boxplot_B([Q10.B' Q10.m'],2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'baseline','stimulus'},[1 2])
ylabel('Q10 repetitive stim')
ylim([0 3])
xlim([.5 3.5])
hold on
ppos = find([AllAni.TNtype] == 1);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/10)+1-1;
scatter(x,Q10.B(1,ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/10)+2-1;
scatter(x,Q10.Bbtf(1,ppos),10,'filled','MarkerFaceColor','c','MarkerFaceAlpha',1)
ppos = find([AllAni.TNtype] == 2);
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/10)+1-1;
scatter(x,Q10.B(1,ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
x = ones(size(ppos)).*(1+(rand(size(ppos))*2-1)/10)+2-1;
scatter(x,Q10.Bftb(1,ppos),10,'filled','MarkerFaceColor','m','MarkerFaceAlpha',1)
set(gcf,'Position',[820 450 320 330])
% text(1,3.5,'ttest(1)')
% [~,ab] = ttest(Q10.B',ones(length(Q10.B),1));
% [~,aftb] = ttest(Q10.Bftb',ones(length(Q10.Bftb),1));
% [~,abtf] = ttest(Q10.Bbtf',ones(length(Q10.Bbtf),1));
% text(.5,3,[num2str(ab),', ',num2str(aftb),', ',num2str(abtf)]);
print('TN_bp_Q10_repetitiveStim_','-dsvg','-r300','-painters')
savefig('TN_bp_Q10_repetitiveStim_.fig')

figure
hold on
ppos = find([AllAni.TNtype] == 1);
plot(1:2,[Q10.B(ppos); Q10.Bbtf(ppos)],'LineWidth',1.2,'Color','c')
ppos = find([AllAni.TNtype] == 2);
plot(1:2,[Q10.B(ppos);Q10.Bftb(ppos)],'LineWidth',1.2,'Color','m')
plot(1:2,nanmedian([Q10.B' Q10.m']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('Q10 repetitive stim')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([0 3])
box on
% ft = friedman([Q10.B(4:end)' Q10.Bftb(4:end)' Q10.Bbtf(4:end)']);
% close
% abftb = signrank(Q10.B(4:end), Q10.Bftb(4:end));
% abbtf = signrank(Q10.B(4:end), Q10.Bbtf(4:end));
% aftbbtf = signrank(Q10.Bftb(4:end), Q10.Bbtf(4:end));
% text(.5,4.5,['friedman(0.01 0.003 0.0003): ',num2str(ft)]);
% text(.5,3.5,'wtest(a-ftb;a-btf;ftb-btf');
% text(.5,3,[num2str(abftb),', ',num2str(abbtf),', ',num2str(aftbbtf)]);
print('TN_lp_Q10_repetitiveStim_','-dsvg','-r300','-painters')
savefig('TN_lp_Q10_repetitiveStim_.fig')

