%% set variables providing pin information and spike information
% usually, Ch14 is providing pin informataion (and Ch13, but we are using
% Ch14 here). The first channel higher than 14 should contain spike
% information from all units. All the following channels provide spike
% information for individual units

% ask for name length of variables -> if length == 4 and position 3 and 4
% are int and int at position 3&4 are higher than 14 -> channel contains
% spike information
% close all

% find stimuli to plot
close all
unit = 5;
T = 1;
ft = 20;
fs = 4;
stimshow = find([rec.width] == fs & [rec.vel] == ft | [rec.width] == fs & [rec.vel] == ft*-1);
startunit = 5;
endunit = startunit;
% startrec = 2;
% endrec = 2;
limup = 40;
% obj.setRLim([0 10])

c = 1; i = []; spikeChannels = []; vars = [];
vars = who;
for i = 1 : length(who)
    if length(vars{i}) == 4 && isempty(str2num(vars{i}(3))) == 0 && isempty(str2num(vars{i}(4))) == 0 && str2num([(vars{i}(3)) (vars{i}(4))]) > 14 % having a number and channel higher 14
        spikeChannels(c) = str2num([(vars{i}(3)) (vars{i}(4))]);
        c = c + 1;
    else
    end
end

spikeChannels = [min(spikeChannels)-1 spikeChannels];

% to analyse information of single units -> remove first channel value as
% this is all units spiking information
spikeChannels = spikeChannels(spikeChannels~=spikeChannels(1));

%% find pins to find time points of stimuli
if exist('Ch1')
    stim = findTimes_ReceptiveField(Ch1.times, Ch1.values, -0.015, rec); % Ch14 or Ch1
else
    stim = findTimes_ReceptiveField(Ch14.times, Ch14.values, -0.015, rec); % Ch14 or Ch1
end

%% calculate background activity

for u = startunit:endunit %length(spikeChannels)
    [baseline, spikesBG, exclude] = testData_ReceptiveField(stim(:,1), stim(:,6), eval(['Ch',num2str(spikeChannels(u))]), rec);
    
    % calculate background acitivty from Hz in counts/bin -> but this only
    % works for 40 Hz so far
%     temp1 = baseline(1) / mean(diff(stim(stimshow(1),:))); % bg divided by duration of stimulus
%     temp2 = temp1 / 36;
%     bg = temp2 * 5;
%     baseline(1)
    bg = baseline(1)/36*5;
    bg0 = mean(spikesBG)/36*5;
    bg1 = (mean(spikesBG)+std(spikesBG))/36*5;
    bg2 = (mean(spikesBG)-std(spikesBG))/36*5;
    % bg = median(spikesBG(stimshow-1));
    duration = (360/calcVelocity(fs,ft))*5; % duration of stimulus for 5 repetitions

    
    %% find spike values
    clear spikes
    
    for j = 1 : size(stim,1)
        for k = 1 : size(stim,2)-1
            temp = find(eval(['Ch',num2str(spikeChannels(u))]).times > stim(j,k) & eval(['Ch',num2str(spikeChannels(u))]).times < stim(j,k+1));
            stim_time(k) = stim(j,k+1)-stim(j,k);
            spike_rate(j,k) = length(temp);%-baseline(1);
            spikes{j,k} = eval(['Ch',num2str(spikeChannels(u))]).times(temp);
            clear temp
        end
        spike_rate_temp(j,:) = spike_rate(j,:);
        spike_rate(j,:) = spike_rate(j,:)/mean(stim_time);
    end
    stimnum = size(spikes,1);
    
    %% calculate spike count per stim
    for s = 1 : stimnum
        clear tAll
        tAll = [];
        for r = 1 : 5
            emptcount = 0;
            clear t
            t = (spikes{s,r}-stim(s,r))*360/diff(stim(s,r:r+1));
            if rec(s).vel < 0
                t = 360-t;
            else
            end
            if isempty(t)
                round{s,r} = NaN;
            else
                round{s,r} = t;
            end
            tAll = [tAll; t];
        end
        if isempty(tAll) == 1 % true
            rec(s).spikes = NaN;
        else
            rec(s).spikes = tAll;
        end
    end
    
    %% plot spike timings (spikes)
    count = 1;
    for s = stimshow % startrec : endrec
        tAll = [];
        if rec(s).vel < 0
            figure('Name',['U',num2str(spikeChannels(u)),' CCW ',num2str(rec(s).vel*-1)],'NumberTitle','off')
            names = num2str(['U',num2str(spikeChannels(u)),'CCW',num2str(rec(s).vel*-1)]);
        else
            figure('Name',['U',num2str(spikeChannels(u)),' CW ',num2str(rec(s).vel)],'NumberTitle','off')
            names = num2str(['U',num2str(spikeChannels(u)),'CW',num2str(rec(s).vel)]);
        end
        for r = 1 : 5
            subp.(names) = subplot(2,3,r,polaraxes);
            if isnan(round{s,r})
            else
                CircHist(round{s,r},36,'areAxialData',false,'parent',subp.(names));
                set(gca,'ThetaZeroLocation', 'top','ThetaDir','clockwise')
            end
        end
        if isempty(rec(s).spikes)
        else
            if rec(s).vel < 0
                figure('Name',['Mean U',num2str(spikeChannels(u)),' CCW ',num2str(rec(s).vel*-1)],'NumberTitle','off')
                CCWData = rec(s).spikes;
                names = num2str(['Mean',names]);
                obj = CircHist(rec(s).spikes,36,'areAxialData',false,'parent',polaraxes,'histType','frequency','binSizeSec',duration/(36));
                set(gca,'ThetaZeroLocation', 'top','ThetaDir','counterclockwise')
            else
                figure('Name',['Mean U',num2str(spikeChannels(u)),' CW ',num2str(rec(s).vel)],'NumberTitle','off')
                CWData = rec(s).spikes;
                names = num2str(['Mean',names]);
                obj = CircHist(rec(s).spikes,36,'areAxialData',false,'parent',polaraxes,'histType','frequency','binSizeSec',duration/(36));
                set(gca,'ThetaZeroLocation', 'top','ThetaDir','clockwise')
            end
            %             figure('Name',['Mean U',num2str(spikeChannels(u)),' CCW ',num2str(rec(s).vel*-1)],'NumberTitle','off')
%             bg = spikesBG(s-1);
            h = obj.drawCirc(bg0, '-g', 'LineWidth', 2,'Color',[5/255 48/255 97/255]);
            h = obj.drawCirc(bg1, '--g', 'LineWidth', 2,'Color',[5/255 48/255 97/255]);
            h = obj.drawCirc(bg2, '--g', 'LineWidth', 2,'Color',[5/255 48/255 97/255]);
            obj.setRLim([-1 limup])
        end
        % save data to plot in new plot (not circular but in catesian coordinate system)
        spikeData(count,:) = obj.histData(:,1);
        count = count + 1;
    end
end

%% plot receptife fields
figure
hold on
plot(obj.edges(1:end-1)+mean(diff(obj.edges))/2,spikeData(1,:),'LineWidth',1.5,'Color',[.3 .5 .8])
plot(obj.edges(1:end-1)+mean(diff(obj.edges))/2,spikeData(2,:),'LineWidth',1.5,'Color',[.8 .5 .7])
xlim([0 360])
set(gca,'XTick',[0 45 90 135 180 225 270 315 360],'XTickLabels',{'-180','-135','-90','-45','0','45','90','135','180'})
upLim = get(gca,'YLim');
plot([180 180],[0 upLim(2)],'--','Color',[.7 .7 .7],'LineWidth',2)
if rec(s).vel < 0
    legend('cw','ccw')
else
    legend('ccw','cw')
end

%%
figure(2)
savefig([char(file.name(41:44)),'_U',num2str(unit),'_RF_ccw_ft',num2str(ft),'_fs',num2str(fs),'_T',num2str(T),'.fig'])
print([char(file.name(41:44)),'_U',num2str(unit),'_RF_ccw_ft',num2str(ft),'_fs',num2str(fs),'_T',num2str(T)],'-dsvg','-r300','-painters')

figure(4)
savefig([char(file.name(41:44)),'_U',num2str(unit),'_RF_cw_ft',num2str(ft),'_fs',num2str(fs),'_T',num2str(T),'.fig'])
print([char(file.name(41:44)),'_U',num2str(unit),'_RF_cw_ft',num2str(ft),'_fs',num2str(fs),'_T',num2str(T)],'-dsvg','-r300','-painters')
