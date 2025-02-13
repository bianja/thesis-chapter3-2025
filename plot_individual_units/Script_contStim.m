% script to extract and visualize response rate during continouous
% stimulation with the same stiumulus

%% set unit information
firstunit = 27-1;                                                
unitnbr = 1; %adapt for each recording
m = 1; % -1 if error pops for amplitude as it can't be -ve
%% find start and end times of stimulation

% number of stimuli
stimNum = size(rec,2);

% find pins
stim_times_all = Ch1.times(find(Ch1.values < -0.02)); 
stim_times_temp = zeros(stimNum*2,1);
j = 1; 
stim_times_temp(1) = stim_times_all(1);
for i = 1 : length(stim_times_all) - 1
    if stim_times_all(i+1) > stim_times_all(i)+1
        stim_times_temp(j+1) = stim_times_all(i+1);
        j = j + 1;
    end
end



stim_times_start = zeros(1,stimNum);
stim_times_end = zeros(1,stimNum);
c = 1;
for i = 0 : 2 : length(stim_times_temp)-2 
    stim_times_start(c) = stim_times_temp(i+1);  
    stim_times_end(c) = stim_times_temp(i+2); 
    c = c + 1;
end

% start_diff = diff(stim_times_start);
% for i = 1 : length(start_diff)
%     if start_diff(i) > 20 || start_diff(i) < 10 %== 1
%         start_diff(i) = NaN;
%     end
% end
% diffVal = nanmedian(diff(stim_times_start(1:stimNum))); 
% 
% stim_times_start = stim_times_start(1) : diffVal : diffVal*stimNum+diffVal;
% stim_times_end = stim_times_end(1) : diffVal : diffVal*stimNum+diffVal;

%% find spikes during stimuli
for k = 1 : unitnbr
    count.(['unit0',num2str(k)]) = zeros(1,stimNum);
    f = 1;
        for i = 1 : stimNum
            temp = eval(['Ch',num2str(firstunit+k)]);
            count.(['unit0',num2str(k)])(i) = (sum((temp.times > stim_times_start(f)...
                & temp.times < stim_times_end(f)) == 1))/(stim_times_end(f) - stim_times_start(f));
            temp2 = find((temp.times > stim_times_start(f) & temp.times < stim_times_end(f)));
            if isempty(temp2) == 1
                temp2 = NaN;            
                delay.(['unit0',num2str(k)])(i) = NaN;
            else
                delay.(['unit0',num2str(k)])(i) = (temp.times(temp2(1)) - stim_times_start(f)) * 1000; % values in ms
            end
            f = f + 1;
        end
        % separate bw and fw response
        bw.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(1:2:end);
        fw.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(2:2:end);
end

%% extract spike shape
times.ex = Ch1.times(1):60:Ch1.times(end)-5; % time points after every 60 secs
for k = 1 : unitnbr
    c = 1;
    temp = eval(['Ch',num2str(firstunit+k)]);
    for j = times.ex 
        d = find(temp.times>j);
        d = d(1); % total time
        for i = 1 : 4 
            shape.(['unit0',num2str(k)]){i,c} = mean(temp.values(d:d+100,:,i));
            shapesd.(['unit0',num2str(k)]){i,c} = std(temp.values(d:d+100,:,i));
        end
        c = c + 1;
    end
end
% mean spike shape of 100 spikes extracted after intervals of 60 secs (look at j & d) 
%% plot

% figure- spike rate v/s time
x = stim_times_start(2:2:end); % bw
y = stim_times_start(1:2:end); % fw
unit = unitnbr; % unit of interest

figure
plot(x-x(1),bw.(['unit0',num2str(unit)]),'-','Color',[.7 .7 .7],'LineWidth',1.5)
hold on
plot(y-y(1),fw.(['unit0',num2str(unit)]),'-','Color',[.3 .3 .3],'LineWidth',1.5)
set(gca,'Box','on')
xlabel('time (s)')
ylabel('spike rate s^-^1')
legend('bw','fw')
set(gcf,'Position',[850 460 390 343])
axis square
xlim([x(1) x(end)-x(1)])
ylim([0 40])

print(['RS_A163_U5_Ch27_T12'],'-dsvg','-r300','-painters')
savefig(['RS_A163_U5_Ch27_T12.fig'])

%% spike shape across all stim for all 4 channels
limity = [-3*10e-5 2.5*10e-5];
x = (0:temp.resolution/2:(temp.items*temp.resolution/2)-temp.resolution/2)*1000;
figure
%k = (['unit0',num2str(unit)])
subplot(2,2,1)
plot(x, shape.(['unit0',num2str(unit)]){1,1},'Color',[.7 .7 .7]) % grey
xlabel('time (ms)')
ylabel('mV')
hold on
plot(x, shape.(['unit0',num2str(unit)]){1,end},'Color',[.3 .3 .3]) % black
set(gca,'YLim',limity)

subplot(2,2,2)
plot(x,shape.(['unit0',num2str(unit)]){2,1},'Color',[.6 .3 .3]) % red
xlabel('time (ms)')
ylabel('mV')
hold on
plot(x,shape.(['unit0',num2str(unit)]){2,end},'Color',[.3 .3 .3])
set(gca,'YLim',limity)

subplot(2,2,3)
plot(x,shape.(['unit0',num2str(unit)]){3,1},'Color',[.3 .3 .8]) % blue
xlabel('time (ms)')
ylabel('mV')
hold on
plot(x,shape.(['unit0',num2str(unit)]){3,end},'Color',[.3 .3 .3])
set(gca,'YLim',limity)

subplot(2,2,4)
plot(x,shape.(['unit0',num2str(unit)]){4,1},'Color',[.3 .5 .3]) % green
xlabel('time (ms)')
ylabel('mV')
hold on
plot(x,shape.(['unit0',num2str(unit)]){4,end},'Color',[.3 .3 .3])
set(gca,'YLim',limity)

% pks= peaks, proms= prominence i.e height, 
% widths= with of spike at half-prominence, locs= location of spike
for i= 1: size(shape.(['unit0',num2str(unit)]), 2)
    for j = 1:4 % channels
[pks,locs,widths,proms]=findpeaks(shape.(['unit0',num2str(unit)]){j,i}*m, 'Annotate','extents');
peak=find(max(shape.(['unit0',num2str(unit)]){j,i}(locs)*m)==pks); %finds position of highest peak
amplitude(j,i)= proms(peak); % finds amplitude of highest peak
    end
end

figure % plot amplitudes for all 4 channels
for i= 1:4 % channels
    subplot(2,2,i)
plot(amplitude(i,:),'-','Color',[.3 .5 .4],'LineWidth',1.5)
set(gca,'Box','on')
xlabel('time(s)') % use time instead
ylabel('amplitude')
end

figure % plot all spike shapes and the mean 
subplot(2,2,1)
hold on
for i = 1:size(shape.(['unit0',num2str(unit)]), 2) % total no.of extracted spikes 
plot(x, shape.(['unit0',num2str(unit)]){1,i},'Color',[.7 .7 .7 .5]) % grey
temp3(i, :) = shape.(['unit0',num2str(unit)]){1,i};
end
plot(x, mean(temp3),'Color',[.3 .3 .3], 'LineWidth',2)
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,2)
hold on 
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)
plot(x,shape.(['unit0',num2str(unit)]){2,i},'Color',[.6 .3 .3 .5]) % red
temp3(i, :) = shape.(['unit0',num2str(unit)]){2,i};
end
plot(x, mean(temp3),'Color',[.3 .3 .3], 'LineWidth',2)
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,3)
hold on 
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)
plot(x,shape.(['unit0',num2str(unit)]){3,i},'Color',[.3 .3 .8 .5]) % blue
temp3(i, :) = shape.(['unit0',num2str(unit)]){3,i};
end
plot(x, mean(temp3),'Color',[.3 .3 .3],'LineWidth',2)
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,4)
hold on 
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)
plot(x,shape.(['unit0',num2str(unit)]){4,i},'Color',[.3 .5 .3 .5]) % green
temp3(i, :) = shape.(['unit0',num2str(unit)]){4,i};
end
plot(x, mean(temp3),'Color',[.3 .3 .3], 'LineWidth',2)
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

%% 3D plot of spike shapes (plot3 command)
figure
subplot(2,2,1)
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)    
plot3(x, shape.(['unit0',num2str(unit)]){1,i},ones(length(x),1)+ i,'Color',[.4+i/500 .3 .3]) % Ch1
hold on
end
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,2)
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)    
plot3(x, shape.(['unit0',num2str(unit)]){2,i},ones(length(x),1)+ i,'Color',[.5+i/500 .4 .8]) % Ch2
hold on
end
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,3)
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)    
plot3(x, shape.(['unit0',num2str(unit)]){3,i},ones(length(x),1)+ i,'Color',[.5+i/500 .7 .7]) % Ch3
hold on
end
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)

subplot(2,2,4)
for i = 1:size(shape.(['unit0',num2str(unit)]), 2)    
plot3(x, shape.(['unit0',num2str(unit)]){4,i},ones(length(x),1)+ i,'Color',[.4+i/500 .5 .3]) % Ch4
hold on
end
xlabel('time (ms)')
ylabel('mV')
set(gca,'YLim',limity)
