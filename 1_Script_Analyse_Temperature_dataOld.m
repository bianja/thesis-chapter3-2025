% analyse Arena data temperature
% this Script puts spiking information

% clear workspace
clear all;  %#ok<CLALL>
close all; clc

% find data file
filename = uigetdir('\\132.187.28.171\home\Data\Ephys');
data = dir([filename,'\*mat*']);
Ani = data(1).folder(39:41);

countR = 0; countRF = 0; c = 1; d = 1;
for i = 1 : size(data,1)
    % if contains(data(i).name, num2str(Ani)) == 1 % yes
        if contains(data(i).name, 'R0') == 1
            countR = countR + 1;
            posR(c) = i;
            c = c + 1;
        elseif contains(data(i).name, 'T0') == 1
            countRF = countRF + 1;
            posRF(d) = i;
            d = d + 1;
        end
    % end
end

% find stim mat file
c = 1;
for i = 1 : size(data,1)
    if contains(data(i).name, 'Arena') && ~contains(data(i).name, 'test') == 1 % yes
        posStim(c) = i;
        c = c + 1;
    end
end

%% for loop over all recordings (different temperature conditions)
for i = 1 : length(posR)
    
    %% load data file and stim file
    load([char(data(posR(i)).folder),'\',char(data(posR(i)).name)])
    load([char(data(posStim(i)).folder),'\',char(data(posStim(i)).name)])
    
    
    %% specify recording channels and stimuli
    if i == 1
        whos
        vars = inputdlg({'Number of units','Channel number of first unit','Number of stimuli repetitions'},'Customer', [1 30; 1 30; 1 30]);
        unitnbr = str2double(vars{1});
        firstunit = str2double(vars{2})-1;
        rep = str2double(vars{3});
        len = 120;
    end
    
       
    %% start and end times of stimuli
    
    stim_times_all = Ch1.times(find(Ch1.values < -0.015));
    stim_times_temp = zeros(len*rep*2,1);
    j = 1;
    for k = 1 : length(stim_times_all) - 1
        if stim_times_all(k+1) > stim_times_all(k)+1
            stim_times_temp(j+1) = stim_times_all(k+1);
            j = j + 1;
        end
    end
    
    stim_times_start = zeros(1,len*rep);
    stim_times_end = zeros(1,len*rep);
    c = 1;
    for k = 0 : 2 : length(stim_times_temp)-2
        stim_times_start(c) = stim_times_temp(k+1);
        stim_times_end(c) = stim_times_temp(k+2);
        c = c + 1;
    end
    
    start_diff = diff(stim_times_start);
    for k = 1 : length(start_diff)
        if start_diff(k) > 20 || start_diff(k) < 10 == 1
            start_diff(k) = NaN;
        end
    end
    
    diffVal = nanmedian(diff(stim_times_start(1:len)));
    
    [stim_times_start, stim_times_end] = findTimes(Ch1.times, Ch1.values, -0.015, rep, diffVal);
    
    
    %% find spike counts for each stimulus - yfreq
    for k = 1 : unitnbr
        count.(['unit0',num2str(k)]) = zeros(1,40*rep);
        f = 1;
        for j = 1 : 120*rep
            temp = eval(['Ch',num2str(firstunit+k)]);
            count.(['unit0',num2str(k)])(j) = (sum((temp.times > stim_times_start(f)...
                & temp.times < stim_times_end(f)) == 1))/(stim_times_end(f) - stim_times_start(f));
            temp2 = find((temp.times > stim_times_start(f) & temp.times < stim_times_end(f)));
            if isempty(temp2) == 1
                temp2 = NaN;
                delay.(['unit0',num2str(k)])(j) = NaN;
            else
                delay.(['unit0',num2str(k)])(j) = (temp.times(temp2(1)) - stim_times_start(f)) * 1000; % values in ms
            end
            f = f + 1;
        end
    end
    
    if stim_times_end(1)-stim_times_start(1) < 6
    else
        error('problem with start times or end times')
    end
    
    
    %% find positions and extract counts
    pos.tp2 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'translational'));
    pos.tp4 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'translational'));
    pos.tp8 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'translational'));
    pos.tn2 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'translational'));
    pos.tn4 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'translational'));
    pos.tn8 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'translational'));
    
    pos.rp2 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'rotational'));
    pos.rp4 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'rotational'));
    pos.rp8 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'rotational'));
    pos.rn2 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'rotational'));
    pos.rn4 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'rotational'));
    pos.rn8 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'rotational'));
    
    % calculate background activity
    for k = 1 : unitnbr
        [background.sum.(['unit0',num2str(k)]), background.raw.(['unit0',num2str(k)]), spikes.(['unit0',num2str(k)]), exclude.(['unit0',num2str(k)])] = testData(stim_times_start, stim_times_end, eval(['Ch',num2str(k+firstunit)]), rep, len);
        [~, posexclude.tp2.(['unit0',num2str(k)])] = intersect(pos.tp2,exclude.(['unit0',num2str(k)]));
        [~, posexclude.tp4.(['unit0',num2str(k)])] = intersect(pos.tp4,exclude.(['unit0',num2str(k)]));
        [~, posexclude.tp8.(['unit0',num2str(k)])] = intersect(pos.tp8,exclude.(['unit0',num2str(k)]));
        %     pos.tp2(posexclude.tp2.(['unit0',num2str(i)])) = NaN;
    end
    
    for k = 1 : unitnbr % this for loop is for the units
        tp2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp2);
        tp4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp4);
        tp8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp8);
        delpos.tp2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp2);
        delpos.tp4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp4);
        delpos.tp8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp8);
        ttpb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp2);
        ttpb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp4);
        ttpb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp8);
        
        tn2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn2);
        tn4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn4);
        tn8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn8);
        delpos.tn2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn2);
        delpos.tn4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn4);
        delpos.tn8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn8);
        ttnb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn2);
        ttnb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn4);
        ttnb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn8);
        
        rp2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp2);
        rp4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp4);
        rp8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp8);
        delpos.rp2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp2);
        delpos.rp4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp4);
        delpos.rp8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp8);
        trpb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp2);
        trpb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp4);
        trpb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp8);
        
        rn2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn2);
        rn4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn4);
        rn8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn8);
        delpos.rn2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn2);
        delpos.rn4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn4);
        delpos.rn8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn8);
        trnb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn2);
        trnb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn4);
        trnb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn8);
    end
    
    xtp2 = [rec(pos.tp2).vel];
    xtp4 = [rec(pos.tp4).vel];
    xtp8 = [rec(pos.tp8).vel];
    xtn2 = [rec(pos.tn2).vel]*(-1);
    xtn4 = [rec(pos.tn4).vel]*(-1);
    xtn8 = [rec(pos.tn8).vel]*(-1);
    xrp2 = [rec(pos.rp2).vel];
    xrp4 = [rec(pos.rp4).vel];
    xrp8 = [rec(pos.rp8).vel];
    xrn2 = [rec(pos.rn2).vel]*(-1);
    xrn4 = [rec(pos.rn4).vel]*(-1);
    xrn8 = [rec(pos.rn8).vel]*(-1);
    
    %% sort stimuli and plot
    
    for k = 1 : unitnbr
        [x_t_p2,temp] = sort(xtp2);
        y_t_p2.(['unit0',num2str(k)]) = tp2.(['unit0',num2str(k)])(temp);
        y_t_p2_delay.(['unit0',num2str(k)]) = delpos.tp2.(['unit0',num2str(k)])(temp);
        tpb2.(['unit0',num2str(k)]) = ttpb2.(['unit0',num2str(k)])(temp);
        
        [x_t_n2,temp] = sort(xtn2);
        y_t_n2.(['unit0',num2str(k)]) = tn2.(['unit0',num2str(k)])(temp);
        y_t_n2_delay.(['unit0',num2str(k)]) = delpos.tn2.(['unit0',num2str(k)])(temp);
        tnb2.(['unit0',num2str(k)]) = ttnb2.(['unit0',num2str(k)])(temp);
        
        [x_r_p2,temp] = sort(xrp2);
        y_r_p2.(['unit0',num2str(k)]) = rp2.(['unit0',num2str(k)])(temp);
        y_r_p2_delay.(['unit0',num2str(k)]) = delpos.rp2.(['unit0',num2str(k)])(temp);
        rpb2.(['unit0',num2str(k)]) = trpb2.(['unit0',num2str(k)])(temp);
        
        [x_r_n2,temp] = sort(xrn2);
        y_r_n2.(['unit0',num2str(k)]) = rn2.(['unit0',num2str(k)])(temp);
        y_r_n2_delay.(['unit0',num2str(k)]) = delpos.rn2.(['unit0',num2str(k)])(temp);
        rnb2.(['unit0',num2str(k)]) = trnb2.(['unit0',num2str(k)])(temp);
    end
    
    num = (len*rep)/rep/12;
    j = 1;
    for m = 1 : rep : num*rep % size(rec,2)/(3*rep)
        x_t_p2_plot(j) = x_t_p2(m); %#ok<*SAGROW>
        x_t_n2_plot(j) = x_t_n2(m);
        x_r_p2_plot(j) = x_r_p2(m);
        x_r_n2_plot(j) = x_r_n2(m);
        for k = 1 : unitnbr
            for c = 1 : rep
                raw.y_t_p2.(['unit0',num2str(k)])(c,j) = y_t_p2.(['unit0',num2str(k)])(m+c-1);
                raw.y_t_n2.(['unit0',num2str(k)])(c,j) = y_t_n2.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_p2.(['unit0',num2str(k)])(c,j) = y_r_p2.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_n2.(['unit0',num2str(k)])(c,j) = y_r_n2.(['unit0',num2str(k)])(m+c-1);
            end
            y_t_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p2.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n2.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p2.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n2.(['unit0',num2str(k)])(m:m+rep-1));            
            
            y_t_p2_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p2.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n2_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n2.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p2_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p2.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n2_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n2.(['unit0',num2str(k)])(m:m+rep-1));
            
            delpos.y_t_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp2.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_t_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn2.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp2.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn2.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawmean.(['unit0',num2str(k)]).translation.bw.w2(j) = nanmedian(tpb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).translation.fw.w2(j) = nanmedian(tnb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.cw.w2(j) = nanmedian(rpb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w2(j) = nanmedian(rnb2.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawsd.(['unit0',num2str(k)]).translation.bw.w2(j) = nanstd(tpb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).translation.fw.w2(j) = nanstd(tnb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.cw.w2(j) = nanstd(rpb2.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w2(j) = nanstd(rnb2.(['unit0',num2str(k)])(m:m+rep-1));
            
        end
        
        j = j + 1;
    end
    
%     for k = 1 : unitnbr
%         unfiltered_t_p2.(['unit0',num2str(k)]) = y_t_p2_plot.(['unit0',num2str(k)]);
%         unfiltered_t_n2.(['unit0',num2str(k)]) = y_t_n2_plot.(['unit0',num2str(k)]);
%         unfiltered_r_p2.(['unit0',num2str(k)]) = y_r_p2_plot.(['unit0',num2str(k)]);
%         unfiltered_r_n2.(['unit0',num2str(k)]) = y_r_n2_plot.(['unit0',num2str(k)]);
%         
%         y_t_p2_plot.(['unit0',num2str(k)]) = smooth(y_t_p2_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_t_n2_plot.(['unit0',num2str(k)]) = smooth(y_t_n2_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_p2_plot.(['unit0',num2str(k)]) = smooth(y_r_p2_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_n2_plot.(['unit0',num2str(k)]) = smooth(y_r_n2_plot.(['unit0',num2str(k)]),0.7,'loess');
%     end
    
    
    
    % sort stimuli
    for k = 1 : unitnbr
        [x_t_p4,temp] = sort(xtp4);
        y_t_p4.(['unit0',num2str(k)]) = tp4.(['unit0',num2str(k)])(temp);
        y_t_p4_delay.(['unit0',num2str(k)]) = delpos.tp4.(['unit0',num2str(k)])(temp);
        tpb4.(['unit0',num2str(k)]) = ttpb4.(['unit0',num2str(k)])(temp);
        
        [x_t_n4,temp] = sort(xtn4);
        y_t_n4.(['unit0',num2str(k)]) = tn4.(['unit0',num2str(k)])(temp);
        y_t_n4_delay.(['unit0',num2str(k)]) = delpos.tn4.(['unit0',num2str(k)])(temp);
        tnb4.(['unit0',num2str(k)]) = ttnb4.(['unit0',num2str(k)])(temp);
        
        [x_r_p4,temp] = sort(xrp4);
        y_r_p4.(['unit0',num2str(k)]) = rp4.(['unit0',num2str(k)])(temp);
        y_r_p4_delay.(['unit0',num2str(k)]) = delpos.rp4.(['unit0',num2str(k)])(temp);
        rpb4.(['unit0',num2str(k)]) = trpb4.(['unit0',num2str(k)])(temp);
        
        [x_r_n4,temp] = sort(xrn4);
        y_r_n4.(['unit0',num2str(k)]) = rn4.(['unit0',num2str(k)])(temp);
        y_r_n4_delay.(['unit0',num2str(k)]) = delpos.rn4.(['unit0',num2str(k)])(temp);
        rnb4.(['unit0',num2str(k)]) = trnb4.(['unit0',num2str(k)])(temp);
    end
    
    j = 1;
    for m = 1 : rep : num*rep % size(rec,2)/(3*rep)
        x_t_p4_plot(j) = x_t_p4(m);
        x_t_n4_plot(j) = x_t_n4(m);
        x_r_p4_plot(j) = x_r_p4(m);
        x_r_n4_plot(j) = x_r_n4(m);
        for k = 1 : unitnbr
            for c = 1 : rep
                raw.y_t_p4.(['unit0',num2str(k)])(c,j) = y_t_p4.(['unit0',num2str(k)])(m+c-1);
                raw.y_t_n4.(['unit0',num2str(k)])(c,j) = y_t_n4.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_p4.(['unit0',num2str(k)])(c,j) = y_r_p4.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_n4.(['unit0',num2str(k)])(c,j) = y_r_n4.(['unit0',num2str(k)])(m+c-1);
            end
            y_t_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p4.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n4.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p4.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n4.(['unit0',num2str(k)])(m:m+rep-1));
            
            y_t_p4_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p4.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n4_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n4.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p4_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p4.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n4_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n4.(['unit0',num2str(k)])(m:m+rep-1));
            
            delpos.y_t_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp4.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_t_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn4.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp4.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn4.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawmean.(['unit0',num2str(k)]).translation.bw.w4(j) = nanmedian(tpb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).translation.fw.w4(j) = nanmedian(tnb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.cw.w4(j) = nanmedian(rpb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w4(j) = nanmedian(rnb4.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawsd.(['unit0',num2str(k)]).translation.bw.w4(j) = nanstd(tpb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).translation.fw.w4(j) = nanstd(tnb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.cw.w4(j) = nanstd(rpb4.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w4(j) = nanstd(rnb4.(['unit0',num2str(k)])(m:m+rep-1));
        end
        j = j + 1;
    end
    
%     for k = 1 : unitnbr
%         unfiltered_t_p4.(['unit0',num2str(k)]) = y_t_p4_plot.(['unit0',num2str(k)]);
%         unfiltered_t_n4.(['unit0',num2str(k)]) = y_t_n4_plot.(['unit0',num2str(k)]);
%         unfiltered_r_p4.(['unit0',num2str(k)]) = y_r_p4_plot.(['unit0',num2str(k)]);
%         unfiltered_r_n4.(['unit0',num2str(k)]) = y_r_n4_plot.(['unit0',num2str(k)]);
%         
%         y_t_p4_plot.(['unit0',num2str(k)]) = smooth(y_t_p4_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_t_n4_plot.(['unit0',num2str(k)]) = smooth(y_t_n4_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_p4_plot.(['unit0',num2str(k)]) = smooth(y_r_p4_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_n4_plot.(['unit0',num2str(k)]) = smooth(y_r_n4_plot.(['unit0',num2str(k)]),0.7,'loess');
%     end
    
    
    % sort stimuli for width 8
    for k = 1 : unitnbr
        [x_t_p8,temp] = sort(xtp8);
        y_t_p8.(['unit0',num2str(k)]) = tp8.(['unit0',num2str(k)])(temp);
        y_t_p8_delay.(['unit0',num2str(k)]) = delpos.tp8.(['unit0',num2str(k)])(temp);
        tpb8.(['unit0',num2str(k)]) = ttpb8.(['unit0',num2str(k)])(temp);
        
        [x_t_n8,temp] = sort(xtn8);
        y_t_n8.(['unit0',num2str(k)]) = tn8.(['unit0',num2str(k)])(temp);
        y_t_n8_delay.(['unit0',num2str(k)]) = delpos.tn8.(['unit0',num2str(k)])(temp);
        tnb8.(['unit0',num2str(k)]) = ttnb8.(['unit0',num2str(k)])(temp);
        
        [x_r_p8,temp] = sort(xrp8);
        y_r_p8.(['unit0',num2str(k)]) = rp8.(['unit0',num2str(k)])(temp);
        y_r_p8_delay.(['unit0',num2str(k)]) = delpos.rp8.(['unit0',num2str(k)])(temp);
        rpb8.(['unit0',num2str(k)]) = trpb8.(['unit0',num2str(k)])(temp);
        
        [x_r_n8,temp] = sort(xrn8);
        y_r_n8.(['unit0',num2str(k)]) = rn8.(['unit0',num2str(k)])(temp);
        y_r_n8_delay.(['unit0',num2str(k)]) = delpos.rn8.(['unit0',num2str(k)])(temp);
        rnb8.(['unit0',num2str(k)]) = trnb8.(['unit0',num2str(k)])(temp);
    end
    
    j = 1;
    for m = 1 : rep : num*rep % size(rec,2)/(3*rep)
        x_t_p8_plot(j) = x_t_p8(m);
        x_t_n8_plot(j) = x_t_n8(m);
        x_r_p8_plot(j) = x_r_p8(m);
        x_r_n8_plot(j) = x_r_n8(m);
        for k = 1 : unitnbr
            for c = 1 : rep
                raw.y_t_p8.(['unit0',num2str(k)])(c,j) = y_t_p8.(['unit0',num2str(k)])(m+c-1);
                raw.y_t_n8.(['unit0',num2str(k)])(c,j) = y_t_n8.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_p8.(['unit0',num2str(k)])(c,j) = y_r_p8.(['unit0',num2str(k)])(m+c-1);
                raw.y_r_n8.(['unit0',num2str(k)])(c,j) = y_r_n8.(['unit0',num2str(k)])(m+c-1);
            end
            y_t_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p8.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n8.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p8.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n8.(['unit0',num2str(k)])(m:m+rep-1));
            
            %         y_t_p8_sd.(['unit0',num2str(k)])(j) = prctile(y_t_p8.(['unit0',num2str(k)])(i:i+rep-1),0.05);
            y_t_p8_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p8.(['unit0',num2str(k)])(m:m+rep-1));
            y_t_n8_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n8.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_p8_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p8.(['unit0',num2str(k)])(m:m+rep-1));
            y_r_n8_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n8.(['unit0',num2str(k)])(m:m+rep-1));
            
            delpos.y_t_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp8.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_t_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn8.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp8.(['unit0',num2str(k)])(m:m+rep-1));
            delpos.y_r_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn8.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawmean.(['unit0',num2str(k)]).translation.bw.w8(j) = nanmedian(tpb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).translation.fw.w8(j) = nanmedian(tnb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.cw.w8(j) = nanmedian(rpb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w8(j) = nanmedian(rnb8.(['unit0',num2str(k)])(m:m+rep-1));
            
            background.rawsd.(['unit0',num2str(k)]).translation.bw.w8(j) = nanstd(tpb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).translation.fw.w8(j) = nanstd(tnb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.cw.w8(j) = nanstd(rpb8.(['unit0',num2str(k)])(m:m+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w8(j) = nanstd(rnb8.(['unit0',num2str(k)])(m:m+rep-1));
        end
        j = j + 1;
    end
    
%     for k = 1 : unitnbr
%         unfiltered_t_p8.(['unit0',num2str(k)]) = y_t_p8_plot.(['unit0',num2str(k)]);
%         unfiltered_t_n8.(['unit0',num2str(k)]) = y_t_n8_plot.(['unit0',num2str(k)]);
%         unfiltered_r_p8.(['unit0',num2str(k)]) = y_r_p8_plot.(['unit0',num2str(k)]);
%         unfiltered_r_n8.(['unit0',num2str(k)]) = y_r_n8_plot.(['unit0',num2str(k)]);
%         
%         y_t_p8_plot.(['unit0',num2str(k)]) = smooth(y_t_p8_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_t_n8_plot.(['unit0',num2str(k)]) = smooth(y_t_n8_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_p8_plot.(['unit0',num2str(k)]) = smooth(y_r_p8_plot.(['unit0',num2str(k)]),0.7,'loess');
%         y_r_n8_plot.(['unit0',num2str(k)]) = smooth(y_r_n8_plot.(['unit0',num2str(k)]),0.7,'loess');
%     end
    
    
    %% analyse freq
    
    figure % tr, bw, freq
    [peak.tp(1,:), width.tp(1,:)] = plotFreq(x_t_p2_plot, y_t_p2_plot,[1/255 102/255 94/255],'x:',y_t_p2_sd,background.sum);
    [peak.tp(2,:), width.tp(2,:)] = plotFreq(x_t_p4_plot, y_t_p4_plot,[1/255 102/255 94/255],'x--',y_t_p4_sd,background.sum);
    [peak.tp(3,:), width.tp(3,:)] = plotFreq(x_t_p8_plot, y_t_p8_plot,[1/255 102/255 94/255],'x-',y_t_p8_sd,background.sum);
    createLegend({'2 bw','4 bw','8 bw'},[1/255 102/255 94/255; 1/255 102/255 94/255; 1/255 102/255 94/255],{':','--','-'},'best',1.5)    
    
    figure % tr, fw, freq
    [peak.tn(1,:), width.tn(1,:)] = plotFreq(x_t_n2_plot, y_t_n2_plot,[140/255 81/255 10/255],'x:',y_t_n2_sd,background.sum);
    [peak.tn(2,:), width.tn(2,:)] = plotFreq(x_t_n4_plot, y_t_n4_plot,[140/255 81/255 10/255],'x--',y_t_n4_sd,background.sum);
    [peak.tn(3,:), width.tn(3,:)] = plotFreq(x_t_n8_plot, y_t_n8_plot,[140/255 81/255 10/255],'x-',y_t_n8_sd,background.sum);
    createLegend({'2 fw','4 fw','8 fw'},[140/255 81/255 10/255; 140/255 81/255 10/255; 140/255 81/255 10/255],{':','--','-'},'best',1.5)
   
    figure % rot, cw, freq
    [peak.rp(1,:), width.rp(1,:)] = plotFreq(x_r_p2_plot, y_r_p2_plot,[1/255 102/255 94/255],'x:',y_r_p2_sd,background.sum);
    [peak.rp(2,:), width.rp(2,:)] = plotFreq(x_r_p4_plot, y_r_p4_plot,[1/255 102/255 94/255],'x--',y_r_p4_sd,background.sum);
    [peak.rp(3,:), width.rp(3,:)] = plotFreq(x_r_p8_plot, y_r_p8_plot,[1/255 102/255 94/255],'x-',y_r_p8_sd,background.sum);
    createLegend({'2 cw','4 cw','8 cw'},[1/255 102/255 94/255; 1/255 102/255 94/255; 1/255 102/255 94/255],{':','--','-'},'best',1.5)
    
    figure % rot, ccw, freq
    [peak.rn(1,:), width.rn(1,:)] = plotFreq(x_r_n2_plot, y_r_n2_plot,[140/255 81/255 10/255],'x:',y_r_n2_sd,background.sum);
    [peak.rn(2,:), width.rn(2,:)] = plotFreq(x_r_n4_plot, y_r_n4_plot,[140/255 81/255 10/255],'x--',y_r_n4_sd,background.sum);
    [peak.rn(3,:), width.rn(3,:)] = plotFreq(x_r_n8_plot, y_r_n8_plot,[140/255 81/255 10/255],'x-',y_r_n8_sd,background.sum);
    createLegend({'2 ccw','4 ccw','8 ccw'},[140/255 81/255 10/255; 140/255 81/255 10/255; 140/255 81/255 10/255],{':','--','-'},'best',1.5)
    
    % adaptAxes(4)
    
    
    %% analyse velocity
    
    % calculate values for x axis
    x = [0.5, 1, 5, 10, 20, 40, 60, 80, 100, 120];
    x2 = calcVelocity(x,2);
    x4 = calcVelocity(x,4);
    x8 = calcVelocity(x,8);
    %
    % figure % tr, bw
    % plotVelocity(x2, y_t_p2_plot,[1/255 102/255 94/255],'x:')
    % plotVelocity(x4, y_t_p4_plot,[1/255 102/255 94/255],'x--')
    % plotVelocity(x8, y_t_p8_plot,[1/255 102/255 94/255],'x-')
    % legend('2 bw','4 bw', '8 bw')
    %
    % figure % tr, fw
    % plotVelocity(x2, y_t_n2_plot,[140/255 81/255 10/255],'x:')
    % plotVelocity(x4, y_t_n4_plot,[140/255 81/255 10/255],'x--')
    % plotVelocity(x8, y_t_n8_plot,[140/255 81/255 10/255],'x-')
    % legend('2 fw', '4 fw', '8 fw')
    %
    % figure % rot, cw
    % plotVelocity(x2, y_r_p2_plot,[1/255 102/255 94/255],'x:')
    % plotVelocity(x4, y_r_p4_plot,[1/255 102/255 94/255],'x--')
    % plotVelocity(x8, y_r_p8_plot,[1/255 102/255 94/255],'x-')
    % legend('2 cw', '4 cw', '8 cw')
    %
    % figure % rot, ccw
    % plotVelocity(x2, y_r_n2_plot,[140/255 81/255 10/255],'x:')
    % plotVelocity(x4, y_r_n4_plot,[140/255 81/255 10/255],'x--')
    % plotVelocity(x8, y_r_n8_plot,[140/255 81/255 10/255],'x-')
    % legend('2 ccw', '4 ccw', '8 ccw')
    %
    % adaptAxes(4)
    
    
    %% put data in table
    % some information are the same for all trials (e.g. xfreq, xvelo, ..),
    % hence extracting those information are sufficient to do once
    % once: folder, animal, unitNbr, rep, xfreq, xvelo
    % for each: temp, yfreq, bg
    
    
    % clear AllAni
    clear AllAni
%     load Results_Temp.mat
%     save Results_Temp_backup.mat AllAni
    load Temperature_TN.mat 
    save Temperature_TN_backup.mat AllAni
    
    if i == 1
        % ask for units which should get saved
        vars = inputdlg({'Which units to save?'},'Customer', [1 30]);
        unitsave = str2num(vars{1});
    end
    
    for j = 1 : length(unitsave)
        
        if i == 1
            k = size(AllAni,2)+1;
            unit = unitsave(j);
            AllAni(k).File = file.name;
            AllAni(k).Animal = str2num(Ani); %str2double(file.name(65:66));
            AllAni(k).UnitNbr = unit;
            AllAni(k).rep = rep;
            AllAni(k).xfreq.w2 = x_t_p2_plot;
            AllAni(k).xfreq.w4 = x_t_p4_plot; %/2;
            AllAni(k).xfreq.w8 = x_t_p8_plot; %/4;
            
            AllAni(k).xvelo.w2 = x2;
            AllAni(k).xvelo.w4 = x4;
            AllAni(k).xvelo.w8 = x8;
        else
            k = size(AllAni,2)-length(unitsave)+j;
            unit = unitsave(j);
        end
        
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.bw.w2 = y_t_p2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.bw.w4 = y_t_p4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.bw.w8 = y_t_p8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.fw.w2 = y_t_n2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.fw.w4 = y_t_n4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.translation.fw.w8 = y_t_n8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.cw.w2 = y_r_p2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.cw.w4 = y_r_p4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.cw.w8 = y_r_p8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.ccw.w2 = y_r_n2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.ccw.w4 = y_r_n4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.mean.rotation.ccw.w8 = y_r_n8_plot.(['unit0',num2str(unit)]);
        
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.bw.w2 = y_t_p2_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.bw.w4 = y_t_p4_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.bw.w8 = y_t_p8_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.fw.w2 = y_t_n2_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.fw.w4 = y_t_n4_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.translation.fw.w8 = y_t_n8_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.cw.w2 = y_r_p2_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.cw.w4 = y_r_p4_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.cw.w8 = y_r_p8_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.ccw.w2 = y_r_n2_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.ccw.w4 = y_r_n4_sd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).yfreq.sd.rotation.ccw.w8 = y_r_n8_sd.(['unit0',num2str(unit)]);
        
        AllAni(k).(['R0',num2str(i)]).delay.translation.bw.w2 = delpos.y_t_p2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.translation.bw.w4 = delpos.y_t_p4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.translation.bw.w8 = delpos.y_t_p8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.translation.fw.w2 = delpos.y_t_n2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.translation.fw.w4 = delpos.y_t_n4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.translation.fw.w8 = delpos.y_t_n8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.cw.w2 = delpos.y_r_p2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.cw.w4 = delpos.y_r_p4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.cw.w8 = delpos.y_r_p8_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.ccw.w2 = delpos.y_r_n2_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.ccw.w4 = delpos.y_r_n4_plot.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).delay.rotation.ccw.w8 = delpos.y_r_n8_plot.(['unit0',num2str(unit)]);
        
        AllAni(k).(['R0',num2str(i)]).background.sum = background.sum.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).background.rawmean = background.rawmean.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).background.rawsd = background.rawsd.(['unit0',num2str(unit)]);
        AllAni(k).(['R0',num2str(i)]).peak.translation.bw = peak.tp(:,unit);
        AllAni(k).(['R0',num2str(i)]).peak.translation.fw = peak.tn(:,unit);
        AllAni(k).(['R0',num2str(i)]).peak.rotation.cw = peak.rp(:,unit);
        AllAni(k).(['R0',num2str(i)]).peak.rotation.ccw = peak.rn(:,unit);
        AllAni(k).(['R0',num2str(i)]).width.translation.bw = width.tp(:,unit);
        AllAni(k).(['R0',num2str(i)]).width.translation.fw = width.tn(:,unit);
        AllAni(k).(['R0',num2str(i)]).width.rotation.cw = width.rp(:,unit);
        AllAni(k).(['R0',num2str(i)]).width.rotation.ccw = width.rn(:,unit);
        
        if exist('Ch17') == 1 %yes
            if contains(Ch17.title, 'Temp') == 1 %yes
                AllAni(k).(['R0',num2str(i)]).Temp = median(Ch17.values(find(Ch17.times > stim_times_start(1) & Ch17.times < stim_times_end(end))));
            elseif exist('Ch19') == 1
                if contains(Ch19.title, 'Temp') == 1 %yes
                    AllAni(k).(['R0',num2str(i)]).Temp = median(Ch19.values(find(Ch19.times > stim_times_start(1) & Ch19.times < stim_times_end(end))));
                else
                    AllAni(k).(['R0',num2str(i)]).Temp = NaN;
                end
            else
                AllAni(k).(['R0',num2str(i)]).Temp = NaN;
            end
        elseif exist('Ch19') == 1
            if contains(Ch19.title, 'Temp') == 1 %yes
                AllAni(k).(['R0',num2str(i)]).Temp = median(Ch19.values(find(Ch19.times > stim_times_start(1) & Ch19.times < stim_times_end(end))));
            else
                AllAni(k).(['R0',num2str(i)]).Temp = NaN;
            end
        else
            AllAni(k).(['R0',num2str(i)]).Temp = NaN;
        end
                    
        
        save Temperature_TN.mat AllAni
%         save AllResultsRF.mat AllAni
%         save Results_Temp.mat AllAni
    end
    
    
    
end
