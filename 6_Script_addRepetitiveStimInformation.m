clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\3_Temperature_TN_correctX_median_BaselineCorrectionFit.mat')
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\AniInfoTN.mat');
load('PhysitempData.mat')


%% go into each folder and load continuous stimulus information and temperature information at start and end of time period
AllAni(1).C01 = []; saveData = 0;
for a = 2 : size(AllAni,2) %A86 does not hold repetitive stimulus information
    if AllAni(a).Animal == AllAni(a-1).Animal
        if ~isempty(AllAni(a-1).C01)
            saveData = 1;
        end
    else
        clearvars -except AllAni a AniInfoTN fitT fitmV saveData savepos
        cd \\132.187.28.171\home\rest\data\Ephys\
        folders = dir; % list all folders and change path to folder of current animal
        for k = 1 : size(folders,1)
            if contains(folders(k).name,['A',num2str(AllAni(a).Animal)])
                cd(folders(k).name)
                cd(AllAni(a).File(strfind(AllAni(a).File,'Rec'):strfind(AllAni(a).File,'Rec')+11))
            end
        end
        files = dir('*mat');
        currInfo = find(AniInfoTN.Ani==AllAni(a).Animal);
        c = 1; b = 1;
        for g = 1 : size(files,1)
            if contains(files(g).name,'C0')
                filesSpike(b) = g;
                b = b + 1;
            elseif contains(files(g).name,'Continuous')
                filesRec(c) = g;
                c = c + 1;
            end
        end
        if exist('filesSpike','var')
            for g = 1 : length(filesSpike)
                load(files(filesSpike(g)).name)
                load(files(filesRec(g)).name)

                %% set unit information
                firstunit = AniInfoTN.first(currInfo)-1; % real first channel number (-1 is because of the for loops)
                unitnbr = AniInfoTN.sum(currInfo);
                Tval = 32; % temperature at the beginning of the stimulus period
                if rec(1).vel > 0 % btf
                    cbtf = 1; cftb = 2;
                else
                    cbtf = 2; cftb = 1;
                end

                %% find start and end times of stimuli
                clear peakTimes stim_times_start stim_times_end
                [~,peakVal] = findpeaks(Ch1.values*-1,'MinPeakHeight',0.02,'MinPeakDistance',500);
                peakTimes = Ch1.times(peakVal);
                stim_times_start = peakTimes(1:2:end);
                stim_times_end = peakTimes(2:2:end);

                %% get temperature information
                clear Tstim bgTstim
                if ~isempty(dir('*xls*')) % temperature tracked by multimeter, excel file exits
                    % load excel file
                    temp = dir('*xls*');
                    TData = readtable(temp.name);
                    startTime = files(filesRec(g)-1).date;
                    if contains(startTime,'Mrz')
                        startTime(strfind(startTime,'Mrz'):strfind(startTime,'Mrz')+2) = 'Mar';
                    elseif contains(startTime,'Mai')
                        startTime(strfind(startTime,'Mai'):strfind(startTime,'Mai')+2) = 'May';
                    elseif contains(startTime,'Okt')
                        startTime(strfind(startTime,'Okt'):strfind(startTime,'Okt')+2) = 'Oct';
                    elseif contains(startTime,'Dez')
                        startTime(strfind(startTime,'Dez'):strfind(startTime,'Dez')+2) = 'Dec';
                    end
                    % find T value for stimuli and save
                    for k = 1 : length(stim_times_start)
                        index = [0 0];
                        [~,index(1)] = ismember(datetime(startTime)+seconds(round(stim_times_start(k))),TData.Date_Time);
                        [~,index(2)] = ismember(datetime(startTime)+seconds(round(stim_times_end(k))),TData.Date_Time);
                        if index(1) == 0
                            in = 1;
                            while index(1) == 0
                                [~,index(1)] = ismember(datetime(startTime)+seconds(round(stim_times_start(k))+in),TData.Date_Time);
                                in = in + 1;
                            end
                        end
                        if index(2) == 0
                            in = 1;
                            while index(2) == 0
                                [~,index(2)] = ismember(datetime(startTime)+seconds(round(stim_times_end(k))+in),TData.Date_Time);
                                in = in + 1;
                            end
                        end
                        Tstim(k) = median(str2num(cell2mat(TData.Value(index(1):index(2)))))+1.0222;
                        % if k < length(stim_times_start)%-1
                            bgTstim(k) = median(str2num(cell2mat(TData.Value(index(2)+5:index(2)+9))))+1.0222;
                        % end
                    end

                elseif exist('Ch17','var') % temperature tracked by Physitemp
                    if AllAni(a).Animal ~= 131 && AllAni(a).Animal ~= 133
                        % find start and end times of stimulus sets
                        for k = 1 : length(stim_times_start)
                            Tstim(k) = round(fitT(find(round(fitmV) == round(median(Ch17.values(find(Ch17.times==stim_times_start(k)):find(Ch17.times==stim_times_end(k))))*10000),1,'first')),1);
                            if k < length(stim_times_start)%-1
                                bgTstim(k) = round(fitT(find(round(fitmV) == round(median(Ch17.values(find(Ch17.times==stim_times_start(k)):find(Ch17.times==stim_times_end(k))))*10000),1,'first')),1);
                            end
                        end
                    else
                        Tstim(1:2) = NaN;
                        if k < length(stim_times_start)%-1
                            bgTstim(k) = NaN;
                        end
                    end
                else % no temperature was tracked
                    Tstim(1:2) = NaN;
                    if k < length(stim_times_start)%-1
                        bgTstim(k) = NaN;
                    end
                end

                %% find spikes during stimuli
                clear shape bw fw count countbg
                for k = 1 : unitnbr
                    count.(['unit0',num2str(k)]) = zeros(1,length(stim_times_start));
                    countbg.(['unit0',num2str(k)]) = zeros(1,length(stim_times_start));
                    f = 1;
                    for i = 1 : length(stim_times_start)
                        temp = eval(['Ch',num2str(firstunit+k)]);
                        count.(['unit0',num2str(k)])(i) = (sum((temp.times > stim_times_start(f)...
                            & temp.times < stim_times_end(f)) == 1))/(stim_times_end(f) - stim_times_start(f));
                        if i <= length(stim_times_start) - 1
                            countbg.(['unit0',num2str(k)])(i) = (sum((temp.times > stim_times_end(f)...
                                & temp.times < stim_times_start(f+1)) == 1))/(stim_times_start(f+1) - stim_times_end(f));
                        end
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
                    bw.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(cbtf:2:end);
                    fw.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(cftb:2:end);
                    bg.(['unit0',num2str(k)]) = countbg.(['unit0',num2str(k)]);
                end

                %% save
                AllAni(a).(['C0',num2str(g)]).bg.T = bgTstim;
                AllAni(a).(['C0',num2str(g)]).bg.spikes = bg.(['unit0',num2str(AllAni(a).UnitNbr)]);
                AllAni(a).(['C0',num2str(g)]).ftb.T = Tstim(cftb:2:end);
                AllAni(a).(['C0',num2str(g)]).btf.T = Tstim(cbtf:2:end);
                AllAni(a).(['C0',num2str(g)]).btf.spikes = bw.(['unit0',num2str(AllAni(a).UnitNbr)]);
                AllAni(a).(['C0',num2str(g)]).ftb.spikes = fw.(['unit0',num2str(AllAni(a).UnitNbr)]);
            end
        end
    end
    if saveData == 1
        for k = 1 : 6
            if ~isempty(AllAni(a-1).(['C0',num2str(k)]))
                AllAni(a).(['C0',num2str(k)]).bg.T = bgTstim;   
                AllAni(a).(['C0',num2str(k)]).bg.spikes = bg.(['unit0',num2str(AllAni(a).UnitNbr)]);
                AllAni(a).(['C0',num2str(k)]).ftb.T = Tstim(cftb:2:end);
                AllAni(a).(['C0',num2str(k)]).btf.T = Tstim(cbtf:2:end);
                AllAni(a).(['C0',num2str(k)]).btf.spikes = bw.(['unit0',num2str(AllAni(a).UnitNbr)]);
                AllAni(a).(['C0',num2str(k)]).ftb.spikes = fw.(['unit0',num2str(AllAni(a).UnitNbr)]);            
            end
        end
    end
    saveData = 0;
end