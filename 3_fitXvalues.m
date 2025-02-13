% Problem: Units were recorded with two different recording sets (Velocity
% and temporal frequency)
% Apply PCHIP on data to get fit that keeps real recording values and fits
% xrange between these values. Afterwards set a x-vector that is the same
% for all units.

clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\2_Temperature_TN_correctX_median.mat')

first = [1 1 1 1 1 1 1 1 1 1];
last = [3 9 7 7 2 2 7 7 7 9];
width = [2 4 8];

% baseline correction
for i = 1 : size(AllAni,2)
    for j = 1 : 2 %first(i) : last(i)
        for w = 1 : 3
            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.(['w',num2str(width(w))]) = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.(['w',num2str(width(w))]) - AllAni(i).(['R0',num2str(j)]).background.sum(1);%rawmean.translation.fw.(['w',num2str(width(w))]);  
            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.(['w',num2str(width(w))]) = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.(['w',num2str(width(w))]) - AllAni(i).(['R0',num2str(j)]).background.sum(1);%rawmean.translation.bw.(['w',num2str(width(w))]);  
        end
    end
end

for i = 1 : size(AllAni,2)
    for w = 1 : 3
        for d = 1 : 2
            if d == 1
                dir = 'fw';
            else
                dir = 'bw';
            end
            for j = 1 : 2 %first(i) : last(i)
                if length(AllAni(i).xvelo.(char(['w',num2str(width(w))]))) < length(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir)).(char(['w',num2str(width(w))])))
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir)).(char(['w',num2str(width(w))]))(length(AllAni(i).xvelo.(char(['w',num2str(width(w))])))+1:end) = []; % correct for not presented stimuli
                end
                f = createFit(AllAni(i).xvelo.(char(['w',num2str(width(w))])),AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir)).(char(['w',num2str(width(w))]))); % Interpolant - shape preserving (PCHIP)
                X = linspace(22.5, 2700, 1000); % get x values
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir)).(char(['w',num2str(width(w))])) = f(X);
            end
        end
        AllAni(i).xvelo.(char(['w',num2str(width(w))])) = [];
        AllAni(i).xvelo.(char(['w',num2str(width(w))])) = X;
    end
end

save Temperature_TN_correctX_median_BaselineCorrectionFit.mat AllAni