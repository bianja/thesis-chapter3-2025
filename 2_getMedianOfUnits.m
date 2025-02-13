% save median of repetitions in new data table

clearvars
close all
load('\\132.187.28.171\home\rest\Manuskript\III_Temperature_TN\data\1_Temperature_TN_correctX.mat')

% vectors to set first and last recording per animal and the position of
% the same temperature conditions
first = [2 4 1 3 1 1 1 1 1 1];
last = [3 9 6 6 2 2 6 6 6 9];

TC{1} = 3;
TW{1} = 2;
TC{2} = [7 8 9];
TW{2} = [4 5 6];
TC{3} = [1 2 3];
TW{3} = [4 5 6];
TC{4} = 6;
TW{4} = 4;
TC{5} = 1;
TW{5} = 2;
TC{6} = 1;
TW{6} = 2;
TC{7} = [2 6];
TW{7} = [1 3];
TC{8} = [4 5 6];
TW{8} = [1 2 3];
TC{9} = [4 5 6];
TW{9} = [1 2 3];
TC{10} = [4 5 6];
TW{10} = [1 2 3 7 8 9];

TemperatureC = [25 26 25.5 26 21 21 24 23 23 23];
TemperatureW = [32 30 30 32 32 32 33 32 32 32];

width = [2 4 8];

for i = 1 : 10 % for all bumblebees
    Data(i).File = AllAni(i).File;
    Data(i).Animal = AllAni(i).Animal;
    Data(i).UnitNbr = AllAni(i).UnitNbr;
    Data(i).rep = AllAni(i).rep;
    Data(i).xvelo = AllAni(i).xvelo;
    for d = 1 : 2 % for ftb and btf direction
        for w = 1 : 3 % for all three spatial frequencies
            if d == 1
                dir = 'fw';
            else
                dir = 'bw';
            end
            % median for room temperature condition
            clear tempC tempCB
            for j = 1 : size(TC{i},2)
                tempC(j,:) = AllAni(i).(['R0',num2str(TC{i}(j))]).yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]);
                tempCB(j) = AllAni(i).(['R0',num2str(TC{i}(j))]).background.sum(1);
                tempCSD(j) = AllAni(i).(['R0',num2str(TC{i}(j))]).background.sum(2);
            end
            if size(tempC,1) > 1 % only median if there are at least two curves
                Data(i).R01.yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]) = median(tempC);
            else
                Data(i).R01.yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]) = tempC;
            end
            Data(i).R01.Temp = TemperatureC(i);
            Data(i).R01.background.sum(1) = median(tempCB);
            Data(i).R01.background.sum(2) = median(tempCSD);
            % median for elevated temperature condition
            clear tempW
            for j = 1 : size(TW{i},2)
                tempW(j,:) = AllAni(i).(['R0',num2str(TW{i}(j))]).yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]);
                tempWB(j) = AllAni(i).(['R0',num2str(TW{i}(j))]).background.sum(1);
                tempWSD(j) = AllAni(i).(['R0',num2str(TW{i}(j))]).background.sum(2);
            end
            if size(tempW,1) > 1 % only median if there are at least two curves
                Data(i).R02.yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]) = median(tempW);
            else
                Data(i).R02.yfreq.mean.translation.(char(dir)).(['w',num2str(width(w))]) = tempW;
            end
            Data(i).R02.Temp = TemperatureW(i);
            Data(i).R02.background.sum(1) = median(tempWB);
            Data(i).R02.background.sum(2) = median(tempWSD);
        end
    end
    Data(i).TNtype = AllAni(i).TNtype;
end

clear AllAni
AllAni = Data;
save Temperature_TN_correctX_median.mat AllAni

