%% Plot the max value of the stimulation artifact as a function of distance
clear all, clc, close all
cd('C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis');
directory = 'C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis';
load('5.1.14_depth.mat');
listing = dir(directory);
% distances = [30.213, 31.19, 31.723, 31.902, 32.268, 33.392, 33.054; ...
%  39, 39, 39, 39, 39, 39, 39];

counter = 0;
minArt = 10^9;
maxArt = 0;
cc=hsv(size(distances,2));

for j = 1:size(listing,1)
    s = listing(j).name;
    if strfind(s,'proc')
        counter = counter + 1;
        load(listing(j).name);
        if max(F.stim(1).meanArt)> maxArt
            maxArt = max(F.stim(1).meanArt);
        end
        if max(F.stim(1).meanArt)<minArt
            minArt = max(F.stim(1).meanArt);
        end
        % determine the electrode site that's stimulated through and thus
        % the color code:
        switch distances(2, counter)
            case 39
                color = 'r';
            case 31
                color = 'g';
            case 15
                color = 'b';
            case 7
                color = 'c';
        end
        
        % plot stim 1
        h(1) = subplot(2,3,1)
        plot([0:1/F.s_rate:(size(F.stim(1).meanArt,2)-1)/F.s_rate], F.stim(1).meanArt, 'color', cc(counter,:)); hold on;  title('350uA stim mean artifacts'); xlabel('Time(ms)'); ylabel('voltage(uV)');
        f(1) = subplot(2,3,4)
        plot(distances(1,counter), max(F.stim(1).meanArt), 'o', 'color', color); hold on;  title('350uA stim artifacts max value vs depth'); xlabel('Depth(mm)'); ylabel('voltage(uV)');
        if size(F.stim,2)>1
            if max(F.stim(2).meanArt)>maxArt
                maxArt = max(F.stim(2).meanArt);
            end
            if max(F.stim(2).meanArt)<minArt
                minArt = max(F.stim(2).meanArt);
            end
            % plot stim 2
            h(2)= subplot(2,3,2)
            plot([0:1/F.s_rate:(size(F.stim(2).meanArt,2)-1)/F.s_rate],F.stim(2).meanArt, 'color', cc(counter,:)); hold on; title('250uA stim mean artifacts'); xlabel('Time(ms)'); ylabel('voltage(uV)');
            f(2)= subplot(2,3,5)
            plot(distances(1,counter), max(F.stim(2).meanArt), 'o', 'color', color); hold on; title('250uA stim artifacts max value vs depth'); xlabel('Depth(mm)'); ylabel('voltage(uV)');
        end
        if size(F.stim,2)>2
            if max(F.stim(3).meanArt)>maxArt
                maxArt = max(F.stim(3).meanArt);
            end
            if max(F.stim(3).meanArt)<minArt
                minArt = max(F.stim(3).meanArt);
            end
            % plot stim 2
            h(3) = subplot(2,3,3)
            plot([0:1/F.s_rate:(size(F.stim(3).meanArt,2)-1)/F.s_rate],F.stim(3).meanArt, 'color', cc(counter,:)); hold on; title('150uA stim mean artifacts'); xlabel('Time(ms)'); ylabel('voltage(uV)');
            f(3) = subplot(2,3,6)
            plot(distances(1,counter), max(F.stim(3).meanArt), 'o', 'color', color); hold on; title('150uA stim artifacts max value vs depth'); xlabel('Depth(mm)'); ylabel('voltage(uV)');
        end
        
        
        
    end
end
xlim(f(1),[floor(distances(1,1))-1,  ceil(distances(1,end))+1]);
ylim(f(1),[minArt-1/10*minArt,  maxArt+1/10*maxArt]);
linkaxes(h,'xy');
linkaxes(f,'xy');

%%  Calculate the max amplitude and then store it back into F

clear all, clc, close all
cd('C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis\Nex_Files');
directory = 'C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis\Nex_Files';
listing = dir(directory);


for j = 1:size(listing,1)
    s = listing(j).name;
    if strfind(s,'proc')
        
        load(listing(j).name); % Now there's 'F'
        for i = 1:length(F.stim)
            meanArt = F.stim(i).meanArt;
            meanArtMaxAmp = max(meanArt); % obtain the max amplitude of the mean artifact
            F.stim(i).meanArtMaxAmp = meanArtMaxAmp; 
            save(s,'F'); % saves matlab output
        end
                
    end
end












