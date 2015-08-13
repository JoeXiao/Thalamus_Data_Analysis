% This code determines if a neuron is modulated by stimulation (firing rate and firing pattern)

%  Load the sorted spike times
clear all, close all, clc
bin=0.1; % (ms) - bin size in PSTH
pplus_std = 6; % number of std above the mean to be considered exicitatory pattern modulation
pplus_bin = 2; % number of consecutive bins that's above (mean+pplus_std*std)of the pre-stim period needed to be considered exci pattern modulation
pminus_std = 2; % number of std below the mean to be considered inhibitory pattern modulation
pminus_bin = 8; % number of consecutive bins that's below (mean-pminus_std*std)of the pre-stim period needed to be considered inhib pattern modulation
rplus_std = 2; % number of std above the mean to be considered exicitatory firing rate modulation
rminus_std = 1.5; % number of std below the mean to be considered inhibitory firing rate modulation
beginChop = 0.5; % (ms) - amount of time at the beginning of the PSTH (DBS-ON) to chop out
endChop = 0.5; % (ms) - amount of time at the end of the PSTH (DBS_ON) to chop out
maxDispFreq = 70; % (Hz) - maximum value on the y-axis when displaying the firing rate of the PSTH 

[fn,pn]=uigetfile('.mat','Pick File'); % Get the proc file
load([pn,fn]);
list=fieldnames(F);

% Extract out the time stamps of all the sorted cells
for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))>=7
        if min(name(1:7)=='d_sub_C')
            try
                ff=readNexFile([pn,fn(1:end-9),sprintf('_%s.nex',name(end-1:end))]); % Read the Nex file
                
            catch
                disp([pn,fn(1:end-9),sprintf('_%s.nex',name(end-1:end)), ' has no nex file'])
            end
            
        end
    end
end

for cellNum = 1: size(ff.neurons,1) % Loop through each sorted cell in the file
    spks = ff.neurons{cellNum}.timestamps; % Extract out the time stamps of the sorted cells
    figure;
    stimNum = length(F.stim); % Number of stimulations. Group the spikes into pre, dur and post periods
    for jj=1:length(F.stim)
        isModulation = [0, 0, 0, 0]; % p+,p-,r+,r-
        % select spikes within the stimulation jj
        pr_spks=spks(spks>F.stim(jj).pr(1) & spks<F.stim(jj).pr(end));
        on_spks=spks(spks>F.stim(jj).on(1) & spks<F.stim(jj).on(end));
        po_spks=spks(spks>F.stim(jj).po(1) & spks<F.stim(jj).po(end));
        % for each spike, finds the time difference to all the pulses in the
        % same epoch
        pr_st_sp_int=pr_spks*ones(1,length(F.stim(jj).pr))-ones(length(pr_spks),1)*F.stim(jj).pr'; % these 3 lines calculates the distance between every spike and every pulse (spike time - pulse time)
        on_st_sp_int=on_spks*ones(1,length(F.stim(jj).on))-ones(length(on_spks),1)*F.stim(jj).on';
        po_st_sp_int=po_spks*ones(1,length(F.stim(jj).po))-ones(length(po_spks),1)*F.stim(jj).po';
        % Find the time difference to the pulse immediately preceding it, i.e.
        % find the smallest non-zero term
        pr_st_sp_int=min(pr_st_sp_int+(pr_st_sp_int<0)*10^6,[],2); % these 3 lines substitutes negative values (corresponding to pulses following the spike) with very large positive numbers so to select the minimum positive distance for every spike (the distance between the spike and the immediately preceding pulse)
        on_st_sp_int=min(on_st_sp_int+(on_st_sp_int<0)*10^6,[],2);
        po_st_sp_int=min(po_st_sp_int+(po_st_sp_int<0)*10^6,[],2);
        % Chop out data in the DBS-on period that are smaller than
        % 'beginChop' and larger than 'endChop'
        pr_st_sp_int((pr_st_sp_int<=beginChop/1000) | (pr_st_sp_int>=(F.stim(jj).stim_int/F.s_rate-endChop/1000))) = [];
        on_st_sp_int((on_st_sp_int<=beginChop/1000) | (on_st_sp_int>=(F.stim(jj).stim_int/F.s_rate-endChop/1000))) = [];
        po_st_sp_int((po_st_sp_int<=beginChop/1000) | (po_st_sp_int>=(F.stim(jj).stim_int/F.s_rate-endChop/1000))) = [];
        % this finds which bin the distance between the spike and the immediately preceding pulse belongs to, sums across
        % bins to find the number of spikes belonging to that bin and
        % normalize to spiking rate, in Hz. Note: because of rounding of
        % F.stim(jj).stim_int/F.s_rate, the last bin will be lost
        psth(1,:)=sum(-diff(pr_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(pr_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).pr)); % length(F.stim(jj).pr) is the number of virtual stim periods
        psth(2,:)=sum(-diff(on_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(on_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).on));
        psth(3,:)=sum(-diff(po_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(po_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).po));
        
        % Determine if there's firing rate modulation
        if mean(psth(2,:)) > (mean(psth(1,:)) + rplus_std*std(psth(1,:))) % There's significant firing rate increase
            isModulation(3) = 1;
        end
        if mean(psth(2,:)) < (mean(psth(1,:)) - rminus_std*std(psth(1,:))) % There's significant firing rate decrease
            isModulation(4) = 1;
        end
        % Find the bins of p+,p-
        psth_excitation=((psth(2,:) > mean(psth(1,:))+ pplus_std*std(psth(1,:))));
        psth_inhibition=((psth(2,:) < mean(psth(1,:))- pminus_std*std(psth(1,:))));
        % Calculate the maximum number of consecutive excitatory or
        % inhibitory bins
        Nmax_excit = maxConsecOnes(psth_excitation);
        Nmax_inhib = maxConsecOnes(psth_inhibition);
        if Nmax_excit >= pplus_bin % consecutive bin criteria for excitation is met
            isModulation(1) = 1;
        end
        if Nmax_inhib >= pminus_bin % consecutive bin criteria for inhibition is met
            isModulation(2) = 1;
        end
        % Plot the results
        h(jj) = subplot(stimNum, 1, jj);
        plot([0:bin:F.stim(jj).stim_int/F.s_rate*1000-bin],psth(1,:),'b', 'linewidth',2);hold on
        plot([0:bin:F.stim(jj).stim_int/F.s_rate*1000-bin],psth(2,:),'r', 'linewidth',2);hold on
        plot([0:bin:F.stim(jj).stim_int/F.s_rate*1000-bin],psth(3,:),'g', 'linewidth',2);hold on
        legend('before', 'during', 'after'); ylabel('Firing Rate (Hz)'); xlabel('Time (msec)');
        % Describe the results in the title of each subplot 
        if sum(isModulation == 0) % no change at all 
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': n'  '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(1) & isModulation(3)
            title([fn,' ', 'cell_', num2str(cellNum)  '_stim_' num2str(jj), ': p+r+' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(1) & isModulation(4)
            title([fn,' ', 'cell_', num2str(cellNum)  '_stim_' num2str(jj), ': p+r-' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(1)
            title([fn,' ', 'cell_', num2str(cellNum)  '_stim_' num2str(jj), ': p+' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(3)
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': r+' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(2) & isModulation(4)
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': p-r-' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(2) & isModulation(3)
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': p-r+' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(2)
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': p-' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        elseif isModulation(4)
            title([fn,' ','cell_', num2str(cellNum), '_stim_' num2str(jj), ': r-' '    cell count' ' ' num2str(length(spks))],'interpreter','none');
        end
        % Plot out lines to indicate mean +/- std
        plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))+ rplus_std*std(psth(1,:)),mean(psth(1,:))+ rplus_std*std(psth(1,:))],'k-.', 'linewidth',2)
        plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))- rminus_std*std(psth(1,:)),mean(psth(1,:))- rminus_std*std(psth(1,:))],'k-.', 'linewidth',2)
        plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))+ pplus_std*std(psth(1,:)),mean(psth(1,:))+ pplus_std*std(psth(1,:))],'r-.', 'linewidth',2)
        plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))- pminus_std*std(psth(1,:)),mean(psth(1,:))- pminus_std*std(psth(1,:))],'r-.', 'linewidth',2)
        % if exitatory pattern modulation is present, locate those that
        % crossed the threshold 
        if isModulation(1)
            plot([(find((psth(2,:) > mean(psth(1,:))+ pplus_std*std(psth(1,:))))-1)*bin], 1.1*max(psth(2,:)),'r*') 
        end
        % if inhibitory pattern modulation is present, locate those that
        % crossed the threshold 
        if  isModulation(2)
            plot([(find((psth(2,:) < mean(psth(1,:))- pminus_std*std(psth(1,:))))-1)*bin], 1.3*max(psth(2,:)),'c*') 
        end
        % Store the relevant data 
        F.Neuron(cellNum).stim(jj).PSTH.pr_psth = psth(1,:);
        F.Neuron(cellNum).stim(jj).PSTH.on_psth = psth(2,:);
        F.Neuron(cellNum).stim(jj).PSTH.po_psth = psth(3,:);
        F.Neuron(cellNum).stim(jj).PSTH.pr_psth_mean = mean(psth(1,:));
        F.Neuron(cellNum).stim(jj).PSTH.pr_psth_std = std(psth(1,:));
        F.Neuron(cellNum).stim(jj).PSTH.on_psth_mean = mean(psth(2,:));
        F.Neuron(cellNum).stim(jj).PSTH.on_psth_std = std(psth(2,:));
        F.Neuron(cellNum).stim(jj).PSTH.po_psth_mean = mean(psth(3,:));
        F.Neuron(cellNum).stim(jj).PSTH.po_psth_std = std(psth(3,:));
        F.Neuron(cellNum).stim(jj).PSTH.Delta_mean = mean(psth(2,:)) - mean(psth(1,:)); % diff in mean FR btw DBS-on and pre-DBS
        F.Neuron(cellNum).cellCount = length(spks);
        disp(['Cell' ' ' num2str(cellNum) '---' ' ' 'Cell count is:' ' ' num2str(length(spks))])
        disp(['Stim' ' ' num2str(jj) ' Max artifact amplitude is:' ' ' num2str(F.stim(jj).splineMeanArtMax)])
        disp(['Stim' ' ' num2str(jj) ' DBS mean firing rate:' ' ' num2str(F.Neuron(cellNum).stim(jj).PSTH.on_psth_mean)])
        disp(['Stim' ' ' num2str(jj) ' Difference in DBS mean firing rate:' ' ' num2str(F.Neuron(cellNum).stim(jj).PSTH.Delta_mean)])
        xlim([0, F.stim(jj).stim_int/F.s_rate*1000-bin]);
    end
    linkaxes(h); ylim([0 maxDispFreq]);
    saveas(gcf,[pn,sprintf([fn '_cell' num2str(cellNum) '_PSTH.jpg'],jj)],'jpg');
    saveas(gcf,[pn,sprintf([fn '_cell' num2str(cellNum) '_PSTH.fig'],jj)],'fig');
    save([pn,fn],'F');
end

