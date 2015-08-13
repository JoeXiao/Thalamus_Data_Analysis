clear,  close all, clc
[fn,pn]=uigetfile('.mat','Pick File','C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis');
load([pn,fn]);% pn: path name, fn: file name

% Parameters that can be changed
channels_to_use =[1]; % change this
thresh=4000; % set treshold
blanklength=30; % blanklength is the number of samples to be blanked for each artifact
SBAB = 15; % number of samples before artifact detection to be blanked
templ_length=10; % number of artifacts to be averaged into a template
SBAD = 20; % This is the number of samples before the artifact crossing the threshold to be included in the artifact waveform
artSeg = [2:12]; % segment of the artifact waveform to compare 
F=[]; % This is where the processed results are stored

for i=1:length(channels_to_use)
    chanNum = channels_to_use(i);
    if chanNum<10
        str=['0', num2str(chanNum)];
    else
        str=[num2str(chanNum)];
    end
    eval(sprintf('dataAll=CSPK_0%s;',str)); % 'dr' is data raw from the channel specified in str, will work only up to 9 channels,
    eval(sprintf('s_rate=round(CSPK_0%s_KHz*10000)/10;',str)); % sampling rate, rounded to one decimal to avoid misaligning with offline sorter (who does not accept too many decimals of Hz)
    % Free up memory
    clearvars -except pn fn channels_to_use thresh blanklength SBAB templ_length SBAD artSeg F chanNum str dataAll s_rate
    dataAll=double(dataAll); % convert to double 
    
    % signal "energy" for stimulation epochs detection
    dd=dataAll([1:floor(0.1*s_rate)]'*ones(1,floor(length(dataAll)/(0.1*s_rate)))...
        +ones(floor(0.1*s_rate),1)*[0:floor(length(dataAll)/(0.1*s_rate))-1]*floor(0.1*s_rate));
    % 0.1*s_rate is grouping data into 0.1s windows
    % This is a rearranging of the array into a matrix format, each column is a non-overlapping window of data that's 0.1s long
    stim_start = find(diff(max(dd)-min(dd)>thresh)==1)*floor(0.1*s_rate); % *floor(0.1*s_rate) gets back to the original sample number
    stim_stops = find(diff(max(dd)-min(dd)>thresh)==-1)*floor(0.1*s_rate); % *floor(0.1*s_rate) gets back to the original sample number
    goodones=(stim_stops-stim_start)/s_rate>20; % there needs to be 20s of stim to be considered 'good'
    stim_start=stim_start(goodones);
    stim_stops=stim_stops(goodones);
    epochs_start=[1,(stim_start(2:end)+stim_stops(1:end-1))/2,length(dataAll)]; % epoch includes the time before stim and the time after stim
    
    for jj=1:length(epochs_start)-1 % length(epochs_start - 1) is the total number of epochs
        sprintf(['processing epoch ', num2str(jj)])
        dataEpoch=dataAll(epochs_start(jj):epochs_start(jj+1)-1); % data within the current epoch
        
        dataAboveThr=find(diff([0,dataEpoch>thresh])==-1); % samples where the artifact first crosses the threshold
        stim_int=median(diff(dataAboveThr)); % inter stimulation period in samples
        
        arts=dataEpoch(-SBAD +dataAboveThr'*ones(1,stim_int)+ones(length(dataAboveThr),1)*[0:stim_int-1]); % extracts all artifacts waveforms
        
        mask=ones(1,length(dataEpoch)); % prepares blanking mask and stimulation timestamps
        m_pre=fliplr([dataAboveThr(1):-stim_int:3])'/s_rate; % marks the starting sample point of the mask
        m_on=dataAboveThr'/s_rate;
        m_post=[dataAboveThr(end)+stim_int:stim_int:epochs_start(jj+1)-epochs_start(jj)-1]'/s_rate;
        m=[fliplr([dataAboveThr(1):-stim_int:3]),dataAboveThr,[dataAboveThr(end)+stim_int:stim_int:length(dataEpoch)-stim_int-1]]; % all the stim start time stamps, including virtual stim
        mask=ones(1,length(dataEpoch));
        if m(1)<5 
            m=m(2:end); 
        end
        if m(end)>length(dataEpoch)-5 
            m=m(1:end-1); 
        end
        
        mask(m'*ones(1,blanklength+1)+ones(length(m),1)*[-SBAB:blanklength-SBAB])=0; % SBAB is the number of samples before art detection to be blanked
        % treatment of data is uniform        
        for ii=1:length(dataAboveThr)
            [Y,I]=sort(sum(abs(ones(length(dataAboveThr),1)*arts(ii,artSeg)-arts(:,artSeg)),2)); % find best artifact matches, artSeg is the segment of the pulse to compare
            t=arts(I(randperm(3*templ_length)+1),:); % create template
            dataEpoch(-SBAD+dataAboveThr(ii):-SBAD+dataAboveThr(ii)+stim_int-1)=dataEpoch(-SBAD+dataAboveThr(ii):-SBAD+dataAboveThr(ii)+stim_int-1)-mean(t); % subtract template from original data
        end
        % applies blanking to data in this epoch
        dataEpoch=dataEpoch.*mask;
        % replace subtracted and blanked data into data array
        dataAll(epochs_start(jj):epochs_start(jj+1)-1)=dataEpoch; 
        % archive pulse timestamps for each stimulation epoch
        F.s_rate=s_rate;
        F.stim(jj).interval=[epochs_start(jj),epochs_start(jj+1)-1];
        F.stim(jj).stim_int=stim_int;
        F.stim(jj).pr=m_pre+epochs_start(jj)/s_rate;
        F.stim(jj).on=m_on+epochs_start(jj)/s_rate;
        F.stim(jj).po=m_post+epochs_start(jj)/s_rate;
        F.stim(jj).meanArt = mean(arts); % This is the mean stimulation artifact in this epoch
        meanArt = mean(arts);
        x = 0:length(meanArt)-1;
        y = meanArt;
        xx = 0:0.01:length(meanArt);
        yy = spline(x,y,xx); % fitting a spline over the stim artifact data
        F.stim(jj).splineMeanArt = yy;
        F.stim(jj).splineMeanArtMax = max(yy);
        
    end
    eval(sprintf('F.d_sub_C%d=dataAll;',chanNum)); % archives subtracted data for Channel chanNum
    dsn=254*255/2*dataAll/400;%(std(dataAll)*5); % rescale subtracted data to better use 16 bit encoding in binary file
    n(1:2:2*length(dsn)-1)=floor(double(dsn)/255); % first 8 bit of every sample (big part of the value)
    n(2:2:2*length(dsn))=rem(double(dsn),255); % second 8 bit of every sample (little endian coding)
    fid=fopen([pn,fn(1:end-4),sprintf('_C%d.2plx',chanNum)],'w+'); % save binary file output to plexon
    fwrite(fid,n,'int8');
    fclose(fid);
    save([pn,fn(1:end-4),'_proc.mat'],'F'); % saves matlab output
end

plot(max(dd)-min(dd))
plot(arts')
figure; plot(dataAll);
figure;
if size(F.stim,2)==1;
    plot(F.stim(1).meanArt,'r'); legend('stim1');
elseif size(F.stim,2) ==2
    plot(F.stim(1).meanArt,'r'); hold on; plot(F.stim(2).meanArt,'b'); legend('stim1','stim2');
else
    plot(F.stim(1).meanArt,'r'); hold on; plot(F.stim(2).meanArt,'b'); plot(F.stim(3).meanArt,'g'); legend('stim1','stim2','stim3');
end

%%  adding timestamps
clear all, clc
[fn,pn]=uigetfile('.mat','Pick File');
load([pn,fn]);
list=fieldnames(F);

for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))>=7
        if min(name(1:7)=='d_sub_C')
            try
                ff=readNexFile([pn,fn(1:end-9),sprintf('_%s.nex',name(end-1:end))]);
                for jj=1:length(ff.neurons)
                    eval(sprintf('F.C%s_n%d=ff.neurons{%d}.timestamps;',name(end-1:end),jj,jj));
                end
                for jj=1:length(F.stim)
                    ff=nexAddEvent(ff,F.stim(jj).pr,sprintf('s_bef%d',jj));
                    ff=nexAddEvent(ff,F.stim(jj).on,sprintf('s_dur%d',jj));
                    ff=nexAddEvent(ff,F.stim(jj).po,sprintf('s_pos%d',jj));
                end
                writeNexFile(ff,[pn,fn(1:end-9),sprintf('_%s.nex',name(end-1:end))]);
            catch
                disp([pn,fn(1:end-9),sprintf('_%s.nex',name(end-1:end)), ' has no nex file'])
            end
            
        end
    end
end
save([pn,fn(1:end-4)],'F');

%%
% PSTH confidence interval analysis
clear
bin=0.2; % ms
[fn,pn]=uigetfile('.mat','Pick File');
load([pn,fn]);

list=fieldnames(F);
for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))==4
        if name(1)=='C' && name(end-2)=='_' && name(end-1)=='n' % these 2 ifs select fields of F corresponding to cell spikes timestamps
            eval(['spks=F.',name,';'])
            for jj=1:length(F.stim)
                %select spikes within the stimulation jj
                pr_spks=spks(spks>F.stim(jj).pr(1) & spks<F.stim(jj).pr(end));
                on_spks=spks(spks>F.stim(jj).on(1) & spks<F.stim(jj).on(end));
                po_spks=spks(spks>F.stim(jj).po(1) & spks<F.stim(jj).po(end));
                
                % this calculates the distance between all spikes and the
                % immediately preceding stimualtion pulse
                pr_st_sp_int=pr_spks*ones(1,length(F.stim(jj).pr))-ones(length(pr_spks),1)*F.stim(jj).pr'; % these 3 lines calculates the distance between every spike and every pulse (spike time - pulse time)
                on_st_sp_int=on_spks*ones(1,length(F.stim(jj).on))-ones(length(on_spks),1)*F.stim(jj).on';
                po_st_sp_int=po_spks*ones(1,length(F.stim(jj).po))-ones(length(po_spks),1)*F.stim(jj).po';
                pr_st_sp_int=min(pr_st_sp_int+(pr_st_sp_int<0)*10^6,[],2); % these 3 lines substitutes negative values (corresponding to pulses following the spike) with very large positive numbers so to select the minimum positive distance for every spike (the distance between the spike and the immediately preceding pulse)
                on_st_sp_int=min(on_st_sp_int+(on_st_sp_int<0)*10^6,[],2);
                po_st_sp_int=min(po_st_sp_int+(po_st_sp_int<0)*10^6,[],2);
                
                % this finds which bin the distance between the spike and
                % the immediately preceding pulse belongs to, sums across
                % bins to find the number of spikes belonging to that bin
                % and normalize to spiking rate
                psth(1,:)=sum(-diff(pr_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(pr_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).pr));
                psth(2,:)=sum(-diff(on_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(on_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).on));
                psth(3,:)=sum(-diff(po_st_sp_int*ones(1,length([0:bin/1000:F.stim(jj).stim_int/F.s_rate]))>ones(length(po_st_sp_int),1)*[0:bin/1000:F.stim(jj).stim_int/F.s_rate],1,2))/(bin/1000*length(F.stim(jj).po));
                %bins during stimulation that cross mean + 3 std of bin of the pre-stim period
                psth_excitation=((psth(2,:)>mean(psth(1,:))+3*std(psth(1,:))));
                psth_inhibition=((psth(2,:)<mean(psth(1,:))-3*std(psth(1,:))));
                
                exc_intervals=[find(diff(psth_excitation)==1)'+1,find(diff(psth_excitation)==-1)'+1];
                psth_maxes=[];
                for kk=1:size(exc_intervals,1)
                    temp_maxes=(exc_intervals(kk,1)+find(psth(2,[exc_intervals(kk,1):exc_intervals(kk,2)])==max(psth(2,[exc_intervals(kk,1):exc_intervals(kk,2)])))-1)*bin;
                    psth_maxes=[psth_maxes,temp_maxes(1)];
                end
                
                ini_intervals=[find(diff(psth_inhibition)==1)',find(diff(psth_inhibition)==-1)'];
                psth_mins=[];
                for kk=1:size(ini_intervals,1)
                    temp_mins=(ini_intervals(kk,1)+find(psth(2,[ini_intervals(kk,1):ini_intervals(kk,2)])==max(psth(2,[ini_intervals(kk,1):ini_intervals(kk,2)])))-1)*bin;
                    psth_mins=[psth_mins,temp_mins(1)];
                end
                
                %place data within the file structure
                eval(sprintf(['F.stim(jj).psth_%s_n',name(end),'=psth;'],name(1:find(name=='_')-1)));
                eval(sprintf(['F.stim(jj).psth_ex_%s_n',name(end),'=exc_intervals;'],name(1:find(name=='_')-1)));
                eval(sprintf(['F.stim(jj).psth_peaks_%s_n',name(end),'=psth_maxes;'],name(1:find(name=='_')-1)));
                eval(sprintf(['F.stim(jj).psth_in_%s_n',name(end),'=ini_intervals;'],name(1:find(name=='_')-1)));
                eval(sprintf(['F.stim(jj).psth_valleys_%s_n',name(end),'=psth_mins;'],name(1:find(name=='_')-1)));
                
                %place data in excel
                
                
                hf=figure;
                plot([0:bin:F.stim(jj).stim_int/F.s_rate*1000-bin],psth(2,:),'r', 'linewidth',2);hold on
                plot([0:bin:F.stim(jj).stim_int/F.s_rate*1000-bin],psth(1,:),'b', 'linewidth',2);hold on
                
                ylabel('Firing Frequency (Hz)'),xlabel('Time (msec)'),title([fn,' ',name]);
                plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:)),mean(psth(1,:))],'k--', 'linewidth',2)
                plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))+3*std(psth(1,:)),mean(psth(1,:))+3*std(psth(1,:))],'k-.', 'linewidth',2)
                plot([0,F.stim(jj).stim_int/F.s_rate*1000-bin],[mean(psth(1,:))-3*std(psth(1,:)),mean(psth(1,:))-3*std(psth(1,:))],'k-.', 'linewidth',2)
                
                
                plot([(find((psth(2,:)>mean(psth(1,:))+3*std(psth(1,:))))-1)*bin], 1.1*max(psth(2,:)),'k*')
                saveas(gcf,[pn,sprintf([fn,'_ch',name(2:find(name=='_')-1),'_cell',name(find(name=='n')+1:end),'_stim%d.jpg'],jj)],'jpg');
            end
        end
    end
end
save([pn,fn(1:end-4)],'F');



% psth plot
clear
cts=['C0';'C1';'C2';'C3';'C4';'C5';'C6']; cts=flipud(cts);
[fn,pn]=uigetfile('.mat','Pick File');
load([pn,fn]);
list=fieldnames(F);
for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))==4
        if sum(name(1:end-1)=='C0_n')==3
            eval(['spks=F.',name,';'])
            for jj=1:length(F.stim)
                pr_spks=spks(spks>F.stim(jj).pr(1) & spks<F.stim(jj).pr(end));
                on_spks=spks(spks>F.stim(jj).on(1) & spks<F.stim(jj).on(end));
                po_spks=spks(spks>F.stim(jj).po(1) & spks<F.stim(jj).po(end));
                
                pr_st_sp_int=pr_spks*ones(1,length(F.stim(jj).pr))-ones(length(pr_spks),1)*F.stim(jj).pr';
                on_st_sp_int=on_spks*ones(1,length(F.stim(jj).on))-ones(length(on_spks),1)*F.stim(jj).on';
                po_st_sp_int=po_spks*ones(1,length(F.stim(jj).po))-ones(length(po_spks),1)*F.stim(jj).po';
                pr_st_sp_int=min(pr_st_sp_int+(pr_st_sp_int<0)*10^6,[],2);
                on_st_sp_int=min(on_st_sp_int+(on_st_sp_int<0)*10^6,[],2);
                po_st_sp_int=min(po_st_sp_int+(po_st_sp_int<0)*10^6,[],2);
                
                psth(1,:)=hist(pr_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).pr)/F.s_rate);
                psth(2,:)=hist(on_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).on)/F.s_rate);
                psth(3,:)=hist(po_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).po)/F.s_rate);
                f1=figure;
                plot([1:100]/100*(F.stim(jj).stim_int/F.s_rate)*1000,psth')
                xlabel('Time (ms)'); ylabel('Firing Frequency (Hz)');
                title(sprintf([fn,' ch %d, cell %d, %d Hz'],str2num(name(2)),str2num(name(4)),round(F.s_rate/F.stim(jj).stim_int)));
                legend ('before', 'during', 'after')
                saveas(gcf,[pn,sprintf([fn,'_ch%d_cell%d_%dHz_%s'],str2num(name(2)),str2num(name(5)),round(F.s_rate/F.stim(jj).stim_int),cts(jj,:)),'.jpg'],'jpg');
                eval(sprintf(['F.stim(jj).psth_n',name(end),'=psth;']));
            end
        end
    end
end
save([pn,fn(1:end-4)],'F');
% close all






clear
[fn,pn]=uigetfile('.mat','Pick File');
load([pn,fn]);
list=fieldnames(F);
for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))==5
        if name(1)=='C' && name(end-2)=='_' && name(end-1)=='n' % sum(name(1:end-1)=='C0_n')==3
            eval(['spks=F.',name,';'])
            
            for jj=1:length(F.stim)
                pr_spks=spks(spks>F.stim(jj).pr(1) & spks<F.stim(jj).pr(end));
                on_spks=spks(spks>F.stim(jj).on(1) & spks<F.stim(jj).on(end));
                po_spks=spks(spks>F.stim(jj).po(1) & spks<F.stim(jj).po(end));
                
                pr_st_sp_int=pr_spks*ones(1,length(F.stim(jj).pr))-ones(length(pr_spks),1)*F.stim(jj).pr';
                on_st_sp_int=on_spks*ones(1,length(F.stim(jj).on))-ones(length(on_spks),1)*F.stim(jj).on';
                po_st_sp_int=po_spks*ones(1,length(F.stim(jj).po))-ones(length(po_spks),1)*F.stim(jj).po';
                pr_st_sp_int=min(pr_st_sp_int+(pr_st_sp_int<0)*10^6,[],2);
                on_st_sp_int=min(on_st_sp_int+(on_st_sp_int<0)*10^6,[],2);
                po_st_sp_int=min(po_st_sp_int+(po_st_sp_int<0)*10^6,[],2);
                
                psth(1,:)=hist(pr_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).pr)/F.s_rate);
                psth(2,:)=hist(on_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).on)/F.s_rate);
                psth(3,:)=hist(po_st_sp_int,100)/(F.stim(jj).stim_int/100*length(F.stim(jj).po)/F.s_rate);
                eval(sprintf(['F.stim(jj).psth_n',name(end),'=psth;']));
                f1=figure;
                plot([1:100]/100*(F.stim(jj).stim_int/F.s_rate)*1000,psth')
                xlabel('Time (ms)'); ylabel('Firing Frequency (Hz)');
                title(sprintf([fn,' ch %d, cell %d, %d Hz'],str2num(name(2)),str2num(name(4)),round(F.s_rate/F.stim(jj).stim_int)));
                legend ('before', 'during', 'after')
                saveas(gcf,[pn,sprintf([fn,'_ch%d_cell%d_%dHz'],str2num(name(2)),str2num(name(5)),round(F.s_rate/F.stim(jj).stim_int)),'.jpg'],'jpg');
                eval(sprintf(['F.stim(jj).psth_n',name(end),'=psth;']));
                
            end
            save([pn,fn(1:end-4)],'F');
        end
    end
end





%% different contact stims plot of psth during dbs for all contacts
clear
cts=['C0';'C1';'C2';'C3';'C4';'C5';'C6'];
order=1;
[fn,pn]=uigetfile('.mat','Pick File','C:\data\Kramer\Recording\13_09_13');
load([pn,fn]);
list=fieldnames(F);
f1=figure;
hold all
for ii=1:length(list)
    name=char(list(ii));
    if length(name(1:end-1))==4
        if sum(name(1:end-1)=='C0_n')==3
            
            
            for jj=1:length(F.stim)
                eval(sprintf(['psth=F.stim(jj).psth_n',name(end),'(2,:);']));
                plot([1:100]/100*(F.stim(jj).stim_int/F.s_rate)*1000,psth')
            end
            
            xlabel('time (msec)'),ylabel('firing rate (Hz)')
            title(fn)
            if order==1
                legend(cts)
            else
                legend(flipud(cts))
            end
            saveas(gcf,[pn,sprintf([fn,'_ch%d_cell%d_%dHz_all contacts'],str2num(name(2)),str2num(name(5)),round(F.s_rate/F.stim(jj).stim_int)),'.jpg'],'jpg');
            
        end
    end
end
