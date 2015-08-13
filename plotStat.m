

Data = [VarName1 pr{:} n{:} n1{:}];

close all
Data = a+b+c;
Data = Data';

bar(Data,'stack')
title ('All_Data_Steering', 'interpreter', 'none')


set(gca,'XTickLabel', {'39', '31', '15','7'});
legend('pr', 'p','r', 'n','r-');


%% calculate percentages of the stacked bar plot
clear all, clc
Data = [12	15	15	16
    4	4	1	1
    9	6	11	6
    9	11	12	12
    8	6	6	4
    
    ];
percent = zeros(5,4);
for i = 1:4
    percent(:,i) = Data(:,i)/sum(Data(:,i));
end
percent


%% Pie Chart for each site
clear all, close all, clc
pieData = [8		5		6		5
1		2		1		1
6		10		11		6
9		7		5		12
7		7		8		4
0		0		0		1
0		0		0		2

    ];
name = 'All_Data';
str = {'pr: ';'p: '; 'r: '; 'n: '; 'r-: ';'p-: '; 'p-r-: '}; % text
for i = 1:4
    
    values = pieData(:,i);
    
    
    emptyIndices = find(values ==0);
    % clean out the labels that do not exist
    str_new = str;
    str_new(emptyIndices)=[];
    subplot(1,4,i)
    h = pie(values);
    
    
    hp = findobj(h, 'Type', 'patch');
    index = find(values); % find the nonzero indices
    for j = 1:length(hp)
        if index(j) == 1
            
            set(hp(j), 'FaceColor', 'g');
        elseif index(j) == 2
            set(hp(j), 'FaceColor', 'r');
        elseif index(j) == 3
            set(hp(j), 'FaceColor', 'b');
        elseif index(j) == 4
            set(hp(j), 'FaceColor', 'm');
        elseif index(j) == 5
            set(hp(j), 'FaceColor', 'c');
        elseif index(j) == 6
            set(hp(j), 'FaceColor', 'y');
        else
            set(hp(j), 'FaceColor', 'k');
        end
    end
    
    
    hText = findobj(h,'Type','text'); % text handles
    percentValues = get(hText,'String'); % percent values
    combinedstrings = strcat(str_new,percentValues);
    %     oldExtents_cell = get(hText,'Extent'); % cell array
    %     oldExtents = cell2mat(oldExtents_cell); % numeric array
    set(hText,{'String'},combinedstrings);
    switch i
        case 1
            title([name '_350uA'],'interpreter', 'none');
        case 2
            title([name '_250uA'],'interpreter', 'none');
        case 3
            title([name '_150uA'],'interpreter', 'none');
    end
end

%% Calculate the average mean and std of average firing rate (DBS period) across 4 contacts
clear all, clc
a = [147.9165
    152.2184
    24.9489
    30.8237
    
    ];
Artifacts = [-29.937/22.296
    -36.0327/134.6642
    -12.974
    -14.4172
    
    ];
CV = std(a)/mean(a)
avgArtifact = mean(Artifacts)

[curve, goodness, output] = fit(artifacts,CV,'poly1');
plot(artifacts, CV,'o'); hold on;
plot(curve)
goodness.rsquare
title(['CV vs Peak Artifact Amplitude', ' r^2 = ', ' ' num2str(goodness.rsquare) ])
xlabel(['Peak Artifact Amplitude (uV)'])
ylabel(['Coefficient of Variability'])
%% Plot the resultant vector
clear all, clc, close all
% Site39 = [1 0];
% Site31 = [0,1];
% Site15 = [-1,0];
% Site7 = [0, -1];
Site_Vecs = [1 0; 0 1; -1 0; 0 -1];


Data = [8.9059	0.021939	3.4451	0.72552	12722.6054	12880.397	13425.1076	13908.9934
-2.2082	-3.1809	-3.6216	-1.6723	13171.4931	13377.493	13625.1772	14660.2923
5.1656	-4.0308	-2.2791	-4.6146	14319.509	13688.8535	13836.5629	14489.5457
8.9997	2.2171	-1.7391	4.2939	14660.6322	14396.8434	14491.4475	15634.6272
-2.0909	-2.1267	-0.12763	5.6572	14083.1231	13578.2576	13606.2403	14715.3761
4.6768	-2.4021	12.9147	40.8415	14695.3087	13920.6572	13850.5005	14804.9805
1.2714	2.2115	-0.30259	2.093	14616.3175	13836.1792	13746.5977	15107.3018
-6.225	-1.4701	-2.8244	6.4171	14970.2994	13966.2677	13961.6097	15421.2955
]; % The first 4 elements are the D-FR for 39,31,15 and 7, the last 4 elements are the peak amplitudes for these 4 sites


Result_Data = [24.7203	-7.289961	8.28978	47.32422	113239.2884	109644.9486	110543.2434	118742.4125
];


% Calculate the resultant vector
Result_D_FR = repmat(Result_Data(1:4)',1,2);
Result_Arts = repmat(Result_Data(5:end)',1,2);
Result_D_FR_Vec = Site_Vecs.*Result_D_FR;
Result_Arts_Vec = Site_Vecs.*Result_Arts;
Result_D_FR_Vec = sum(Result_D_FR_Vec,1);
Result_Arts_Vec = sum(Result_Arts_Vec,1);
Resultant_Angle_Degrees  = acos(dot(Result_D_FR_Vec/norm(Result_D_FR_Vec),Result_Arts_Vec/norm(Result_Arts_Vec)))*180/pi

% Randomize the D_FR data for each row, and then obtain the resultant
% vector
% Randomized_Data = Data(:,1:4);
% for j = 1:size(Randomized_Data,1)
%     idx = randperm(4);
%     Randomized_Data(j,:) = Randomized_Data(j,idx);
%
% end
% Randomized_D_FR = sum(Randomized_Data,1);
% Randomized_D_FR = repmat(Randomized_D_FR',1,2);
% Randomized_Vec = Site_Vecs.*Randomized_D_FR;
% Randomized_Vec = sum(Randomized_Vec,1);
%
% Randomized_Angle_Degrees  = acos(dot(Randomized_Vec/norm(Randomized_Vec),Result_Arts_Vec/norm(Result_Arts_Vec)))*180/pi
%
%
% Randomized_Vec =  [Randomized_Vec/norm(Randomized_Vec) 0];
Result_D_FR_Vec = [Result_D_FR_Vec/norm(Result_D_FR_Vec) 0];
Result_Arts_Vec = [Result_Arts_Vec/norm(Result_Arts_Vec) 0];
% Result_D_FR_Vec = [Result_D_FR_Vec 0];
% Result_Arts_Vec = [Result_Arts_Vec 0];


% find the minimum of the mean amplitudes of 4 contacts
minAmp = min(mean(Data(:,5:end),2));
for i = 1:size(Data, 1)
    
    D_FR = repmat(Data(i,1:4)',1,2); % D_FR for all 4 sites
    Arts = repmat(Data(i,5:end)',1,2); % the peak artifact for all 4 sites
    FR_Vec = Site_Vecs.*(D_FR);
    Art_Vec =Site_Vecs.*(Arts);
    
    FR_Vec = sum(FR_Vec, 1);
    Art_Vec = sum(Art_Vec, 1);
    
    FR_Vec = [FR_Vec/norm(FR_Vec) 0];
    Art_Vec = [Art_Vec/norm(Art_Vec) 0];
    %     FR_Vec = [FR_Vec 0];
    %     Art_Vec = [Art_Vec 0];
    
    
    % a = [0 2];
    % b = [2 0];
    % c = a+b;
    % c = c/norm(c)
    % starts = zeros(5,5);
    % ends = [a;b;c];
    
    
    
    starts = repmat([0,0,mean(Data(i,5:end))/minAmp],2,1); % The z start point is the mean of the amplitudes for the 4 sites divided by the smallest mean amplitude from all the data
    ends = [FR_Vec; Art_Vec];
    
    quiver3(starts(1,1), starts(1,2), starts(1,3), ends(1,1), ends(1,2), ends(1,3),'color','b'); hold on
    quiver3(starts(2,1), starts(2,2), starts(2,3), ends(2,1), ends(2,2), ends(2,3),'color','r');
    
    
    l = legend('dFR_Resultant_Direction', 'Artifact_Result_Direction');
    set(l, 'Interpreter', 'none')
    % axis equal
    
    angle_Degrees  = acos(dot(FR_Vec/norm(FR_Vec),Art_Vec/norm(Art_Vec)))*180/pi
    
end


Result_Start = [ 0 0 0.99; 0,0, 0.99];
Result_End = [Result_D_FR_Vec; Result_Arts_Vec];
% Randomized_End = [Randomized_Vec];
quiver3(Result_Start(1,1), Result_Start(1,2), Result_Start(1,3),Result_End(1,1), Result_End(1,2), Result_End(1,3),'color','b', 'LineWidth', 5); hold on
quiver3(Result_Start(2,1), Result_Start(2,2), Result_Start(2,3),Result_End(2,1), Result_End(2,2), Result_End(2,3),'color','r', 'LineWidth', 5);
% quiver3(Result_Start(2,1), Result_Start(2,2), Result_Start(2,3),Randomized_End(1,1), Randomized_End(1,2), Randomized_End(1,3),'color','y', 'LineWidth', 5);
axis fill
title(['4.29_All_Data_Abs(dFR)_Artifact_Vector', ' n = ', num2str(4*size(Data,1)), ' Resultant_Angle = ' num2str(Resultant_Angle_Degrees)],'interpreter', 'none')
view(-40.5,12)
%% Plot the resultant vector in 3D, the DFR and Amplitude vectors are not unit vectors, but normalized to the smallest vector. 
clear all, clc, close all
Site_Vecs = [1 0; 0 1; -1 0; 0 -1];


Data = [-0.55485	-0.1101	0.041568	0.11164	6519.1874	6813.8728	6772.2209	7398.8161
-1.9377	5.0343	-1.3988	-3.1807	7076.235	6829.5423	7226.7423	7667.0715
1.244	-1.107	-2.7152	-1.9529	7472.5731	7303.1193	7407.7906	7228.5683
-0.2361	1.6402	-0.064872	3.9769	7835.4802	7677.9692	7715.1365	8349.1985
0.93346	0.084175	1.051	-2.3804	7685.2436	7203.2085	7416.591	7467.2232
0.69814	-2.0278	-9.063	-10.9605	7711.6431	7380.6361	7123.3598	8135.8516
3.162	-2.9238	0.51191	-1.8363	7352.0065	7582.1476	7038.649	7767.0063
-5.7397	0.94185	-11.7277	-2.2608	8076.072	7434.3015	7475.502	7983.5962
]; % The first 4 elements are the D-FR for 39,31,15 and 7, the last 4 elements are the peak amplitudes for these 4 sites


Result_Data = [3.30895	0.589975	-11.637394	-16.22226	59728.4409	58224.7973	58175.9921	61997.3317
];


% Calculate the resultant vector
Result_D_FR = repmat(Result_Data(1:4)',1,2);
Result_Arts = repmat(Result_Data(5:end)',1,2);
Result_D_FR_Vec = Site_Vecs.*Result_D_FR;
Result_Arts_Vec = Site_Vecs.*Result_Arts;
Result_D_FR_Vec = sum(Result_D_FR_Vec,1);
Result_Arts_Vec = sum(Result_Arts_Vec,1);
Resultant_Angle_Degrees  = acos(dot(Result_D_FR_Vec/norm(Result_D_FR_Vec),Result_Arts_Vec/norm(Result_Arts_Vec)))*180/pi


Result_D_FR_Vec = [Result_D_FR_Vec/norm(Result_D_FR_Vec) 0];
Result_Arts_Vec = [Result_Arts_Vec/norm(Result_Arts_Vec) 0];

% find the minimum of the mean amplitudes of 4 contacts
minAmp = min(mean(Data(:,5:end),2));

FR_Vec_starts = zeros(size(Data,1),3); 
FR_Vec_ends = zeros(size(Data,1),3); 
Art_Vec_starts = zeros(size(Data,1),3); 
Art_Vec_ends = zeros(size(Data,1),3); 

for i = 1:size(Data, 1)
    
    D_FR = repmat(Data(i,1:4)',1,2); % D_FR for all 4 sites
    Arts = repmat(Data(i,5:end)',1,2); % the peak artifact for all 4 sites
    FR_Vec = Site_Vecs.*(D_FR);
    Art_Vec =Site_Vecs.*(Arts);
    
    FR_Vec = sum(FR_Vec, 1);
    Art_Vec = sum(Art_Vec, 1);
    
    FR_Vec = [FR_Vec 0]; % not normalized
    Art_Vec = [Art_Vec 0];
    FR_Vec_ends(i,:) = FR_Vec;
    Art_Vec_ends(i,:) = Art_Vec;
    
    FR_Vec_starts(i,:) = [0,0,mean(Data(i,5:end))/minAmp];
    Art_Vec_starts(i,:) = [0,0,mean(Data(i,5:end))/minAmp]; 
    
end

% Find the smallest of the resultant vectors in each category

min_FR_Vec = min(sqrt(sum(abs(FR_Vec_ends).^2,2)));
min_Art_Vec = min(sqrt(sum(abs(Art_Vec_ends).^2,2)));
% normalize all the vectors to the smallest resultant vector in each
% category
FR_Vec_ends = FR_Vec_ends/min_FR_Vec;
Art_Vec_ends = Art_Vec_ends/min_Art_Vec;
% Plot everything out 
for j = 1: size(Data,1)
    
    quiver3(FR_Vec_starts(j,1), FR_Vec_starts(i,2), FR_Vec_starts(j,3), FR_Vec_ends(j,1), FR_Vec_ends(j,2), FR_Vec_ends(j,3),'color','b'); hold on
    quiver3(Art_Vec_starts(j,1), Art_Vec_starts(j,2), Art_Vec_starts(j,3), Art_Vec_ends(j,1), Art_Vec_ends(j,2), Art_Vec_ends(j,3),'color','r'); hold on; 
    l = legend('dFR_Resultant_Direction', 'Artifact_Result_Direction');
    set(l, 'Interpreter', 'none')
    
end


Result_Start = [ 0 0 0.99; 0,0, 0.99];
Result_End = [Result_D_FR_Vec; Result_Arts_Vec];
% Randomized_End = [Randomized_Vec];
quiver3(Result_Start(1,1), Result_Start(1,2), Result_Start(1,3),Result_End(1,1), Result_End(1,2), Result_End(1,3),'color','b', 'LineWidth', 5); hold on
quiver3(Result_Start(2,1), Result_Start(2,2), Result_Start(2,3),Result_End(2,1), Result_End(2,2), Result_End(2,3),'color','r', 'LineWidth', 5);
% quiver3(Result_Start(2,1), Result_Start(2,2), Result_Start(2,3),Randomized_End(1,1), Randomized_End(1,2), Randomized_End(1,3),'color','y', 'LineWidth', 5);
axis fill
title(['4.25_150uA_dFR_Artifact_Vector', ' n = ', num2str(4*size(Data,1)), ' Resultant_Angle = ' num2str(Resultant_Angle_Degrees)],'interpreter', 'none')
view(-40.5,12)

%% Plot angle vs mean amplitude (of 4 contacts)

[curve, goodness, output] = fit(Data(:,1),Data(:,2),'poly1');
plot(curve); hold on;
plot(Data(:,1), Data(:,2), 'o');

goodness.rsquare
title(['Angle vs Mean Artifact Amplitude', ' r^2 = ', ' ' num2str(goodness.rsquare) ])
xlabel(['Peak Artifact Amplitude (uV)'])
ylabel(['Angle between D_FR and Mean Amplitude'], 'interpreter', 'none')




%% Plot the resultant vectors using D_FR
clear all, clc, close all
Site39 = [1 0];
Site31 = [0,1];
Site15 = [-1,0];
Site7 = [0, -1];

Site39_FR = Site39*3.30895;
Site31_FR = Site31*-7.622211;
Site15_FR = Site15*35.467356;
Site7_FR = Site7*102.627843;

Site39_Art = Site39*422774.5104;
Site31_Art = Site31*417611.7509;
Site15_Art = Site15*428375.1714;
Site7_Art = Site7*439756.4849;
FR_Vec = Site39_FR+Site31_FR+Site15_FR+Site7_FR
Art_Vec = Site39_Art + Site31_Art+ Site15_Art+ Site7_Art


% a = [0 2];
% b = [2 0];
% c = a+b;
% c = c/norm(c)
% starts = zeros(5,5);
% ends = [a;b;c];



starts = zeros(2,2);
ends = [FR_Vec/norm(FR_Vec); Art_Vec/norm(Art_Vec)];

quiver(starts(1,1), starts(1,2),  ends(1,1), ends(1,2), 'color', 'b'); hold on
quiver(starts(2,1), starts(2,2),  ends(2,1), ends(2,2), 'color', 'r');
l = legend('dFR_Resultant_Direction', 'Artifact_Result_Direction');
set(l, 'Interpreter', 'none')
axis equal

angle_Degrees  = acos(dot(FR_Vec/norm(FR_Vec),Art_Vec/norm(Art_Vec)))*180/pi

title(['4.25.14&4.29.14_All_Data_dFR_Artifact_Vector', ' n = ', num2str(180), ' angle = ' num2str(angle_Degrees)],'interpreter', 'none')


%% Calculat the mean percent increase in peak artifact amplitude
clear all, clc, close all
artifacts_peaks = [6864.8333
    6591.761
    6696.5457
    7197.8087
    
    ];
artifacts_peaks_sorted = sort(artifacts_peaks,'ascend');
percent_increase = ((artifacts_peaks_sorted(2:end)./artifacts_peaks_sorted(1:end-1))-1)*100;

mean_percent_increase = mean(percent_increase);
max_percent_increase = (artifacts_peaks_sorted(end)./artifacts_peaks_sorted(1)-1)*100;


disp(num2str(percent_increase))
disp(num2str(mean_percent_increase))
disp(num2str(max_percent_increase))


%% For every row of amplitudes, plot them out
% The x-axis is the smallest amplitude of each row, the y-axis plots out
% all 4 amplitudes in that row
Min_Amps = zeros(1,length(Amplitudes)/4);
Max_Amps = zeros(1,length(Amplitudes)/4);
counter = 1;

for i = 1:4:length(Amplitudes)
    Sorted = sort(Amplitudes(i:i+3),'ascend');
    plot(Sorted(1), (Sorted(4)./Sorted(1)-1)*100,'ro'); hold on;
    %     plot(Sorted(1), Sorted(3)./Sorted(1)-1,'ko'); hold on;
    %     plot(Sorted(1), Sorted(2)./Sorted(1)-1,'co'); hold on;
    %     plot(Sorted(1), Sorted(1)./Sorted(1)-1,'mo'); hold on;
    Min_Amps(counter) = Sorted(1);
    Max_Amps(counter) = (Sorted(4)./Sorted(1)-1)*100;
    counter = counter + 1;
end

[curve, goodness, output] = fit(Min_Amps',Max_Amps','poly1');
plot(curve)
goodness.rsquare
title(['Max percent change vs min row amplitude', ' r^2 = ', ' ' num2str(goodness.rsquare) ])
xlabel(['Min Artifact Amplitude of Row (uV)'])
ylabel(['Max Percent Increase'])


%% Plot the mean firing rate vs peak artifact amplitude
clc, close all
[Peak_Art_Amp_sorted sorted_indices] = sort(Peak_Art_Amp,'ascend');
mean_FR_sorted = mean_FR(sorted_indices);
duplicate_ind = [];
for i = 1:length(Peak_Art_Amp_sorted)-1
    if Peak_Art_Amp_sorted(i) == Peak_Art_Amp_sorted(i+1);
        duplicate_ind = [duplicate_ind, i];
    end
end
Peak_Art_Amp_sorted(duplicate_ind) = [];
mean_FR_sorted(duplicate_ind) = [];

[curve, goodness, output] = fit(Peak_Art_Amp_sorted,mean_FR_sorted,'poly1');
plot(Peak_Art_Amp_sorted,mean_FR_sorted,'o'); hold on;
plot(curve)
goodness.rsquare
title(['CV vs Peak Artifact Amplitude', ' r^2 = ', ' ' num2str(goodness.rsquare) ])

% calculate the average change in firing rate when increasing by 1% of max stim
% amplitude

Mean_Change = mean(diff(mean_FR_sorted)./((Peak_Art_Amp_sorted(2:end)./Peak_Art_Amp_sorted(1:end-1)-1)*100));
std_Mean_Change = std(diff(mean_FR_sorted)./((Peak_Art_Amp_sorted(2:end)./Peak_Art_Amp_sorted(1:end-1)-1)*100));


%% %% Plot the change in mean firing rate vs peak artifact amplitude
% first import the excel data
clc, close all


plot(PAA, D_FR,'o'); hold on;


[curve, goodness, output] = fit(PAA,D_FR,'poly1');
plot(PAA,D_FR,'o'); hold on;
plot(curve)
goodness.rsquare
title(['%Diff in FR vs Peak Artifact Amplitude', ' r^2 = ', ' ' num2str(goodness.rsquare) ])
xlabel(['Peak Artifact Amplitude (uV)'])
ylabel(['%Difference in Firing Rate (Hz)'])

%% Plot type of change vs amplitude.
% The y axis is the type of  change: pr, p, r, n, r-,p-,p-r-, from top to
% bottom
for i = 1: length(Type_350)
    
    if strfind(Type_350{i},'/') % means there's more than one cell in there
        index = strfind(Type_350{i},'/');
        first = Type_350{i}(1:index-1);
        second = Type_350{i}(index+1:end);
        switch first
            case 'pr'
                plot(PAA_350(i),7,'ro'); hold on;
            case 'p'
                plot(PAA_350(i),6,'ro'); hold on;
            case 'r'
                plot(PAA_350(i),5,'ro'); hold on;
            case 'n'
                plot(PAA_350(i),4,'ro'); hold on;
            case 'r-'
                plot(PAA_350(i),3,'ro'); hold on;
            case 'p-'
                plot(PAA_350(i),2,'ro'); hold on;
            case 'p-r-'
                plot(PAA_350(i),1,'ro'); hold on;
        end
        switch second
            case 'pr'
                plot(PAA_350(i),7,'ro'); hold on;
            case 'p'
                plot(PAA_350(i),6,'ro'); hold on;
            case 'r'
                plot(PAA_350(i),5,'ro'); hold on;
            case 'n'
                plot(PAA_350(i),4,'ro'); hold on;
            case 'r-'
                plot(PAA_350(i),3,'ro'); hold on;
            case 'p-'
                plot(PAA_350(i),2,'ro'); hold on;
            case 'p-r-'
                plot(PAA_350(i),1,'ro'); hold on;
        end
    else
        content = Type_350{i};
        switch content
            case 'pr'
                plot(PAA_350(i),7,'ro'); hold on;
            case 'p'
                plot(PAA_350(i),6,'ro'); hold on;
            case 'r'
                plot(PAA_350(i),5,'ro'); hold on;
            case 'n'
                plot(PAA_350(i),4,'ro'); hold on;
            case 'r-'
                plot(PAA_350(i),3,'ro'); hold on;
            case 'p-'
                plot(PAA_350(i),2,'ro'); hold on;
            case 'p-r-'
                plot(PAA_350(i),1,'ro'); hold on;
        end
        
    end
    
end

for i = 1: length(Type_250)
    
    if strfind(Type_250{i},'/') % means there's more than one cell in there
        index = strfind(Type_250{i},'/');
        first = Type_250{i}(1:index-1);
        second = Type_250{i}(index+1:end);
        switch first
            case 'pr'
                plot(PAA_250(i),7,'co'); hold on;
            case 'p'
                plot(PAA_250(i),6,'co'); hold on;
            case 'r'
                plot(PAA_250(i),5,'co'); hold on;
            case 'n'
                plot(PAA_250(i),4,'co'); hold on;
            case 'r-'
                plot(PAA_250(i),3,'co'); hold on;
            case 'p-'
                plot(PAA_250(i),2,'co'); hold on;
            case 'p-r-'
                plot(PAA_250(i),1,'co'); hold on;
        end
        switch second
            case 'pr'
                plot(PAA_250(i),7,'co'); hold on;
            case 'p'
                plot(PAA_250(i),6,'co'); hold on;
            case 'r'
                plot(PAA_250(i),5,'co'); hold on;
            case 'n'
                plot(PAA_250(i),4,'co'); hold on;
            case 'r-'
                plot(PAA_250(i),3,'co'); hold on;
            case 'p-'
                plot(PAA_250(i),2,'co'); hold on;
            case 'p-r-'
                plot(PAA_250(i),1,'co'); hold on;
        end
    else
        content = Type_250{i};
        switch content
            case 'pr'
                plot(PAA_250(i),7,'co'); hold on;
            case 'p'
                plot(PAA_250(i),6,'co'); hold on;
            case 'r'
                plot(PAA_250(i),5,'co'); hold on;
            case 'n'
                plot(PAA_250(i),4,'co'); hold on;
            case 'r-'
                plot(PAA_250(i),3,'co'); hold on;
            case 'p-'
                plot(PAA_250(i),2,'co'); hold on;
            case 'p-r-'
                plot(PAA_250(i),1,'co'); hold on;
        end
        
    end
    
end

for i = 1: length(Type_150)
    
    if strfind(Type_150{i},'/') % means there's more than one cell in there
        index = strfind(Type_150{i},'/');
        first = Type_150{i}(1:index-1);
        second = Type_150{i}(index+1:end);
        switch first
            case 'pr'
                plot(PAA_150(i),7,'bo'); hold on;
            case 'p'
                plot(PAA_150(i),6,'bo'); hold on;
            case 'r'
                plot(PAA_150(i),5,'bo'); hold on;
            case 'n'
                plot(PAA_150(i),4,'bo'); hold on;
            case 'r-'
                plot(PAA_150(i),3,'bo'); hold on;
            case 'p-'
                plot(PAA_150(i),2,'bo'); hold on;
            case 'p-r-'
                plot(PAA_150(i),1,'bo'); hold on;
        end
        switch second
            case 'pr'
                plot(PAA_150(i),7,'bo'); hold on;
            case 'p'
                plot(PAA_150(i),6,'bo'); hold on;
            case 'r'
                plot(PAA_150(i),5,'bo'); hold on;
            case 'n'
                plot(PAA_150(i),4,'bo'); hold on;
            case 'r-'
                plot(PAA_150(i),3,'bo'); hold on;
            case 'p-'
                plot(PAA_150(i),2,'bo'); hold on;
            case 'p-r-'
                plot(PAA_150(i),1,'bo'); hold on;
        end
    else
        content = Type_150{i};
        switch content
            case 'pr'
                plot(PAA_150(i),7,'bo'); hold on;
            case 'p'
                plot(PAA_150(i),6,'bo'); hold on;
            case 'r'
                plot(PAA_150(i),5,'bo'); hold on;
            case 'n'
                plot(PAA_150(i),4,'bo'); hold on;
            case 'r-'
                plot(PAA_150(i),3,'bo'); hold on;
            case 'p-'
                plot(PAA_150(i),2,'bo'); hold on;
            case 'p-r-'
                plot(PAA_150(i),1,'bo'); hold on;
        end
        
    end
    
end

ylim([0 8])
legend_handle = legend(['350uA'], ['250uA'], ['150uA']);
legend_markers = findobj(get(legend_handle, 'Children'), ...
    'Type', 'line', '-and', '-not', 'Marker', 'none');

set(legend_markers(3), 'Color', 'r');
set(legend_markers(2), 'Color', 'c');
set(legend_markers(1), 'Color', 'b');

title('Site39_Type_vs_Amplitude','interpreter','none')
xlabel('Amplitude (uV)')
ylabel('Change Type')


%% 3D bar plot showing for each category, count vs amplitude
% [Amplitudes_Sorted,sorted_indices] = sort(Amplitudes, 'ascend');
% Type_Sorted = Type(sorted_indices)
change_types = {'pr','p','r','n','r-','p-','p-r-'};

count = zeros(length(Type_Sorted),7);


for i = 1: length(Type_Sorted)
    
    if strfind(Type_Sorted{i},'/') % means there's more than one cell in there
        index = strfind(Type_Sorted{i},'/');
        first = Type_Sorted{i}(1:index-1);
        second = Type_Sorted{i}(index+1:end);
        switch first
            case 'pr'
                count(i,1) = count(i,1)+1;
            case 'p'
                count(i,2) = count(i,2)+1;
            case 'r'
                count(i,3) = count(i,3)+1;
            case 'n'
                count(i,4) = count(i,4)+1;
            case 'r-'
                count(i,5) = count(i,5)+1;
            case 'p-'
                count(i,6) = count(i,6)+1;
            case 'p-r-'
                count(i,7) = count(i,7)+1;
        end
        switch second
            case 'pr'
                count(i,1) = count(i,1)+1;
            case 'p'
                count(i,2) = count(i,2)+1;
            case 'r'
                count(i,3) = count(i,3)+1;
            case 'n'
                count(i,4) = count(i,4)+1;
            case 'r-'
                count(i,5) = count(i,5)+1;
            case 'p-'
                count(i,6) = count(i,6)+1;
            case 'p-r-'
                count(i,7) = count(i,7)+1;
                
        end
    else
        content = Type_Sorted{i};
        switch content
            case 'pr'
                count(i,1) = count(i,1)+1;
            case 'p'
                count(i,2) = count(i,2)+1;
            case 'r'
                count(i,3) = count(i,3)+1;
            case 'n'
                count(i,4) = count(i,4)+1;
            case 'r-'
                count(i,5) = count(i,5)+1;
            case 'p-'
                count(i,6) = count(i,6)+1;
            case 'p-r-'
                count(i,7) = count(i,7)+1;
                
        end
        
    end
    
end

% Turn each column of the matrix into a running tally

counts_sum = cumsum(counts, 1);

Amplitudes_Sorted = num2cell(Amplitudes_Sorted);

% Create the 3D bar chart
figure;
width = 0.5;
bar3(counts_sum,width);

axis([10 70 0 length(Amplitudes_Sorted)+1 0 max(max(counts_sum))] )

% Add title and axis labels
title('All Data Change vs Amplitude');
xlabel('Type of Change');
ylabel('Amplitudes (uV)');
zlabel('Counts');

set(gca,'XLim',[1 7])
set(gca,'XTick',[1:7])

set(gca,'YTick',[1:5:length(Amplitudes_Sorted)])
% Change the x and y axis tick labels
set(gca, 'XTickLabel', change_types);
set(gca, 'YTickLabel', Amplitudes_Sorted(1:5:end));

%% 3D line plot showing for each category, count vs amplitude
% [Amplitudes_Sorted,sorted_indices] = sort(Amplitudes, 'ascend');
% Type_Sorted = Type(sorted_indices)
change_types = {'pr','p','r','n','r-','p-','p-r-'};

% Turn each column of the matrix into a running tally

counts_sum = cumsum(counts, 1);
Amplitudes_Sorted = repmat(Amplitudes_Sorted,1,size(counts,2));
change_Types = repmat([1:7],size(counts,1),1);
% Create the 3D bar chart
figure;
width = 0.5;
plot3(Amplitudes_Sorted, change_Types,counts_sum, '-');

axis([1 7 0 length(Amplitudes_Sorted)+1 0 max(max(counts_sum))] )

% Add title and axis labels
title('Site39 Cumulative Count vs Amplitude');
ylabel('Type of Change');
xlabel('Amplitudes (uV)');
zlabel('Counts');
set(gca,'XTick',[1:5:size(Amplitudes_Sorted,1)])
% Change the x and y axis tick labels
set(gca, 'YTickLabel', change_types);
set(gca, 'XTickLabel', Amplitudes_Sorted(1:5:end));


%% 2D line plot showing for each category, count vs amplitude
% [Amplitudes_Sorted,sorted_indices] = sort(Amplitudes, 'ascend');
% Type_Sorted = Type(sorted_indices)
change_types = {'pr','p','r','n','r-','p-','p-r-'};

% Turn each column of the matrix into a running tally

counts_sum = cumsum(counts, 1);
Amplitudes_Sorted = repmat(Amplitudes_Sorted,1,size(counts,2));
figure;

plot(Amplitudes_Sorted, counts_sum, 'o-'); hold on


% Add title and axis labels
title('Cumulative Count vs Amplitude');
xlabel('Amplitudes (uV)');
ylabel('Cumulative Count');
legend({'pr','p','r','n','r-','p-','p-r-'});


%% Plot the Type count vs amplitude, binning the amplitudes.
change_types = {'pr','p','r','n','r-','p-','p-r-'};
binSize = 60;
for i = 1:size(counts,2)
    subplot(7,1,i)
    binranges = linspace(min(Amplitudes),max(Amplitudes),binSize);
    activity1_Amplitudes = Amplitudes(find(counts(:,i)==1)); % These are the amplitudes that has 1 cell modulated for the current category
    activity2_Amplitudes = Amplitudes(find(counts(:,i)==2)); % These are the amplitudes that has 2 cell modulated for the current category
    activity2_Amplitudes = repmat(activity2_Amplitudes, 2,1); % have the duplicate these amplitudes since there are 2 cells for each
    activity_Amplitudes = horzcat(activity1_Amplitudes', activity2_Amplitudes');
    bincounts = histc(activity_Amplitudes,binranges);
    sum(bincounts)
    bar(binranges,bincounts,'histc')
    set(gca, 'YLim', [0 15])
    ylabel(change_types{i})
    if i ==1
        title(['Change Type Count vs. Peak Stimulus Amplitude bin size = ', num2str(binSize)]);
    end
    
end
xlabel('Peak Stimulus Artifact Amplitude (uV)')


%% Calculate percent change


for i = 1:length(mean)
    mean_content = mean{i};
    diff_content = diff{i};
    if strfind(mean_content,'/') % there are two cells
        mean_index = strfind(mean_content,'/');
        diff_index = strfind(diff_content,'/');
        C1 = (str2num(mean_content(1:mean_index-1))/(str2num(mean_content(1:mean_index-1)) - str2num(diff_content(1:diff_index-1)))-1)*100;
        C2 = (str2num(mean_content(mean_index+1:end))/(str2num(mean_content(mean_index+1:end)) - str2num(diff_content(diff_index+1:end))))*100;
        disp([num2str(C1) '/' num2str(C2)])
    else
        C = (str2num(mean_content)/(str2num(mean_content) - str2num(diff_content)))*100;
        disp([num2str(C)])
    end
end

C = (mean./(mean-D)-1)*100


