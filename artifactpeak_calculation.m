function [artifactwave_maxavg, artifactwave_minavg, peak2peak] = artifactpeak_calculation(CSPK_001, CSPK_001_KHz)
%% 
% Define an arbitrary threshold
thr = 1000;

% Find threshold crossing timestamps
thr_cross = find(CSPK_001>thr);

% Calculate a threshold crossing index for each stimulus pulse
thr_begin = thr_cross(1);
ctr_begin = 2;
for i = 2:length(thr_cross)-1,
    if thr_cross(i-1) < thr_cross(i)-1,
        thr_begin(ctr_begin) = thr_cross(i);
        ctr_begin = ctr_begin + 1;
    end
end

% Define waveforms for each threshold crossing
artifactwin = [CSPK_001_KHz*2, CSPK_001_KHz*3];
for i=1:length(thr_begin),
    artifactwave(i,:) = CSPK_001(thr_begin(i)-artifactwin(1):thr_begin(i)+artifactwin(2));
end

% Calculate the average minimum and maximum artifact points
artifactwave_minavg = mean(min(artifactwave,[],2));
artifactwave_maxavg = mean(max(artifactwave,[],2));
peak2peak = artifactwave_maxavg - artifactwave_minavg;