% This file gathers the depth vs. amplitude data
clear all, clc 
% List the recording files to be analyzed, in order
fileFormat = ['15_07_09k_00'];
recNo = [2, 3, 4, 5, 6, 7, 10,11,12,13,14,15,16,17,18,  19,20,  22,23,  24,25,  26];
depth = [15,17,19,21,23,25,27,29,31,32,33,34,35,36,36.5,37,37.5,38,38.5,39,39.5,40];
dataTable = zeros(size(numel(recNo, 5)));

for i = 1:numel(recNo)
    if recNo(i)<10
        load([fileFormat,'0', num2str(recNo(i)), '.mat']);
    else
        load([fileFormat, num2str(recNo(i)), '.mat']);
    end
    [artifactwave_maxavg, artifactwave_minavg, peak2peak] = artifactpeak_calculation(CSPK_001, CSPK_001_KHz);
    dataTable(i,1) = recNo(i);
    dataTable(i,2) = depth(i);
    dataTable(i,3) = artifactwave_maxavg;
    dataTable(i,4) = artifactwave_minavg;
    dataTable(i,5) = peak2peak;  
    clearvars -except fileFormat recNo depth dataTable i
end
plot(dataTable(:,2), dataTable(:,5),'r.-', 'markersize', 10); title('peak2peak_vs_distance'), xlabel('distance'),ylabel('amplitude');
plot(dataTable(:,2), dataTable(:,3),'b.-', 'markersize', 10); title('avgMaxAmp_vs_distance'), xlabel('distance'),ylabel('amplitude');
plot(dataTable(:,2), dataTable(:,4),'g.-', 'markersize', 10); title('avgMinAmp_vs_distance'), xlabel('distance'),ylabel('amplitude');
