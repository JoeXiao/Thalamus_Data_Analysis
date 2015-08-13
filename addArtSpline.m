clear, close all, clc
[fn,pn]=uigetfile('.mat','Pick File','C:\Users\Joe Xiao\Desktop\Lab\Research\Projects\Thalamus_3D_Map_Project\Data_Analysis');
load([pn,fn]);

for jj = 1: size(F.stim,2)
    meanArt = F.stim(jj).meanArt;
    
    x = 0:length(meanArt)-1;
    y = meanArt;
    xx = 0:0.01:length(meanArt);
    yy = spline(x,y,xx);
    F.stim(jj).splineMeanArt = yy;
    F.stim(jj).splineMeanArtMax = max(yy);
    disp(['stim ' num2str(jj) ' max spline artifact amplitude is ' num2str(max(yy))  ])
end

save([pn,fn],'F'); % saves matlab output