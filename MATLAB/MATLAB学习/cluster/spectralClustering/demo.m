% clc
% clear

%twoCircles
%load('twoCircles.mat');
% load('Iris.mat');
% dataSetOri=Iris;
% dataSet=dataSetOri/(max(max(abs(dataSetOri)))-min(min(abs(dataSetOri))));
% num_clusters=3;
% sigma=0.1;

%XOR
% load('XOR.mat');
% dataSetOri=XOR;
dataSet=dataSetOri/(max(max(abs(dataSetOri)))-min(min(abs(dataSetOri))));
num_clusters=4;
sigma=0.1;

Z=pdist(dataSet);
W=squareform(Z);

C = spectral(W,sigma, num_clusters);

plot(dataSetOri(C==1,1),dataSetOri(C==1,2),'ro', dataSetOri(C==2,1),dataSetOri(C==2,2),'bo', dataSetOri(C==3,1),dataSetOri(C==3,2),'go', dataSetOri(C==4,1),dataSetOri(C==4,2),'mo');


