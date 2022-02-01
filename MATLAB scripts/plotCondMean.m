%MATLAB script to plot conditional means

%clear all vars and plots
clear all;
close all; 
clc;

%load data
pltfile = "plt_02490"
dat=dlmread(pltfile+'_CM.dat')
data = dlmread(pltfile+'_scatter.dat',' ',1,0)

names = ["Temperature(K)","Y_{CH3OCH2O2}", "Y_{OCH2OCHO}","Y_{CH2O}","heatRelease","density","\chi"]
colors = [0.9290 0.6940 0.1250;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0 0 1;...
    1 1 0;...
    0 1 0;...
    1 0 1
    ]
%initialize vars
nVars=7;
oVar=0; % offset (counting from 0)

iZ=1;

%Conditional means
clf()
r=tiledlayout(2,4);
title(r,'Mixture fraction (\xi)-conditioned averages of variables');
xlabel(r,'Mixture Fraction (\xi)');
for i=2:nVars+1
    iSum=i;
    iSSq=iSum+nVars;
    iAvg=iSSq+nVars;
    iStd=iAvg+nVars;
    iN=iStd+nVars;
    ip=iN+1;

    c=[colors(i-1,1) colors(i-1,2) colors(i-1,3)];

    x=dat(:,iZ);
    a=dat(:,iAvg+oVar);
    s=dat(:,iStd+oVar);
    %plot stddev band
    nexttile
    fillStdDev(x,a-s,a+s,c,0);
    hold on;
    plot(x,a,'k-','LineWidth',2);
    xlim([0 inf])
    ylim([0 inf])
    title(names(i-1))
    xline(0.1024,'r-',{'\xi_{st}=0.1024'},'LabelOrientation','horizontal','LabelVerticalAlignment','top','LineWidth',1)
end
saveas(r,'conditionalMean_t1.png')
