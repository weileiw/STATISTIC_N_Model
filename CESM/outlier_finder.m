clc 
clear all
close all

load ~/DATA/CESM_output.mat
load NF_wo_outliers.mat
load ~/DATA/tempobs_90x180x24.mat
load ~/DATA/radiation_90x180.mat 

%load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat
spy = 365*24*60^2;
%NF2d = nansum(NF,3)*10*spy; % the upper 16 layers have the same thickness (10m);
SI = ocean_solar;

O2 = O2(:,:,1:33); % select the upper ~500m.
O2min = nanmin(O2,[],3);

SST = nanmean(tempobs(:,:,1:2),3);
DFe = nanmean(Fe(:,:,1:10),3);
DIP = nanmean(PO4(:,:,1:10),3);
DIN = nanmean(NO3(:,:,1:10),3);
DOP = nanmean(DOP(:,:,1:10),3);
DIP = DIP+DOP;

iposi = all([NF2d(:)>0,DIP(:)>0,DIN(:)>0,DFe(:)>0,SST(:)>0,SI(:)>0],2); 

% use only positive numbers.

NF  = log10(NF2d(iposi));
DIP = log10(DIP(iposi));
DIN = log10(DIN(iposi));
DFe = log10(DFe(iposi));
SST   = SST(iposi);
SI    = SI(iposi);
O2min = O2min(iposi);

% get z-score
NF = (NF-mean(NF))/std(NF);
DFe   = (DFe-mean(DFe))/std(DFe);
DIP  = (DIP-mean(DIP))/std(DIP);
DIN  = (DIN-mean(DIN))/std(DIN);
SST  = (SST-mean(SST))/std(SST);
SI  = (SI-mean(SI))/std(SI);
O2min  = (O2min-mean(O2min))/std(O2min);

fprintf('Global ocean regression\n')
tbl = [DFe,DIP,DIN,SST,SI,O2min];
mdl = stepwiselm(tbl,NF,'quadratic','Criterion','Rsquared','PRemove',0.05)

% get rid of outliers based on Cook's distance.
ibad = find((mdl.Diagnostics.CooksDistance)>5*mean(mdl.Diagnostics.CooksDistance));
NF(ibad) = nan;
NF2d_tmp = NF2d;
NF2d_tmp(iposi) = NF;
inan = find(isnan(NF2d_tmp(:)));
NF2d(inan) = -999;
%fname = sprintf('NF_wo_outliers');
%save(fname,'NF2d')

