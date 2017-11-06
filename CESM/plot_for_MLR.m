clc
clear all
close all
load ~/DATA/CESM_output.mat PO4 NO3 DOP Fe NH4 NF dz O2 
load ~/DATA/tempobs_90x180x24.mat
load ~/DATA/radiation_90x180.mat
load NF_wo_outliers.mat 
load Fe_lim.mat 
load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat
spy = 365*24*60^2;
SI  = ocean_solar;

O2 = O2(:,:,1:33); % select the upper ~500m.
O2min = nanmin(O2,[],3);

Fe_lim = nanmean(Fe_lim(:,:,1:10),3);
SST = nanmean(tempobs(:,:,1:2),3);
DFe = nanmean(Fe(:,:,1:10),3);
PO4 = nanmean(PO4(:,:,1:10),3);
NO3 = nanmean(NO3(:,:,1:10),3);
NH4 = nanmean(NH4(:,:,1:10),3);
DOP = nanmean(DOP(:,:,1:10),3);

DIP = PO4+DOP;
DIN = NO3+NH4;

Pstar = DIP-DIN/16;
iposi = all([NF2d(:)>0,DIP(:)>0,DIN(:)>0,DFe(:)>0,SST(:)>0,SI(:)>0],2); 
% use only positive numbers.

NF = log10(NF2d(iposi));
DIP = log10(DIP(iposi));
DIN = log10(DIN(iposi));
DFe = log10(DFe(iposi));
Fe_lim = Fe_lim(iposi);
Pstar = Pstar(iposi);
SST   = SST(iposi);
SI    = SI(iposi);
O2min = O2min(iposi);

% get z-score
NFz   = (NF-mean(NF))/std(NF);
%DFe   = (DFe-mean(DFe))/std(DFe);
Fe_lim   = (Fe_lim-mean(Fe_lim))/std(Fe_lim);
DIP   = (DIP-mean(DIP))/std(DIP);
DIN   = (DIN-mean(DIN))/std(DIN);
Pstar = (Pstar-mean(Pstar))/std(Pstar);
SST   = (SST-mean(SST))/std(SST);
SI    = (SI-mean(SI))/std(SI);
O2min = (O2min-mean(O2min))/std(O2min);
keyboard
T = [0 0 0 0 0 0 0; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0;...
     0 0 0 0 0 1 0; 2 0 0 0 0 0 0; 0 2 0 0 0 0 0; 0 0 2 0 0 0 0; 0 0 0 2 0 0 0;...
     0 0 0 0 2 0 0; 0 0 0 0 0 2 0];
tbl = [Fe_lim,DIP,DIN,SST,SI,O2min];
mdl = stepwiselm(tbl,NFz,'quadratic','Criterion','Rsquared','PRemove',0.05)
keyboard
Intercept = mdl.Coefficients.Estimate(1); 
x2        = mdl.Coefficients.Estimate(2); 
x5        = mdl.Coefficients.Estimate(3); 
x6        = mdl.Coefficients.Estimate(4); 
x2square  = mdl.Coefficients.Estimate(5); 
x5square  = mdl.Coefficients.Estimate(6); 


NF_fit = nanmean(M3d(:,:,1:2),3)+nan;
tmp = x2*DIP+x5*SI+x6*O2min+x2square*DIP.^2+x5square*SI.^2+Intercept;
NF_fit(iposi) = 10.^(tmp*std(NF)+mean(NF));
figure()
plot(NFz,tmp,'+')
fname = sprintf('Glob_pridiction');
save(fname,'NF_fit')
figure()
contourf(NF_fit);colorbar
