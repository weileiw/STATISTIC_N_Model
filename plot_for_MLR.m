clc
clear all
close all

load NF2d_OPT_N2P_tanh.mat
load DFe_from_Ben.mat
load ~/DATA/no3obs_90x180x24.mat
load ~/DATA/po4obs_90x180x24.mat
load ~/DATA/tempobs_90x180x24.mat
load ~/DATA/o2obs_90x180x24.mat
load ~/DATA/radiation_90x180.mat

load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat

o2obs = o2obs*44.661;      % convert unit form [ml/l] to [umol/l].
o2obs = o2obs*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l].
o2obs = o2obs(:,:,1:8); % select the upper 447m.
O2min = nanmin(o2obs,[],3);

DFe = nanmean(DFe_3d(:,:,1:2),3);
DIP = nanmean(po4obs(:,:,1:2),3);
DIN = nanmean(no3obs(:,:,1:2),3);
Pstar = DIP-DIN/16;
SST = nanmean(tempobs(:,:,1:2),3);
SI  = ocean_solar;

iposi = all([NF2d(:)>0,DIP(:)>0,DIN(:)>0,DFe(:)>0],2); 
% use only positive numbers.
NF = log10(NF2d(iposi));
DIP = log10(DIP(iposi));
DIN = log10(DIN(iposi));
DFe = log10(DFe(iposi));
Pstar = Pstar(iposi);
SST   = SST(iposi);
SI    = SI(iposi);
O2min = O2min(iposi);

% get z-score
NFz = (NF-mean(NF))/std(NF);
DFe = (DFe-mean(DFe))/std(DFe);
DIP = (DIP-mean(DIP))/std(DIP);
DIN = (DIN-mean(DIN))/std(DIN);
Pstar = (Pstar-mean(Pstar))/std(Pstar);
SST = (SST-mean(SST))/std(SST);
SI = (SI-mean(SI))/std(SI);
O2min = (O2min-mean(O2min))/std(O2min);

T = [0 0 0 0 0 0 0; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0;...
      0 0 0 0 0 1 0; 2 0 0 0 0 0 0; 0 2 0 0 0 0 0; 0 0 2 0 0 0 0; 0 0 0 2 0 0 0;...
      0 0 0 0 2 0 0; 0 0 0 0 0 2 0];
tbl = [DFe,DIP,DIN,SST,SI,O2min];
mdl = stepwiselm(tbl,NFz,'quadratic','Criterion','Rsquared','PEnter',0.1)
%mdl = stepwiselm(tbl,NFz,'quadratic','Criterion','Rsquared','PRemove',0.02)
Intercept = -0.06989;    
x3 =       -0.75705;  
x4 =        0.22136; 
x3x4 =      0.71076;
x4square = 0.49398;

NF_fit = nanmean(M3d(:,:,1:2),3)+nan;
tmp = x3*DIN+x4*SST+x3x4*DIN.*SST+x4square*SST.^2+Intercept;
NF_fit(iposi) = tmp*std(NF)+mean(NF);
figure()
plot(NF,tmp,'+')
fname = sprintf('Glob_pridiction');
save(fname,'NF_fit')
figure()
contourf(NF_fit);colorbar
