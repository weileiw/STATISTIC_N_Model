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

%load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat

SI = ocean_solar;

o2obs = o2obs*44.661;      % convert unit form [ml/l] to [umol/l].
o2obs = o2obs*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l].
o2obs = o2obs(:,:,1:8); % select the upper 447m.
O2min = nanmin(o2obs,[],3);

SST = nanmean(tempobs(:,:,1:2),3);
DFe = nanmean(DFe_3d(:,:,1:2),3);
DIP = nanmean(po4obs(:,:,1:2),3);
DIN = nanmean(no3obs(:,:,1:2),3);
Pstar = DIP-DIN/16;

iposi = all([NF2d(:)>0,DIP(:)>0,DIN(:)>0,DFe(:)>0],2); 
% use only positive numbers.

NF  = log10(NF2d(iposi));
DIP = log10(DIP(iposi));
DIN = log10(DIN(iposi));
DFe = log10(DFe(iposi));
Pstar = Pstar(iposi);
SST   = SST(iposi);
SI    = SI(iposi);
O2min = O2min(iposi);

% get z-score
NF  = (NF-mean(NF))/std(NF);
DFe = (DFe-mean(DFe))/std(DFe);
DIP = (DIP-mean(DIP))/std(DIP);
DIN = (DIN-mean(DIN))/std(DIN);
Pstar = (Pstar-mean(Pstar))/std(Pstar);
SST = (SST-mean(SST))/std(SST);
SI = (SI-mean(SI))/std(SI);
O2min = (O2min-mean(O2min))/std(O2min);

fprintf('Global ocean regression\n')
tbl = [DFe,DIP,DIN,SST,SI,O2min];
%mdl = stepwiselm(tbl,NF,'quadratic')
mdl = stepwiselm(tbl,NF,'quadratic','Criterion','Rsquared','PRemove',0.02)

keyboard
ATL_srf = MSKS.ATL(:,:,1);
ATL_srf = ATL_srf(iposi);
PAC_srf = MSKS.PAC(:,:,1);
PAC_srf = PAC_srf(iposi);
IND_srf = MSKS.IND(:,:,1);
IND_srf = IND_srf(iposi);
MED_srf = MSKS.MED(:,:,1);
MED_srf = MED_srf(iposi);
ARC_srf = MSKS.ARC(:,:,1);
ARC_srf = ARC_srf(iposi);

% apply masks to fixation.
NF_ATL    = NF.*ATL_srf;
DFe_ATL   = DFe.*ATL_srf;
DIP_ATL   = DIP.*ATL_srf;
DIN_ATL   = DIN.*ATL_srf;
Pstar_ATL = Pstar.*ATL_srf;
SST_ATL   = SST.*ATL_srf;
SI_ATL    = SI.*ATL_srf;
O2min_ATL = O2min.*ATL_srf;

NF_PAC    = NF.*PAC_srf;
DFe_PAC   = DFe.*PAC_srf;
DIP_PAC   = DIP.*PAC_srf;
DIN_PAC   = DIN.*PAC_srf;
Pstar_PAC = Pstar.*PAC_srf;
SST_PAC   = SST.*PAC_srf;
SI_PAC    = SI.*PAC_srf;
O2min_PAC = O2min.*PAC_srf;
SI_PAC    = SI.*PAC_srf;

NF_IND    = NF.*IND_srf;
SI_IND    = SI.*IND_srf;
DFe_IND   = DFe.*IND_srf;
DIP_IND   = DIP.*IND_srf;
DIN_IND   = DIN.*IND_srf;
Pstar_IND = Pstar.*IND_srf;
SST_IND   = SST.*IND_srf;
SI_IND    = SI.*IND_srf;
O2min_IND = O2min.*IND_srf;

NF_ARC = NF.*ARC_srf;
SI_ARC = SI.*ARC_srf;
% ========================================================================
fprintf('\n\n')
fprintf('Atlantic ocean regression \n')
% select data only in the Atlantic ocean for regression. 
ido = find(NF_ATL~=0);
PAR_ATL = [DFe_ATL(ido),DIP_ATL(ido),DIN_ATL(ido),Pstar_ATL(ido),...
          SST_ATL(ido),SI_ATL(ido),O2min_ATL(ido)];
mdl_ATL = stepwiselm(PAR_ATL,NF_ATL(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)
keyboard
%==========================================================================
fprintf('\n\n')
fprintf('Pacific ocean regression \n')
% select data only in the Pacific ocean for regression. 
ido   = find(NF_PAC~=0);

PAR_PAC = [DFe_PAC(ido),DIP_PAC(ido),DIN_PAC(ido),Pstar_PAC(ido),...
            SST_PAC(ido),SI_PAC(ido),O2min_PAC(ido)];
mdl_PAC = stepwiselm(PAR_PAC,NF_PAC(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)
% =========================================================================
fprintf('\n\n')
fprintf('Indian ocean regression \n')
% select data only in the Indian ocean for regression. 
ido = find(NF_IND~=0);
PAR_IND = [DFe_IND(ido),DIP_IND(ido),DIN_IND(ido),Pstar_IND(ido),...
           SST_IND(ido),SI_IND(ido),O2min_IND(ido)];
mdl_IND = stepwiselm(PAR_IND,NF_IND(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)

%===========================================================================
% ============================================================================
% create masks for north and south ocean sections.
NATL = MSKS.ATL;
NATL(1:45,:,:) = 0;
SATL = MSKS.ATL;
SATL(46:end,:,:) = 0;

NPAC = MSKS.PAC;
NPAC(1:45,:,:) = 0;
SPAC = MSKS.PAC;
SPAC(46:end,:,:) = 0;

NATL_srf = NATL(:,:,1);
NATL_srf = NATL_srf(iposi);
SATL_srf = SATL(:,:,1);
SATL_srf = SATL_srf(iposi);
NPAC_srf = NPAC(:,:,1);
NPAC_srf = NPAC_srf(iposi);
SPAC_srf = SPAC(:,:,1);
SPAC_srf = SPAC_srf(iposi);

% apply masks to fixation.
NF_NATL    = NF.*NATL_srf;
SI_NATL    = SI.*NATL_srf;
DFe_NATL   = DFe.*NATL_srf;
DIP_NATL   = DIP.*NATL_srf;
DIN_NATL   = DIN.*NATL_srf;
Pstar_NATL = Pstar.*NATL_srf;
SST_NATL   = SST.*NATL_srf;
SI_NATL    = SI.*NATL_srf;
O2min_NATL = O2min.*NATL_srf;

NF_SATL    = NF.*SATL_srf;
SI_SATL    = SI.*SATL_srf;
DFe_SATL   = DFe.*SATL_srf;
DIP_SATL   = DIP.*SATL_srf;
DIN_SATL   = DIN.*SATL_srf;
Pstar_SATL = Pstar.*SATL_srf;
SST_SATL   = SST.*SATL_srf;
SI_SATL    = SI.*SATL_srf;
O2min_SATL = O2min.*SATL_srf;

NF_NPAC    = NF.*NPAC_srf;
SI_NPAC    = SI.*NPAC_srf;
DFe_NPAC   = DFe.*NPAC_srf;
DIP_NPAC   = DIP.*NPAC_srf;
DIN_NPAC   = DIN.*NPAC_srf;
Pstar_NPAC = Pstar.*NPAC_srf;
SST_NPAC   = SST.*NPAC_srf;
SI_NPAC    = SI.*NPAC_srf;
O2min_NPAC = O2min.*NPAC_srf;

NF_SPAC    = NF.*SPAC_srf;
SI_SPAC    = SI.*SPAC_srf;
DFe_SPAC   = DFe.*SPAC_srf;
DIP_SPAC   = DIP.*SPAC_srf;
DIN_SPAC   = DIN.*SPAC_srf;
Pstar_SPAC = Pstar.*SPAC_srf;
SST_SPAC   = SST.*SPAC_srf;
SI_SPAC    = SI.*SPAC_srf;
O2min_SPAC = O2min.*SPAC_srf;

% =====================================================================
fprintf('\n\n')
fprintf('N. Atlantic ocean regression \n')
% select data only in the N. Atlantic ocean for regression. 
ido = find(NF_NATL~=0);
PAR_NATL = [DFe_NATL(ido),DIP_NATL(ido),DIN_NATL(ido),Pstar_NATL(ido),...
            SST_NATL(ido),SI_NATL(ido),O2min_NATL(ido)];
mdl_NATL = stepwiselm(PAR_NATL,NF_NATL(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)

% ======================================================================
fprintf('\n\n')
fprintf('S. Atlantic ocean regression \n')
% select data only in the S. Atlantic ocean for regression. 
ido = find(NF_SATL~=0);
RG_SATL  = fitlm(SI_SATL(ido),NF_SATL(ido));
PAR_SATL = [DFe_SATL(ido),DIP_SATL(ido),DIN_SATL(ido),Pstar_SATL(ido),...
            SST_SATL(ido),SI_SATL(ido),O2min_SATL(ido)];
mdl_SATL = stepwiselm(PAR_SATL,NF_SATL(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)

% =======================================================================
fprintf('\n\n')
fprintf('N. Pacific ocean regression \n')
% select data only in the N. Pacific ocean for regression. 
ido = find(NF_NPAC~=0);
PAR_NPAC = [DFe_NPAC(ido),DIP_NPAC(ido),DIN_NPAC(ido),Pstar_NPAC(ido),...
            SST_NPAC(ido),SI_NPAC(ido),O2min_NPAC(ido)];
mdl_NPAC = stepwiselm(PAR_NPAC,NF_NPAC(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)

% =======================================================================
fprintf('\n\n')
fprintf('S. Pacific ocean regression \n')
% select data only in the S. Pacific ocean for regression. 
ido = find(NF_SPAC~=0);
PAR_SPAC = [DFe_SPAC(ido),DIP_SPAC(ido),DIN_SPAC(ido),Pstar_SPAC(ido),...
            SST_SPAC(ido),SI_SPAC(ido),O2min_SPAC(ido)];
mdl_SPAC = stepwiselm(PAR_SPAC,NF_SPAC(ido),'quadratic','Criterion','Rsquared','PRemove',0.02)
