clc 
clear all
close all

load NF_wo_outliers.mat 
load ~/DATA/no3obs_90x180x24.mat 
load ~/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat

DIN = nanmean(no3obs(:,:,1:2),3);

iposi = and(NF2d(:)>0,DIN(:)>0); % use only positive numbers.

NF  = log10(NF2d(iposi));
DIN = log10(DIN(iposi));

NF(imag(NF)~=0) = nan;    % change imaginary number to nan.
DIN(imag(DIN)~=0) = nan;  % change imaginary number to nan.
% get rid of nan and infinity numbers.
ivalue = and(~or(isnan(NF),isnan(DIN)),~or(isinf(NF), isinf(DIN)));

NF  = NF(ivalue);
DIN = DIN(ivalue);

%NF_z  = (NF-mean(NF))/std(NF);
%DIN_z = (DIN-mean(DIN))/std(DIN);
%NF    = NF_z;
%DIN   = DIN_z;

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
NF_ATL = NF.*ATL_srf(ivalue);
DIN_ATL = DIN.*ATL_srf(ivalue);

NF_PAC = NF.*PAC_srf(ivalue);
DIN_PAC = DIN.*PAC_srf(ivalue);

NF_IND = NF.*IND_srf(ivalue);
DIN_IND = DIN.*IND_srf(ivalue);

NF_ARC = NF.*ARC_srf(ivalue);
DIN_ARC = DIN.*ARC_srf(ivalue);

% ==========================================================================
fprintf('global ocean regression \n')
% linear regression 
RG_GLOB = fitlm(DIN,NF,'linear');
fprintf('Linear Regression R^2 = %0.3f; \n',  RG_GLOB.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_GLOB.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_GLOB.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_GLOB.Coefficients.SE(2))

txt.Basion = 'Global';
txt.parm   = 'DIN';
txt.target = 'NF';
txt.target_longname = 'Nitrogen fixation of OCIM prediction';
txt.parm_longname = 'Surface water inorganic nitrogen concnetration';
txt.xunit  = texlabel(sprintf('(mmol m^{-3})'));
txt.yunit  = texlabel(sprintf('(mmol N m^{-2} y^{-1})')); 

% quadratic regression
figure()
plt_conf(DIN,NF,txt,1);

% ========================================================================
fprintf('Atlantic ocean regression \n')
% select data only in the Atlantic ocean for regression. 
ido = find(NF_ATL~=0);
RG_ATL = fitlm(DIN_ATL(ido),NF_ATL(ido));

fprintf('Linear Regression R^2 = %0.3f; \n',  RG_ATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_ATL.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_ATL.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_ATL.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'Atlantic';
plt_conf(DIN_ATL(ido),NF_ATL(ido),txt,1)
SB.ATL_x  = DIN_ATL(ido);
SB.ATL_y = NF_ATL(ido);

%==========================================================================
fprintf('Pacific ocean regression \n')
% select data only in the Pacific ocean for regression. 
ido = find(NF_PAC~=0);
RG_PAC = fitlm(DIN_PAC(ido),NF_PAC(ido));

fprintf('Linear Regression R^2 = %0.3f; \n',  RG_PAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_PAC.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_PAC.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_PAC.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'Pacific';
plt_conf(DIN_PAC(ido),NF_PAC(ido),txt,1)
SB.PAC_x  = DIN_PAC(ido);
SB.PAC_y  = NF_PAC(ido);

%===========================================================================
fprintf('Indian ocean regression \n')
% select data only in the Indian ocean for regression. 
ido = find(NF_IND~=0);
RG_IND = fitlm(DIN_IND(ido),NF_IND(ido));

fprintf('Linear Regression R^2 = %0.3f; \n',  RG_IND.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_IND.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_IND.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_IND.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'Indian';
plt_conf(DIN_IND(ido),NF_IND(ido),txt,1)
SB.IND_x  = DIN_IND(ido);
SB.IND_y  = NF_IND(ido);

%===========================================================================
fprintf('Arctic ocean regression \n')
% select data only in the Arctic ocean for regression. 
ido = find(NF_ARC~=0);
RG_ARC = fitlm(DIN_ARC(ido),NF_ARC(ido));

fprintf('Linear Regression R^2 = %0.3f; \n',  RG_ATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n \n \n',RG_ATL.Coefficients.pValue(2))

% quadratic regression
figure()
txt.Basion = 'Arctic';
plt_conf(DIN_ARC(ido),NF_ARC(ido),txt,1)
SB.ARC_x  = DIN_ARC(ido);
SB.ARC_y  = NF_ARC(ido);

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
NF_NATL = NF.*NATL_srf(ivalue);
DIN_NATL = DIN.*NATL_srf(ivalue);
NF_SATL = NF.*SATL_srf(ivalue);
DIN_SATL = DIN.*SATL_srf(ivalue);

NF_NPAC = NF.*NPAC_srf(ivalue);
DIN_NPAC = DIN.*NPAC_srf(ivalue);
NF_SPAC = NF.*SPAC_srf(ivalue);
DIN_SPAC = DIN.*SPAC_srf(ivalue);

% =====================================================================
fprintf('N. Atlantic ocean regression \n')
% select data only in the N. Atlantic ocean for regression. 
ido = find(NF_NATL~=0);
RG_NATL = fitlm(DIN_NATL(ido),NF_NATL(ido));
fprintf('Linear Regression R^2 = %0.3f; \n',  RG_NATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NATL.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_NATL.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_NATL.Coefficients.SE(2))

[p,ErrorEst] = polyfit(DIN_NATL(ido),NF_NATL(ido),1);
[NF_fit,delta] = polyval(p,DIN_NATL(ido),ErrorEst);

% quadratic regression
figure()
txt.Basion = 'N. Atlantic';
plt_conf(DIN_NATL(ido),NF_NATL(ido),txt,1)
SB.NATL_x  = DIN_NATL(ido);
SB.NATL_y  = NF_NATL(ido);

% ======================================================================
fprintf('S. Atlantic ocean regression \n')
% select data only in the S. Atlantic ocean for regression. 
ido = find(NF_SATL~=0);
RG_SATL = fitlm(DIN_SATL(ido),NF_SATL(ido));
fprintf('Linear Regression R^2 = %0.3f; \n',  RG_SATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SATL.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_SATL.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_SATL.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'S. Atlantic';
plt_conf(DIN_SATL(ido),NF_SATL(ido),txt,1)
SB.SATL_x  = DIN_SATL(ido);
SB.SATL_y  = NF_SATL(ido);

% =======================================================================
fprintf('N. Pacific ocean regression \n')
% select data only in the N. Pacific ocean for regression. 
ido = find(NF_NPAC~=0);
RG_NPAC = fitlm(DIN_NPAC(ido),NF_NPAC(ido));
fprintf('Linear Regression R^2 = %0.3f; \n',  RG_NPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NPAC.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_NPAC.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_NPAC.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'N. Pacific';
plt_conf(DIN_NPAC(ido),NF_NPAC(ido),txt,1)
SB.NPAC_x  = DIN_NPAC(ido);
SB.NPAC_y  = NF_NPAC(ido);

% =======================================================================
fprintf('S. Pacific ocean regression \n')
% select data only in the S. Pacific ocean for regression. 
ido = find(NF_SPAC~=0);
RG_SPAC = fitlm(DIN_SPAC(ido),NF_SPAC(ido));
fprintf('Linear Regression R^2 = %0.3f; \n',  RG_SPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SPAC.Coefficients.pValue(2))
fprintf('Slope   = %5.2e \n',RG_SPAC.Coefficients.Estimate(2))
fprintf('SlopeSE = %5.2e \n\n',RG_SPAC.Coefficients.SE(2))

% quadratic regression
figure()
txt.Basion = 'S. Pacific';
plt_conf(DIN_SPAC(ido),NF_SPAC(ido),txt,1)
SB.SPAC_x  = DIN_SPAC(ido);
SB.SPAC_y  = NF_SPAC(ido);

figure()
txt.Basion = 'Global';
plt_conf(1,1,txt,1,SB)

